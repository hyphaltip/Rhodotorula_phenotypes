#!/usr/bin/env bash
# LOCO (leave-one-chromosome-out) sensitivity check, 3 traits (chroma, AUC_10,
# resilience_30 -- matching the ORIGINAL run's actual scope: PROGRESS.md prose
# says "6 traits" but the real output directory only has these 3, x 2 panel
# sets = 6 (trait,panel) combinations, x ~20 chromosomes = 120 scans) x both
# panels (gwas=213, gwasc=182-culled).
#
# Reuses run_loco_shared.sh + plink_arch.sh (both ported verbatim from
# analysis/ideas/2026-08-15-color-phenotype-space/results/gwas/loco/ -- this
# machinery WAS saved as a script originally, unlike the culling algorithm).
#
# Run from the repo root, inside a SLURM job (needs $SCRATCH):
#   sbatch --wrap="bash analysis/gwas/scripts/run_loco.sh"
set -euo pipefail

REPO_ROOT="$PWD"
[ -f "$REPO_ROOT/pixi.toml" ] || { echo "ERROR: run from repo root." >&2; exit 1; }

: "${SCRATCH:?SCRATCH env var not set -- run inside a SLURM job}"
WORK="$SCRATCH/gwas_loco"
mkdir -p "$WORK"
cp analysis/gwas/scripts/plink_arch.sh "$WORK/plink_arch.sh"

PHENO_CSV="analysis/gwas/results/gwas_next_phenotypes_fam_order.csv"
OUTDIR="analysis/gwas/results/gwas/loco"
mkdir -p "$OUTDIR/output"

TRAITS="chroma AUC_10 resilience_30"

run_panel() {
  local panel="$1" bfile="$2" kinroot="$3"
  for trait in $TRAITS; do
    echo "=== LOCO [$panel]: $trait ==="
    TWORK="$WORK/${panel}_${trait}"
    mkdir -p "$TWORK"
    cp "${bfile}.bed" "$TWORK/g.bed"
    cp "${bfile}.bim" "$TWORK/g.bim"
    pixi run python3 -c "
import pandas as pd
fam = pd.read_csv('${bfile}.fam', sep=r'\s+', header=None, names=['fid','iid','pid','mid','sex','pheno'])
pheno = pd.read_csv('$PHENO_CSV')
pheno_map = dict(zip(pheno['strain_code'], pheno['$trait']))
vals = fam['iid'].map(pheno_map)
fam['pheno'] = vals.apply(lambda v: 'NA' if pd.isna(v) else v)
fam.to_csv('$TWORK/g.fam', sep=' ', header=False, index=False)
"
    cp "${kinroot}.bed" "$TWORK/kinroot.bed"
    cp "${kinroot}.bim" "$TWORK/kinroot.bim"
    cp "${kinroot}.fam" "$TWORK/kinroot.fam"

    bash analysis/gwas/scripts/run_loco_shared.sh \
      "$trait" "$TWORK/g" "$TWORK/kinroot" "${panel}" "1" "" "$WORK"

    mkdir -p "$OUTDIR/output"
    cp "$WORK/output/loco_${panel}_${trait}_"*.assoc.txt "$OUTDIR/output/" 2>/dev/null || true
    rm -rf "$TWORK"
  done
}

run_panel gwas \
  "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas" \
  "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas.pruned"
run_panel gwasc \
  "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship_culled/gwasc" \
  "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship_culled/gwasc.pruned"

echo "LOCO complete: 3 traits x 2 panels x ~20-23 chromosomes."
