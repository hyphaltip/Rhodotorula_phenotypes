#!/usr/bin/env bash
# Tier A add-on: raw CIELAB components (L*, a*, b* medians, clone-mean over
# plates) as 3 additional traits, requested 2026-08-26 after the initial 12-trait
# Tier A rebuild only covered derived chroma/sat/bright. Same full unpruned
# SNP set + per-panel kinship as run_tiera_gemma.sh / resume_tiera_gemma.sh --
# skip-if-present so this is safe to rerun.
#
# Run from the repo root, inside a SLURM job (needs $SCRATCH):
#   sbatch --partition=stajichlab --time=3:00:00 --wrap="bash analysis/gwas/scripts/run_tiera_lab_traits.sh"
set -euo pipefail

REPO_ROOT="$PWD"
[ -f "$REPO_ROOT/pixi.toml" ] || { echo "ERROR: run from repo root." >&2; exit 1; }

: "${SCRATCH:?SCRATCH env var not set -- run inside a SLURM job}"
WORK="$SCRATCH/gwas_tiera_lab"
mkdir -p "$WORK"

PHENO_CSV="analysis/gwas/results/gwas_next_phenotypes_fam_order.csv"
TOOLCHAIN="/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/GWAS/vc2gwas_env/bin"
OUTDIR="analysis/gwas/results/gwas/tierA_summary"
mkdir -p "$OUTDIR/gemma_output"

TRAITS="lab_L lab_a lab_b"

run_tiera_panel() {
  local panel="$1" bfile="$2" kinship="$3"
  for trait in $TRAITS; do
    outfile="$OUTDIR/gemma_output/${panel}_${trait}.assoc.txt"
    if [ -s "$outfile" ]; then
      echo "skip [$panel/$trait]: already present ($outfile)"
      continue
    fi
    echo "=== Tier A [$panel]: $trait ==="
    TWORK="$WORK/${panel}_${trait}"
    mkdir -p "$TWORK"
    cp "$bfile.bed" "$TWORK/g.bed"
    cp "$bfile.bim" "$TWORK/g.bim"
    pixi run python3 -c "
import pandas as pd
fam = pd.read_csv('$bfile.fam', sep=r'\s+', header=None, names=['fid','iid','pid','mid','sex','pheno'])
pheno = pd.read_csv('$PHENO_CSV')
pheno_map = dict(zip(pheno['strain_code'], pheno['$trait']))
vals = fam['iid'].map(pheno_map)
n_missing = vals.isna().sum()
fam['pheno'] = vals.apply(lambda v: 'NA' if pd.isna(v) else v)
fam.to_csv('$TWORK/g.fam', sep=' ', header=False, index=False)
print(f'  [$panel] $trait: {len(fam) - n_missing}/{len(fam)} with phenotype, {n_missing} NA')
"
    "$TOOLCHAIN/gemma" -bfile "$TWORK/g" -k "$kinship" -lmm 4 \
      -o "${panel}_${trait}" -outdir "$OUTDIR/gemma_output" \
      > "$OUTDIR/gemma_output_${panel}_${trait}.log" 2>&1 \
      || { echo "ERROR: GEMMA failed for $panel/$trait" >&2; exit 1; }
    rm -rf "$TWORK"
  done
}

run_tiera_panel gwas \
  "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas" \
  "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/kins.cXX.txt"
run_tiera_panel gwasc \
  "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship_culled/gwasc" \
  "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship_culled/kinsc.cXX.txt"

echo "Lab-trait Tier A scan complete: 3 traits x 2 panels."
