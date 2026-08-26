#!/usr/bin/env bash
# Tier A: per-trait kinship-only GEMMA LMM scan on the rebuilt 213-strain kinship
# (Task 7's rebuild_tiers_abc.sh output). Reconstructed from
# analysis/ideas/2026-08-15-color-phenotype-space/PROGRESS.md section 6 (N4/N5) --
# those commands ran ad hoc on $SCRATCH and were never saved as a script.
#
# GEMMA reads the phenotype from .fam column 6, NOT from -p, when both -bfile and
# -p are given (learning L-16) -- so this bakes each trait into its own per-trait
# .fam copy rather than passing -p.
#
# Run from the repo root: bash analysis/gwas/scripts/run_tiera_gemma.sh
set -euo pipefail

REPO_ROOT="$PWD"
if [ ! -f "$REPO_ROOT/pixi.toml" ]; then
  echo "ERROR: run this from the repo root." >&2
  exit 1
fi

: "${SCRATCH:?SCRATCH env var not set -- run this inside a SLURM job}"
WORK="$SCRATCH/gwas_tiera"
mkdir -p "$WORK"

KIN_DIR="analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship"
BED="$KIN_DIR/gwas.pruned.bed"
BIM="$KIN_DIR/gwas.pruned.bim"
FAM="$KIN_DIR/gwas.pruned.fam"
KINSHIP="$KIN_DIR/kins.cXX.txt"
PHENO_CSV="analysis/gwas/results/gwas_next_phenotypes_fam_order.csv"
TOOLCHAIN="/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/GWAS/vc2gwas_env/bin"

for f in "$BED" "$BIM" "$FAM" "$KINSHIP" "$PHENO_CSV"; do
  [ -f "$f" ] || { echo "ERROR: required input missing: $f" >&2; exit 1; }
done

OUTDIR="analysis/gwas/results/gwas/tierA_summary"
mkdir -p "$OUTDIR"

TRAITS="chroma sat bright clone_mean_area AUC_0 AUC_10 AUC_20 AUC_30 AUC_ratio_10 resilience_30 cu_dose_slope IC50_est"

for trait in $TRAITS; do
  echo "=== Tier A: $trait ==="
  TWORK="$WORK/$trait"
  mkdir -p "$TWORK"
  cp "$BED" "$TWORK/gwas.bed"
  cp "$BIM" "$TWORK/gwas.bim"

  # Bake this trait's value into .fam col 6 (missing -> NA, GEMMA's missing marker).
  pixi run python3 -c "
import pandas as pd
fam = pd.read_csv('$FAM', sep=r'\s+', header=None,
                   names=['fid','iid','pid','mid','sex','pheno'])
pheno = pd.read_csv('$PHENO_CSV')
pheno_map = dict(zip(pheno['strain_code'], pheno['$trait']))
vals = fam['iid'].map(pheno_map)
n_missing = vals.isna().sum()
fam['pheno'] = vals.apply(lambda v: 'NA' if pd.isna(v) else v)
fam.to_csv('$TWORK/gwas.fam', sep=' ', header=False, index=False)
print(f'  $trait: {len(fam) - n_missing}/{len(fam)} strains with phenotype, {n_missing} NA')
"

  "$TOOLCHAIN/gemma" -bfile "$TWORK/gwas" -k "$KINSHIP" -lmm 4 \
    -o "gwas_${trait}" -outdir "$OUTDIR/gemma_output" \
    > "$OUTDIR/gemma_output_${trait}.log" 2>&1 \
    || { echo "ERROR: GEMMA failed for $trait -- see $OUTDIR/gemma_output_${trait}.log" >&2; exit 1; }
  echo "  wrote $OUTDIR/gemma_output/gwas_${trait}.assoc.txt"
done

echo "Tier A scan complete: 12 traits x 213 strains. Assoc files in $OUTDIR/gemma_output/."
