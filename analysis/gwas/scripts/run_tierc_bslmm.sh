#!/usr/bin/env bash
# Tier C: BSLMM (Bayesian sparse linear mixed model) on 5 representative traits,
# all-213 panel only (matches the prior run's scope: tierC_summary/ there has
# bslmm_{chroma,AUC_10,AUC_30,clone_mean_area,resilience_30} and nothing for the
# culled set). GEMMA's `-bslmm 1` estimates its own relatedness structure from
# genotypes directly -- it does not take a `-k` kinship argument.
#
# MCMC params reconstructed from PROGRESS.md section 7 ("100k MCMC, 20% burn-in
# discard") and the prior run's .hyp.txt row count (100,000 recorded samples,
# i.e. "100k MCMC" = 100k RETAINED samples, not 100k total iterations). GEMMA's
# default -rpace (record pace) is 10, so -s 1000000 records exactly 100,000
# rows (matches the original file's row count precisely, verified directly:
# `zcat bslmm_chroma.hyp.txt.gz | wc -l` = 100001 in both the original run and
# this reconstruction). -w 200000 (burn-in, ~20% of the 1,000,000 sampling
# phase, matching the "20% burn-in" prose) is discarded, not recorded.
# An earlier attempt at this script used -w 20000 -s 100000 (a literal but
# wrong reading of "100k MCMC" as total iterations) -- that produced a 10x
# SHORTER chain (10,000 recorded rows), caught by directly diffing row counts
# against the original before trusting the posterior estimates.
#
# Run from the repo root, inside a SLURM job:
#   sbatch --wrap="bash analysis/gwas/scripts/run_tierc_bslmm.sh"
set -euo pipefail

REPO_ROOT="$PWD"
[ -f "$REPO_ROOT/pixi.toml" ] || { echo "ERROR: run from repo root." >&2; exit 1; }

: "${SCRATCH:?SCRATCH env var not set -- run inside a SLURM job}"
WORK="$SCRATCH/gwas_tierc"
mkdir -p "$WORK"

BFILE_BASE="analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas"
PHENO_CSV="analysis/gwas/results/gwas_next_phenotypes_fam_order.csv"
TOOLCHAIN="/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/GWAS/vc2gwas_env/bin"
OUTDIR="analysis/gwas/results/gwas/tierC_summary"
mkdir -p "$OUTDIR"

for f in "$BFILE_BASE.bed" "$BFILE_BASE.bim" "$BFILE_BASE.fam" "$PHENO_CSV"; do
  [ -f "$f" ] || { echo "ERROR: required input missing: $f" >&2; exit 1; }
done

TRAITS="chroma AUC_10 AUC_30 clone_mean_area resilience_30"

for trait in $TRAITS; do
  echo "=== Tier C BSLMM: $trait ==="
  TWORK="$WORK/$trait"
  mkdir -p "$TWORK"
  cp "$BFILE_BASE.bed" "$TWORK/g.bed"
  cp "$BFILE_BASE.bim" "$TWORK/g.bim"
  pixi run python3 -c "
import pandas as pd
fam = pd.read_csv('$BFILE_BASE.fam', sep=r'\s+', header=None, names=['fid','iid','pid','mid','sex','pheno'])
pheno = pd.read_csv('$PHENO_CSV')
pheno_map = dict(zip(pheno['strain_code'], pheno['$trait']))
vals = fam['iid'].map(pheno_map)
n_missing = vals.isna().sum()
fam['pheno'] = vals.apply(lambda v: 'NA' if pd.isna(v) else v)
fam.to_csv('$TWORK/g.fam', sep=' ', header=False, index=False)
print(f'  $trait: {len(fam) - n_missing}/{len(fam)} with phenotype, {n_missing} NA')
"
  "$TOOLCHAIN/gemma" -bfile "$TWORK/g" -bslmm 1 -w 200000 -s 1000000 \
    -o "bslmm_${trait}" -outdir "$OUTDIR" \
    > "$OUTDIR/bslmm_${trait}.runlog.txt" 2>&1 \
    || { echo "ERROR: BSLMM failed for $trait -- see $OUTDIR/bslmm_${trait}.runlog.txt" >&2; exit 1; }
  gzip -f "$OUTDIR/bslmm_${trait}.hyp.txt" "$OUTDIR/bslmm_${trait}.gamma.txt" 2>/dev/null || true
  rm -rf "$TWORK"
  echo "  wrote $OUTDIR/bslmm_${trait}.{hyp,gamma,param}.txt(.gz)"
done

echo "Tier C BSLMM complete: 5 traits x 213 strains."
