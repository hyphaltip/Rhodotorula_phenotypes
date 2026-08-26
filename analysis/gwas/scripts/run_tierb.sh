#!/usr/bin/env bash
# Tier B: SKAT/burden/min-p set tests on pixy high-dxy windows, both panels
# (gwas=213, gwasc=182-culled). Ported tierb_set_tests.py (unmodified logic,
# analysis/ideas/2026-08-15-color-phenotype-space/scripts/tierb_set_tests.py)
# against this port's rebuilt Tier A assoc files and kinship's full genotype
# bfile (LD for windows computed from the FULL unpruned set, matching the
# original invocation in tierB_submit.sh).
#
# CAVEAT (documented in GWAS.md): --pixy-dir points at the ORIGINAL run's pixy
# output, computed on the prior 201-strain cohort.all.vcf.gz -- pixy was NOT
# recomputed for the 213-strain panel (a 6h35m job) in this session. The
# high-dxy WINDOW DEFINITIONS this selects are a population-genetic property
# of the genome, not expected to be strongly sensitive to 12 more strains
# within existing populations, but this is an approximation, not a rebuild.
#
# Run from the repo root, inside a SLURM job (needs $SCRATCH-scale memory,
# ~40G, per the original tierB_submit.sh):
#   sbatch --wrap="bash analysis/gwas/scripts/run_tierb.sh"
set -euo pipefail

REPO_ROOT="$PWD"
[ -f "$REPO_ROOT/pixi.toml" ] || { echo "ERROR: run from repo root." >&2; exit 1; }

ASSOC_CSV_DIR="analysis/gwas/results/gwas/tierA_summary/assoc_csv"
PIXY_DIR="analysis/ideas/2026-08-15-color-phenotype-space/results/gwas/pixy"
WORK_ARCH="analysis/gwas/scripts"   # holds plink_arch.sh
OUTDIR="analysis/gwas/results/gwas/tierB"

for panel_bfile in \
  "gwas:analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas" \
  "gwasc:analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship_culled/gwasc"
do
  panel="${panel_bfile%%:*}"
  bfile="${panel_bfile#*:}"
  [ -f "${bfile}.bed" ] || { echo "ERROR: missing ${bfile}.bed -- run rebuild_full_genotypes_and_tiera.sh first" >&2; exit 1; }
  echo "=== Tier B: $panel ==="
  export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
  pixi run python3 analysis/gwas/scripts/tierb_set_tests.py \
    --assoc-dir "$ASSOC_CSV_DIR" \
    --pixy-dir "$PIXY_DIR" \
    --bfile-base "$bfile" \
    --work "$WORK_ARCH" \
    --out "$OUTDIR" \
    --prefix "$panel" --mc 50000 \
    > "$OUTDIR/tierb_run_${panel}.log" 2>&1 \
    || { echo "ERROR: Tier B failed for $panel -- see $OUTDIR/tierb_run_${panel}.log" >&2; exit 1; }
  echo "  wrote $OUTDIR/tierb_settests_${panel}.csv"
done

echo "Tier B complete for both panels."
