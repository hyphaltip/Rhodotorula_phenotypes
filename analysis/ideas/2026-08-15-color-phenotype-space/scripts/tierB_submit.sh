#!/bin/bash
#SBATCH --job-name=tierB_gwas
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=2:00:00
#SBATCH --output=slurm_tierB_gwas_%j.out
#SBATCH --error=slurm_tierB_gwas_%j.err
#SBATCH --nodes=1
set -euo pipefail

HERE=/bigdata/stajichlab/jstajich/projects/Rhodotorula_phenotypes/analysis/ideas/2026-08-15-color-phenotype-space
PY=/bigdata/stajichlab/jstajich/projects/Rhodotorula_phenotypes/.pixi/envs/default/bin/python3

cd "$HERE"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
"$PY" scripts/tierb_set_tests.py \
  --assoc-dir results/gwas/tierA_summary \
  --pixy-dir results/gwas/pixy \
  --bfile-base /scratch/jstajich/27384933/gwas/work/gwas \
  --work /scratch/jstajich/27384933/gwas/work \
  --out results/gwas/tierB \
  --prefix gwas --mc 50000 \
  > tierb_run_gwas.log 2>&1
echo "TIERB_DONE rc=$?" >> tierb_run_gwas.log
