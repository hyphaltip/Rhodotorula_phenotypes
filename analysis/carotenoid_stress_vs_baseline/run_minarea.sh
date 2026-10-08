#!/usr/bin/env bash
# Sensitivity to a minimum colony area (pixels). Run from the repo root as a SLURM job.
set -euo pipefail
export SKIP_ZINC=1
S=analysis/carotenoid_stress_vs_baseline/scripts; R=analysis/carotenoid_stress_vs_baseline/results
for A in 1000 2000 3000 5000; do
  pixi run python $S/prepare_wells.py --min-area $A
  mkdir -p $R/minarea_$A
  pixi run Rscript $S/fit_mixed_models.R $R/wells_min$A.csv $R/minarea_$A
done
