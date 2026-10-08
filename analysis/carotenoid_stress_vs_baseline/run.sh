#!/usr/bin/env bash
# Run from the repo root as a SLURM job:
#   sbatch -p short -c 4 --mem=24G -t 60 --wrap="cd $PWD && bash analysis/carotenoid_stress_vs_baseline/run.sh"
set -euo pipefail
S=analysis/carotenoid_stress_vs_baseline/scripts
pixi run python $S/prepare_wells.py
pixi run python $S/prepare_wells.py --all-objects          # sensitivity table; compare with: pixi run Rscript $S/fit_mixed_models.R <wells_allobj.csv> <outdir>
pixi run Rscript $S/fit_mixed_models.R
pixi run python $S/plot_figures.py
