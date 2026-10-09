#!/usr/bin/env bash
# Run from the repo root as a SLURM job:
#   sbatch -p short -c 4 --mem=24G -t 60 --wrap="cd $PWD && bash analysis/carotenoid_stress_vs_baseline/run.sh"
set -euo pipefail
S=analysis/carotenoid_stress_vs_baseline/scripts
pixi run python $S/prepare_wells.py --min-area 2000      # main dataset: area >= 2000 px, Zinc window ends 80 h
cp analysis/carotenoid_stress_vs_baseline/results/wells_min2000.csv analysis/carotenoid_stress_vs_baseline/results/wells.csv
pixi run python $S/prepare_wells.py --min-area 0         # unfiltered, for the sensitivity section
pixi run python $S/prepare_wells.py --all-objects          # sensitivity table; compare with: pixi run Rscript $S/fit_mixed_models.R <wells_allobj.csv> <outdir>
pixi run Rscript $S/fit_mixed_models.R
pixi run python $S/plot_figures.py
