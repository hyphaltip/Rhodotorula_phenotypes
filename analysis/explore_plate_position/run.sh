#!/usr/bin/env bash
# Reproduce the full plate-position analysis end to end.
# Run from the repo root: bash analysis/explore_plate_position/run.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

pixi run Rscript analysis/explore_plate_position/scripts/00_build_dataset.R
pixi run Rscript analysis/explore_plate_position/scripts/01_variance_partition.R
pixi run Rscript analysis/explore_plate_position/scripts/02_adjacency_effect.R
pixi run Rscript analysis/explore_plate_position/scripts/03_plots.R
pixi run Rscript -e 'rmarkdown::render("analysis/explore_plate_position/explore_plate_position.Rmd")'

echo "Done. Report: analysis/explore_plate_position/explore_plate_position.html"
