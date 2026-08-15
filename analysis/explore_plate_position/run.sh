#!/usr/bin/env bash
# Reproduce the full plate-position analysis end to end.
# Run from the repo root: bash analysis/explore_plate_position/run.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

pixi run Rscript analysis/explore_plate_position/scripts/00_build_dataset.R
pixi run Rscript analysis/explore_plate_position/scripts/01_variance_partition.R
pixi run Rscript analysis/explore_plate_position/scripts/02_adjacency_effect.R
pixi run Rscript analysis/explore_plate_position/scripts/03_plots.R
# Render HTML + GitHub-flavored markdown + PDF. The PDF needs a complete TeX
# distribution: the default pdflatex on PATH is a minimal TinyTeX missing
# common packages (e.g. xcolor.sty), so prefer the full HPCC TeX Live 2022
# install when present. options(tinytex.tlmgr.path = "") stops the tinytex R
# package from prepending its own minimal TinyTeX bin back onto PATH.
TEX_BIN="/opt/linux/rocky/8.x/x86_64/pkgs/texlive/20220403/bin/x86_64-linux"
if [ -x "$TEX_BIN/pdflatex" ]; then
  export PATH="$TEX_BIN:$PATH"
fi
pixi run Rscript -e 'options(tinytex.tlmgr.path = ""); rmarkdown::render("analysis/explore_plate_position/explore_plate_position.Rmd", output_format = "all")'

echo "Done. Report: analysis/explore_plate_position/explore_plate_position.{html,pdf,md} (+ figure/ folder for GitHub)"
