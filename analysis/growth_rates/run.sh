#!/usr/bin/env bash
# Reproduce the full growth-rate analysis end to end.
# Run from the repo root: bash analysis/growth_rates/run.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

pixi run Rscript analysis/growth_rates/scripts/00_build_series.R
pixi run Rscript analysis/growth_rates/scripts/01_fit_growth_models.R
pixi run Rscript analysis/growth_rates/scripts/02_species_cu_rates.R
pixi run Rscript analysis/growth_rates/scripts/03_color_interaction.R
pixi run Rscript analysis/growth_rates/scripts/04_doubling_time.R
pixi run Rscript analysis/growth_rates/scripts/05_species_cu_sensitivity.R
# Render HTML + GitHub-flavored markdown + PDF. Prefer the full HPCC TeX Live
# 2022 install over the minimal TinyTeX on PATH; options(tinytex.tlmgr.path="")
# stops tinytex from prepending its minimal TeX bin back onto PATH.
TEX_BIN="/opt/linux/rocky/8.x/x86_64/pkgs/texlive/20220403/bin/x86_64-linux"
if [ -x "$TEX_BIN/pdflatex" ]; then
  export PATH="$TEX_BIN:$PATH"
fi
pixi run Rscript -e 'options(tinytex.tlmgr.path = ""); rmarkdown::render("analysis/growth_rates/growth_rates.Rmd", output_format = "all")'

echo "Done. Report: analysis/growth_rates/growth_rates.{html,pdf,md}"
