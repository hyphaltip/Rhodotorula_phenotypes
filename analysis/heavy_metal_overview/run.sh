#!/usr/bin/env bash
# Rebuild the heavy-metal overview figures. Run from the repo root:
#   sbatch --wait -p short -c 4 --mem=24G -t 60 -o analysis/heavy_metal_overview/results/logs/make_overview.log \
#     --wrap="cd $PWD && bash analysis/heavy_metal_overview/run.sh"
# Do not run on the head node. No BASH_SOURCE (breaks on SLURM).
set -euo pipefail
cd "${HM_REPO_ROOT:-$PWD}"
test -f analysis/heavy_metal_overview/run.sh || { echo "run from repo root or set HM_REPO_ROOT" >&2; exit 1; }
pixi run python analysis/heavy_metal_overview/scripts/make_overview.py
