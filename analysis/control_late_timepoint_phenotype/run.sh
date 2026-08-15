#!/usr/bin/env bash
# Reproduce the control-only late-timepoint strain phenotype tables (3 windows).
# Run from the repo root: bash analysis/control_late_timepoint_phenotype/run.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

S=analysis/control_late_timepoint_phenotype/scripts/build_phenotype_table.py
pixi run python "$S" --tmin 70 --tmax 80
pixi run python "$S" --tmin 80 --tmax 90
pixi run python "$S" --tmin 90 --tmax 110

echo "Done. Tables: analysis/control_late_timepoint_phenotype/results/phenotype_control_timepoint_{70_80,80_90,90_110}.csv"
