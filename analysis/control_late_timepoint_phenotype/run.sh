#!/usr/bin/env bash
# Reproduce the control-only late-timepoint strain phenotype table.
# Run from the repo root: bash analysis/control_late_timepoint_phenotype/run.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

pixi run python analysis/control_late_timepoint_phenotype/scripts/build_phenotype_table.py

echo "Done. Table: analysis/control_late_timepoint_phenotype/results/phenotype_control_late_timepoint.csv"
