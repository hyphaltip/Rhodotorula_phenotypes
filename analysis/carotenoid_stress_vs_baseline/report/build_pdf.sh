#!/usr/bin/env bash
# Build REPORT.pdf from REPORT.md (images are referenced relative to this directory). Run from the repo root inside a SLURM job.
set -euo pipefail
module load texlive
export PATH="/opt/linux/rocky/8.x/x86_64/pkgs/texlive/20220403/bin/x86_64-linux:$PATH"   # ~/bin holds an older partial TeX Live that would otherwise win
export TEXMFHOME=/nonexistent   # ~/texmf holds an old xcolor.sty that breaks TeX Live 2022
cd analysis/carotenoid_stress_vs_baseline/report
pandoc REPORT.md -o REPORT.pdf --pdf-engine=pdflatex --resource-path=. -V geometry:margin=1.6cm -V fontsize=9pt -V colorlinks=true \
  --toc --toc-depth=2 -V title="a* (carotenoid proxy) under metal stress: stratified analyses"
