#!/usr/bin/env bash
# Genome-wide KEGG Orthology (KO) assignment via KofamScan -- the foundation layer of
# this repo's functional-annotation practice (per user request 2026-08-26), letting any
# future gene query check KEGG pathway membership (e.g. map00906 carotenoid
# biosynthesis) directly, not just via HMM/keyword matching against one curated library.
#
# Invocation pattern matches the lab's existing usage
# (/bigdata/stajichlab/shared/projects/Black_yeasts/KEGG_KOFAM_Scan/pipeline/03_KOALA_scan.sh):
# exec_annotation with the eukaryote profile subset (eukaryote.hal, ~15,700 of the full
# 27,233 KO profiles -- restricting to eukaryote-relevant KOs cuts runtime substantially
# without losing recall for a fungal proteome) and the standard E<=1e-4 KofamScan default.
#
# Usage: sbatch --partition=stajichlab --mem=64G --time=12:00:00 --cpus-per-task=8 \
#          --wrap="bash run_kofamscan.sh <proteome.fa> <out_prefix>"
set -euo pipefail

PROTEOME="${1:?need a proteome FASTA}"
OUT_PREFIX="${2:?need an output prefix}"
source /etc/profile.d/modules.sh
module load kofamscan
module load workspace/scratch

KOFAM_DB=/bigdata/stajichlab/shared/lib/gtotree/kofamscan_data
PROFILES="$KOFAM_DB/profiles/eukaryote.hal"
KOLIST="$KOFAM_DB/ko_list"
CPU="${SLURM_CPUS_ON_NODE:-8}"

mkdir -p "$(dirname "$OUT_PREFIX")"
: "${SCRATCH:?SCRATCH env var not set -- run inside a SLURM job}"
mkdir -p "$SCRATCH/tmp"

exec_annotation "$PROTEOME" -p "$PROFILES" -k "$KOLIST" \
  -o "${OUT_PREFIX}.detail.txt" -E 0.0001 --cpu "$CPU" --tmp-dir "$SCRATCH/tmp"
exec_annotation "$PROTEOME" -p "$PROFILES" -k "$KOLIST" \
  -o "${OUT_PREFIX}.mapper.txt" -f mapper -E 0.0001 --cpu "$CPU" --tmp-dir "$SCRATCH/tmp2"

echo "KofamScan complete: ${OUT_PREFIX}.detail.txt (annotated hits, * = best), ${OUT_PREFIX}.mapper.txt (protein<TAB>KO)"
