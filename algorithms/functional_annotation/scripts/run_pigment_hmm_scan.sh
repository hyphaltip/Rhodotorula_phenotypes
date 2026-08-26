#!/usr/bin/env bash
# Scan a proteome against the pigmentation-pathway HMM library
# (data/raw/pigmentation-pathway-hmms/) -- 28 profiles across 7 pigment biosynthesis
# pathways (fungal + cyanobacterial), imported from an external mining project
# (data/DATA_MANIFEST.md "pigmentation-pathway-hmms" entry).
#
# t3hnr/t4hnr use their DISCRIMINATIVE variants instead of the originals, per
# data/raw/pigmentation-pathway-hmms/specific/RECOMMENDED_CONFIG.md (28x/5x fewer false
# positives against the broad bacterial SDR superfamily): t4hnr_disc_top15_f3.hmm at
# E<=1.0, t3hnr_disc_top20_f4.hmm at E<=1e-5. The other 26 profiles have no calibrated
# per-profile threshold in this repo (the source project's search_genomes.py pipeline,
# which had them, was not copied over -- only its output HMMs) -- run at a permissive
# E<=1e-5 and report full-sequence + best-domain E-values for manual review, per the
# report's own guidance ("full sequence score used as primary criterion").
#
# Usage: bash run_pigment_hmm_scan.sh <proteome.fa> <out_dir>
set -euo pipefail

PROTEOME="${1:?need a proteome FASTA}"
OUTDIR="${2:?need an output directory}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HMM_DIR="$HERE/../../../data/raw/pigmentation-pathway-hmms"

mkdir -p "$OUTDIR"
# NOTE: `module load ... | grep ...` runs module load in a subshell, so the
# environment changes never reach the parent shell -- must be a bare statement.
module load hmmer/3.4

# --- Main 28-profile library (excludes t3hnr/t4hnr originals -- superseded below) ---
COMBINED="$OUTDIR/pigmentation_profiles.hmm"
cp "$HMM_DIR/pigmentation_profiles.hmm" "$COMBINED"
hmmpress -f "$COMBINED"

hmmscan --domtblout "$OUTDIR/pigment_scan.domtbl" -E 1e-5 --cpu "${SLURM_CPUS_ON_NODE:-4}" \
  "$COMBINED" "$PROTEOME" > "$OUTDIR/pigment_scan.hmmscan.log"

# --- t3hnr/t4hnr discriminative variants (RECOMMENDED_CONFIG.md thresholds) ---
for pair in "t4hnr_disc_top15_f3.hmm:1.0" "t3hnr_disc_top20_f4.hmm:1e-5"; do
  hmmfile="${pair%%:*}"; evalue="${pair##*:}"
  src="$HMM_DIR/specific/$hmmfile"
  if [ ! -f "$src" ]; then
    echo "WARNING: discriminative HMM not found: $src -- check data/raw/pigmentation-pathway-hmms/specific/ for the current filename" >&2
    continue
  fi
  cp "$src" "$OUTDIR/$hmmfile"
  hmmpress -f "$OUTDIR/$hmmfile"
  hmmscan --domtblout "$OUTDIR/${hmmfile%.hmm}.domtbl" -E "$evalue" --cpu "${SLURM_CPUS_ON_NODE:-4}" \
    "$OUTDIR/$hmmfile" "$PROTEOME" > "$OUTDIR/${hmmfile%.hmm}.hmmscan.log"
done

# Clean up regenerable intermediates -- only the small .domtbl files (the actual
# parseable output) are worth keeping; the copied .hmm files, pressed .h3* indexes, and
# verbose per-domain-alignment .hmmscan.log text (~11MB each) are not.
rm -f "$OUTDIR"/*.hmm "$OUTDIR"/*.h3f "$OUTDIR"/*.h3i "$OUTDIR"/*.h3m "$OUTDIR"/*.h3p "$OUTDIR"/*.hmmscan.log

echo "Pigment HMM scan complete. Outputs in $OUTDIR/*.domtbl"
