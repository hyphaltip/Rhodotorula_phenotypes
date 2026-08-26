#!/usr/bin/env bash
# Screen the (separately-called, excluded-from-this-pipeline's-sequence-construction)
# INDEL VCF for any indel falling within a candidate gene's CDS +/- flank -- per Opus's
# review: the design's "SNP-only VCF" claim only holds for the file chosen, and this
# same panel has a real INDEL VCF (data/raw/genotypes/RmucY2510_v2/*.INDEL.*) that could
# contain a real coding indel invisible to the substitution-based pipeline. This does
# NOT build indels into any sequence -- it only reports presence/absence per gene so a
# reader knows whether "no indels considered" is actually true for that gene, or is a
# known gap.
#
# Usage: bash screen_indels.sh <gene_id> <scaffold> <cds_start> <cds_end> <flank_bp> <out_csv>
set -euo pipefail
GENE="${1:?gene id}"; SCAFFOLD="${2:?scaffold}"; START="${3:?start}"; END="${4:?end}"
FLANK="${5:?flank bp}"; OUT="${6:?out csv}"

INDEL_VCF="data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.INDEL.combined_selected.vcf.gz"
REGION="${SCAFFOLD}:$((START - FLANK))-$((END + FLANK))"

# Fail loudly if the toolchain isn't available: a silent 0 here means "no indels",
# which is the exact wrong answer if bcftools simply isn't on PATH (this bit the pilot
# -- discovered 2026-08-26: pilot genes reported 0 indels when the real INDEL VCF has
# hundreds; the pipe + `|| true` swallowed the failure). Check the binary and that the
# VCF/tabix are present and readable BEFORE counting.
command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not on PATH -- can't screen indels." >&2; exit 1; }
[ -s "$INDEL_VCF" ] || { echo "ERROR: INDEL VCF missing: $INDEL_VCF" >&2; exit 1; }
if grep -q '^[1-9]' <(bcftools view -r "$REGION" -H "$INDEL_VCF" 2>/dev/null); then :; fi
RC=$(bcftools view -r "$REGION" -H "$INDEL_VCF" 2>/dev/null | wc -l)
# Sanity: bcftools failing (e.g. missing tabix index) produces empty output = RC 0.
[ "$RC" -gt 0 ] || { echo "ERROR: bcftools returned 0 records for $REGION -- likely missing tabix index or tool failure, refusing to report 0 indels." >&2; exit 1; }

echo "gene,scaffold,region,n_indels" > "$OUT.header"
echo "${GENE},${SCAFFOLD},${REGION},${RC}" >> "$OUT.header"
mv "$OUT.header" "$OUT"
echo "Indel screen for $GENE ($REGION): $RC indel sites -- see $OUT"
