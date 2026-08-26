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

echo "gene,scaffold,region,n_indels" > "$OUT.header"
N=$(bcftools view -r "$REGION" "$INDEL_VCF" 2>/dev/null | grep -vc "^#" || true)
echo "${GENE},${SCAFFOLD},${REGION},${N}" >> "$OUT.header"
mv "$OUT.header" "$OUT"
echo "Indel screen for $GENE ($REGION): $N indel sites -- see $OUT"
