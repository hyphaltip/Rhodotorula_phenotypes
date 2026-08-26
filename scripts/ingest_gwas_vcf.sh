#!/usr/bin/env bash
# Symlink a versioned genotype VCF pair into data/raw/genotypes/<version>/
# with a checksum manifest. Usage:
#   bash scripts/ingest_gwas_vcf.sh <version> <source_dir> <basename_prefix>
# Example:
#   bash scripts/ingest_gwas_vcf.sh RmucY2510_v2 \
#     /bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rmuc_popgen_NRRLY2510/vcf \
#     RmucY2510_v2.All
set -euo pipefail
VERSION="$1"
SRC_DIR="$2"
PREFIX="$3"
cd "$(dirname "${BASH_SOURCE[0]}")/.."

DEST="data/raw/genotypes/${VERSION}"
mkdir -p "$DEST"

MANIFEST="$DEST/MANIFEST.yaml"
echo "version: ${VERSION}" > "$MANIFEST"
echo "ingested: $(date -I)" >> "$MANIFEST"
echo "source_dir: ${SRC_DIR}" >> "$MANIFEST"
echo "files:" >> "$MANIFEST"

for suffix in All.SNP.combined_selected.vcf.gz All.SNP.combined_selected.vcf.gz.tbi \
              All.INDEL.combined_selected.vcf.gz All.INDEL.combined_selected.vcf.gz.tbi; do
  base="${PREFIX%.All}.${suffix}"
  src="${SRC_DIR}/${base}"
  if [ ! -e "$src" ]; then
    echo "ERROR: source file missing: $src" >&2
    exit 1
  fi
  real="$(readlink -f "$src")"
  ln -sf "$real" "$DEST/$base"
  size=$(stat -c%s "$real")
  sum=$(sha256sum "$real" | cut -d' ' -f1)
  {
    echo "  - name: ${base}"
    echo "    resolved_path: ${real}"
    echo "    size_bytes: ${size}"
    echo "    sha256: ${sum}"
  } >> "$MANIFEST"
  echo "ingested ${base} (${size} bytes, sha256 ${sum:0:12}...)"
done

echo "Manifest written: $MANIFEST"
