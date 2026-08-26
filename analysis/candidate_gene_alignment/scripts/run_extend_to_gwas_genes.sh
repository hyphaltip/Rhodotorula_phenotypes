#!/usr/bin/env bash
# Extend the candidate-gene alignment pipeline to the 8 GWAS-locus nearest genes
# (pilot = OM429_003333/003336 carotenoid genes, done). Per gene:
#   1. extract_gene_sequences.py  -- per-strain spliced CDS + protein (asserts
#      REF-vs-genome, equal CDS length; renders N for missing calls)
#   2. screen_indels.sh           -- read-only check of the INDEL VCF for the gene's
#      CDS +/- flank (OPUS review: SNP-only premise only holds for the SNP file)
#   3. snpEff (4.3m + RmucNRRLY2510 DB) on the CDS+/-flank region
#   4. build_variant_table.py     -- consequence x genotype x population x phenotype x
#      lead-SNP genotype + each site's Tier-A p-value
#   5. render_alignment_image.py  -- static PNG, strains sorted by lead-SNP genotype.
#
# Must be run from the repo root with samtools/bcftools/snpEff loaded:
#   module load samtools bcftools snpEff/4.3m
#   bash analysis/candidate_gene_alignment/scripts/run_extend_to_gwas_genes.sh
set -euo pipefail

# Resolve module-provided binaries explicitly: they are NOT on the pixi-env PATH
# (pixi prepends its own bin/ and the module dirs drop out), so hardcode the resolved
# absolute paths. Sub-shells (extract_gene_sequences.py subprocess calls) use these.
# NOTE: `module load snpEff/4.3m` resets PATH and silently drops java/21.0.7/bin, so
# JAVA is hardcoded too (discovered 2026-08-26 -- java: command not found at snpEff step).
# Also: bcftools/samtools are dynamically linked against htslib/libdeflate which live
# OUTSIDE the rpath -- run without the module's LD_LIBRARY_PATH and they exit 127.
BCFTOOLS="${BCFTOOLS:-/opt/linux/rocky/8.x/x86_64/pkgs/bcftools/1.22/bin/bcftools}"
SAMTOOLS="${SAMTOOLS:-/opt/linux/rocky/8.x/x86_64/pkgs/samtools/1.22.1/bin/samtools}"
JAVA="${JAVA:-/opt/linux/rocky/8.x/x86_64/pkgs/java/21.0.7/bin/java}"
SNPEFFJAR="${SNPEFFJAR:-/opt/linux/rocky/8.x/x86_64/pkgs/snpEff/4.3m/snpEff/snpEff.jar}"
export PATH="$(dirname "$BCFTOOLS"):$(dirname "$SAMTOOLS"):/opt/linux/rocky/8.x/x86_64/pkgs/htslib/1.22.1/bin:$(dirname "$JAVA"):$PATH"
export LD_LIBRARY_PATH="/opt/linux/rocky/8.x/x86_64/pkgs/htslib/1.22.1/lib:/opt/linux/rocky/8.x/x86_64/pkgs/libdeflate/1.17/lib64:/opt/linux/rocky/8.x/x86_64/pkgs/java/21.0.7/lib:/opt/linux/rocky/8.x/x86_64/pkgs/java/21.0.7/lib/server:${LD_LIBRARY_PATH:-}"

SRC=analysis/candidate_gene_alignment/scripts
OUT=analysis/candidate_gene_alignment/results
GFF=/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/genome/Rhodotorula_mucilaginosa_NRRL_Y-2510.gff3.gz
GENOME=/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/genome/Rhodotorula_mucilaginosa_NRRL_Y-2510.scaffolds.fa
VCF=data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.SNP.combined_selected.vcf.gz
STRAINS=analysis/gwas/results/strain_reconciliation/accepted_vcf_ids.txt
POP=analysis/gwas/data/prior_run_state/pop_assignment_at_run.csv
PHENO=analysis/gwas/results/gwas_next_phenotypes.csv
SNPEFF_CONFIG=/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/snpEff/snpEff.config
SNPEFF_DIR=/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/snpEff
FLANK=2000
# gene    scaffold  cds_start cds_end trait      panel  lead_snp
GENES=(
"OM429_001533 scaffold_3 544610 546394 cu_dose_slope gwas scaffold_3:546065"
"OM429_001415 scaffold_3 227541 229892 lab_L         gwas scaffold_3:229210"
"OM429_001430 scaffold_3 264409 266988 lab_L         gwasc scaffold_3:265939"
"OM429_003729 scaffold_8 38699  41014  lab_a         gwas scaffold_8:38068"
"OM429_002663 scaffold_5 882726 884304 lab_a         gwasc scaffold_5:882396"
"OM429_005034 scaffold_11 608768 612742 lab_b         gwas scaffold_11:608491"
"OM429_001521 scaffold_3 502667 504489 lab_b         gwasc scaffold_3:503756"
"OM429_000065 scaffold_1 207031 212775 sat           gwas scaffold_1:208569"
)

mkdir -p "$OUT"

for SPEC in "${GENES[@]}"; do
  read -r GENE SCAFFOLD START END TRAIT PANEL LEAD <<< "$SPEC"
  GDIR="$OUT/$GENE"
  echo ""
  echo "############################################################"
  echo "# $GENE  ($SCAFFOLD:$START-$END, $TRAIT/$PANEL, lead=$LEAD)"
  echo "############################################################"

  # ---- 1. extract per-strain CDS/protein ----
  echo "--- step 1/5: extract gene sequences ---"
  pixi run python "$SRC/extract_gene_sequences.py" \
    --gene-ids "$GENE" --gff3 "$GFF" --genome-fasta "$GENOME" \
    --vcf "$VCF" --strain-list "$STRAINS" --samtools "$SAMTOOLS" --bcftools "$BCFTOOLS" \
    --out-dir "$OUT"
  echo "--- step 2/5: indel screen ---"
  bash "$SRC/screen_indels.sh" "$GENE" "$SCAFFOLD" "$START" "$END" "$FLANK" "$GDIR/indel_screen.csv"

  # ---- 3. snpEff region annotation ----
  echo "--- step 3/5: snpEff annotation ---"
  REGION="$SCAFFOLD:$((START - FLANK))-$((END + FLANK))"
  mkdir -p "$GDIR/snpeff"
  bcftools view -r "$REGION" -Ou "$VCF" \
    | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o "$GDIR/snpeff/region.vcf.gz" -
  tabix -p vcf "$GDIR/snpeff/region.vcf.gz"
  "$JAVA" -Xmx8g -jar "$SNPEFFJAR" eff -c "$SNPEFF_CONFIG" -dataDir "$SNPEFF_DIR/data" \
    -noLog RmucNRRLY2510 "$GDIR/snpeff/region.vcf.gz" > "$GDIR/snpeff/region.ann.vcf" 2> "$GDIR/snpeff/snpeff.log"
  echo "   snpEff log tail:"; tail -3 "$GDIR/snpeff/snpeff.log"
  # ---- 4. variant table ----
  echo "--- step 4/5: build variant table ---"
  ASSOC="analysis/gwas/results/gwas/tierA_summary/gemma_output/${PANEL}_${TRAIT}.assoc.txt"
  pixi run python "$SRC/build_variant_table.py" \
    --gene-id "$GENE" --ann-vcf "$GDIR/snpeff/region.ann.vcf" \
    --strain-list "$STRAINS" --pop-csv "$POP" --pheno-csv "$PHENO" \
    --lead-snp "$LEAD" --lead-snp-trait "$TRAIT" --assoc-file "$ASSOC" \
    --out-csv "$GDIR/variant_table.csv"

  # ---- 5. alignment image ----
  echo "--- step 5/5: render alignment image ---"
  pixi run python "$SRC/render_alignment_image.py" \
    --dna-polymorphic-csv "$GDIR/dna_polymorphic_positions.csv" \
    --strain-context-csv "$GDIR/variant_table_strain_context.csv" \
    --sort-by lead_snp --phenotype-col "$TRAIT" --gene-id "$GENE" \
    --out-png "$GDIR/alignment.png"

  echo "completed $GENE"
done

echo ""
echo "All 8 GWAS-locus genes extended."
