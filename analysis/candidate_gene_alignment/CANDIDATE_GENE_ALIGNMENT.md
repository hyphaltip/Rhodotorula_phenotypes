# Candidate-gene DNA/protein sequence alignment across the 213-strain panel

**Design doc**: `docs/superpowers/specs/2026-08-26-candidate-gene-alignment-design.md`
(reviewed by an Opus general-review pass and a bioinformatics-expert/fable consult
before implementation; both reviews' findings are addressed below).

**Status**: pilot (2/10 genes) complete and validated. Remaining 8 GWAS-locus genes
not yet run.

## What this is (and isn't)

For each candidate gene, this builds a per-strain spliced CDS + protein sequence
across all 213 panel strains, a positional variant table (snpEff consequence/codon/AA
change per site, joined with each strain's genotype, population, phenotype value, and
— for GWAS-locus genes — the locus's own lead-SNP genotype and each coding variant's
own Tier A marginal p-value), and a static color-block alignment image.

**This is descriptive, not a second association test.** A coding variant's own Tier A
p-value is reported for context (to distinguish a plausible driver variant from an
LD-linked passenger within one gene), but no new hypothesis test is performed here —
statistical validation of each locus was already done in
`analysis/gwas/scripts/check_population_vs_locus.py` and `check_mas_gates.py`.

## Pipeline

1. `scripts/extract_gene_sequences.py` — GFF3 CDS exons + genome FASTA + SNP VCF →
   per-strain spliced CDS FASTA, translated protein FASTA, and a polymorphic-positions
   CSV (invariant CDS columns dropped). No MAFFT: sequences are asserted equal-length
   (see below) and built as a plain positional character matrix, not a general aligner
   output — safe given (1) all target genes are single-exon-transcript loci with no
   isoform ambiguity, (2) confirmed no indels in the CDS+/-2kb region (step 2), so
   every strain's CDS really is a pure per-site substitution of the reference.
2. `scripts/screen_indels.sh` — read-only screen of the **separate INDEL VCF**
   (`RmucY2510_v2.All.INDEL.combined_selected.vcf.gz`) against each gene's CDS+/-2kb.
   This VCF is never used to build any sequence; it exists only to confirm the
   SNP-only assumption in step 1 actually holds for that gene (Opus's review flagged
   this as the single most likely blind spot — a coding indel would be invisible to a
   substitution-only pipeline).
3. snpEff (`snpEff/4.3m` module + the pre-built `RmucNRRLY2510` database, format-version
   matched — the default `snpEff` module loads 5.1 and refuses this database) annotates
   the SNP VCF restricted to each gene's CDS+flank region.
4. `scripts/build_variant_table.py` — parses the snpEff `ANN` field for the entry
   matching the target gene, joins per-strain genotype, population
   (`pop_assignment_at_run.csv`), phenotype snapshot, and (for GWAS-locus genes) the
   lead SNP's own genotype and this variant's Tier A p-value.
5. `scripts/render_alignment_image.py` — static PNG, one per gene, strains sorted by
   phenotype (pathway genes, no locus) or by lead-SNP genotype then phenotype
   (GWAS-locus genes), CDS-position ruler on the x-axis, A/C/G/T/N color-coded.

**Follow-up flagged, not built this round**: an interactive viewer (user's explicit
request when scoping output format — static images + CSV first, interactive version
next).

## Pilot genes (carotenoid pathway, no GWAS locus)

| gene | product | CDS length | variant sites (region) | polymorphic (panel) | premature stops | indels in CDS+/-2kb |
|---|---|---|---|---|---|---|
| `OM429_003333` | phytoene synthase / lycopene beta-cyclase | 1293 bp | 58 | 49 | 0/213 | 0 |
| `OM429_003336` | phytoene desaturase (misannotated as "transcription factor" in GFF3 — see D-20) | 1929 bp | 81 | 72 | 0/213 | 0 |

Consequence breakdown (snpEff, CDS+flank region, both genes combined): mostly
upstream/downstream-gene variants from neighboring genes in the same dense scaffold_7
cluster; within each gene's own CDS, synonymous variants dominate, with a handful of
missense calls (7 in `OM429_003333`, 23 in `OM429_003336`) and one
`splice_donor_variant&intron_variant` (HIGH impact) in `OM429_003333`.

### Reference-strain sanity check (resolves design-doc open question 4)

`NRRL_Y-2510` — the reference-genome individual itself — is in the accepted 213-strain
panel. Queried directly: **0/0 (homozygous reference) at all 139 CDS variant sites
across both pilot genes, 0 alt calls.** This is exactly the expected result and
confirms the substitution pipeline has no systematic REF/ALT-swap or off-by-one bug —
if there were one, this strain would show spurious heterozygous/alt calls at every
site since it's genetically identical to the reference assembly.

## Open items before extending to the remaining 8 GWAS-locus genes

- Confirm none of the 8 GWAS-locus genes is mitochondrially encoded (bioinformatics
  reviewer's flag; quick check against the GFF3/genome scaffold naming, not yet done).
- Add an explicit hard assertion of 213/213 strain coverage in
  `extract_gene_sequences.py` (currently implicit via the strain-list filter in
  `get_variants_in_region`) — Opus's recommendation.
- Pull each GWAS-locus gene's lead SNP region (may be outside the CDS — intergenic/
  upstream) so `build_variant_table.py`'s `--lead-snp` join actually has a genotype to
  report; the pilot genes have no lead SNP by design.
- Extend the flanking-window screen and snpEff region to each gene's own coordinates
  (currently these scripts were run ad hoc per-gene; batch this into `run.sh`).

## Files

- `scripts/extract_gene_sequences.py`, `scripts/screen_indels.sh`,
  `scripts/build_variant_table.py`, `scripts/render_alignment_image.py`
- `results/<gene_id>/{dna.fasta, protein.fasta, dna_polymorphic_positions.csv,
  indel_screen.csv, variant_table.csv, variant_table_strain_context.csv, alignment.png}`
- `results/snpeff_pilot/` — region-restricted VCF + snpEff-annotated VCF for the pilot
- `results/PROVENANCE.json` — VCF/genome/database checksums and versions, sanity-check
  result, indel-screen result
