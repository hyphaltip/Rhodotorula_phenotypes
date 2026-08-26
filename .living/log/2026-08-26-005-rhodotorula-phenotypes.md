---
session_id: 2026-08-26-005
project: rhodotorula-phenotypes
branch: "main"
started: 2026-08-26T15:20:00-0700
ended: 2026-08-26T16:10:00-0700
duration_minutes: 50
files_changed: 30
---

## Session Log

### 15:20 — Session started; picked up candidate_gene_alignment extension + color-phenotype association
- Resumed from 2026-08-26-004 (candidate-gene alignment pilot). Focus: user asked whether any candidate gene variation can be associated with color phenotypes.

### Correction found: pilot's "0 indels" was a silent toolchain failure
- `screen_indels.sh` reported 0 indels for both pilot genes, but with bcftools on PATH the real INDEL VCF has 433 (OM429_003333) and 821 (OM429_003336) records in those CDS±2kb windows. Root cause: `bcftools ... 2>/dev/null | grep -vc "^#" || true` reported an empty-pipe count of 0 when the tool was missing. Rewrote script to fail loudly; added `check_gene_coding_indels.py` (CDS-exon-exact screen) — **10/10 target genes have ≥2 segregating in-CDS indels** (`results/gene_coding_indel_screen.csv`). Learned L-31/L-32.

### Extended pipeline to all 8 GWAS-locus genes (`run_extend_to_gwas_genes.sh`)
- Fixed batch driver: absolute BCFTOOLS/SAMTOOLS/JAVA paths, htslib/bin on PATH, LD_LIBRARY_PATH for htslib/libdeflate/java, and explicit snpEff `-dataDir <.../snpEff/data>` (DB config's `data.dir=./data/` is relative). All 8 genes produced dna/protein FASTA, indel screen, snpEff region, variant table, alignment PNG.

### Phenotype-association step (population-aware battery) — the core deliverable
- `associate_variants_with_phenotype.py`: 954 tests across 10 genes × 6 color traits → **39 FDR-sig / 40 likely_real** SNP coding variants.
- `associate_indels_with_phenotype.py`: 2046 tests (1110 frameshift) → **32 FDR-sig** in-CDS indel associations.
- **Headline**: OM429_000065 (sat-locus gene) c.839C>T p.Ala280Val (lead SNP 208569) AND a +6bp in-frame GAGCGG-repeat insertion at 208398 in 100% LD with it — both associate with lab_a/chroma/sat (meta_p≈8.8e-11, partial R² 0.086–0.10, replicated 3/3 testable pops). The insertion is invisible to SNP-only pipelines (learned L-33).
- OM429_003333 splice_donor is monomorphic (0/213 alt) — resolves that TODO as N.A.
- Carotenoid genes show no replicated color association (consistent with no co-localization).

### Post-action protocol updates
- Updated `CANDIDATE_GENE_ALIGNMENT.md`, `ANALYSIS_MANIFEST.md`, `results/PROVENANCE.json`. Logged D-22 (decisions), L-31/L-32/L-33 (learnings), finding `candidate-gene-color-association`, TODO registry (3 done, 2 new open items), session log.

## Files Modified
- analysis/ANALYSIS_MANIFEST.md
- analysis/candidate_gene_alignment/CANDIDATE_GENE_ALIGNMENT.md
- analysis/candidate_gene_alignment/results/PROVENANCE.json
- analysis/candidate_gene_alignment/results/{gene_coding_indel_screen,candidate_gene_phenotype_assoc_all10,candidate_gene_phenotype_assoc_all10_within_pop_details,candidate_indel_phenotype_assoc_all10,candidate_indel_phenotype_assoc_all10_within_pop_details,candidate_indel_phenotype_assoc_pilot,candidate_indel_phenotype_assoc_pilot_within_pop_details,candidate_gene_phenotype_assoc,candidate_gene_phenotype_assoc_within_pop_details}.csv + extend_batch.log
- analysis/candidate_gene_alignment/results/OM429_{001533,001415,001430,003729,002663,005034,001521,000065}/** (extension outputs)
- analysis/candidate_gene_alignment/scripts/{screen_indels.sh,check_gene_coding_indels.py,run_extend_to_gwas_genes.sh,associate_variants_with_phenotype.py,associate_indels_with_phenotype.py}
- .living/decisions.md, .living/learnings.md, .living/findings/{FINDINGS_REGISTRY.md,gwas-colocalization.md}
- todo/TODO_REGISTRY.md
