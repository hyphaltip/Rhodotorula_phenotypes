---
session_id: 2026-08-15-005
project: rhodotorula-phenotypes
branch: main
started: 2026-08-15T22:30:00-0700
ended: 2026-08-15T23:15:00-0700
duration_minutes: 45
files_changed: ~12
---

## Session Log

### 22:30 — Idea 11: GWAS power & variant feasibility (part of the color-phenotype campaign)
- Established the target set: 278 phenotyped strains -> 200 fully-tip strains; 178 effective independent haplotypes (redundancy 1.12, 22 near-clones collapsed). Min detectable per-SNP R² at n=202: 0.164 GW (5e-8) / 0.140 exome (1e-6) / 0.100 candidate (1e-4); detectable allele ≈0.71 SD (≈22-25% of add. genetic variance); Cu-slope NOT mappable (ICC 0.19) -> GWAS targets size/chroma/pace.

### 22:40 — Located existing SNP panel (key breakthrough)
- Found `/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/` — completed pop-genomics + GWAS project for these strains. `vcf/RmucY2510_v2.All.SNP.combined_selected.vcf.gz` = 728,581 SNP sites x 422 samples, haploid GTs, GATK-hard-filtered (PASS-only). 201/278 strains genotype-covered; 200/201 complete for all 4 GWAS traits. Verified v2/v1 sample sets identical; GEMMA 0.98.3 + plink2 + bcftools in their `GWAS/vc2gwas_env/bin`. Prior-art: 218-strain growth-rate GWAS (T4C-T37C/Salt6) found few hits (chr13 13_30149 p=1.68e-11) - consistent with idea-11 large-effect-only limit.

### 23:05 — Decision & write-up
- D-8: reuse existing SNP panel + GEMMA framework (no de-novo variant calling), with extra QC (biallelic, MAC>=5, missingness<=10%, depth-censor aneuploid scaffolds). L-15: always survey shared Population_Genomics projects before building variant pipelines. Written idea 11 into FINDINGS.md (section 11), 00_index (row 11), ANALYSIS_MANIFEST (new YAML block), `.living/findings/phenotype-color-space.md` (F8/F9/F10), learnings L-15, decisions D-8.

### Next
- GWAS dataset prep: `bcftools view -S` our 201 strains -> filter (biallelic, MAC>=5, geno<=0.1, mind<=0.1, censor high-depth scaffolds) -> 4-trait phenotype TSV -> PLINK LD-pruned kinship -> GEMMA -lmm 4 per trait -> manhattan/QQ figures + significant-hits summary.
