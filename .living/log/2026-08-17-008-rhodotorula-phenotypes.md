# Session Log 2026-08-17-008 — Tier D (locus→gene), Tier E (fine-mapping), Tier G (prior-locus replication) COMPLETE

## What was done
- **Tier D** — `scripts/annotate_gwas_loci.py` (bug-fixed: FDR05_sig filter, BSLMM `PIP<num>`
  regex, int-chr→`scaffold_N` mapping) run → `results/gwas/tierD/`:
  `tierD_fdr_snps_annotated.csv` (12,348 FDR-sig rows), `tierD_independent_loci.csv`
  (5,286 loci; 250 kb clump; AUC_0 2832 / AUC_10 1636 / resilience_30 404 / chroma 282 /
  rest <60), `tierD_bslmm_loci_annotated.csv` (8 BSLMM top-PIP loci, 7 overlap a gene).
  Notable genes: chroma scaffold_8 → telomerase RT (OM429_004009); AUC_10 scaffold_10 →
  DBP3 RNA-dependent ATPase (OM429_004640), RNA-pol-I TF, endodeoxyribonuclease;
  BSLMM chroma scaffold_3 → methionine aminopeptidase 1 (OM429_001379).
- **Tier E** — `scripts/finemap_credible_sets.py` (Wakefield ABF in z-space, NCP prior SD=0.2,
  logsumexp PIPs, candidate filter p<1e-3) run → `results/gwas/tierE/tierE_credible_sets.csv`
  (42 sets × 14 anchors × 90/95/99%). chroma_lead_10 95% CS n=24 (lead pp=0.054,
  β=1.11±0.19, common AF 0.81) well-bounded; auc10 DBP3 and resil scaffold_13 rare-driven
  (resil 95% CS n=1, pp=0.615); ic50 95% CS n=10 (pp=0.183).
- **Tier G** — prior lab's growth-rate top hit chr13:13_30149 (p=1.68e-11, snpID 100002,
  from `SNP_GWAS_T37C_linear.tsv`) **replicates in our AUC_10** via nearest proxy
  scaffold_13_30134 (15 bp away): p_wald=4.03e-6, FDR-sig, β=804,778, af=0.015. Exact pos
  30149 MAF-filtered out of our panel. 30134 sits inside the single chr13 rare-haplotype
  block (217 FDR-sig SNPs → one lead 13_791853 p=2.4e-7, block spans 11,701–800,664).
  Gene at locus: **OM429_005439** (scaffold_13:30,096–30,670), hypothetical protein, no
  GO/InterPro/PFAM; flanked by Ark1-family Ser/Thr kinase (OM429_005441, +6.7 kb).
- Documentation: GWAS_REPORT.md §9 (Tier D/E/G + §9.4 outputs), PROGRESS.md §14,
  FINDINGS.md Tier D/E/G section, living decisions D-14, learning L-24, new finding
  block in `.living/findings/gwas-colocalization.md`, FINDINGS_REGISTRY.md entry.

## Key numbers / corrections
- FDR sig file = 12,348 rows (ALL FDR-sig; earlier assumed partial — corrected across docs).
- GEMMA `-lma` gives multiple p columns (Wald vs score test): chr13:30134 p_wald 4.03e-6 vs
  score-test p 7.90e-6 — both tiny; always cite the column used (p_wald).
- Tier D 250 kb clump keeps only the single best SNP/chromosome as lead → second distant
  block on the same small scaffold collapses (chr13 case). Use
  `tierD_fdr_snps_annotated.csv` for per-window top-loci enumeration.

## Decisions
- D-14: Tier D/E/G executed with z-space ABF fine-mapping, 250 kb single-lead-per-chr clump,
  BSLMM CSV defensive parsing, nearest-proxy replication testing (15 bp offset stated).
- Replication verdict: prior growth-rate locus replicates in the growth-like trait AUC_10
  only (other traits null: best nominal clone_mean_area p=0.044, bright p=0.082).

## Next steps
- Functional follow-up targeting: OM429_005439 (replicated locus, unannotated), OM429_004640
  DBP3 (rare-EF), OM429_004009 telomerase RT (chroma scaffold_8).
- Rare-EF loci (auc10 DBP3, chr13 block) need imputation/denser markers for resolution.
- Consult TODO_REGISTRY.md for the open dataviz item; finalize report/dataviz.

## Tags
gwas, tierd, annotation, tierte, finemapping, credible-sets, abf, tierg, replication,
prior-locus, chr13, AUC_10, DBP3, telomerase, methionine, aminopeptidase, OM429_005439,
near-clone, rare-EF, decision, learning, finding
