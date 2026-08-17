# Last session: 2026-08-17 (session 009, continued) — GWAS Tier D/E/G COMPLETE + figure & methods

## Done
- **Tier D — locus→gene annotation** (`scripts/annotate_gwas_loci.py`, bug-fixed FDR filter /
  BSLMM `PIP<num>` regex / int-chr→scaffold mapping): `results/gwas/tierD/tierD_fdr_snps_annotated.csv`
  (12,348 FDR-sig rows), `tierD_independent_loci.csv` (5,286 loci; 250 kb clump, single lead/chr),
  `tierD_bslmm_loci_annotated.csv` (8 top-PIP loci). Notable genes: chroma scaffold_8 →
  **telomerase RT (OM429_004009)**; AUC_10 scaffold_10 → **DBP3 RNA-dependent ATPase (OM429_004640)**,
  RNA-pol-I TF, endodeoxyribonuclease; BSLMM chroma scaffold_3 → **methionine aminopeptidase 1 (OM429_001379)**.
- **Tier E — fine-mapping** (`scripts/finemap_credible_sets.py`, Wakefield ABF in z-space,
  NCP prior SD=0.2, logsumexp, candidate p<1e-3): `results/gwas/tierE/tierE_credible_sets.csv`
  (42 sets). chroma_lead_10 95% CS n=24 (lead pp=0.054, β=1.11±0.19, common AF 0.81) =
  paper-grade well-bounded anchor. Rare-EF loci (af≈0.015) wide (auc10 DBP3 CS n=67,
  rare_driven) or fully-resolved singleton (resil scaffold_13 95% CS n=1, pp=0.615).
- **Tier G — prior-locus replication**: prior lab's growth-rate chr13:13_30149 (p=1.68e-11)
  **replicates in our AUC_10** via nearest proxy scaffold_13_30134 (15 bp): p_wald=4.03e-6,
  FDR-sig, β=804,778, af=0.015 — inside the single chr13 rare-haplotype block (217 FDR SNPs →
  one lead 13_791853 p=2.4e-7). Gene at locus = **OM429_005439** (hypothetical, unannotated,
  flanked by Ark1-family kinase OM429_005441). Other traits null → replication is
  growth-phenotype-specific.
- **Tier D/E figure** (`scripts/make_tierde_figure.py`, 2-panel, PNG+PDF per convention):
  `results/gwas/figures/tierde_gene_finemap.{png,pdf}`. Panel A = gene map of the scaffold_10
  anchor region (chroma lead 384,905 & AUC_10/DBP3 lead 396,172 ~11 kb apart): -log10 p of
  FDR-sig SNPs per trait + schematic gene track. Panel B = the 14 anchors' 95% credible-set
  size (log x) vs lead posterior probability, coloured by rare-driven flag.
- **Publication methods paragraph** added as GWAS_REPORT.md §10 — Tier A (kinship-only GEMMA
  `-lmm 4 -k`, all-201 + culled-173, BH-FDR q<0.05) → Tier B (burden/SKAT/min-p on pixy
  high-dxy windows, MC-verified) → Tier C (BSLMM 100k, 20% burn-in) → LOCO → Tier D (6,799-gene
  annotation, 250 kb clump) → Tier E (z-space ABF credible sets) → Tier G (nearest-proxy
  replication of chr13:13_30149).
- **Docs**: GWAS_REPORT.md §9 (9.1–9.4) & §10 (methods), PROGRESS.md §14, FINDINGS.md Tier D/E/G
  section, living: D-14, L-24, new finding block in `gwas-colocalization.md`, FINDINGS_REGISTRY,
  ANALYSIS_MANIFEST (`gwas-tierdeg` + `make_tierde_figure.py`/figure outputs), session logs
  `2026-08-17-008` & `2026-08-17-009`, TODO_REGISTRY: closed dataviz note updated, new items
  `rare-ef-imputation.md` (high) + `functional-followup.md` (medium).

## Key gotchas (L-24)
- GEMMA `-lma` has multiple p columns (Wald vs score): chr13:30134 p_wald 4.03e-6 vs 7.90e-6 —
  always cite the column used (p_wald).
- 250 kb clump keeps only the single best SNP/chromosome as lead → second distant block on
  the same small scaffold collapses (the chr13 case). Enumerate per-window top loci from
  `tierD_fdr_snps_annotated.csv` instead.
- BSLMM summary CSV is unquoted-comma malformed → parse with `csv.reader`.
- Wakefield ABF must be in z-space (trait-scale betas span ~1e5 AUC_10 vs ~1 chroma).

## Next steps
1. **Functional follow-up** targets: OM429_005439 (replicated locus, unannotated),
   OM429_004640 DBP3 (rare-EF), OM429_004009 telomerase RT (chroma scaffold_8) — see
   `todo/functional-followup.md`.
2. Rare-EF loci (auc10 DBP3, chr13 block) need imputation/denser markers for resolution —
   see `todo/rare-ef-imputation.md`.
3. Dataviz consult for the new `tierde_gene_finemap.{png,pdf}` (TODO_REGISTRY open item);
   figure not yet visually QA'd by a human.
