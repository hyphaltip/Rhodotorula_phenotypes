# Last session: 2026-08-17 (session 008, continued) — GWAS post-Tier-A follow-up COMPLETE

## Done
- **LOCO sensitivity COMPLETE**: 120/120 GEMMA per-chromosome-kinship scans (6 traits × 20 chrs, gwas + culled gwasc).
  All anchors reproduce at unchanged p (chroma 2.46e-8, AUC_10 1.44e-8, resilience_30 6.04e-9 all-chr).
  LOCO λ medians 0.34–0.73 track Tier-A λ → the λ<1 deflation is stable near-clonal structure, not per-chromosome rescue.
  → `results/gwas/loco/loco_merged_{gwas,gwasc}.csv`, figures `results/gwas/figures/loco_sensitivity.{png,pdf}`.
- **Tier B set tests**: burden/SKAT/min-p on 178 high-dxy windows (gwas + gwasc). NO set-level signal survives FDR(q<0.05).
  ...burden over-conservative (λ~0.5 sign-cancellation), SKAT mildly deflated (λ~0.6–0.8), min-p recovers Tier-A hits.
  SKAT three-moment approx verified vs exact MC: gwas r=0.984 (n=50, worst drift 0.469), gwasc r=0.997 (n=47, worst 0.871.
  at p~1e-3 — still non-significant). Jobs 27511559 (gwas), 27511892 (gwasc, after L-22 eigdec fix).
  → `results/gwas/tierB/tierb_settests_{gwas,gwasc}.csv`, `tierb_skat_mcver_{gwas,gwasc}.csv`, `tierB_settests.{png,pdf}`.
- **dxy/Fst co-localization**: 12,348 FDR loci → pixy windows. GWAS loci NOT enriched in high-divergence windows
  (dxy OR=0.81 p=0.61; Fst OR=0.41 p=0.065). Phenotypes segregate within near-clonal focal clade (standing variation).
  scaffold_20:100001 (dxy 0.070) flagged as mis-assembly artifact, excluded.
  → `results/gwas/tierB/coloc/*`, `coloc_dxy_fst.{png,pdf}` (scripts `coloc_dxy_fst.py`, `make_coloc_figure.py`).
- **GWAS_REPORT.md §8** written (8.1 LOCO, 8.2 Tier B calibration, 8.3 co-localization, 8.4 outputs).
- **Living repo**: D-13 (Python Tier B + LOCO/coloc finalized), L-22 (`_safe_eigvals` for eigvalsh LinAlgError),
  new finding `gwas-colocalization.md`; FINDINGS_REGISTRY + TODO_REGISTRY + ANALYSIS_MANIFEST + PROGRESS.md updated.

## Next steps
1. Dataviz consult on the three new figure sets (LOCO sensitivity, Tier B calibration, coloc scatter/bar) before publication.
2. Finalize publication-editable methods text for LOCO + Tier B (SKAT three-moment + MC verify) + co-localization.
3. Decision: does Tier-B null + no coloc justify a short "post-Tier-A follow-up" methods paragraph only, or a standalone section?
