# Session 2026-08-16-007 — Rhodotorula phenotypes

**Date**: 2026-08-16 (afternoon)
**Project**: rhodotorula-phenotypes
**Branch**: main
**Focus**: Next-gen expert-guided GWAS — Tier A (24 LMM scans) + Tier C (BSLMM architecture); full-genome pixy job running

## Summary
Completed the next-gen GWAS designed per expert QG review (08b): engineered 12 clone-mean
phenotypes (color chroma/sat/bright, copper AUC block, resilience, dose-slope, IC50, size),
built a near-clone-culled 173-strain set (IBS0<0.005 greedy, 28 dropped), ran all 24 GEMMA
LMM scans (12 traits × all-201 + culled-173), summarized (lambda/PVE/top hits + multi-trait
Fisher for color & copper blocks), then ran Tier-C BSLMM on 5 representative traits.
Completed full-genome pixy (job 27502734, COMPLETED 6h35m). Generated 20 manhattan+QQ
figures (png+pdf) and embedded them in the report. Finished Tier-C significance with FDR
after rejecting the maxT permutation null (uncalibrated on near-clonal panel, L-21/D-11).

## Key outputs
- `09-NEXT-GWAS-DESIGN.md` — design + confirmed decisions (D10 in `.living/decisions.md`).
- `results/gwas/tierA_summary/` (29 files) — `tiera_summary.csv`, per-trait assoc, multi-trait Fisher top500.
- `results/gwas/tierC_summary/` (11 files) — `tierc_bslmm_summary.csv` + hyp/gamma MCMC traces.
- `results/gwas/fdr/` (13 files) — `tiera_fdr_summary.csv` + per-trait FDR-significant SNP lists.
- `results/gwas/figures/` (20 png + 20 pdf) — manhattan+QQ, embedded in `GWAS_REPORT.md`.
- `GWAS_REPORT.md` Stages 6 & 7; `PROGRESS.md` §6-9.
- `results/gwas/pixy/` — full-genome `genome_{pi,dxy,fst,watterson_theta}.txt` (job COMPLETED).
- `culled_keep.list`, `gwasc.*`, `kinsc.cXX.txt` (culled kinship) in scratch workdir.

## Tier-A results (all-201, kinship-only LMM)
- Color multi-Fisher p=7.5e-16 @ scaffold_10:384905 (chroma 2.4e-8, sat, bright co-localize).
- Copper block p=5.9e-12 @ scaffold_16:421208; AUC_30 scaffold_4:508257 repeats both sets.
- Lambdas mostly <1 (conservative); clone_mean_area λ=2.12, AUC_0 λ=2.74 inflated; IC50 culled λ=0.19 (structure inflation).
- Culled-173 top color hit shifts to scaffold_3:570085 → rank instability flagged for LOCO.

## Tier-C BSLMM results (100k MCMC, 20% burn-in, all-201)
- **AUC_10 near-oligogenic**: PVE 0.40, PGE 0.96, only ~3 active LD-clusters; top loci
  scaffold_5:533603 + scaffold_2:939277 (PIP 1.0). Strongest interpretable copper signal.
- chroma: PVE 0.24, PGE 0.39, ~138 clusters, scaffold_3-enriched loci (26959/137161/132903/7440 PIP 1.0).
- resilience_30: PVE 0.25, PGE 0.62, sparse; scaffold_16:563722 (PIP 0.95) ≈ Tier-A multi-copper block.
- clone_mean_area: PVE 0.07 (weakly heritable — matches inflated λ).
- BSLMM converges with Tier A (scaffold_16 copper, scaffold_3 color).

## Learnings
- **GEMMA 0.98.3 BSLMM syntax**: `-bslmm 1` (option 1 = BSLMM; 2 = ridge/GBLUP; 3 = probit;
  4/5 = DAP). 100k samples default fine. Output: `output/bslmm_<o>.hyp.txt`
  (cols h pve rho pge pi n_gamma), `bslmm_<o>.gamma.txt` (per-sample × 300-LD-cluster
  matrix, value = representative SNP index, 0 = inactive). **PIP must be computed from
  gamma-nonzero fraction, NOT by averaging the raw index values.**
- **Meff proviso**: PCA-based Meff on n=201 is bounded by sample size (~197), NOT the
  independent-test count — needs SNP×SNP LD (intractable at n≪M). Use FDR or np permutation
  null; nominal 5e-8 conservative (~5k–50k effective tests).
- **maxT phenotype-permutation is unusable on near-clonal panels (L-21)**: fixed kinship +
  permuted phenotype → min-p collapses to 1e-17..1e-18 under the null (expected ~2e-6).
  Median stays uniform but the tail (which a maxT threshold uses) is garbage. Use FDR instead.
  FDR doubles as an inflation detector (AUC_0 λ=2.74→7340 "sig", clone_mean_area λ=2.12→0).
- Near-clone-relatedness panels: BSLMM sparse clusters collapse to ~300 LD clusters,
  giving clean near-oligogenic architecture calls (AUC_10: 3 loci carry ~all h²).

## Decisions
- Tier A → Tier B/C first; LOCO as committed sensitivity follow-up (D-10).
- Tier C BSLMM on 5 representative traits (color flag + copper reps + size) — all-201 set.
- **Tier-C significance: BH FDR, not maxT permutation (D-11)** — permutation null uncalibrated
  on the near-clonal panel; FDR robust and endorsed by expert rec #4.

## Next steps
- **Tier B**: SKAT/set tests on pixy high-dxy windows (full-genome pixy tables now ready).
- LOCO sensitivity run (runbook `10-LOCO-RUNBOOK.md` ready, wrapper validated).
- dxy/fst trait co-localization against Tier-A/FDR loci.
- Manhattan/QQ already done; consider FDR-flagged manhattan overlays.

## Tags
gwas, gemma, LMM, bslmm, PVE, PGE, near-oligogenic, copper, color, chroma, AUC, kinship,
kinship-only, multi-trait, fisher, culled, near-clone, pixy, idea11
