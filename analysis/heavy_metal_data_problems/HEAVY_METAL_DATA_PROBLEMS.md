# Heavy-metal data: known problems (2026-10-08)

Scope: the five-metal intermediate data (`data/raw/heavy-metal-array-intermediate/`, table `heavy_metal_measurement`, 517,371 colony rows) and the analyses built on it. Every number below was computed in this session, unless marked "subagent-reported" (not re-run by the parent) or "not checked".

## A. Data design problems

1. **Zinc layout.** 9 plates in 2 runs, treated as one experiment. Two disjoint strain sets: set A (run d000388) has 80 strains at doses 0, 5, 10, 15; set B (run d000390) has about 70 strains at doses 10-30. One plate per dose and 1 well per strain per dose, so no replicates. Set B has no dose-0 baseline.
   - **Rescue (2026-10-08).** The late window used for the other metals (end 107.6 h) dropped plate d000388/41 (dose 15, stops at 83.7 h) and cut plate d000390/100 (dose 30, stops at 89.7 h). A window ending at 80 h keeps both. Plate d000390/96 (dose 10; 7 images, 35.8 h) cannot be rescued. Result: 8 plates, 588 wells, 152 strains (502 wells after the 2,000 px cutoff).
   - **What is valid.** Strain comparisons within a plate (ranks, between-strain differences, correlations across strains) are not confounded with the plate effect. Strain a* ranks agree between dose plates (Spearman 0.40-0.84). The two sets agree at the shared dose 15 (a* 17.7 vs 18.1, b* 9.8 vs 10.2).
   - **What is not valid.** The population-mean dose effect cannot be separated from the plate effect, and mixed-model standard errors and p-values for Zinc dose effects are understated (the plate variance collapses to zero). See `analysis/carotenoid_stress_vs_baseline/report/REPORT.md` section 13.
2. **Unnamed rows.** 10,921 rows have NULL `strain_id` (Cr 5,923; Cu 3,489; Pb 1,509; none in Fe or Zn), all in the last two runs of those metals. 4,399 rows are `Control-N` (13 IDs, no sample_name or species). Both are excluded from the analyses.
3. **Uneven imaging spans.** Plate span (h): Cr 114-119; Cu 96-117 (one 96 h plate); Fe 36-96; Pb 101-113; Zn 36-114. A plate-relative window would compare different growth stages, so an absolute window [T-24, T] per metal is used (T = median span; Zinc T = 80 h).
4. **Missing wells and survivor selection.** Wells with no growth have no objects and are absent. Wells dropped for fewer than 2 window images are listed in `wells_lost_by_dose.csv` (largest at the top doses, for example Cr 1.2: 158 of 1,210). High-dose estimates describe wells that grew.
5. **Iron coverage.** 146 of 219 strains have both dose 0 and dose 30. Run d000408 has a single short plate at dose 25 (83 wells lost at dose 25).
6. **Several objects per well.** Up to 28 objects in one well-image. The largest object is kept; a sensitivity table using all objects (`wells_allobj.csv`) exists but is not modelled.
7. **Dose units.** Not stated in the source. Cr spans 0-1.2; the other metals 0-30. Do not pool Cr.
8. **Different column sets.** Pb and Zn have 32 extra columns (later pipeline version). Cross-metal comparison of those columns is not possible.

## B. Strain identity problems

9. **Species placeholder.** 16 strains still read `Species Not Found` after 3 user-confirmed overrides. Triage: 5 can take species from the old table (not confirmed), strain 65 (`223D-8` vs `223A-8`) is unresolved, and 10 have no record in the old table, the VCF or the tree. See `data/metadata/heavy-metal-array-intermediate/strains_needing_species_id.tsv`.
10. **Duplicate sample name.** `strain_id` 165 and 269 share `TFCN_17-332Y-1`.
11. **Old vs new names.** 16 of 320 strain_ids differ between the old strain table and the new `sample_name`.
12. **Old vs new assignment (subagent-reported, not re-run).** Re-derived Copper traits agree poorly with the old ones (Pearson r 0.49 log10 area, 0.37 AUC_0, 0.10 AUC_ratio_10, 0.32-0.53 texture) although image-level median area agrees (r=0.9997) and 2,397 image names match. Object-level strain assignment agreed only 41%. Unresolved. Do not regenerate GWAS traits from DuckDB until this is understood.
13. **Reconcile on the new table (subagent-reported).** 213 exact, 74 fuzzy, 33 unmatched, 54 collisions against the VCF.

## C. Analysis limitations (a* study, `analysis/carotenoid_stress_vs_baseline/`)

14. **a* and size.** a* correlates with ln(area) at about 0.7 in every metal. Size is partly an effect of stress, so the model at fixed size is a direct effect, not the whole induction. One size slope is used for all strains.
15. **Small colonies.** a* depends on colony area even in unstressed colonies (flat below about 1,500-2,000 px, then rising). The main analyses drop colonies below 2,000 px; this was tested for Cr and Cu and assumed for Pb, Fe and Zn (Zinc has no unstressed colonies below 3,000 px). After the cutoff few Pb colonies remain at doses 20-30. See REPORT.md section 10.
16. **Convergence artifact.** In Cr, Cu and Pb strains converge to a common a* at high dose (strain share of variance falls to 0.24, 0.28, 0.03), which forces a strongly negative baseline-vs-induction correlation (-0.8 to -0.99). Do not read it as biology. Pb M1 is near the boundary.
17. **Non-linear dose response.** Cr and Pb are biphasic. The linear-dose models hide this; use the dose-factor models.
18. **Species tests.** R. mucilaginosa has 215 strains; the other species have 5-18. Kruskal-Wallis on strain means ignores phylogenetic relatedness. A mixed model with species fixed and strain random has not been run.
19. **Trait correlations.** Many traits are near-duplicates (radius, Feret diameter, minor axis; `ColorLab_a*Medoid` is a* itself). Iron shows opposite-sign size correlations.
20. **Not inspected.** Fig 2, 4, 5 and 6 of the a* study; 5 of 7 overview figures (subagent).

## D. Pipeline and bookkeeping problems

21. **Scripts reading removed inputs.** `analysis/gwas/scripts/common.py`, `build_gwas_phenotypes.py`, `check_mas_gates.py`, `analysis/ideas/...`, `analysis/control_late_timepoint_phenotype/...` will fail if run. The `cu_doseauc_v0151` source dataset is gone.
22. **Stale docs.** README.md, DATABASE_DESIGN.md, SCHEMA.md still describe `colony_measurement`.
23. **Gitignored data not covered by the tag.** The old DB and Parquets are not in `data-v1-pre-metal-replace`. A DB backup is in `db/`.
24. **DB redundancy.** The 6 strain text columns hold 41 MB raw in `heavy_metal_measurement`; the on-disk saving from dropping them is not measured.
25. **Bugs found and fixed this session.** A datetime unit error (window included all rows) and a plate-relative window flaw. Earlier numbers from those versions were discarded.

## Update 2026-10-08: Fe and Zn tables are incomplete (D-56)
- Run-by-run comparison with the 0.15.1 tables (`old/0.15.1_Analysis/Results`): ingested Zinc has 2 of 5 runs (9 of 119 plates); ingested Iron has 3 of 5 runs (62 of 120 plates, run 408 only 6 of 28). Cr, Cu and Pb match the 0.15.1 row counts.
- The data owner's note (`Data/Temp/zinc_sources/ZINC_SOURCES.md`, 2026-09-30) says the `Data/Interm` Zinc table comes from a run where 2,121 of 2,275 images failed.
- Problem 1 (Zinc layout) above is explained by this truncation. Fe and Zn are held out of the analyses until complete tables exist.
