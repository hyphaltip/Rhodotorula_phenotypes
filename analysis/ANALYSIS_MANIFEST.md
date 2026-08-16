# Analysis Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

### explore-plate-position
```yaml
name: explore-plate-position
question: How much of strain-replicate variance in color (L*a*b*) and morphology is attributable to plate identity / within-plate grid position, stratified by Copper concentration?
input: db/rhodotorula_phenotypes.duckdb (v_phenotype), endpoint timepoint per plate
scripts:
  - scripts/00_build_dataset.R      # last-image-per-plate endpoint colonies -> endpoint_colonies.rds/csv
  - scripts/01_variance_partition.R # lmer variance partition, STRATIFIED by copper_mm (plate nested in Cu)
  - scripts/02_adjacency_effect.R   # 4-connected grid-neighbor mean predicts own trait (controls copper_mm)
  - scripts/03_plots.R              # per-Cu variance bars, plate-% heatmap, plate heatmap, edge/adjacency plots
  - scripts/04_build_timecourse.R   # per-colony x timepoint table (runs 353-356, 6-hourly) -> colony_timecourse.rds + detection stats
  - scripts/05_time_variance_partition.R # day-block (d1-d5) x Cu variance partition rerun of 01's structure
  - scripts/06_growth_curves.R      # L*/log(area) growth mixed models (quadratic time x Cu, random strain slope/plate/colony)
outputs:
  - results/tables/variance_components_by_copper.csv
  - results/tables/fixed_effects_by_copper.csv
  - results/tables/categorical_position_anova_by_copper.csv
  - results/tables/adjacency_fixed_effects.csv
  - results/tables/colony_timecourse.rds/csv
  - results/tables/variance_components_by_day.csv
  - results/tables/growth_curve_fixed_effects.csv
  - results/tables/detection_curve_by_copper.csv
  - results/figures/variance_components_by_copper.png
  - results/figures/plate_pct_heatmap.png
  - results/figures/day_plate_pct_heatmap.png
  - results/figures/day_variance_components_d1_d5.png
  - results/figures/growth_L_by_cu.png
  - results/figures/growth_area_by_cu.png
  - results/figures/detection_curve_by_cu.png
  - explore_plate_position.html
reproduce: bash analysis/explore_plate_position/run.sh
status: complete
key_findings:
  - Plate-explained variance is SMALL within each Cu concentration (L* ~0.1-4.7%, b* 0-3.2%, area/solidity ~0.4-5.3%) after stratifying by Cu.
  - Earlier pooled (across-Cu) model over-estimated plate variance (23% L*, 33% b*) because plate_id is nested in Cu and the linear copper_mm covariate left non-linear Cu response to be absorbed by the random plate term.
  - Neighbor adjacency signal for L* persists (p~6e-117) with no a*/b*/area/solidity neighbor effect.
  - Plate % declines from day 1 to day 5 e.g. L* at 30 mM 19.1% -> 1.6%: plate position explains least variance at the mature stage captured by the endpoint analysis.
  - Growth curves confirm expectation: L* rises ~44->74 plateauing ~day 3-3.5, log area rises monotonically; Cu slows/limits both. At 25-30 mM ~15-18% of (strain x plate) never reach late timepoints (missingness is the phenotype).
  - run 357 (plates 113-120) excluded from the time course; all plates share the same relative 6-hourly clock except run 353's offset schedule (~3h grid to 117h).
tags: [copper, plate-position, variance-partition, mixed-model, lme4, rhodotorula, timecourse, growth-curve, detection]
```

### growth-rates
```yaml
name: growth-rates
question: How do per-strain colony growth-rate parameters (data-derived peak slopes of log area and of consumer-light intensity) vary with copper concentration and with species, and how do they interact with endpoint color (CIELAB L*, a*, b*) where L* is the light-intensity readout?
input: db/rhodotorula_phenotypes.duckdb (v_phenotype), runs 353-356, per (plate, well) x timepoint
scripts:
  - scripts/00_build_series.R      # per-colony x timepoint table + species join -> colony_growth_series + coverage
  - scripts/01_fit_growth_models.R # Gompertz + Logistic per colony/trait (mu/lambda/A, AIC-preferred), fallback log-linear; primary rate = peak 6 h slope of log(area) [rate_area] and of intensity [rate_int]
  - scripts/02_species_cu_rates.R  # strain x Cu aggregation (up to 4 run-replicates), mixed models rate ~ Cu(factor) + (1|species/strain), Cu x species interaction (well-sampled spp >= 8 strains, unnamed 'sp. clade I' excluded), contrasts vs 0 mM; rate_by_cu_spp / rate_int_by_cu_spp facet top-16 species (4x4)
  - scripts/03_color_interaction.R # endpoint L*/a*/b* as outcomes of rate_area/rate_int x copper (+ species/strain random); documents L*-Intensity ~0.999 collinearity; also endpoint color (L*/a*/b*) vs Cu plots
  - scripts/04_doubling_time.R     # phase-independent estimator: exponential-region specific rate (best-R2 4-pt window) -> doubling time, t50, saturation fraction; tests whether peak-slope Cu trend persists
  - scripts/05_species_cu_sensitivity.R # extent/yield-based Cu-sensitivity index per strain & species (log2 max-area ratio 25-30 vs 0-5 mM, saturation drop, dbl fold), species extent model (Cu x species); sensitivity figures facet top-16 species (4x4)
outputs:
  - results/tables/colony_growth_series.rds/csv
  - results/tables/series_coverage.csv
  - results/tables/growth_model_fits.csv
  - results/tables/growth_model_preferred.csv
  - results/tables/strain_x_cu_rates.csv
  - results/tables/rate_models_species_cu.csv/txt
  - results/tables/rate_area_diff_0mM.csv
  - results/tables/rate_int_diff_0mM.csv
  - results/tables/color_growth_models.csv/txt
  - results/tables/color_growth_anova.csv
  - results/tables/colony_doubling.csv
  - results/tables/dbl_model_cu.csv
  - results/tables/dbl_by_cu.csv
  - results/tables/strain_sensitivity.csv
  - results/tables/species_sensitivity.csv
  - results/tables/species_extent_model.txt
  - results/tables/species_extent_anova.csv
  - results/figures/rate_by_cu_overall.png
  - results/figures/rate_by_cu_spp.png
  - results/figures/rate_int_by_cu_spp.png
  - results/figures/color_L_vs_ratearea_by_cu.png
  - results/figures/color_a_vs_ratearea_by_cu.png
  - results/figures/color_b_vs_ratearea_by_cu.png
  - results/figures/ratearea_vs_rateint.png
  - results/figures/color_L_by_cu.png
  - results/figures/color_a_by_cu.png
  - results/figures/color_b_by_cu.png
  - results/figures/color_L_by_cu_spp.png
  - results/figures/color_a_by_cu_spp.png
  - results/figures/color_b_by_cu_spp.png
  - results/figures/dbl_region_by_cu.png
  - results/figures/saturation_by_cu.png
  - results/figures/sensitivity_species_rank.png
  - results/figures/sensitivity_extent_by_cu_spp.png
  - results/figures/sensitivity_saturation_by_cu_spp.png
  - growth_rates.html
reproduce: bash analysis/growth_rates/run.sh
status: complete
key_findings:
  - Peak growth/brightening rate increases monotonically with Cu (area F=79 p~6e-96; intensity F=540 p~0) - contrary to naive toxicity expectation; consistent with high-Cu colonies staying in log-growth (unsaturating) phase longer, i.e. peak-slope is phase/length sensitive, not a phase-independent rate constant.
  - Gompertz/Logistic lag is structurally unidentifiable in-window (~94% boundary lambda, ~82% area asymptote extrapolated; unbounded GN gives degenerate lags -297..-20 h) -> primary rate estimator is data-derived max 6 h slope.
  - rate_area vs rate_int only r=0.21: area expansion and consumer-light-intensity rise are distinct growth axes.
  - L* (endpoint) vs Intensity_MeanIntensity r=0.999: same axis; L* is "light intensity"; intensity only used as growth readout.
  - Growth-rate x Cu interaction on endpoint L* significant (F=8.16 p~8e-9): faster-growing colonies end lighter, slope largest at low Cu (5mM +3.7 L*/unit rate) smallest at 30mM (+1.3); Cu dominates chromatic outcome (~1000x F on L*/a*).
  - 8,446 colonies (97.6% of 8,652 design slots; 309 strains, 112 plates); 130/2183 strain-Cu aggregates are single-culture, 282 colonies lack species label (-> unknown level).
  - Doubling-time test (04): exponential doubling time is flat ~37 h across 0-30 mM (fold-change 0.99 at 30 mM; Kendall tau between peak-slope rate_area and doubling time = 0.002) -> the rising peak-slope with Cu IS A PHASE ARTIFACT. Cu's real effect is extent-limited: saturation fraction 95.4% -> 73.9%, lower max area, median t50 72.6 -> 64.1 h.
  - Species Cu sensitivity (05, extent-based log2 ratio 25-30/0-5 mM): most tolerant R. glutinis (+1.01) and, among well-sampled, R. mucilaginosa (-0.72); most sensitive R. kratochvilovae (-4.13), R. araucariae (-3.68), R. taiwanensis (-3.06); well-sampled R. diobovata (-2.00). Cu x species interaction on log(max_area) F = 4.34 (p ~ 6.0e-4).
  - Species mislabels in data/metadata/Copper.Strain_info.csv (paludigenum -> paludigena, evergladiensis -> evergladensis) corrected and re-imported via scripts/db/10_import_experiment.py --strain-only (strain 84 & 327 regrouped).
tags: [copper, growth-rate, growth-curve, species, mixed-model, lme4, rhodotorula, color, light-intensity, gompertz, logistic, timecourse]
```

### control-late-timepoint-phenotype
```yaml
name: control-late-timepoint-phenotype
question: What are the per-strain colony size and CIELAB color (L*, a*, b*) summary statistics on control media (Cu=0, YPD) at a late timepoint, aggregated across all replicate colonies?
input: db/rhodotorula_phenotypes.duckdb (v_phenotype + condition_plate_factor Cu=0 + strain)
scripts:
  - scripts/build_phenotype_table.py # latest imaging pass (rounded hour) per strain within a [tmin,tmax] window on Cu=0; per-strain median/mean/var/sd of Shape_Area and per-colony ColorLab_{L*,a*,b*}Median; run per window (70-80, 80-90, 90-110)
outputs:
  - results/phenotype_control_timepoint_70_80.csv
  - results/phenotype_control_timepoint_80_90.csv
  - results/phenotype_control_timepoint_90_110.csv
  - CONTROL_LATE_PHENOTYPE.md
reproduce: bash analysis/control_late_timepoint_phenotype/run.sh
status: complete
key_findings:
  - 314/320 strains have a late Cu=0 image in each window (70-80, 80-90, 90-110 h); 286 have >=3 replicate colonies. Latest pass is the alternate-cadence (75/87/105 h, 84 strains) or main-cadence (78/90/108 h, 230 strains) pass; using an exact max-hour per strain would fragment a single imaging pass (1-2 colonies/strain) and was rejected in favor of rounding to the imaging pass.
  - Strain sets are identical across the three windows; only the sampled pass differs. Colony size highly stable across windows (area_median r ~0.98 for 70-80 vs 90-110).
  - Color drives the table: L* medians 70.2-79.8, a* 0.6-14.2 (carotenoid red/orange), b* -1.4-9.1 on YPD control; colony area 0.8-70 k px.
  - Per-strain stats reproduced exactly by independent manual aggregation (strain 185).
tags: [copper, control, YPD, color, CIELAB, L*a*b*, colony-size, phenotype-table, duckdb, strain, rhodotorula]
```

### 2026-08-15-color-phenotype-space
```yaml
name: 2026-08-15-color-phenotype-space
question: What else should we explore for phenotypic data space or color space in this dataset? (persona-driven ideation campaign, 8 personas x 2 ideas = 16 ideas, all implemented)
input: db/rhodotorula_phenotypes.duckdb (v_phenotype, 211,800 rows) via shared extract data/db_extract.parquet (116 cols)
scripts:
  - scripts/build_series.py        # shared extract: v_phenotype -> data/db_extract.tsv.gz/.parquet + strain_metadata.tsv
  - scripts/common.py              # read_extract/read_meta/save/boot_ci/circular_stats/hue asat helpers
  - scripts/idea_01_arrest.py      # Weibull arrest collapse + Gibrat dispersion (statistical-physicist)
  - scripts/idea_01b_gibrat_within.py # within-strain + size-controlled resolution of the Gibrat dispersion ambiguity
  - scripts/idea_02_color.py       # hue/chroma/L* triple, morphs, onset, run calibration (color-imaging)
  - scripts/idea_03_information.py # species info, redundancy, forward selection (information-theorist)
  - scripts/idea_04_qg.py          # ICC/repeatability, reaction norms, GxE (quantitative-geneticist)
  - scripts/idea_05_repl.py        # rank-PCA/varimax factors + UMAP/HDBSCAN atlas (representation-learning)
  - scripts/idea_06_mediation.py   # Cu -> growth -> pigment endpoint mediation (causal-inference)
  - scripts/idea_07_ecology.py     # environment stratification, species-conditional permutation (trait-ecology)
  - scripts/idea_08_onset.py       # onset threshold sweep + pigment-area coupling (temporal-dynamics)
  - scripts/idea_09_phylogeny.R    # phylogeny join + Mantel perm signal (phylogeneticist, post-campaign) on data/raw/rhodotorula-phyling-protein-tree
outputs:
  - results/idea01_{arrest_collapse,dispersion}.csv; figures/fig01_{master_collapse,arrest_exponent,gibrat_dispersion}.png
  - results/idea02_{species_color_triple,heterogeneity,pigment_morphs,run_calibration,onset_times}.csv; fig02{a,b}_*.png
  - results/idea03_{feature_species_info,redundancy_summary,top_correlated_pairs,forward_selection}.csv; fig03_{mi_top,redundancy_heatmap}.png
  - results/idea04_{icc_audit,heterogeneity_repeatability,reaction_norms,gxe_interaction_mucilaginosa}.csv; fig04_{icc_audit,reaction_norms}.png
  - results/idea05_{pca_variance,varimax_loadings,factor_block_top,atlas,cluster_species_matrix}.csv; fig05_{atlas_umap,cluster_species}.png
  - results/idea06_mediation_decomposition.csv; fig06_growth_pigment_decouple.png
  - results/idea07_{trait_environment,env_trait_means,env_trait_z}.csv; fig07_trait_environment.png
  - results/idea08_{onset_threshold_sweep,pigment_pace,logistic_onset,well_growth_onset}.csv; fig08_{onset_vs_capture,example_growth_pigment}.png
  - results/idea09_{phylo_signal,tip_traits}.csv; results/idea09_pruned_trees.rds; fig09_trait_on_tree_{cu_slope,baseline_chroma,dispersion}.png
  - FINDINGS.md  # master synthesis of all 8 results
reproduce: pixi run python scripts/idea_XX_name.py per script
status: complete
key_findings:
  - sd(log10 Asat) "dispersion widens with Cu" (idea 01 pooled) is largely a SIZE/FLOOR ARTIFACT: sd(log10 Asat) ~ -0.63 correlated with colony size; within-strain raw slopes +0.0033/Mm (p~5e-16) collapse to ~0 size-controlled (p=0.46); only R. taiwanensis genuinely widens (+0.014, p=0.03) and R. paludigena narrows (-0.025, p=0.001) (idea 01b).
  - Cu acts mostly THROUGH growth: log(chroma) loss at 90-110 h is ~fully mediated by colony area (mediation frac ~1.2, direct b~0 per well, boot CI [0.49,2.69]); only the a* (redness) channel keeps a ~40% direct, growth-independent Cu effect.
  - Species/stress/environment signal rides on intra-colony heterogeneity and shape, NOT on L*/intensity: L*Median has the least species MI (~0.07 bit vs best feature 0.19 bit of H=3.0); environment associates only with b*CoeffVar (p_strat~0.006) and L*Median after species-blocking. Camouflage caveat: a*/b*CoeffVar are colony-noise (ICC_strain 0.19/~0) while shape/chroma are repeatable (ICC 0.82-0.93).
  - The strain atlas is a CONTINUUM: varimax factors F1 color-level / F2 heterogeneity / F6-F8 shape are block-diagonal, but UMAP+HDBSCAN finds no discrete species clusters; GxE on log(chroma) within R. mucilaginosa is real but modest (F=1.66, p<<1e-6).
  - Define onset with chroma>7-10: basal chroma noise ~1.5-2 in late Cu=0 colonies; at thr7 median onset 48 h, 91% ever pigment, while >=20 mM Cu gives 0% ever pigment. Colony size gates onset only weakly (rho_onset,size 0.30-0.42 across thresholds).
  - Strain phenotypes carry phylogenetic signal when measured tree-robustly (idea09): Mantel permutation, all-strains scope, l10med_fixed r=0.43, partial_slope_sd_cu r=0.31, intercept_logchroma r=0.19, slope_logchroma_per_mM r=0.10 (all p<=0.003, n~266-272); pace_loglog n.s. Within R. mucilaginosa only baseline chroma (p=0.002) and size (p=0.027) retain signal -> between-species structure, consistent with within-species continuum (idea05).
  - IMPORTANT methodological caveat (idea09): the PHYling protein tree is a near-comet topology (giant polytomy of near-duplicate R. mucilaginosa; 167/541 edges <=1e-7, 22 zero-length pendant tips). Blomberg's K collapses to ~1e-7 here regardless of truth (BM power check K=2.25) and likelihood lambda is numerically unstable (geiger lnL inconsistent) -> use rank-based Mantel permutation for phylogenetic signal on near-comet/genome-cluster trees, never raw K.
  - Species monophyly on this tree (idea09): dairenensis, diobovata, graminis, kratochvilovae, sphaerocarpa, sp. clade I monophyletic; mucilaginosa/paludigena/taiwanensis/toruloides NOT.
tags: [ideas, ideation, persona, phenotype-space, color-space, CIELAB, mediation, information-theory, heritability, GxE, UMAP, HDBSCAN, varimax, ecology, onset, rhodotorula, copper, duckdb, pixi, phylogenetics, mantel, phylogenetic-signal]
```
