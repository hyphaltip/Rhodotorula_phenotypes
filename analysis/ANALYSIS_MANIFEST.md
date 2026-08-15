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
  - scripts/02_species_cu_rates.R  # strain x Cu aggregation (up to 4 run-replicates), mixed models rate ~ Cu(factor) + (1|species/strain), Cu x species interaction (well-sampled spp), contrasts vs 0 mM
  - scripts/03_color_interaction.R # endpoint L*/a*/b* as outcomes of rate_area/rate_int x copper (+ species/strain random); documents L*-Intensity ~0.999 collinearity; also endpoint color (L*/a*/b*) vs Cu plots
  - scripts/04_doubling_time.R     # phase-independent estimator: exponential-region specific rate (best-R2 4-pt window) -> doubling time, t50, saturation fraction; tests whether peak-slope Cu trend persists
  - scripts/05_species_cu_sensitivity.R # extent/yield-based Cu-sensitivity index per strain & species (log2 max-area ratio 25-30 vs 0-5 mM, saturation drop, dbl fold), species extent model (Cu x species)
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
tags: [copper, growth-rate, growth-curve, species, mixed-model, lme4, rhodotorula, color, light-intensity, gompertz, logistic, timecourse]
```
