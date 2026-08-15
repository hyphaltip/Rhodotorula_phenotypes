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
