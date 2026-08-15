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
outputs:
  - results/tables/variance_components_by_copper.csv
  - results/tables/fixed_effects_by_copper.csv
  - results/tables/categorical_position_anova_by_copper.csv
  - results/tables/adjacency_fixed_effects.csv
  - results/figures/variance_components_by_copper.png
  - results/figures/plate_pct_heatmap.png
  - explore_plate_position.html
reproduce: bash analysis/explore_plate_position/run.sh
status: complete
key_findings:
  - Plate-explained variance is SMALL within each Cu concentration (L* ~0.1-4.7%, b* 0-3.2%, area/solidity ~0.4-5.3%) after stratifying by Cu.
  - Earlier pooled (across-Cu) model over-estimated plate variance (23% L*, 33% b*) because plate_id is nested in Cu and the linear copper_mm covariate left non-linear Cu response to be absorbed by the random plate term.
  - Neighbor adjacency signal for L* persists (p~6e-117) with no a*/b*/area/solidity neighbor effect.
tags: [copper, plate-position, variance-partition, mixed-model, lme4, rhodotorula]
```
