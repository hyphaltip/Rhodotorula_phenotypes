#!/usr/bin/env Rscript
# Primary question: does plate identity or within-plate position explain
# strain-to-strain-replicate variance in color (L*, a*, b*) or colony
# size/morphology?
#
# Each strain is measured once per run (= "replicate", but see README quirk
# #3: replicate_label is really an imager-run proxy, not a true technical
# replicate) at a given Copper concentration, so a strain's 4 observations at
# one concentration sit on 4 different physical plates/wells. That gives a
# genuine (if small, n=4) within-strain sample to ask: how much of that
# spread is plate-of-origin / grid position, vs strain identity, vs residual
# noise?
#
# DESIGN NOTE: every plate was imaged at exactly one Copper concentration
# (plate_id is nested inside copper_mm -- no plate carries multiple Cu
# levels). A single pooled model would therefore average the plate variance
# component across Cu conditions that have very different biology (0 mM control
# vs 30 mM stress). So the analysis is STRATIFIED BY COPPER CONCENTRATION: the
# variance partition is fit independently within each copper_mm level, with
# the copper term dropped (invariant within a stratum).
#
# Variance-components approach (crossed/nested random effects via lme4):
#   trait ~ grid_row_c + grid_col_c + is_edge
#           + (1 | strain_code) + (1 | run_number) + (1 | plate_id)
# strain_code and plate_id are non-nested (strain moves plate-to-plate across
# runs), run_number is a batch effect, plate_id is nested in run_number.
# grid_row / grid_col enter as fixed effects (a monotonic edge/center
# gradient, e.g. evaporation, is a more plausible mechanism than an effect
# that's constant within a row/col level, so numeric position is the primary
# model; a fully categorical version is fit alongside as a robustness check).
#
# Run from the repo root:
#   pixi run Rscript analysis/explore_plate_position/scripts/01_variance_partition.R

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
})

tab_dir <- "analysis/explore_plate_position/results/tables"
fig_dir <- "analysis/explore_plate_position/results/figures"
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

endpoint <- readRDS(file.path(tab_dir, "endpoint_colonies.rds"))

traits <- c(
  lab_l_mean     = "L* (lightness)",
  lab_a_mean     = "a* (green-red)",
  lab_b_mean     = "b* (blue-yellow)",
  shape_area     = "Colony area",
  shape_solidity = "Solidity (morphology)",
  shape_eccentricity = "Eccentricity (morphology)"
)

copper_levels <- sort(unique(endpoint$copper_mm))

# Restrict to strains observed on >= 2 plates. A strain seen once can't
# inform a within-strain variance-component estimate, and including it just
# adds noise to the (1 | strain_code) term. The filter is applied *within each
# copper stratum* so it adapts to each condition's realized design.
n_plates_per_strain <- function(d) d %>%
  distinct(strain_code, plate_id) %>%
  count(strain_code, name = "n_plates")
usable_strains <- function(d) n_plates_per_strain(d) %>% filter(n_plates >= 2) %>% pull(strain_code)

dat_by_copper <- lapply(copper_levels, function(cu) {
  d <- endpoint %>%
    filter(copper_mm == cu, strain_code %in% usable_strains(.)) %>%
    mutate(
      grid_row_c = grid_row - mean(grid_row),
      grid_col_c = grid_col - mean(grid_col)
    )
  d
})
names(dat_by_copper) <- copper_levels

for (cu in copper_levels) {
  d <- dat_by_copper[[as.character(cu)]]
  cat(sprintf(
    "Stratum Cu=%d: %d colonies, %d strains (each on >=2 plates %s), %d plates, %d runs\n",
    cu, nrow(d), n_distinct(d$strain_code),
    sprintf("within Cu=%d", cu),
    n_distinct(d$plate_id), n_distinct(d$run_number)
  ))
}

fit_one <- function(trait, dat) {
  form <- as.formula(sprintf(
    "%s ~ grid_row_c + grid_col_c + is_edge + (1 | strain_code) + (1 | run_number) + (1 | plate_id)",
    trait
  ))
  lmer(form, data = dat, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
}

fit_categorical <- function(trait, dat) {
  form <- as.formula(sprintf(
    "%s ~ grid_row_f + grid_col_f + (1 | strain_code) + (1 | run_number) + (1 | plate_id)",
    trait
  ))
  lmer(form, data = dat, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
}

# models[[cu]][[trait]]; ordered for legible output
models    <- lapply(copper_levels, function(cu) {
  setNames(lapply(names(traits), fit_one, dat = dat_by_copper[[as.character(cu)]]), names(traits))
})
names(models) <- copper_levels
cat_models <- lapply(copper_levels, function(cu) {
  setNames(lapply(names(traits), fit_categorical, dat = dat_by_copper[[as.character(cu)]]), names(traits))
})
names(cat_models) <- copper_levels

# --- Variance components table: per Cu concentration, how much of the total
# variance sits at the strain level vs run (batch) vs plate
# (position/microenvironment) vs residual (within-plate, within-strain noise)?
vc_table <- bind_rows(lapply(copper_levels, function(cu) {
  bind_rows(lapply(names(models[[as.character(cu)]]), function(tr) {
    vc <- as.data.frame(VarCorr(models[[as.character(cu)]][[tr]]))
    vc %>%
      transmute(
        trait = traits[[tr]],
        copper_mm = cu,
        group = grp,
        variance = vcov,
        sd = sdcor
      ) %>%
      mutate(pct_of_total = 100 * variance / sum(variance))
  }))
}))
write.csv(vc_table, file.path(tab_dir, "variance_components_by_copper.csv"), row.names = FALSE)

# --- Fixed-effects table: does grid position (row/col/edge) have a detectable
# effect once strain/run/plate are accounted for, within each Cu concentration,
# with lmerTest's Satterthwaite-approximated p-values?
fe_table <- bind_rows(lapply(copper_levels, function(cu) {
  bind_rows(lapply(names(models[[as.character(cu)]]), function(tr) {
    tidy(models[[as.character(cu)]][[tr]], effects = "fixed") %>%
      mutate(trait = traits[[tr]], copper_mm = cu, .before = 1)
  }))
}))
write.csv(fe_table, file.path(tab_dir, "fixed_effects_by_copper.csv"), row.names = FALSE)

# --- Categorical-position robustness check: row/col as factors instead of
# linear terms, to catch a non-monotonic (e.g. single problem row) effect
# that the linear model above would average away. Fit per Cu concentration.
anova_table <- bind_rows(lapply(copper_levels, function(cu) {
  bind_rows(lapply(names(cat_models[[as.character(cu)]]), function(tr) {
    a <- anova(cat_models[[as.character(cu)]][[tr]])
    as.data.frame(a) %>%
      tibble::rownames_to_column("term") %>%
      mutate(trait = traits[[tr]], copper_mm = cu, .before = 1)
  }))
}))
write.csv(anova_table, file.path(tab_dir, "categorical_position_anova_by_copper.csv"), row.names = FALSE)

# Old pooled-model files are no longer produced; remove any leftover from a
# previous run so results/tables stays self-consistent.
unlink(file.path(tab_dir, c(
  "variance_components.csv", "fixed_effects.csv", "categorical_position_anova.csv",
  "lmer_models_linear.rds", "lmer_models_categorical.rds"
)))

saveRDS(models, file.path(tab_dir, "lmer_models_by_copper_linear.rds"))
saveRDS(cat_models, file.path(tab_dir, "lmer_models_by_copper_categorical.rds"))

cat("\nPlate variance (% of total) per trait x Cu concentration:\n")
print(vc_table %>%
        filter(group == "plate_id") %>%
        select(trait, copper_mm, pct_of_total) %>%
        tidyr::pivot_wider(names_from = copper_mm, values_from = pct_of_total))

cat("\nWrote variance_components_by_copper.csv, fixed_effects_by_copper.csv, categorical_position_anova_by_copper.csv to", tab_dir, "\n")
