#!/usr/bin/env Rscript
# Secondary question: beyond a colony's own grid position, do its immediate
# neighbors on the same plate/image predict its trait value? A significant
# neighbor effect after controlling for the focal colony's own strain, run,
# plate, and grid position would point to a local microenvironment or
# colony-colony interaction effect (e.g. shared moisture/copper diffusion,
# optical bleed into neighboring segmentation masks) distinct from a
# monotonic plate-position gradient.
#
# Neighbors are defined by grid adjacency (von Neumann: up/down/left/right)
# on the SAME plate_id (same physical plate, same endpoint image) -- adjacent
# grid_row/grid_col, not adjacent well_position (well_position is column-major
# and numerically adjacent positions are not spatially adjacent).
#
# Run from the repo root:
#   pixi run Rscript analysis/explore_plate_position/scripts/02_adjacency_effect.R

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
})

tab_dir <- "analysis/explore_plate_position/results/tables"
endpoint <- readRDS(file.path(tab_dir, "endpoint_colonies.rds"))

traits <- c(
  lab_l_mean     = "L* (lightness)",
  lab_a_mean     = "a* (green-red)",
  lab_b_mean     = "b* (blue-yellow)",
  shape_area     = "Colony area",
  shape_solidity = "Solidity (morphology)"
)

# Build the neighbor-mean covariate per (plate_id, colony) for each trait.
# Self join within plate_id on |drow|+|dcol| == 1 (4-connected grid).
neighbor_means <- endpoint %>%
  select(plate_id, well_position, grid_row, grid_col, all_of(names(traits))) %>%
  {
    focal <- rename_with(., ~ paste0(.x, "_focal"), c(well_position, grid_row, grid_col))
    cand  <- rename_with(., ~ paste0(.x, "_nb"), c(well_position, grid_row, grid_col, all_of(names(traits))))
    inner_join(focal, cand, by = "plate_id", relationship = "many-to-many")
  } %>%
  filter(
    (abs(grid_row_focal - grid_row_nb) + abs(grid_col_focal - grid_col_nb)) == 1
  ) %>%
  group_by(plate_id, well_position_focal) %>%
  summarise(
    n_neighbors = n(),
    across(paste0(names(traits), "_nb"), mean, .names = "{.col}_mean"),
    .groups = "drop"
  ) %>%
  rename(well_position = well_position_focal) %>%
  rename_with(~ sub("_nb_mean$", "_neighbor_mean", .x), ends_with("_nb_mean"))

dat <- endpoint %>%
  inner_join(neighbor_means, by = c("plate_id", "well_position")) %>%
  filter(n_neighbors >= 2)  # corner/edge wells with <2 neighbors give a noisy 1-colony average

cat(sprintf(
  "Adjacency dataset: %d colonies with >=2 grid-neighbors on their plate (of %d total endpoint colonies)\n",
  nrow(dat), nrow(endpoint)
))

fit_one <- function(trait) {
  nb_col <- paste0(trait, "_neighbor_mean")
  form <- as.formula(sprintf(
    "%s ~ copper_mm + is_edge + %s + (1 | strain_code) + (1 | run_number) + (1 | plate_id)",
    trait, nb_col
  ))
  lmer(form, data = dat, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
}

models <- setNames(lapply(names(traits), fit_one), names(traits))

fe_table <- bind_rows(lapply(names(models), function(tr) {
  tidy(models[[tr]], effects = "fixed") %>% mutate(trait = traits[[tr]], .before = 1)
}))
write.csv(fe_table, file.path(tab_dir, "adjacency_fixed_effects.csv"), row.names = FALSE)

saveRDS(models, file.path(tab_dir, "lmer_models_adjacency.rds"))
saveRDS(dat, file.path(tab_dir, "adjacency_dataset.rds"))

cat("\nNeighbor-mean coefficient (per trait; the row 'traitname_neighbor_mean'):\n")
print(fe_table %>%
  filter(grepl("_neighbor_mean$", term)) %>%
  select(trait, term, estimate, std.error, statistic, df, p.value))

cat("\nWrote adjacency_fixed_effects.csv to", tab_dir, "\n")
