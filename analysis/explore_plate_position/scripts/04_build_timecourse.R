#!/usr/bin/env Rscript
# Build the per-colony x per-timepoint growth table used by the time-component
# analyses (05_time_variance_partition.R, 06_growth_curves.R).
#
# The endpoint analysis (00_build_dataset.R) keeps only each plate's *last*
# image so growth stage cannot confound plate-position effects. Here we go the
# other way and keep the FULL per-plate time course:
#
#   - every plate is imaged on its own 6-hourly relative clock (t = 0 at plate
#     start, then ~6, 12, 18, ... up to ~114 h; run 353 extends to ~117 h), so
#     hours_since_plate_start is the plate's own age at each image, NOT wall
#     clock;
#   - a "day block" (d1 = 0-24 h, ..., d5 = 96-120 h) is therefore the same
#     developmental stage for every plate, which is what makes block-stratified
#     or growth-curve modelling meaningful;
#   - runs 353-356 only (run 357, plates 113-120, is excluded for consistency
#     with the endpoint table until its provenance is pinned down).
#
# Same conventions as 00/01/03: strain_id IS NOT NULL drops unmatched wells;
# when a well yields >1 segmented object at a timepoint keep the larger by
# shape_area; position columns are re-computed identically.
#
# Run from the repo root:
#   pixi run Rscript analysis/explore_plate_position/scripts/04_build_timecourse.R

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(dplyr)
})

out_dir <- "analysis/explore_plate_position/results/tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

con <- dbConnect(duckdb(), dbdir = "db/rhodotorula_phenotypes.duckdb", read_only = TRUE)

# One row per (run, plate, well, timepoint) for the Copper experiment, all
# imaging rounds of each plate (not just the last one).
tc <- dbGetQuery(con, "
  SELECT
      p.run_number,
      p.replicate_label,
      p.plate_number,
      p.factors['Copper concentration'] AS copper_mm,
      p.strain_id,
      p.strain_code,
      p.well_position,
      p.\"Grid_RowNum\"          AS grid_row,
      p.\"Grid_ColNum\"          AS grid_col,
      p.hours_since_plate_start,
      p.\"Shape_Area\"            AS shape_area,
      p.\"Shape_Perimeter\"       AS shape_perimeter,
      p.\"Shape_Eccentricity\"    AS shape_eccentricity,
      p.\"Shape_Solidity\"        AS shape_solidity,
      p.\"Shape_Circularity\"     AS shape_circularity,
      p.\"Shape_MajorAxisLength\" AS shape_major_axis,
      p.\"Shape_MinorAxisLength\" AS shape_minor_axis,
      p.\"ColorLab_L*Mean\"       AS lab_l_mean,
      p.\"ColorLab_a*Mean\"       AS lab_a_mean,
      p.\"ColorLab_b*Mean\"       AS lab_b_mean,
      p.\"ColorLab_L*StdDev\"     AS lab_l_sd,
      p.\"ColorLab_a*StdDev\"     AS lab_a_sd,
      p.\"ColorLab_b*StdDev\"     AS lab_b_sd
  FROM v_phenotype p
  WHERE p.experiment_name = 'Copper'
    AND p.run_number NOT IN (357)       -- excluded: provenance not yet pinned down
    AND p.strain_id IS NOT NULL
    AND p.hours_since_plate_start > 0   -- t=0 is the plating image, no growth yet
")
dbDisconnect(con, shutdown = TRUE)

tc <- tc %>%
  mutate(
    row_edge  = grid_row %in% c(0, 7),
    col_edge  = grid_col %in% c(0, 11),
    is_edge   = row_edge | col_edge,
    # Chebyshev distance from the plate border (0 = on the border).
    edge_dist = pmin(grid_row, 7 - grid_row, grid_col, 11 - grid_col),
    grid_row_f = factor(grid_row),
    grid_col_f = factor(grid_col),
    run_number  = factor(run_number),
    plate_id    = paste(run_number, plate_number, sep = "_"),
    strain_code = factor(strain_code),
    # Snap each image to the nearest 6 h ring so the (small, ~0.01-0.1 h)
    # per-plate jitter and run 353's 3 h offset map onto the shared grid.
    time_idx    = round(hours_since_plate_start / 6),
    hours_snap  = round(hours_since_plate_start / 6) * 6,
    hours_d     = hours_snap / 24,
    # Same-day-block rule for every plate (0 < t <= 24h is "day 1").
    day_block   = cut(hours_snap,
                      breaks = c(-Inf, 24, 48, 72, 96, Inf),
                      labels = c("d1 0-24h", "d2 24-48h", "d3 48-72h",
                                 "d4 72-96h", "d5 96+h")),
    across(c(lab_l_mean, lab_a_mean, lab_b_mean, shape_area),
           ~ suppressWarnings(as.numeric(.x)))
  )

# A well can contain >1 segmented object at a given timepoint (one colony plus
# a small fragmentation/aggregation artifact). Keep the larger object per
# (plate, well, timepoint), matching the endpoint-table convention.
dup_keys <- tc %>% count(plate_id, well_position, time_idx) %>% filter(n > 1)
if (nrow(dup_keys) > 0) {
  cat(sprintf(
    "NOTE: %d (plate, well, timepoint) keys had >1 segmented object; keeping the larger (by shape_area) at each\n",
    nrow(dup_keys)
  ))
}
tc <- tc %>%
  group_by(plate_id, well_position, time_idx) %>%
  slice_max(shape_area, n = 1, with_ties = FALSE) %>%
  ungroup()

saveRDS(tc, file.path(out_dir, "colony_timecourse.rds"))
write.csv(tc, file.path(out_dir, "colony_timecourse.csv"), row.names = FALSE)

# --- Missingness summary. Not every strain grows on every plate / Cu: some
# colonies are never detected at a given concentration, and at high Cu a real
# (smaller) share never reach the late timepoints at all. That dropout is
# biologically informative (growth failure under Cu stress is itself a
# phenotype) but means the growth/variance models see a selected subset at
# high Cu. Report the realized detection so downstream strata know their n.
miss <- tc %>%
  group_by(strain_code, copper_mm, plate_id) %>%
  summarise(tmax = max(hours_snap), .groups = "drop") %>%
  mutate(reached_late = tmax >= 100)
miss_sum <- miss %>%
  group_by(copper_mm) %>%
  summarise(
    n_strain_plates      = n(),
    n_reached_late       = sum(reached_late),
    pct_reached_late     = round(100 * mean(reached_late), 1),
    n_strains_seen       = n_distinct(strain_code),
    n_strains_reached_late = n_distinct(strain_code[reached_late]),
    .groups = "drop"
  ) %>%
  arrange(copper_mm)
write.csv(miss_sum, file.path(out_dir, "growth_missingness_by_copper.csv"), row.names = FALSE)

cat(sprintf(
  "Wrote %d colony-timepoint rows (%d plates, %d strains, %d runs) to %s\n",
  nrow(tc),
  n_distinct(tc$plate_id),
  n_distinct(tc$strain_code),
  n_distinct(tc$run_number),
  out_dir
))
cat("Detection by copper concentration (%% of strain-plates reaching >=100 h):\n")
print(miss_sum)
cat("Day blocks:\n")
print(tc %>% count(day_block, time_idx))
