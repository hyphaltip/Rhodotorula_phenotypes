#!/usr/bin/env Rscript
# Build the analysis-ready colony-endpoint table used by every downstream
# script in this analysis.
#
# For every colony_measurement row we keep only the *last* image per plate
# (the most-developed timepoint, ~117h) -- colony size/color are still
# changing earlier in the time course, so mixing timepoints would confound
# growth stage with plate-position effects.
#
# Run from the repo root:
#   pixi run Rscript analysis/explore_plate_position/scripts/00_build_dataset.R

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(dplyr)
})

out_dir <- "analysis/explore_plate_position/results/tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

con <- dbConnect(duckdb(), dbdir = "db/rhodotorula_phenotypes.duckdb", read_only = TRUE)

# One row per (run_number, plate_number) giving the last imaged timepoint,
# then join back to v_phenotype restricted to that timepoint. Done in SQL so
# duckdb -- not R -- does the heavy lifting over the ~212k-row fact table.
endpoint <- dbGetQuery(con, "
  WITH last_image AS (
      SELECT run_number, plate_number, MAX(hours_since_plate_start) AS t_end
      FROM v_phenotype
      WHERE experiment_name = 'Copper'
      GROUP BY run_number, plate_number
  )
  SELECT
      p.run_number,
      p.replicate_label,      -- imager-run proxy, NOT a true technical replicate (README quirk #3)
      p.plate_number,
      p.factors['Copper concentration'] AS copper_mm,
      p.strain_id,
      p.strain_code,
      p.well_position,        -- 1..96, column-major (Grid_ColMajorIdx + 1)
      p.\"Grid_RowNum\"          AS grid_row,   -- 0..7
      p.\"Grid_ColNum\"          AS grid_col,   -- 0..11
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
  JOIN last_image li
    ON li.run_number = p.run_number
   AND li.plate_number = p.plate_number
   AND li.t_end = p.hours_since_plate_start
  WHERE p.experiment_name = 'Copper'
    AND p.strain_id IS NOT NULL   -- drop unmatched wells (README quirk #6 gaps)
")

dbDisconnect(con, shutdown = TRUE)

# Plate-position derived columns, all computed once here so every downstream
# script uses the same definitions.
endpoint <- endpoint %>%
  mutate(
    row_edge  = grid_row %in% c(0, 7),
    col_edge  = grid_col %in% c(0, 11),
    is_edge   = row_edge | col_edge,
    # Chebyshev distance from the plate border (0 = on the border).
    edge_dist = pmin(grid_row, 7 - grid_row, grid_col, 11 - grid_col),
    grid_row_f = factor(grid_row),
    grid_col_f = factor(grid_col),
    run_number = factor(run_number),
    plate_id   = paste(run_number, plate_number, sep = "_"),
    strain_code = factor(strain_code)
  )

# A handful of wells (3 in the full Copper set) have two CellProfiler objects
# at the endpoint timepoint -- almost certainly one real colony plus a small
# segmentation fragment, not two independent colonies (a 96-well plate has
# one strain per well by design; see well_placement). Keep the larger object
# and drop the fragment, rather than silently letting duplicate rows inflate
# n or averaging a real colony together with a segmentation artifact.
dup_keys <- endpoint %>%
  count(plate_id, well_position) %>%
  filter(n > 1)
if (nrow(dup_keys) > 0) {
  cat(sprintf(
    "NOTE: %d well(s) had >1 segmented object at the endpoint timepoint; keeping the larger object (by shape_area) at each:\n",
    nrow(dup_keys)
  ))
  print(dup_keys)
}
endpoint <- endpoint %>%
  group_by(plate_id, well_position) %>%
  slice_max(shape_area, n = 1, with_ties = FALSE) %>%
  ungroup()

saveRDS(endpoint, file.path(out_dir, "endpoint_colonies.rds"))
write.csv(endpoint, file.path(out_dir, "endpoint_colonies.csv"), row.names = FALSE)

cat(sprintf(
  "Wrote %d endpoint colony rows (%d plates, %d strains, %d runs) to %s\n",
  nrow(endpoint),
  n_distinct(endpoint$plate_id),
  n_distinct(endpoint$strain_code),
  n_distinct(endpoint$run_number),
  out_dir
))
