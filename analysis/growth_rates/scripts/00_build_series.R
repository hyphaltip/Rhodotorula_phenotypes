#!/usr/bin/env Rscript
# Build the per-colony x per-timepoint growth series for the growth-rate
# analysis, with TWO size readouts sampled from v_phenotype:
#
#   shape_area              -- CellProfiler Shape_Area, biomass/expansion proxy
#   intensity_mean          -- CellProfiler Intensity_MeanIntensity, colony
#                              brightness/bulk density. Note it saturates early
#                              (~48-54 h) while area keeps rising to 120 h, so
#                              the two readouts cover different growth windows.
#
# Replicates the exact conventions of the plate-position time course
# (explore_plate_position/scripts/04_build_timecourse.R):
#   - experiment 'Copper', runs 353-356 (run 357 excluded),
#   - strain_id IS NOT NULL,
#   - hours_since_plate_start > 0 (skip the t=0 plating image),
#   - one row per (run, plate, well, timepoint); when a well has >1 segmented
#     object keep the larger by shape_area,
#   - plate_id = "<run>_<plate_number>" so colonies are physically distinct.
#
# Additionally joins a species label per strain (from the `strain` table) so
# the between-species growth analysis has a grouping key.
#
# Run from the repo root:
#   pixi run Rscript analysis/growth_rates/scripts/00_build_series.R

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(dplyr)
})

out_dir <- "analysis/growth_rates/results/tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

con <- dbConnect(duckdb(), dbdir = "db/rhodotorula_phenotypes.duckdb", read_only = TRUE)

series <- dbGetQuery(con, "
  SELECT
      p.run_number,
      p.plate_number,
      p.factors['Copper concentration'] AS copper_mm,
      p.strain_id,
      p.strain_code,
      p.well_position,
      p.hours_since_plate_start,
      p.\"Shape_Area\"              AS shape_area,
      p.\"Intensity_MeanIntensity\" AS intensity_mean,
      p.\"ColorLab_L*Mean\"         AS lab_l_mean,
      p.\"ColorLab_a*Mean\"         AS lab_a_mean,
      p.\"ColorLab_b*Mean\"         AS lab_b_mean
  FROM v_phenotype p
  WHERE p.experiment_name = 'Copper'
    AND p.run_number NOT IN (357)
    AND p.strain_id IS NOT NULL
    AND p.hours_since_plate_start > 0
")
# Species label: best-effort from the strain table.
sp <- dbGetQuery(con, "
  SELECT DISTINCT strain_id, strain_code, species FROM strain
")
dbDisconnect(con, shutdown = TRUE)

series <- series %>%
  mutate(
    run_number   = factor(run_number),
    plate_id     = paste(run_number, plate_number, sep = "_"),
    strain_code  = factor(strain_code),
    time_idx     = round(hours_since_plate_start / 6),
    hours_snap   = round(hours_since_plate_start / 6) * 6,
    across(c(shape_area, intensity_mean, lab_l_mean, lab_a_mean, lab_b_mean),
           ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  left_join(sp, by = c("strain_id", "strain_code"))

# Keep the largest segmented object per (plate, well, timepoint).
dup <- series %>% count(plate_id, well_position, time_idx) %>% filter(n > 1)
if (nrow(dup) > 0) {
  cat(sprintf("NOTE: %d (plate, well, timepoint) keys had >1 object; keeping larger by shape_area\n", nrow(dup)))
}
series <- series %>%
  group_by(plate_id, well_position, time_idx) %>%
  slice_max(shape_area, n = 1, with_ties = FALSE) %>%
  ungroup()

# Drop any row where the primary growth readout is unusable.
series <- series %>% filter(!is.na(shape_area) & shape_area > 0)

saveRDS(series, file.path(out_dir, "colony_growth_series.rds"))
write.csv(series, file.path(out_dir, "colony_growth_series.csv"), row.names = FALSE)

# --- Coverage summary: how many timepoints / colonies are usable per Cu.
cov <- series %>%
  count(plate_id, well_position, copper_mm, name = "n_tp") %>%
  count(copper_mm, n_tp) %>%
  rename(n_colonies = n)
write.csv(cov, file.path(out_dir, "series_coverage.csv"), row.names = FALSE)

cat(sprintf(
  "Wrote %d colony-timepoint rows, %d colonies, %d strains, %d species, %d plates to %s\n",
  nrow(series),
  n_distinct(series$plate_id, series$well_position),
  n_distinct(series$strain_code),
  n_distinct(series$species[!is.na(series$species)]),
  n_distinct(series$plate_id),
  out_dir
))
cat("\nColonies by usable-timepoint count per Cu concentration:\n")
print(as.data.frame(cov))
cat("\nSpecies N (strains):\n")
print(as.data.frame(series %>%
  distinct(strain_code, species) %>%
  count(species, name = "n_strains") %>%
  arrange(desc(n_strains))), max = 40)
