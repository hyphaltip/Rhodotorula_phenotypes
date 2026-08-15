#!/usr/bin/env Rscript
# Minimal example: query the phenotype database from R.
#
# Run from the repo root:
#   pixi run Rscript scripts/db/query_examples/query_example.R

library(DBI)
library(duckdb)

con <- dbConnect(duckdb(), dbdir = "db/rhodotorula_phenotypes.duckdb", read_only = TRUE)

strain_code <- "TFCN_152C-3"

# 1. Every colony observation for one strain, across every condition and
#    replicate it appears in.
obs <- dbGetQuery(con, "
    SELECT strain_code, replicate_label, hours_since_plate_start,
           Shape_Area, Intensity_MeanIntensity
    FROM v_phenotype
    WHERE strain_code = ?
    ORDER BY hours_since_plate_start
", params = list(strain_code))
cat(sprintf("%s: %d colony observations\n", strain_code, nrow(obs)))
print(head(obs))

# 2. Mean colony size per Copper concentration for that strain (pull the
#    level straight out of the `factors` MAP).
summary_df <- dbGetQuery(con, "
    SELECT factors['Copper concentration'] AS concentration_mm,
           COUNT(*) AS n, AVG(Shape_Area) AS mean_area
    FROM v_phenotype
    WHERE strain_code = ? AND experiment_name = 'Copper'
    GROUP BY 1
    ORDER BY 1
", params = list(strain_code))
cat("\nMean colony area by Copper concentration:\n")
print(summary_df)

dbDisconnect(con, shutdown = TRUE)
