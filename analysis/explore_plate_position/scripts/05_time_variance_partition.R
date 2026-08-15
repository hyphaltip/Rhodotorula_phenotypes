#!/usr/bin/env Rscript
# Time-component part 1: does the share of trait variance explained by plate
# identity / position change across developmental stage?
#
# Each plate is imaged on its own relative 6-hourly clock from t=0 to ~114 h
# (run 353 to ~117 h), so the same "day block" (0-24h = d1, 24-48h = d2, ...,
# 96-120h = d5) is the same growth stage on every plate. We re-run the exact
# endpoint variance partition (01_variance_partition.R) on each (day block x
# copper concentration) stratum independently, dropping the copper term
# (invariant within a stratum) exactly as the endpoint analysis does.
#
# Missingness caveat (see 04_build_timecourse.R): strains that do not grow at
# high Cu drop out of the later strata. The >=2-plates-per-strain filter runs
# within each stratum, and every output row carries the realized n_colonies /
# n_strains / n_plates / n_runs so the report can state how big (and how
# selected) each stratum actually was.
#
# Run from the repo root:
#   pixi run Rscript analysis/explore_plate_position/scripts/05_time_variance_partition.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(ggplot2)
})

tab_dir <- "analysis/explore_plate_position/results/tables"
fig_dir <- "analysis/explore_plate_position/results/figures"
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

tc <- readRDS(file.path(tab_dir, "colony_timecourse.rds"))

traits <- c(
  lab_l_mean     = "L* (lightness)",
  lab_a_mean     = "a* (green-red)",
  lab_b_mean     = "b* (blue-yellow)",
  shape_area     = "Colony area",
  shape_solidity = "Solidity (morphology)",
  shape_eccentricity = "Eccentricity (morphology)"
)

day_levels <- levels(tc$day_block)
copper_levels <- sort(unique(tc$copper_mm))

n_plates_per_strain <- function(d) d %>%
  distinct(strain_code, plate_id) %>%
  count(strain_code, name = "n_plates")
usable_strains <- function(d) n_plates_per_strain(d) %>% filter(n_plates >= 2) %>% pull(strain_code)

# Stratified by (day_block x copper); within each stratum repeat the endpoint
# fixed-effect structure (linear position gradient + edge indicator + crossed
# strain/run/plate random effects).
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

fit_strata <- function(fitter) {
  out <- vector("list", length(day_levels) * length(copper_levels))
  k <- 1L
  for (db in day_levels) {
    for (cu in copper_levels) {
      d <- tc %>%
        filter(day_block == db, copper_mm == cu,
               strain_code %in% usable_strains(.)) %>%
        mutate(grid_row_c = grid_row - mean(grid_row),
               grid_col_c = grid_col - mean(grid_col))
      tag <- sprintf("day=%s Cu=%g", db, cu)
      if (nrow(d) < 50) {
        cat(sprintf("  SKIP %s: only %d rows\n", tag, nrow(d)))
        next
      }
      # force plate/run to be factors present in the slice
      d <- d %>% mutate(run_number = factor(run_number), plate_id = factor(plate_id))
      nd <- list(day_block = db, copper_mm = cu,
                 n_colonies = nrow(d), n_strains = n_distinct(d$strain_code),
                 n_plates = n_distinct(d$plate_id), n_runs = n_distinct(d$run_number),
                 models = setNames(lapply(names(traits), fitter, dat = d), names(traits)))
      out[[k]] <- nd; k <- k + 1L
    }
  }
  out[seq_len(k - 1L)]
}

cat("Fitting linear-position models per (day block x copper)...\n")
linear <- fit_strata(fit_one)
cat("Fitting categorical-position robustness models per (day block x copper)...\n")
categorical <- fit_strata(fit_categorical)

# --- Variance components with realized stratum sizes attached.
vc_table <- bind_rows(lapply(linear, function(st) {
  bind_rows(lapply(st$models, function(m) {
    as.data.frame(VarCorr(m)) %>%
      transmute(group = grp, variance = vcov, sd = sdcor) %>%
      mutate(pct_of_total = 100 * variance / sum(variance))
  }), .id = "trait") %>%
    mutate(trait = unname(traits[.$trait]), .before = 1) %>%
    mutate(day_block = st$day_block, copper_mm = st$copper_mm,
           n_colonies = st$n_colonies, n_strains = st$n_strains,
           n_plates = st$n_plates, n_runs = st$n_runs, .after = "trait")
}))
write.csv(vc_table, file.path(tab_dir, "variance_components_by_day.csv"), row.names = FALSE)

# --- Fixed effects (position gradient) per stratum.
fe_table <- bind_rows(lapply(linear, function(st) {
  bind_rows(lapply(st$models, function(m) tidy(m, effects = "fixed")), .id = "trait") %>%
    mutate(trait = unname(traits[.$trait]), day_block = st$day_block,
           copper_mm = st$copper_mm, n_strains = st$n_strains,
           n_plates = st$n_plates, .before = 1)
}))
write.csv(fe_table, file.path(tab_dir, "fixed_effects_by_day.csv"), row.names = FALSE)

# --- Categorical-position robustness check per stratum.
anova_table <- bind_rows(lapply(categorical, function(st) {
  bind_rows(lapply(st$models, function(m) {
    a <- as.data.frame(anova(m)) %>% tibble::rownames_to_column("term")
  }), .id = "trait") %>%
    mutate(trait = unname(traits[.$trait]), day_block = st$day_block,
           copper_mm = st$copper_mm, n_strains = st$n_strains,
           n_plates = st$n_plates, .before = 1)
}))
write.csv(anova_table, file.path(tab_dir, "categorical_position_anova_by_day.csv"), row.names = FALSE)

# --- Figures.
ggsave_ <- function(name, plot, width = 7, height = 5) {
  ggsave(file.path(fig_dir, name), plot, width = width, height = height, dpi = 150)
}

# Plate-explained variance (% of total) heatmap: day block x trait, one panel
# per Cu concentration. Colors are comparable because the axis is %.
plate_pct <- vc_table %>%
  filter(group == "plate_id") %>%
  mutate(trait_short = sub(" \\(.*", "", trait),
         copper_mm = factor(copper_mm, levels = sort(unique(copper_mm)),
                            labels = paste0(sort(unique(copper_mm)), " mM Cu")))
p_day <- plate_pct %>%
  ggplot(aes(x = day_block, y = trait_short, fill = pct_of_total)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", pct_of_total)), color = "grey15", size = 2.8) +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~ copper_mm, ncol = 4) +
  labs(x = "Day block", y = NULL,
       title = "Plate-explained variance (%) by day block x Cu concentration",
       fill = "% of total variance") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave_("day_plate_pct_heatmap.png", p_day, width = 13, height = 9)

# Full partition (stacked bars) for the first and last day blocks, per Cu.
p_part <- vc_table %>%
  filter(day_block %in% c("d1 0-24h", "d5 96+h")) %>%
  mutate(group = factor(group,
                        levels = c("strain_code", "run_number", "plate_id", "Residual"),
                        labels = c("Strain", "Run (batch)", "Plate", "Residual")),
         copper_mm = factor(copper_mm, levels = sort(unique(copper_mm)),
                            labels = paste0(sort(unique(copper_mm)), " mM Cu"))) %>%
  ggplot(aes(x = trait, y = pct_of_total, fill = group)) +
  geom_col() +
  coord_flip() +
  facet_grid(day_block ~ copper_mm) +
  labs(x = NULL, y = "% of total variance", fill = NULL,
       title = "Variance partition: day block 1 vs 5, per Cu concentration") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave_("day_variance_components_d1_d5.png", p_part, width = 15, height = 8)

cat("\nPlate variance (% of total) by day block and Cu concentration:\n")
print(plate_pct %>%
        select(trait_short, copper_mm, day_block, pct_of_total) %>%
        pivot_wider(names_from = day_block, values_from = pct_of_total) %>%
        as.data.frame() %>%
        arrange(trait_short, copper_mm))

cat("\nWrote variance_components_by_day.csv, fixed_effects_by_day.csv, categorical_position_anova_by_day.csv to", tab_dir, "\n")
