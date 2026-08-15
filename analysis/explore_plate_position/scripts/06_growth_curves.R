#!/usr/bin/env Rscript
# Time-component part 2: growth-curve mixed models.
#
# Question: do colonies get brighter (L* up) and bigger (area up) over the
# plate time course, does the rate depend on Cu, and do strains differ in
# their growth trajectories?
#
# Every colony is tracked across its plate's full 6-hourly relative clock
# (t=0 .. ~114 h, run 353 to ~117 h). We treat time as continuous (days) and
# fit, per trait:
#
#   trait ~ (time_c + time_c^2) * copper
#          + (1 + time_c | strain_code)      # strain curve offsets & rate diffs
#          + (1 | plate_id) + (1 | colony_id) # plate batch + repeated colony
#
# where colony_id = plate_id_well. copper is a factor for the 7 tested levels.
# time_c is centered at day 4 to stabilize the quadratic.
#
# Missingness caveat: at high Cu, strains that fail to grow are simply absent,
# so late-time estimates come from the surviving (tolerant) subset -- we report
# realized n per (Cu, timepoint) and show it as a detection curve.
#
# Run from the repo root:
#   pixi run Rscript analysis/explore_plate_position/scripts/06_growth_curves.R

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

tc <- tc %>%
  mutate(
    copper = factor(copper_mm, levels = sort(unique(copper_mm))),
    time_c = hours_d - 4,                 # centered days (median ~ day 4)
    colony_id = paste(plate_id, well_position, sep = "_"),
    log_area = log(shape_area)
  )

# Detection curve: fraction of (strain x plate) combos present at each age per
# Cu, versus the max seen for that Cu. Declines at late ages = strains not
# growing / dropping out (informative missingness).
det_curve <- tc %>%
  group_by(copper, hours_snap) %>%
  summarise(n_colonies = n(), .groups = "drop") %>%
  group_by(copper) %>%
  mutate(n_max = max(n_colonies), frac_detected = n_colonies / n_max) %>%
  ungroup()
write.csv(det_curve, file.path(tab_dir, "detection_curve_by_copper.csv"), row.names = FALSE)

fit_growth <- function(trait, formula_suffix) {
  form <- as.formula(sprintf(
    "%s ~ %s + (1 + time_c | strain_code) + (1 | plate_id) + (1 | colony_id)",
    trait, formula_suffix
  ))
  lmer(form, data = tc, REML = TRUE,
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
}

cat("Fitting L* growth curve (ranks may take a minute)...\n")
m_L   <- fit_growth("lab_l_mean", "(time_c + I(time_c^2)) * copper")
cat("Fitting log(area) growth curve...\n")
m_A   <- fit_growth("log_area", "(time_c + I(time_c^2)) * copper")

growth_fe <- bind_rows(
  tidy(m_L, effects = "fixed") %>% mutate(trait = "L* (lightness)", .before = 1),
  tidy(m_A, effects = "fixed") %>% mutate(trait = "log Colony area", .before = 1)
)
write.csv(growth_fe, file.path(tab_dir, "growth_curve_fixed_effects.csv"), row.names = FALSE)

growth_re <- bind_rows(
  as.data.frame(VarCorr(m_L)) %>% transmute(group = grp, variance = vcov, sd = sdcor) %>%
    mutate(trait = "L* (lightness)", .before = 1),
  as.data.frame(VarCorr(m_A)) %>% transmute(group = grp, variance = vcov, sd = sdcor) %>%
    mutate(trait = "log Colony area", .before = 1)
)
write.csv(growth_re, file.path(tab_dir, "growth_curve_random_effects.csv"), row.names = FALSE)

cat("\nGrowth-curve fixed effects:\n")
print(growth_fe %>% select(-statistic) %>% as.data.frame(), digits = 4)

# --- Figures ---
ggsave_ <- function(name, plot, width = 7, height = 5) {
  ggsave(file.path(fig_dir, name), plot, width = width, height = height, dpi = 150)
}

# Fixed-effect prediction across the time grid for each Cu level (mean strain,
# zero plate/colony random effects) -- the "population" growth trajectory.
predict_fixed <- function(model, trait_name) {
  b <- fixef(model)
  cu_levels <- levels(tc$copper)
  grid_time <- seq(6, 114, by = 6) / 24
  # intercept = b["(Intercept)"]; copper k effect = b[paste0("copper", lvl)];
  # time effects are centered at day 4.
  preds <- lapply(cu_levels, function(cu) {
    time_c <- grid_time - 4
    t1 <- b["time_c"]; t2 <- b["I(time_c^2)"]
    eff_cu <- if (cu == cu_levels[1]) 0 else b[paste0("copper", cu)]
    eff_t1 <- b["time_c"]
    eff_cut1 <- if (cu == cu_levels[1]) 0 else b[paste0("time_c:copper", cu)]
    eff_t2 <- b["I(time_c^2)"]
    eff_cut2 <- if (cu == cu_levels[1]) 0 else b[paste0("I(time_c^2):copper", cu)]
    y <- b["(Intercept)"] + b["time_c"] * time_c + b["I(time_c^2)"] * time_c^2 +
      eff_cu + eff_cut1 * time_c + eff_cut2 * time_c^2
    data.frame(copper = cu, hours_d = grid_time, pred = y)
  })
  bind_rows(preds)
}

obs_summ <- tc %>%
  group_by(copper, hours_snap) %>%
  summarise(L_med = median(lab_l_mean), area_med = median(shape_area),
            logA_med = median(log_area), .groups = "drop")

pred_L <- predict_fixed(m_L, "L*")
pL <- ggplot() +
  geom_point(data = obs_summ, aes(hours_snap / 24, L_med, color = copper),
             alpha = 0.5, size = 1) +
  geom_line(data = pred_L, aes(hours_d, pred, color = copper)) +
  scale_color_viridis_d() +
  labs(x = "Hours since plate start / 24 (days)", y = "Median L* (lightness)",
       color = "Cu (mM)", title = "Colony brightness over time, by Cu concentration") +
  theme_minimal(base_size = 11)
ggsave_("growth_L_by_cu.png", pL, width = 9, height = 6)

pred_A <- predict_fixed(m_A, "log area")
pA <- ggplot() +
  geom_point(data = obs_summ, aes(hours_snap / 24, logA_med, color = copper),
             alpha = 0.5, size = 1) +
  geom_line(data = pred_A, aes(hours_d, pred, color = copper)) +
  scale_color_viridis_d() +
  labs(x = "Hours since plate start / 24 (days)", y = "Median log(colony area)",
       color = "Cu (mM)", title = "Colony size over time, by Cu concentration") +
  theme_minimal(base_size = 11)
ggsave_("growth_area_by_cu.png", pA, width = 9, height = 6)

# Detection (informative missingness) curve.
pD <- det_curve %>%
  mutate(copper = factor(copper, labels = paste0(levels(copper), " mM Cu"))) %>%
  ggplot(aes(hours_snap / 24, frac_detected, color = copper)) +
  geom_line() + geom_point(size = 1) +
  scale_color_viridis_d() +
  labs(x = "Hours since plate start / 24 (days)", y = "Fraction of (strain x plate) combos detected",
       color = "Cu (mM)",
       title = "Colony detection over time -- strains drop out of high-Cu late strata") +
  theme_minimal(base_size = 11)
ggsave_("detection_curve_by_cu.png", pD, width = 9, height = 6)

cat("\nWrote growth_curve_fixed_effects.csv, growth_curve_random_effects.csv, detection_curve_by_copper.csv to", tab_dir, "\n")
cat("Wrote growth_L_by_cu.png, growth_area_by_cu.png, detection_curve_by_cu.png to", fig_dir, "\n")
