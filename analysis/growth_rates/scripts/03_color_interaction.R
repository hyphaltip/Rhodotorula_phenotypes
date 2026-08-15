#!/usr/bin/env Rscript
# 03_color_interaction.R
# Interaction between growth-rate parameters and color / light-intensity.
#
# Outcomes (endpoint = late image, matching the plate-position analysis):
#   lab_l_mean  -- CIELAB L* (lightness) == the "light intensity" readout
#   lab_a_mean  -- CIELAB a* (red/green)
#   lab_b_mean  -- CIELAB b* (yellow/blue)
# Predictors: per-colony growth-rate params (rate_area = peak d(log area)/dt,
#             rate_int = peak d(intensity)/dt) and copper (factor).
#
# Design notes:
#  - Intensity_MeanIntensity is r ~ 0.99 collinear with L* (they are the same
#    axis), so intensity is used as a *growth* readout (rate_int) and L* as an
#    *outcome*; the two are never both predictors.
#  - L* (and a*, b*) are modeled AS OUTCOMES of growth rate x copper, not as
#    predictors (per project decision).
#  - Species + strain enter as random effects (shrinkage).
#
# Run from repo root:
#   pixi run Rscript analysis/growth_rates/scripts/03_color_interaction.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
})

out_tab <- "analysis/growth_rates/results/tables"
out_fig <- "analysis/growth_rates/results/figures"

notes <- character()
note <- function(...) {
  msg <- sprintf(...); cat(msg, "\n"); notes <<- c(notes, msg)
}

series <- readRDS(file.path(out_tab, "colony_growth_series.rds"))
pref   <- read.csv(file.path(out_tab, "growth_model_preferred.csv"), stringsAsFactors = FALSE)

# --- endpoint color per colony (last non-NA reading) ---
last_non_na <- function(x) { y <- x[!is.na(x)]; if (length(y) == 0) NA else y[length(y)] }
endpoint <- series %>%
  group_by(plate_id, well_position) %>%
  arrange(hours_since_plate_start) %>%
  summarise(lab_l_late     = last_non_na(lab_l_mean),
            lab_a_late     = last_non_na(lab_a_mean),
            lab_b_late     = last_non_na(lab_b_mean),
            intensity_late = last_non_na(intensity_mean),
            area_late      = last_non_na(shape_area),
            .groups = "drop")

colony <- pref %>%
  select(plate_id, well_position, copper_mm, species, strain_code, trait, rate_area, rate_int) %>%
  tidyr::pivot_wider(names_from = trait,
                     values_from = c(rate_area, rate_int),
                     names_glue = "{.value}_{trait}", names_vary = "fastest")

d <- inner_join(colony, endpoint, by = c("plate_id", "well_position"))
d$species[is.na(d$species) | d$species == ""] <- "unknown"
note("Joined colonies w/ endpoint color: %d; dropped %d (no late image)",
     nrow(d), n_distinct(colony$plate_id, colony$well_position) - nrow(d))
stopifnot(nrow(d) == n_distinct(d$plate_id, d$well_position))

# Documented collinearity check (L* vs intensity at endpoint)
r_li <- cor(d$lab_l_late, d$intensity_late, use = "complete.obs")
note("Endpoint L* vs Intensity correlation r = %.3f (expected ~0.99)", r_li)

# Two growth modes: are they one axis or two?
r_modes <- cor(d$rate_area_ln_area, d$rate_int_intensity, use = "complete.obs")
note("Correlation of growth modes rate_area vs rate_int: r = %.3f", r_modes)

# --- Models: color/lightness ~ growth rate x copper (+ species/strain random) ---
fit <- function(fml, data, label) {
  m <- tryCatch(lmer(fml, data = data, REML = TRUE), error = function(e) e)
  if (inherits(m, "error")) { note("  MODEL FAILED [%s]: %s", label, conditionMessage(m)); return(NULL) }
  list(cf = coef(summary(m)), an = anova(m), model = m, label = label)
}

run_model <- function(outcome, rate_col, label) {
  fml <- as.formula(sprintf("%s ~ %s * factor(copper_mm) + (1 | species / strain_code)", outcome, rate_col))
  m <- fit(fml, d, label)
  if (is.null(m)) return(NULL)
  cf <- as.data.frame(m$cf); cf$term <- rownames(cf); cf$model <- label
  an <- as.data.frame(m$an);  an$term <- rownames(an);  an$model <- label
  list(cf = cf, an = an, model = m$model, label = label)
}

mods <- list(
  L_star_area  = run_model("lab_l_late", "rate_area_ln_area",    "L* ~ rate_area x Cu"),
  a_star_area  = run_model("lab_a_late", "rate_area_ln_area",    "a* ~ rate_area x Cu"),
  b_star_area  = run_model("lab_b_late", "rate_area_ln_area",    "b* ~ rate_area x Cu"),
  L_star_int   = run_model("lab_l_late", "rate_int_intensity",   "L* ~ rate_int x Cu"),
  a_star_int   = run_model("lab_a_late", "rate_int_intensity",   "a* ~ rate_int x Cu"),
  b_star_int   = run_model("lab_b_late", "rate_int_intensity",   "b* ~ rate_int x Cu")
)

cf_all <- bind_rows(lapply(mods, function(m) m$cf))
an_all <- bind_rows(lapply(mods, function(m) m$an))
write.csv(cf_all, file.path(out_tab, "color_growth_models.csv"), row.names = FALSE)
write.csv(an_all, file.path(out_tab, "color_growth_anova.csv"), row.names = FALSE)

sink(file.path(out_tab, "color_growth_models.txt"))
cat("=== Color / light-intensity as outcome of growth rate x copper ===\n")
cat(sprintf("Endpoint L* vs Intensity r = %.3f\n", r_li))
cat(sprintf("Growth modes rate_area vs rate_int r = %.3f\n", r_modes))
for (nm in names(mods)) {
  m <- mods[[nm]]
  if (is.null(m)) next
  cat("\n\n##### ", nm, ": ", m$label, " #####\n", sep = "")
  cat("-- anova --\n"); print(m$an)
  cat("\n-- fixed effects (interaction rows) --\n")
  ix <- grepl(":", rownames(m$cf))
  print(m$cf[ix, , drop = FALSE], digits = 3)
}
sink()

# --- Figures ---
theme_set(theme_bw(base_size = 11))
lvl <- sort(unique(d$copper_mm))
d$copper_mm <- factor(d$copper_mm, levels = lvl)

mk <- function(df, x, y, lab_x, lab_y, tt) {
  ggplot(df, aes(.data[[x]], .data[[y]])) +
    geom_point(aes(colour = copper_mm), size = 0.8, alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, colour = "grey30") +
    labs(x = lab_x, y = lab_y, title = tt, colour = "Cu (mM)") +
    theme(legend.position = "bottom")
}

ggsave(file.path(out_fig, "color_L_vs_ratearea_by_cu.png"),
       mk(d, "rate_area_ln_area", "lab_l_late", "peak d(log area)/dt (per h)", "endpoint L*",
          "L* (lightness) vs growth rate, by copper"), width = 7, height = 5.5, dpi = 150)
ggsave(file.path(out_fig, "color_a_vs_ratearea_by_cu.png"),
       mk(d, "rate_area_ln_area", "lab_a_late", "peak d(log area)/dt (per h)", "endpoint a*",
          "a* (red-green) vs growth rate, by copper"), width = 7, height = 5.5, dpi = 150)
ggsave(file.path(out_fig, "color_b_vs_ratearea_by_cu.png"),
       mk(d, "rate_area_ln_area", "lab_b_late", "peak d(log area)/dt (per h)", "endpoint b*",
          "b* (yellow-blue) vs growth rate, by copper"), width = 7, height = 5.5, dpi = 150)
ggsave(file.path(out_fig, "ratearea_vs_rateint.png"),
       mk(d, "rate_area_ln_area", "rate_int_intensity",
          "peak d(log area)/dt (per h)", "peak d(intensity)/dt (per h)",
          "Two growth modes: area expansion vs brightness rise, by copper"),
       width = 7, height = 5.5, dpi = 150)

writeLines(notes, file.path(out_tab, "color_growth_notes.txt"))
cat("\nDone: 03_color_interaction.R\n")
