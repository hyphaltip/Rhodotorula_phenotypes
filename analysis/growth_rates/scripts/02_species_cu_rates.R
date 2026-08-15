#!/usr/bin/env Rscript
# 02_species_cu_rates.R
# Growth-rate analysis: how do per-strain growth rates (data-derived peak slopes)
# vary with copper concentration and species?
#
# Input : results/tables/growth_model_preferred.csv  (per-colony fits from 01)
# Output: results/tables/strain_x_cu_rates.csv       (strain x copper aggregate)
#         results/tables/rate_models_species_cu.csv  (mixed-model coefficients)
#         results/tables/rate_models_species_cu.txt  (human-readable summaries)
#         results/figures/rate_by_cu_{overall,spp}.png
#         results/tables/rate_fit_notes.txt          (validation log)
#
# Conventions: robust-analysis (asserts, row counts, sensitivity checks),
# Cu as factor (nonlinearity across levels), species = random effect (shrinkage),
# sensitivity: colony-level vs strain x Cu-level model agreement.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
})

# Run from the repo root (same convention as 01):
#   pixi run Rscript analysis/growth_rates/scripts/02_species_cu_rates.R
out_tab <- "analysis/growth_rates/results/tables"
out_fig <- "analysis/growth_rates/results/figures"
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

notes <- character()
note <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  notes <<- c(notes, msg)
}

pref <- read.csv(file.path(out_tab, "growth_model_preferred.csv"), stringsAsFactors = FALSE)
note("Loaded growth_model_preferred.csv: %d rows (colonies x traits)", nrow(pref))

stopifnot(all(c("copper_mm", "species", "strain_code", "rate_area", "rate_int", "trait") %in% names(pref)))

# --- wide per-colony table: one row per (plate_id, well_position) ---
colony <- pref %>%
  select(plate_id, well_position, strain_id, strain_code, species,
         run_number, copper_mm, pref_model, converged, trait,
         mu, rate_area, rate90_area, rate_int, rate90_int, sigma) %>%
  tidyr::pivot_wider(names_from = trait,
                     values_from = c(pref_model, converged, mu,
                                     rate_area, rate90_area, rate_int, rate90_int, sigma),
                     names_glue = "{.value}_{trait}", names_vary = "fastest")
note("Colony-level wide table: %d colonies", nrow(colony))
stopifnot(nrow(colony) == n_distinct(colony$plate_id, colony$well_position))

# Sanity: every colony maps to exactly one Cu / species / strain
stopifnot(all(!is.na(colony$copper_mm)))
colony$species[is.na(colony$species) | colony$species == ""] <- "unknown"
note("Species assignment: %d colonies unmatched -> 'unknown'", sum(colony$species == "unknown"))

# --- Primary rate readouts ---
# rate_area = peak positive 6 h slope of log(shape_area) per h
# rate_int  = peak positive 6 h slope of intensity_mean per h
for (t in c("ln_area", "intensity")) {
  colt_rate <- paste0("rate_area_", t)
  colori    <- paste0("rate_int_", t)
  na_rate <- sum(is.na(colony[[colt_rate]]))
  na_int  <- sum(is.na(colony[[colori]]))
  note("Trait %-10s members with rate_area: %d (%d NA), rate_int: %d (%d NA)",
       t, sum(!is.na(colony[[colt_rate]])), na_rate, sum(!is.na(colony[[colori]])), na_int)
}
colony_area <- colony %>% filter(!is.na(rate_area_ln_area))
colony_int  <- colony %>% filter(!is.na(rate_int_intensity))

# --- Strain x Cu aggregation (plate replicates = up to 4 runs per strain x Cu) ---
agg <- function(d, rate_col) {
  d %>%
    group_by(strain_id, strain_code, species, copper_mm) %>%
    summarise(n_repl   = sum(!is.na(.data[[rate_col]])),
              mean_rate = mean(.data[[rate_col]], na.rm = TRUE),
              sd_rate   = sd(.data[[rate_col]], na.rm = TRUE),
              se_rate   = sd_rate / sqrt(n_repl),
              .groups = "drop")
}
agg_area <- agg(colony_area, "rate_area_ln_area")
agg_int  <- agg(colony_int,  "rate_int_intensity")

agg_area <- agg_area %>% rename(n_repl_area = n_repl,
                                mean_rate_area = mean_rate, sd_rate_area = sd_rate, se_rate_area = se_rate)
agg_int  <- agg_int  %>% rename(n_repl_int = n_repl,
                                mean_rate_int  = mean_rate, sd_rate_int  = sd_rate,  se_rate_int  = se_rate)
sxc <- full_join(agg_area, agg_int, by = c("strain_id", "strain_code", "species", "copper_mm"))
sxc <- sxc %>% arrange(strain_code, copper_mm)
note("Strain x Cu aggregate: %d rows (%d strains x up to 7 Cu)", nrow(sxc), n_distinct(sxc$strain_id))
write.csv(sxc, file.path(out_tab, "strain_x_cu_rates.csv"), row.names = FALSE)

note("Replicate histogram (n_repl per strain x Cu):")
nhist <- table(sxc$n_repl_area)
note("  %s", paste(sprintf("n=%d:%d", as.integer(names(nhist)), as.integer(nhist)), collapse = "  "))
note("Strain x Cu rows with n_repl_area < 2: %d (single-culture estimates, flagged)",
     sum(sxc$n_repl_area < 2))

# --- Mixed models (primary: colony-level, Cu as factor, species + strain random) ---
fit_model <- function(f, data, label) {
  m <- tryCatch(lmer(f, data = data, REML = TRUE), error = function(e) e)
  if (inherits(m, "error")) {
    note("  MODEL FAILED [%s]: %s", label, conditionMessage(m))
    return(NULL)
  }
  a <- as.data.frame(anova(m))
  cf <- as.data.frame(coef(summary(m)))
  cf$term <- rownames(cf)
  cf$model <- label
  list(anova = a %>% mutate(term = rownames(a), model = label), cf = cf, model = m)
}

models <- list()
# area growth rate
models$area_all <- fit_model(
  rate_area_ln_area ~ factor(copper_mm) + (1 | species / strain_code),
  colony_area, "rate_area ~ Cu(factor) + species/strain (colony-level)")
# intensity rise rate
models$int_all <- fit_model(
  rate_int_intensity ~ factor(copper_mm) + (1 | species / strain_code),
  colony_int, "rate_int ~ Cu(factor) + species/strain (colony-level)")

# linear trend (numeric Cu) for direction/sensitivity
models$area_trend <- fit_model(
  rate_area_ln_area ~ copper_mm + (1 | species / strain_code),
  colony_area, "rate_area ~ Cu(numeric) + species/strain (sensitivity)")
models$int_trend <- fit_model(
  rate_int_intensity ~ copper_mm + (1 | species / strain_code),
  colony_int, "rate_int ~ Cu(numeric) + species/strain (sensitivity)")

# Interaction among well-sampled species (>= 8 strains), Cu x species fixed
well_spp <- sxc %>%
  filter(!species %in% c("unknown", "R. sp. clade I")) %>%
  count(species)  # strain-Cu rows; better count distinct strains:
well_spp <- colony_area %>%
  filter(!species %in% c("unknown", "R. sp. clade I")) %>%
  group_by(species) %>% summarise(n_strains = n_distinct(strain_id), .groups = "drop") %>%
  filter(n_strains >= 8)
note("Well-sampled species (>= 8 strains): %s",
     paste(well_spp$species, "(", well_spp$n_strains, ")", collapse = ", "))
well_spp_areas <- colony_area %>% filter(species %in% well_spp$species)
well_spp_ints  <- colony_int  %>% filter(species %in% well_spp$species)
well_spp_areas$species <- factor(well_spp_areas$species)
well_spp_ints$species  <- factor(well_spp_ints$species)
models$area_xspp <- fit_model(
  rate_area_ln_area ~ factor(copper_mm) * species + (1 | strain_code),
  well_spp_areas, "rate_area ~ Cu x species (well-sampled only)")
models$int_xspp <- fit_model(
  rate_int_intensity ~ factor(copper_mm) * species + (1 | strain_code),
  well_spp_ints, "rate_int ~ Cu x species (well-sampled only)")

# --- Collect and write model outputs ---
cf_all <- bind_rows(lapply(models, function(m) if (!is.null(m)) m$cf))
an_all <- bind_rows(lapply(models, function(m) if (!is.null(m)) m$anova))
write.csv(cf_all, file.path(out_tab, "rate_models_species_cu.csv"), row.names = FALSE)

sink(file.path(out_tab, "rate_models_species_cu.txt"))
cat("=== Growth-rate models: species x copper ===\n")
for (nm in names(models)) {
  m <- models[[nm]]
  if (is.null(m)) next
  cat("\n\n##### ", nm, ": ", attr(m, "label"), " #####\n", sep = "")
  cat("\n-- anova --\n"); print(m$anova)
  cat("\n-- fixed effects --\n"); print(m$cf[!grepl("Intercept|strain|species", rownames(m$cf), ignore.case = TRUE), , drop = FALSE], digits = 3)
}
sink()

# Marginal means vs control (0 mM) with 95% CI via lmerTest difflsmeans
# Contrasts vs control (0 mM = reference level of factor(copper_mm)):
# the factor coefficients ARE differences vs 0 mM; add Wald 95% CI.
write_diff_vs_control <- function(m, label, file) {
  if (is.null(m)) return(invisible(NULL))
  f <- fixef(m); V <- vcov(m)
  met <- grepl("factor\\(copper_mm\\)", names(f)) | grepl("^factor\\(copper_mm\\)", names(f))
  cn <- names(f)[met]
  if (length(cn) == 0) return(invisible(NULL))
  seb <- sqrt(diag(V)[met])
  d <- data.frame(copper_mm = sub("factor\\(copper_mm\\)", "", cn),
                  diff = f[met], se = seb,
                  lo = f[met] - 1.96 * seb, hi = f[met] + 1.96 * seb,
                  tvalue = f[met] / seb,
                  pvalue = 2 * pnorm(-abs(f[met] / seb)))
  d$model <- label
  rownames(d) <- NULL
  write.csv(d, file, row.names = FALSE)
}
write_diff_vs_control(models$area_all$model, "rate_area", file.path(out_tab, "rate_area_diff_0mM.csv"))
write_diff_vs_control(models$int_all$model,  "rate_int",  file.path(out_tab, "rate_int_diff_0mM.csv"))

# --- Figures ---
theme_set(theme_bw(base_size = 11))
cu_order <- sort(unique(sxc$copper_mm))
sxc$copper_mm <- factor(sxc$copper_mm, levels = cu_order)

# overall: strain-Cu means with species-colored points
p_overall <- ggplot(sxc, aes(copper_mm, mean_rate_area)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_jitter(aes(colour = species), width = 0.12, size = 0.9, alpha = 0.55) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 1.1, colour = "black") +
  stat_summary(aes(group = 1), fun = mean, geom = "point", size = 2.3, colour = "black") +
  labs(x = "Copper (mM)", y = "peak d(log area)/dt per h",
       title = "Growth rate (area) vs copper, strain x Cu means",
       colour = "Species") +
  theme(legend.position = "bottom")
ggsave(file.path(out_fig, "rate_by_cu_overall.png"), p_overall,
       width = 7.5, height = 5.5, dpi = 150)

# per well-sampled species + others
sxc$grp <- ifelse(sxc$species %in% well_spp$species, sxc$species, "rare / other")
sxc$grp <- factor(sxc$grp, levels = c(well_spp$species, "rare / other"))
p_spp <- ggplot(sxc, aes(copper_mm, mean_rate_area)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_jitter(width = 0.10, size = 0.7, alpha = 0.45, colour = "grey40") +
  stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 1.0, colour = "darkred") +
  stat_summary(aes(group = 1), fun = mean, geom = "point", size = 1.8, colour = "darkred") +
  facet_wrap(~ grp, scales = "free_y") +
  labs(x = "Copper (mM)", y = "peak d(log area)/dt per h",
       title = "Growth rate (area) by species and copper") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_fig, "rate_by_cu_spp.png"), p_spp,
       width = 9, height = 6.5, dpi = 150)

# intensity-rate version (well-sampled spp)
sxc$grp2 <- ifelse(sxc$species %in% well_spp$species, sxc$species, "rare / other")
sxc$grp2 <- factor(sxc$grp2, levels = c(well_spp$species, "rare / other"))
p_int <- ggplot(sxc, aes(copper_mm, mean_rate_int)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_jitter(width = 0.10, size = 0.7, alpha = 0.45, colour = "grey40") +
  stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 1.0, colour = "darkblue") +
  stat_summary(aes(group = 1), fun = mean, geom = "point", size = 1.8, colour = "darkblue") +
  facet_wrap(~ grp2, scales = "free_y") +
  labs(x = "Copper (mM)", y = "peak d(intensity)/dt per h",
       title = "Consumer-light-intensity rise rate by species and copper") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_fig, "rate_int_by_cu_spp.png"), p_int,
       width = 9, height = 6.5, dpi = 150)

sink(file.path(out_tab, "rate_fit_notes.txt"))
cat(paste(notes, collapse = "\n"), "\n")
sink()

cat("\nDone: 02_species_cu_rates.R\n")
