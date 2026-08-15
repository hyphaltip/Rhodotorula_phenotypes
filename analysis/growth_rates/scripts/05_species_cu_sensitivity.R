#!/usr/bin/env Rscript
# 05_species_cu_sensitivity.R
# Which strains / species are most and least sensitive to copper?
#
# Rationale (from 04_doubling_time.R): the exponential doubling time is flat
# across 0-30 mM (~37 h, fold-change 0.99), so growth-RATE is not the readout
# of Cu sensitivity. The Cu phenotype is EXTENT-limited: saturation fraction
# falls 95% -> 74% and colonies stay in expansion through the window. The
# sensitivity index is therefore defined on final extent/yield and on the
# plateau-reachability (missingness-as-phenotype) axis, with run-rate and
# doubling-time fold-changes reported for completeness:
#
#   extent_ratio = median max(area) @ high Cu (25,30 pooled)
#                  / median max(area) @ low Cu (0,5 pooled)      (<1 = sensitive)
#   log2_extent  = log2(extent_ratio)
#   sat_drop     = saturation% @ low Cu - saturation% @ high Cu  (>0 = sensitive)
#   dbl_fold     = median dbl @ high / median dbl @ low          (~1 = rate-neutral)
#
# Outputs:
#   strain_sensitivity.csv        (per-strain metrics; all strains)
#   species_sensitivity.csv       (per-species metrics; ranked)
#   species_extent_model.txt/.csv (well-sampled species: log(max_area) ~ Cu x species)
#   results/figures/sensitivity_*.png
#
# Run from repo root: pixi run Rscript analysis/growth_rates/scripts/05_species_cu_sensitivity.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
})

out_tab <- "analysis/growth_rates/results/tables"
out_fig <- "analysis/growth_rates/results/figures"

d <- read.csv(file.path(out_tab, "colony_doubling.csv"), stringsAsFactors = FALSE)
d$species[is.na(d$species) | d$species == ""] <- "unknown"
d$copper_mm <- as.numeric(d$copper_mm)
d$grp_cu <- ifelse(d$copper_mm >= 25, "high", ifelse(d$copper_mm <= 5, "low", "mid"))

cat("Loaded colony_doubling.csv: ", nrow(d), " colonies\n", sep = "")

# ---- per-strain and per-species metrics ------------------------------------
strain_sens <- d %>%
  filter(!is.na(max_area), grp_cu %in% c("low", "high")) %>%
  group_by(strain_id, strain_code, species) %>%
  summarise(n_low = sum(grp_cu == "low"), n_high = sum(grp_cu == "high"),
            max_a_low = median(max_area[grp_cu == "low"], na.rm = TRUE),
            max_a_high = median(max_area[grp_cu == "high"], na.rm = TRUE),
            dbl_low = median(dbl_region_area[grp_cu == "low"], na.rm = TRUE),
            dbl_high = median(dbl_region_area[grp_cu == "high"], na.rm = TRUE),
            sat_low = 100 * mean(saturated[grp_cu == "low"], na.rm = TRUE),
            sat_high = 100 * mean(saturated[grp_cu == "high"], na.rm = TRUE),
            .groups = "drop") %>%
  mutate(extent_ratio = max_a_high / max_a_low,
         log2_extent = log2(extent_ratio),
         sat_drop = sat_low - sat_high,
         dbl_fold = dbl_high / dbl_low,
         n_group = pmin(n_low, n_high))
write.csv(strain_sens, file.path(out_tab, "strain_sensitivity.csv"), row.names = FALSE)

species_sens <- strain_sens %>%
  group_by(species) %>%
  summarise(n_strains = n_distinct(strain_id),
            n_colonies_low = sum(n_low, na.rm = TRUE),
            n_colonies_high = sum(n_high, na.rm = TRUE),
            median_extent_ratio = median(extent_ratio, na.rm = TRUE),
            median_log2_extent = median(log2_extent, na.rm = TRUE),
            median_sat_drop = median(sat_drop, na.rm = TRUE),
            median_dbl_fold = median(dbl_fold, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(min_group = pmin(n_colonies_low, n_colonies_high)) %>%
  arrange(median_log2_extent)
write.csv(species_sens, file.path(out_tab, "species_sensitivity.csv"), row.names = FALSE)
cat("\n-- Species ranked by log2 extent sensitivity (most sensitive first; <1 = extent shrinks) --\n")
print(species_sens %>% select(species, n_strains, min_group, median_extent_ratio,
                              median_log2_extent, median_sat_drop, median_dbl_fold), digits = 3)

# ---- extent mixed model: does the Cu slope differ among well-sampled species? ----
well <- d %>%
  filter(species %in% c("Rhodotorula mucilaginosa", "Rhodotorula paludigena",
                        "Rhodotorula diobovata", "Rhodotorula toruloides",
                        "Rhodotorula dairenensis", "Rhodotorula sphaerocarpa") &
           species != "unknown" & !is.na(max_area))
well <- well %>% mutate(species = factor(species), Cu = copper_mm)
cat("\nExtent model (log max_area ~ Cu x species, well-sampled): ", nrow(well), " colonies\n", sep = "")

m_ext <- tryCatch(
  lmer(log(max_area) ~ Cu * species + (1 | strain_code), data = well, REML = TRUE),
  error = function(e) { cat("  model failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(m_ext)) {
  an <- as.data.frame(anova(m_ext)); an$term <- rownames(an)
  write.csv(an, file.path(out_tab, "species_extent_anova.csv"), row.names = FALSE)
  sink(file.path(out_tab, "species_extent_model.txt"))
  print(an)
  cat("\n-- fixed effects (Cu term + Cu:species interactions vs mucilaginosa ref) --\n")
  cf <- as.data.frame(coef(summary(m_ext))); cf$term <- rownames(cf)
  ix <- grepl("Cu", rownames(cf))
  print(cf[ix, ], digits = 3)
  sink()
}

# ---- figures ---------------------------------------------------------------
theme_set(theme_bw(base_size = 11))
rank10 <- species_sens %>% filter(!species %in% c("unknown", "Rhodotorula sp. clade I")) %>%
  top_n(12, n_strains)

p_rank <- rank10 %>%
  mutate(species = reorder(species, median_log2_extent)) %>%
  ggplot(aes(median_log2_extent, species)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(size = 3, colour = "darkred") +
  labs(x = "log2(extent at 25-30 mM / extent at 0-5 mM)  [<0 = Cu-sensitive]",
       y = "", title = "Species Cu sensitivity by final-extent ratio") +
  scale_x_continuous(labels = function(x) sprintf("%.2f", x))
ggsave(file.path(out_fig, "sensitivity_species_rank.png"), p_rank, width = 6.5, height = 5, dpi = 150)

# per well-sampled species, extent by Cu
spp_med <- d %>%
  filter(species %in% c("Rhodotorula mucilaginosa", "Rhodotorula paludigena",
                        "Rhodotorula diobovata", "Rhodotorula toruloides",
                        "Rhodotorula dairenensis", "Rhodotorula sphaerocarpa")) %>%
  mutate(species = factor(species),
         copper_mm = factor(copper_mm, levels = sort(unique(copper_mm)))) %>%
  group_by(species, copper_mm) %>%
  summarise(m = median(max_area, na.rm = TRUE), .groups = "drop")

p_ext <- spp_med %>% ggplot(aes(copper_mm, m, group = 1)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ species, scales = "free_y") +
  scale_y_continuous(labels = scales::scientific) +
  labs(x = "Copper (mM)", y = "median max colony area (px)",
       title = "Final extent by copper per species (sensitivity of yield)")
ggsave(file.path(out_fig, "sensitivity_extent_by_cu_spp.png"), p_ext, width = 9, height = 6, dpi = 150)

p_sat <- d %>%
  filter(species %in% c("Rhodotorula mucilaginosa", "Rhodotorula paludigena",
                        "Rhodotorula diobovata", "Rhodotorula toruloides",
                        "Rhodotorula dairenensis", "Rhodotorula sphaerocarpa")) %>%
  mutate(species = factor(species),
         copper_mm = factor(copper_mm, levels = sort(unique(copper_mm)))) %>%
  group_by(species, copper_mm) %>%
  summarise(sat = 100 * mean(saturated, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(copper_mm, sat, group = 1)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ species) +
  labs(x = "Copper (mM)", y = "% colonies saturated in-window",
       title = "Plateau reachability by copper per species") +
  scale_y_continuous(limits = c(0, 100))
ggsave(file.path(out_fig, "sensitivity_saturation_by_cu_spp.png"), p_sat, width = 9, height = 6, dpi = 150)

cat("\nDone: 05_species_cu_sensitivity.R\n")
