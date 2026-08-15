#!/usr/bin/env Rscript
# Diagnostic and summary plots for the plate-position analysis. Figures are
# written as PNG to results/figures/ and used directly by the Rmd report.
#
# Run from the repo root:
#   pixi run Rscript analysis/explore_plate_position/scripts/03_plots.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

tab_dir <- "analysis/explore_plate_position/results/tables"
fig_dir <- "analysis/explore_plate_position/results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

endpoint <- readRDS(file.path(tab_dir, "endpoint_colonies.rds"))
vc_table <- read.csv(file.path(tab_dir, "variance_components_by_copper.csv"))

ggsave_ <- function(name, plot, width = 7, height = 5) {
  ggsave(file.path(fig_dir, name), plot, width = width, height = height, dpi = 150)
}

# 1. Variance-components stacked bars: strain vs run vs plate vs residual, per
# trait, one panel per Cu concentration (analysis is stratified by Cu).
p1 <- vc_table %>%
  mutate(group = factor(group, levels = c("strain_code", "run_number", "plate_id", "Residual"),
                         labels = c("Strain", "Run (batch)", "Plate", "Residual")),
         copper_mm = factor(copper_mm, levels = sort(unique(copper_mm)),
                            labels = paste0(sort(unique(copper_mm)), " mM Cu"))) %>%
  ggplot(aes(x = trait, y = pct_of_total, fill = group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ copper_mm, ncol = 4) +
  labs(x = NULL, y = "% of total variance", fill = NULL,
       title = "Variance partition by Cu concentration (stratified)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave_("variance_components_by_copper.png", p1, width = 13, height = 9)

# 1b. Plate-explained variance (%) heatmap across trait x Cu concentration.
p1b <- vc_table %>%
  filter(group == "plate_id") %>%
  mutate(trait_short = sub(" \\(.*", "", trait),
         copper_mm = factor(copper_mm)) %>%
  ggplot(aes(x = copper_mm, y = trait_short, fill = pct_of_total)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", pct_of_total)), color = "grey15", size = 3) +
  scale_fill_viridis_c(option = "C") +
  labs(x = "Copper concentration (mM)", y = NULL,
       title = "Plate-explained variance (%) by Cu concentration",
       fill = "% of total variance") +
  theme_minimal(base_size = 11)
ggsave_("plate_pct_heatmap.png", p1b, width = 8, height = 5)

# 2. Plate heatmap of mean L*/a*/b* residual-from-strain-mean, one representative plate.
# Center each trait on its own strain mean so plate structure isn't swamped by
# strain-to-strain differences.
plate_example <- endpoint %>%
  group_by(strain_code) %>%
  mutate(across(c(lab_l_mean, lab_a_mean, lab_b_mean, shape_area),
                ~ .x - mean(.x), .names = "{.col}_resid")) %>%
  ungroup()

one_plate <- plate_example %>% count(plate_id, sort = TRUE) %>% slice(1) %>% pull(plate_id)

p2 <- plate_example %>%
  filter(plate_id == one_plate) %>%
  ggplot(aes(x = grid_col, y = grid_row, fill = lab_l_mean_resid)) +
  geom_tile(color = "grey90") +
  scale_y_reverse() +
  scale_fill_gradient2(midpoint = 0, low = "#2166ac", mid = "white", high = "#b2182b") +
  coord_fixed() +
  labs(title = sprintf("L* (strain-centered) across plate %s", one_plate),
       x = "grid_col", y = "grid_row", fill = "L* - strain mean") +
  theme_minimal(base_size = 12)
ggsave_("plate_heatmap_L_example.png", p2, width = 8, height = 5)

# 3. Strain-centered trait vs edge distance, all plates pooled.
p3 <- plate_example %>%
  select(edge_dist, lab_l_mean_resid, lab_a_mean_resid, lab_b_mean_resid, shape_area_resid) %>%
  pivot_longer(-edge_dist, names_to = "trait", values_to = "resid") %>%
  mutate(trait = recode(trait,
                         lab_l_mean_resid = "L*", lab_a_mean_resid = "a*",
                         lab_b_mean_resid = "b*", shape_area_resid = "Colony area")) %>%
  ggplot(aes(x = factor(edge_dist), y = resid)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ trait, scales = "free_y") +
  labs(x = "Distance from plate edge (0 = border well)", y = "Trait - strain mean",
       title = "Strain-centered trait value by distance from plate edge") +
  theme_minimal(base_size = 12)
ggsave_("edge_distance_boxplots.png", p3, width = 9, height = 6)

# 4. Neighbor-mean scatter for L* (the trait with the strongest adjacency signal).
adj_dat <- readRDS(file.path(tab_dir, "adjacency_dataset.rds"))
p4 <- adj_dat %>%
  ggplot(aes(x = lab_l_mean_neighbor_mean, y = lab_l_mean)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_smooth(method = "lm", color = "#b2182b") +
  labs(x = "Mean L* of grid-adjacent neighbors (same plate)", y = "Focal colony L*",
       title = "Focal colony L* vs neighboring-colony L* (same plate, 4-connected)") +
  theme_minimal(base_size = 12)
ggsave_("adjacency_scatter_L.png", p4)

cat("Wrote figures to", fig_dir, "\n")
