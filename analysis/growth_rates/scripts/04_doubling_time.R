#!/usr/bin/env Rscript
# 04_doubling_time.R
# Test the doubling-time question: is the copper trend robust to a
# phase-independent growth-rate estimator?
#
# The headline rate from 01/02 (rate_area = max 6 h slope of log(area)) is
# phase/length sensitive: it can reflect "still growing when observed" rather
# than "growing faster". This script computes per-colony phase-independent
# descriptors and asks whether the monotone-Cu pattern persists, flips, or
# vanishes:
#
#   k_region  -- slope of log(area)~t over the best-fitting (highest R^2)
#                sliding 4-point window = exponential-region specific rate
#   dbl_region-- ln(2)/k_region  (exponential doubling time, h)
#   k_max6    -- max 6 h (2-point) slope of log(area)  [same as 01's rate_area]
#   t50       -- interpolated hours to reach 50% of in-window max area
#   saturated -- end area >= 0.97 * in-window max area (plateaued in-window)
#   end_ratio -- end area / in-window max area
#
# Same descriptors for the bounded intensity trait are computed with the
# caveat that log-linear growth breaks near saturation (kept for reference).
#
# Run from repo root:
#   pixi run Rscript analysis/growth_rates/scripts/04_doubling_time.R

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
note <- function(...) { msg <- sprintf(...); cat(msg, "\n"); notes <<- c(notes, msg) }

series <- readRDS(file.path(out_tab, "colony_growth_series.rds"))
note("Series: %d rows, %d colonies (%d trays %s)",
     nrow(series), n_distinct(series$plate_id, series$well_position),
     n_distinct(series$plate_id), names(table(series$copper_mm))[1] %||% "x")

# --- per-colony growth descriptors -----------------------------------------
W2 <- 2  # ~6 h window -> max slope
W4 <- 4  # ~18 h window -> exponential region (best linear fit in log space)

desc_colony <- function(plate_id, well_position, h, area, intensity) {
  one <- function(y, kind) {
    n <- length(y)
    if (n < max(W2, W4) || any(is.na(y)) || min(y, na.rm = TRUE) <= 0) {
      return(list(k_max = NA_real_, k_region = NA_real_, dbl_region = NA_real_,
                  r2_region = NA_real_))
    }
    ly <- log(y)
    # max 6 h (2-point) slope
    k_max <- NA_real_
    if (n >= W2 + 1) {
      k_max <- max(-Inf, sapply(seq_len(n - W2), function(i)
        (ly[i + W2] - ly[i]) / (h[i + W2] - h[i])))
    }
    # best 4-point window in log space
    if (n >= W4) {
      gr <- sapply(seq_len(n - W4 + 1), function(i) {
        j <- i:(i + W4 - 1)
        m <- stats::lm(ly[j] ~ h[j])
        c(k = unname(coef(m)[2]), r2 = summary(m)$r.squared)
      })
      best <- which.max(gr["r2", ])
      r2b <- gr["r2", best]
      kb <- gr["k", best]
      if (!is.na(kb) && r2b >= 0.97) {
        list(k_max = k_max, k_region = kb, dbl_region = log(2) / kb, r2_region = r2b)
      } else {
        list(k_max = k_max, k_region = NA_real_, dbl_region = NA_real_, r2_region = r2b)
      }
    } else {
      list(k_max = k_max, k_region = NA_real_, dbl_region = NA_real_, r2_region = NA_real_)
    }
  }

  a <- one(area, "area")
  i <- one(intensity, "intensity")

  # t50 (area): first crossing of 0.5 * in-window max, linear interpolation
  t50 <- NA_real_
  mmax <- max(area, na.rm = TRUE)
  thr <- 0.5 * mmax
  idx <- which(area >= thr)
  if (length(idx) > 0 && idx[1] > 1) {
    p <- idx[1]
    h0 <- h[p - 1]; a0 <- area[p - 1]; a1 <- area[p]; h1 <- h[p]
    t50 <- h0 + (thr - a0) * (h1 - h0) / (a1 - a0)
  } else if (length(idx) > 0) {
    t50 <- h[idx[1]]
  }
  saturated <- !is.na(mmax) && last_non_na(area) >= 0.97 * mmax
  end_ratio <- last_non_na(area) / mmax

  tibble(plate_id = plate_id, well_position = well_position,
         n_tp = length(area), k_max6_area = a$k_max, k_region_area = a$k_region,
         dbl_region_area = a$dbl_region, r2_region_area = a$r2_region,
         k_max6_int = i$k_max, k_region_int = i$k_region, dbl_region_int = i$dbl_region,
         r2_region_int = i$r2_region,
         t50_h = t50, saturated = saturated, end_ratio = end_ratio,
         max_area = mmax)
}

last_non_na <- function(x) { y <- x[!is.na(x)]; if (length(y) == 0) NA else y[length(y)] }
`%||%` <- function(a, b) if (length(a)) a else b

d <- series %>%
  filter(shape_area > 0, !is.na(shape_area), !is.na(intensity_mean)) %>%
  group_by(plate_id, well_position) %>%
  arrange(hours_since_plate_start) %>%
  do(desc_colony(.$plate_id[1], .$well_position[1], .$hours_since_plate_start,
                 .$shape_area, .$intensity_mean)) %>%
  ungroup()

# join colony metadata
meta <- series %>%
  group_by(plate_id, well_position) %>%
  summarise(species = last(species), strain_id = last(strain_id),
            strain_code = last(strain_code), copper_mm = last(copper_mm),
            run_number = last(run_number), .groups = "drop")
d <- left_join(d, meta, by = c("plate_id", "well_position"))
d$species[is.na(d$species) | d$species == ""] <- "unknown"
note("Colony descriptors: %d colonies (n_tp >= 4 required)", nrow(d))
write.csv(d, file.path(out_tab, "colony_doubling.csv"), row.names = FALSE)

# coverage of each estimator
note("rate_area defined: %d, k_region_area defined: %d, t50 defined: %d, saturated: %d",
     sum(!is.na(d$k_max6_area)), sum(!is.na(d$k_region_area)),
     sum(!is.na(d$t50_h)), sum(d$saturated))

# --- does a phase-independent doubling time change the Cu conclusion? ------
# relationship between the two estimators
tau <- cor(d$k_max6_area, d$dbl_region_area, method = "kendall", use = "complete.obs")
note("Kendall tau (rate_area vs dbl_region): %.3f (expect strongly negative if same axis)", tau)

# mixed model on log(doubling time) - phase-independent rate
byn <- d %>% filter(!is.na(dbl_region_area), dbl_region_area > 1)
note("log-doubling-time model on %d colonies", nrow(byn))
m_dbl <- tryCatch(
  lmer(log(dbl_region_area) ~ factor(copper_mm) + (1 | species / strain_code),
       data = byn, REML = TRUE),
  error = function(e) { note("  dbl model failed: %s", conditionMessage(e)); NULL })
if (!is.null(m_dbl)) {
  an <- anova(m_dbl); cat("-- log(dbl) ~ Cu(factor)+species/strain --\n"); print(an)
  cf <- as.data.frame(coef(summary(m_dbl)))
  cf$term <- rownames(cf)
  write.csv(cf, file.path(out_tab, "dbl_model_cu.csv"), row.names = FALSE)
  # back-transform: fold-change in doubling time vs 0 mM at 30 mM
  f30 <- exp(cf["factor(copper_mm)30", "Estimate"])
  note("Fold-change in doubling time at 30 mM vs 0 mM: %.2fx", f30)
}

# saturation fraction by Cu: extent-based view
sat_cu <- d %>%
  mutate(copper_mm = factor(copper_mm, levels = sort(unique(copper_mm)))) %>%
  group_by(copper_mm) %>%
  summarise(n = n(), saturated_pct = 100 * mean(saturated, na.rm = TRUE),
            median_t50 = median(t50_h, na.rm = TRUE),
            median_dbl = median(dbl_region_area, na.rm = TRUE), .groups = "drop")
note("Saturation / t50 / doubling time by Cu:")
print(sat_cu, digits = 3)
write.csv(sat_cu, file.path(out_tab, "dbl_by_cu.csv"), row.names = FALSE)

# --- figures ---------------------------------------------------------------
theme_set(theme_bw(base_size = 11))
p_dbl <- d %>%
  filter(!is.na(dbl_region_area)) %>%
  mutate(copper_mm = factor(copper_mm, levels = sort(unique(copper_mm)))) %>%
  ggplot(aes(copper_mm, dbl_region_area)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  geom_boxplot(outlier.size = 0.4, width = 0.6, fill = "grey90") +
  stat_summary(aes(group = 1), fun = median, geom = "line", colour = "darkred", linewidth = 1) +
  labs(x = "Copper (mM)", y = "exponential doubling time (h)",
       title = "Phase-independent doubling time vs copper (shorter = faster growth)") +
  scale_y_log10()
ggsave(file.path(out_fig, "dbl_region_by_cu.png"), p_dbl, width = 7, height = 5, dpi = 150)

p_sat <- sat_cu %>% ggplot(aes(copper_mm, saturated_pct)) +
  geom_col(fill = "steelblue", alpha = 0.85) +
  geom_text(aes(label = sprintf("%.0f%%", saturated_pct)), vjust = -0.4, size = 3) +
  labs(x = "Copper (mM)", y = "% colonies saturated in-window",
       title = "Fraction of colonies reaching plateau within 120 h window") +
  ylim(0, 105)
ggsave(file.path(out_fig, "saturation_by_cu.png"), p_sat, width = 7, height = 5, dpi = 150)

writeLines(notes, file.path(out_tab, "dbl_fit_notes.txt"))
cat("\nDone: 04_doubling_time.R\n")
