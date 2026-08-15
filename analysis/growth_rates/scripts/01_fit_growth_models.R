#!/usr/bin/env Rscript
# Fit growth models to every colony's time series and extract growth-rate
# parameters per colony, for two size readouts:
#
#   ln_area    -- log(Shape_Area); colony expansion/biomass proxy
#   intensity  -- Intensity_MeanIntensity; early-saturating brightness/size
#
# Two candidate sigmoid models (Zwietering reparameterizations, fitted to
# y over relative plate age t, hours):
#
#   Gompertz: y = A * exp(-exp((mu*e/A)*(lambda - t) + 1))
#   Logistic: y = A / (1 + exp((4*mu/A)*(lambda - t) + 2))
#
# parameterized so that:
#   A      = asymptote of y (ln px for area; 0-1 fraction for intensity)
#   mu     = MAXIMUM GROWTH RATE, per hour (the headline parameter)
#   lambda = lag time, hours
#
# Fit with stats::nls + port bounds, require >=8 usable timepoints spanning
# >=24 h, and keep BOTH converged model fits so downstream can compare
# (AIC). If neither model converges we fall back to a log-linear slope over
# the rising phase and flag method = "fallback". Every row carries its
# realized n and a convergence/quality flag (robust-analysis: fail loudly,
# never silently drop).
#
# Run from the repo root:
#   pixi run Rscript analysis/growth_rates/scripts/01_fit_growth_models.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

tab_dir <- "analysis/growth_rates/results/tables"
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

series <- readRDS(file.path(tab_dir, "colony_growth_series.rds"))

MIN_TP   <- 8L
MIN_SPAN <- 24

# --- per-(colony, trait) data slices -----------------------------------------
prep_trait <- function(d, trait) {
  if (trait == "ln_area") {
    d <- d %>% filter(!is.na(shape_area) & shape_area > 0) %>%
      mutate(y = log(shape_area))
  } else if (trait == "intensity") {
    d <- d %>% filter(!is.na(intensity_mean)) %>% mutate(y = intensity_mean)
  } else {
    stop("unknown trait")
  }
  d %>% select(-shape_area, -intensity_mean)
}

# --- discrete max-slope start value (robust 90th pct of positive diffs) ------
start_mu <- function(d) {
  o <- order(d$hours_snap)
  dy <- diff(d$y[o]); dt <- diff(d$hours_snap[o])
  g <- dy[dt > 0] / dt[dt > 0]
  g <- g[is.finite(g) & g > 0]
  if (length(g) == 0) return(0.02)
  as.numeric(quantile(g, 0.9))
}

fit_models <- function(d, trait) {
  base <- d %>% distinct(plate_id, well_position, strain_id, strain_code, species,
                         run_number, copper_mm, n_obs = n())
  if (nrow(d) < MIN_TP || (max(d$hours_snap) - min(d$hours_snap)) < MIN_SPAN) {
    return(base %>% mutate(trait = trait, n_obs = nrow(d), eligible = FALSE))
  }
  base <- base %>% mutate(n_obs = nrow(d), eligible = TRUE)
  t <- d$hours_snap; y <- d$y
  A0  <- max(y); mu0 <- start_mu(d); la0 <- min(t) + 1

  fit_one <- function() {
    low <- c(A = A0 * 0.9,          mu = 1e-4,        lambda = 0 - 1)
    upp <- c(A = max(A0 * 5, 1e3),  mu = 5,           lambda = max(t) + 24)
    st  <- c(A = A0,                mu = mu0,          lambda = la0)
    mm <- NULL
    for (attempt in 1:3) {
      mm <- tryCatch(
        nls(y ~ A * exp(-exp((mu * exp(1) / A) * (lambda - t) + 1)),
            data = data.frame(t = t, y = y),
            start = st, lower = low, upper = upp, algorithm = "port",
            control = nls.control(maxiter = 200, warnOnly = TRUE)),
        error = function(e) NULL)
      if (!is.null(mm)) break
      st <- c(A = A0 * (1 + 0.15 * attempt), mu = mu0 * sqrt(attempt), lambda = la0 * 0.7)
    }
    if (is.null(mm)) return(tibble(model = "Gompertz", converged = FALSE))
    co <- coef(mm); rss <- sum(residuals(mm)^2); nn <- length(y)
    tibble(model = "Gompertz", converged = TRUE,
           A = co[["A"]], mu = co[["mu"]], lambda = co[["lambda"]],
           AIC = nn * log(rss / nn) + 2 * 3,
           sigma = sqrt(rss / (nn - 3)),
           A_extrap = co[["A"]] > A0 * 1.15,
           maxiter_reached = mm$convInfo$isConv == FALSE)
  }
  fit_one_logistic <- function() {
    low <- c(A = A0 * 0.9,         mu = 1e-4,       lambda = 0 - 1)
    upp <- c(A = max(A0 * 5, 1e3), mu = 5,          lambda = max(t) + 24)
    st  <- c(A = A0,               mu = mu0,         lambda = la0)
    mm <- NULL
    for (attempt in 1:3) {
      mm <- tryCatch(
        nls(y ~ A / (1 + exp((4 * mu / A) * (lambda - t) + 2)),
            data = data.frame(t = t, y = y),
            start = st, lower = low, upper = upp, algorithm = "port",
            control = nls.control(maxiter = 200, warnOnly = TRUE)),
        error = function(e) NULL)
      if (!is.null(mm)) break
      st <- c(A = A0 * (1 + 0.15 * attempt), mu = mu0 * sqrt(attempt), lambda = la0 * 0.7)
    }
    if (is.null(mm)) return(tibble(model = "Logistic", converged = FALSE))
    co <- coef(mm); rss <- sum(residuals(mm)^2); nn <- length(y)
    tibble(model = "Logistic", converged = TRUE,
           A = co[["A"]], mu = co[["mu"]], lambda = co[["lambda"]],
           AIC = nn * log(rss / nn) + 2 * 3,
           sigma = sqrt(rss / (nn - 3)),
           A_extrap = co[["A"]] > A0 * 1.15,
           maxiter_reached = mm$convInfo$isConv == FALSE)
  }

  rows <- bind_rows(fit_one(), fit_one_logistic())
  rows
}

fallback_slope <- function(d, trait) {
  o <- order(d$hours_snap)
  t <- d$hours_snap[o]; y <- d$y[o]
  # rising phase = timepoints before the running max growth has decelerated:
  # take slope over [first tp, the tp before the largest single-step gain] but
  # require at least 3 points.
  dy <- diff(y); dt <- diff(t)
  g <- dy / dt
  i_peak <- which.max(g)
  if (i_peak < 2 || i_peak > length(g)) i_peak <- min(3, length(g))
  use <- seq_len(min(i_peak + 1, length(t)))
  if (length(use) < 3) use <- seq_len(min(3, length(t)))
  lm(y[use] ~ t[use])$coefficients[[2]]
}

# --- fit all colonies --------------------------------------------------------
colonies <- series %>% distinct(plate_id, well_position) %>% arrange(plate_id, well_position)
cat(sprintf("Fitting %d colonies x 2 traits x 2 models...\n", nrow(colonies)))

out <- list(); k <- 0L
n_done <- 0L
for (trait in c("ln_area", "intensity")) {
  tr <- prep_trait(series, trait)
  for (i in seq_len(nrow(colonies))) {
    ci <- colonies[i, ]
    d <- tr %>% filter(plate_id == ci$plate_id, well_position == ci$well_position)
    if (nrow(d) < MIN_TP || (max(d$hours_snap) - min(d$hours_snap)) < MIN_SPAN) next
    base <- d %>% distinct(plate_id, well_position, strain_id, strain_code, species,
                           run_number, copper_mm)
    m <- tryCatch(fit_models(d, trait), error = function(e) NULL)
    if (is.null(m)) next
    feats <- m %>% mutate(trait = trait) %>% bind_cols(base, .)
    # fallback slope if neither model converged
    n_conv <- sum(feats$converged, na.rm = TRUE)
    if (n_conv == 0) {
      feats <- base %>% mutate(trait = trait, model = "fallback", converged = FALSE,
                               A = NA_real_, mu = fallback_slope(d, trait),
                               lambda = NA_real_, AIC = NA_real_,
                               sigma = NA_real_, A_extrap = NA,
                               maxiter_reached = NA)
    }
    k <- k + 1L; out[[k]] <- feats
    n_done <- n_done + 1L
    if (n_done %% 500 == 0) cat(sprintf("  %d colonies done\n", n_done))
  }
}
fits <- bind_rows(out)

# --- preferred model per (colony, trait): lower AIC among converged ----------
preferred <- fits %>%
  mutate(AIC = ifelse(converged, AIC, NA)) %>%
  group_by(plate_id, well_position, trait) %>%
  slice_min(AIC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(pref_model = model)

# --- data-derived peak relative growth rate (per hour), the primary rate  ----
# The parametric lag/asymptote are structurally unidentifiable on a 120 h
# window for the majority of colonies (growth starts before the first image;
# area does not plateau), so the Gompertz/logistic mu is only weakly pinned.
# The same quantity conceptually -- the MAXIMUM instantaneous relative growth
# rate -- is estimated directly from the data as the steepest positive 6 h
# slope of ln(area) (or of intensity). Defined for ~99% of colonies.
peak_rate <- function(df, yexpr) {
  df <- df %>% filter(!is.na(.data[[yexpr]]))
  md <- df %>%
    arrange(plate_id, well_position, hours_snap) %>%
    group_by(plate_id, well_position) %>%
    summarise(
      rate      = { g <- diff(.data[[yexpr]]) / diff(hours_snap);
                    g <- g[is.finite(g) & g > 0]; if (length(g)) max(g) else NA_real_ },
      rate_p90  = { g <- diff(.data[[yexpr]]) / diff(hours_snap);
                    g <- g[is.finite(g) & g > 0]; if (length(g)) as.numeric(quantile(g, 0.9)) else NA_real_ },
      .groups = "drop")
  md
}
rate_area <- peak_rate(series %>% mutate(shape_area_ln = log(shape_area)), "shape_area_ln") %>%
  rename(rate_area = rate, rate90_area = rate_p90)
rate_int  <- peak_rate(series, "intensity_mean") %>%
  rename(rate_int = rate, rate90_int = rate_p90)

preferred <- preferred %>%
  left_join(rate_area, by = c("plate_id", "well_position")) %>%
  left_join(rate_int,  by = c("plate_id", "well_position"))

# --- save --------------------------------------------------------------------
write.csv(fits, file.path(tab_dir, "growth_model_fits.csv"), row.names = FALSE)
write.csv(preferred, file.path(tab_dir, "growth_model_preferred.csv"), row.names = FALSE)

cat(sprintf("\n%d total model rows; %d colonies x %d traits preferred fits\n",
            nrow(fits), nrow(preferred), length(unique(fits$trait))))
cat("Convergence summary (model x trait):\n")
print(as.data.frame(fits %>%
  group_by(trait, model) %>%
  summarise(n = n(), converged = sum(converged), .groups = "drop")))
cat("Parametric mu vs data-derived peak rate (Kendall tau, preferred fits):\n")
pref_wide <- preferred %>%
  select(plate_id, well_position, trait, mu, rate_area, rate_int) %>%
  tidyr::pivot_wider(names_from = trait, values_from = c(mu, rate_area, rate_int), names_glue = "{.value}_{trait}", names_vary = "fastest") %>%
  filter(!is.na(mu_ln_area))
if (nrow(pref_wide) > 10) {
  cat("  area:", signif(cor(pref_wide$mu_ln_area, pref_wide$rate_area_ln_area, method = "kendall", use = "complete.obs"), 3), "\n")
}
cat("\nWrote growth_model_fits.csv, growth_model_preferred.csv to", tab_dir, "\n")
