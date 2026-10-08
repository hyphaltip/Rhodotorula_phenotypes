#!/usr/bin/env Rscript
# s3: dose regimes (sub-inhibitory vs inhibitory by colony-size ratio) ; s4: size-matched dose effect ; s6: batch (run) heterogeneity + variance components.
# All models: strain, run, plate random intercepts unless stated. Zinc skipped (one plate per dose, runs differ in strain set and dose range).
suppressMessages({library(lme4); library(lmerTest)})
R <- "analysis/carotenoid_stress_vs_baseline/results"; T <- "analysis/carotenoid_stress_vs_baseline/report/tables"
w <- read.csv(file.path(R, "wells.csv"), stringsAsFactors = FALSE)
ctl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5), check.conv.singular = "ignore")
RATIO_CUT <- 0.5   # regime threshold: median over strains of (colony area at dose / area at dose 0)
res3 <- list(); res3b <- list(); ratio <- list(); res4 <- list(); res6 <- list(); res6b <- list()
for (m in c("Chromium", "Copper", "Iron", "Lead")) {
  d <- w[w$Metal == m, ]; d$strain_id <- factor(d$strain_id); d$run_number <- factor(d$run_number); d$plate_id <- factor(d$plate_id)
  # ---- s3 size ratio per dose
  sa <- aggregate(lnA ~ strain_id + conc, d, mean); s0 <- sa[sa$conc == 0, c("strain_id", "lnA")]; names(s0)[2] <- "lnA0"
  sa <- merge(sa, s0); sa$ratio <- exp(sa$lnA - sa$lnA0)
  rr <- aggregate(ratio ~ conc, sa[sa$conc > 0, ], function(x) c(med = median(x), q25 = quantile(x, .25), q75 = quantile(x, .75), n = length(x)))
  rr <- data.frame(Metal = m, conc = rr$conc, median_ratio = rr$ratio[, 1], q25 = rr$ratio[, 2], q75 = rr$ratio[, 3], n_strains = rr$ratio[, 4])
  rr$regime <- ifelse(rr$median_ratio >= RATIO_CUT, "sub-inhibitory", "inhibitory"); ratio[[m]] <- rr
  cat(sprintf("== %s regimes: %s\n", m, paste(sprintf("%g:%s(%.2f)", rr$conc, substr(rr$regime, 1, 3), rr$median_ratio), collapse = " ")))
  for (reg in unique(rr$regime)) {
    doses <- c(0, rr$conc[rr$regime == reg]); x <- d[d$conc %in% doses, ]
    if (length(unique(x$conc)) < 3) { cat("   skip", reg, "(fewer than 3 doses)\n"); next }
    for (adj in c(FALSE, TRUE)) {
      f <- if (adj) a ~ dose_s + lnA_c + (1|strain_id) + (1|run_number) + (1|plate_id) else a ~ dose_s + (1|strain_id) + (1|run_number) + (1|plate_id)
      fit <- lmer(f, x, REML = TRUE, control = ctl); co <- summary(fit)$coefficients["dose_s", ]
      res3[[length(res3) + 1]] <- data.frame(Metal = m, regime = reg, doses = paste(doses, collapse = ","), n_wells = nrow(x), size_adjusted = adj,
        slope_per_0.1_of_max_dose = co["Estimate"] / 10, se = co["Std. Error"] / 10, p = co["Pr(>|t|)"])
    }
  }
  # dose-factor effects (vs 0), with and without size
  d$cf <- relevel(factor(d$conc), ref = "0")
  for (adj in c(FALSE, TRUE)) {
    f <- if (adj) a ~ cf + lnA_c + (1|strain_id) + (1|run_number) + (1|plate_id) else a ~ cf + (1|strain_id) + (1|run_number) + (1|plate_id)
    fit <- lmer(f, d, REML = TRUE, control = ctl); co <- summary(fit)$coefficients; i <- grep("^cf", rownames(co))
    res3b[[length(res3b) + 1]] <- data.frame(Metal = m, conc = as.numeric(sub("^cf", "", rownames(co)[i])), size_adjusted = adj, effect_vs_0 = co[i, "Estimate"], se = co[i, "Std. Error"])
  }
  # ---- s4 size-matched: bins on pooled ln area quantiles
  q <- quantile(d$lnA, seq(0, 1, .2)); q[1] <- q[1] - 1e-6; d$sbin <- cut(d$lnA, q, labels = paste0("S", 1:5))
  for (b in levels(d$sbin)) {
    x <- d[d$sbin == b, ]; nd <- length(unique(x$conc))
    row <- data.frame(Metal = m, size_bin = b, ln_area_lo = unname(q[as.integer(sub("S", "", b))]), ln_area_hi = unname(q[as.integer(sub("S", "", b)) + 1]),
      n_wells = nrow(x), n_doses = nd, doses = paste(sort(unique(x$conc)), collapse = ","), mean_a = mean(x$a), sd_a = sd(x$a), slope_full_range = NA, se = NA, p = NA)
    if (nd >= 3 && nrow(x) >= 150) {
      fit <- lmer(a ~ dose_s + (1|strain_id) + (1|run_number) + (1|plate_id), x, REML = TRUE, control = ctl); co <- summary(fit)$coefficients["dose_s", ]
      row$slope_full_range <- co["Estimate"]; row$se <- co["Std. Error"]; row$p <- co["Pr(>|t|)"]
    }
    res4[[length(res4) + 1]] <- row
  }
  # ---- s6 per-run dose effect at fixed size + heterogeneity; variance components
  for (r in levels(d$run_number)) {
    x <- d[d$run_number == r, ]
    if (length(unique(x$conc)) < 5) next
    fit <- lmer(a ~ dose_s + lnA_c + (1|strain_id) + (1|plate_id), x, REML = TRUE, control = ctl); co <- summary(fit)$coefficients["dose_s", ]
    res6[[length(res6) + 1]] <- data.frame(Metal = m, run = r, n_wells = nrow(x), n_strains = length(unique(x$strain_id)), n_doses = length(unique(x$conc)),
      dose_effect = co["Estimate"], se = co["Std. Error"], p = co["Pr(>|t|)"])
  }
  fit <- lmer(a ~ dose_s + lnA_c + (1|strain_id) + (1|run_number) + (1|plate_id), d, REML = TRUE, control = ctl); vc <- as.data.frame(VarCorr(fit))
  v <- setNames(vc$vcov, vc$grp); tot <- sum(v)
  res6b[[length(res6b) + 1]] <- data.frame(Metal = m, n_runs = nlevels(d$run_number), var_strain = v["strain_id"], var_run = v["run_number"], var_plate = v["plate_id"], var_resid = v["Residual"],
      share_strain = v["strain_id"] / tot, share_run = v["run_number"] / tot, share_plate = v["plate_id"] / tot, share_resid = v["Residual"] / tot, singular = isSingular(fit))
  cat("   done", m, "\n")
}
bind <- function(l) { x <- do.call(rbind, l); rownames(x) <- NULL; x }
write.csv(bind(ratio), file.path(T, "s3_size_ratio_by_dose.csv"), row.names = FALSE)
write.csv(bind(res3), file.path(T, "s3_regime_slopes.csv"), row.names = FALSE)
write.csv(bind(res3b), file.path(T, "s3_dose_effects.csv"), row.names = FALSE)
write.csv(bind(res4), file.path(T, "s4_size_matched.csv"), row.names = FALSE)
r6 <- bind(res6); r6$weight <- 1 / r6$se^2
het <- do.call(rbind, lapply(split(r6, r6$Metal), function(z) { mu <- sum(z$weight * z$dose_effect) / sum(z$weight); Q <- sum(z$weight * (z$dose_effect - mu)^2); k <- nrow(z)
  data.frame(Metal = z$Metal[1], k_runs = k, pooled_effect = mu, Q = Q, Q_p = pchisq(Q, k - 1, lower.tail = FALSE), I2 = max(0, (Q - (k - 1)) / Q)) }))
write.csv(r6[, setdiff(names(r6), "weight")], file.path(T, "s6_per_run_dose_effect.csv"), row.names = FALSE); write.csv(het, file.path(T, "s6_run_heterogeneity.csv"), row.names = FALSE)
write.csv(bind(res6b), file.path(T, "s6_variance_components.csv"), row.names = FALSE)
cat("done\n")
