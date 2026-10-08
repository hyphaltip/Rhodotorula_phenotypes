#!/usr/bin/env Rscript
# Mixed models per metal (Zinc included; treated as ONE experiment, see below) on the well table (replicates = wells on different plates within a run).
#  M0: a ~ dose_s                 + (1 + dose_s | strain) + (1 | plate)   total dose effect on a*
#  M1: a ~ dose_s + lnA_c         + (1 + dose_s | strain) + (1 | plate)   dose effect at the same colony size
# dose_s = conc / max conc (0..1) so slopes are "change in a* from no metal to the top dose".
# strain intercept variance = inherent (baseline) differences; strain dose-slope variance = strain-specific induction.
suppressMessages({library(lme4); library(lmerTest)})
args <- commandArgs(trailingOnly = TRUE)
inp <- if (length(args) >= 1) args[1] else "analysis/carotenoid_stress_vs_baseline/results/wells.csv"
out <- if (length(args) >= 2) args[2] else "analysis/carotenoid_stress_vs_baseline/results"
w <- read.csv(inp, stringsAsFactors = FALSE)
if (Sys.getenv("SKIP_ZINC") == "1") w <- w[w$Metal != "Zinc", ]   # Zinc random-slope fits are slow and uninformative (1 plate per dose)
w$strain_id <- factor(w$strain_id); w$plate_id <- factor(w$plate_id)
res <- list(); blups <- list()
# Zinc: user states runs d000388 and d000390 are one experiment. No run term. Two disjoint strain sets cover doses 0-15 and 10-30
# (1 well per strain x dose, 9 plates), so the dose slope is learned across strain sets and per-strain slopes are weakly identified.
for (m in unique(w$Metal)) {
  d <- w[w$Metal == m, ]
  cat(sprintf("== %s: %d wells, %d strains, %d plates\n", m, nrow(d), nlevels(droplevels(d$strain_id)), nlevels(droplevels(d$plate_id))))
  for (mod in c("M0", "M1")) {
    f <- if (mod == "M0") a ~ dose_s + (1 + dose_s | strain_id) + (1 | plate_id)
         else             a ~ dose_s + lnA_c + (1 + dose_s | strain_id) + (1 | plate_id)
    fit <- tryCatch(lmer(f, data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
                    error = function(e) { cat("  FAILED", mod, conditionMessage(e), "\n"); NULL })
    if (is.null(fit)) next
    vc <- as.data.frame(VarCorr(fit)); co <- summary(fit)$coefficients
    g <- function(grp, v1, v2 = NA) { r <- vc[vc$grp == grp & vc$var1 == v1 & (is.na(v2) & is.na(vc$var2) | !is.na(v2) & vc$var2 %in% v2), "vcov"]; if (length(r)) r[1] else NA }
    cv <- vc[vc$grp == "strain_id" & vc$var1 == "(Intercept)" & !is.na(vc$var2), "vcov"]
    sI <- g("strain_id", "(Intercept)"); sS <- g("strain_id", "dose_s"); sP <- g("plate_id", "(Intercept)"); sR <- vc[vc$grp == "Residual", "vcov"]
    cr <- attr(VarCorr(fit)$strain_id, "correlation")[1, 2]
    one <- data.frame(Metal = m, model = mod, n_wells = nrow(d), n_strains = nlevels(droplevels(d$strain_id)),
      dose_effect = co["dose_s", "Estimate"], dose_se = co["dose_s", "Std. Error"], dose_p = co["dose_s", "Pr(>|t|)"],
      size_effect_per_lnA = if (mod == "M1") co["lnA_c", "Estimate"] else NA,
      var_strain_intercept = sI, var_strain_dose_slope = sS, cor_intercept_slope = cr, var_plate = sP, var_resid = sR,
      cov_intercept_slope = if (length(cv)) cv[1] else NA,
      repeatability_dose0 = sI / (sI + sP + sR),
      repeatability_dose1 = (sI + 2 * (if (length(cv)) cv[1] else 0) + sS) / ((sI + 2 * (if (length(cv)) cv[1] else 0) + sS) + sP + sR),
      singular = isSingular(fit))
    res[[length(res) + 1]] <- one
    {
      b <- ranef(fit)$strain_id
      blups[[paste(m, mod)]] <- data.frame(Metal = m, model = mod, strain_id = rownames(b), blup_intercept = b[["(Intercept)"]], blup_dose_slope = b[["dose_s"]],
                               fixed_intercept = fixef(fit)[["(Intercept)"]], fixed_dose = fixef(fit)[["dose_s"]])
    }
    cat(sprintf("  %s dose effect %.2f (SE %.2f, p=%.2g); var strain-int %.2f, strain-slope %.2f, plate %.2f, resid %.2f; singular=%s\n",
        mod, one$dose_effect, one$dose_se, one$dose_p, sI, sS, sP, sR, one$singular))
  }
}
# M2: dose as a factor (no linearity assumption); per-dose effects at fixed size
fx <- list()
for (m in unique(w$Metal)) {
  d <- w[w$Metal == m, ]; d$cf <- relevel(factor(d$conc), ref = "0")
  f2 <- tryCatch(lmer(a ~ cf + lnA_c + (1 | strain_id) + (1 | plate_id), data = d, REML = TRUE), error = function(e) NULL)
  if (!is.null(f2)) { co <- summary(f2)$coefficients; i <- grep("^cf", rownames(co))
    fx[[m]] <- data.frame(Metal = m, conc = as.numeric(sub("^cf", "", rownames(co)[i])), effect_at_fixed_size = co[i, "Estimate"], se = co[i, "Std. Error"]) }
}
if (length(fx)) write.csv(do.call(rbind, fx), file.path(out, "dose_factor_effects_M2.csv"), row.names = FALSE)
stopifnot(length(res) > 0)
write.csv(do.call(rbind, res), file.path(out, "mixed_model_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, blups), file.path(out, "strain_blups.csv"), row.names = FALSE)
cat("done\n")
