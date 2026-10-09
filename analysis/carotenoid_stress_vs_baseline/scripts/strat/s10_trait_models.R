#!/usr/bin/env Rscript
# Dose effect on colour and morphology traits. Each trait is z-scored per metal on the dose-0 wells (mean 0, SD 1), so effects are in
# SD units of the unstressed wells over the full dose range. Model: z ~ dose_s [+ lnA_c] + (1|strain) + (1|plate); plate nested in run.
# Run on the main wells (largest object >= 2000 px) and on all wells. Zinc skipped.
suppressMessages({library(lme4); library(lmerTest)})
R <- "analysis/carotenoid_stress_vs_baseline/results"; T <- "analysis/carotenoid_stress_vs_baseline/report/tables"
TR <- c("L", "a", "b", "chroma", "hue_deg", "sat", "val", "circ", "solid", "ecc", "compact", "extent", "aspect")
ctl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5), check.conv.singular = "ignore")
out <- list()
for (ds in c("wells_traits", "wells_traits_nofilter")) {
  w <- read.csv(file.path(R, paste0(ds, ".csv")), stringsAsFactors = FALSE); w <- w[w$Metal != "Zinc", ]
  for (m in c("Chromium", "Copper", "Lead")) {
    d <- w[w$Metal == m, ]; d$strain_id <- factor(d$strain_id); d$plate_id <- factor(d$plate_id)
    for (tr in TR) {
      d$z <- d[[tr]]; z0 <- d$z[d$conc == 0]; d$z <- (d$z - mean(z0, na.rm = TRUE)) / sd(z0, na.rm = TRUE); x <- d[!is.na(d$z), ]
      for (adj in c(FALSE, TRUE)) {
        f <- if (adj) z ~ dose_s + lnA_c + (1|strain_id) + (1|plate_id) else z ~ dose_s + (1|strain_id) + (1|plate_id)
        fit <- tryCatch(lmer(f, x, REML = TRUE, control = ctl), error = function(e) NULL); if (is.null(fit)) next
        co <- summary(fit)$coefficients["dose_s", ]
        out[[length(out) + 1]] <- data.frame(dataset = ifelse(ds == "wells_traits", "main: largest object >= 2000 px", "no area filter"), Metal = m, trait = tr, model = ifelse(adj, "at fixed size", "total"),
          n_wells = nrow(x), effect_sd = co["Estimate"], se = co["Std. Error"], p = co["Pr(>|t|)"])
      }
    }
    cat(ds, m, "done\n")
  }
}
x <- do.call(rbind, out); rownames(x) <- NULL; write.csv(x, file.path(T, "s10_trait_dose_effects.csv"), row.names = FALSE); cat("done\n")
