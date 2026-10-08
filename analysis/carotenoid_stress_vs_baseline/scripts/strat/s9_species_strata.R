#!/usr/bin/env Rscript
# Stratified dose models: R. mucilaginosa alone, every other species with >= 5 strains alone, and the other species pooled.
# Same model in every stratum (so strata are comparable): a ~ dose_s [+ lnA_c] + (1|strain) + (1|plate), plate nested in run.
# Small strata (5-18 strains) cannot support random dose slopes, so strain intercept and plate only. dose_s = conc / max conc.
suppressMessages({library(lme4); library(lmerTest)})
R <- "analysis/carotenoid_stress_vs_baseline/results"; T <- "analysis/carotenoid_stress_vs_baseline/report/tables"
w  <- read.csv(file.path(R, "wells.csv"), stringsAsFactors = FALSE)
st <- read.csv(file.path(R, "strat/strain_table.csv"), stringsAsFactors = FALSE)
w  <- merge(w, st[, c("strain_id", "species", "pop")], by = "strain_id", all.x = TRUE)
w  <- w[!is.na(w$species) & w$species != "Species Not Found" & w$Metal != "Zinc", ]
REF <- "Rhodotorula mucilaginosa"; MINS <- 5
ctl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5), check.conv.singular = "ignore")
mod <- list(); dfc <- list(); cnt <- list()
for (m in c("Chromium", "Copper", "Iron", "Lead")) {
  d <- w[w$Metal == m, ]; ns <- tapply(d$strain_id, d$species, function(x) length(unique(x)))
  keep <- names(ns)[ns >= MINS]; strata <- c(REF, setdiff(keep, REF))
  d$stratum <- ifelse(d$species %in% keep, d$species, NA)
  sets <- lapply(strata, function(s) d[d$species == s, ]); names(sets) <- strata
  sets[["Other species (>= 5 strains, pooled)"]] <- d[d$species %in% setdiff(keep, REF), ]
  for (s in names(sets)) {
    x <- sets[[s]]; x$strain_id <- factor(x$strain_id); x$plate_id <- factor(x$plate_id)
    cat(sprintf("%s | %s: %d strains, %d wells, %d doses\n", m, s, nlevels(x$strain_id), nrow(x), length(unique(x$conc))))
    cnt[[length(cnt) + 1]] <- data.frame(Metal = m, stratum = s, n_strains = nlevels(x$strain_id), n_wells = nrow(x), n_doses = length(unique(x$conc)),
      mean_a_dose0 = mean(x$a[x$conc == 0]), mean_lnA_dose0 = mean(x$lnA[x$conc == 0]))
    f <- a ~ dose_s + (1|strain_id) + (1|plate_id)   # total effect (no size term); a size slope per small stratum is not identifiable (dose and size confounded)
    fit <- tryCatch(lmer(f, x, REML = TRUE, control = ctl), error = function(e) NULL)
    if (!is.null(fit)) { co <- summary(fit)$coefficients["dose_s", ]
      mod[[length(mod) + 1]] <- data.frame(Metal = m, stratum = s, model = "total", n_strains = nlevels(x$strain_id), n_wells = nrow(x), dose_effect = co["Estimate"], se = co["Std. Error"], p = co["Pr(>|t|)"], singular = isSingular(fit)) }
  }
}

# Per-species dose-factor effects at ONE shared size slope (same size treatment as Table 5): a ~ cf*species + lnA_c + (1|strain) + (1|plate)
for (m in c("Chromium", "Copper", "Iron", "Lead")) {
  d <- w[w$Metal == m, ]; ns <- tapply(d$strain_id, d$species, function(x) length(unique(x))); keep <- names(ns)[ns >= MINS]
  d <- d[d$species %in% keep, ]; d$species <- relevel(factor(d$species), REF); d$cf <- relevel(factor(d$conc), ref = "0"); d$strain_id <- factor(d$strain_id); d$plate_id <- factor(d$plate_id)
  fit <- tryCatch(lmer(a ~ cf * species + lnA_c + (1|strain_id) + (1|plate_id), d, REML = TRUE, control = ctl), error = function(e) NULL); if (is.null(fit)) next
  b0 <- fixef(fit); V <- as.matrix(vcov(fit)); nm <- names(b0)
  for (s in levels(d$species)) for (dz in setdiff(levels(d$cf), "0")) {
    L <- setNames(rep(0, length(b0)), nm); k1 <- paste0("cf", dz); if (!(k1 %in% nm)) next; L[k1] <- 1
    if (s != REF) { k2 <- paste0("cf", dz, ":species", s); if (k2 %in% nm) L[k2] <- 1 else next }
    dfc[[length(dfc) + 1]] <- data.frame(Metal = m, stratum = s, conc = as.numeric(dz), effect_at_fixed_size = sum(L * b0), se = sqrt(as.numeric(t(L) %*% V %*% L)))
  }
  cat(m, "shared-slope dose-factor model done\n")
}
b <- function(l) { x <- do.call(rbind, l); rownames(x) <- NULL; x }
write.csv(b(mod), file.path(T, "s9_stratum_models.csv"), row.names = FALSE)
write.csv(b(dfc), file.path(T, "s9_stratum_dose_factor.csv"), row.names = FALSE)
write.csv(b(cnt), file.path(T, "s9_stratum_counts.csv"), row.names = FALSE)
cat("done\n")
