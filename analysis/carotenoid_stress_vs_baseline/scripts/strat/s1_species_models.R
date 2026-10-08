#!/usr/bin/env Rscript
# Species (fixed) mixed models; strain, run and plate random. Reference species = R. mucilaginosa.
#  B : dose-0 wells:  a ~ species [+ lnA_c] + (1|strain) + (1|run) + (1|plate)          -> baseline differences
#  D : all wells:     a ~ dose_s*species + lnA_c + (1|strain) + (0+dose_s|strain) + (1|run) + (1|plate)  -> species-specific dose slope
#  P : same within R. mucilaginosa with population (pop1 = reference)
#  G : mucilaginosa vs all other species pooled (contrast only)
# dose_s = conc / max conc (0..1). Slope tables are per full dose range and per 0.1 of it. Zinc skipped (only mucilaginosa has >= 5 strains).
suppressMessages({library(lme4); library(lmerTest)})
R <- "analysis/carotenoid_stress_vs_baseline/results"; O <- file.path(R, "strat"); T <- "analysis/carotenoid_stress_vs_baseline/report/tables"
w  <- read.csv(file.path(R, "wells.csv"), stringsAsFactors = FALSE)
st <- read.csv(file.path(O, "strain_table.csv"), stringsAsFactors = FALSE)
w  <- merge(w, st[, c("strain_id", "species", "pop")], by = "strain_id", all.x = TRUE)
w  <- w[!is.na(w$species) & w$species != "Species Not Found", ]
REF <- "Rhodotorula mucilaginosa"; MINS <- 5
ctl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5), check.conv.singular = "ignore")
bh <- function(p) p.adjust(p, "BH")
out <- list(base = list(), slope = list(), pbase = list(), pslope = list(), pooled = list(), omni = list())
for (m in c("Chromium", "Copper", "Iron", "Lead")) {
  d <- w[w$Metal == m, ]
  ns <- tapply(d$strain_id, d$species, function(x) length(unique(x))); keep <- names(ns)[ns >= MINS]
  d <- d[d$species %in% keep, ]; d$species <- relevel(factor(d$species), REF)
  d$strain_id <- factor(d$strain_id); d$run_number <- factor(d$run_number); d$plate_id <- factor(d$plate_id)
  ns <- tapply(d$strain_id, d$species, function(x) length(unique(x)))
  cat(sprintf("== %s: %d wells, species (>= %d strains): %s\n", m, nrow(d), MINS, paste(names(ns), ns, sep = "=", collapse = "; ")))
  b0 <- d[d$conc == 0, ]
  for (adj in c(FALSE, TRUE)) {
    f <- if (adj) a ~ species + lnA_c + (1|strain_id) + (1|run_number) + (1|plate_id) else a ~ species + (1|strain_id) + (1|run_number) + (1|plate_id)
    fit <- lmer(f, b0, REML = TRUE, control = ctl); co <- summary(fit)$coefficients; an <- anova(fit)
    i <- grep("^species", rownames(co))
    out$base[[length(out$base) + 1]] <- data.frame(Metal = m, size_adjusted = adj, species = sub("^species", "", rownames(co)[i]), n_strains = ns[sub("^species", "", rownames(co)[i])],
        estimate = co[i, "Estimate"], se = co[i, "Std. Error"], p = co[i, "Pr(>|t|)"], singular = isSingular(fit))
    out$omni[[length(out$omni) + 1]] <- data.frame(Metal = m, model = ifelse(adj, "baseline_size_adjusted", "baseline"), term = "species", F = an["species", "F value"], p = an["species", "Pr(>F)"])
  }
  fd <- lmer(a ~ dose_s * species + lnA_c + (1|strain_id) + (0 + dose_s|strain_id) + (1|run_number) + (1|plate_id), d, REML = TRUE, control = ctl)
  co <- summary(fd)$coefficients; V <- as.matrix(vcov(fd)); an <- anova(fd); b <- fixef(fd)
  sp <- levels(d$species)
  for (s in sp) {
    L <- setNames(rep(0, length(b)), names(b)); L["dose_s"] <- 1; nm <- paste0("dose_s:species", s); if (s != REF) L[nm] <- 1
    est <- sum(L * b); se <- sqrt(as.numeric(t(L) %*% V %*% L))
    cp <- if (s == REF) NA else co[nm, "Pr(>|t|)"]
    out$slope[[length(out$slope) + 1]] <- data.frame(Metal = m, species = s, n_strains = ns[s], slope_full_range = est, se = se, z_p = 2 * pnorm(-abs(est / se)),
        slope_per_0.1 = est / 10, vs_reference_p = cp, singular = isSingular(fd))
  }
  out$omni[[length(out$omni) + 1]] <- data.frame(Metal = m, model = "dose_response_size_adjusted", term = "dose_s:species", F = an["dose_s:species", "F value"], p = an["dose_s:species", "Pr(>F)"])
  # pooled "other species" vs mucilaginosa (contrast only)
  d$grp <- factor(ifelse(d$species == REF, "mucilaginosa", "other_pooled"), levels = c("mucilaginosa", "other_pooled"))
  fb <- lmer(a ~ grp + lnA_c + (1|strain_id) + (1|run_number) + (1|plate_id), d[d$conc == 0, ], REML = TRUE, control = ctl)
  fs <- lmer(a ~ dose_s * grp + lnA_c + (1|strain_id) + (0 + dose_s|strain_id) + (1|run_number) + (1|plate_id), d, REML = TRUE, control = ctl)
  cb <- summary(fb)$coefficients["grpother_pooled", ]; cs <- summary(fs)$coefficients["dose_s:grpother_pooled", ]
  out$pooled[[length(out$pooled) + 1]] <- data.frame(Metal = m, n_other_strains = length(unique(d$strain_id[d$grp == "other_pooled"])), baseline_diff = cb["Estimate"], baseline_se = cb["Std. Error"], baseline_p = cb["Pr(>|t|)"],
        slope_diff = cs["Estimate"], slope_se = cs["Std. Error"], slope_p = cs["Pr(>|t|)"])
  # population within mucilaginosa
  p <- d[d$species == REF & !is.na(d$pop), ]; p$pop <- factor(paste0("pop", p$pop))
  npop <- tapply(p$strain_id, p$pop, function(x) length(unique(x)))
  fpb <- lmer(a ~ pop + lnA_c + (1|strain_id) + (1|run_number) + (1|plate_id), p[p$conc == 0, ], REML = TRUE, control = ctl)
  fps <- lmer(a ~ dose_s * pop + lnA_c + (1|strain_id) + (0 + dose_s|strain_id) + (1|run_number) + (1|plate_id), p, REML = TRUE, control = ctl)
  cb <- summary(fpb)$coefficients; i <- grep("^pop", rownames(cb))
  out$pbase[[length(out$pbase) + 1]] <- data.frame(Metal = m, population = sub("^pop", "pop", rownames(cb)[i]), n_strains = npop[sub("^pop", "", rownames(cb)[i])], estimate = cb[i, "Estimate"], se = cb[i, "Std. Error"], p = cb[i, "Pr(>|t|)"])
  out$omni[[length(out$omni) + 1]] <- data.frame(Metal = m, model = "population_baseline_size_adjusted", term = "pop", F = anova(fpb)["pop", "F value"], p = anova(fpb)["pop", "Pr(>F)"])
  out$omni[[length(out$omni) + 1]] <- data.frame(Metal = m, model = "population_dose_response_size_adjusted", term = "dose_s:pop", F = anova(fps)["dose_s:pop", "F value"], p = anova(fps)["dose_s:pop", "Pr(>F)"])
  cs <- summary(fps)$coefficients; V <- as.matrix(vcov(fps)); b <- fixef(fps)
  for (s in levels(p$pop)) {
    L <- setNames(rep(0, length(b)), names(b)); L["dose_s"] <- 1; nm <- paste0("dose_s:pop", s); if (nm %in% names(b)) L[nm] <- 1
    est <- sum(L * b); se <- sqrt(as.numeric(t(L) %*% V %*% L))
    out$pslope[[length(out$pslope) + 1]] <- data.frame(Metal = m, population = s, n_strains = npop[s], slope_full_range = est, se = se, z_p = 2 * pnorm(-abs(est / se)), vs_pop1_p = if (nm %in% rownames(cs)) cs[nm, "Pr(>|t|)"] else NA)
  }
  cat(sprintf("   done %s\n", m))
}
fin <- function(l, nm) { x <- do.call(rbind, l); rownames(x) <- NULL; write.csv(x, file.path(T, nm), row.names = FALSE); x }
b <- fin(out$base, "s1_species_baseline.csv"); b$p_BH <- NA
for (a in c(FALSE, TRUE)) { i <- b$size_adjusted == a; b$p_BH[i] <- bh(b$p[i]) }; write.csv(b, file.path(T, "s1_species_baseline.csv"), row.names = FALSE)
s <- fin(out$slope, "s1_species_dose_slope.csv"); fin(out$pbase, "s1_population_baseline.csv"); fin(out$pslope, "s1_population_dose_slope.csv")
fin(out$pooled, "s1_pooled_other_vs_mucilaginosa.csv"); o <- fin(out$omni, "s1_omnibus_tests.csv"); o$p_BH <- bh(o$p); write.csv(o, file.path(T, "s1_omnibus_tests.csv"), row.names = FALSE)
cat("done\n")
