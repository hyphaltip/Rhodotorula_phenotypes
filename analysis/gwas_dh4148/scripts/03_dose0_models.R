#!/usr/bin/env Rscript
# Dose-0 (YPD control) consistency: variance components of a*, ln area and growth rate across strains, experiments (metal screens), runs and plates.
suppressMessages({library(lme4)})
R <- "analysis/gwas_dh4148/results"; T <- "analysis/gwas_dh4148/report/tables"
w <- read.csv(file.path(R, "dose0_wells.csv.gz"), stringsAsFactors = FALSE)
st <- read.csv(file.path(R, "strain_traits.csv"), stringsAsFactors = FALSE)[, c("strain_id", "gwas_panel")]; st <- unique(st)
w$strain_id <- factor(w$strain_id); w$Metal <- factor(w$Metal); w$plate_id <- factor(w$plate_id)
w$run <- factor(paste(w$Metal, w$run_number)); w$lnR <- ifelse(is.na(w$rgr), NA, log(pmax(w$rgr, 1e-4)))
ctl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5), check.conv.singular = "ignore")
traits <- c(a = "a", lnA = "lnA", rgr = "rgr")
out <- list()
fit_vc <- function(d, tr, label) {
  d <- d[!is.na(d[[tr]]), ]
  f <- if (nlevels(droplevels(d$Metal)) > 1) as.formula(paste(tr, "~ 1 + (1|strain_id) + (1|Metal) + (1|run) + (1|plate_id)")) else as.formula(paste(tr, "~ 1 + (1|strain_id) + (1|run) + (1|plate_id)"))
  m <- lmer(f, d, REML = TRUE, control = ctl); v <- as.data.frame(VarCorr(m))
  vv <- setNames(v$vcov, v$grp); tot <- sum(vv)
  data.frame(scope = label, trait = tr, n_wells = nrow(d), n_strains = nlevels(droplevels(d$strain_id)), grand_mean = fixef(m)[[1]],
             strain = vv["strain_id"], metal_screen = if ("Metal" %in% names(vv)) vv["Metal"] else NA, run = vv["run"], plate = vv["plate_id"], residual = vv["Residual"],
             share_strain = vv["strain_id"] / tot, share_metal = if ("Metal" %in% names(vv)) vv["Metal"] / tot else NA, share_run = vv["run"] / tot, share_plate = vv["plate_id"] / tot, share_resid = vv["Residual"] / tot,
             singular = isSingular(m), row.names = NULL)
}
main3 <- w[w$Metal %in% c("Chromium", "Copper", "Lead"), ]
for (tr in names(traits)) {
  out[[length(out) + 1]] <- fit_vc(main3, tr, "Cr+Cu+Pb pooled")
  out[[length(out) + 1]] <- fit_vc(w, tr, "all five metals pooled (Fe and Zn partial plates)")
  for (m in c("Chromium", "Copper", "Lead")) out[[length(out) + 1]] <- fit_vc(w[w$Metal == m, ], tr, m)
  out[[length(out) + 1]] <- fit_vc(main3[as.character(main3$strain_id) %in% as.character(st$strain_id[which(as.character(st$gwas_panel) %in% c("True","TRUE"))]), ], tr, "Cr+Cu+Pb pooled, GWAS panel strains")
}
res <- do.call(rbind, out); write.csv(res, file.path(T, "dose0_variance_components.csv"), row.names = FALSE)
# strain best linear unbiased predictions for the pooled model (strain effect shared across screens)
bl <- list()
for (tr in names(traits)) {
  d <- main3[!is.na(main3[[tr]]), ]
  m <- lmer(as.formula(paste(tr, "~ 1 + (1|strain_id) + (1|Metal) + (1|run) + (1|plate_id)")), d, REML = TRUE, control = ctl)
  re <- ranef(m, condVar = TRUE)$strain_id; pv <- attr(re, "postVar")
  bl[[tr]] <- data.frame(strain_id = rownames(re), trait = tr, blup = re[, 1], se = sqrt(as.numeric(pv)), n_wells = as.vector(table(d$strain_id)[rownames(re)]))
}
write.csv(do.call(rbind, bl), file.path(T, "dose0_strain_blups.csv"), row.names = FALSE)
cat("done\n")
