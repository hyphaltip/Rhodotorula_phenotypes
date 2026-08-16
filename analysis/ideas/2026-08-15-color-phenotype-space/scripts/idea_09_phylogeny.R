#!/usr/bin/env Rscript
# Idea 09 (Phylogeneticist) — is any of the strain-level color/growth phenotype
# phylogenetically structured, or (consistent with idea 05) a continuum that
# spreads across clades?
#
# Tests Blomberg's K and Pagel's lambda (with permutation/LRT p-values) for
# strain-level traits with a tree + trait value:
#   - Cu reaction-norm slope     (slope_logchroma_per_mM; idea04)  sensitivity
#   - baseline chroma            (intercept_logchroma; idea04)      pigmentation
#   - colony size                (l10med_fixed; idea01b)            growth
#   - within-strain dispersion   (partial_slope_sd_cu; idea01b)     heterogeneity
#   - pigment pace               (median pace_loglog; idea08)
# Two scopes: all 277 matched tree strains, and within R. mucilaginosa only
# (the 69%-of-culture dominant lineage).
#
# Outputs (in results/):
#   idea09_phylo_signal.csv          per trait x scope: K, K p, lambda, lambda p
#   idea09_tip_traits.csv            tip-label -> strain -> trait matrix used
#   idea09_pruned_trees.rds          phylo objects for all + mucilaginosa subtrees
#   fig09_trait_on_tree_*.png        trait mapped on the pruned tree for each trait
suppressMessages({
  library(ape)
  library(phytools)
  library(ggtree)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

set.seed(7)
BASE <- here::here()
DATA_DIR <- file.path(BASE, "data/raw/rhodotorula-phyling-protein-tree")
RES <- file.path(BASE, "analysis/ideas/2026-08-15-color-phenotype-space/results")
META <- file.path(BASE, "analysis/ideas/2026-08-15-color-phenotype-space/data/strain_metadata.tsv")
N_PERM <- 1000

treefile <- file.path(DATA_DIR, "protein-Rhodotorula-taxa_278.fungi_odb10.fasttree.support.treefile")
tr <- read.tree(treefile)

# some pendant edges are 0 length (identical/duplicate seqs) -> singular covariance;
# give them a tiny positive length so phylosig's covariance solve is well-posed
z <- which(tr$edge.length == 0)
if (length(z)) {
  cat(sprintf("Jittering %d zero-length branches to %g\n", length(z), max(1e-6, min(tr$edge.length[tr$edge.length > 0])) * 1e-3))
  tr$edge.length[z] <- max(1e-6, min(tr$edge.length[tr$edge.length > 0])) * 1e-3
}
# FastTree yields extremely small branch lengths (down to ~5e-9); rescale the tree
# so mean root-to-tip = 1 to keep phylosig's covariance solve numerically stable
# (K is scale-invariant, lambda is unaffected by uniform scaling).
rescale_tree <- function(t) {
  d <- dist.nodes(t)
  dtt <- d[length(t$tip.label) + 1, 1:length(t$tip.label)]
  t$edge.length <- t$edge.length / mean(dtt)
  t
}
tr <- rescale_tree(tr)

## ---- 1. Normalize tip labels ------------------------------------------------
strip_suffix <- function(x) {
  x <- sub("\\.proteins\\.fa$", "", x)
  x <- sub("\\.proteins$", "", x)
  x
}
tiplab <- strip_suffix(tr$tip.label)

## ---- 2. Map tips -> strain_code ---------------------------------------------
meta <- read.delim(META, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
meta$strain_code <- as.character(meta$strain_code)
meta <- meta[!is.na(meta$strain_code) & nzchar(meta$strain_code), ]
# some physical strains occur under 2 strain_id rows (same code); dedupe to one
meta <- meta[order(!is.na(meta$species), -as.numeric(meta$strain_id)), ]
meta <- meta[!duplicated(meta$strain_code), ]
codes <- unique(meta$strain_code)

# match each tip to the LONGEST strain_code that is a suffix (preceded by _ or start)
tip_code <- vapply(tiplab, function(t) {
  hits <- codes[vapply(codes, function(c) grepl(paste0("_?", c, "$"), t), logical(1))]
  if (length(hits) == 0) return(NA_character_)
  hits[which.max(nchar(hits))]
}, character(1))

tips_df <- data.frame(tip_label = tr$tip.label, tiplab = tiplab,
                      strain_code = unname(tip_code), stringsAsFactors = FALSE)
tips_df <- merge(tips_df, meta[, c("strain_code", "strain_id", "species", "origin", "environment")],
                 by = "strain_code", all.x = TRUE)

cat(sprintf("Tips: %d matched strain_code, %d unmatched\n",
            sum(!is.na(tips_df$strain_code)), sum(is.na(tips_df$strain_code))))
cat("Unmatched tips:", paste(tips_df$tip_label[is.na(tips_df$strain_code)], collapse = "; "), "\n")

## ---- 3. Load per-strain traits ----------------------------------------------
rn  <- read.csv(file.path(RES, "idea04_reaction_norms.csv"), stringsAsFactors = FALSE)
inb <- read.csv(file.path(RES, "idea01b_within_strain.csv"), stringsAsFactors = FALSE)
ppc <- read.csv(file.path(RES, "idea08_pigment_pace.csv"), stringsAsFactors = FALSE)

ppc_med <- ppc %>% group_by(strain_id) %>% summarise(pace_loglog = median(pace_loglog, na.rm = TRUE))

traits <- rn[, c("strain_id", "species", "slope_logchroma_per_mM", "intercept_logchroma", "r2")]
traits <- merge(traits, inb[, c("strain_id", "l10med_fixed", "partial_slope_sd_cu")], by = "strain_id", all.x = TRUE)
traits <- merge(traits, ppc_med, by = "strain_id", all.x = TRUE)

## ---- 4. Align traits to tree tips --------------------------------------------
trait_cols <- c("slope_logchroma_per_mM", "intercept_logchroma", "l10med_fixed",
                "partial_slope_sd_cu", "pace_loglog")
tip <- tips_df[!is.na(tips_df$strain_id), ]
tip$strain_id <- as.numeric(tip$strain_id)
tip <- merge(tip, traits, by = "strain_id", all.x = TRUE, suffixes = c("", ".sp"))

cat(sprintf("Tips with strain_id in trait table: %d\n", nrow(tip)))

# build a named trait matrix keyed by tip_label
mk_mat <- function(td) {
  rownames(td) <- td$tip_label
  td[, trait_cols, drop = FALSE]
}
mat_all <- mk_mat(tip)

drop_tips <- setdiff(tr$tip.label, rownames(mat_all))
tr_pr <- drop.tip(tr, drop_tips)
mat_pr <- mat_all[tr_pr$tip.label, , drop = FALSE]

## ---- 5. Phylogenetic signal per trait ----------------------------------------
# Primary statistic: Mantel (Spearman) permutation test between phylogenetic
# (patristic) distance and absolute trait distance. Rank-based + permutation null ->
# robust to the near-comet shape of this tree (167/541 edges near-zero because of a
# huge polytomy of near-duplicate R. mucilaginosa genomes), where Blomberg's K is
# known to collapse (degenerate ~1e-7) and likelihood lambda is numerically unstable.
# Positive mantel r = close relatives are phenotypically more similar (signal).
mantel_test <- function(trok, x, nsim = 999, seed = 7) {
  x <- x[complete.cases(x)]
  D <- cophenetic(trok)[names(x), names(x)]
  dT <- as.matrix(dist(x))
  iu <- upper.tri(D)
  r.obs <- cor(D[iu], dT[iu], method = "spearman")
  set.seed(seed)
  rr <- numeric(nsim)
  for (k in seq_len(nsim)) {
    sh <- sample(length(x))
    Ds <- D[sh, sh]
    rr[k] <- cor(Ds[iu], dT[iu], method = "spearman")
  }
  list(r = r.obs, p = (1 + sum(rr >= r.obs)) / (nsim + 1))
}

signal_row <- function(name, trx, matx) {
  vals <- matx[, name, drop = TRUE]
  x <- vals; names(x) <- rownames(matx)
  x <- x[!is.na(x)]
  trok <- drop.tip(trx, setdiff(trx$tip.label, names(x)))
  x <- x[trok$tip.label]
  mt <- mantel_test(trok, x)
  # secondary (tree-shape-sensitive) estimates for reference:
  k  <- tryCatch(phylosig(trok, x, method = "K", test = FALSE), error = function(e) NA_real_)
  la <- tryCatch(phylosig(trok, x, method = "lambda", test = TRUE), error = function(e) NULL)
  data.frame(
    trait = name, n_tips = length(x),
    mantel_r_spearman = mt$r, mantel_p_perm = mt$p,
    blomberg_K_ref = as.numeric(k),
    pagel_lambda_ref = if (!is.null(la)) unname(la$lambda) else NA_real_,
    pagel_lambda_p_ref = if (!is.null(la)) unname(la$P) else NA_real_,
    stringsAsFactors = FALSE
  )
}

# --- all strains scope ---
res_all <- do.call(rbind, lapply(trait_cols, signal_row, trx = tr_pr, matx = mat_pr))
res_all$scope <- "all"

# --- within R. mucilaginosa scope ---
mu_ids <- tip$tip_label[tip$species == "Rhodotorula mucilaginosa"]
tr_mu <- drop.tip(tr_pr, setdiff(tr_pr$tip.label, mu_ids))
mat_mu <- mat_pr[tr_mu$tip.label, , drop = FALSE]
res_mu <- do.call(rbind, lapply(trait_cols, signal_row, trx = tr_mu, matx = mat_mu))
res_mu$scope <- "Rhodotorula mucilaginosa"

res <- rbind(res_all, res_mu)
res <- res[, c("scope", "trait", setdiff(names(res), c("scope", "trait")))]

cat("\n=== Phylogenetic signal (Mantel permutation primary) ===\n")
print(res, digits = 3)
cat("\nMantel r > 0 means close relatives are MORE phenotypically similar.\n")
cat("Blomberg K is ~1e-7 for all traits -> DEGENERATE (near-comet tree), not informative.\n")

## ---- 6. Outputs ---------------------------------------------------------------
write.csv(res, file.path(RES, "idea09_phylo_signal.csv"), row.names = FALSE)
tmp_tip <- tips_df
tmp_tip$strain_id <- as.numeric(tmp_tip$strain_id)
out_tip <- merge(tmp_tip, traits, by = "strain_id", all.x = TRUE, suffixes = c("", ".sp"))
write.csv(out_tip, file.path(RES, "idea09_tip_traits.csv"), row.names = FALSE)
saveRDS(list(all = tr_pr, mucilaginosa = tr_mu), file.path(RES, "idea09_pruned_trees.rds"))

## ---- 7. Figures: trait mapped on tree ------------------------------------------
plot_tree_trait <- function(trx, matx, col, title, fname) {
  vals <- matx[, col, drop = TRUE]
  x <- vals; names(x) <- rownames(matx)
  x <- x[!is.na(x)]
  trok <- drop.tip(trx, setdiff(trx$tip.label, names(x)))
  if (length(unique(x)) < 3) return(NULL)
  p <- ggtree(trok, layout = "rectangular") +
    geom_tippoint(aes(color = x), size = 1.6) +
    scale_color_gradient2(low = "#2c7fb8", mid = "#f7f7f7", high = "#d95f0e",
                          name = col) +
    ggtitle(title) +
    theme_tree2() + theme(plot.title = element_text(size = 10))
  ggsave(file.path(RES, fname), p, width = 9, height = 12, dpi = 150)
}

plot_tree_trait(tr_pr, mat_pr, "slope_logchroma_per_mM",
                "Idea 09 | Cu reaction-norm slope (log chroma / mM)", "fig09_trait_on_tree_cu_slope.png")
plot_tree_trait(tr_pr, mat_pr, "intercept_logchroma",
                "Idea 09 | Baseline chroma (log)", "fig09_trait_on_tree_baseline_chroma.png")
plot_tree_trait(tr_pr, mat_pr, "partial_slope_sd_cu",
                "Idea 09 | Within-strain dispersion widening", "fig09_trait_on_tree_dispersion.png")

cat("\nDone. Wrote idea09_phylo_signal.csv, idea09_tip_traits.csv, idea09_pruned_trees.rds, fig09_*.png\n")
