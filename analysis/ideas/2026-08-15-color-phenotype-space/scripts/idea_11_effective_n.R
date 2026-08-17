#!/usr/bin/env Rscript
# idea_11_effective_n.R
# Quantify effective number of independent genomes within R. mucilaginosa
# (the mapping population) from the PHYling protein tree.
#
# Rationale: GWAS power scales with effective independent haplotypes, not
# nominal strain count. Near-clonal genomes (patristic distance ~ 0) collapse
# into one effective unit and inflate LD, reducing power and raising FPR.
#
# Inputs:
#   - protein tree: BFD/results/phyling_pep/protein/tree/fungi_odb10/final_tree.nw
#   - phenotype tip set: idea09_tip_traits.csv
# Outputs:
#   - results/idea11_genome_redundancy.csv   (per-species tip/effective counts)
#   - results/idea11_effective_n.csv         (mapping-relevant effective n)

suppressMessages({
  library(ape)
})

HERE <- dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)
                                  [grep("--file=", commandArgs(trailingOnly = FALSE))])))
ROOT <- dirname(HERE)
TREE_PATH <- Sys.getenv(
  "PHYLING_TREE",
  "/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS/BFD/results/phyling_pep/protein/tree/fungi_odb10/final_tree.nw"
)
OUT <- file.path(ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("Reading protein tree:", TREE_PATH, "\n")
tr <- read.tree(TREE_PATH)
cat("Tips:", length(tr$tip.label), "\n")

# phenotype tip set
tips <- read.csv(file.path(ROOT, "results", "idea09_tip_traits.csv"))
tips <- tips[!is.na(tips$species.sp) & !is.na(tips$tip_label), ]
tip_to_species <- setNames(tips$species.sp, tips$tip_label)
tip_to_species <- tip_to_species[tr$tip.label]
tip_to_species <- tip_to_species[!is.na(tip_to_species)]

# distance matrix (patristic) for redundancy clustering
D <- cophenetic(tr)
# restrict to tips that have phenotype data for species with enough samples
species_counts <- table(tip_to_species)
species_counts <- species_counts[species_counts >= 3]

threshold <- 1e-7   # near-zero distance: same effective haplotype

res <- lapply(names(species_counts), function(sp) {
  lab <- names(tip_to_species)[tip_to_species == sp]
  if (length(lab) < 2) {
    return(data.frame(species = sp, n_tips = length(lab), n_distinct_haplotypes = length(lab),
                      redundancy_ratio = 1, n_independent_singletons = length(lab)))
  }
  Ds <- D[lab, lab]
  # single-linkage clustering at threshold
  h <- hclust(as.dist(Ds), method = "single")
  clusters <- cutree(h, h = threshold)
  n_eff <- length(unique(clusters))
  # count how many tips are in multi-tip (redundant) clusters vs singletons
  tab <- table(clusters)
  n_indep <- sum(tab == 1)
  data.frame(species = sp, n_tips = length(lab), n_distinct_haplotypes = n_eff,
             redundancy_ratio = length(lab) / n_eff,
             n_independent_singletons = n_indep)
})
red <- do.call(rbind, res)
red <- red[order(-red$n_tips), ]
rownames(red) <- NULL
cat("\n--- Genotype redundancy by species (threshold =", threshold, ") ---\n")
print(red)
write.csv(red, file.path(OUT, "idea11_genome_redundancy.csv"), row.names = FALSE)

# mapping-relevant effective n for R. mucilaginosa
muc <- red[red$species == "Rhodotorula mucilaginosa", ]
stopifnot(nrow(muc) == 1)
cat("\nR. mucilaginosa: nominal tips =", muc$n_tips,
    "| effective independent haplotypes =", muc$n_distinct_haplotypes, "\n")

eff_df <- data.frame(
  population = c("R. mucilaginosa (tips)", "R. mucilaginosa (effective)",
                 "All strains (tips)", "All strains (effective, est)"),
  sample_size = c(muc$n_tips, muc$n_distinct_haplotypes,
                  nrow(tips), round(nrow(tips) * min(1, muc$n_distinct_haplotypes / muc$n_tips), 0)),
  note = c("phenotype tips in protein tree", "independent haplotypes after collapsing near-clones",
           "idea09 mapping tip set", "same redundancy ratio applied to all-strain set")
)
write.csv(eff_df, file.path(OUT, "idea11_effective_n.csv"), row.names = FALSE)
cat("Wrote:", file.path(OUT, "idea11_effective_n.csv"), "\n")
