# GWAS, dose-0 consistency, growth rate and IC50 on the DH4148 reference

Generated 2026-10-09. Phenotypes: Chromium, Copper and Lead screens (Iron and Zinc are held out, D-56). Genotypes: popgen callset on the DH4148 assembly (GCA_058775505.1). Panel: 126 haploid *R. mucilaginosa* strains (D-49). All numbers come from `report/tables/*.csv`. Scripts: `scripts/01` to `scripts/12`.

## Summary

1. **The 126-strain panel is 14 clonal lineages.** Strains that differ by at most 2,000 SNPs form 14 groups (sizes 42, 34, 16, 6, 6, 6, 6, 3, 2 and five single strains). Within a lineage the median difference is 557 SNPs. Between lineages it is at least 48,641 SNPs. The 5-SNP de-cloning of D-49 does not remove this structure. The effective sample size for a genome-wide scan is about 14, not 126.
2. **Genome-wide hits are lineage markers.** With the relatedness matrix only, 64 loci pass the Bonferroni threshold (unadjusted phenotypes) and 44 loci pass it (run-adjusted phenotypes). **All of them are SNPs that do not vary inside any lineage.** They mark which lineage a strain belongs to. They cannot be mapped to a gene.
3. **Within lineages, nothing passes the threshold.** With run adjustment and lineage covariates the scans find 0 Bonferroni loci in 20 traits. 17 SNPs reach the suggestive level, against about 94 expected by chance. The 12 suggestive loci are listed for reference (Table 13). They are hypotheses, not results.
4. **Lineage explains part of some traits.** Intraclass correlation by lineage: Cr relative area at the top dose 0.29, Cr area AUC 0.28, Cu a\* AUC 0.26, dose-0 ln area 0.24, dose-0 a\* 0.23, Cu relative area at the top dose 0.20. Pb traits and all growth-rate traits are at most 0.10. With 14 lineages these estimates are uncertain.
5. **Non-additive variance is not detectable.** The strains are haploid, so there is no dominance. For epistasis, the additive-by-additive variance component is not significant in any of 20 traits (smallest p = 0.095). A scan of about 65,000 SNP pairs per trait found no pair below p = 1e-5 (about 0.65 expected by chance).
6. **Dose-0 (YPD) controls are consistent for a\*, less so for size and not for growth rate.** Strain explains 38% of dose-0 a\* variance across the three screens, and the screen itself explains about 0%. For ln area, strain explains 18% and the screen 17%. For the maximum growth rate, strain explains 3%.
7. **Growth rate and IC50.** Relative growth rate falls with dose in Cr (0.41 of control at dose 1.2) and Pb (0.12 at dose 30). Cu barely changes it (0.93 at dose 30). IC50 of relative area could be fitted for almost all strains in Cr and Pb. In Cu 138 of 317 strains are right-censored, because relative area stays at or above 0.5 at the top dose.

8. **No recombination inside the large lineages.** In the three largest lineages (42, 34 and 16 strains) the fraction of SNP pairs showing all four gametes does not change with distance: it is the same for pairs less than 1 kb apart, 500 kb apart and on different contigs (0.0% in L01, 1.7% to 3.0% in L02, 2.9% to 5.0% in L03). Recombination would make close pairs compatible and unlinked pairs incompatible. The pattern fits strict clonality (section 7).
9. **The Cr tree signal is partly a run effect.** Pagel's lambda of baseline a\* in pure haploid *R. mucilaginosa* is 0.91 (p = 0.008, 163 strains). After removing run effects it is 0.89 (p = 0.12). The lineages with the lowest Cr a\* (L05, L07) were assayed mostly in the two later runs, and the later runs have lower a\* (section 8).
10. **Within lineages there is not enough signal to link Cr traits to loci.** In per-lineage permutation scans, 1 of 18 Cr tests (3 lineages, 6 traits) has a family-wise p below 0.05, which is what chance gives (0.9 expected). No candidate survives correction for the 60 lineage-by-trait tests. Candidate SNP and gene panels with window plots are in section 9. The metal and stress genes in the panels are not enriched (permutation p = 0.78).

## 1. Data

- **Phenotypes** come from the DuckDB `heavy_metal_measurement` table with strain, species and ploidy from the `strain_info` view (D-58). a\* is `ColorLab_a*GeoMedian`. The largest object per well and image is used. The late window and the 2,000 px colony cutoff follow the a\* report (`analysis/carotenoid_stress_vs_baseline`).
- **Growth rate** of a well is the largest slope of ln area against time over windows of 5 consecutive images. A window must start at 300 px or more and span at least 4 h. A well needs at least 6 images.
- **Relative traits** compare a strain with itself at dose 0: area AUC (relative area against dose scaled 0 to 1), relative area at the top dose, relative growth-rate AUC, a\* change at the top dose and a\* AUC.
- **IC50** is fitted for relative area with a two-parameter log-logistic curve. It is reported only when the curve crosses 0.5 inside the tested dose range.
- **GWAS traits** are 20 per strain (Table 1). Each is rank-normalised before the scan. Copper IC50 is not used (46% of panel strains are censored).
- **Genotypes.** Biallelic SNPs from popgen's `rmuc_core` file, kept at minor allele frequency 0.05 or more and at most 10% missing genotypes inside the panel: 244,956 SNPs. After collapsing identical genotype patterns there are 7,199 distinct patterns (the effective number of tests for the lineage-free scans). With lineage covariates there are 2,940.
- **Annotation.** snpEff effects come from the VCF. Gene IDs and product names come from the DH4148 GFF. No GO or KEGG annotation exists for DH4148 yet.

{{T_PHENO}}

*Table 1. GWAS traits in the 126-strain panel (raw values; 119-125 strains per trait).*

{{IMG:gwas_trait_correlations|Figure 1. Spearman correlation among the 20 traits (panel strains).}}

## 2. Dose-0 (YPD) consistency across experiments

Dose-0 wells of the three screens are pooled. Model: trait ~ 1 + (1 | strain) + (1 | screen) + (1 | run) + (1 | plate). Wells are replicate wells on dose-0 plates. Fe and Zn dose-0 plates exist but are partial, so the main model uses Cr, Cu and Pb. A version with all five screens is in `report/tables/dose0_variance_components.csv`.

{{IMG:d0_variance_components|Figure 2. Shares of dose-0 variance: strain, screen (metal), run, plate and residual.}}

{{T_D0_VC}}

*Table 2. Variance shares of dose-0 wells. `share_metal` is empty for single-screen models. Fits marked singular in the table have a variance component at zero.*

{{IMG:d0_strain_means_between_screens|Figure 3. Dose-0 strain means in one screen against another (dashed = equal values).}}

{{T_D0_COR}}

*Table 3. Agreement of dose-0 strain means between screens.*

{{IMG:d0_plate_means|Figure 4. Mean of each dose-0 plate, by screen.}}

{{IMG:d0_strain_heterogeneity|Figure 5. Spread of a\*, ln area and growth rate across the dose-0 wells of each strain (red = GWAS panel).}}

{{T_D0_HET}}

*Table 4. Strains whose dose-0 wells vary most in a\* (robust z of the SD across wells; strains with at least 6 wells).*

**Reading.**
- **a\* at dose 0 is a reproducible strain trait.** Strain explains 38% of the variance in the pooled model, and the screen explains about 0%. Strain means agree between screens (Spearman 0.47 to 0.55). Run and plate explain 3% and 4%. Half of the variance (55%) is replicate-well noise.
- **Colony size at dose 0 differs between screens.** The screen explains 17% of ln area variance, against 18% for strain. Plate means of ln area are 10.60 (Cr), 10.33 (Cu) and 10.37 (Pb). Strain means agree only weakly between screens (0.24 to 0.38).
- **Growth rate is not a reproducible strain trait here.** Strain explains 3% of its variance and between-screen correlations are 0.00 to 0.15. It is not suitable for mapping.
- **Cr dose-0 plates vary most in a\*.** In the Cr-only model run explains 15% of the variance, against 2% in Cu and 0% in Pb. The plate-mean SD of a\* is 1.98 in Cr, 0.87 in Cu and 1.10 in Pb.
- **Within-strain spread.** The median SD of a\* across a strain's dose-0 wells is 2.98. 13 of 312 strains have a robust z above 3 (Table 4). Spread alone does not show mixed cultures. The identity columns of `strain_curation.csv` (`identity_flag_categories`) should be checked for these strains; this has not been done here.
- **Panel strains.** In the panel strain subset the strain share of a\* is 28% (all strains 38%).

## 3. Growth rate and IC50

{{IMG:growth_curves_panel|Figure 6. Median ln colony area of panel strains against time, by dose.}}

{{IMG:relative_growth_rate|Figure 7. Maximum growth rate of each strain relative to the same strain at dose 0, by dose (panel strains).}}

{{T_RGR}}

*Table 5. Relative maximum growth rate by dose (median and quartiles across strains). Pb at doses 20 to 30 has only 67 to 81 strains with enough images for a growth rate.*

{{IMG:ic50_distribution|Figure 8. IC50 of relative area (all strains).}}

{{T_IC50_STATUS}}

*Table 6. IC50 fit status by metal (all strains).*

{{IMG:ic50_by_group|Figure 9. IC50 by species and ploidy group (groups with at least 5 fitted strains; Cr and Pb).}}

{{T_IC50_GRP}}

*Table 7. IC50 (dose units of the screen) by group. Copper IC50 values are close to the top dose (30) and half of the strains are censored, so Copper IC50 is coarse.*

**Reading.**
- **Cr:** growth rate rises at low doses (1.15 at 0.2) and falls to 0.64 at dose 1.0 and 0.41 at dose 1.2. Median IC50 is about 0.6 dose units for pure haploid *R. mucilaginosa* (IQR 0.58 to 0.66).
- **Cu:** growth rate stays at 0.93 to 1.16 of control through dose 30. Relative area falls only at the top doses. Most Cu IC50 values are outside the tested range.
- **Pb:** growth rate is unchanged to dose 10, 0.73 at dose 15, and 0.12 at dose 30. Median IC50 is about 13.3 for pure haploid strains.
- **Group differences.** *R. sphaerocarpa* and *R. taiwanensis* have the lowest Cr IC50 (0.52 and 0.55, 6 strains each). Hybrid diploids have a higher Pb IC50 (14.4) than pure haploids (13.3). Group sizes are small (5 to 12 strains), so these are descriptive.
- **IC50 against other traits.** Cr IC50 correlates with Cr area AUC (rho 0.74) but not with the a\* change (0.10). Pb IC50 correlates weakly with Pb area AUC (0.31).

## 4. Genetic structure: 14 lineages

{{IMG:gwas_panel_structure|Figure 10. First two principal components of the centred relatedness matrix and the variance of the first 20 components. Each visible point is a cluster of near-identical strains.}}

{{T_LIN_THR}}

*Table 8a. Number of single-linkage groups of panel strains by the largest allowed SNP difference. The count is 14 for any cut from 2,000 to 20,000 differences.*

{{T_LIN_SIZES}}

*Table 8b. Lineage sizes (lineage 1 is the largest).*

- The first two components explain 40% and 20% of the relatedness. The panel has visible clusters, not a continuum.
- Lineage and run are associated in all three screens (chi-square p < 0.001). Strains of one lineage were often on the same runs. Removing run effects therefore removes part of the lineage signal as well.

### Lineage as a source of trait variance

{{T_ICC}}

*Table 9. Intraclass correlation of strain-level traits by lineage (one-way ANOVA estimator, 14 lineages), before and after removing run effects (fixed run effects, lineage as a random group).*

- Lineage explains up to 29% of a trait (Cr relative area). With 14 lineages each ICC has a wide interval. This is a broad-sense measure: it mixes additive and non-additive genetic effects and anything else that differs between lineages.
- Pb traits and the growth-rate traits show no lineage component (ICC at most 0.10).

## 5. GWAS

Three scans per trait (GEMMA, univariate LMM, Wald test, centred relatedness matrix from all 244,956 SNPs):
- **(a)** relatedness only, unadjusted phenotypes;
- **(b)** relatedness, phenotypes with run effects removed;
- **(c)** as (b) plus lineage indicator covariates, so that only variation inside lineages is tested.

Bonferroni thresholds use the number of distinct genotype patterns (7,199 for (a) and (b), p < 6.95e-6; 2,940 for (c), p < 1.7e-5). The suggestive level is 1 / that number.

{{T_INFL}}

*Table 10. Genomic inflation (lambda) of the scans, and the number of SNPs below the Bonferroni threshold. "OLS, no relatedness" is the inflation of a plain regression on the run-adjusted phenotypes.*

{{T_LINCHECK}}

*Table 11. Lead SNPs of the loci (SNPs within 50 kb are one locus) that are lineage markers: SNPs without variation inside any lineage.*

{{IMG:gwas_inflation_runadj|Figure 11. Inflation with and without the relatedness matrix (run-adjusted phenotypes).}}

### Relatedness only (scans a and b)

{{IMG:gwas_manhattan_Cu_runadj|Figure 12. Cu traits, run-adjusted, relatedness only. The wide vertical blocks are sets of perfectly linked SNPs that mark lineages.}}

{{IMG:gwas_manhattan_Cr_runadj|Figure 13. Cr traits, run-adjusted, relatedness only.}}

{{T_LOCI_UNADJ_TOP}}

*Table 12. The 12 strongest Bonferroni loci of scan (a). Every one is a lineage marker: it carries the alternate allele in a set of whole lineages. The gene is the one at or nearest the lead SNP, and the gene is not the cause.*

- A locus here is a statement about which lineages differ in the trait, and not about a gene. All SNPs on a lineage's backbone have the same p-value. The lead SNP sits at an arbitrary position along the genome.
- The relatedness matrix does not solve this. With 14 lineages, a lineage-level trait difference is perfectly explained by thousands of SNPs.

### Within lineages (scan c)

{{IMG:gwas_manhattan_Cr_runadj_lineage|Figure 14. Cr traits, run-adjusted, lineage covariates.}}

{{IMG:gwas_manhattan_Cu_runadj_lineage|Figure 15. Cu traits, run-adjusted, lineage covariates.}}

{{IMG:gwas_manhattan_Pb_runadj_lineage|Figure 16. Pb traits, run-adjusted, lineage covariates.}}

{{IMG:gwas_manhattan_dose0_runadj_lineage|Figure 17. Dose-0 traits, run-adjusted, lineage covariates.}}

{{IMG:gwas_qq_runadj_lineage|Figure 18. QQ plots of scan (c).}}

{{T_LOCI_LIN}}

*Table 13. Suggestive loci of scan (c) (p below 1 / 2,940). None passes the Bonferroni threshold (p < 1.7e-5). `effect` is the snpEff effect of the lead SNP; the gene is the one overlapping or nearest the lead SNP, with distance in bp.*

- None of the 12 loci passes the threshold, and 17 SNPs reach the suggestive level across 20 traits, against about 94 expected by chance. The scan shows no excess of signal.
- Two loci lie 2 kb apart on contig CM179495.1 (Pb a\* AUC, ACY3AU_002286, and Pb relative area, ACY3AU_002287, both missense variants; p = 2.9e-4 each). The two traits are different summaries of the same Pb data, and the result is not replicated in an independent experiment.
- Candidates with a named product are: ACY3AU_005580 (Nuclear control of ATPase 2, Cr a\* AUC, missense), ACY3AU_003010 (a mannosyltransferase of the dolichol pathway, Cr a\* AUC, splice-region variant) and ACY3AU_004087 (cis-prenyltransferase, dose-0 a\*). Prenyl and dolichol pathway genes are of interest for a\* biology, but this is a hypothesis for later work.
- The SNP heritability estimates that GEMMA reports with lineage covariates are at the boundary (0 or 1) and are not used.

## 6. Non-additive architecture

Strains are haploid, so there is no dominance. Non-additive variance here means epistasis. Two analyses:
- **Variance components (REML).** Model: additive kernel (centred relatedness) plus additive-by-additive kernel (its element-wise square) plus error. The epistatic component is tested by a likelihood ratio test with the boundary correction.
- **Pairwise scan.** All pairs of the 1,401 LD-pruned SNPs (r² < 0.5), interaction term in an LMM that is adjusted for relatedness. Pairs need at least 5 strains in each of the four two-locus classes. Between 51,888 and 65,081 of 980,700 pairs per trait meet this (the others are collinear or rare because of the lineage structure).

{{T_VC_EPI}}

*Table 14. Variance shares from the REML models (RINT traits, panel strains). `h2_additive_only` is the SNP heritability with the additive kernel alone.*

{{T_EPI}}

*Table 15. Pairwise interaction scan. No pair passes the Bonferroni threshold. `n_expected_at_1e_5` is the number expected below 1e-5 by chance.*

- SNP heritability from the additive kernel alone is 0 to 0.26: d0 a\* 0.26, Cu a\* AUC 0.23, Cr a\* AUC 0.18, Cu a\* change 0.18, Cr and Cu relative area at the top dose 0.16. Standard errors are 0.08 to 0.17 (GEMMA, `gwas_inflation_pve_unadjusted.csv`).
- The epistatic component is zero or small in all traits. The smallest p is 0.095 (d0 a\*). When it is non-zero it takes variance away from the additive component (d0 a\*: 0.26 additive only, 0.00 additive and 0.22 epistatic). The two kernels are highly correlated in a clonal panel, so the split between them is not identifiable.
- The pairwise scan has no hit and no excess (lambda 0.5 to 1.7; median 1.02).
- With 14 lineages, no analysis can tell an additive lineage effect from an epistatic one. Both are a difference between whole genotypes.

## 7. Recombination inside the large lineages

Informative SNPs are those with minor allele count of at least 2 inside the lineage (all SNPs of the QC VCF, not only common ones): 1,691 in L01 (42 strains), 2,756 in L02 (34) and 1,661 in L03 (16). Two tests:
- **Four-gamete test and r² by distance.** Under strict clonality every pair of sites is equally linked, so the fraction of pairs with all four two-locus classes does not depend on distance. Under recombination close pairs stay compatible and unlinked pairs (other contigs) become incompatible.
- **Windowed distance matrices.** For 300 kb windows the strain-by-strain SNP distance matrix is compared between windows. With recombination, distant windows would show different relationships between strains than adjacent windows.

{{IMG:recombination_four_gamete|Figure 19. Fraction of SNP pairs with all four gametes, by physical distance (and pairs on different contigs), for the three largest lineages.}}

{{T_RECOMB}}

*Table 16. Pair classes: mean r², fraction of pairs with all four gametes (each gamete in at least 1 strain, and in at least 2 strains).*

{{T_RECOMB_WIN}}

*Table 17. Spearman correlation of strain-by-strain distance matrices between 300 kb windows.*

**Reading.**
- **L01 has no incompatible pair at any distance** (0 of about 190,000 sampled pairs). Its r² is 0.27 below 1 kb and 0.14 to 0.15 at every other distance, including other contigs.
- **L02 and L03 have a few percent incompatible pairs, with no trend with distance** (L02: 3.0% below 1 kb, 1.7% to 2.3% beyond; 2.0% on other contigs). Requiring each gamete in at least 2 strains lowers this to 0.1% to 0.3% beyond 1 kb in L02 (1.5% below 1 kb) and 0.5% to 0.8% in L03 (0% below 1 kb). The rate on other contigs equals the rate within contigs. This is the signature of recurrent mutation or genotype errors, not of recombination.
- **The distance-matrix correlation does not depend on window distance** (L01 0.64 adjacent, 0.58 distant, 0.68 other contigs; L02 0.58, 0.59, 0.60; L03 0.83, 0.86, 0.81). The strains have the same relationships at every position.
- **A small excess of r² below 1 kb** (0.27 in L01, 0.29 in L02) is visible. This may come from clustered mutations or mapping errors. It does not extend beyond 1 kb.
- **Conclusion.** There is no evidence of recombination in these lineages. They behave as clones that accumulated mutations. Rare events (for example gene conversion tracts shorter than 1 kb or a few recombinants among many strains) cannot be excluded with these tests.

## 8. The Cr tree signal

{{T_CR_TREE}}

*Table 18. Pagel's lambda of baseline a\* in pure haploid *R. mucilaginosa* (outgroup-rooted PHYling tree), before and after removing run effects (run mean subtracted).*

{{T_RUNLIN}}

*Table 19. Share of variance in strain-mean baseline a\* explained by run and by lineage (adjusted R²). Lineage needs a lineage label, so n differs.*

- **Cr:** lambda is 0.91 (p = 0.008, n = 163). After removing run means it is 0.89 (p = 0.12). The point estimate hardly changes, but the evidence for lambda above 0 weakens, so the signal is not robust. **Cu and Pb:** lambda is 0 in both versions.
- **Run matters for Cr.** Mean baseline a\* by run: d000320 20.9, d000321 21.2, d000322 19.6 (SD 3.9), d000323 18.1. Run explains 14% of the variance in strain means (Cu 2%, Pb 1%). Lineage explains 15% (Cu 13%, Pb 12%).
- **Low-a\* lineages sit in late runs.** In Cr, L05 (mean a\* 15.6) has 5 of 6 strains in run d000322, and L07 (16.0) has 6 of 6 in d000322 or d000323. The large lineages L01 (20.1) and L02 (20.7) are spread over all four runs. So a few lineages with low a\* coincide with the runs with low a\*. Two explanations fit: a genetic difference of those lineages, or a run effect. The data cannot separate them.
- **Cr lineage structure on the tree.** The tree signal comes mostly from the separation between a small number of lineages and clades, consistent with section 4. It does not show a continuous relationship between relatedness and trait inside the species.

## 9. Candidate SNP and gene panels, with window context (Cr)

**Panel A: lineage-level partitions.** The Cr Bonferroni loci of the unadjusted scan (the run-adjusted scan has none for Cr) fall into a few genotype partitions, each carried by a set of whole lineages. For each partition the table counts the SNPs that share its genotype pattern and the genes with high or moderate snpEff effects among them.

{{T_PART}}

*Table 20. Cr lineage partitions (lead trait and scan in brackets). `n_perfectly_linked_snps` counts SNPs with the same genotype pattern.*

**Panel B: within-lineage candidates.** SNPs with p < 1e-3 in scan (c) for the six Cr traits, and SNPs with p < 0.01 in the per-lineage permutation scans (L01, L02, L03). The full table is `report/tables/cr_candidate_snp_panel.csv` (91 rows, 74 SNPs, 175 genes within 5 kb).

{{T_PANEL_B}}

*Table 21. The 15 strongest rows of panel B. Genes within 5 kb of the SNP are listed with their products.*

{{T_WL_CR}}

*Table 22. Per-lineage permutation scans for Cr traits (family-wise p from 2,000 permutations of the trait across strains).*

{{T_ENRICH}}

*Table 23. Are genes with a metal, stress or pigment-related product keyword enriched among the Cr scan (c) signal? Mean -log10 p of the best SNP in or within 2 kb of each gene, keyword genes against all genes with SNPs; p from 10,000 random gene sets of the same size.*

{{IMG:zoom_cr_locus1|Figure 20. Locus 1: Cr relative growth-rate AUC near CM179490.1:904,292. Top: scan (c). Middle: scan (b). Bottom: gene models in the window (orange: product matches the keyword list).}}

{{IMG:zoom_cr_locus3|Figure 21. Locus 3: Cr a\* AUC near CM179487.1:957,384 (missense variant in ACY3AU_005580, Nuclear control of ATPase 2).}}

{{IMG:zoom_cr_locus4|Figure 22. Locus 4: Cr a\* AUC near CM179497.1:588,056 (splice-region variants in ACY3AU_003010).}}

{{IMG:zoom_cr_locus5|Figure 23. Locus 5: Cr relative area at the top dose near CM179491.1:380,349.}}

{{IMG:zoom_cr_locus7|Figure 24. Locus 7: the only family-wise hit of the Cr within-lineage scans, in L02 (Cr relative area at the top dose, CM179495.1:102,060; 15 of 33 strains carry the alternate allele).}}

{{IMG:zoom_cr_locus8|Figure 25. Locus 8: best Cr SNP of the L03 scan (16 strains).}}

**Reading.**
- **Lineage level.** The Cr traits differ among lineages (ICC 0.28 for relative area). The partitions that track the difference are carried by 16 to 50 strains from 5 to 8 lineages. Each has 7 to 619 perfectly linked SNPs and up to 12 genes with high or moderate effects (partition P8: 3 high, 92 moderate). The data cannot say which SNP matters, and cannot tell a genetic effect from a run effect (section 8).
- **Within lineages.** Only one Cr test reaches a family-wise p below 0.05: L02, relative area at the top dose, 2 SNPs 5 bp apart (p = 4.6e-4, family-wise 0.0085) upstream of ACY3AU_002196 (sister chromatid cohesion protein 2). Over the 60 lineage-by-trait tests, 3 have a family-wise p below 0.05 (the other two are Cu and Pb traits in L02), which equals the 3 expected by chance. With Bonferroni over 60 tests the best p is 0.51. L01 and L03 have none for Cr.
- **Window context.** In the zoom plots the SNPs at the top of each locus are single SNPs or a small cluster with no supporting neighbours. A true causal locus would show a peak of SNPs in linkage with the lead SNP. Here the lead SNPs are mostly isolated points, which is the pattern of noise.
- **Keyword genes.** 107 of 3,605 genes with SNPs have a metal, stress or pigment-related product. Their mean signal is not higher than for other genes (1.17 against 1.21; permutation p = 0.78).
- **Is there enough signal?** No, not in this panel: Cr traits clearly differ among lineages, but no locus can be assigned, and clone-corrected analysis inside Rmuc shows no result that survives correction.

## 10. Limits and next steps

- **Panels in section 9 are hypotheses.** None survives correction.
- **The panel cannot support gene mapping.** Every genome-wide signal tracks lineage. A usable panel needs many more independent lineages or recombinant strains (crosses, or populations with sexual recombination). D-49 (5-SNP cutoff) should be revisited: at 2,000 SNPs only 14 strains remain.
- **Run is confounded with lineage** (all screens). Run adjustment is therefore partial. A design with lineage spread over runs, or repeated runs, would separate them.
- **Replicate noise.** The strain share of the dose-0 variance is 38% for a\* and lower for size and growth. Strain-level traits from 1 to 4 wells per dose are noisy, which lowers power further.
- **Cu IC50** is censored for 46% of panel strains; only coarse.
- **Fe and Zn** are held out (D-56). When complete tables exist, add them to `METALS` in `01_phenotypes.py`.
- **Population analyses** (DH4148 SNP clusters, D-51) are not part of this report. The lineage table (`results/lineages.csv`) is a first version of that grouping.
- **Annotation.** Gene products come from the assembly GFF only. A GO and KEGG annotation of DH4148 would help to interpret candidates.
- **Scope of claims.** Candidate genes in section 5 are suggestive only. They have not been replicated.
