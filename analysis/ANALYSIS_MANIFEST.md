# Analysis Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

### explore-plate-position
```yaml
name: explore-plate-position
question: How much of strain-replicate variance in color (L*a*b*) and morphology is attributable to plate identity / within-plate grid position, stratified by Copper concentration?
input: db/rhodotorula_phenotypes.duckdb (v_phenotype), endpoint timepoint per plate
scripts:
  - scripts/00_build_dataset.R      # last-image-per-plate endpoint colonies -> endpoint_colonies.rds/csv
  - scripts/01_variance_partition.R # lmer variance partition, STRATIFIED by copper_mm (plate nested in Cu)
  - scripts/02_adjacency_effect.R   # 4-connected grid-neighbor mean predicts own trait (controls copper_mm)
  - scripts/03_plots.R              # per-Cu variance bars, plate-% heatmap, plate heatmap, edge/adjacency plots
  - scripts/04_build_timecourse.R   # per-colony x timepoint table (runs 353-356, 6-hourly) -> colony_timecourse.rds + detection stats
  - scripts/05_time_variance_partition.R # day-block (d1-d5) x Cu variance partition rerun of 01's structure
  - scripts/06_growth_curves.R      # L*/log(area) growth mixed models (quadratic time x Cu, random strain slope/plate/colony)
outputs:
  - results/tables/variance_components_by_copper.csv
  - results/tables/fixed_effects_by_copper.csv
  - results/tables/categorical_position_anova_by_copper.csv
  - results/tables/adjacency_fixed_effects.csv
  - results/tables/colony_timecourse.rds/csv
  - results/tables/variance_components_by_day.csv
  - results/tables/growth_curve_fixed_effects.csv
  - results/tables/detection_curve_by_copper.csv
  - results/figures/variance_components_by_copper.png
  - results/figures/plate_pct_heatmap.png
  - results/figures/day_plate_pct_heatmap.png
  - results/figures/day_variance_components_d1_d5.png
  - results/figures/growth_L_by_cu.png
  - results/figures/growth_area_by_cu.png
  - results/figures/detection_curve_by_cu.png
  - explore_plate_position.html
reproduce: bash analysis/explore_plate_position/run.sh
status: complete
key_findings:
  - Plate-explained variance is SMALL within each Cu concentration (L* ~0.1-4.7%, b* 0-3.2%, area/solidity ~0.4-5.3%) after stratifying by Cu.
  - Earlier pooled (across-Cu) model over-estimated plate variance (23% L*, 33% b*) because plate_id is nested in Cu and the linear copper_mm covariate left non-linear Cu response to be absorbed by the random plate term.
  - Neighbor adjacency signal for L* persists (p~6e-117) with no a*/b*/area/solidity neighbor effect.
  - Plate % declines from day 1 to day 5 e.g. L* at 30 mM 19.1% -> 1.6%: plate position explains least variance at the mature stage captured by the endpoint analysis.
  - Growth curves confirm expectation: L* rises ~44->74 plateauing ~day 3-3.5, log area rises monotonically; Cu slows/limits both. At 25-30 mM ~15-18% of (strain x plate) never reach late timepoints (missingness is the phenotype).
  - run 357 (plates 113-120) excluded from the time course; all plates share the same relative 6-hourly clock except run 353's offset schedule (~3h grid to 117h).
tags: [copper, plate-position, variance-partition, mixed-model, lme4, rhodotorula, timecourse, growth-curve, detection]
```

### growth-rates
```yaml
name: growth-rates
question: How do per-strain colony growth-rate parameters (data-derived peak slopes of log area and of consumer-light intensity) vary with copper concentration and with species, and how do they interact with endpoint color (CIELAB L*, a*, b*) where L* is the light-intensity readout?
input: db/rhodotorula_phenotypes.duckdb (v_phenotype), runs 353-356, per (plate, well) x timepoint
scripts:
  - scripts/00_build_series.R      # per-colony x timepoint table + species join -> colony_growth_series + coverage
  - scripts/01_fit_growth_models.R # Gompertz + Logistic per colony/trait (mu/lambda/A, AIC-preferred), fallback log-linear; primary rate = peak 6 h slope of log(area) [rate_area] and of intensity [rate_int]
  - scripts/02_species_cu_rates.R  # strain x Cu aggregation (up to 4 run-replicates), mixed models rate ~ Cu(factor) + (1|species/strain), Cu x species interaction (well-sampled spp >= 8 strains, unnamed 'sp. clade I' excluded), contrasts vs 0 mM; rate_by_cu_spp / rate_int_by_cu_spp facet top-16 species (4x4)
  - scripts/03_color_interaction.R # endpoint L*/a*/b* as outcomes of rate_area/rate_int x copper (+ species/strain random); documents L*-Intensity ~0.999 collinearity; also endpoint color (L*/a*/b*) vs Cu plots
  - scripts/04_doubling_time.R     # phase-independent estimator: exponential-region specific rate (best-R2 4-pt window) -> doubling time, t50, saturation fraction; tests whether peak-slope Cu trend persists
  - scripts/05_species_cu_sensitivity.R # extent/yield-based Cu-sensitivity index per strain & species (log2 max-area ratio 25-30 vs 0-5 mM, saturation drop, dbl fold), species extent model (Cu x species); sensitivity figures facet top-16 species (4x4)
outputs:
  - results/tables/colony_growth_series.rds/csv
  - results/tables/series_coverage.csv
  - results/tables/growth_model_fits.csv
  - results/tables/growth_model_preferred.csv
  - results/tables/strain_x_cu_rates.csv
  - results/tables/rate_models_species_cu.csv/txt
  - results/tables/rate_area_diff_0mM.csv
  - results/tables/rate_int_diff_0mM.csv
  - results/tables/color_growth_models.csv/txt
  - results/tables/color_growth_anova.csv
  - results/tables/colony_doubling.csv
  - results/tables/dbl_model_cu.csv
  - results/tables/dbl_by_cu.csv
  - results/tables/strain_sensitivity.csv
  - results/tables/species_sensitivity.csv
  - results/tables/species_extent_model.txt
  - results/tables/species_extent_anova.csv
  - results/figures/rate_by_cu_overall.png
  - results/figures/rate_by_cu_spp.png
  - results/figures/rate_int_by_cu_spp.png
  - results/figures/color_L_vs_ratearea_by_cu.png
  - results/figures/color_a_vs_ratearea_by_cu.png
  - results/figures/color_b_vs_ratearea_by_cu.png
  - results/figures/ratearea_vs_rateint.png
  - results/figures/color_L_by_cu.png
  - results/figures/color_a_by_cu.png
  - results/figures/color_b_by_cu.png
  - results/figures/color_L_by_cu_spp.png
  - results/figures/color_a_by_cu_spp.png
  - results/figures/color_b_by_cu_spp.png
  - results/figures/dbl_region_by_cu.png
  - results/figures/saturation_by_cu.png
  - results/figures/sensitivity_species_rank.png
  - results/figures/sensitivity_extent_by_cu_spp.png
  - results/figures/sensitivity_saturation_by_cu_spp.png
  - growth_rates.html
reproduce: bash analysis/growth_rates/run.sh
status: complete
key_findings:
  - Peak growth/brightening rate increases monotonically with Cu (area F=79 p~6e-96; intensity F=540 p~0) - contrary to naive toxicity expectation; consistent with high-Cu colonies staying in log-growth (unsaturating) phase longer, i.e. peak-slope is phase/length sensitive, not a phase-independent rate constant.
  - Gompertz/Logistic lag is structurally unidentifiable in-window (~94% boundary lambda, ~82% area asymptote extrapolated; unbounded GN gives degenerate lags -297..-20 h) -> primary rate estimator is data-derived max 6 h slope.
  - rate_area vs rate_int only r=0.21: area expansion and consumer-light-intensity rise are distinct growth axes.
  - L* (endpoint) vs Intensity_MeanIntensity r=0.999: same axis; L* is "light intensity"; intensity only used as growth readout.
  - Growth-rate x Cu interaction on endpoint L* significant (F=8.16 p~8e-9): faster-growing colonies end lighter, slope largest at low Cu (5mM +3.7 L*/unit rate) smallest at 30mM (+1.3); Cu dominates chromatic outcome (~1000x F on L*/a*).
  - 8,446 colonies (97.6% of 8,652 design slots; 309 strains, 112 plates); 130/2183 strain-Cu aggregates are single-culture, 282 colonies lack species label (-> unknown level).
  - Doubling-time test (04): exponential doubling time is flat ~37 h across 0-30 mM (fold-change 0.99 at 30 mM; Kendall tau between peak-slope rate_area and doubling time = 0.002) -> the rising peak-slope with Cu IS A PHASE ARTIFACT. Cu's real effect is extent-limited: saturation fraction 95.4% -> 73.9%, lower max area, median t50 72.6 -> 64.1 h.
  - Species Cu sensitivity (05, extent-based log2 ratio 25-30/0-5 mM): most tolerant R. glutinis (+1.01) and, among well-sampled, R. mucilaginosa (-0.72); most sensitive R. kratochvilovae (-4.13), R. araucariae (-3.68), R. taiwanensis (-3.06); well-sampled R. diobovata (-2.00). Cu x species interaction on log(max_area) F = 4.34 (p ~ 6.0e-4).
  - Species mislabels in data/metadata/Copper.Strain_info.csv (paludigenum -> paludigena, evergladiensis -> evergladensis) corrected and re-imported via scripts/db/10_import_experiment.py --strain-only (strain 84 & 327 regrouped).
tags: [copper, growth-rate, growth-curve, species, mixed-model, lme4, rhodotorula, color, light-intensity, gompertz, logistic, timecourse]
```

### control-late-timepoint-phenotype
```yaml
name: control-late-timepoint-phenotype
question: What are the per-strain colony size and CIELAB color (L*, a*, b*) summary statistics on control media (Cu=0, YPD) at a late timepoint, aggregated across all replicate colonies?
input: db/rhodotorula_phenotypes.duckdb (v_phenotype + condition_plate_factor Cu=0 + strain)
scripts:
  - scripts/build_phenotype_table.py # latest imaging pass (rounded hour) per strain within a [tmin,tmax] window on Cu=0; per-strain median/mean/var/sd of Shape_Area and per-colony ColorLab_{L*,a*,b*}Median; run per window (70-80, 80-90, 90-110)
outputs:
  - results/phenotype_control_timepoint_70_80.csv
  - results/phenotype_control_timepoint_80_90.csv
  - results/phenotype_control_timepoint_90_110.csv
  - CONTROL_LATE_PHENOTYPE.md
reproduce: bash analysis/control_late_timepoint_phenotype/run.sh
status: complete
key_findings:
  - 314/320 strains have a late Cu=0 image in each window (70-80, 80-90, 90-110 h); 286 have >=3 replicate colonies. Latest pass is the alternate-cadence (75/87/105 h, 84 strains) or main-cadence (78/90/108 h, 230 strains) pass; using an exact max-hour per strain would fragment a single imaging pass (1-2 colonies/strain) and was rejected in favor of rounding to the imaging pass.
  - Strain sets are identical across the three windows; only the sampled pass differs. Colony size highly stable across windows (area_median r ~0.98 for 70-80 vs 90-110).
  - Color drives the table: L* medians 70.2-79.8, a* 0.6-14.2 (carotenoid red/orange), b* -1.4-9.1 on YPD control; colony area 0.8-70 k px.
  - Per-strain stats reproduced exactly by independent manual aggregation (strain 185).
tags: [copper, control, YPD, color, CIELAB, L*a*b*, colony-size, phenotype-table, duckdb, strain, rhodotorula]
```

### 2026-08-15-color-phenotype-space
```yaml
name: 2026-08-15-color-phenotype-space
question: What else should we explore for phenotypic data space or color space in this dataset? (persona-driven ideation campaign, 8 personas x 2 ideas = 16 ideas, all implemented)
input: db/rhodotorula_phenotypes.duckdb (v_phenotype, 211,800 rows) via shared extract data/db_extract.parquet (116 cols)
scripts:
  - scripts/build_series.py        # shared extract: v_phenotype -> data/db_extract.tsv.gz/.parquet + strain_metadata.tsv
  - scripts/common.py              # read_extract/read_meta/save/boot_ci/circular_stats/hue asat helpers
  - scripts/idea_01_arrest.py      # Weibull arrest collapse + Gibrat dispersion (statistical-physicist)
  - scripts/idea_01b_gibrat_within.py # within-strain + size-controlled resolution of the Gibrat dispersion ambiguity
  - scripts/idea_02_color.py       # hue/chroma/L* triple, morphs, onset, run calibration (color-imaging)
  - scripts/idea_03_information.py # species info, redundancy, forward selection (information-theorist)
  - scripts/idea_04_qg.py          # ICC/repeatability, reaction norms, GxE (quantitative-geneticist)
  - scripts/idea_05_repl.py        # rank-PCA/varimax factors + UMAP/HDBSCAN atlas (representation-learning)
  - scripts/idea_06_mediation.py   # Cu -> growth -> pigment endpoint mediation (causal-inference)
  - scripts/idea_07_ecology.py     # environment stratification, species-conditional permutation (trait-ecology)
  - scripts/idea_08_onset.py       # onset threshold sweep + pigment-area coupling (temporal-dynamics)
  - scripts/idea_09_phylogeny.R    # phylogeny join + Mantel perm signal (phylogeneticist, post-campaign) on data/raw/rhodotorula-phyling-protein-tree
  - scripts/idea_10_species_variation.py # per-species boxplots + between/within-species variance decomposition (scale of variation) on the 5 key traits
outputs:
  - results/idea01_{arrest_collapse,dispersion}.csv; figures/fig01_{master_collapse,arrest_exponent,gibrat_dispersion}.png
  - results/idea02_{species_color_triple,heterogeneity,pigment_morphs,run_calibration,onset_times}.csv; fig02{a,b}_*.png
  - results/idea03_{feature_species_info,redundancy_summary,top_correlated_pairs,forward_selection}.csv; fig03_{mi_top,redundancy_heatmap}.png
  - results/idea04_{icc_audit,heterogeneity_repeatability,reaction_norms,gxe_interaction_mucilaginosa}.csv; fig04_{icc_audit,reaction_norms}.png
  - results/idea05_{pca_variance,varimax_loadings,factor_block_top,atlas,cluster_species_matrix}.csv; fig05_{atlas_umap,cluster_species}.png
  - results/idea06_mediation_decomposition.csv; fig06_growth_pigment_decouple.png
  - results/idea07_{trait_environment,env_trait_means,env_trait_z}.csv; fig07_trait_environment.png
  - results/idea08_{onset_threshold_sweep,pigment_pace,logistic_onset,well_growth_onset}.csv; fig08_{onset_vs_capture,example_growth_pigment}.png
  - results/idea09_{phylo_signal,tip_traits}.csv; results/idea09_pruned_trees.rds; fig09_trait_on_tree_{cu_slope,baseline_chroma,dispersion}.png
  - results/idea10_{variance_decomposition,species_variation}.csv; figures/fig10_boxplot_{slope_logchroma_per_mM,intercept_logchroma,l10med_fixed,partial_slope_sd_cu,pace_loglog}.png
  - FINDINGS.md  # master synthesis of all 8 results
reproduce: pixi run python scripts/idea_XX_name.py per script
status: complete
key_findings:
  - sd(log10 Asat) "dispersion widens with Cu" (idea 01 pooled) is largely a SIZE/FLOOR ARTIFACT: sd(log10 Asat) ~ -0.63 correlated with colony size; within-strain raw slopes +0.0033/Mm (p~5e-16) collapse to ~0 size-controlled (p=0.46); only R. taiwanensis genuinely widens (+0.014, p=0.03) and R. paludigena narrows (-0.025, p=0.001) (idea 01b).
  - Cu acts mostly THROUGH growth: log(chroma) loss at 90-110 h is ~fully mediated by colony area (mediation frac ~1.2, direct b~0 per well, boot CI [0.49,2.69]); only the a* (redness) channel keeps a ~40% direct, growth-independent Cu effect.
  - Species/stress/environment signal rides on intra-colony heterogeneity and shape, NOT on L*/intensity: L*Median has the least species MI (~0.07 bit vs best feature 0.19 bit of H=3.0); environment associates only with b*CoeffVar (p_strat~0.006) and L*Median after species-blocking. Camouflage caveat: a*/b*CoeffVar are colony-noise (ICC_strain 0.19/~0) while shape/chroma are repeatable (ICC 0.82-0.93).
  - The strain atlas is a CONTINUUM: varimax factors F1 color-level / F2 heterogeneity / F6-F8 shape are block-diagonal, but UMAP+HDBSCAN finds no discrete species clusters; GxE on log(chroma) within R. mucilaginosa is real but modest (F=1.66, p<<1e-6).
  - Define onset with chroma>7-10: basal chroma noise ~1.5-2 in late Cu=0 colonies; at thr7 median onset 48 h, 91% ever pigment, while >=20 mM Cu gives 0% ever pigment. Colony size gates onset only weakly (rho_onset,size 0.30-0.42 across thresholds).
  - Strain phenotypes carry phylogenetic signal when measured tree-robustly (idea09): Mantel permutation, all-strains scope, l10med_fixed r=0.43, partial_slope_sd_cu r=0.31, intercept_logchroma r=0.19, slope_logchroma_per_mM r=0.10 (all p<=0.003, n~266-272); pace_loglog n.s. Within R. mucilaginosa only baseline chroma (p=0.002) and size (p=0.027) retain signal -> between-species structure, consistent with within-species continuum (idea05).
  - IMPORTANT methodological caveat (idea09): the PHYling protein tree is a near-comet topology (giant polytomy of near-duplicate R. mucilaginosa; 167/541 edges <=1e-7, 22 zero-length pendant tips). Blomberg's K collapses to ~1e-7 here regardless of truth (BM power check K=2.25) and likelihood lambda is numerically unstable (geiger lnL inconsistent) -> use rank-based Mantel permutation for phylogenetic signal on near-comet/genome-cluster trees, never raw K.
  - Species monophyly on this tree (idea09): dairenensis, diobovata, graminis, kratochvilovae, sphaerocarpa, sp. clade I monophyletic; mucilaginosa/paludigena/taiwanensis/toruloides NOT.
  - Within-species (among-strain) variation is the DOMINANT scale of variation for all 5 key traits (idea10): exact SS decomposition, 11 species n>=3, fraction of variance WITHIN species = Cu slope 87%, baseline chroma 92%, colony size 62%, within-strain heterogeneity 67%, pigment pace 84% (ANOVA F 2.6-17.5). Colony size is the most species-structured (38% between); chroma/Cu-slope least (~8-13% between). Consistent with idea09 (between-species Mantel signal in the smaller component) + idea05 (within-species continuum). CV% unreliable for near-zero-mean traits (pace diobovata).
tags: [ideas, ideation, persona, phenotype-space, color-space, CIELAB, mediation, information-theory, heritability, GxE, UMAP, HDBSCAN, varimax, ecology, onset, rhodotorula, copper, duckdb, pixi, phylogenetics, mantel, phylogenetic-signal, within-species, variance-partition, boxplot]
```

### 2026-08-15-color-phenotype-space-idea11-power
```yaml
name: 2026-08-15-color-phenotype-space (idea 11) - GWAS power & variant feasibility
question: Given the within-species dominance (idea09/10), is a GWAS feasible for our color/growth traits, and what effect sizes are detectable? (Also: do we even need variant calling?)
input: analysis/ideas/2026-08-15-color-phenotype-space/results/idea09_tip_traits.csv (278 strains x 5 traits); PHYling protein tree final_tree.nw (278 tips); external SNP panel /bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510
scripts:
  - scripts/idea_11_effective_n.R  # effective independent haplotypes from protein tree (collapse near-clones)
  - scripts/idea_11_power.py       # noncentral-chi2 power table + detectable effect sizes
outputs:
  - results/idea11_effective_n.csv, results/idea11_genome_redundancy.csv  # 200 tips -> 178 effective for R. mucilaginosa
  - results/idea11_power.csv  # min detectable R2 per n at alpha 5e-8/1e-6/1e-4
  - results/idea11_detectable_effect.csv  # beta in SD & trait units, frac of additive genetic variance (h2=ICC)
  - figures/fig11_power_curves.png
reproduce: pixi run Rscript scripts/idea_11_effective_n.R; pixi run python scripts/idea_11_power.py
status: complete
key_findings:
  - Effective independent genomes: R. mucilaginosa 200 tips -> 178 effective haplotypes (redundancy 1.12, 22 near-clones collapsed at dist 1e-7); everyone else redundancy 1.0.
  - Min detectable per-SNP R2 @80% power at n=202: 0.164 genome-wide (5e-8), 0.140 exome (1e-6), 0.100 candidate (1e-4); at real eff n=178 GW ~0.21. Only large-effect loci (~0.71 SD/allele @ MAF 0.3, ~22-25% of additive genetic var assuming h2=ICC) are detectable.
  - Copper slope unmappable (ICC 0.19 < required heritability); colony size/chroma/heterogeneity/pace are the GWAS-able traits.
  - DISCOVERY: no variant calling needed. Existing population-genomics project for R. mucilaginosa (ref NRRL Y-2510, /bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510) has vcf/RmucY2510_v2.All.SNP.combined_selected.vcf.gz = 728,581 SNPs x 422 haploid strains (GATK hard-filtered). 201 of our 278 phenotyped strains (all R. mucilaginosa) carry genotypes; 200/201 have complete data for all 4 GWAS traits -> n matches the n=202 power row.
  - Prior art: their 218-strain growth-rate GWAS (GEMMA+LMM+linear, T4C-T37C/Salt6) found only a handful of hits (e.g. chr13, p~1.7e-11 T37C linear) - consistent with idea-11 large-effect-only limit; our color traits untested (gap).
tags: [power, gwas, variant, effective-n, redundancy, noncentral-chi2, snp, r-mucilaginosa, NRRL-Y2510, idea11, within-species]
```

### 2026-08-15-color-phenotype-space-gwas-loco-tierb-coloc
```yaml
name: 2026-08-15-color-phenotype-space GWAS follow-up (LOCO sensitivity + Tier B set tests + dxy/Fst co-localization)
question: Do Tier-A GWAS hits (chroma/AUC_10/resilience_30) survive LOCO kinship sensitivity? Do pixy high-dxy windows carry multi-SNP (set-level) signal beyond single-SNP Tier-A? Are GWAS loci co-localized with population-divergence (dxy/Fst) outliers?
input: results/gwas/tierA_summary/{gwas,gwasc}_*_assoc.csv.gz; results/gwas/pixy/genome_{dxy,fst}.txt; results/gwas/loco/output/loco_*.assoc.txt (120 GEMMA LOCO scans); LD cache results/gwas/tierB/ld_cache/
scripts:
  - scripts/merge_loco.py        # per-chr lambda + top hits vs Tier-A -> loco_merged_{gwas,gwasc}.csv
  - scripts/make_loco_figure.py  # two-panel LOCO sensitivity PNG/PDF
  - scripts/tierb_set_tests.py   # burden/SKAT(min-p) set tests on 178 high-dxy windows (12 traits x both sets), MC-verified top-50
  - scripts/tierB_submit.sh / tierB2_submit.sh  # SLURM runners (gwas/gwasc)
  - scripts/coloc_dxy_fst.py     # FDR-locus to pixy-window join + Fisher enrichment
  - scripts/make_tierb_figure.py / make_coloc_figure.py  # result figures
outputs:
  - results/gwas/loco/loco_merged_{gwas,gwasc}.csv
  - results/gwas/tierB/tierb_settests_{gwas,gwasc}.csv; tierb_skat_mcver_{gwas,gwasc}.csv
  - results/gwas/tierB/coloc/coloc_final.csv.gz, coloc_anchors.csv, coloc_enrichment.txt
  - results/gwas/figures/{loco_sensitivity,tierB_settests,coloc_dxy_fst}.{png,pdf}
  - GWAS_REPORT.md section 8
reproduce: sbatch scripts/tierB_submit.sh; python scripts/{merge_loco,make_loco_figure,coloc_dxy_fst,make_tierb_figure,make_coloc_figure}.py
status: complete
key_findings:
  - LOCO reproduces every Tier-A anchor at unchanged p (chroma 2.46e-8, AUC_10 1.44e-8, resilience_30 6.04e-9 on all-201); signals are not kinship-absorption artifacts. LOCO lambda (medians 0.34-0.73) tracks Tier-A lambda - the lambda<1 deflation is stable near-clonal structure, not per-chromosome rescue.
  - Tier B: no set-level signal survives FDR(q<0.05) in either set. burden over-conservative (lambda~0.5; + and - z cancel with no direction prior), SKAT mildly deflated (lambda~0.6-0.8; top p~1.4e-3 not consistent across sets), min_p recovers Tier-A single-SNP hits (383/384 & 408/408 high-dxy windows). MC-verified moment-approx SKAT (r=0.984 log10, n=50 windows).
  - dxy/Fst co-localization: FDR GWAS loci NOT enriched in high-divergence windows (OR=0.81 dxy p=0.61; OR=0.41 Fst p=0.065). Phenotype alleles segregate within the near-clonal focal clade (standing variation), decoupled from the deep pop splits driving dxy/Fst extremes.
  - scaffold_20:100001 (dxy 0.070 genome max, ~0 Fst, 12 SNPs) = assembly/collapsed-repeat artifact, excluded, not a hotspot.
tags: [gwas, loco, sensitivity, skat, burden, set-test, pixy, dxy, fst, co-localization, tierb, highdxy, near-clone, standing-variation, rhodotorula, gemma, mcverification]
```

### 2026-08-15-color-phenotype-space-gwas-tierdeg
```yaml
name: 2026-08-15-color-phenotype-space GWAS Tier D/E/G (locus->gene mapping + ABF fine-mapping + prior-locus replication)
question: Which genes do the Tier-A/BSLMM FDR-significant loci map to (Tier D)? Can anchors be resolved into credible sets with effect-size bounds (Tier E)? Does the prior lab's growth-rate locus chr13:13_30149 (p=1.68e-11) replicate in our color/copper panel (Tier G)?
input: results/gwas/fdr/*_fdr05_sig.txt; results/gwas/tierA_summary/gwas_*_assoc.csv.gz (p_wald column); results/gwas/tierC_summary/tierc_bslmm_summary.csv; gene index (mRNA-derived, 6,799 genes) results/gwas/tierD/scripts/gene_index.json; prior-GWAS file /bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/GWAS/SNP_GWAS_T37C_linear.tsv
scripts:
  - scripts/annotate_gwas_loci.py       # FDR-sig + BSLMM loci -> gene overlap/nearest; 250kb positional clump (single lead per chr)
  - scripts/finemap_credible_sets.py    # Wakefield ABF in z-space (NCP prior SD=0.2, logsumexp, candidate p<1e-3) -> 90/95/99% credible sets
  - scripts/make_tierde_figure.py       # Tier D gene-map (scaffold_10 anchor region) + Tier E credible-set resolution figure
outputs:
  - results/gwas/tierD/tierD_fdr_snps_annotated.csv.gz (12,348 rows)
  - results/gwas/tierD/tierD_independent_loci.csv.gz (5,286 loci; 250kb clump caveat)
  - results/gwas/tierD/tierD_bslmm_loci_annotated.csv (8 top-PIP loci)
  - results/gwas/tierD/scripts/gene_index.json
  - results/gwas/tierE/tierE_credible_sets.csv (42 sets x 14 anchors)
  - results/gwas/figures/tierde_gene_finemap.{png,pdf} (panel A: scaffold_10 chroma+AUC_10 gene map; panel B: 95% CS size vs lead pp, rare-driven flag)
  - GWAS_REPORT.md section 9 (results) + section 10 (draft publication methods, Tier A->G chain)
reproduce: PY .pixi/envs/default/bin/python; $PY scripts/annotate_gwas_loci.py --fdr results/gwas/fdr --index results/gwas/tierD/scripts/gene_index.json --outdir results/gwas/tierD; $PY scripts/finemap_credible_sets.py; $PY scripts/make_tierde_figure.py --sig results/gwas/tierD/tierD_fdr_snps_annotated.csv.gz --credset results/gwas/tierE/tierE_credible_sets.csv --out-dir results/gwas/figures
status: complete
key_findings:
  - Tier D: 12,348 FDR-sig SNPs annotated; 5,286 independent loci across 9 traits. Notable genes: chroma scaffold_8 -> telomerase RT (OM429_004009); AUC_10 scaffold_10 -> DBP3 RNA-dependent ATPase (OM429_004640), RNA-pol-I TF, endodeoxyribonuclease; BSLMM chroma scaffold_3 -> methionine aminopeptidase 1 (OM429_001379). Caveat: 250kb clump keeps only the single most-significant SNP per chromosome as lead, so a second distant block on the same small scaffold collapses (chr13's 217 FDR SNPs = ONE block, single lead 13_791853 p=2.4e-7).
  - Tier E: chroma scaffold_10:384905 resolves well (95% CS n=24, lead pp=0.054, beta=1.11+/-0.19 SE, common AF 0.81) - the paper-grade common-variant anchor. Rare-EF loci (AF~0.015) give wide sets (auc10 DBP3 CS n=67, rare_driven) or fully-resolved singletons (resil scaffold_13 95% CS n=1, pp=0.615). ABF must be in z-space (trait-scale priors break across ~1e5 unit spans); candidate filter p<1e-3 required.
  - Tier G: prior growth-rate locus chr13:13_30149 REPLICATES in our AUC_10 via nearest proxy scaffold_13_30134 (15 bp): p_wald=4.03e-6, FDR-sig, beta=804,778, af=0.015, inside the single chr13 rare-haplotype block. Gene at locus = OM429_005439, a hypothetical protein with no GO/InterPro/PFAM annotation (flanked by Ark1-family Ser/Thr kinase OM429_005441). Other traits null -> replication is growth-phenotype-specific. Causal gene under a p~1e-11 locus is functionally unknown - priority validation target.
tags: [gwas, tierd, annotation, gene-mapping, tierte, finemapping, credible-sets, abf, wakefield, tierg, replication, prior-locus, chr13, AUC_10, DBP3, telomerase, methionine-aminopeptidase, OM429_005439, rare-EF, near-clone, rhodotorula, gemma]
```

### candidate_gene_significance_plots_and_homology
```yaml
name: candidate_gene_significance_plots_and_homology
question: For the 10 GWAS candidate genes (analysis/candidate_gene_alignment/), what do the phenotype-association results look like visually (Manhattan-style), and are the significant missense residues conserved in other Rhodotorula species and S. cerevisiae?
input: analysis/candidate_gene_alignment/results/candidate_{gene,indel}_phenotype_assoc_all10.csv; analysis/gwas/results/gwas/tierA_summary/assoc_csv/gwas_{lab_a,chroma,bright}_assoc.csv.gz; analysis/candidate_gene_alignment/results/<gene>/protein.fasta (reference strain); 11 local Rhodotorula/Cystobasidium proteomes + S. cerevisiae SGD proteome (both external, see script headers for paths)
scripts:
  - scripts/make_candidate_gene_manhattan.py   # per-gene variant significance plot + genome-wide Manhattan with candidate loci highlighted
  - scripts/build_homolog_impact_table.py      # diamond blastp of each gene's significant missense variant against 12 proteomes; maps mutated residue through the local alignment to call conserved/diverged/outside-aligned-region
outputs:
  - figures/{candidate_gene_association_manhattan,genome_wide_manhattan_with_candidates}.{png,pdf}
  - results/candidate_genes/sequences/candidate_genes_{cds,protein}.fa (reference-strain sequences, quick-lookup copy)
  - results/candidate_genes/homology/{blast/*.tsv, homolog_impact_table.{csv,md}, homolog_impact_summary.csv}
  - GWAS.md section 17
reproduce: module load diamond/2.1.7 (NOT piped through `| tail` etc — see .living/conventions.md L-30, re-hit a 4th time building this); python3 scripts/make_candidate_gene_manhattan.py; python3 scripts/build_homolog_impact_table.py (after building diamond dbs, see script/section 17 for db source paths)
status: complete
key_findings:
  - Genome-wide Manhattan overlay confirms the candidate genes are not an arbitrary pick -- several (OM429_000065, OM429_001430/001533, OM429_003729, OM429_005034) sit directly at or beside genome-wide lab_a/chroma/bright peaks.
  - Cross-species homolog scan: every one of the 5 genes with FDR-sig missense hits has a clear same-genome paralog (distinct gene ID) and 61-97%-identity homologs in all 10 other local Rhodotorula/Cystobasidium genomes; S. cerevisiae homologs found for most but not all. Residue-level conservation calls are noisy (single-best-hit diamond, not a curated MSA) -- many N/M-aligned denominators are small because local alignments don't extend to N/C-terminal variant positions (e.g. OM429_005034's variants at residues 5/12).
  - OM429_002663 p.Asn4Lys (LYS1, replicated color association) is conserved (Asn) in 3/4 aligned cross-genus homologs, suggesting the phenotype-associated Lys4 allele is the derived state in this panel.
tags: [gwas, candidate-gene, manhattan, visualization, homology, diamond, blastp, cross-species, conservation, saccharomyces, rhodotorula, orthology, complete]
```

### gwas
```yaml
name: gwas
question: Can Rhodotorula mucilaginosa color and copper-response phenotypes be mapped with GWAS (GEMMA kinship-only LMM), after an audited strain-name reconciliation and per-strain ploidy validation?
input: data/raw/genotypes/RmucY2510_v2/ (versioned symlinks), data/metadata/Copper.Strain_info.csv
scripts:
  - scripts/ingest_gwas_vcf.sh                          # repo-root: versioned VCF ingestion
  - analysis/gwas/scripts/reconcile_strains.py
  - analysis/gwas/scripts/check_ploidy.py
  - analysis/gwas/scripts/diff_strain_state.py
  - analysis/gwas/scripts/check_grm_conditioning.py
  - analysis/gwas/scripts/cull_near_clones.py
  - analysis/gwas/scripts/rebuild_full_genotypes_and_tiera.sh    # supersedes rebuild_tiers_abc.sh's Tier A scope
  - analysis/gwas/scripts/resume_tiera_gemma.sh
  - analysis/gwas/scripts/build_gwas_phenotypes.py
  - analysis/gwas/scripts/summarize_tiera.py
  - analysis/gwas/scripts/convert_assoc_for_tierb.py
  - analysis/gwas/scripts/check_population_confounding.py
  - analysis/gwas/scripts/run_tierb.sh + tierb_set_tests.py
  - analysis/gwas/scripts/run_tierc_bslmm.sh + summarize_tierc.py
  - analysis/gwas/scripts/run_loco.sh + run_loco_shared.sh + merge_loco.py
  - analysis/gwas/scripts/add_copper_v0151_trait.py     # reconciles + merges copper-heavy-metal-screen-v0.15.1's mean_auc_rate into the fam-order phenotype table
  - analysis/gwas/scripts/compare_copper_v0151_trait.py # correlates cu_doseauc_v0151 vs. existing copper traits
  - analysis/gwas/scripts/summarize_cu_doseauc_v0151_gemma.py  # FDR + nearest-gene + cross-trait replication table for the new trait's GEMMA scan
  - analysis/gwas/scripts/run_population_vs_locus_cu_v0151.sh  # drives check_population_vs_locus.py once per locus (5 loci) and merges results
  - analysis/gwas/scripts/qc_scaffold16_investigation.py  # per-scaffold VCF QC + focal-SNP LD + pixy Fst summary (D-26)
  - analysis/gwas/scripts/finemap_candidate_genes.py  # reused (unmodified) for scaffold_13:810026, D-27
  - analysis/gwas/scripts/rare_variant_carrier_permutation_test.py  # population-stratified exact carrier permutation, D-29
outputs:
  - analysis/gwas/results/strain_reconciliation/strain_match_table.reviewed.csv
  - analysis/gwas/results/ploidy_check/ploidy_flags.csv
  - analysis/gwas/results/strain_state_diff/state_diff_report.json
  - analysis/gwas/results/gwas/near_clone_culling/culled_keep.txt
  - analysis/gwas/results/gwas/grm_conditioning/{grm_diagnostic.json,grm_diagnostic_culled.json,rebuilt_kinship/,rebuilt_kinship_culled/}
  - analysis/gwas/results/gwas/tierA_summary/{tiera_summary_gwas.csv,tiera_summary_gwasc.csv,assoc_csv/,fdr/}
  - analysis/gwas/results/gwas/tierA_summary/population_confounding_{gwas,gwasc}.csv
  - analysis/gwas/results/gwas/tierB/tierb_settests_{gwas,gwasc}.csv
  - analysis/gwas/results/gwas/tierC_summary/{tierc_bslmm_summary.csv,pip/}
  - analysis/gwas/results/gwas/loco/output/loco_merged_{gwas,gwasc}.csv
  - analysis/gwas/results/strain_reconciliation/copper_v0151_match_table.csv
  - analysis/gwas/results/gwas/tierA_summary/gemma_output/gwas_cu_doseauc_v0151.assoc.txt
  - analysis/gwas/results/gwas/tierA_summary/cu_doseauc_v0151_{correlations,top_hits_annotated}.csv
  - analysis/gwas/results/gwas/population_vs_locus/population_vs_locus_cu_doseauc_v0151.csv
  - analysis/gwas/results/gwas/population_vs_locus/population_vs_locus_scaffold13_{resilience_30,AUC_30}.csv
  - analysis/gwas/results/gwas/scaffold16_investigation/{per_scaffold_qc,focal_snp_ld,per_scaffold_fst}.csv
  - analysis/gwas/results/gwas/tierE/scaffold13_{credible_sets,window_genes}.csv
  - analysis/gwas/results/gwas/rare_variant_validation/{cu_doseauc_v0151_scaffold_9_704260,resilience_30_scaffold_13_810026,AUC_30_scaffold_13_810026}.csv
reproduce: bash analysis/gwas/run.sh (human-review gates at reconciliation and ploidy steps; SLURM submission steps for the kinship rebuild, Tier A/B/C, and LOCO documented inline in each script's header). cu_doseauc_v0151 (§18, added 2026-08-27) reproduced separately via add_copper_v0151_trait.py -> compare_copper_v0151_trait.py -> the same run_tiera_gemma.sh GEMMA invocation pattern for one trait -> summarize_cu_doseauc_v0151_gemma.py -> run_population_vs_locus_cu_v0151.sh.
status: in-progress (Tier A [corrected, full unpruned SNP set], Tier B, Tier C, and LOCO all complete on both the 213-strain and 182-strain near-clone-culled panels; pixy reused from the prior 201-strain run as a documented approximation, not recomputed). cu_doseauc_v0151 (new copper-v0.15.1 dose-response-AUC trait) added, GEMMA-scanned, and population-vs-locus validated 2026-08-27 -- all 5 candidate loci verdict likely_population_artifact, see GWAS.md §18. Tier D/E/G (fine-mapping) run on the resilience_30/AUC_30 scaffold_13:810026 anchor 2026-08-27 -- it too fails population-vs-locus validation (D-27, rare-variant blind spot); scaffold_16's recurring-hit pattern investigated and explained (D-26, not a supergene/artifact); chroma/lab_L/a/b instability resolved (D-28). Orthogonal rare-variant carrier-permutation test built and run 2026-08-27 (D-29): scaffold_13:810026 passes decisively (p=0.001-0.003, best-supported rare hit in the project), scaffold_9:704260 does not (p=0.199, unconfirmed on two independent grounds).
key_findings:
  - Audited reconciliation accepted 213/308 phenotype strains (12 more than the prior informal our200.txt match), all 62 fuzzy candidates hand-reviewed and rejected as coincidental string similarity, 0 removed vs. prior run.
  - Ploidy check confirms haploid GT encoding genome-wide (0 het in every strain); 46/213 strains flag watch_contamination (depth >2.5x panel median) and 14 watch (1.5-2.5x) -- kept per user decision, informational only for now.
  - GRM conditioning diagnostic flagged singular_risk (condition_number ~2e11) on both panels, but investigation traced it to a small near-clone triplet (DBVPG_5757/5758/5759) plus one benign GRM-centering null eigenvalue -- consistent with D-9's precedent, not new pathology.
  - Population-structure check (GRM PC1=33%/PC2=10% var explained, cleanly separating 6 subpopulations, Fst~0.45) confirms real structure but kinship-only LMM is already the correct correction (D-9); no additional PC/population covariate is warranted or would even run (collinearity with near-clone eigenvectors).
  - Near-clone IBS0 culling reconstructed as reusable code (never saved originally), validated at 93.6% membership overlap against the prior run's known 173-strain result; 182/213 kept on the new panel.
  - **resilience_30/AUC_30's scaffold_13:810026 anchor** replicates across every panel and method tested (prior 201-strain run, full-213, full-182-culled, LOCO with chr13 excluded from kinship) at p ranging 6.4e-9 to 8.7e-10, independently corroborated at lower confidence by Tier C BSLMM's nearby top locus (scaffold_13:793374, ~17kb away), and is not in a high-dxy window (Tier B). **REVISED 2026-08-27 (D-27): population-vs-locus tested for the first time and comes back `likely_population_artifact`** -- af=0.014 (~3 alt carriers total), 0/6 populations have enough carriers to test, so none of its cross-panel/cross-method "replications" are actually independent evidence (same handful of carrier strains drive all of them). No longer citable as a confirmed §13-verdict finding; revealed a structural rare-variant blind spot in the disambiguation battery itself (can fail rare hits, cannot validate them). Fine-mapped anyway for reference: nearest gene `OM429_005716`, unannotated. **FURTHER UPDATE (D-29): a purpose-built orthogonal test (population-stratified exact carrier permutation) DOES support this locus** -- p=0.0013 (resilience_30) / 0.0025 (AUC_30), carriers span 2 real populations (pop1, pop6), not confounded by the "same carrier strains every time" critique above. Now the best-supported rare hit in the project, though still short of a full §13 `likely_real` verdict (can't rule out within-population relatedness).
  - **scaffold_16 investigated (D-26)**: recurs as a top hit for many unrelated traits, but is NOT one non-recombining haploblock or an assembly artifact -- per-scaffold QC (depth/density/rare-fraction) is unremarkable, LD between the recurring hit positions decays normally (only truly-close hits, e.g. cu_doseauc_v0151's own 39kb cluster, are in strong LD), and mean Fst is only modestly elevated (5th of 23 scaffolds). Several independent, moderately-differentiated loci, not one phenomenon.
  - **chroma/lab_L/a/b instability resolved (D-28)**: lab_L/a/b were already `likely_real` per the existing D-17/D-18 battery (a documentation cross-link gap, not open science). chroma genuinely has no single-locus signal -- 0 FDR-sig SNPs in most marker-set variants AND 0 SNPs above any BSLMM posterior-inclusion threshold (PGE=9.3% of PVE=0.213, i.e. nearly fully polygenic) -- its shifting "top hit" is single-SNP-scan noise, not a bug or hidden confound.
  - chroma's top hit is unstable across every SNP-set/panel choice tried (5 different top loci across 5 variants) and is NOT population-confounded -- cause still open, flagged not to be cited as a replicated finding.
  - Tier B: no set-level (burden/SKAT) signal survives BH-FDR in either panel, matching the original run's pattern exactly.
  - Tier C BSLMM architecture looks substantively more polygenic here than the original run (no locus reaches PIP>=0.5 in any of the 5 traits tested, vs. the original's several near-oligogenic PIP~1.0 loci) -- not yet explained, candidate causes include the larger unpruned marker set diluting PIP and/or the different near-clone structure.
  - **cu_doseauc_v0151 (dose-response AUC from copper-heavy-metal-screen-v0.15.1) is only moderately correlated with existing copper traits (Spearman rho 0.04-0.46)** and its GEMMA scan surfaces the single strongest single-SNP signal in the whole GWAS port so far -- scaffold_9:704260 (p=5.9e-14, FDR q=1.7e-9), nominally significant in all 12 existing traits but never before FDR-significant. **BUT all 5 candidate loci (incl. scaffold_9:704260 and the recurring scaffold_16 region) FAIL population-vs-locus validation, verdict `likely_population_artifact`** (D-25): 4 rare SNPs (af 0.015-0.030) are essentially private to 0-1 of 6 populations (too rare to ever test); scaffold_16:455499 is common but near-completely population-fixed (fst_proxy=0.84, only 1/6 populations testable). None citeable as real findings -- the cleaner trait mainly surfaced classic population-stratification false positives more clearly than the noisier existing traits did.
tags: [gwas, gemma, rhodotorula, mucilaginosa, copper, color, kinship, ploidy, strain-reconciliation, near-clone-culling, tierb, tierc, bslmm, loco, population-structure, ported, in-progress]
```

### candidate_gene_alignment
```yaml
name: candidate_gene_alignment
question: For candidate genes identified during the GWAS port (2 confirmed carotenoid pathway genes + 8 GWAS-locus nearest genes), what do the actual per-strain DNA/protein sequence changes look like across the 213-strain panel -- moving from "this SNP is statistically associated" to "here is the specific allelic/amino-acid change", and does any candidate-gene coding variation (SNP or indel) associate with the 6 color traits under the population-aware battery?
input: data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.SNP.combined_selected.vcf.gz, .../RmucY2510_v2.All.INDEL.combined_selected.vcf.gz, .../genome/Rhodotorula_mucilaginosa_NRRL_Y-2510.{gff3.gz,scaffolds.fa}, pre-built snpEff database RmucNRRLY2510
scripts:
  - scripts/extract_gene_sequences.py       # GFF3 + VCF + genome -> per-strain spliced CDS/protein FASTA + polymorphic-positions CSV
  - scripts/screen_indels.sh                # read-only screen of the separate INDEL VCF against CDS+/-2kb (fixed 8/26 to fail loudly; pilot's silent-0 was wrong)
  - scripts/check_gene_coding_indels.py     # CDS-exon-exact INDEL screen + segregating-in-panel count (results/gene_coding_indel_screen.csv)
  - scripts/run_extend_to_gwas_genes.sh     # batch driver: 8 GWAS-locus genes through steps 1-5 (abs-path toolchain + snpEff dataDir fixes)
  - scripts/build_variant_table.py          # snpEff ANN parsing + per-strain genotype/population/phenotype/lead-SNP/Tier-A-p-value join
  - scripts/render_alignment_image.py       # static color-block alignment PNG per gene
  - scripts/associate_variants_with_phenotype.py  # population-aware battery (within-pop meta + covariate + FDR) for coding SNP variants x 6 color traits
  - scripts/associate_indels_with_phenotype.py    # same battery for segregating in-CDS INDEL genotypes (INDEL VCF direct)
outputs:
  - results/<gene_id>/{dna.fasta, protein.fasta, dna_polymorphic_positions.csv, indel_screen.csv, variant_table.csv, variant_table_strain_context.csv, alignment.png} for all 10 genes
  - results/gene_coding_indel_screen.csv
  - results/candidate_gene_phenotype_assoc_all10.csv (+_within_pop_details.csv), results/candidate_indel_phenotype_assoc_all10.csv (+_within_pop_details.csv)
  - results/snpeff_pilot/, results/PROVENANCE.json
reproduce: bash scripts/run_extend_to_gwas_genes.sh (8 GWAS-locus genes; pilot run gene-by-gene per CANDIDATE_GENE_ALIGNMENT.md); then run the two associate_* scripts; visualization/homology follow-up lives in ../gwas/scripts/{make_candidate_gene_manhattan.py,build_homolog_impact_table.py} (see gwas manifest entry below)
status: complete for all 10 target genes (2 carotenoid + 8 GWAS-locus). Phenotype-association step (SNP + indel) complete. Significance plots, reference CDS/protein FASTAs, and cross-species homolog conservation table added under analysis/gwas/results/candidate_genes/ + GWAS.md section 17 (2026-08-26). Interactive viewer still deferred.
key_findings:
  - CORRECTION (8/26): pilot's screen_indels.sh silently reported 0 indels (bcftools not on PATH; pipe+`|| true` swallowed failure). Real INDEL VCF has 362-821 records in every candidate gene's CDS+/-2kb, and all 10 genes have >=2 segregating indels INSIDE their CDS exons -> the SNP-only sequence premise does NOT hold for any target gene; indel genotypes must be tested separately (now done).
  - Association battery (within-pop replication + meta + covariate + FDR) across 10 genes x 6 color traits: 39 FDR-sig SNP-coding tests / 32 FDR-sig indel tests.
  - Strongest finding: OM429_000065 (sat-locus gene) c.839C>T p.Ala280Val (lead SNP 208569) + a +6bp in-frame GAGCGG-repeat insertion at 208398 in perfect LD with it -- both associate with lab_a/chroma/sat (meta_p~9e-11, partial R2 0.09-0.10) replicated 3/3 testable pops. The insertion is invisible to SNP-only pipelines.
  - Other replicated coding hits: OM429_001521 p.Gln245His (bright), OM429_005034 p.Tyr12His (lab_a/chroma/sat), OM429_002663 p.Asn4Lys (lab_a/chroma/sat), OM429_001533 p.Tyr427Phe.
  - OM429_003333 splice_donor (HIGH) is monomorphic in the 213 panel (0/213 alt) -> cannot be tested / cannot drive within-panel color variation.
  - Neither carotenoid gene (OM429_003333/003336) shows an FDR-sig + replicated color association (SNP or indel), consistent with no co-localization with validated color loci.
tags: [gwas, candidate-gene, sequence-alignment, snpeff, carotenoid, pigment, rhodotorula, mucilaginosa, indels, association, population-structure, complete]
```
