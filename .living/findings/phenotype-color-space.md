# Finding: Copper acts mostly through growth; phenotype signal lives in heterogeneity/shape, not lightness

- **Date**: 2026-08-15
- **Source**: `analysis/ideas/2026-08-15-color-phenotype-space/` (8 scripts, ideas 01–08; summary in FINDINGS.md) + idea 09 phylogeny script and tree asset `data/raw/rhodotorula-phyling-protein-tree/`
- **Topic**: phenotype-color-space, copper-pigment-mediation

## Result

0. **[robust] sd(log₁₀ Asat) is a size/floor-coupled metric — the v1 "dispersion
   widens with Cu" (idea 01) claim was largely confounded.**
   sd(log Asat) correlates ~−0.63 with colony size; within-strain (paired) raw
   Cu slopes are positive overall (median +0.0033/Mm, Wilcoxon p≈5e-16) but
   collapse to ~0 when colony size is controlled (+0.0002/Mm, p=0.46). Only
   lineage-specific exceptions survive: R. taiwanensis genuinely widens
   (+0.014/Mm, p=0.03) and R. paludigena narrows (−0.025/Mm, p=0.001);
   species differ in within-strain slope (Kruskal H=40.5, p≈1e-6).
1. **[robust] Cu suppresses pigmentation largely via growth arrest, not a direct
   pigment block.** At 90–110 h well level, total Cu effect on log(chroma)
   (Δb≈−0.022) vanishes once colony area + species are controlled
   (direct b≈+0.004; mediation fraction ≈1.2, bootstrap CI [0.49, 2.69]).
   The a* (redness) channel is different: total −0.32 → direct −0.13, i.e.
   ~60% growth-mediated but a real ~40% growth-independent direct component.
2. **[robust] The strain phenotype space is a continuum, and species/interactions
   are distributed across it, not localized.** Varimax-factorized axes form a
   near-block-diagonal tripartition: (F1) color level — chroma/saturation/L*,
   (F2) intra-colony heterogeneity — texture entropy + L*CoeffVar, (F6/F8)
   shape/area. UMAP + HDBSCAN finds no discrete species clusters (all clusters
   mucilaginosa-dominated). GxE on log(chroma) within R. mucilaginosa is real
   but modest (F=1.66).
3. **[robust] Lightness/intensity carry almost no species or stress information;
   within-colony heterogeneity and shape carry the signal others miss.**
   - Species MI: best single feature ≈0.19 bit of H(species)=3.0 bit;
     L*Median and intensity are bottom (≈0.07 bit); top features are
     b*CoeffVar, shape (Area/Perimeter/Solidity).
   - Environment (marsh_tidalflat/soil/plant/food), species-blocked
     permutation: only b*CoeffVar (p≈0.006) and L*Median (p≈0.027) associate;
     mean pigmentation (a*, chroma) and size do NOT.
   - Repeatability: shape ICC_strain 0.86–0.93 and chroma 0.82 are heritable-
     style axes, but a*CoeffVar (0.19) and b*CoeffVar (~0) are colony-noise —
     the "stress heterogeneity" in a*/b* variance is not genotypic.
4. **Colony size only weakly gates the timing of pigmentation.** Spearman
   (time-to-size, onset-time) ≈ 0.30–0.42 across chroma thresholds 2→10;
   and pigment accrues sub-linearly vs area after onset, with strong
   species-specific exponents (pigment~area^β, β 0.19–0.53).
5. **[robust] Strain phenotypes DO carry phylogenetic signal (idea 09, Mantel
   permutation on the PHYling protein tree), but mostly BETWEEN species, not a
   gradation within R. mucilaginosa.** All-strains scope: l10med_fixed (colony
   size) r=0.43 (p=0.001), partial_slope_sd_cu (within-strain heterogeneity
   broadening) r=0.31 (p=0.001), intercept_logchroma (baseline chroma) r=0.19
   (p=0.001), slope_logchroma_per_mM (Cu sensitivity) r=0.10 (p=0.003);
   pace_loglog not structured (r=0.06, p=0.07). Within R. mucilaginosa only
   baseline chroma (r=0.13, p=0.002) and colony size (p=0.027) survive.
   Reconciles the idea-05 "continuum" with real phylogenetic structure: it is a
   within-species continuum on top of between-species differences.
6. **[methodological, important] Blomberg's K is degenerate (~1e-7 = "no signal")
   on near-comet trees with a giant polytomy of near-duplicate genomes.**
   Here 167/541 edges ≤1e-7 (near-duplicate R. mucilaginosa clade); a BM power
   check reproduced K≈2.25 on simulated data on this same tree (should be ~1),
   and likelihood lambda was numerically unstable (geiger white-vs-lambda lnL
   internally inconsistent). Use rank-based Mantel permutation on such
   topologies, never raw K.
7. **[robust] Within-species (among-strain) variation is the DOMINANT scale of
   variation for every trait (idea 10).** Exact SS decomposition (11 species,
   n≥3): fraction of variance WITHIN species = Cu-sensitivity slope 87%,
   baseline chroma 92%, colony size 62%, within-strain heterogeneity 67%,
   pigment pace 84% (ANOVA F 2.6–17.5; n≈293–301). Colony size is most
   species-structured (38% between-species), chroma & Cu-slope least
   (~8–13%). Reconciles idea 09 (detectable between-species Mantel signal in
   the smaller component) with idea 05 (within-species continuum).

## Caveats

- Onset-time results are threshold-sensitive; use chroma>7 (basal noise ~1.5–2).
- Mediation fraction is imprecisely estimated (CI [0.49, 2.69]); "full
  mediation" is directionally robust, not exactly 1.0.
- Do not interpret "no direct Cu pigment effect" as "no evolutionary response":
  the residual a* direct channel and GxE indicate a second, pigment-specific
  mechanism.
- Environment signal is small (η²≈0.06) and depends on coarse environment
  categories; finer habitat fields are not available in the strain table.
