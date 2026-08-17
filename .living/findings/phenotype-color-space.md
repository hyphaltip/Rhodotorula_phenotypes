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

8. **[robust] GWAS feasibility (idea 11): R. mucilaginosa has 178 effective
   independent haplotypes (200 tree tips, redundancy ratio 1.12 from 22
   near-clone strains), so only LARGE-EFFECT loci are detectable — min per-SNP
   R² ≈0.16 genome-wide (5e-8) / 0.10 candidate (1e-4) at 80% power, n=202.**
   A detectable allele moves the phenotype ≈0.71 SD (≈22–25% of additive genetic
   variance assuming h²=ICC upper bound); colony size / chroma / pace / 
   heterogeneity are mappable, copper-slope is NOT (ICC 0.19, β≈0.008 trait
   units). Direct power check (noncentral χ²) matches the table.

## (idea 11 continued) Variant resources already exist — no de-novo calling needed
9. **A completed population-genomics + GWAS project for these strains exists**
   at `/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/`
   (reference NRRL Y-2510): `vcf/RmucY2510_v2.All.SNP.combined_selected.vcf.gz`,
   728,581 SNP sites × 422 strains, **haploid** GT format (no het, no HWE
   filtering), GATK-hard-filtered (QD<2, MQ<40, SOR>3, FS>60, PASS-only).
   **201 of our 278 phenotyped strains (all R. mucilaginosa) carry genotypes**;
   200/201 complete for all 4 GWAS traits. So the variant pipeline blocker is
   resolved entirely by reuse; our contribution is the colour-trait phenotype
   side (never tested in their 218-strain growth-rate GWAS, which found only a
   handful of loci — consistent with the idea-11 effective-n power limit).
10. **QC signal for the GWAS run:** heterogeneous, high-depth scaffolds
    (scaffold_6 max 25,373×, scaffold_18 10,702×, scaffold_21 mean 550×)
    look like aneuploidy/collapsed-repeats and should be depth-QA'd; strain
    heterozygosity varies 0.014→10.5 (they genotyped many near-clones), so
    kinship-based LMM + LD-pruned relatedness is the right model; 6 population
    assignments cover our 201 strains for optional structure conditioning.

## Caveats

- Onset-time results are threshold-sensitive; use chroma>7 (basal noise ~1.5–2).
- Mediation fraction is imprecisely estimated (CI [0.49, 2.69]); "full
  mediation" is directionally robust, not exactly 1.0.
- Do not interpret "no direct Cu pigment effect" as "no evolutionary response":
  the residual a* direct channel and GxE indicate a second, pigment-specific
  mechanism.
- Environment signal is small (η²≈0.06) and depends on coarse environment
  categories; finer habitat fields are not available in the strain table.

## (GWAS + pixy run) Population-structure correction numerically unsupported; scaffold_1 diversity is strongly population-structured

11. **GEMMA LMM + structure covariates degenerates in this near-clonal panel.**
    Adding 10 genotype-PCs (or 3, or population dummies) to the kinship-corrected
    LMM collapses the model: `pve = 0.99997`, `se(pve) = NaN`,
    `GSL: matrix is singular`; ~55 min/trait vs ~3 min kinship-only. Root cause is
    structural — the GRM from 22 near-clone strains is singular, and top genotype
    PCs are nearly collinear with its eigenvectors, so any fixed covariate is
    redundant and the joint model breaks. The **kinship-only LMM is the valid
    analysis**: PVE = 0.180 ± 0.077 (sane), lambdas 0.357–0.638, top size hit
    scaffold_16:121473 p_wald = 6.2e-11 (13 independent loci after clumping).
    This empirically documents that explicit structure covariates are unsupported
    by the near-clonal relatedness, not that structure is absent.
12. **Scaffold_1 diversity (pixy, 201 strains, 6 populations, 22×100 kb windows) is
    strongly population-structured.** Mean π varies >600× across populations:
    pop5 = 0.044, pop6 = 0.016, pop3 = 8.9e-3 vs pop1 = 5.1e-5, pop4 = 7.0e-5
    (near-monomorphic). Mean Hudson Fst across the 15 pop-pairs = **0.494**;
    most divergent pop pairs are pop4–pop2 (0.876) and pop4–pop1 (0.773), least
    pop3–pop6 (0.251). dxy ranges 2.9e-5 → 0.053. Pop architecture aligns with the
    near-clone continuum and explains the conservative GWAS lambda (<1): the
    kinship over-absorbs variance concentrated in a few divergent clones/populations.
    Full-genome pixy (all 23 scaffolds) is the natural next step to map where
    trait-relevant divergence lives.

13. **Next-gen GWAS (expert-guided, 12 clone-mean phenotypes; kinship-only LMM):**
    Color (chroma/sat/bright) multi-trait Fisher p = **7.5e-16 at scaffold_10:384905**
    (chroma p=2.4e-8 co-localizes sat/bright) — a coherent pigmentation locus.
    Copper block multi-trait Fisher p = **5.9e-12 at scaffold_16:421208**; AUC_30
    scaffold_4:508257 repeats on all-201 and culled-173. Lambdas mostly <1
    (conservative; near-clone structure), but clone_mean_area λ=2.12 and AUC_0
    λ=2.74 are inflated; IC50_est-culled λ=0.19 with 673 SNPs <5e-8 is structure
    inflation, not signal. Culled-173 (IBS0<0.005, 173 strains) shifts the top color
    hit to scaffold_3:570085 — rank instability flagged for LOCO.
14. **Tier-C BSLMM architecture (all-201, 100k MCMC):** AUC_10 (copper tolerance)
    is **near-oligogenic**: PVE 0.40, PGE 0.96, only ~3 active LD-clusters, top loci
    scaffold_5:533603 + scaffold_2:939277 (PIP 1.0) — a handful of variants carry
    essentially all heritability. chroma is moderately heritable (PVE 0.24), partly
    sparse (PGE 0.39), scaffold_3-enriched. resilience_30 sparse with scaffold_16:563722
    (PIP 0.95) matching the Tier-A multi-copper block; clone_mean_area weakly heritable
    (PVE 0.07). BSLMM clusters converge with Tier-A signals, giving trait-specific
    candidate loci for follow-up annotation.
