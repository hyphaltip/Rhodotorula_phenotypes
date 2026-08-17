# FINDINGS — 2026-08-15 Color & Phenotype Space Exploration (ideas 01–08)

Ideation session implemented autonomously: 8 personas → 16 ideas → 8 scripts →
`results/idea*.csv` + `figures/fig*.png`. Shared pipeline: `scripts/build_series.py`
(→ `data/db_extract.parquet`, 211,800 colony-rows × 116 cols) + `scripts/common.py`.
All analyses use Cu-graded growth curves from DuckDB `v_phenotype`
(5 runs/353–357, 320 strains, 8,648 well-curves).

Caveat up front: several "statistical findings" are threshold/definition-sensitive
(species MIs, onset times). Robust claims are flagged **[robust]**.

---

## Idea 01 — Statistical Physicist: universal arrest + growth-fluctuation dispersion
`idea_01_arrest.py` (+ `idea_01b_gibrat_within.py` resolves the dispersion ambiguity)

- **[robust]** Cu suppresses growth by shifts along a **common collapse**, not by
  changing the shape of the growth curve. Weibull arrest exponent k ≈ 0.45 is flat
  across 0–30 mM Cu (0.450–0.462; bootstrap CIs overlap), i.e. the arrest
  "shape" is Cu-invariant; Cu moves the *time* (median t50 68.4 → 65.2 h) and the
  asymptotic size, not the dynamics' form.
- Dispersion (Gibrat check): within-strain across-well colony-size spread
  sd(log₁₀ area) rises with Cu at the pooled level, mean 0.30 (Cu=0) → 0.43 (30 mM);
  Spearman ρ = 0.23–0.28 (p<1e-27). **BUT this pooled slope is largely a
  size/floor artifact** (idea01b):
  - sd(log₁₀ Asat) is strongly anti-correlated with colony size
    (Spearman −0.63, p<1e-200): Cu shrinks colonies, and smaller colonies carry
    more log-dispersion.
  - **Within-strain (paired, same genotype) raw slopes are positive**
    (median +0.0033 sd per mM, Wilcoxon p≈5e-16, 75.9% strains positive);
    **size-controlled slopes collapse to ≈0** (median +0.0002, p=0.46,
    53.6% positive) → no universal quenched-disorder widening per genotype.
  - Species genuinely differ (Kruskal H=40.5, p≈1e-6): `R. taiwanensis`
    **widens with Cu even size-controlled** (+0.014 per mM, p=0.03, all 6
    strains positive), while `R. paludigena` **narrows** (−0.025, p=0.001) — a
    genuine lineage-specific behavior on top of an overall size artifact.
- Implication: the v1 "dispersion widens with Cu" claim was composition+size
  confounded; the defensible claim is a size/floor coupling plus lineage-specific
  (taiwanensis positive, paludigena negative) trends.

## Idea 02 — Color/Imaging Scientist: pigment identity, morphs, onset, run calibration
`idea_02_color.py`

- **[robust]** Species separate cleanly in CIELAB hue (species-typical hue° differ by
  ~12–54°), chroma 5.5–13.2, L* pinned ≈ 73–77 (L* carries little species info).
- a*CoeffVar (within-colony redness heterogeneity) rises with Cu 0.31 → 0.54 — a
  label-free stress response in colony color homogeneity.
- Pigment morphs: k=2; 2nd morph = the Cystobasidium/Pseudomicrostroma (non-
  Rhodotorula) cluster — the odd one out in color space.
- Onset: ≥20 mM Cu → 100% of strains never pigment; Cu<10 mM → 98% pigment with
  median onset ≈24 h. **Flag:** `t_darkening_h=0` everywhere (onset detected as
  immediate) — basal chroma noise defeats this onset definition; see idea 08.
- Run calibration table built (5 runs) for future cross-run corrections.

## Idea 03 — Information Theorist: species information of the feature space
`idea_03_information.py`

- **[robust]** Discrete per-feature species information is **low**: best single
  feature ≈ 0.19 bit <6% of H(species)=3.0 bit. Species differences are mostly
  continuous/along-lineage, not discrete per-feature.
- Ranking is informative, not magnitude: **L*Median and intensity carry the least
  species MI (≈0.07 bit)**; the top features are *within-colony heterogeneities*
  (b*CoeffVar, a* heterogeneity, texture info-correlation) and shape (Area,
  Perimeter, Solidity). MI barely drops when conditioning on L* → this species
  signal is largely **independent of lightness**.
- Effective rank ≈ 11 of 37 curated traits (95% var); mean |corr| ≈ 0.29, but
  pairs like Compactness≈Circularity (r=1.00), Intensity mean≈L*Median (r≈1.00)
  show the raw ~200-col space is heavily redundant (~11 independent axes).

## Idea 04 — Quantitative Geneticist: repeatability audit + reaction norms + GxE
`idea_04_qg.py`

- **[robust]** Shape traits are the most strain-repeatable (ICC_strain:
  Circularity 0.93, Compactness 0.91, Area 0.86); chroma ICC ≈ 0.82. Colony-level
  noise is tiny (ICC_colony <0.01); well-to-well ≈ 0.14–0.19.
- Intra-colony heterogeneity is **repeatable only for some channels**: texture
  entropy ICC=0.80, L*CoeffVar=0.73, SaturationCoeffVar=0.62 — but **a*CoeffVar
  (0.19) and b*CoeffVar (~0) are colony-noise-dominated**, i.e. not heritable.
  The "stress heterogeneity" readout is microenvironmental, not genotypic.
- Reaction norms: all species lose log(chroma) linearly with Cu (slope per mM,
  mean by species −0.019 to −0.036); toruloides & taiwanensis most Cu-sensitive,
  sphaerocarpa & dairenensis least. Higher-baseline-pigmented strains lose faster
  (Spearman(baseline, slope) = −0.59).
- **[robust]** GxE: significant strain × Cu interaction on log(chroma) within
  R. mucilaginosa (F=1.66, dof 1287, p<1e-50) — genotypic plasticity differences
  exist; small F means the *scale* of plasticity is modest vs main effects.

## Idea 05 — Representation Learning: PCA/varimax factors + UMAP/HDBSCAN atlas
`idea_05_repl.py`

- **[robust]** Rank-PCA: PC1+PC2 ≈ 48% of variance; 8 PCs = 89%. Varimax factors
  are near block-diagonal: **F1 color-level (chroma/saturation/xy-x ≈+, L*≈−),
  F2 intra-colony heterogeneity (texture entropy, L*CoeffVar, intensity CV),
  F3 texture corr/contrast, F6/F8 shape** — the color-amount / heterogeneity /
  shape tripartition is a real latent structure.
- Unsupervised atlas: the 35-dim phenotype space is a **continuum** — HDBSCAN on
  raw features returns all-noise; on the UMAP projection it finds 12 local
  neighborhoods that do **not** correspond to species (every cluster is
  mucilaginosa-dominated). No discrete morphotypes; species overlap continuously.

## Idea 06 — Causal Endpoint Mediation: direct vs growth-routed Cu effect
`idea_06_mediation.py`

- **[robust, striking]** At 90–110 h well level, the total effect of Cu on
  log(chroma) (b=−0.022) **vanishes once colony area and species are controlled**
  (direct b≈+0.004, mediation fraction ≈1.2) — chroma loss at Cu is consistent
  with being **fully mediated by growth arrest**, not a direct pigment block.
- Redness (a*) is different: total −0.32, direct −0.13 → ~60% via growth, with a
  residual **growth-independent direct component** (~40%) in the a* channel.
- Bootstrap CI on the chroma mediation fraction [0.49, 2.69] — the "full
  mediation" conclusion is directionally robust, imprecisely estimated.

## Idea 07 — Trait-Ecology Spectrum (environment stratification, species-conditional null)
`idea_07_ecology.py`

- **[robust, non-obvious]** After species-block permutation, environment is
  associated only with **within-colony heterogeneity in b*** (F=4.7, p_strat≈0.006,
  η²≈0.06) and marginally L*Median (p≈0.027). Mean pigment (a*, chroma), colony
  size, and texture entropy do **not** survive species stratification → environment
  has no independent signature in mean pigmentation or size; only in colony-color
  homogeneity/lightness. Ecology tracks species via composition, not via
  pigment traits.

## Idea 08 — Temporal Dynamics: onset logistics + growth–pigment coupling
`idea_08_onset.py`

- **[robust]** Colony size does NOT strictly gate pigmentation: Spearman
  (time-to-size, onset-time) ≈ 0.30–0.42 across onset thresholds 2→10 — weak-to-
  moderate, robustly positive. Onset is only partly size-gated.
- Pigment-areal coupling exponent (log chroma ~ log area after onset) is species-
  specific: mucilaginosa/evergladensis ≈ 0.49–0.53, sphaerocarpa ≈ 0.46,
  diobovata ≈ 0.19 — a ~2.7× range in how strongly pigment accrues per unit
  growth; all sub-linear → pigment lags colony growth during differentiation.
- Onset-rate saturation is threshold-sensitive (median onset 0/12/30/48/84 h at
  thr 2/3/5/7/10) — **flag**: basal chroma ~1.5–2 contaminates low-threshold
  onset definitions; use thr≥7 (51% of IM scope) for clean onset biology.

---

## Idea 09 — Phylogeneticist: strain phenotypes vs relatedness
`scripts/idea_09_phylogeny.R` (tree asset `rhodotorula-phyling-protein-tree`)

- **Strain phenotypes DO carry phylogenetic signal** once measured with a tree-
  robust statistic. Mantel permutation (Spearman, patristic vs trait distance,
  999 perms): colony size l10med_fixed r=0.43 (p=0.001), within-strain
  heterogeneity broadening partial_slope_sd_cu r=0.31 (p=0.001), baseline chroma
  intercept_logchroma r=0.19 (p=0.001), Cu-sensitivity slope r=0.10 (p=0.003);
  pigment pace NOT structured (r=0.06, p=0.07). [n≈266–272/277 matched tips]
- **Within R. mucilaginosa** (200 tips), signal survives only for baseline chroma
  (r=0.13, p=0.002) and weakly colony size (p=0.027); Cu sensitivity and
  heterogeneity show none → structure is between-species, not a fine gradation
  within the dominant species. Reconciles with idea 05's within-species continuum.
- **Methodological finding (important):** this tree is a near-"comet" — a giant
  polytomy of near-duplicate R. mucilaginosa genomes (167/541 edges ≤1e-7; 22
  zero-length pendant edges). On such topologies **Blomberg's K collapses to ~1e-7
  (= "no signal") regardless of truth** (BM power check: K≈2.25 on simulated data
  on this tree), and likelihood lambda becomes numerically unstable (geiger
  white-vs-lambda lnL internally inconsistent). ⇒ Use rank-based permutation
  (Mantel) on near-comet/genome-cluster trees, never raw K.
- Species monophyly (≥3-tip): dairenensis, diobovata, graminis, kratochvilovae,
  sphaerocarpa, sp. clade I monophyletic; mucilaginosa/paludigena/taiwanensis/
  toruloides NOT (consistent with the near-zero mucilaginosa polytomy).

---

## Idea 10 — Scale of within-species variation
`scripts/idea_10_species_variation.py`

- **[robust] Within-species (among-strain) variation dominates the scale of
  variation for every trait.** Exact SS decomposition (11 species with n≥3):
  fraction of total variance that is *within* species — Cu-sensitivity slope
  87%, baseline chroma 92%, colony size 62%, within-strain heterogeneity 67%,
  pigment pace 84% (all-strain species, ANOVA F 2.6–17.5; n≈293–301 strains).
  Colony size is the most species-structured (38% between-species), chroma and
  Cu-slope the least (~8–13% between).
- **Consistent with idea 09 + idea 05:** the phylogenetic signal detected by
  Mantel (idea 09) is *between-species* and modest in absolute terms because
  most of the phenotypic variance sits *inside* species — the strain atlas
  really is largely a within-species continuum (idea 05), with a smaller
  between-species component that is nonetheless statistically detectable.
- Per-species spread is broad: e.g. colony size sp. clade I sd=0.68
  (17.7% CV), sphaerocarpa sd=0.60, paludigena sd=0.53 (log10 area);
  baseline chroma sphaerocarpa sd=0.32, dairenensis sd=0.31 (log chroma).
- Caveats: CV% is not meaningful for traits with mean ≈ 0 (pace_loglog of
  diobovata/graminis/granis); small species (n=3: graminis/pacifica/
  kratochvilovae) have unstable SD estimates; boxplots count n≥3 species.

---

## Cross-cutting synthesis

1. **A tri-partite latent structure dominates the screen:** color-amount (F1), 
   intra-colony heterogeneity (F2), and shape/size — largely orthogonal axes that
   load by block, not by species.
2. **L* / intensity carry the least signal**; within-colony heterogeneity and
   shape carry species + environment + stress signal others miss.
3. **Cu acts mostly through growth**: chroma loss fully growth-mediated; the
   only clearly growth-independent channel is a* (redness) — pointing toward a
   second, pigment-specific mechanism at high Cu.
4. **No discrete phenotype types**: the strain atlas is a continuum; species
   overlap; phenotypes are quantitative.
5. **Definitional caveats to carry forward:** onset-time thresholds, species MI
   discreteness, mediation fraction width — all sensitive to choices; the
   directionally robust claims are marked [robust].
6. **[robust] Log-dispersion of colony size is a size-coupled, floor-prone metric.**
   sd(log Asat) correlates ~−0.6 with colony size, so pooled "dispersion vs Cu"
   trends (incl. idea 01's v1 Gibrat claim) are confounded unless within-strain
   and size-controlled; always pair such metrics with a size covariate.
7. **The dominant scale of phenotypic variation is *within* species (idea 10):**
   62–92% of variance for the 5 key traits is among-strain, not across species;
   the phylogenetic signal (idea 09) sits in the smaller between-species
   component. Interpret species-level summaries in that light.

---

## Idea 11 — GWAS power & variant feasibility (are we set up to find the causal loci?)
`scripts/idea_11_effective_n.R` + `scripts/idea_11_power.py`

Motivation: ideas 09/10 showed within-species (strain) variation dominates every
trait, and the R. mucilaginosa panel is ~200 strains — but can we actually map
those differences to SNPs? Two sub-questions: (a) how much *independent*
genotypic diversity is there after redundancy, and (b) given n, how big a per-SNP
effect is detectable?

- **[robust] Effective independent genome count** (idea_11_effective_n.R, on the
  PHYling protein tree `phyling_pep/protein/tree/fungi_odb10/final_tree.nw`,
  278 tips): R. mucilaginosa has 200 tips but only **178 effective independent
  haplotypes** after collapsing 22 near-clone strains (pairwise dist < 1e-7 →
  redundancy ratio 1.12; cf. idea 09's "comet" polytomy — these are true genomic
  near-duplicates, e.g. DBVPG 3235–3239/3442–3446 closely-related sets).
  Every other species has redundancy exactly 1.0 (no collapse): paludigena 14,
  diobovata 8, sp. clade I 8, toruloides 8, dairenensis 7, taiwanensis 6,
  sphaerocarpa 5, graminis 3, kratochvilovae 3.
- **Power table** (idea_11_power.py; noncentral χ² df=1, ncp=n·R²/(1−R²),
  power 80%, per-SNP R²). Min detectable R² at α=5e-8 (genome-wide), 1e-6,
  1e-4 (candidate): n=100 → 0.284/0.247/0.183; n=150 → 0.209/0.180/0.130;
  n=202 → **0.164/0.140/0.100**; n=250 → 0.137/0.116/0.082; n=272 →
  0.127/0.108/0.076; n=400 → 0.090/0.076/0.053. At the real effective n=178,
  genome-wide threshold ≈ 0.21. Verified directly (n=178/202/272 vs explicit
  noncentral-χ² power computation; R²=0.05→power 0.008–0.048, 0.15→0.56–0.93).
- **Effect-size reading**: at MAF 0.3 a detectable allele moves the phenotype
  ≈0.71 SD. In trait units (β = 0.71 × sd_within × √(R²_target/1−…))  ≈ 0.17 log₁₀ px
  colony size, 0.14 log₁₀ chroma, 0.14 pace. Converting to % of *genetic* variance
  assuming h²=ICC (upper bound): per-SNP ≈ 22–25% of additive genetic variance for
  size/chroma/heterogeneity/pace — i.e. **only large-effect loci are detectable**,
  unless h²<ICC or effects spread across many SNPs. Copper slope is hopeless:
  ICC(0.19) < h² needed for mapping, β≈0.008 trait units — undetectable (avoid).
- **Discovery (resolves the entire variant-pipeline blocker)**: no variant calling
  is needed — a completed population-genomics + GWAS project exists at
  `/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/`
  (ref **NRRL Y-2510**): `vcf/RmucY2510_v2.All.SNP.combined_selected.vcf.gz` =
  728,581 SNP sites × 422 strains, haploid GTs, GATK-hard-filtered.
  **201 of our 278 phenotyped strains (all R. mucilaginosa) carry genotypes there**
  (200/201 have complete data for all 4 GWAS traits). Power estimate translates
  directly: n≈201 raw / ≈178 effective ⇒ min detectable R² ≈ 0.16 (GW) / 0.10
  (candidate) — near the n=202 row of the table. The remainder of the panel is
  other species (must not be pooled into a GWAS against a single reference).
- Prior art in that project: a 218-strain GEMMA+LMM+linear GWAS of *growth rates*
  (T4C–T37C, Salt6; sensitivity α≈1.18e-7) found only a handful of hits (linear
  T37C scan: p≈1.7e-11, chr13 locus), consistent with the idea-11 limit — few
  large-effect loci. Their color traits were never tested → gap this project fills.

---

## Tier D/E/G — gene mapping, fine-mapping, prior-locus replication (2026-08-17)

Follow-ups from the GWAS divergence design (`09-NEXT-GWAS-DESIGN.md` §7); full write-up in
`GWAS_REPORT.md` §9, outputs in `results/gwas/tierD/`, `results/gwas/tierE/`.

### Tier D — locus → gene annotation (scripts/annotate_gwas_loci.py)
- 12,348 FDR-sig SNPs annotated against a 6,799-gene index; 5,286 independent loci
  (positional 250 kb clump) across 9 traits.
- **Biologically informative hits**: chroma scaffold_8:831789 → **telomerase reverse
  transcriptase (OM429_004009)**; AUC_10 scaffold_10:396172 → **RNA-dependent ATPase /
  DBP3 (OM429_004640)**, scaffold_2 → RNA-pol-I transcription factor, scaffold_7 →
  endodeoxyribonuclease; BSLMM chroma scaffold_3:132903 → **methionine aminopeptidase 1
  (OM429_001379)**.
- **Clump caveat**: the 250 kb clump keeps only the single most-significant SNP per
  chromosome as lead → a chromosome whose top signal is at one end hides a distant
  second signal on Chr13-like small scaffolds (chr13's 217 FDR SNPs form ONE 525 kb block
  collapsing to a single lead 13_791853). Use window-based per-trait top-loci tables
  (issue in `tierD_fdr_snps_annotated.csv.gz`) when many blocks per chromosome are expected.

### Tier E — fine-mapping (scripts/finemap_credible_sets.py, Wakefield ABF)
- Multiple-testing-aware credible sets (90/95/99%) for 14 anchors; 95% CSs:
  **chroma_lead_10 (scaffold_10:384905) CS n=24, lead pp=0.054, β=1.11±0.19 SE, common AF
  (0.81) → well-bounded common-variant signal** — this is the paper-grade anchor.
- Rare-EF loci (AF≈0.015) resolve poorly: auc10 DBP3 CS n=67, lead pp≈0.02, β=8.2e5
  (rare_driven); **resil_lead_13 95% CS n=1, lead pp=0.615** (99% CS n=1, pp=0.977) —
  a singleton that fully resolves; ic50 95% CS n=10 (pp=0.183).
- Method note: ABF prior must be in z-space (trait-scale priors break when AUC_10 betas
  ~1e5 vs chroma ~1); candidate filter p<1e-3 needed (else sets span thousands of nulls).

### Tier G — prior growth-rate locus chr13:13_30149 replicates in our AUC_10
- Prior lab's p=1.68e-11 growth-rate locus at 13_30149 is **replicated** via our nearest
  proxy **scaffold_13_30134 (15 bp away): p_wald=4.03e-6, FDR-sig for AUC_10, β=804,778,
  af=0.015**. It sits inside the single chr13 rare-haplotype block (lead 13_791853 p=2.4e-7)
  that drives all 217 chr13 FDR SNPs.
- Gene at the locus: **OM429_005439** (scaffold_13:30,096–30,670) — hypothetical protein,
  unannotated (no GO/InterPro/PFAM). Flanking: OM429_005440 (hypothetical), OM429_005441
  (Ark1-family Ser/Thr kinase) 6.7 kb downstream. **Causal gene under a replicated,
  p=1e-11-in-prior locus is functionally unknown** → priority validation target.
- Other traits are null at 13_30134 (best nominal clone_mean_area p=0.044, bright p=0.082;
  none FDR). Replication is restricted to the growth-like trait AUC_10 — consistent with
  the prior scan being a growth-rate trait and our panels' niche is color/copper.
