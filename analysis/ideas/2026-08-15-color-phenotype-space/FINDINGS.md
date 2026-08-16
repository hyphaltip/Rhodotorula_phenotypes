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
