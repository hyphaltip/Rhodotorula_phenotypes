# Ideation Session — 2026-08-15 — Phenotype Space & Color Space Exploration

Question: *"What else should we explore for phenotypic data space or color space in
this dataset?"*

Session following the `idea-generator` convention (persona-based, data-grounded).
8 personas × 2 ideas = **16 ideas**. Shared context in `_context.md`.

---

## Idea Index

| # | Persona | Idea | Ambitious/Tractable | Effort | File |
|---|---------|------|---------------------|--------|------|
| 01a | Statistical Physicist | Universality of the arrest: scaling-collapse of colony growth curves + quenched-disorder phase diagram across Cu | Ambitious | High | `01-statistical-physicist.md` |
| 01b | Statistical Physicist | Roughness & ordering transition: intra-colony fluctuation exponents (CV_L, ΔC) + edge-vs-core asymmetry as early label-free Cu/strain discriminators | Tractable | Low–Med | `01-statistical-physicist.md` |
| 02a | Color/Imaging Scientist | Chroma-weighted hue space: separating pigment identity (hue angle) / amount (chroma) / density (L*) + origin of intra-colony heterogeneity (a*CoeffVar rises 0.31→0.54 with Cu) | Ambitious | High | `02-color-imaging-scientist.md` |
| 02b | Color/Imaging Scientist | Pigment appearance time: chroma/hue onset trajectories + within-run color calibration across the 5 imager runs | Tractable | Low–Med | `02-color-imaging-scientist.md` |
| 03a | Information Theorist | Minimal sufficient statistic of the phenotype screen: entropy/redundancy decomposition of the ~200-col feature space | Ambitious | High | `03-information-theorist.md` |
| 03b | Information Theorist | Does Hue/Texture carry species information that L* does not? (conditional MI, PID) | Tractable | Low–Med | `03-information-theorist.md` |
| 04a | Quantitative Geneticist | Heritability & repeatability audit of the ~200-ccol feature space + strain × Cu GxE for color traits | Ambitious | Med–High | `04-quantitative-geneticist.md` |
| 04b | Quantitative Geneticist | Strain-level repeatability of intra-colony pigment heterogeneity (StdDev/CoeffVar) | Tractable | Low–Med | `04-quantitative-geneticist.md` |
| 05a | Representation Learning | Intrinsic phenotype atlas: bootstrap-stable latent map (rank-PCA/varimax, UMAP+HDBSCAN) of the full feature space | Ambitious | High | `05-representation-learning.md` |
| 05b | Representation Learning | Parts-based "pigment parts" NMF on the chromatic block (Saturation/Chroma/a*/b*/xy + heterogeneity) | Tractable | Low–Med | `05-representation-learning.md` |
| 06a | Causal Inference | Decomposing the copper phenotype: growth-density mediation, latent-tolerance structure, well/plate interference | Ambitious | High | `06-causal-inference.md` |
| 06b | Causal Inference | Endpoint mediation: is Cu's chromatic effect direct, or routed through growth density? (NDE vs NIE at 90–110 h) | Tractable | Low–Med | `06-causal-inference.md` |
| 07a | Trait Ecology | Colony trait spectrum: pigment-investment vs growth-defense tradeoffs, where Cu tolerance lives in trait space | Ambitious | High | `07-trait-ecology.md` |
| 07b | Trait Ecology | Do marsh_tidalflat isolates carry a shared pigment + Cu-tolerance syndrome? (two-group trait filter, species-conditional permutation null) | Tractable | Low–Med | `07-trait-ecology.md` |
| 08a | Temporal Dynamics | Pigment ontogeny atlas: full trajectory reconstruction, functional clustering, area–color coupling over time | Ambitious | High | `08-temporal-dynamics.md` |
| 08b | Temporal Dynamics | Onset, not endpoint: when does each strain pigment, and how does Cu move that moment? | Tractable | Low–Med | `08-temporal-dynamics.md` |

---

## Cross-cutting themes (synthesis)

1. **The ~170 unused feature columns are the frontier.** Every prior analysis used only
   Shape_Area, per-colony L*/a*/b* Medians, and Intensity Mean. The unused space —
   intra-colony heterogeneity (StdDev/CoeffVar, Q1/Q3), HSV, xy chromaticity, chroma,
   ~52 TextureGray Haralick features, morphology shape descriptors — recurs as the raw
   material in 14 of 16 ideas.

2. **Hue, chroma, and lightness must be treated as three distinct axes** (color-science
   insight): L* ≈ Intensity (density, collinear ≥0.99); chroma = pigment amount;
   CIELAB hue angle = pigment identity. HSV Hue is likely **degenerate** in this pipeline
   (verified live: ColorHSV_HueMedian ≈ 0.1 everywhere) — use CIELAB-derived hue with
   chroma-weighted circular statistics, never arithmetic hue means.

3. **Timing is under-exploited** — pigment onset time (t_chroma_onset), darkening onset,
   hue-stability time all build directly on the existing 0–117 h time-course and the
   control-table pass-rounding trick; several ideas agree Cu likely acts on *amount/yield*
   and *timing* more than on *identity*.

4. **The replicate design is a strength in reserve** — 16 plates/dose, replicate wells
   per strain: unused muscle for ICC/heritability (04a), bootstrap collapse residual (01a),
   and interference tests (06a).

5. **Risks to respect** (from existing findings): late-time colony confluence inflates
   roughness/Q3−Q1 (filter on Solidity/Extent); plate identity nested in Cu (need run/plate
   reference centering); R. mucilaginosa dominates the panel (216 strains — but that is a
   strength for within-species tests); 15–18% high-Cu dropout censors late curves
   (survival-style treatment for onset-∞).

## Suggested sequencing (tractable-first)

- **T0 (now):** 02b pigment onset + run calibration → directly reusable normalization.
- **T0:** 01b early roughness discriminator (Low–Med, pure v_phenotype).
- **T1:** 02a chroma-weighted hue space (needs 02b's calibration constants).
- **T1:** 03b conditional MI on unused features (cheap, sets which columns matter).
- **T2:** 04b heterogeneity repeatability, 05b pigment-parts NMF, 08b onset models.
- **T3 (thesis-scale):** 01a scaling collapse, 04a heritability+GxE, 06a mediation.

Everything is implementable on existing `v_phenotype` data via the pixi env (DuckDB +
R lme4/circular/changepoint); no new imaging required.
