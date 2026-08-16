# Trait-Ecology Ideas — Two Proposals

Two research ideas from a Trait-Ecology / Macro-Trait-Phylogeneticist lens on the
Rhodotorula copper screen. Both use the same mental model: a colony is a
functionally expressed phenotype drawn from a trait distribution whose *axes* we
can estimate, and whose *structure* (who clusters, who filters by environment,
who trades off growth for defense) is what trait ecology is about.

---

# Idea 1 (AMBITIOUS): The Colony Trait Spectrum — pigment-investment, growth-defense tradeoffs, and where Cu tolerance lives in trait space

## Persona
**Trait-Ecology / Macro-Trait-Phylogeneticist** — the lens is functional-trait
spectra and null models for trait structure: do species occupy distinct regions
of phenotype space (niche separation) or converge on shared optima (trait
convergence)? Is tolerance a trait that rides a trait axis, i.e. correlated with
carotenoid investment and with growth, along a fast–slow / grow-fast-vs-invest-in-defense
continuum?

## Motivation
In macro-ecology we rarely get to measure ~300 individuals of a clade for a
hundred traits under a controlled stress gradient. Here we do: ~320 yeast
strains (Rhodotorula, Cystobasidium, Pseudomicrostroma), ~200 colony features,
and a designed 0–30 mM Cu gradient, all in one image system. The obvious
ecological question is then the *trait-spectrum* one: is there a dominant
low-dimensional axis (a "colony fast–slow spectrum") that orders strains from
"big, pale, fast expanding, Cu-sensitive" to "small, deeply pigmented, slow,
Cu-tolerant"? Carotenoid pigmentation is costly and protective (membrane
antioxidant); a negative growth–chroma correlation with a positive
chroma–tolerance correlation is exactly the growth-defense tradeoff morphology
would predict for the genus (Rhodotorula as industrial carotenoid producers).
The existing single-trait analyses (endpoint L*a*b* tables, extent growth-rate,
species log2 sensitivity) each see one column of the matrix; the un-tested claim
is that these are **columns of one or two latent trait axes**.

## Connection to Existing Data
- The ~170 unused `v_phenotype` columns make up nearly the whole trait matrix —
  the existing analyses use only `Shape_Area`, `ColorLab_{L,a,b} Median`, and
  `Intensity_MeanIntensity`. The morphology block (`Shape_Solidity, Eccentricity,
  Extent, BboxArea`), the intra-colony heterogeneity statistics (`StdDev`, `CoeffVar`
  of every color space), the `ColorHSV` Hue/Saturation and `ColorLab chroma`
  (`ChromaEstimatedMedian/Mean`) descriptors, and all 52 `TextureGray` Haralick
  features are available in the same DuckDB view with zero new acquisition.
- Existing outputs to pull (not recompute): per-strain Cu tolerance from
  `analysis/growth_rates/results/tables/strain_sensitivity.csv` (log2 extent
  ratio, 25–30 vs 0–5 mM), the flat ~37 h doubling time (knowing `dbl` is
  Cu-neutral tells us growth-defense axes must be *extent*- and
  *pigment*-anchored, not rate-anchored), and `rate_area` vs `rate_int` being
  only r = 0.21 (two genuinely distinct growth axes to feed the spectrum).
- Metadata `strain.origin` (China 180 / Italy 119 / America 9) and
  `strain.environment` (marsh_tidalflat 162 / soil 45 / plant 24 / food 24 /
  unknown 12 / marine 4 / cave 3 / air_cloud 2 / insect 2 / sand 2 / snow_ice 2 /
  rock 2 / built env 1 / water 1) give the environment and geography labels for
  filtering tests.

## Approach
1. **Define the trait set (strain-level, control Cu = 0 only first).** From the
   ~170 unused columns + 2 existing indices, curate ~12 non-redundant traits:
   (i) expansion: `Shape_Area` at plateau, `Shape_Solidity`/`Extent` (packing);
   (ii) pigment amount: control `ColorLab_a* Median`, `ChromaEstimatedMedian`,
   `Saturation (HSV) Median`; (iii) hue: circular mean of `ColorLab_b* Median`
   (keep b* separate from hue to avoid circular-mean artifacts; context warns
   against arithmetic means of hue); (iv) intra-colony heterogeneity:
   mean `CoeffVar` across `ColorLab {a*,b*}` (edge-core pigment structure — the
   "is the colony uniform or patterned" axis); (v) texture: `TextureGray`
   `InverseDifferenceMoment` + `Entropy` angle-averaged (granularity);
   (vi) existing rate/tolerance indices: `rate_area`, `rate_int`, log2-extent Cu
   sensitivity, saturation-drop. Centroids across replicate colonies (≤4
   run-plates) per strain; flag n-rep < 2 rows as in the growth-rates analysis.
2. **Rank-reduce and name the spectrum.** PCA (and, for robustness, a
   correlations/IMFS redundancy audit following the context's "how much of the
   200-col space is redundant" direction) on centered-scaled traits; retain PCs
   with broken-stick-eigenvalue criterion. Test the prediction that the dominant
   PC is interpretable as *pigment-investment* (PC loads chroma/a*/CoeffVar
   positively, Area or rate negatively) and that **PC2** separates `rate_area`
   vs `rate_int` (two growth axes). Report whether Cu-tolerance vectors (log2
   extent, saturation drop) align with PC1 — the claim "tolerance is a trait
   living on the pigment-investment axis".
3. **Species niche separation vs convergence (null model 1, permutation of
   species labels).** Compute species trait-space centroids (20 named species +
   `unknown` handled per step 5). Measure pairwise Mahalanobis D² between
   species centroids and the fraction of species whose 95% trait-space ellipses
   overlap the global centroid — under a null where the species label is
   permuted across strains *within* genus. If realized species separation
   exceeds the permutation null, species occupy distinct trait niches
   (separation); if not, the colony range is a single convergent spectrum
   (continuous trait variation, not discrete morphs). Also report whether the
   "pigment morph" clusters the color-space directions describe are species
   collapses or cuts across species.
4. **Environment → trait filtering (null model 2, permutation of environment
   labels).** Collapse environment to a binary marsh_tidalflat (162) vs
   non-marsh (~124 named), given every other environment holds 1–45 strains
   (small effective sample; keep them as verbose posterior points, not testable
   bins). Test PC1 (pigment-investment) and log2-extent tolerance by environment
   with a **species-conditional** permutation: shuffle environment labels only
   across strains of the same species (so species identity is held fixed and the
   test asks "net of who-you-are, does the isolation *habitat* predict the trait
   axis?"), 10,000 perms, report observed minus null-median effect size. Because
   the marsh bin dominates, this is really "are marsh isolates
   off-spectrum relative to everything else". Origin (China/Italy/America) is the
   same fixed-factor logic at coarser resolution — phylogenetic confound
   acknowledged (origin ≈ lineage in a sample of ~320).
5. **Handle non-species-typed strains honestly.** Strains with `species =
   'unknown'` (12) and `NULL` species (12) are excluded from species-niche tests
   (step 3) but kept as passive points in the PCA once loading directions are
   fixed (project into spectrum, never rotate it with them). In the
   environment-filter test all ~24 are kept under the marsh/non-marsh collapse
   except the 12 `environment = unknown`, which are handled separately as a
   "where would trait-space nearest-neighbors place these?" sanity check (see
   Idea 2 step 4).
6. **Validate.** (a) Repeat the PCA and separation tests with the log2-extent
   index recomputed strain-level on control vs 30 mM only (robustness variant of
   the existing table); (b) jackknife by run_number (drop one imager run) to
   confirm PC loadings are not run artifacts, given the plate-position analysis
   showed plate variance is small (L* 0.1–4.7%) but nested in Cu; (c) report all
   effect sizes (Cohen's d for environment contrasts, D² for niches) not just p.

## Expected Insights
We'd learn whether the ~200-column colony feature space collapses to a small
named trait spectrum (pigment-investment × growth-mode), whether that spectrum
has structure (species niches) or is a single convergent packing axis, whether
Cu tolerance *lives on* the pigment axis (prediction: chroma-positive PC1
correlates with higher log2 extent, i.e. more tolerant) — evidence for or
against "carotenoid pigment is a Cu-protective trait", — and whether the
marsh_tidalflat habitat leaves a detectible trait signature beyond species
identity. If the spectrum claim fails (e.g., pigment and growth are orthogonal),
that itself is the finding: tolerance would then be an idiosyncratic species
trait, not an axis coordinate.

## Feasibility
- **Effort**: High — multivariate trait pipeline (PCA, permutation engines,
  species-ellipse D²), multiple null models, run-jackknife validation.
- **Data ready**: Yes — all feature columns already sit in `v_phenotype`;
  tolerance indices exist in `growth_rates/results/tables`.
- **Methods available**: Standard (R: prcomp/FactoMineR, ade4-style permutation;
  all base stats); needs a modest custom permutation loop, no new packages.
- **Key risk**: The huge strain-count skew toward one species (R. mucilaginosa,
  216) means species-niche tests are dominated by one ellipse; trait-spectrum
  axes could be driven by the mucilaginosa-vs-everything split, and small
  per-environment effective samples (1–45) restrict the filter test to a single
  contrast (marsh vs rest). Mitigate with within-species conditional
  permutations and by separating "between-species" and "within-mucilaginosa"
  spectra.

---

# Idea 2 (TRACTABLE): Do marsh_tidalflat isolates carry a shared pigment + Cu-tolerance syndrome? A two-group trait-filter test with a species-conditional permutation null

## Persona
**Trait-Ecology / Macro-Trait-Phylogeneticist** — the lens is environmental
filtering: when one habitat dominates a collection (marsh_tidalflat 162/286
named isolates), the cheapest ecologically meaningful question is whether the
dominant habitat's isolates express a shared trait syndrome (pigmented →
defensive → tolerant) that distinguishes them from all other habitats pooled.
This is the "do strains from the pigmented-origin habitat share traits" claim,
made falsifiable.

## Motivation
The collection is overwhelmingly coastal: marsh_tidalflat holds 57% of named
isolates, far ahead of soil (45), plant (24), food (24), and all other
micro-environments (≤4 each). Coastal tidal flats are carotenoid-culture
grounds (UV + osmotic + metal stress selecting for pigmented,
antioxidant-equipped microbes), and the same pigment that protects them there is
the industrial carotenoid product (torulene/β-carotene). A trait-ecology prior:
marsh isolates should be *more chromatically invested* and *more Cu-tolerant* at
once, because both are on the defensive end of a pigment-investment axis. The
12 `environment = unknown` isolates are a ready-made out-of-sample: if the marsh
syndrome is real, trait-space nearest neighbors should assign them to marsh more
than to non-marsh by trait chance alone — a light-touch imputation sanity check.

## Connection to Existing Data
- Metadata gives `environment` counts directly (marsh 162 vs soil 45 vs plant 24
  vs food 24 vs unknown 12 vs marine 4 vs cave 3 vs air_cloud 2 vs insect 2 vs
  sand 2 vs snow_ice 2 vs rock 2 vs built env 1 vs water 1) and `origin`
  (China 180 / Italy 119 / America 9).
- Existing per-strain Cu tolerance: `analysis/growth_rates/results/tables/
  strain_sensitivity.csv` (log2 extent, saturation drop, doubling fold) and the
  species table showing R. mucilaginosa the tolerant species (−0.72) vs
  R. kratochvilovae (−4.13) — we reuse, not recompute.
- Control-endpoint color from `control_late_timepoint_phenotype` per-strain
  tables: `ColorLab_a* Median`, `ColorLab_b* Median`, plus `ChromaEstimatedMedian`
  and `CoeffVar` columns from `v_phenotype` (Cu = 0 only) — the pigment trait set.
- No new imaging or extraction: everything is in the DuckDB `v_phenotype` view +
  2 existing result tables.

## Approach
1. **Build a strain-level trait matrix at Cu = 0.** For the marsh-vs-non-marsh
   contrast only, use a small trait set (≤6): pigment amount (control
   `ColorLab_a* Median`, `ChromaEstimatedMedian`), intra-colony pigment
   heterogeneity (mean of `ColorLab_b* CoeffVar` across replicate colonies),
   colony size (`Shape_Area` plateau median), and the two existing Cu-response
   reads from `strain_sensitivity.csv` (log2 extent, saturation drop). Aggregate
   replicate colonies to strain medians; drop `n < 2` replicate strains from the
   contrast and report how many that removes.
2. **Two-group test with a crisp species-conditional null.** Compare marsh (162)
   vs named non-marsh (~124) on each trait with a Wilcoxon + scaled mean
   difference (Cohen's d). Null: permute the marsh/non-marsh label 10,000 times
   **within species only** (strain labels shuffled within each species), so the
   null preserves both the species structure and the marsh skew; the observed d
   is significant only if it exceeds the within-species permutation distribution.
   This directly answers "marsh predicts trait *net of species*, or is this just
   R. mucilaginosa being over-sampled in the marsh?"
3. **Confound control, the decisive sub-test.** Because R. mucilaginosa (216
   strains) is both the collection's majority and its most tolerant species, run
   the entire step-2 contrast **inside R. mucilaginosa only** (the one species
   with enough marsh/non-marsh spread). If marsh vs non-marsh still separates on
   chroma and/or log2-extent within this single species, the environment-filter
   signal is trait-level, not taxonomy-level. Any other species with ≥4 strains
   on both habitat sides joins a small per-species replication table.
4. **Use the 12 `environment = unknown` strains as out-of-sample.** Fix the
   discrimination (e.g., logistic on the ≤6-trait set fit on marsh vs named
   non-marsh), then predict the 12 unknowns and compare their predicted habitat
   posterior to the 50% chance expectation — under the null that the unknown
   isolates are drawn from the non-marsh pool (permutation of which 12 strains
   are "unknown"). Species-NULL (12) and `species = unknown` (12) strains are
   excluded from step 3's within-species test but included in step 2's
   permutation (species identity then just constrains shuffle structure).
5. **Validate and report.** (a) jackknife strains (drop 10% at a time) for
   stable d; (b) repeat with a*b*-Median endpoint tables from
   `control_late_timepoint_phenotype` instead of `v_phenotype` averages to check
   the aggregation window doesn't flip the sign; (c) state effect sizes (d, and
   the within-mucilaginosa d) with CIs, not p-only.

## Expected Insights
We'd learn whether the dominant coastal habitat prints a measurable trait
signature — both the headline (does marsh → more chroma/higher tolerance hold
net of species?) and the mechanism slice (is the marsh association actually
carried entirely by R. mucilaginosa's ubiquity, i.e., a sampling artifact rather
than filtering?). A null result at step 2–3 is still a result: it would say the
Cu-tolerance and pigment traits assort by *lineage*, not by *habitat*, in this
collection — the cleanest trait-vs-phylogeny statement this dataset can make.
The unknown-environment imputation gives a free, honest preview of whether
habitat is even predictable from the colony phenotype.

## Feasibility
- **Effort**: Low–Medium — a small R script (Wilkoxon + within-species
  permutation + within-mucilaginosa subset + 12-unknown classification); hours,
  not days.
- **Data ready**: Yes — all columns in `v_phenotype`; tolerance indices already
  in `growth_rates/results/tables/`.
- **Methods available**: Standard (base R + any permutation helper; no new
  tooling).
- **Key risk**: The contrast is marsh vs a composite "everything else" whose
  internal variance (soil vs plant vs food) is baked into non-marsh; and within
  R. mucilaginosa the marsh/non-marsh split may be so skewed that the
  within-species comparison has few effectively independent isolates. Small
  per-habitat samples (≤4 for 8 of 15 habitats) cap us at a single binary
  contrast — high-resolution per-habitat trait filtering is not feasible.
