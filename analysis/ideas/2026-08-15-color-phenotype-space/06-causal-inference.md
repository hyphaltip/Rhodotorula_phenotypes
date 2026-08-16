# Causal-Inference Methodologist — Ideas

Two ideas from a causal-inference lens: one ambitious (full structural model with
mediation + interference + latent tolerance), one immediately tractable (endpoint
mediation decomposition). Both are grounded in the `v_phenotype` abundant layer,
the growth-rates findings (extent-limited Cu effect; doubling time ~37 h flat
across 0–30 mM; L* collinear with intensity >0.99), and the plate-position
finding (plate identity nested in Cu; neighbor adjacency signal only for L*).

---

# Ambitious Idea 1 — Decomposing the copper phenotype: growth-density mediation, latent-tolerance structure, and well/plate interference

## Persona
**Causal-Inference Methodologist** — identifies effects of a randomized dose in a
nested observational-within-experiment design, worries about mediation through
post-treatment variables, interference between units, collider selection from
Cu-dependent dropout, and whether the "tolerance" construct is a single latent
trait or a dose threshold.

## Motivation
The dose is randomized at the **plate** level (each plate carries exactly one
`copper_mm`, 16 plates/dose, interleaved across 5 imager runs), which is the one
clean causal lever in this dataset — so the interesting identification problems
are not confounding, they are (a) **composition**: Cu acts partly *through*
growth extent and partly directly on per-cell pigment chemistry, and the two
readouts share a measurement axis (L* is 0.99-collinear with intensity, so "dark"
means "denser" at the per-pixel optics level); (b) **selection**: the 15–18% of
(strain, plate) cultures that drop out before the late timepoints at 25–30 mM are
exactly the sensitive ones — restricting the endpoint sample to survivors closes
a collider that sits on the mediation path; and (c) **interference**: the
documented plate-position neighbor-adjacency signal for L* (and its built-in
"similar strains were placed near each other by design" confound) means colony
outcomes may depend on their well-neighbors, which changes what the 8,446 colonies
actually replicate. A causal diagram of these layers changes *which* estimator is
defensible at all.

## Connection to Existing Data
- Dose gradient `copper_mm` 0–30 (7 levels), 16 replicate plates/dose, 112
  experimental plates, 5 runs, strain reassigned to a *fresh plate per (run,dose)*
  → up to 4 physically independent (plate, well) colonies per strain/dose (the
  plate-position analysis's own design note: `plate_id` nested in `copper_mm`).
- growth-rates already established the **extent-limited** picture: exponential
  doubling time flat ~37 h, but saturation fraction 95.4%→73.9%, t50 72.6→64.1 h —
  the clean *growth-density mediator* M. `rate_area` vs `rate_int` r = 0.21 tells
  us extent and intensity/lightness are quasi-independent axes, so a two-mediator
  structure is warranted.
- explore-plate-position gives us the nuisance and interference geometry:
  within-Cu plate variance small (L* 0.1–4.7%), naive pooled model overestimates
  plate variance (plate ≈ dose → not estimable as a confounder), and an adjacency
  effect on **L* only** (4-connected grid neighbors), which we treat as a
  hypothesis, not a fact, because the design clusters similar strains spatially.
- The unused ~170 columns are exactly what a per-cell pigment decomposition needs:
  intra-colony pixel StdDev/CoeffVar and Min/Q1 vs Q3/Max (surface vs core),
  HSV/xy/chroma, and the 52 Haralick texture columns — i.e., descriptors of *how*
  pigment is distributed inside a colony, separable from *how much* colony there is.
- Strain metadata (`origin`: China 180/Italy 119/America 9; `environment`:
  marsh_tidalflat 162, soil 45, plant 24, food 24, …) give effect-modifier contrasts
  (e.g., marsh_tidalflat strains as an a priori tolerant group to check whether the
  tolerance latent trait tracks ecology).

## Approach
Causal diagrams. Layer 1 (per-colony, mediation + selection):

```
  Cu dose (plate-rand., interleaved) ──────────────────────────────► chromatic Y(t)
       │                                                                 (a*,b*,chroma,hue,
       │  ┌── indirect via M ──────────────┐                             not L*: L*=intensity)
       │  v                                │
       ├─► extent/density M(t) ──► survival S(t) ─► measured Y(t)│S=1
       │     (Shape_Area, saturation,        (S: dropdown 15-18%    ^
       │      t50; extent-limited effect)     at 25-30 mM =            collider:
       │                                      collider, DAG edge        do NOT condition
       │                                      from latent θ too)        on S alone
       └─► direct per-cell chemistry (oxidative stress → carotenogenesis)
       latent tolerance θ (strain) ──────────► M(t), Y(t), S(t)   [unmeasured
       species/origin/environment ──────────► θ, Y                  Messier "odd" link]
       neighbor well j on same plate ────────► Y_i  (interference; documented L*)
       plate row×col segment / edge ────────► Y    (nuisance; motivates stratification)
       run (batch) ── stratified nuisance ──► Y    (instrument NOT needed: dose is random)
```

Concrete steps:

1. **Build the analysis-ready longitudinal + graph layer.** Reuse the growth-rates
   series construction (`hours_since_plate_start` 0–117 h, 6-hourly, runs 353–356)
   but retain, per colony-timepoint: Shape_Area, saturation/t50, per-colony
   `a*/b*` Median + StdDev + CoeffVar, `ChromaEstimatedMedian`, hue (circular — use
   directional stats, never arithmetic mean), and 2–3 Haralick entropy/difference
   columns. Construct the **96-well 4-connected adjacency graph** per plate (the
   plate-position analysis already keyed adjacency to the grid, not to
   `well_position` column-major order) and a **plate-segment** factor (row×col
   block + `is_edge`) as the design-stratification nuisance.

2. **Plate-as-unit design-based estimator (population-level, blend-invariant).**
   Estimate the plate-level ITT for dose → (plate-mean outcome) over 112 plates
   with cluster-robust SEs. This estimate is *immune to within-plate interference*
   and to within-well mis-replication (it never uses well-within-plate as
   replication), and it is unbiased under the interleaved randomization with no
   instruments needed. This is the anchor number every model below must reproduce.

3. **Longitudinal mediation with two mediators.** Fit a g-computation/sequential
  -regression decomposition where Cu → final chromatic endpoint Y (at fixed window,
   e.g., 90–110 h matching the control-late-window cadence) is split into:
   - NIE via M1 (extent/density trajectory: Shape_Area@endpoint, saturation),
   - NIE via M2 (intra-colony heterogeneity: a*/b* StdDev/CoeffVar @endpoint —
        the "pixel optics / surface-vs-core" pathway),
   - NDE (direct per-cell chromatic response at fixed density).
   Identify M-path effects *within strain* by differencing across doses for a fixed
   strain (strain fixed effects kill between-strain latent-θ confounding of M→Y,
   the main unmeasured-mediator-outcome threat), then average over strains. Report
   both the weighted and cluster-bootstrapped CIs (cluster = plate for dose terms,
   well for within-plate terms).

4. **Correct the collider: Cu-dependent dropout.** Fit a survival/reach-endpoint
   model `P(reach 100 h | Cu, species, strain, early M(0-24h))` and IPW the
   mediation/decomposition, and re-run the plain decomposition on the complete-case;
   the gap between the two is the *selection bias from conditioning on survivors*.
   Also refit with outcome = "last-observed pigment carried to plateau + dropout
   indicator" so missingness is part of the outcome vector, not deleted.

5. **Interference, formally.** (a) *Within-well*: code wells with >1 segmented
   object (`object_label` collisions; the plate-position analysis logged 3 such
   wells at endpoint but the time-course has more). Compare within-well
   across-object variance of the endpoints vs across-plate same-strain same-dose
   variance; run a nested mixed model `Y ~ (1|plate/well/object_label)`. If
   within-well spread ≈ across-plate spread, the extra objects are near-pure
   interference/mis-replication — drop them and show effective-n inflation in a
   sensitivity column. (b) *Between-well (plate) neighborhoods*: after strain-fixed
   + plate-segment residualization, fit a conditional autoregressive (CAR)/spillover
   model of `residual Y_i` on mean residual of 4-connected neighbors, *separately
   per trait* (test whether adjacency extends beyond the documented L*, e.g., to
   a*/b*/Haralick). Run a **placebo** that re-permutes well identities within each
   plate (1000×) to break the "similar strains placed nearby by design" confound —
   only keep the adjacency claim if the observed neighbor coefficient exceeds the
   placebo distribution. (c) Because dose is plate-constant there is no
   *cross-dose* spillover; the interference channel to model is local
   microenvironment (Co-deposition shadows, drying) — proxy with plate-segment +
   neighbor residual traits, and verify with the population ITT from step 2 as the
   hard reference.

6. **Latent-tolerance structure: one trait or a threshold?** Take the strain-level
   dose→extent/pigment response surfaces and fit three nested models: (i) a single
   latent tolerance dimension with dose as "difficulty" (Graded-response/IRT or a
   1-factor mixed model on strain×dose responses); (ii) a threshold model where
   extent collapses past a strain-specific Cu breakpoint (change-point per
   strain); (iii) dose-specific traits (unrestricted). Compare via LOO-CV and
   report, per strain, whether the threshold model is selected and where the
   breakpoint sits (expect it low for R. kratochvilovae/araucariae/taiwanensis and
   absent for R. mucilaginosa given their −4.13 vs −0.72 log2-extent indices).
   Then test whether the *chromatic* direct effect and the *extent* effect are the
   same latent variable (confirmatory factor/multitrait model) — i.e., is
   "tolerance" one construct or two?

7. **Heterogeneity and batch.** Strain-level NDE/NIE estimates → regress on
   species, `origin` (China/Italy/America), `environment` (marsh_tidalflat vs
   others), checking whether marsh/China strains show both smaller extent effects
   and smaller direct pigment effects (a shared driver) or only one. Keep `run` as
   a 5-level fixed/strata nuisance throughout (it is exchangeable by design, so it
   is precision, not a confounder); as a robustness framing, re-estimate the dose
   contrast within-run only (removes any residual inter-run drift bias).

## Expected Insights
- An honest split of the Cu→final-color effect into **direct per-cell chromatic**
  (chemically regulated carotenogenesis) vs **through growth-density** (extent +
  intra-colony pixel structure/optics). If NDE ≫ NIE, Cu rewrites pigment
  composition per cell (a real "pigment phenotype"); if NIE dominates, most of the
  color shift is a screened re-scaling of colony density — a claim none of the
  existing endpoint/L*a*b* tables can make because they regress color on area as
  a correlational covariate, not through a mediation decomposition.
- A verdict on whether "copper tolerance" is a **single latent trait** (strain
  ordering consistent across 0–30 mM, graded) versus a **threshold phenotype**
  (per-strain collapse dose). This changes how the −4.13/−0.72 species sensitivity
  ranking should be interpreted mechanically, and whether a single IC-style dose
  can ever summarize a strain.
- Whether the documented **L*-adjacency** is biological interference, a
  segmentation/lighting artifact, or a design artifact (placebo permutation
  adjudicates); and whether interference is truly trait-specific (L* but not
  a*/b*).
- A quantitative **replication-vs-interference** ruling on multi-object wells —
  correcting effective sample sizes for the whole project, not just this analysis.
- Whether within-colony pixel heterogeneity (StdDev/CoeffVar, core-vs-edge) is a
  *mediator* that carries part of the direct effect (pigment packing under stress)
  — connecting the "phenotypic data space" goal to a mechanism.

## Feasibility
- **Effort**: High (longitudinal g-computation, survival IPW, CAR + placebo
  permutation, IRT/threshold model comparison — a small research pipeline over the
  time-course, not just the endpoint).
- **Data ready**: Yes — all inputs live in `v_phenotype` (runs 353–356 have full
  time courses; strain/origin/environment metadata available). Needs reprocessing
  to keep per-colony pixel-stat columns and to build the well-adjacency graph
  (not a new experiment).
- **Methods available**: Research-grade — mediation/censoring/g-computation are
  standard but the CAR-with-placebo and the threshold-vs-graded model comparison
  need custom implementation (R `lavaan`/`mediation`, survival/IPW, custom
  permutation; DuckDB for the wide table).
- **Key risk**: no sensitivity analysis can rescue the within-strain mediation from
  residual M→Y confounding if latent tolerance also drives *intra-colony* pigment
  packing; and plate-level power is limited (16 plates/dose) for the interference
  and ITT components — so the mediation claims rest on assumption checks that the
  design cannot fully verify.

---

# Tractable Idea 2 — Endpoint mediation decomposition: is Cu's chromatic effect direct, or routed through growth density?

## Persona
**Causal-Inference Methodologist** — applies the same causal mental model but to
the *endpoint* layer only: a small, immediately runnable decomposition with honest
clustering and an explicit censoring correction, delivering a decision-relevant
number rather than a full structural model.

## Motivation
The growth-rates analysis proved Cu is extent-limited on the *biomass* axis
(flat ~37 h doubling; saturation 95.4%→73.9%; t50 72.6→64.1 h). A disciplined
causal reading of the *color* axis then asks the same question: whenever Cu
shifts endpoint a*/b*/chroma, is that because (i) the colony is smaller/less dense
so the same per-cell pigment reads differently at the pixel level (and total
pigment amount trivially scales with biomass), or (ii) Cu actually changes
per-cell pigment composition/loading directly? The existing color v Cu tables
(Cu F ~ 900 on L*, ~1100 on a*) are *total* effects; they cannot answer this, and
L* itself cannot be a chromatic outcome at all (0.99-collinear with intensity).
At the endpoint this is a two-mediator decomposition that can be estimated by
sequential regression on an already-buildable table — tractable now.

## Connection to Existing Data
- Endpoint datasets already exist from two analyses: the plate-position endpoint
  table (last timepoint ~117 h, 8,446 colonies, runs 353–356, `strain_id`,
  `copper_mm`, `plate_id`, `well_position`, grid position, `is_edge`) and the
  growth-rates per-colony summary (max Shape_Area, saturation fraction, t50,
  `rate_area`/`rate_int`) — plus the control-late-window cadence (70–80/80–90/
  90–110 h) showing the project's notion of "endpoint."
- The per-colony pixel statistics (a*/b* StdDev, CoeffVar; `ChromaEstimatedMean/
  Median`; hue Quartiles) and Haralick columns are **already in `v_phenotype`** —
  the "unused ~170 columns" — so the mediator and outcome variables need no new
  imaging.
- The documented dropout (15–18% of strain-plates never reach late timepoints at
  25–30 mM) is the selection channel this idea corrects with IPW; the
  plate-position "naive pooled model overestimates plate variance" result tells me
  to cluster at the right level (`plate_id` for dose terms, `well` within plate for
  M terms).

## Approach
DAG (endpoint, one time window; L* excluded as outcome):

```
  Cu dose (plate-rand.; one dose/plate) ──► M = extent/density       ─┐
       │                                    (max Shape_Area, sat.,     │ indirect
       │ direct per-cell chromatic          intra-colony a*/b* Stud.) ─┤→  Y = endpoint
       │ (oxidative stress → carotenogenesis)                          │    a*/b*/chroma Median
       └───────────────────────────────────────────────────────────────┘    + a*/b* CoeffVar
   strain (random; species/origin/env fixed) ──► M, Y      run ── stratified nuisance
   S = reached endpoint;  Cu ─► S ─► observed Y|S=1  (collider; handle by IPW, not deletion)
```

1. **Assemble the endpoint decomposition table.** For runs 353–356, one row per
   (plate, well) at the ~90–110 h window (matches the control cadence so pigment
   is mature): exposure `copper_mm`; mediators M = {max Shape_Area (extent),
   saturation fraction, a*/b* within-colony CoeffVar (pixel heterogeneity)};
   outcomes Y = {a* Median, b* Median, `ChromaEstimatedMedian`, hue Median
   (directional), a*/b* CoeffVar}. Drop the L*/intensity axis explicitly (it is
   not chromatic). Record S = reached window.

2. **Correct survival-selection first.** Logistic `P(S|Cu, species, strain)`; build
   IPW weights; verify balance. If weights are extreme at 25–30 mM (heavy loss),
   report the complete-case and weighted decompositions side by side — the
   difference is the collider bias magnitude.

3. **Mediation by sequential regression / g-computation** (a two-mediator,
   unexposed-value-parameterized version, estimable in base R):
   - M1/M2 models: `M ~ Cu(factor) + (1|strain) + run` (mediators as functions of
     dose), plus `Y ~ Cu + M1 + M2 + (1|strain) + run`.
   - Compute NDE = contrast in Y when Cu changes but M is fixed at its 0 mM
     potential (via coefficients/hold-out prediction); NIE = contrast in Y when Cu
     is fixed and M shifts by Cu's effect on M. Bootstrap over **plates** (not
     wells) for CIs; report NDE, NIE_via_extent, NIE_via_pixelStruct, and their
     fractions of the total effect, per species and pooled.
   - Run the same on 30-mM-vs-0-mM and on 5-mM-vs-0-mM (low-dose) separately — the
     decomposition at a mild dose is where a "direct pigment regulation" signal
     would be cleanest, since extent is barely limited at 5 mM.

4. **Honesty checks.**
   - *Within-strain version*: redo M2 and Y models with strain fixed effects
     (within-strain dose contrast) so between-strain latent tolerance can't
     confound M→Y; if NDE/NIE flip sign or magnitude, report the discrepancy.
   - *Replication-vs-interference diagnostic* (protects #3): code multi-object
     wells; fit `Y ~ (1|plate/well/object_label)` and compare within-well vs
     across-plate spread for the endpoint features. Drop true duplicates and
     report how much the mediation CIs tighten/spread when the effective n is
     corrected.
   - *Placebo dose*: verify no residual run-cadence gradient by splitting on `run`
     and confirming the NDE/NIE ordering is stable across runs.

5. **Deliverable.** One table + figure: stacked bar of NDE / NIE_via_extent /
   NIE_via_pixel vs total Cu effect on a*/b*/chroma, per species (tolerant
   R. mucilaginosa vs sensitive R. diobovata as bookends), with the censoring-
   corrected and cluster-corrected CIs.

## Expected Insights
- A direct answer: **"Cu alters pigment per cell" vs "Cu just changes how much
  colony you read."** Concretely, if NIE_via_extent ≈ total effect, the chromatic
  phenotype is a screening of the density effect — the "color space" story is
  mostly geometry; if a large NDE survives (especially at 5 mM where extent is
  barely affected), Cu genuinely modulates carotenogenesis, and per-colony
  a*/b* CoeffVar (pigment packing) is the phenotype worth mining next.
- Magnitude of the collider bias from Cu-dependent dropout on endpoint color
  (likely concentrated in sensitive species, where the survivor subset is most
  unrepresentative — a caveat that should attach to *every* existing endpoint
  table at 25–30 mM).
- A quick ruling on whether multi-object wells inflate replicate confidence
  (the mis-replication question) that all future endpoint analyses can adopt
  as a default filter.
- Whether tolerant and sensitive species differ in the *pathway* (decomposition
  type) rather than just the magnitude (total effect) — a sharper statement than
  the existing log2-extent ranking.

## Feasibility
- **Effort**: Low–Medium (one new table assembled from `v_phenotype` + existing
  endpoint and growth-rate outputs; g-computation in base R; bootstrap over 112
  plates; IPW via `glm`).
- **Data ready**: Mostly — the endpoint rows and per-colony pixel-stat columns
  already exist; needs the ~90–110 h window keyed on `hours_since_plate_start`
  and the well-adjacency/multi-object coding (small preprocessing).
- **Methods available**: Standard tools (sequential regression / g-computation,
  IPW, cluster bootstrap in R/pixi); no research-grade machinery required.
- **Key risk**: sequential ignorability for the mediator fails if strain-level
  latent tolerance drives both pigment packing and extent — the within-strain
  fixed-effect check will reveal it, but if it does, the NDE/NIE split cannot be
  trusted beyond a descriptive statement (the ambitious idea is then the honest
  upgrade path). Censoring+mediation interaction is the second fragile point.
