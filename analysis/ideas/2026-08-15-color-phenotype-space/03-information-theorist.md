# Idea 1 (AMBITIOUS): The Minimal Sufficient Statistic of the Phenotype Screen

## Persona
**Information-Theorist** — I think in entropy, mutual information, redundancy, and compression. The ~200-column feature vector per colony is a *code*; I want to know how many independent bits it genuinely carries, how redundant it is, and what the smallest feature set is that preserves all the information this screen can possibly contain about the biological labels.

## Motivation
Every screen with a fat feature table is really asking a compression question: how much of what we measure is new information vs. a deterministic transform of something already measured? This dataset is almost purpose-built for that. L* is ~>0.99-collinear with Intensity (two names, one axis). The 24 ColorLab columns are eight statistics (Min,Q1,Mean,Median,Q3,Max,StdDev,CoeffVar) of just three channels — order statistics of one histogram, so highly dependent. TextureGray has the same 13 Haralick features at 4 angles plus angle-average (near-rotational copies). A correlation matrix will *tell* you this; only an entropy/MI decomposition tells you the **information-theoretic cost**: how many independent bits are left after the redundancy, and what a minimal sufficient feature set — the "skeleton" that loses no species information, no copper information, and no batch confound — looks like. This directly answers the user's "how much of the ~200-col feature space is redundant and what are the few independent axes" direction, quantitively.

## Connection to Existing Data
- `v_phenotype` gives the full ~200-col vector per colony-timepoint: Morphology (Shape_Area, Circularity, Solidity, Eccentricity...), Intensity, ColorLab{Min..CoeffVar} for L*/a*/b*, ColorHSV{...} incl. Hue, Colorxy{...} incl. x/y, ChromaEstimatedMedian, TextureGray (13 features × 4 angles + avg).
- `strain` was metadata gives species (3 genera, ~10 spp), genus, origin (China 180 / Italy 119 / America 9), environment (marsh_tidalflat 162, soil 45, ...); `copper_mm`, `run_number`, `plate_number`, `well_position`, `hours_since_plate_start` key every row.
- growth-rates already established L* vs Intensity collinearity (>0.99) — the canonical redundancy seed; plausibly the *limit* of redundancy, not the whole story.
- explore-plate-position established plate identity is nested in Cu and contributes small but nonzero variance — the natural seed for the batch-MI half of this analysis.

## Approach
1. **Collapse to per-colony endpoint objects.** Use the same selection/dedup as growth-rates (runs 353-356, strain_id NOT NULL, hours>0, largest object per plate/well/time, plate_id = run_plate) and take the window 70–80 h; aggregate each of the ~200 features per colony (median across timepoints for location stats). Treat Hue circularly (embed as sin/cos of hue angle — never arithmetic-mean the circular axis). This yields ~8k colonies × ~200 features; **explicitly treat this as the unit** so timepoint pseudo-replication cannot inflate MI.
2. **Entropy/MI estimation.** Standardize continuous features; estimate pairwise I(Xi;Xj) and I(Xi;Species) with the continuous k-NN (Kraskov–Stögbauer–Grassberger) estimator, and cross-check with a discrete-binned version (per-feature quantile bins, e.g. 5–10 levels) for interpretability and bias robustness. Report bootstrap CIs; this is the parameter-sensitivity sweep the robust-analysis conventions demand.
3. **Redundancy structure.** Build the N×N mutual-information matrix; cluster features by information overlap (hierarchical clustering on 1−normalized-I). Identify the near-deterministic families (the L*–Intensity pair; the Min..Max order-statistic chains; the 4-angle texture copies). Report the estimated **effective rank of the code**: how many roughly independent bits total, H(full feature vector), vs. H of the L*a*b* median subset.
4. **Minimal sufficient statistic.** Forward-greedy feature selection that maximizes I(feature set; Species) while discounting features already covered (conditional-MI / MIFS criterion: add Xj maximizing I(feature_set ∪ Xj; Species) gain per added bit). Report the smallest feature set recovering the full I(features; Species) and the fraction of that MI recoverable from shapes+colors+texture at each step. Repeat the target as Cu-tolerance (the log2 sensitivity index from species-Cu-sensitivity) and origin/environment to see whether different labels need different skeletons.
5. **Information yield of the ~170 unused columns.** Compute ΔI = I(full feature set; Species) − I(L*a*b* Median subset; Species) *after conditioning out batch*: the number of *net new* bits about species carried by IntraColony StdDev/CoeffVar, HSV/xy/chroma, TextureGray, and morphology. This is the headline number: do the unused columns pay rent?
6. **Confounding MI test.** Estimate I(features; run_number) and I(features; plate_number | run) to bound the batch-channel capacity that could be mis-credited to species; then compare I(features; Species) marginal vs I(features; Species | run). The gap is the spurious species-MI that plate/run batch could inject — the information-theoretic twin of the plate-position variance finding, now in bits rather than variance components.
7. **Validate:** permutation null (shuffle species labels *within plate/run strata* to preserve batch structure) for every headline MI; report only features whose MI exceeds the null's 95% quantile.

## Expected Insights
- A quantified redundancy budget: e.g., "the 200 columns carry ≈H bits of entropy; >85% of the species information is in ≤8 features; the L*a*b* medians alone capture ≈X% of it, and the remaining Y bits come from texture/heterogeneity." This is the minimal-sufficient-statistic answer the user asked for, in their own vocabulary.
- The single largest redundancy families ranked by wasted bits (likely: L*–Intensity, the order-statistic chains, and the 4-angle texture copies) — a concrete recommendation of which columns to keep in any downstream model/DB table.
- Whether batch (run/plate) injects enough MI to change species rankings — i.e., whether naive downstream MI/ML would be learning the imager, not the yeast.

## Feasibility
- **Effort**: High (entropy estimation, greedy selection, batch-MI, bootstrap CIs; needs custom code beyond off-the-shelf tables)
- **Data ready**: Yes — `v_phenotype` view and the growth-rates selection logic exist; only aggregation/estimation code is new
- **Methods available**: Needs custom implementation (KSG estimator, conditional-MI greedy selection) atop scipy/scikit-bio-style primitives — write it in R (existing stack) or Python
- **Key risk**: k-NN entropy/MI estimators are biased in 200 dims and with ~8k rows; the per-species imbalance (R. mucilaginosa dominant) and circular hue could inflate or deflate MI — the discrete/continuous cross-check and stratum-preserving permutations are the defense, but a wrong estimator choice could produce confidently wrong headline numbers.

---

# Idea 2 (TRACTABLE): Does Hue/Texture Carry Species Information that L* Does Not?

## Persona
**Information-Theorist** — edge cases: marginal MI can be large and still useless after conditioning. The existing analysis only ever looks at ColorLab {L*,a*,b*}.Median. I want the *conditional* and *partial-decomposition* view: how much of the species signal in the ~170 unused columns is **unique** (independent of the median-vector channel) vs. **redundant** (already carried by L*a*b* medians) vs. **synergistic** (only visible as a conjunction)? That is exactly the Williams–Beer partial-information-decomposition (PID) question, and it is the honest formulation of "do the unused columns add anything."

## Connection to Existing Data
- Specific candidate features from `v_phenotype`: `ColorHSV_Hue.Median` (circular — embed as sin/cos of hue angle), `ColorLab_a.StdDev` & `ColorLab_b.CoeffVar` (intra-colony chromatic heterogeneity), `TextureGray_Entropy_avg` / `TextureGray_Contrast_avg` / `TextureGray_AngularSecondMoment_avg` (pigment spatial arrangement), `Colorxy_x.Median`/`Colorxy_y.Median`, `ChromaEstimatedMedian`, alongside the already-used `ColorLab_L.Median`, `ColorLab_a.Median`, `ColorLab_b.Median`, `Shape_Area`.
- Note Hue ≈ atan2(b*, a*) in CIELAB — so a naive "Hue carries species info" claim may be pure redundancy with a*/b* medians; the interesting quantity is **I(Hue; Species | a*, b*)**.
- Labels available per colony: species (from strain), origin, environment (`strain` metadata), `copper_mm`, `run_number`, `plate_number`.
- Species-Cu-sensitivity analysis provides the discrete target list and confirms species-level structure exists to explain.

## Approach
1. **Build the endpoint table** (same selection/collapse as Idea 1: runs 353-356, 70–80 h window, per-colony median-of-timepoints), Cu=0 plates only so copper does not enter as a confound on first pass; ~8k colonies × a hand-picked ~10 features above.
2. **Discretize**: quantile-bin each continuous feature (5–7 levels); encode HueMedian as {sin(H), cos(H)} bins or a wrapped-hue discretization. This makes exact-marginal MI computable by contingency tables and gives interpretable "phenotype levels."
3. **Marginal + conditional MI vs Species.** For each candidate feature X, compute I(X; Species) and the conditional I(X; Species | Z) for Z = {L*.Median}, then Z = {L*,a*,b*}.Median, using discrete exact counts first, cross-checked with the continuous KSG estimator. Rank features by the *conditional* MI — that ranks "rent-paying" new axes. Interpret: hue conditional-on-a*/b* should be ≈0 (redundant transform); a*.StdDev and Texture Entropy may be the true additions.
4. **PID for the key pairs.** For inspired pairs (X1, X2) = (L*.Median, TextureGray_Entropy_avg), (a*.Median, a*.StdDev / b*.CoeffVar), (a*.Median, b*.Median), (a*.Median, Hue.Median-embedded), decompose I(X1,X2; Species) into {Unique1, Unique2, Redundant, Synergistic} via the Williams–Beer lattice on the discrete marginals. The Redundant term for (a*, b*) shows how much of the chromatic channel was already 'double-counted'; the Unique2 term for (median, StdDev) is the direct estimate of how much intra-colony heterogeneity reveals about species beyond its mean color.
5. **Null / batch control:** permutation null shuffling species within (plate, run) strata; recompute conditional MI conditioning out `run_number` (I(X; Species | L*, run)) to show the estimate is not an artifact of strain-×-plate layout.
6. **Extend target:** repeat the top-ranked features against origin and against Cu-tolerance class, since species is not the only label of interest.

## Expected Insights
- A concrete, ranked answer to "do the ~170 columns add information beyond L*a*b* medians?": likely *intra-colony heterogeneity* (a*/b* StdDev/CoeffVar) and texture entropy carry a small but significant **unique** species channel, whereas Hue is shown to be ~fully redundant with (a*, b*) — an information-theoretic proof that no new color axis is hiding in pure hue.
- The PID numbers quantify synergy too: e.g., "the conjunction (dark AND high-inside-texture) is what separates R. kratochvilovae" would be a Species×feature insight no median table can see.
- Feasibility demo: an I(X; Species | medians) ranking table and a small heatmap of the PID terms — immediately interpretable, directly reusable by any follow-on modeling.

## Feasibility
- **Effort**: Low–Medium (one endpoint table + contingency/MI library; PID for bivariate pairs is a small routine)
- **Data ready**: Yes — `v_phenotype` view and selection logic from growth-rates/control-late-timepoint exist
- **Methods available**: Mostly standard tools (scipy/MI, or a ~50-line PID implementation); conditional KSG needs a small custom routine
- **Key risk**: conditional-MI estimators are bias-prone on small conditional strata (rare species × rare hue bin), and the circular-hue discretization is easy to botch; the exact-discrete contingency cross-check plus stratum-preserving permutations are the main safeguards.
