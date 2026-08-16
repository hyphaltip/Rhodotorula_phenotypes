# Temporal-Dynamics Persona — Two Ideas

Persona folder: `2026-08-15-color-phenotype-space`
File id: `08-temporal-dynamics`

---

# Idea 1 (Ambitious) — Pigment Ontogeny Atlas: Full trajectory reconstruction, functional clustering, and the area–color coupling over time

## Persona
**Temporal-Dynamics Specialist** — thinks in trajectories, not snapshots: which cells/colonies change *when*, *how fast*, and *whether they keep changing* after the biomass curve flattens.

## Motivation
Every existing analysis in this project reads a colony at one or a few hours: endpoint CIELAB medians (control-late-timepoint windows 70–80/80–90/90–110 h), endpoint `L* ~ rate_area × Cu` interaction, and extent/t50/saturation summaries. But a colony time-lapse with `hours_since_plate_start` spanning **0–117 h** and ~20 imaging passes per colony window is a *movie*, and carotenoid biology is inherently a timing phenotype — pigment is synthesized through the growth cycle, and industrially it matters *when* it accumulates (Rhodotorula are beta-carotene/torulene producers). The central temporal question nobody has asked yet: **do strains differ more in WHEN they pigment than in their final color, and how does Cu move each strain on the timing-vs-level plane?** This is the "pigment ontogeny" of the screen.

## Connection to Existing Data
- The full time-course lives in `v_phenotype`: keyed by `strain_id`, `strain_code`, `species`, `run_number`, `plate_number`, `well_position`, `object_label`, `image_name`, `hours_since_plate_start` (0–117 h), `copper_mm`, `imaged_at`; per row we have `Shape_Area` plus per-colony `ColorLab_L*Median`, `ColorLab_a*Median`, `ColorLab_b*Median` (and `StdDev`/`CoeffVar` for intra-colony spread).
- The series-construction pattern already exists and is reproducible: `analysis/growth_rates/scripts/00_build_series.R` (runs 353–356, `strain_id IS NOT NULL`, `hours_since_plate_start > 0`, largest object per plate/well/time, `plate_id = run_plate`) gives per-colony time series with up to 20 timepoints.
- The growth-rates analysis documents the two cruxes this idea must engineer around: (a) **the phase artifact** — the rising peak 6 h slope under Cu was shown to be a length-of-observation artifact (doubling time flat ~37 h across 0–30 mM, `rate_area` vs doubling time τ = 0.002); (b) **the plateau/floor limits** — `ln(area)` rarely plateaus in-window (asymptote extrapolated in ~82% of fits, Gompertz lags structurally unidentifiable in ~94%), while `Intensity`/`L*` *does* saturate (~48–54 h). So pigment curves are more "completable" than area curves; area curve shape must be treated cautiously.
- The imaging cadence grounding: the rig images on two interleaved ~3 h cadences, so a fixed hour drops ~85 strains; the control table's **rounding trick** (round `hours_since_plate_start` to the nearest integer "pass"; passes land at e.g. 75/78, 87/90, 105/108 h) is the alignment primitive to generalize for *all* passes, not just the late windows.
- Two growth axes already shown decorrelated (`rate_area` vs `rate_int` r = 0.21) — a hint that color-rate and extent-rate trajectories can be decoupled, which is exactly what an area–color coupling analysis quantifies.

## Approach
1. **Build per-colony long series.** Generalize `00_build_series.R` to retain *all* timepoints (not just peak windows): one row per (colony = `run_number` × `well_position`; Cu dose fixed at plate level) per pass, with `t = hours_since_plate_start`, `Shape_Area`, `ColorLab_{L,a,b}*Median`, and `ColorLab_{a*,b*}*StdDev`/`CoeffVar`. Round `t` to the nearest integer pass (the control-table trick) and record both raw `t` and rounded `pass_h`.
2. **Align onto a common grid.** Two per-platform options, run in parallel and compared: (i) *imputation/interpolation* — for each colony, monotonically interpolate each trait onto a common grid of rounded pass hours (e.g., every 3 h from 3 to 117 using PCHIP over the colony's own passes), so both cadences contribute and strains are comparable at identical `pass_h`; (ii) *model-based* — fit a small nonlinear mixed model per trait, `trait ~ f(t; theta) + (1|strain_code) + (1|colony)`, on raw irregular times, so cadence offset is absorbed by the model rather than by rounding. Prefer (ii) for inference, (i) for visualization/alignment of curves.
3. **Extract an onset/rate/level signature per strain × Cu.** For each of `a*` (red=carotenoid), `b*` (yellow), `L*` (darkening; recall `L*` ⇄ `Intensity_MeanIntensity` at r > 0.99), and `log(Shape_Area)`, define per-colony descriptors computed on *changes* to defeat the phase artifact: **onset time** `t_onset` = first pass where the trait exceeds colony baseline + k·(background SD resampled from early passes); **maximal sustained slope** `rate` over a fixed post-onset window (not the global 6 h max that produced the artifact); **saturation level** `level` and **time to saturation** `t_sat`. Aggregate per strain × Cu (median over replicate colonies; ~3–4 per strain–Cu). This produces a compact trajectory fingerprint per (strain, Cu): `(t_onset, rate, level)` for each of L*, a*, b*, area.
4. **Trajectory typing and the timing-vs-level plane.** Standardize the signature vectors across strains (isolate timing vs level axes via PCA of the fingerprint matrix, or functional PCA on the aligned `a*`/`b*` curves). Cluster strains into trajectory "morphs" (k-means on fingerprints; silhouette-validated). Then plot per morph, per Cu: timing axis (`t_onset`) on x, level axis (endpoint `a*`) on y — ask *how Cu moves each strain on that plane* (delays onset vs suppresses level vs both), and cross-tabulate morph membership against `species`, `origin`, `environment`, and the known `log2_extent` Cu-sensitivity ranks.
5. **Area–color coupling over time (does pigment accumulate after growth arrest?).** Per colony, compute a time-shifted coupling: cross-correlation between `d(log Shape_Area)/dt` and `d(a*)/dt`; and a **pigment-after-arrest index** = rate of `a*`/`b*` rise over the window after the colony's own area-saturation pass (using the saturation/missingness behavior from `04_doubling_time.R`: at 0 mM 95% of colonies saturate by ~72 h, at 30 mM only ~74%). Test whether high-Cu-sensitive strains pigment *late but still* vs *never*, and whether pigment rate couples to the area axis (`rate_area`) or the intensity axis (`rate_int`) per strain.

## Expected Insights
- A **per-strain × Cu "pigment ontogeny" atlas**: onset/rate/level triplets for L*, a*, b*, area — the first view of *when* each strain colors up, not just how dark it ends.
- Whether strain diversity in color is dominated by **timing** (t_onset, t_sat) or **level** (saturation `a*`) — e.g., strains might fall into "early-rapid-pigmenting" vs "late-slow" morphs that look similar at 90 h but are biologically distinct (→ production relevance: harvest-timing).
- How **Cu shifts each strain on the timing-vs-level plane**: does the documented level response (a*/b* drop with Cu; `a*` F ~ 1100 on endpoint) actually arise as a *delay* (onset pushed past the observation window / past growth arrest) or a true *level suppression*? Distinguishing these has a direct industrial implication (whether low-Cu media only buys time, not more pigment).
- Whether pigment is **growth-coupled or arrest-coupled**: if the pigment-after-arrest index is high, pigment accumulates on a schedule independent of extent (matching `rate_area` vs `rate_int` r = 0.21); if onset tracks area saturation, pigment is growth-phase-locked and the Cu timing shift is a downstream phase effect.

## Feasibility
- **Effort**: High
- **Data ready**: Mostly — `v_phenotype` has all timepoints; the series builder exists (`00_build_series.R`); needs new aggregation scripts and curve-fitting infrastructure (DuckDB query + R/Python NLME).
- **Methods available**: Needs custom implementation — nonlinear mixed models / functional PCA are research-grade, and the identifiability traps are documented (unidentifiable lags, extrapolated area asymptotes) so the model-based leg carries real risk.
- **Key risk**: Sparse irregular sampling (only ~19–20 passes per colony; two cadences ±3 h) plus high-Cu dropout (~15–18% of strain-plates never reach late timepoints) makes per-(strain × Cu) trajectory fitting underpowered, and the same **phase-artifact trap** that caught the growth-rate headline could corrupt any raw slope under Cu — mitigated only by change-based descriptors and the doubling-time-style validation (e.g., a no-timing-shift null where Cu is shuffled within species).

---

# Idea 2 (Immediately Tractable) — Onset, not endpoint: when does each strain pigment, and how does copper move that moment?

## Persona
**Temporal-Dynamics Specialist** — the minimal temporal claim that costs one table and changes every downstream reading of the color space: replace "final color" with "color-at-a-fixed *phase*" and "onset time."

## Motivation
The cheapest temporal win is a discrete, hypothesis-testable summary that sidesteps curve fitting entirely: **pigment onset time**. The dataset already forces a pass-alignment discipline (two interleaved ~3 h cadences), so onset can be defined cleanly as "first rounded pass where the chromatic median crosses a per-colony noise threshold." This turns the color phenotype from a number (endpoint `a*`) into a moment (when red/yellow appeared) — and Cu's effect on color can then be re-expressed as **which strains get their pigment delayed by Cu, by how many hours, before we even talk about amplitude**. It directly serves the industrial question (carotenoid accumulation timing) and is buildable from the control-table's own rounding machinery.

## Connection to Existing Data
- Same `v_phenotype` rows: per colony per image, `hours_since_plate_start` (0–117 h), `copper_mm`, `run_number`, `well_position`, `ColorLab_a*Median`, `ColorLab_b*Median`, `ColorLab_L*Median`, `Shape_Area` — and the intra-colony spread columns (`ColorLab_a*CoeffVar`, etc.) for the noise-threshold baseline.
- **The rounding trick** is already implemented and validated in `analysis/control_late_timepoint_phenotype/scripts/build_phenotype_table.py`: rounds `hours_since_plate_start` to the nearest integer pass (passes at 75/78, 87/90, 105/108 h; 314/320 strains covered, spot-checked on strain 185). Generalizing it to *all* passes (0–117 h) is a small edit, not new machinery.
- The phase-artifact caution is quantified: `growth_rates/04_doubling_time.R` shows raw slopes under Cu are unreliable (peak-slope rises with Cu yet doubling time is flat ~37 h); so this idea deliberately uses **threshold-crossing times**, not slopes, as the primary readout.
- Plateau/floor grounding: `Intensity`/`L*` saturates ~48–54 h (intensity axis is "completable") while `Shape_Area` rarely plateaus in-window — so onset is far better defined on `a*`/`b*`/`L*` than on area, and area enters only as a co-variate/phase-marker.

## Approach
1. **Rebuild all-pass series for the chromatic axes.** Extend the control-table's rounding filter to runs 353–356 at every pass (not just [70,110]; include early/late): per (colony, `pass_h`) keep `ColorLab_{a*,b*}*Median`, `ColorLab_a*CoeffVar`, `Shape_Area`. Keep both `t` and rounded `pass_h`.
2. **Define per-colony onset.** For each colony, use the early passes (passes where `a*` is at its stable baseline, e.g., first 3–4 passes) to get background median and SD of `a*`; define `t_onset(a*)` = first `pass_h` where `Median` exceeds background + 3× background SD *and* stays above it for the next 2 passes (hysteresis to kill noise). Same for `b*` and (as a sanity, since `L*`⇄Intensity) `L*`. Record `t_onset`, `t_sat` (first pass where `a*` is within 95% of its plateau), and the colony's final `a*`.
3. **Per-strain × Cu aggregation with cadence fairness.** Median `t_onset` and fraction-pigmenting across replicate colonies (≥3 per strain–Cu, flagged otherwise — same rule as the control table's `n_colonies < 3` flag). Use the rounded `pass_h` and, where a strain sits on the alternate cadence (3 h off), bin to the *other* strain's pass grid ±3 h or explicitly record both—the ±3 h is below the granularity of the question (we care about 10–30 h shifts, not 3 h).
4. **Test the timing shift.** Mixed model `log(t_onset(a*)) ~ factor(copper_mm) + (1|species/strain_code)` on the well-sampled subset (like `02_species_cu_rates.R`), plus a **Cu × species interaction on onset** — the temporal analogue of the extent analysis (`F = 4.34, p ~ 6e-4` on log max_area). Report effect sizes in hours and 95% CIs, and confirm against a permutation null (shuffle Cu labels within strain) so the phase artifact is ruled out by construction.
5. **Contrast onset-vs-level.** Join `t_onset` to the existing endpoint/level readouts (`control_late_timepoint_phenotype` tables; species-sensitivity `log2_extent`). Rank species both ways (by onset delay under Cu, by level loss under Cu) and ask whether the two rankings disagree — i.e., species that keep pigmenting (late onset is fine) vs species that stop (onset never reached).

## Expected Insights
- A **per-strain, per-Cu pigment-onset time in hours** — the first *when* answer in the project; whether Cu delays pigment onset monotonically (e.g., +5–30 h as 0→30 mM) or abolishes it outright for sensitive taxa (`R. kratochvilovae`, `R. araucariae`, `R. taiwanensis`).
- **Whether species differ more in onset time than in final color** — the strain-level temporal axis orthogonal to the endpoint tables; joint onset×level ranking may reshuffle the current species color/order story.
- A **Cu × species interaction on onset timing** comparable to the known interaction on extent — showing whether Cu sensitivity is a *growth* phenomenon, a *color-timing* phenomenon, or both.
- A production-relevant rule of thumb: for tolerant strains (e.g., `R. mucilaginosa`), optimal pigment harvest occurs at a *time* that shifts with Cu even when final pigment is similar — the onset table is the argument.

## Feasibility
- **Effort**: Low–Medium
- **Data ready**: Yes — all inputs are already in `v_phenotype`; the rounding/aggregation pattern ships in `build_phenotype_table.py` and the series builder in `00_build_series.R`; no new imaging or features.
- **Methods available**: Standard tools — mixed models (as in `02_species_cu_rates.R`), threshold-crossing, permutation null; only the onset rule is mildly custom.
- **Key risk**: Per-strain × Cu replicate colonies are few (~3–4), so onset medians are noisy; and a strain that *never* pigments is censored (no `t_onset`) — must be treated as an informative "onset > 117 h" outcome via survival-style analysis, not dropped, or the high-Cu timing shift will be underestimated.

---

## Persona meta-block (for the session harness)
- **Persona name**: Temporal-Dynamics / Time-Series Specialist
- **Idea 1 title**: Pigment Ontogeny Atlas: Full trajectory reconstruction, functional clustering, and the area–color coupling over time
- **Idea 1 feasibility**: High effort; series data is ready but irregular-cadence + high-Cu-missingness makes per-(strain × Cu) curve fitting statistically tight — the same phase-artifact trap that caught the growth-rate headline is the main thing that could kill it.
- **Idea 2 title**: Onset, not endpoint: when does each strain pigment, and how does copper move that moment?
- **Idea 2 feasibility**: Low–Medium effort; buildable immediately from the existing rounding/series machinery with standard mixed models; main risk is censoring — strains that never cross the pigment threshold must be modeled as onset-∞ (survival), not dropped.
