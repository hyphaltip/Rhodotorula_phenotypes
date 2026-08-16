# Idea 1 — Chroma-Weighted Hue Space: Separating Pigment *Identity* from Pigment *Amount* from Tissue *Density*, and the Origin of Intra-Colony Heterogeneity

## Persona
**Color/Imaging Scientist** — A colorimetric-rigor lens: RGB/CIELAB/HSV interrelationships, hue as a circular variate that must never be arithmetically averaged, chroma as pigment amount, L* as tissue density, and pixel-histogram *shape* (Min/Q1/Median/Q3/Max, StdDev, CoeffVar) as a phenotype in its own right.

## Motivation
Every downstream analysis so far treats L*, a*, b* as three interchangeable scalar readouts and reports their *medians*. But colorimetry tells us these are not three independent axes: they decompose into **hue angle** (what pigment — identity), **chroma/saturation** (how much pigment — amount), and **lightness/density** (how much biomass per pixel, a growth axis, since L*~Intensity r≥0.99 is already documented). Mixing these three in a single median table conflates "this strain makes a different carotene" with "this strain makes more carotene" with "this colony is denser". The open direction the dataset keeps pointing at — "color space / hue vs chroma vs lightness / discrete pigment morphs" — is precisely the question a color scientist is trained not to answer with arithmetic means. And there is a distinct, unused phenotype sitting in the *shape* of the per-colony pixel histogram: **intra-colony heterogeneity**, which I can already see rising with copper stress (a*CoeffVar: 0.31 @ 0 mM → 0.54 @ 30 mM in R. mucilaginosa at 90–100 h), i.e. heterogeneity looks biologically structured, not noise.

## Connection to Existing Data
- `v_phenotype` gives per-colony **ColorLab** L*/a*/b* with Min, Q1, Mean, Median, Q3, Max, StdDev, CoeffVar; **ColorLab_ChromaEstimatedMean/Median**; **ColorHSV** Hue/Saturation/Brightness (same 8 stats); and **Colorxy** x,y chromaticity (same 8 stats) — all fully populated (verified: COUNT = COUNT(Colorxy_xMedian) = COUNT(ColorHSV_SaturationMedian) = COUNT(ColorLab_ChromaEstimatedMedian) = 176,269 rows). These are the **~96 color-descriptor columns that all prior analyses ignore** (they only use ColorLab_{L,a,b}*Mean/Median).
- Live query shows the core colorimetric fact: species differ ~4× in chroma (Color Lab_ChromaEstimatedMedian 2.93 for R. sp. clade XIII → 9.82–10.4 for dairenensis/evergladensis/paludigena at 70–80 h) while L*Median is pinned at ≈70–71 for *all* species — i.e. L* carries no species signal at this window but chroma and hue angle (atan2 of a*/b* medians, ≈23–59° by species) do. This is a concrete, already-visible answer to "a*/b* after controlling L*".
- Copper response is amount-centric at late time and trajectory-with-chroma-collapse early: at 90–100 h chroma drops 10.3 → 5.9 with Cu and a*CoeffVar rises — pure amount + heterogeneity shift with near-constant hue; but at t≈0 the same strain shows hue flipping to ≈−160° at high Cu because chroma collapses into the noise floor (a*≈b*≈0 makes hue undefined).
- Known confounds to build in (from existing analyses): L*~Intensity collinearity (0.99+), plate identity nested within Cu, 5 imager runs / 8 control plates for run-level normalization, run 357 exclusion precedent.

## Approach
1. Derive the colorimetric triple from `v_phenotype` in DuckDB. Per colony compute chroma `C = "ColorLab_ChromaEstimatedMedian"` (amount), CIELAB hue `h° = (atan2("ColorLab_b*Median", "ColorLab_a*Median") * 180/pi()) % 360.0` folded to the non-negative half-open unit circle (identity), and density `D = "ColorLab_L*Median"` (keep L* for amount×density coupling since C*L* ≈ α-carotene *optical* density logic, but *always* partition D out when examining C). Add heterogeneity descriptors: `H_aRange = "ColorLab_a*Q3" - "ColorLab_a*Q1"`, `H_CRange = "ColorLab_L*Q3" - "ColorLab_L*Q1"`, `CV_a = "ColorLab_a*CoeffVar"`, `Hue_Disp = 1 - R` (circular resultant length across the colony's 8 hue stats is NOT derivable pixel-wise — use across-colony dispersion instead).
2. **Chroma-weighted circular statistics** (the rigor the data demands): never average hue arithmetically. For each strain × Cu group, compute resultant vector:
   `R̄x = SUM(C·cos(h°))`, `R̄y = SUM(C·sin(h°))`, then circular mean `μ = atan2(R̄y,R̄x)` and concentration `R̄ = hypot(R̄x,R̄y)/SUM(C)` in DuckDB aggregates. Low-chroma colonies (C < noise floor ≈ 1, and/or where max(a*,b*) < sensor floor) enter with weight ≈ 0, automatically killing the −160° t≈0 artifact instead of contaminating the mean. Filter t < 12 h separately (imperfect segmentation at colony onset).
3. Ask the **pigment-identity vs pigment-amount** question with mixed-effects models in R (as in `growth_rates/scripts/03_color_interaction.R`):
   - `str"ColorLab_a*Median" ~ Cu + (1|species/strain) ` already trending in existing code — extend to the *angle* via circular-linear regression of `h°` on Cu **weighted by C**, and to the chromium axis: `log(ChromaEstimatedMedian) ~ Cu × species + (1|strain)` to test whether Cu acts by *amount* (chroma slope, identity constant) or by *identity* (hue–Cu interaction). Species-level hue existing spread (≈23–59°) sets the null magnitude.
   - `log(a*CoeffVar) ~ Cu × species + (1|strain)`: is stress-induced intra-colony heterogeneity generic or species-specific (visible: a*CV 0.31→0.54 with Cu)?
4. Test **discrete pigment morphs**: cluster strain-Cu-mean hue angles on the circle (von Mises mixture / circular k-means with 2–6 components) and compare against species, origin, environment, and the existing extent-based Cu sensitivity index (−0.72 to −4.13, `growth_rates`). Quantify redundancy of hand-derived axes vs the ~96 raw color columns via variance-explained / mutual-information (open direction: "how much of the color space is independent").
5. Validate adversarially: (a) confirm L*, Interior Q1→Q3 core-vs-edge contrast, and TextureGray_* do not duplicate the hue axis (corrmatrix; partial out L* and replot); (b) permutation/nulling: shuffle strain labels within Cu and Cu within plate to assert the chroma/hue/CV effects exceed between-plate spread (plate-position analysis shows within-Cu plate variance is small, 0.1–4.7%, for L* — recheck for chroma); (c) replicate the whole pipeline using `ColorHSV_HueMedian` to expose whether the HSV pathway is degenerate (see Key risk).

## Expected Insights
- Whether copper stress is an *amount* phenotype (chroma drops, hue stable — the late-time data already hints yes at 90–100 h) or an *identity* phenotype (hue angle rotates), and whether that differs by species — a direct read on whether the screen measures carotene *yield* (industrial relevance) or carotene *quality*.
- Whether strains form discrete circular pigment morphs (e.g. 3–4 hue clusters) that map onto species/elades, vs a smooth continuum — intervening variable is chroma, not L*.
- Whether within-colony heterogeneity (a*CoeffVar rising 0.31→0.54 with Cu) is a *stress-response phenotype* independent of both area and mean hue — a candidate second axis beyond the growth axis (recall rate_area vs rate_int only r=0.21: color-space axes may be similarly independent).
- Quantitative answer to "does a*/b* carry signal after L*": we predict yes, because chroma/hue are near-orthogonal to the L*~Intensity axis by construction, but the analysis will state how much variance survives.

## Feasibility
- **Effort**: High
- **Data ready**: Yes — 176,269 rows, all HSV/xy/chroma/Texture columns populated; only light preprocessing (chroma floor + hue circular weights + t<12 h filter).
- **Methods available**: Standard — DuckDB (built-in math aggregates for circular resultants), R `circular`/`circular`-compatible lm for circular-linear, `movMF` for von Mises mixtures; all reproducible under pixi.
- **Key risk**: (1) The HSV-derived Hue is likely **degenerate** in this pipeline (verified live: ColorHSV_HueMedian ≈ 0.1 for essentially every species — CellProfiler's HSV conversion or storage wraps the angle unusably), so the HSV pathway may contribute nothing and the circular machinery must be run on the CIELAB-derived hue; (2) chroma in low-pigment strains at high Cu sits in the sensor noise floor where hue is undefined — mitigable only by the chroma-floor weighting, which is exactly the design above; (3) mixing of plate/run "identity" with Cu (nested design) needs the reference-centering from Idea 2.

---

# Idea 2 — Pigment Appearance Time: Temporal Chroma/Hue Onset Trajectories with a Within-Run Color Calibration

## Persona
**Color/Imaging Scientist** — Same colorimetric-rigor lens, but applied to the *time domain*: a color scientist treats the pixel not as a static value but as a signal with an **onset** (when pigment crosses the detection floor), a **saturation rise** (chroma ramp), and a **stable hue lobe** (identity confirmed) — and knows that any multi-run acquisition needs a color-normalization pass before trajectories are comparable (commonest silent killer in imaging screens).

## Motivation
The existing temporal analysis only follows **L*Mean over time** (`growth_rates`, `explore_plate_position`/05_variance + 06_growth_curves), which collapses the pigment story: L* is collinear with Intensity (r≥0.99), so L*-curves are *density* curves, not *pigment* curves. The untold story is **when pigment appears** — i.e. chroma/hue onset. Live data shows the raw material: at Cu=0, R. mucilaginosa chroma ramps 2.1 → 4.97 → 6.22 → 7.3 → 8.17 across timepoints while hue stabilizes ≈40–80° (yellow-orange carotene sector) — a clean, slow, extractable onset. And the b* pathway is *non-monotonic* in Cu (median b* 2.67 @0 / 7.7 @5 / 7.3 @10 / … / 4.9 @30 at 90–100 h) — a hue-chroma trajectory that a median-endpoint table structurally cannot see. Finally, with 5 imager runs and replicate plates, run-to-run illumination drift will smear any onset estimate unless reference-normalized — a small, data-grounded calibration is possible right now.

## Connection to Existing Data
- Temporal grain is already supported: `hours_since_plate_start` (0–117 h) on every `v_phenotype` row, plus `run_number` (5 runs), `plate_number` (112 exp + 8 control), `is_control`, and strain keying. R. mucilaginosa dominates the strain set (9,614 colony-rows in a single 70–80 h window) giving ~24+ replicate colonies/cu for stable per-strain–Cu trajectories.
- A control-plate + common-reference normalization is grounded because: (a) plate identity is nested within Cu (explore-plate-position finding), so a naive across-runs pass would alias run with Cu; (b) 8 `is_control` plates exist across the 5 runs to estimate run-level bias; (c) the same reference strain/species recurs in every run (mucliaginosa), enabling a per-run shift/scale anchor.
- Duplicates nothing existing: prior time-series used only `Shape_Area` and L*/intensity Means; no one has estimated chroma/hue **onset** times or run-level color offsets.

## Approach
1. **Build the per-strain–Cu color trajectory** from `v_phenotype` in DuckDB (mirroring `growth_rates/scripts/00_build_series.R`):
   ```
   SELECT strain_id, strain_code, species,
          CAST(factors['Copper concentration'] AS DOUBLE) AS copper_mm,
          hours_since_plate_start,
          "ColorLab_ChromaEstimatedMedian"           AS chroma,
          atan2("ColorLab_b*Median", "ColorLab_a*Median")*180/pi() mod 360 AS hue,
          "ColorLab_L*Median"                        AS L,
          "ColorLab_a*CoeffVar"                      AS aCV,
          "Colorxy_xMedian", "Colorxy_yMedian"       AS chromaticity
   FROM v_phenotype
   WHERE species='Rhodotorula mucilaginosa'  -- and any subset enabled by the panel
     AND run_number != 357  -- existing exclusion precedent
   ```
2. **Run-level color calibration** (reference-centering, not whitening): for each `(run_number, plate_number)` estimate the drift of the reference species' (or control plate's) `(chroma, hue, L)` against its global median; then per-colony residualize: `chroma_adj = chroma - run_ref_delta_chroma`, `hue_adj = hue_adj rotated so run_ref hue maps to 0` (circular shift). This removes the run/plate offset that the nested design otherwise smuggles into every Cu comparison, and makes the 5 imager runs a *strength* (5 independent bias estimates) instead of a confound.
3. **Onset/extract per-colony timing phenotypes** on the adjusted series — per Colony (well_position × Cu × strain), define in DuckDB window functions:
   - `t_chroma_onset` = first hours_since_plate_start where `chroma_adj > noise_floor` (floor = 2× median chroma of that strain's t<12 h window, i.e. the pigment-free colony margin) for 2 consecutive images — a "pigment appearance time".
   - `t_darkening` = crossing of first derivative of L* (peak negative slope of LOWESS-smoothed L*) — the density/checkpoint onset.
   - `t_hue_stable` = first time the chroma-weighted circular angle enters and *stays within* ±30° of the strain's own late (70–110 h) stable hue lobe (protects against the −160° artifact — only trust hue where chroma is above floor, the idea-1 rule).
   - Summarize as per-strain–Cu medians (matching the existing 70–80/80–90/90–110 window convention).
4. **Model the timing phenotypes**: `t_chroma_onset ~ Cu × species + (1|strain)` and `t_hue_stable ~ Cu + (1|strain)` in R (lme4, reuse the pipeline skeleton); test whether onset *precedes* or *follows* darkening (`t_chroma_onset` vs `t_darkening` paired within colony) — is pigment induced before or after growth stalls, and does Cu advance or delay onset?
5. Validate: (a) permutation of the reference run-shifts to show calibration changes effect estimates materially (i.e. that not calibrating would have changed conclusions); (b) cross-check onset with TextureGray_* scale-05 (Entropy/Contrast across angles) as an independent pigment-texture readout of surface maturing; (c) compare `t_chroma_onset` to the extent-based sensitivity index — is Cu-sensitive species' pigment simply *slower*, or *never*?

## Expected Insights
- A new, easy-to-communicate per-strain phenotype: **pigment appearance time** (e.g., "xx hours until first chroma") and darkening onset — orthogonal to both growth rate (extent) and endpoint color (yield), and directly relevant to BRET/carotenoid timing questions (when is the culture worth harvesting).
- Whether copper *delays* pigment onset (toxicity) or *cuts* final chroma (yield) while keeping onset timing — discriminating the two mechanisms Cu screens invariably conflate; the current endpoint median tables cannot separate them, but onset + late chroma jointly can.
- A run-level calibration constant table (5 runs × reference deltas) that benefits *every* subsequent color analysis in this project — reusable normalization, not a dead-end analysis.

## Feasibility
- **Effort**: Low–Medium
- **Data ready**: Yes — full 0–117 h time-course per colony in `v_phenotype`; control plates and repeated reference strain in every run; only aggregation + window math needed.
- **Methods available**: Standard — DuckDB window functions + R lme4/circular; no new measurement required.
- **Key risk**: (1) Time sampling may be coarse/irregular at late hours (0–117 h with variable cadence), so "onset" resolution is limited by image cadence — mitigate by reporting onset in floor-of-bin hours and requiring 2 consecutive above-floor images; (2) colonies are undetected/small at early t (segmentation onset), so early-time chroma is measured on tiny pixel sets — must require a minimum Shape_Area before trusting a timepoint; (3) if HSV Hue is degenerate (idea-1 caveat), the `t_hue_stable` score must lean on the CIELAB-derived hue only.
