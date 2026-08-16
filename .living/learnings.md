# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### [2026-08-15] L-1 — mycelium scripts require Python >= 3.10
- **Category**: tooling
- **What happened**: `install_convention.py` (and other mycelium core scripts) use PEP 604 union syntax (`Path | None`), which throws `TypeError: unsupported operand type(s) for |` under the system default `python3` (3.9.18 on this host).
- **Why it matters**: Home-page python3 is too old; running the scripts fails loudly but misleadingly during mycelium init.
- **Resolution**: Run with `/usr/bin/python3.12` (available) or under the repo pixi env; also pass `--network-dir` explicitly because the script's auto-detect may miss the plugin marketplace path.
- **mitigation_type**: structural — add a version guard/`sys.version_info` check or a shebang in the scripts.
- **structural_mitigation_candidate**: assert `sys.version_info >= (3,10)` with a clear message at script start.
- **Tags**: tooling, mycelium, python, setup

### [2026-08-15] L-2 — skillpacks/Autonomous-Science repo is no longer public
- **Category**: tooling
- **What happened**: `git clone https://github.com/arjunrajlaboratory/Autonomous-Science.git` fails (`Repository not found`) — the repo is absent from the org's public repos. Only `scientific-agent-skills` and `bioSkills` cloned successfully.
- **Why it matters**: The `skill-bridge` convention's `skill-sources.yaml` lists this as a verified persona source (`personas/library/arjun_raj.json`); the persona routing feature is therefore currently unavailable.
- **Resolution**: Left `skillpacks/` with the two working repos. TODO: re-point `skill-bridge` `skill-sources.yaml` to a live personas source or drop the source entry.
- **mitigation_type**: ambient-awareness — monitor whether the org restores/renames the repo.
- **Tags**: tooling, skill-bridge, skillpacks, third-party

### [2026-08-15] L-3 — dplyr mutate using `named_vec[[.$col[1]]]` broadcasts row 1's value to all rows
- **Category**: R / data-wrangling
- **What happened**: In `05_time_variance_partition.R` I relabeled a bound list's `.id` column with `mutate(trait = traits[[.$trait[1]]])`. Because `.[1]` takes the *first element of the whole tibble*, every row got the first trait's label ("L* (lightness)") instead of its own — 6 distinct model fits collapsed to 6 identical labels, and the `pivot_wider` in validation then threw "values not uniquely identified".
- **Why it matters**: This exact anti-pattern (indexing a lookup vector with `[[...]]` + `[1]` inside a single mutate) silently labels every row with the first group's name; it only surfaces later as duplicate keys. Recovery was easy here because the models had already been fit correctly.
- **Resolution**: Use vectorized name lookup `unname(traits[.$trait])` so each row maps through `traits[name]`. Validate by checking `n_distinct(trait)` or that a `pivot_wider` on the relabeled key doesn't warn.
- **mitigation_type**: structural — prefer `lookup[.$col]` over `lookup[[.$col[1]]]` when relabeling within a mutate.
- **Tags**: R, dplyr, mutate, pivot_wider, label-map, debugging

### [2026-08-15] L-4 — dplyr `summarise(across(...))` with `.names` throws "subscript out of bounds" when a group's lambda returns length-0
- **Category**: R / data-wrangling
- **What happened**: `summarise(across(c(a,b,c), ~ .x[which(!is.na(.x))[...]], .names = "last_{.col}"))` crashed with `Error in names(dots)[[i]] : subscript out of bounds` when a plate/well had an all-NA column (empty index vector). Repro was minimal and isolated to the `across`-with-`.names` path.
- **Why it matters**: The "take last non-NA per group" idiom is common for endpoints; the failure is opaque and can't be debugged from the message.
- **Resolution**: Replace with explicit per-column `summarise(x_late = last_non_na(x), ...)` using a small helper `last_non_na <- function(x){y <- x[!is.na(x)]; if(length(y)==0) NA else y[length(y)]}`.
- **mitigation_type**: structural — avoid variadic `across` when the lambda can return length-0 inside `summarise`.
- **Tags**: R, dplyr, summarise, across, last, endpoint, gotcha

### [2026-08-15] L-5 — lmerTest `difflsmeans()` needs the fitting environment, not a list wrapper; realize such "diff vs ref" tables from fixef/vcov instead
- **Category**: R / mixed-models
- **What happened**: Stashed lmer fits in a list and passed the *list* to `difflsmeans()`, which then failed with `no applicable method` (class "list"); passing the model object itself failed with `object 'copper_mm' not found` because the fitting data.frame was a local that difflsmeans re-evaluates via model.frame.
- **Why it matters**: Post-hoc contrast helpers in lmerTest/lme4 are brittle about env/data retention; the failure wastes a full model rerun cycle.
- **Resolution**: Since `factor(copper_mm)` coefficients are already differences vs the reference level, compute contrasts directly as `fixef(m)` + `sqrt(diag(vcov(m)))` Wald 95% CI (no emmeans needed — the pixi env lacks emmeans/AICcmodavg/minpack.lm).
- **Tags**: R, lme4, lmerTest, difflsmeans, contrasts, mixed-model, gotcha

### [2026-08-15] L-6 — rmarkdown chunk CWD = directory of the .Rmd, not the render() call; pdflatex chokes on Unicode ≥/≤/→/⇒/−/…
- **Category**: R / reporting
- **What happened**: Merged growth-rates report failed reading `"analysis/growth_rates/results/tables/..."` because knitr executes chunks with working directory set to the document's folder (so paths must be `results/tables`). Separately, pdflatex compiled em-dashes/× fine but aborted on U+2265/≤/→/⇒/−/…; xelatex in the shared TeX Live 2022 failed on a stale l3names kernel instead.
- **Why it matters**: Two silent portability traps cost several render cycles. The repo's existing plate-position report already follows the document-relative path convention.
- **Resolution**: In Rmds under this repo use doc-relative `results/tables` and keep Unicode to `—` and `×` (what the plate-position report proves pdflatex accepts); else switch the engine. Detect remaining non-ASCII with a quick python one-liner printing U+ codepoints before rendering.
- **mitigation_type**: structural — a report template/lint for allowed Unicode + doc-relative paths.
- **Tags**: R, rmarkdown, knitr, working-directory, pdflatex, unicode, reporting

### [2026-08-15] L-7 — a metadata-CSV label fix does NOT propagate to DuckDB: `strain` used `ON CONFLICT DO NOTHING` and full re-import is blocked by FK from measurements to imager_run
- **Category**: data / DB import
- **What happened**: Corrected two species-name misspellings in `data/metadata/Copper.Strain_info.csv`, then re-ran `10_import_experiment.py`: (1) the global `strain` upsert was `DO NOTHING`, so existing rows kept the old labels; (2) a full re-import crashed with `Constraint error: key run_number 353 is still referenced by a foreign key` — measurements reference `imager_run.run_number`, so the imager_run upsert cannot mid-flight replace rows already referenced.
- **Why it matters**: Analyses join species from the DB `strain` table, not from the CSV, so the CSV fix silently had no effect on any downstream output until the DB was updated. And the only "official" propagation path (full re-import) is structurally impossible in place once measurements load.
- **Resolution**: Switched the `strain` upsert to `DO UPDATE SET` the metadata columns (CSV authoritative) and added `--strain-only` to `10_import_experiment.py` to sync just the strain table without touching imager_run/well_placement; made the previously-required experiment/plate args optional with a fail-loud guard under `--strain-only`. Running `pixi run python scripts/db/10_import_experiment.py --strain-only --strain-info-csv data/metadata/Copper.Strain_info.csv` fixed strains 84 and 327.
- **mitigation_type**: structural — any metadata edit workflow should use `--strain-only`; the full import is only safe on a fresh DB.
- **Tags**: data, duckdb, import, metadata, strain, FK, copper

### [2026-08-15] L-8 — choosing a "timepoint" in this DB: bin hours_since_plate_start to the imaging pass, do NOT take per-strain max hour
- **Category**: data / phenotype extraction
- **What happened**: For a control-only late-timepoint strain phenotype table, taking each strain's exact maximum `hours_since_plate_start` (<=110 h) to define its "latest timepoint" yielded only 1-2 colonies/strain for 84 strains: within a single imaging pass, colonies of one strain image a few minutes apart, so an exact-max filter split one pass into fragments. Rounding to the nearest integer hour and taking the max *rounded* pass restored the full pass (median 4 colonies/strain).
- **Why it matters**: The imaging rig interleaves two ~3 h cadences; a fixed single hour covers only ~231/320 strains, while pass-of-latest-image-in-window covers 314/320 (286 with >=3 colonies). Any per-strain/plate timepoint pick should use the rounded-hour pass.
- **Resolution**: `with c0 as (select *, round(hours_since_plate_start,0) tp_h ...), latest as (select *, max(tp_h) over (partition by strain_id) tp_hmax from c0) ... where tp_h = tp_hmax` — i.e. round, then window-max the rounded hour.
- **mitigation_type**: structural — a reusable timepoint-binning helper / convention for time-series extractions on this pipeline data.
- **Tags**: data, duckdb, timepoint, imaging-pass, cadence, phenotype, gotcha, copper

### [2026-08-15] L-9 — user-site Python packages shadow the pixi conda env unless PYTHONNOUSERSITE=1
- **Category**: environment / pixi
- **What happened**: `numba` refused to import inside the pixi env with `TypeError: numpy.core.multiarray failed to import` because the conda env's numpy 2.4.x was being shadowed by a newer `~/.local/lib/python3.12/site-packages` numpy (2.5.2) on `sys.path`.
- **Why it matters**: User-site packages are a permanent source of version skew on this host; even after pinning `numpy = ">=2.4,<2.5"` in `pixi.toml`, the shadowing can still occur depending on `PYTHONUSERBASE`/`.pth` behavior.
- **Resolution**: Add to `pixi.toml` `[activation.env] PYTHONNOUSERSITE = "1"` so the env never reads the user-site directory; verify with `python -c "import numpy, sys; print(sys.path)"`.
- **mitigation_type**: structural — always set PYTHONNOUSERSITE in repo pixi.toml activation env.
- **Tags**: pixi, environment, numpy, numba, user-site, activation-env, gotcha

### [2026-08-15] L-10 — statsmodels `anova_lm` cannot compare multiple fitted models; compute the F-test by hand
- **Category**: statistics / pandas / statsmodels
- **What happened**: Passing several model refits to `anova_lm(m1, m2, ...)` for type-II tests of a strain×Cu interaction failed (multi-model comparisons not supported in this api).
- **Why it matters**: Hypothesis tests for GxE / added variables are a routine need; blocking on the statsmodels API was wasting time.
- **Resolution**: Fit the reduced and full model separately and compute `F = ((RSS_reduced - RSS_full) / (df_reduced - df_full)) / (RSS_full / df_full)`, p from `scipy.stats.f.sf`.
- **mitigation_type**: procedural — keep a small helper for hand-rolled nested-model F-tests.
- **Tags**: statsmodels, anova, F-test, nested-model, GxE, stats, gotcha

### [2026-08-15] L-11 — basal chroma noise (~1.5-2) in late Cu=0 colonies forces a minimum onset pigment threshold
- **Category**: analysis / phenotype
- **What happened**: With onset defined as first time chroma crosses a low threshold (2-3), nearly every colony was "onset at t=0" (`t_darkening_h=0` everywhere), and median onset collapsed to 0-12 h; threshold sweep showed median onset jumps 0/12/30/48/84 h as threshold rises 2→3→5→7→10.
- **Why it matters**: Baseline colony color has intrinsic chroma spread (~1.5-2) unrelated to pigment onset; any "time-to-pigment" derived from low thresholds reports imaging noise, not biology.
- **Resolution**: Use chroma threshold >=7 (`ColourLab`-scale chroma) for onset definitions; report the sweep as a robustness figure. Only ~51% of colonies reach chroma >= 7 under the IM scope window, so onset-definition must be a documented parameter.
- **mitigation_type**: analysis, procedural — threshold-sweep before committing to an onset definition.
- **Tags**: onset, chroma, threshold, pigment, timecourse, imaging, noise, gotcha

### [2026-08-15] L-12 — median one-time: naive residual-ON-residual mediation blows up numerically
- **Category**: statistics / modeling
- **What happened**: A rough "regress predictor on mediator residuals, then fit outcome on those residuals" mediation attempt produced wildly unstable / non-interpretable coefficients (collinearity of constructed regressors).
- **Why it matters**: Stats courses illustrate mediation with regressions of residuals, but those models are numerically pathological with correlated phenotype measures; the resulting "effects" are not meaningful.
- **Resolution**: Use a proper structural decomposition: outcome ~ species + area + Cu (total) versus outcome ~ species + Cu (direct effect with the mediator in the model), i.e. compare coefficients of Cu between models with/without the mediator; bootstrap the a*b indirect component.
- **mitigation_type**: analysis, procedural — use product-of-coefficients / model-difference decomposition, never residual-residual hacks.
- **Tags**: mediation, collinearity, regression, stats, bootstrap, gotcha

### [2026-08-15] L-13 — sd(log₁₀ colony area) is a size/floor-coupled metric; pool across strains and Cu and you will fabricate a "dispersion grows with stress" trend
- **Category**: analysis / statistics
- **What happened**: Idea 01's Gibrat check found sd(log₁₀ Asat) across replicate wells rising with Cu (0.30→0.43; pooled Spearman 0.23–0.28, p<1e-27). A within-strain follow-up showed raw paired slopes are indeed positive (+0.0033 per mM, Wilcoxon p≈5e-16) but that sd(log₁₀ Asat) is anti-correlated with colony size (Spearman −0.63, p<1e-200): Cu shrinks colonies, and small/floor-adjacent colonies carry inflated log-dispersion. Once colony size was included as a covariate, the within-strain Cu slope collapsed to ~0 (+0.0002, p=0.46; 53.6% positive).
- **Why it matters**: Any "heterogeneity / disorder / dispersion vs treatment" claim built on log-size (or log-intensity, log-chroma) dispersion across replicates will be confounded by the confound of treatment moving the mean near the detection floor. This subsumes the idea 01b diversity:higher-Cu => more log-variance => uncertainty story — it was an artifact.
- **Resolution**: For dispersion-type statistics, always (1) compute within-unit slopes (paired design), (2) control for the metric's own level (median size) as a covariate, and (3) report the pooled trend only alongside the size-controlled one. Here this revealed genuine lineage exceptions (R. taiwanensis widens, R. paludigena narrows) rather than a universal effect.
- **mitigation_type**: analysis, procedural — dispersion metrics need within-unit + size-controlled design, never naive pooled trends.
- **Tags**: dispersion, gibrat, noise, size-artifact, log-variance, variance, stats, floor-effect, gotcha
