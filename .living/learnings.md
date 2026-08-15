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
