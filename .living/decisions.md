# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

### [2026-08-15] D-3 — Stratify plate-effect variance partition by Copper concentration
- **Context**: `analysis/explore_plate_position/01_variance_partition.R` estimated a pooled mixed model (`trait ~ copper_mm + grid_row_c + grid_col_c + is_edge + (1|strain_code) + (1|run_number) + (1|plate_id)`) reporting 23% (L*) and 33% (b*) plate-explained variance. Inspection showed `plate_id` is **nested** inside `copper_mm`: every plate was imaged at exactly one Cu level (0 of 112 plates have >1 concentration), so a random slope of Cu within plate is not estimable and the pooled plate term can absorb any non-linear Cu response.
- **Decision**: Replaced the pooled model with a per-Cu-concentration stratified analysis: fit the variance partition independently within each of the 7 Cu levels (0-30 mM), dropping `copper_mm` (invariant within a stratum), applying the ≥2-plates-per-strain filter within each stratum, plus the categorical row/col robustness check per stratum.
- **Alternatives considered**: Keep pooled + add per-stratum as detail (rejected: user chose stratified-only); nest `(1|copper_mm:plate_id)` in one model (redundant given Cu is a fixed effect and plate is already nested in Cu).
- **Rationale**: Stratification removes the plate↔Cu confound directly and yields per-condition plate variance estimates, which are the quantities of interest given the stress gradient (0 mM control vs 30 mM).
- **Consequences**: Stratified plate variance is much smaller (L* ~0.1-4.7%, b* 0-3.2%) than the pooled 23%/33% — the pooled numbers were substantially an artifact of pooling across concentrations with a non-linear Cu response. Outputs renamed to `*_by_copper.csv`; pooled files (`variance_components.csv`, `fixed_effects.csv`, lmer_models_linear/categorical.rds) deleted. Adjacency test (02) unchanged (already controls `copper_mm`).
- **Tags**: analysis, copper, plate-position, variance-partition, mixed-model

### [2026-08-15] D-2 — Recorded data metadata for the Copper phenotype dataset (mycelium ingest)
- **Context**: The mycelium init installed the metadata framework but the Copper dataset (2,398 parquet files, 211,800 measurements in DuckDB `colony_measurement`) had no per-dataset `schema.yaml`/`provenance.md`/`summary_stats.md` and the `DATA_MANIFEST.md` was empty.
- **Decision**: Registered the dataset as `copper-colony-measurements` with metadata generated from the populated database (not a file sample) so counts/missingness are exact. Added `scripts/db/05_generate_metadata.py` as a reproducible generator. The two annotation CSVs remain at `data/metadata/` and are referenced as annotation sources in provenance/manifest.
- **Alternatives considered**: Hand-written metadata; sampling a few parquets for stats; registering the two annotation CSVs as separate datasets.
- **Rationale**: Computing from DuckDB gives exact, reproducible stats and avoids drift between the DB and metadata. One primary dataset entry keeps the manifest readable; annotation CSVs are lookup tables rather than independent measurement datasets.
- **Consequences**: Metadata now reflects DB reality (181 cols = 180 parquet - 2 dropped + 3 added; per-run coverage; no nulls in probed numeric columns). Maint. = re-run `05_generate_metadata.py` on new data.
- **Tags**: data, metadata, ingest, copper

### [2026-08-15] D-1 — Completed mycelium init: install all convention packs; skip broken skillpack
- **Context**: Finishing a partially-completed mycelium `init` for this repo. Core packs are mandatory; domain packs were user-selected.
- **Decision**: Installed core (`robust-analysis`, `report-generator`, `idea-generator`) plus domain (`bioinformatics`, `image-analysis`, `skill-bridge`) convention packs. Cloned `scientific-agent-skills` and `bioSkills` into `skillpacks/`; did **not** clone `Autonomous-Science` because `github.com/arjunrajlaboratory/Autonomous-Science` no longer exists publicly.
- **Alternatives considered**: Install only core packs; skip domain packs until needed; wait for `Autonomous-Science` to be restored.
- **Rationale**: User explicitly requested all three domain packs. The dead repo cannot be cloned, and blocking on it stalls init; the `skill-bridge` pack still functions for the two successfully-cloned skill repos.
- **Consequences**: `skill-bridge` persona-routing (personas/library) is currently unavailable (see L-2). Official CLAUDE.md generated; `.living/INDEX.md` built; `validate_structure.py` passes.
- **Tags**: setup, mycelium-init, conventions, skillpacks
