# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

### [2026-08-15] D-4 — Make strain metadata authoritative over CSV + add `--strain-only` sync to the DB import
- **Context**: Two species-name misspellings were corrected in `data/metadata/Copper.Strain_info.csv` (`Rhodotorula paludigenum` -> `paludigena`, `Rhodotorula evergladiensis` -> `evergladensis` for strains 84 and 327). The growth-rate analysis joins species from DuckDB `strain`, which had been populated with `ON CONFLICT (strain_id) DO NOTHING`, so the CSV fix did not propagate. A full re-import is impossible in place: measurements reference `imager_run.run_number` via FK, so the imager_run upsert on a populated DB raises a constraint violation.
- **Decision**: (1) Changed the `strain` upsert to `DO UPDATE SET` for the metadata columns (strain_code, strain_name, species, origin, environment), making the CSV authoritative for strain metadata on re-import; (2) added a `--strain-only` mode to `scripts/db/10_import_experiment.py` that syncs only the global `strain` table from the CSV and exits, avoiding imager_run/well_placement/condition_plate rewrites that FK constraints preclude after measurements load; made experiment-name/factor-name/plate-info-csv optional with a fail-loud guard when not in --strain-only.
- **Alternatives considered**: Direct one-off `UPDATE` SQL against the two rows (simplest but not reproducible); full DB rebuild from scratch (correct but requires re-importing all 211,800 measurements).
- **Rationale**: The repository treats the metadata CSVs as source of truth and the DuckDB as a derived store, so sync-on-conflict is the natural semantic; `--strain-only` gives a safe, repeatable path for exactly this class of metadata-correction fix without a full reload.
- **Consequences**: `strain` rows now refresh from the CSV on any re-import or `--strain-only` run; the analysis pipeline (00-05 + render) was regenerated cleanly with the corrected labels. Minor: curated strain metadata no longer survives a conflicting CSV value (documented trade-off; CSV is authoritative).
- **Tags**: data, metadata, import, duckdb, strain, copper

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
