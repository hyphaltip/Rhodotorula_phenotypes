# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

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
