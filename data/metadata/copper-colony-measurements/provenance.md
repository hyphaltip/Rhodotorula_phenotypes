<!-- Provenance for the Copper colony-measurement dataset -->
# Provenance: copper-colony-measurements

## Source

**Type**: instrument-derived (automated time-course colony imaging + segmentation)

**Origin**:
- Instrument: Automated plate imager (imager runs `d000353`-`d000357`, temperature token `300`). Exact make/model not recorded.
- Annotation sources (design/plate-layout lookup tables, live at `data/metadata/`):
  - `Copper.Strain_info.csv` — strain/well-placement layout per run (strain names, species, origin, environment, media, incubation temp).
  - `Copper.Plate_info.csv` — plate -> batch (Copper concentration mM) + replicate mapping.
- Derived from: Raw plate images segmented per well; per-colony feature export written to `data/preprocessed/*.parquet` (one file per image). Feature families are MorphoLibJ/CellProfiler-style (`Shape_`, `Intensity_`, `TextureGray_` Haralick, `Bbox_`, `Colorxy_`/`ColorLab_`/`ColorHSV_`).

**Citation / accession**: N/A

## Acquisition details

**Date acquired**: Preprocessed files span 2026-02-20 .. 2026-02-26 (imaged_at timestamps); placed in this repo around 2026-08-14.

**Obtained by**: Jason Stajich lab / instrumentation pipeline (not recorded in detail).

**Method**: Automated plate imaging on a fixed cadence (~6h interval) per run; segmentation + feature extraction per colony; CSV/Parquet export. Ingestion/load into DuckDB is scripted (`scripts/db/10_import_experiment.py`, `20_import_measurements.py`).

**Checksum**: Not recorded for the raw images; aggregate Parquet size = 428,465,532 bytes across 2,398 files. Source images are gitignored (see `data/preprocessed/`).

## Access restrictions

**Restriction level**: none (repo is internal; no data-use agreement in place)

**Details**: None.

## Known issues

- **Variable object counts per image** — not every well is fully detected at every timepoint; object counts vary by image. Treat per-image well occupancy as informative, not complete.
- **`Metadata_Dataset` / `Metadata_ImageName` dropped on import** — per-run temp identifiers / already encoded in the image filename; see `scripts/db/20_import_measurements.py` and `DATABASE_DESIGN.md`.
- **Fully-blank strain rows dropped** — `10_import_experiment.py` drops padding rows where `Strain ID` is null.
- **Plate-number -> run/configuration arithmetic** is an implicit Copper convention (28 plates/run, 4 configs) applied once at import time — see the fragile-convention note in `scripts/db/10_import_experiment.py`.

## Contact

**Primary contact**: Jason Stajich (jstajich)

**Backup contact**: None.

## Version history

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-08-15 | Initial metadata ingestion (schema, summary, provenance) |
| <!-- 1.1 --> | <!-- YYYY-MM-DD --> | <!-- e.g., "Added samples from cohort B." --> |