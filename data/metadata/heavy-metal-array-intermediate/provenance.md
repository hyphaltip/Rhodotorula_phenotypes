# Provenance: heavy-metal-array-intermediate

- **Source**: shared lab project `ArrayedHeavyMetalScreen`, directory `/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Phenotyping/Heavy_Metals/ArrayedHeavyMetalScreen/Data/Interm/`
- **Files**: `CrArrayRun/chromiummeasurementmeta.csv`, `CuArrayRun/coppermeasurementmeta.csv`, `FeArrayRun/ironmeasurementmeta.csv`, `PbArrayRun/leadmeasurementmeta.csv`, `ZnArrayRun/zincmeasurementmeta.csv`
- **Source owner**: file owner on the shared drive is `cona002`. Contact not otherwise recorded.
- **Source file dates**: Cr 2026-09-22, Cu 2026-09-22, Fe 2026-09-30, Pb 2026-09-22, Zn 2026-09-24. `about.md` dated 2026-10-05.
- **Date acquired**: 2026-10-07
- **Acquisition method**: symlink of the CSVs into `data/raw/heavy-metal-array-intermediate/`. SHA256 in `SHA256SUMS.txt`.
- **Transformation**: CSV to Parquet (zstd) with DuckDB. `Concentration` cast to DOUBLE. Row counts checked equal to the CSV for every file. Script: `scripts/ingest/heavy_metal_csv_to_parquet.py`.
- **Imaging runs**: Cr d000320-324, Cu d000353-357, Fe d000406-408, Pb d000399-403, Zn d000388 and d000390.
- **Access restrictions**: shared lab storage, not redistributable.
- **Known issues**: see `HEAVY_METAL_ARRAY_INTERMEDIATE.md` (Caveats).
