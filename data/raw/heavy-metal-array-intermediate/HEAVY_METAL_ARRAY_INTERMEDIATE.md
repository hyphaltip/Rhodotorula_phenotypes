# heavy-metal-array-intermediate

Per-colony intermediate measurements for five arrayed heavy-metal screens (Cr, Cu, Fe, Pb, Zn).
The source author calls this the "intermediate phase": compiled from primary tool output, before summary statistics (see `source_about.md`).

## Layout

- `<Metal>_<metal>measurementmeta.csv`: symlinks to the shared-lab CSVs (not copied; 1.2 GB total).
- `SHA256SUMS.txt`: checksums of the CSVs at ingest (2026-10-07). Verify with `sha256sum -c SHA256SUMS.txt`.
- `source_about.md`: copy of the source `about.md`.
- Parquet copies (one per metal): `data/preprocessed/heavy_metal_array/<Metal>.parquet`. These are gitignored. Rebuild with `scripts/ingest/heavy_metal_csv_to_parquet.py`.

## Rebuild

```bash
srun -p short -c 4 --mem=16G -t 30 pixi run python scripts/ingest/heavy_metal_csv_to_parquet.py \
  --src-dir data/raw/heavy-metal-array-intermediate --out-dir data/preprocessed/heavy_metal_array
pixi run python scripts/ingest/heavy_metal_make_metadata.py
```

## Single-table model in DuckDB

All five files can be read as one table. Use `union_by_name=true`:

```sql
SELECT * FROM read_parquet('data/preprocessed/heavy_metal_array/*.parquet', union_by_name=true);
```

- 517,371 rows, 179 columns. `Metal` is already a column in every file.
- 147 columns are shared by all metals. 32 columns exist only in Lead and Zinc and are NULL for Cr/Cu/Fe.
- Run numbers do not overlap between metals.

## Caveats

- `Concentration` has no unit in the source. Cr spans 0-1.2. The other metals span 0-30. Do not pool Cr with the others until the unit is confirmed.
- `strain_id` is text. It holds numeric IDs, `Control-N` labels, and NULL (10,921 rows: Cr 5,923; Cu 3,489; Pb 1,509).
- Pb and Zn have 32 extra columns, so they come from a later version of the image-analysis pipeline. Cross-metal comparison of the extra columns is not possible.
- Copper run numbers (d000353-d000357) match the existing `copper-colony-measurements`. Values were not compared at ingest.
- `Concentration` was cast to DOUBLE in all Parquet files. No other values were changed.
- The source files may be updated by their owner. Check `SHA256SUMS.txt` before relying on the symlinks.
