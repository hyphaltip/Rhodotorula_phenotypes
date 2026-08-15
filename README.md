# Rhodotorula copper-stress colony phenotyping

Reference documentation for the colony-based phenotyping dataset: per-image colony
measurements from an arrayed *Rhodotorula* screen under copper stress, plus the
strain and experiment metadata needed to interpret them.

This document is the entry point for anyone writing an analysis script in `scripts/`.
It describes the data format, the directory layout, and — most importantly — exactly
how the files join together.

---

## Directory layout

```
data/
  metadata/
    Copper.Plate_info.csv        # 120 rows: Batch Number, Concentration (mM), Replicate
    Copper.Strain_info.csv       # 2400 data rows, 14 columns: strain identity per Run/Config/Well
  preprocessed/
    d000353_300_028_2026-02-21_10-31-05.parquet    # one file per imaged plate/timepoint
    ...                                            # 2398 parquet files total, Parquet only

ignore/                          # background only -- not read by any analysis script
  About.MD                       # human-written notes on file naming and the join process
  d000355_300_079_2026-02-24_17-14-27.csv   # leftover one-off test export of one plate

scripts/                         # new analysis scripts go here; read from data/ and its subfolders
lib/                             # currently empty
```

All analysis work reads from `data/`. `ignore/` holds background material on how this
dataset used to be organized (the original `About.MD` naming notes, and a one-off test
CSV export of a single plate) — useful for context, not part of the pipeline going
forward.

---

## Per-image measurement files

`data/preprocessed/*.parquet`. One file per imaged plate at one timepoint; one row per
detected colony object within that image. Roughly 50–90 rows per plate — not all 96 grid
positions have a detected colony, so occupancy varies.

**Parquet is the only supported input format.** All measurement data should be read from
the `*.parquet` files in `data/preprocessed/` (2398 files).

Columns, in order, by group:

| Group | Columns | Meaning |
| --- | --- | --- |
| Image identifiers | `Metadata_Dataset`, `Metadata_FileSuffix`, `Metadata_BitDepth`, `Metadata_ImageType`, `Metadata_ImageName` | Image-level identifiers. `Metadata_ImageName` is the key that gets decoded into run/temperature/plate/datetime. |
| Object ID | `ObjectLabel` | Per-colony object ID within the image. |
| Bounding box | `Bbox_*` | Bounding box center/min/max coordinates, in pixel space plus intensity- and distance-weighted variants. |
| Grid position | `Grid_RowNum` (0–7), `Grid_ColNum` (0–11) | 0-indexed position of the colony in the 8x12 (96-position) plate grid. |
| Grid index | `Grid_RowMajorIdx` = `Grid_RowNum * 12 + Grid_ColNum` | 0-indexed row-major numbering. |
| Grid index | `Grid_ColMajorIdx` = `Grid_ColNum * 8 + Grid_RowNum` | 0-indexed column-major numbering. **This is the field used for the well-position join**, not `Grid_RowMajorIdx`. |
| Shape | `Shape_*` | Colony morphology: area, perimeter, circularity, radii, Feret diameters, eccentricity, solidity, axis lengths, etc. |
| Intensity | `Intensity_*` | Grayscale intensity stats: integrated, min/max/mean/median, stddev, coefficient of variation, quartiles, density. |
| Texture | `TextureGray_*` | Haralick texture features (angular second moment, contrast, correlation, variance, entropy, etc.) at 4 angles (`deg000`/`deg045`/`deg090`/`deg135`) plus an `-avg-` variant, all at `scale05`. |
| Color | `Colorxy_*`, `ColorLab_*`, `ColorHSV_*` | Per-channel summary stats (min/Q1/mean/median/Q3/max/stddev/coeffvar) in xy, CIE L\*a\*b\*, and HSV color spaces, plus `ColorLab_ChromaEstimatedMean` / `ColorLab_ChromaEstimatedMedian`. |

### Filename convention

From `ignore/About.MD` (authoritative background notes). Example:
`d000320_300_001_2026-02-01_15-44-58`

| Field | Example | Meaning |
| --- | --- | --- |
| Imager run number | `d000320` | Which imager run produced the image. Parsed as `Batch_Number` in code; called "Run Number" in `About.MD` and in the strain metadata. |
| Temperature | `300` | Temperature in Celsius x10, i.e. `300` = 30.0 °C. Constant `300` across every file in this dataset. |
| Plate/batch position | `001` | Plate position within the imager. Valid range 1–160; this dataset uses 1–120. Parsed as `Plate_Number`. |
| Datetime | `2026-02-01_15-44-58` | `year-month-day_hour-min-sec` the image was taken. |

Standard parse (Polars):

```python
df = df.with_columns(
    pl.col("Metadata_ImageName").str.split("_").list.get(0).alias("Batch_Number"),
    pl.col("Metadata_ImageName").str.split("_").list.get(1).alias("Temperature"),
    pl.col("Metadata_ImageName").str.split("_").list.get(2).alias("Plate_Number"),
    pl.col("Metadata_ImageName").str.split("_").list.slice(3, 4).list.join("_").alias("Datetime"),
)
```

---

## Metadata files

### `data/metadata/Copper.Plate_info.csv` — experiment conditions

120 rows. Columns: `Batch Number`, `Concentration (mM)`, `Replicate`.

- `Batch Number` (1–120) is the *global* plate/batch position. It is exactly
  `Plate_Number` from the filename — joined directly, no transformation.
- `Concentration (mM)` cycles 0, 5, 10, 15, 20, 25, 30 mM (7 steps) within contiguous
  blocks of 28 `Batch Number`s.
- `Replicate` values are `1`, `2`, `3`, `4`, `Control` — 28 rows each for 1–4, and 8 rows
  for `Control` (batches 113–120).
- Verified empirically: `Batch Number` blocks 1–28, 29–56, 57–84, 85–112, 113–120
  correspond 1:1 to imager run numbers 353, 354, 355, 356, 357 respectively. This was
  confirmed by cross-tabulating every `Plate_Number` that appears against its file's
  `Batch_Number` across all 2399 preprocessed files, with no exceptions. `Replicate` is
  therefore a proxy for imager run, not a technical-replicate count (see quirks).

### `data/metadata/Copper.Strain_info.csv` — strain identity

2400 data rows (2401 lines including the header), 14 columns:
`Index`, `Strain ID`, `Strain Name`, `Strain`, `Species`, `Origin`, `Environment`,
`Run Number`, `Library Plate`, `Configuration`, `Well position`,
`Incubation Temp (°C)`, `Time stamped`, `Media`.

- `Run Number` values present: 353, 354, 355, 356. **357 never appears** — it is the
  Control-only imager run and has no strain rows.
- Each `Run Number` maps 1:1 to a `Library Plate`: 353→1, 354→3, 355→2, 356→4.
- Each `Run Number` has 336 strain rows split across `Configuration` values 1–4
  (84 rows each), i.e. up to 84 of the 96 possible well positions are populated. The
  skipped positions follow a consistent pattern (for example, Run 355 / Configuration 1
  is missing positions 1, 11, 21, 31, 34, 44, 54, 64, 65, 75, 85, 95), consistent with
  intentional edge/corner exclusion in the rearray design rather than a data gap.

---

## How the files join

Two sequential left joins onto the per-image measurement table.

```mermaid
flowchart TD
    P["data/preprocessed/*.parquet<br/>one row per colony object<br/>Metadata_ImageName, Grid_ColMajorIdx, Shape_*, Intensity_*, ..."]

    D["Decode Metadata_ImageName<br/>split on '_'<br/>→ Batch_Number, Temperature,<br/>Plate_Number, Datetime"]

    H["Derive Hours<br/>Datetime − min(Datetime)<br/>over Plate_Number"]

    PI["data/metadata/Copper.Plate_info.csv<br/>Batch Number, Concentration (mM), Replicate<br/>(filter out Replicate == 'Control' first)"]

    J1["JOIN 1 — experiment condition (left)<br/>Plate_Number (int) == Batch Number"]

    DER["Derive strain join keys<br/>run_number = 353 + (Plate_Number−1)//28<br/>configuration = (Plate_Number−1)//7 %% 4 + 1<br/>well_position = Grid_ColMajorIdx + 1"]

    SI["data/metadata/Copper.Strain_info.csv<br/>Run Number, Configuration, Well position<br/>Strain ID, Strain Name, Species, Origin, ..."]

    J2["JOIN 2 — strain identity (left)<br/>(run_number, configuration, well_position)<br/>== (Run Number, Configuration, Well position)"]

    OUT["Combined strain + experiment + measurement table<br/>one row per colony per image"]

    P --> D --> H --> J1
    PI --> J1
    J1 --> DER --> J2
    SI --> J2
    J2 --> OUT
```

### Join 1 — experiment condition

`Plate_Number` (cast to integer, parsed from the filename) is matched against
`Plate_info.Batch Number`, after first dropping `Replicate == "Control"` rows from the
plate metadata. This attaches `Concentration (mM)` and `Replicate`.

### Join 2 — strain identity

Three keys are derived from `Plate_Number` and the per-colony `Grid_ColMajorIdx`:

```python
run_number    = 353 + (Plate_Number - 1) // 28     # 353 = min non-null Run Number in Strain_info
configuration = (Plate_Number - 1) // 7 % 4 + 1
well_position = Grid_ColMajorIdx + 1               # column-major, 1-indexed
```

These are matched against `Strain_info`'s `(Run Number, Configuration, Well position)`.

Because it is a left join on a unique key triple, it should never fan out or drop rows —
assert the row count is unchanged after this join as a sanity check. Rows with no strain
match (Control batches 113–120, and any unoccupied well position) simply carry null
strain columns and should be filtered out downstream.

### The `Hours` column

`Hours = Datetime - min(Datetime) grouped by Plate_Number`, as whole hours. This gives
elapsed imaging time per plate. Each plate is imaged repeatedly over time, so `Hours` is
the time axis for growth-curve and kinetics analysis.

---

## Quirks and inconsistencies

Each item below was verified directly against the files, not inferred.

1. **`Temperature` is not actually variable.** All 2399 files in `data/preprocessed/`
   have literally `300` in the temperature filename slot, and `Strain_info`'s own
   `Incubation Temp (°C)` column is always 30. It is a fixed value baked into the
   filename, not a per-experiment variable — despite being modeled as a "Temperature"
   field, which invites treating it as if it varies.

2. **The well-position join uses `Grid_ColMajorIdx`, not `Grid_RowMajorIdx`, despite the
   row-major field's more intuitive name.** Both fields exist in every measurement row:
   `Grid_RowMajorIdx = row*12 + col`, `Grid_ColMajorIdx = col*8 + row`. Anyone
   re-deriving this join from first principles (as happened during initial exploration
   of this repo) will naturally reach for `RowMajorIdx` and get a silently wrong join —
   it will not error, it will just mismatch strains to wells. Easy to make, hard to
   detect.

3. **The `Replicate` column name in `Plate_info.csv` is misleading.** Its values
   (1, 2, 3, 4, Control) are not independent technical replicates of the same condition
   — they are a 1:1 proxy for which imager run (353, 354, 355, 356, 357) the plate
   belongs to. Analyzing "replicate effects" without knowing this conflates imager-run
   batch effects with true biological or technical replication.

4. **The run-number/configuration derivation is not stored anywhere as data — it is
   arithmetic inferred from plate-number block boundaries.**
   `run_number = 353 + (Plate_Number-1)//28` and
   `configuration = (Plate_Number-1)//7 % 4 + 1` are documented in neither metadata CSV;
   they were reverse-engineered from the plate numbering (per `ignore/About.MD`), and
   are only correct because plate numbering happens to be laid out in perfectly regular
   contiguous blocks of 28. Any future batch that does not follow this exact block size
   and order (a run with a different number of configurations or concentration steps,
   for instance) will silently break the formula, with no validation in place to catch
   it.

5. **The Control-only imager run (`d000357`, `Plate_Number` 113–120) has no
   corresponding rows in `Strain_info.csv` at all** — `Run Number` 357 never appears
   there — and it is also explicitly filtered out of `Plate_info.csv` by
   `Replicate == "Control"` before the first join. These plates are dropped from the
   pipeline for two independent, overlapping reasons. Worth confirming that is
   intentional, and that the same filter logic could not silently exclude a *non-control*
   row in the future.

6. **Not every well position (1–96) is populated for a given Run/Configuration.**
   `Strain_info` has only 84 of 96 possible positions filled per configuration, in a
   consistent skipped pattern (for example, positions 1, 11, 21, 31, 34, 44, 54, 64, 65,
   75, 85, 95 are missing for Run 355 / Configuration 1). This looks like intentional
   plate-edge exclusion in the rearray design rather than a data error, but the rule is
   documented nowhere — it is only visible empirically.

7. **`Copper.Strain_info.csv` has 2400 data rows by line count, but only 1284 of them
   are real.** The other 1116 are fully-blank padding rows — every column empty except
   a populated `Index` value (e.g. `Index=1285` through `2400` with every other field
   blank). A naive `wc -l`-style row count, or any load that doesn't explicitly filter
   on `Strain ID IS NOT NULL`, silently carries 1116 all-null strain rows into any
   downstream join. Real per-run strain counts are uneven, not a flat 336 each as a
   single-run sample might suggest: Run 353→336, 354→332, 355→336, 356→280.

### Resolved since initial exploration (no longer concerns)

- Duplicate copies of both metadata CSVs (previously under a `PhenotypicMeasurements/`
  symlink) no longer exist in the working tree — `data/metadata/` is the single source.
- The stray one-off test CSV (`d000355_300_079_2026-02-24_17-14-27.csv`) has been moved
  out of `data/preprocessed/` into `ignore/` and is no longer a candidate input.
- `Copper.Strain_info.csv`'s `Well position` header no longer has a trailing space.
- `Strain Name` values no longer contain embedded newlines inside quoted CSV fields;
  the file now parses to its true row count with any tool, not just a real CSV parser.
