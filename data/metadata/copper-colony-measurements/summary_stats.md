# Summary Statistics: copper-colony-measurements

Computed by scripts/db/05_generate_metadata.py from DuckDB colony_measurement (all imported rows).
Generated: 2026-08-15

## Overview

| Property | Value |
|----------|-------|
| Image files (Parquet) | 2,398 |
| Measurement rows | 211,800 |
| Database columns | 181 (180 parquet - 2 dropped + 3 added) |
| Distinct well positions | 96 |
| Imager runs | 5 (d000353-d000357) |
| Date range | 2026-02-20 to 2026-02-26 |
| Aggregate source size | 428,465,532 bytes |
| Format | Parquet -> DuckDB `colony_measurement` |

## Per-run coverage

| Run | Images | Measurement rows |
|-----|--------|------------------|
| 353 | 585 | 57,385 |
| 354 | 558 | 44,904 |
| 355 | 558 | 45,794 |
| 356 | 556 | 37,188 |
| 357 | 141 | 26,529 |

## Column family layout

| Family | Meaning |
|--------|---------|
| `Metadata_` | 5 col(s) — Image/run metadata stamped onto every object by the segmentation export. |
| `ObjectLabel` | 1 col(s) — Sequential label of the detected object (colony) within its image. |
| `Bbox_` | 10 col(s) — Bounding-box geometry of the object in pixel coordinates. |
| `Grid_` | 4 col(s) — Row/column grid position the object was laid out in. |
| `Shape_` | 17 col(s) — Morphometric properties of the object (MorphoLibJ). |
| `Intensity_` | 12 col(s) — Pixel-intensity statistics of the object in the (gray) channel. |
| `TextureGray_` | 65 col(s) — Haralick co-occurrence texture metrics computed from the gray-level co-occurrence matrix at four orientations and their average (MorphoLibJ). |
| `Colorxy_` | 16 col(s) — Median-filtered chromaticity (x,y) statistics of the object. |
| `ColorLab_` | 26 col(s) — Perceptual CIE L*a*b* color statistics of the object. |
| `ColorHSV_` | 24 col(s) — HSV color statistics of the object (Hue, Saturation, Brightness). |

## Missing data summary

Probe of representative numeric columns (computed over all imported rows):

| Column | Missing count | Missing % | Pattern / notes |
|--------|---------------|-----------|------------------|
| None | 0 | 0.00% | No numeric probe columns with nulls |

## Most frequent categorical values

- **Metadata_FileSuffix**: .tif (211800)
- **Metadata_ImageType**: GridImage (211800)
- **Grid_RowNum**: 7 (32128), 1 (27257), 5 (26449)
- **Grid_ColNum**: 6 (19308), 5 (19056), 7 (18507)

## Quality flags

- **Object count can vary by image/timepoint** — wells not always fully detected;
  per-image object counts differ (see variable file sizes in preprocessed/).
- Metadata_Dataset / Metadata_ImageName source columns are dropped on import
  (per-dataset temp identifiers / already encoded in image_name). See DATABASE_DESIGN.md.

## Notes

- Image filenames encode run_number, temperature_token, plate_number, imaged_at  (scripts/db/lib/imagename.py).
- well_position is derived at import from Grid row/col (scripts/db/20_import_measurements.py).

