# Data Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

### copper-colony-measurements
```yaml
name: copper-colony-measurements
type: imaging
source: automated time-course colony imaging (imager runs d000353-d000357) + segmentation export
annotation_sources: data/metadata/Copper.Strain_info.csv, data/metadata/Copper.Plate_info.csv
date_acquired: 2026-08-14
format: Parquet (2,398 files) -> DuckDB colony_measurement (181 cols)
rows: 211800
columns: 181
size: 428 MB (preprocessed parquets)
raw_path: data/preprocessed/
metadata_path: data/metadata/copper-colony-measurements/
db_table: colony_measurement
status: processed
known_issues:
  - Variable object count per image (wells not always fully detected)
  - Metadata_Dataset / Metadata_ImageName source columns dropped on import
  - Plate->run/configuration map is an implicit Copper convention applied at import
access_restrictions: none
tags: [copper, colony, rhodotorula, morphology, time-course, segmentation]
```

Per-colony measurements for the Copper exposure phenotyping experiment (runs 353-357,
temperature token 300). One row per segmented colony object per image; 211,800 rows across
2,398 images, keyed by `image_name` + `object_label`. Feature families: Shape, Intensity,
Haralick Texture (gray), Bbox, and color (xy / CIELAB / HSV). Strain/plate annotation lives
in `Copper.Strain_info.csv` and `Copper.Plate_info.csv` (see provenance.md). Fully loaded into
the DuckDB `colony_measurement` table by `scripts/db/`.
