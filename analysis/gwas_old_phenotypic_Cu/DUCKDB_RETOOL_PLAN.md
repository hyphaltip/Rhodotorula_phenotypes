# GWAS retool to DuckDB: plan and status (2026-10-08)

Goal: take strain and phenotype inputs for `analysis/gwas/` from `db/rhodotorula_phenotypes.duckdb`
(`strain_info` view, `heavy_metal_measurement` table). The removed inputs were
`data/metadata/Copper.Strain_info.csv`, `v_phenotype`, and `analysis/gwas/data/db_extract.parquet`
(an extract of `v_phenotype`).

All numbers below were computed in SLURM jobs on 2026-10-08. Logs are in
`analysis/gwas/results/duckdb_retool/`. Nothing in `analysis/gwas/results/` that existed before was changed.

## 1. Where the removed inputs are used

| File | Line(s) | What it reads | Needs |
|---|---|---|---|
| `analysis/gwas/run.sh` | 12 (old) | `data/metadata/Copper.Strain_info.csv`, `--pheno-col Strain` | one column of strain IDs |
| `analysis/gwas/scripts/reconcile_strains.py` | 6, 55-58 | the CSV passed as `--pheno` | column named by `--pheno-col`; reads only that column |
| `analysis/gwas/scripts/common.py` | 18-19, 47-48 | `data/db_extract.parquet` (`read_extract`), `data/strain_metadata.tsv` (`read_meta`) | see below |
| `analysis/gwas/scripts/build_gwas_phenotypes.py` | 74, 115 | `CM.read_extract()` | `species, strain_code, Shape_Area, ColorLab_ChromaEstimatedMedian, ColorHSV_Saturation/BrightnessMedian, ColorLab_{L*,a*,b*}Median, TextureGray_*-avg-scale05 (13), tp_h, copper_mm, run_number, plate_number, well_position` |
| `analysis/gwas/scripts/check_mas_gates.py` | 205, 222-228 | `--extract-parquet` (db_extract), uses `species, strain_code, copper_mm, tp_h, run_number` + trait column | timepoint-robustness of Cu=0 traits at windows |
| `analysis/gwas/scripts/add_copper_v0151_trait.py`, `compare_copper_v0151_trait.py`, `summarize_cu_doseauc_v0151_gemma.py`, `rare_variant_carrier_permutation_test.py`, `run_population_vs_locus_cu_v0151.sh` | - | `data/raw/copper-heavy-metal-screen-v0.15.1/copper_auc_mean_by_strain.csv` (dataset removed, D-35) | cannot be rebuilt; the existing column `cu_doseauc_v0151` in `gwas_next_phenotypes_fam_order.csv` stays as it is |
| `analysis/gwas/GWAS.md` | 13, 22, 741, 1026 | prose only | no code change |
| `analysis/ideas/2026-08-15-color-phenotype-space/scripts/build_series.py` | 91-100 | `v_phenotype` -> `data/db_extract.tsv.gz`; `strain` table | asserts 211,800 rows; cannot run |
| `analysis/ideas/.../scripts/common.py` | 18 | `data/db_extract.parquet` | same as above |
| `analysis/control_late_timepoint_phenotype/scripts/build_phenotype_table.py` | 27, 92 | `v_phenotype`, `imager_run`, `condition_plate_factor`, `strain` | needs rewrite (see 5) |
| `scripts/db/00_init_schema.sql`, `10_import_experiment.py`, `20_*`, `05_generate_metadata.py`, `query_examples/*` | - | `colony_measurement`, `well_placement`, `Copper.*_info.csv`, `v_phenotype` | legacy; D-35 says not to run them |

Existing GWAS results (GEMMA outputs, `gwas_next_phenotypes*.csv`, `strain_match_table.reviewed.csv`) do not read the DB.
They stay valid for the old data. They are not regenerated here.

## 2. Strain ID mapping

`reconcile_strains.py` reads only the column given by `--pheno-col`. The old column was `Strain`
(e.g. `TFCN_17-291Y-1`). The same string is `strain_info.sample_name` in the new data.

| strain_info | exported column | Use |
|---|---|---|
| `sample_name` | `Strain` | matched to VCF sample IDs |
| `strain` | `Strain_name` | lab name, e.g. `17-291Y-1 BY126-B8`; not used for matching |
| `species` | `Species` | panel filter (`Rhodotorula mucilaginosa`) |
| `strain_id` | `strain_id` | key to `heavy_metal_measurement.strain_id` |

`analysis/gwas/scripts/duckdb_inputs.py` exports this table (non-control strains only).
Output: `analysis/gwas/results/strain_reconciliation_duckdb/strain_info_from_duckdb.csv` (321 rows).

## 3. Coverage (computed)

`bcftools query -l` on the VCF gives 422 sample IDs (`duckdb_retool/vcf_sample_ids.txt`).

| Quantity | Value |
|---|---|
| Old `Copper.Strain_info.csv`: rows / distinct `Strain` | 2,400 / 308 |
| New `strain_info`: non-control strains / distinct `sample_name` | 321 / 320 |
| Old and new share the same string | 296 |
| Only in old / only in new | 12 / 24 |
| New `sample_name` exactly equal to a VCF ID | 213 |
| Old `Strain` exactly equal to a VCF ID | 213 |
| GWAS panel (reviewed `accept` list = .fam strains) | 213 |
| Panel strains found by exact `sample_name` in the new table | 212 of 213 |
| Panel strain not found | `TFCN_223A-8` |

Details:

- The one missing panel strain: new `sample_name` is `TFCN_223D-8` (strain_id 65, lab name `223D-8 BY35-B9`, species 'Species Not Found').
  The VCF ID is `TFCN_223A-8`. `reconcile_strains.py` gives a fuzzy match (score 0.909). This is probably the same strain
  with a letter difference, but I could not verify that. A person must decide.
- Spelling differs for 12 old/new pairs (old `TFCN_152A_3`, new `TFCN_152A-3`; underscores became hyphens). None is in the panel: the exact-match
  sets differ only by `TFCN_223A-8` (old) and `DBVPG_10656` (new).
- Re-running `reconcile_strains.py` on the new table (output in `strain_reconciliation_duckdb/`): exact 213, fuzzy 74, unmatched 33
  (old run: exact 213, normalized 1, fuzzy 62, unmatched 32). Collisions: 54 (old run: 48). The script exits with status 1 when there
  are collisions. It did so before too, so `run.sh` stops at that step by design.
- 23 fuzzy/unmatched strains in the new table have no row in the old `strain_match_table.reviewed.csv`. They need human review.
  Exact matches: the set is the same except `TFCN_223A-8` (old) and `DBVPG_10656` (new). `DBVPG_10656` is `Rhodotorula mucilaginosa`
  in the new table. It is exactly in the VCF, but not in the 213 panel. Whether it was excluded earlier (ploidy, clone culling)
  or is new, I did not check.
- `strain_id` 165 and 269 share `sample_name` `TFCN_17-332Y-1` (species 'Species Not Found'). The export keeps both rows and sets
  `duplicate_sample_name=True`. `reconcile_strains.py` collapses them to one ID. Neither is in the panel.
- Species: of the 212 panel strains found, all 212 are `Rhodotorula mucilaginosa` in the DB now. Species text changed vs the old CSV for 18 of
  296 shared strains (spelling such as `sp. clade I` vs `sp_clade_I`, plus some re-assignments, e.g. `TFCN_1A-1-3` mucilaginosa -> pacifica).
  The parent task is changing species values; re-check before use.
- Tested metals for the 212 panel strains: all 5 metals 121; Cr,Cu,Fe,Pb 53; Cr,Cu,Pb 32; Cr,Cu,Pb,Zn 5; Cr,Cu 1.

## 4. What was changed

1. `analysis/gwas/scripts/duckdb_inputs.py` (new): opens the DB read-only. If another process holds the write lock, it copies the file to
   `$SCRATCH` and opens the copy. (The DB was write-locked during this work.) Exports the strain table.
2. `analysis/gwas/run.sh`: exports the strain table, then reconciles with `--pheno .../strain_info_from_duckdb.csv --pheno-col Strain`.
   Output goes to `strain_reconciliation_duckdb/`, so the reviewed table in `strain_reconciliation/` is not overwritten.
   `BASH_SOURCE[0]` was removed; the script uses `$PWD` or `$GWAS_REPO_ROOT`. Step 6 text points to the new builder.
3. `analysis/gwas/scripts/check_duckdb_strain_coverage.py` (new): produced the coverage numbers above.
4. `analysis/gwas/scripts/build_gwas_phenotypes_duckdb.py` (new): re-derives the Copper traits from `heavy_metal_measurement`
   into `results/duckdb_retool/`. It does not touch `gwas_next_phenotypes*.csv`.
5. `analysis/gwas/scripts/check_old_new_object_match.py` (new): attempt to match old and new objects (see 6).

Not changed: `common.py`, `build_gwas_phenotypes.py`, `check_mas_gates.py`, `analysis/ideas/*`, `analysis/control_late_timepoint_phenotype/*`.
They still read removed inputs and will fail if run. See 7.

## 5. Trait re-derivation from `heavy_metal_measurement`

Column mapping (new data has no `hours_since_plate_start`, `copper_mm`, `plate_number`, `well_position`, `strain_code`):

| Old | New | Status |
|---|---|---|
| `strain_code` | `strain_info.sample_name` via `strain_id` | exact string match for the panel (212/213) |
| `copper_mm` | `Concentration` where `Metal='Copper'` | same dose grid 0,5,...,30 |
| `hours_since_plate_start` | `capture_datetime - min(capture_datetime) per (run_number, plate_position)` | the old rule was `min(imaged_at) OVER (PARTITION BY run_number, plate_number)` (read from the old DB backup `v_image`). Same images: old and new image names agree for 2,397 of the 2,398 old Cu images. |
| `plate_number` | `plate_position` | |
| `well_position` | `Grid_RowMajorIdx` | well identity inside a plate; the old-to-new well index mapping was not verified |
| `Shape_Area` | `Shape_Area` | same name |
| `ColorLab_{L*,a*,b*}Median` | `ColorLab_{L*,a*,b*}GeoMedian` | different statistic |
| `ColorLab_ChromaEstimatedMedian` | `sqrt(a*^2+b*^2)` of GeoMedian | approximation |
| `ColorHSV_SaturationMedian`, `...BrightnessMedian` | `ColorHSV_SaturationRobustMean`, `ColorHSV_ValueRobustMean` | approximation |
| `TextureGray_<M>-avg-scale05` | `Texture_<M>-avg-scale05` | same 13 metrics, name changed; same definition not verified |

Not possible without the removed dataset: `cu_doseauc_v0151` (source `copper-heavy-metal-screen-v0.15.1` was removed).

Implementation (`build_gwas_phenotypes_duckdb.py`), same algorithm as `build_gwas_phenotypes.py`: Cu=0, 85-110 h window, median per plate then mean over plates;
AUC by trapezoid over rounded hour per well, averaged over wells; slope, IC50 as before. Added step: 129 (well, rounded hour) duplicates are collapsed to their median.
Result: 215 strains in the new table; 211 of the 213 panel strains aligned. Missing: `TFCN_223A-8` (name issue above) and `DBVPG_10619`
(strain_id 317 has only 47 colony observations in total, versus about 1,300-2,200 for typical strains, so it has too few Cu rows).

### Old vs new values (211 panel strains, `trait_old_vs_new_comparison.csv`)

The comparison shows poor agreement, including traits whose definition is unchanged:

| Trait | n | Pearson | Spearman |
|---|---|---|---|
| clone_mean_area (log10 area, Cu=0, 85-110 h) | 210 | 0.486 | 0.477 |
| AUC_0 | 210 | 0.372 | 0.451 |
| AUC_10 / AUC_20 / AUC_30 | 210/211/211 | 0.175 / 0.363 / 0.175 | 0.272 / 0.342 / 0.340 |
| AUC_ratio_10 | 210 | 0.099 | 0.204 |
| cu_dose_slope | 209 | 0.432 | 0.488 |
| lab_L / lab_a / lab_b | 210 | 0.506 / 0.550 / 0.318 | 0.452 / 0.485 / 0.357 |
| chroma / sat / bright (approximations) | 210 | 0.403 / 0.342 / 0.614 | 0.388 / 0.368 / 0.608 |
| 13 texture traits | 210 | 0.316-0.529 | 0.323-0.476 |

Other checks on old vs new Copper data:

- Image level: old 2,398 images, new 2,397, common 2,397. Per-image median `Shape_Area` correlates r = 0.9997. Per-image object count equals in 26% of images;
  the new/old count ratio has median 0.97, mean 0.90, min 0.02. Old total objects 211,800 vs new 155,470.
- Median colour values differ in level (`lab_b` median 2.7 old, 18.3 new; chroma 10.2 vs 26.3). The statistic and possibly the pipeline changed.
- Object-level and well-level joins (image + `Object_Label`, or image + old `well_position` vs new grid index) did NOT give a reliable one-to-one match
  (best well join: area equal within 0.1% for 5.5% of matched wells; same strain string for 41%). I could not tell whether this is a wrong
  join key, re-segmentation, or a changed well-to-strain assignment. This is unresolved. The mismatches are systematic, not random: e.g. old `JES-025-001`
  vs new `DBVPG_7539` (517 wells), old `DBVPG_10619` vs new `DBVPG_6741` (493), old `TFCN_17-332C-2` vs new `DBVPG_6093`/`DBVPG_3772`. That pattern points to
  strain re-assignment of some wells (or a wrong well index in my join), not noise. Full list: `duckdb_retool/object_match.log`.

Conclusion: the new traits are NOT drop-in replacements. Do not mix the new table with GEMMA results from the old table.
Before any GWAS rerun, someone who knows the new pipeline should confirm the well-to-strain assignment and the colour statistic.
Reported correlations are between per-strain summaries with few replicate wells per strain and dose, so some disagreement is expected
from noise alone. I did not measure replicate reliability, so I cannot say how much.

## 6. Open items for a decision

1. Is `TFCN_223D-8` (new) = `TFCN_223A-8` (VCF)? Needed to keep the 213-strain panel.
2. Review the 23 new fuzzy/unmatched strains and the 54 collisions (`strain_reconciliation_duckdb/NEEDS_REVIEW.md`).
3. Confirm how the new data assigns strains to wells, and why old and new per-strain Cu summaries differ (r 0.2-0.6).
4. Decide whether `ColorLab_*GeoMedian` (and the RobustMean HSV columns) are acceptable replacements for the old `*Median` traits.
5. `common.py`, `build_gwas_phenotypes.py`, `check_mas_gates.py` still need `db_extract.parquet`. A DuckDB-based `read_extract()` would need
   the column mapping in section 5. I did not make it, because the comparison above says the trait values change.
6. The Cr unit (0-1.2) and the Zn design (doses 15-30 were run on different strains than 0/5/10) matter if the GWAS is extended to other metals.
