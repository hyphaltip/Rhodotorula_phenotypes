# Summary Statistics: copper-heavy-metal-screen-v0.15.1

<!-- Generated: 2026-08-27 -->
<!-- Script: manual (python3 csv.DictReader inline, see session) -->

## Overview

| Property | Value |
|----------|-------|
| Rows (copied files, 16 total) | see per-file table below |
| Large symlinked files | Copper_terminal_by_concentration/ (115 MB, 40 files), coppermeasurementmeta.csv (211,800 rows / 643 MB), copper_measurements_combined.csv (639 MB) |
| File size (copied files only) | ~925 KB |
| Date range | Underlying imaging: 2026-02-25 onward (run d000353-d000357, matches copper-colony-measurements); analysis outputs generated/exported 2026-08-27 |
| Format | CSV |

## Per-file row counts (logical CSV records, quoted-newline-safe)

| File | Rows | Unit of observation |
|------|------|----------------------|
| copper_auc_mean_by_strain.csv | 298 strains | strain (averaged across configurations) |
| copper_auc_all_strains_by_configuration.csv | 1081 | strain x plate configuration |
| copper_auc_strain254_by_configuration.csv | 4 | single strain (254) x configuration -- QC/prototype check |
| copper_radial_growth_rates_all_concentrations.csv | 28 | strain 254 x concentration x configuration only (NOT all-strain, see known issues) |
| copper_radial_growth_rates_30mM.csv | 8 | strain 254 x configuration, 30 mM only |
| copper_config_pairwise_correlation.csv | 7 | one row per concentration, 6 pairwise-config correlations + mean |
| copper_texture_distinct_strains_summary.csv | ~26 | strain (outlier distance summary) |
| copper_texture_outlier_strains.csv | ~28 | strain (full texture feature vector, outliers only) |
| copper_texture_outliers_all_runs.csv | 317 | per-run outlier calls |
| comparative_tolerance_all_strains_all_vs_all.csv | 306 strains | strain, 4-metal wide join |
| comparative_tolerance_extreme_strains_all_vs_all.csv | subset of above | strain (extreme only) |
| comparative_tolerance_top_bottom_quartile_membership.csv | 306 x up to 4 metals | strain (wide, one row/strain) |
| comparative_tolerance_top_bottom_quartile_membership_v2.csv | long format | strain x metal (one row per metal per strain) |
| comparative_tolerance_top_quartile_membership.csv | subset | strain, top-quartile only |
| comparative_tolerance_top_quartile_membership_v2.csv | subset | strain x metal, top-quartile only |
| copperterminal_dropped.csv | 1197 | well (QC-dropped wells, one row per dropped well) |

## Column summaries (key strain-level AUC columns)

| Column | Type | Non-null | Min | Max | Mean | SD |
|--------|------|----------|-----|-----|------|----|
| copper_auc_mean_by_strain.csv: mean_auc_rate | numeric | 298 | 0.806 | 29.835 | 22.159 | 4.116 |
| copper_auc_all_strains_by_configuration.csv: auc_rate | numeric | 1081 | 0.806 | 31.190 | 22.330 | 4.780 |
| comparative_tolerance_all_strains_all_vs_all.csv: mean_auc_rate_copper | numeric | 298/306 | 0.806 | 29.835 | 22.159 | 4.116 |

`profile_class` distribution (comparative_tolerance_all_strains_all_vs_all.csv, n=306):
Mixed (top in some, bottom in others)=100, Polysensitive (bottom in >=2)=54,
Single-condition extreme (top)=42, Polytolerant (top in >=2)=41,
Single-condition extreme (bottom)=31, Not extreme (middle 50% in all)=38.

## Missing data summary

- `comparative_tolerance_all_strains_all_vs_all.csv`: 306 strains total but only
  298 have non-null `mean_auc_rate_copper` -- 8 strains present in the cross-metal
  join lack a copper AUC (likely screened for another metal but not copper, or
  failed all-configuration QC for copper). Not further diagnosed.
- No other missingness quantified in this pass; copied files are small enough to
  inspect directly if needed.

## Quality flags

- `copper_radial_growth_rates_all_concentrations.csv` / `_30mM.csv` cover only
  strain 254, not all ~298-436 strains -- do not use these two files as an
  all-strain rate table (see Known issues in `COPPER_HEAVY_METAL_SCREEN_V0_15_1.md`).
- `copper_auc_mean_by_strain.csv` (298 rows) and `copper_auc_all_strains_by_configuration.csv`
  (1081 rows / up to 4 configs each) are internally consistent (identical
  min/mean-scale for mean_auc_rate vs auc_rate) -- no cross-file discrepancy found.
- Row-count mismatch between `wc -l` (437) and logical CSV records (298) for
  `copper_auc_mean_by_strain.csv` is expected and not a data-quality issue: the
  `Strain Name` field contains embedded newlines inside quoted CSV fields.

## Notes

Strain-level counts here (~298-306) are smaller than this project's existing
213-strain GWAS panel and the 422-sample genotype VCF, and use a different ID
system (`Strain ID` integer + `Strain` code) -- reconciliation against
`strain_code` via `analysis/gwas/results/strain_reconciliation/` is required
before any cross-dataset merge, and the overlap size is not yet known.
