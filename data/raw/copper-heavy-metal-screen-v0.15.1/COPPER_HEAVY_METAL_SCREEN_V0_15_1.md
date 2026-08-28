# copper-heavy-metal-screen-v0.15.1

Updated copper phenotyping outputs from the shared lab's `0.15.1_Analysis` Quarto
project (`ArrayedHeavyMetalScreen`), copied 2026-08-27. Same underlying imaging
runs as the existing `copper-colony-measurements` dataset (runs d000353-d000357,
211,800 per-colony rows match exactly) but reprocessed through a newer pipeline
version that adds: (1) a "terminal" per-colony extract split by plate arrangement
and concentration, (2) a **linear radial-growth-rate model** (Shape_MaxRadius ~
Hours, OLS) per (strain, concentration, plate configuration) superseding this
project's earlier logistic/power-law fits, (3) a **dose-response AUC** — trapezoidal
area under growth-rate-vs-concentration (0-30 mM), per strain and per configuration
plus a configuration-mean, and (4) a **multi-metal comparative tolerance** join
(copper/chromium/iron/pH mean AUC-rate per strain with top/bottom-quartile flags
and a polytolerant/single-condition/susceptible classification).

**This is a different AUC axis than this project's existing `AUC_0/AUC_10/AUC_20/
AUC_30` GWAS phenotypes** (`analysis/gwas/scripts/build_gwas_phenotypes.py`), which
integrate colony area over *time* at a *fixed* copper dose. The AUC here integrates
growth *rate* over *dose* (0-30 mM) at a fixed endpoint — i.e. a single-number
dose-response-curve summary per strain, closer in spirit to the existing `cu_dose_slope`
trait but computed from an independent (OLS-linear, not per-timepoint-difference)
growth model. Treat as a complementary/candidate replacement trait, not a
drop-in substitute — see `PROVENANCE.md`-equivalent notes in `data/metadata/
copper-heavy-metal-screen-v0.15.1/provenance.md` for the exact caveats before
using in GWAS.

## Layout

```
copper_auc_all_strains_by_configuration.csv     # AUC(rate vs Cu dose) per (strain, plate config); 1577 rows
copper_auc_mean_by_strain.csv                   # AUC averaged across configs; 436 strains
copper_auc_strain254_by_configuration.csv       # single-strain AUC prototype/QC check (D-019 in source repo)
copper_radial_growth_rates_all_concentrations.csv  # linear-fit y0/mumax/K-style params per strain x conc x config (small sample; see known issues)
copper_radial_growth_rates_30mM.csv             # same, 30 mM only
copper_config_pairwise_correlation.csv          # correlation of rate-vs-conc curves across the 4 plate replicate configurations
copper_texture_distinct_strains_summary.csv     # Haralick-texture outlier distance-from-centroid, one row/strain
copper_texture_outlier_strains.csv              # full texture feature vector for outlier strains
copper_texture_outliers_all_runs.csv            # outlier calls per run
comparative_tolerance_*.csv                     # copper vs chromium/iron/pH mean-AUC-rate cross-metal comparison + quartile membership (2 schema versions, v2 adds `source`/quartile_direction long format)
copperterminal_dropped.csv                      # QC: wells dropped at each run's last measured hour (per plate/concentration)
large_source/                                   # symlinks to shared storage; NOT copied into git (see below)
  Copper_terminal_by_concentration/  -> .../Results/Copper/   (115 MB: per-colony terminal extract, all + split by arrangement x concentration)
  coppermeasurementmeta.csv          -> .../Results/coppermeasurementmeta.csv   (643 MB: full time-course, metadata-joined)
  copper_measurements_combined.csv   -> .../Results/copper_measurements_combined.csv  (639 MB: raw input to the R growth-model scripts)
```

## Known issues

- `copper_measurements_combined.csv` and `copper_measurements_combined(1).csv` exist
  side-by-side in the source directory and are **not byte-identical** (first
  difference at line 43); only the non-`(1)` file is symlinked here. Unclear from
  the source project which is canonical — flagged to the source team, not resolved.
- `copper_radial_growth_rates_all_concentrations.csv` / `_30mM.csv` have only 28/8
  data rows respectively (single-strain, strain 254 only) despite the AUC tables
  covering ~436 strains — these two files look like the strain-254 rate-model
  prototype (source repo decision D-019), not the all-strain rate table. The
  all-strain linear rates that actually feed the AUC integration
  (`outputs/radial_growth_linear_all_strains.csv` per `CopperRadialGrowthLinear.R`)
  were **not present** in the `Results/Copper/` or `Results/` directories examined;
  only the downstream `copper_auc_*` tables (which already fold in per-strain rate)
  were available. If per-concentration rate (not just dose-integrated AUC) is
  needed, that intermediate file should be requested from the source team.
- Source repo's `CopperAUCTesting.r` documents 3 successive growth-model choices
  (logistic mumax -> power-law alpha -> linear rate, decisions D-015/D-016/D-018)
  before settling on linear rate as canonical; only the linear-rate-based AUC is
  provided here.
- `copperterminalmeta_arrangement{1..4}_concentration{0,5,10,15,20,25,30}.csv` and
  the non-split `copperterminalmeta_arrangement{1..4}.csv` / `copperterminalmeta.csv`
  in `large_source/Copper_terminal_by_concentration/` are heavily overlapping views
  of the same terminal-timepoint rows sliced different ways (by arrangement, by
  concentration, or both) — not independent data.
- Large per-colony files (`large_source/`) are symlinks into shared lab storage
  (`/bigdata/stajichlab/shared/.../0.15.1_Analysis/Results/`), consistent with how
  `RmucY2510_v2-genotypes` was handled — not duplicated, not git-tracked, and only
  as durable as that shared path.

## Suggested next step

Compare `copper_auc_mean_by_strain.csv` (dose-response AUC, this dataset) against
this project's existing `cu_dose_slope` / `AUC_ratio_10` / `resilience_30` GWAS
traits (`analysis/gwas/scripts/build_gwas_phenotypes.py`) for the strain subset in
common, to decide whether to add this as a new GWAS trait or use it as an
independent validation/replication check on the existing copper-response GWAS hits.
Strain identity here uses `Strain ID` (shared-project integer, e.g. `254`) and a
free-text `Strain`/`Strain Name` — reconcile against this project's `strain_code`
via the same join used in `analysis/gwas/results/strain_reconciliation/` before
comparing.
