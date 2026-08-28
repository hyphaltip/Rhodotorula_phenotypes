# Provenance: copper-heavy-metal-screen-v0.15.1

## Source

**Type**: collaborator (shared lab project, same lab)

**Origin**:
- Derived from: `ArrayedHeavyMetalScreen` shared project's `0.15.1_Analysis` Quarto
  pipeline (R/tidyverse growth-curve fitting + AUC integration), at
  `/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Phenotyping/Heavy_Metals/ArrayedHeavyMetalScreen/analysis/0.15.1_Analysis/`.
  Underlying raw images/segmentation are the same imaging runs (d000353-d000357)
  already ingested into this project as `copper-colony-measurements`
  (`data/metadata/copper-colony-measurements/`); this dataset is a *re-processing*
  of that same experiment by the source lab's own pipeline, not new wet-lab data.

**Citation / accession**: N/A (internal shared lab storage)

## Acquisition details

**Date acquired**: 2026-08-27

**Obtained by**: Jason Stajich (via Claude Code, at user request)

**Method**: Files copied directly from
`.../0.15.1_Analysis/Results/{Copper subfolder small files, and sibling Results/
comparative_tolerance_*.csv / copper_auc_*.csv / copper_radial_growth_rates_*.csv /
copper_texture_*.csv / copper_config_pairwise_correlation.csv}`. Large per-colony
files (`Copper/` subfolder, `coppermeasurementmeta.csv`, `copper_measurements_combined.csv`,
totaling ~1.4 GB) were **symlinked, not copied**, following the precedent set by
`RmucY2510_v2-genotypes` (see its `MANIFEST.yaml`) -- see
`data/raw/copper-heavy-metal-screen-v0.15.1/large_source/`.

**Checksum**: SHA256 recorded for the 4 primary strain-level tables:
- `copper_auc_mean_by_strain.csv`: `db28b3d40eb70cd00f10aee0e1a034a813ac4b873e914e7eef7d9f76b79842fb`
- `copper_auc_all_strains_by_configuration.csv`: `436cf3b1bf9a12563bce1ac96e78f70ca991dea419a8ac95e64def4be2a6ec4a`
- `comparative_tolerance_all_strains_all_vs_all.csv`: `892703e5f4dca0ee20c014b7d964bf97f0fd10fea97d1023d379210a34175cf3`
- `copperterminal_dropped.csv`: `983f58894e495ebc0e22d6db1c82cdef0a43709f9bed462f8afcf1482e90fe3c`

(Remaining small files not individually hashed; re-derivable by re-copying from
the source path above, which the source lab controls and may update.)

## Access restrictions

**Restriction level**: institutional-only

**Details**: Shared lab storage (`stajichlab` group), same access tier as this
project's other shared-storage datasets (e.g. `RmucY2510_v2-genotypes`).

## Known issues

- Source directory (`Results/`) has duplicate `copper_measurements_combined.csv`
  and `copper_measurements_combined(1).csv` that are NOT byte-identical (differ
  at line 43); only the non-`(1)` file was symlinked -- provenance of the `(1)`
  variant is unclear and was not resolved with the source lab.
- The all-strain per-concentration linear growth-rate table that the source
  pipeline's `CopperAUCTesting.r` reads as input
  (`outputs/radial_growth_linear_all_strains.csv`) was not found alongside the
  files inspected in `Results/`; only the downstream, already-dose-integrated AUC
  tables were available. `copper_radial_growth_rates_all_concentrations.csv` /
  `_30mM.csv` (28 / 8 rows) appear to be an earlier single-strain (strain 254)
  prototype, not the all-strain rate table -- see source repo decision D-019
  referenced in `CopperAUCTesting.r`'s header comment.
- Source pipeline tried 3 growth-curve model families (logistic -> power-law ->
  linear-rate; source decisions D-015/D-016/D-018) before settling on linear
  rate as canonical; this dataset only reflects the final (linear-rate) choice.
- This dataset's AUC is **area under growth-rate-vs-concentration** (a
  dose-response-curve summary), which differs conceptually from this project's
  existing `AUC_0/AUC_10/AUC_20/AUC_30` GWAS phenotypes (area under
  colony-area-vs-time at one fixed dose). Do not treat as directly comparable
  without checking definitions -- see `COPPER_HEAVY_METAL_SCREEN_V0_15_1.md`.
- Strain identity uses the shared project's integer `Strain ID` / free-text
  `Strain` code, not this project's `strain_code` -- must be reconciled via
  the join used in `analysis/gwas/results/strain_reconciliation/` before use
  in any GWAS-adjacent analysis.

## Contact

**Primary contact**: ArrayedHeavyMetalScreen shared-project maintainer (see
`.claude/` / `.qwen/` project metadata at the source path for current owner;
not independently confirmed here).

**Backup contact**: None.

## Version history

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-08-27 | Initial ingestion of 0.15.1_Analysis copper outputs (AUC, radial growth rate, texture-outlier, comparative-tolerance, terminal per-colony QC-dropped tables). |
