# Heavy-metal array overview plots

Source: `heavy_metal_measurement` (Parquet in `data/preprocessed/heavy_metal_array/`, 517,371 rows; row counts per metal asserted)
and `strain_info` (DB, read-only). Run: `analysis/heavy_metal_overview/run.sh` inside a SLURM job. Log: `results/logs/make_overview.log`.
Figures: `results/figures/`. Tables: `results/*.csv`.

## Decisions

- **Time.** The new data has no `hours_since_plate_start`. Hours = `capture_datetime` minus the first capture of the same `(Metal, run_number)`.
  Run-start is used, not plate-start, because plate 14 of Cu run d000353 has a missed first capture (3 h late). (The GWAS builder uses plate-start to match the old definition.)
  Max hours per run: Cu 114-118, Fe 84-96, Pb 108-114, Zn 108-114, Cr 114-119.
- **Matched timepoint.** T = 90 h. For each plate, take the capture nearest 90 h; keep the plate only if within +/-3 h.
  Why 90: Pb has no capture between about 69 and 87 h, and Fe runs end at 84-96 h. Result at 90 h: 410 of 429 plates kept
  (excluded: Fe 16, Pb 2, Zn 1). Max deviation kept 2.85 h, median 0.26 h. Sensitivity: T = 66 h (425 of 429 plates).
- **Replicates.** Several objects in one well at one capture (242 wells of 26,135) are collapsed to their median. Strain value at a dose = median over wells.
  Most strain-dose cells have one well; replicate reliability was not measured.
- **Exclusions.** Rows with NULL `strain_id` (544 at the 90 h capture; 10,921 overall) cannot be assigned to a strain: excluded from strain plots, shown in the coverage plot title.
  Control-N strains are excluded from (b)-(e). Species panels (e) also exclude 'Species Not Found'
  (1,057 tolerance rows -> 1,002 after exclusions). Species labels are read from `strain_info` at run time and may change as the species overrides are updated.
- **Cr** has its own concentration axis (0-1.2); Cr is never pooled. Units of all concentrations are not documented. Area is `Shape_Area` (pixels assumed).
- **Colour.** All 5 metals have non-null `ColorLab_{L*,a*,b*}GeoMedian` (checked: counts equal row counts). Plotted for all five.
- **Tolerance** = log2(median area at top dose / median area at 0 dose), same strain, at T. Top dose = highest dose paired with 0 dose in at least 30 strains:
  Cu 30, Fe 30, Pb 30, Cr 1.2, **Zn 10** (Zn doses 15-30 were run on different strains than 0/5/10; at T=66 the top dose for Zn is 15).
  Strains with no detected colony (area missing or 0) at either dose are dropped, not set to zero. This biases the ratio upward at toxic doses (censoring).
  Strains used: Cu 314, Fe 141, Pb 294, Zn 80, Cr 308.

## Results (computed)

- Median area falls with dose in every metal (b). At the top dose the median relative area is about 0.5 (Cu), 0.66 (Fe), 0.05 (Pb), 0.41 (Zn, dose 10), 0.06 (Cr).
  Cu and Fe show a rise at dose 5 (about +20% Cu, +5% Fe).
- Median log2 tolerance at 90 h: Cu -0.97, Fe -0.59, Pb -4.23, Zn -1.28, Cr -4.15.
- Cross-metal Spearman (d1; pairwise n in `tolerance_pairwise_n.csv`): all |rho| below 0.25 except Cu-Zn 0.41 (n=80). Cu-Fe -0.22, Cu-Pb -0.01, Cu-Cr 0.10, Pb-Cr 0.14.
  No p-values or CIs were computed. Weak correlation can be real or come from single-well noise; not separated here.
- Time sensitivity (T=66 vs 90): Spearman Cu 0.95, Fe 0.84, Pb 0.77, Zn 0.82 (different top dose), Cr 0.85 (`tolerance_timepoint_sensitivity.csv`).

## Figures

| File | Content |
|---|---|
| a1_coverage_strain_by_metal.png | strain x metal colony counts; strains and runs per metal |
| a2_coverage_runs_by_metal.png | distinct strains per run x concentration |
| b_dose_response_area.png | area vs dose, median and IQR across strains; absolute and relative to 0 |
| c_dose_response_colour.png | L*, a*, b* (GeoMedian) vs dose |
| d1_tolerance_correlation.png, d2_tolerance_scatter_matrix.png | cross-metal tolerance |
| e_tolerance_by_species.png | tolerance by species (boxes for >=3 strains per metal) |

Not verified: figures d2 and e and c were generated but only a1 and b were inspected by eye.
