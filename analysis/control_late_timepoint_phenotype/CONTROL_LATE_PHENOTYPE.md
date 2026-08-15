# Control-Media Late-Timepoint Strain Phenotype Table

**Purpose.** A per-strain phenotype table for the control media only (Copper
concentration = 0 mM, i.e. growth on plain YPD), taken at a late timepoint
(hours_since_plate_start in [80, 110]) that has data for essentially all
strains. The phenotype emphasis is colony color (CIELAB L\*, a\*, b\*) plus
colony size, aggregated across every replicate colony of a strain.

## Inputs
- `db/rhodotorula_phenotypes.duckdb`
  - `v_phenotype` — colony measurements joined with strain, timepoint, and
    condition (96 replicate control plates at Cu = 0).
  - `condition_plate_factor` — factor `'Copper concentration'` (mM), filter to
    `factor_value = 0`.
  - `strain` — strain metadata: `strain_id`, `strain_code`, `strain_name`,
    `species`, `origin`, `environment`.

## Method
1. **Control filter** — only plates with `Copper concentration = 0 mM`.
2. **Late timepoint** — the imaging rig runs on two interleaved ~3 h cadences,
   so a fixed single hour would drop ~85 strains imaged only on the alternate
   cadence. Instead, hours are rounded to the nearest integer (the imaging
   "pass"), and each strain gets its **latest pass within the chosen time
   window**; all colonies from that pass are used as the strain's replicates.
   Three windows are produced:
   - `phenotype_control_timepoint_70_80.csv` — latest pass within [70, 80] h
     (passes 75 h / 78 h)
   - `phenotype_control_timepoint_80_90.csv` — latest pass within [80, 90] h
     (passes 87 h / 90 h)
   - `phenotype_control_timepoint_90_110.csv` — latest pass within [90, 110] h
     (passes 105 h / 108 h)
3. **Per-colony trait values** (one row of `v_phenotype` = one colony at one
   timepoint):
   | trait | column | unit |
   |---|---|---|
   | colony size | `Shape_Area` | px |
   | L\* | `ColorLab_L*Median` | CIELAB L\* (0–100) |
   | a\* | `ColorLab_a*Median` | CIELAB a\* |
   | b\* | `ColorLab_b*Median` | CIELAB b\* |
   Per-colony color values are the pixel-median of the colony's Lab histogram
   (robust to intra-colony outliers).
4. **Per-strain aggregation** across all replicate colonies at the strain's
   chosen pass: `median`, `mean`, sample `variance`, and sample `sd` for each of
   area / L\* / a\* / b\*; plus `n_colonies` and `n_replicate_wells`
   (distinct run:well spots) for transparency.

## Output — `results/phenotype_control_timepoint_{tmin}_{tmax}.csv`
One table per window `[tmin, tmax]` (70–80, 80–90, 90–110 h).
Columns: `strain_id, strain_code, strain_name, genus, species, origin,
environment, n_colonies, n_replicate_wells, timepoint_h` +
(for each of `area, l, a, b`) `_median, _mean, _var, _sd`.
- `genus` = first token of `species`.
- `timepoint_h` = rounded hour of the strain's latest pass in the window.
- `n_colonies` = number of colonies aggregated (typical 3–4; strains with
  `n_colonies < 3` have unreliable variance/SD and are flagged).

## Coverage (2026-08-15 re-run)
- Each window covers **314/320** strains; 286 strains with ≥3 colonies at
  their latest pass, 28 flagged (`n_colonies < 3`); 1 strain has no control
  (Cu = 0) image in any window and 5 have no late control image at all.
- Timepoint passes used per window: 75/78 h, 87/90 h, 105/108 h (84 alternate
  cadence + 230 main cadence strains each window). The strain sets are
  identical across windows; only the sampled pass differs.
- Range sanity: L\* 70–80 (yeast colony lightness), a\* ~0.6–14 (carotenoid
  red/orange), b\* −1.4–9.1 (yellowish), colony area ~0.8 k–70 k px.
  Colony size across windows is highly stable (area_median correlation
  70–80 vs 90–110 ≈ 0.98).
- Spot-check: per-strain stats reproduced exactly by an independent manual
  aggregation (strain 185).

## Caveats
- Variance/SD are computed on small samples (median 4 colonies) — treat as
  within-strain replicate spread, not population variance.
- Strains on the alternate cadence are measured ~3 h earlier than the main
  group (105 h vs 108 h); negligible for a late, plateau-stage timepoint but
  recorded per strain in `timepoint_h`.
- `strain_id = NULL` rows (unidentified spots) are excluded.
- Sample variance (`var_samp`/`stddev_samp`); use population versions if the
  colonies are considered the full population.

## Reproduce
```
bash analysis/control_late_timepoint_phenotype/run.sh
```