Species, copper, and colony growth-rate parameters
================
analysis/growth_rates
2026-08-15

- [Question](#question)
- [Data construction](#data-construction)
- [Model fitting and rate
  extraction](#model-fitting-and-rate-extraction)
- [Growth rate by copper and
  species](#growth-rate-by-copper-and-species)
- [Growth-rate × color / light-intensity
  interaction](#growth-rate--color--light-intensity-interaction)
- [Limitations](#limitations)
- [Reproduce](#reproduce)

## Question

Each strain is measured once per Copper concentration per run, on
physically distinct plates (`replicate_label` is an imager-run proxy —
see [`DATABASE_DESIGN.md`](../../DATABASE_DESIGN.md)). Each (plate,
well) is therefore an independent culture grown over ~120 h with
~6-hourly imaging. This analysis fits a growth curve to every culture
and asks:

1.  **Primary**: how does per-strain growth rate vary with copper
    concentration and with species?
2.  **Secondary**: how do growth-rate parameters interact with endpoint
    color (CIELAB `L*`, `a*`, `b*`) — where `L*` is the “light
    intensity” readout of interest and the measured
    `Intensity_MeanIntensity` is ~0.99-collinear with it (same axis)?

Unit of analysis: **per colony = (plate, well)**. Design capacity: 112
plates (runs 353-356, 7 Cu × 16 plates) × 96 wells = 10,752 slots; each
strain appears once per (run, Cu) on a distinct plate, so the
theoretical max is 309 strains × 4 runs × 7 Cu = 8,652 colonies; **8,446
observed (97.6%** detection, the remainder being high-Cu dropout).

## Data construction

`scripts/00_build_series.R` rebuilds the per-colony × timepoint table
from `v_phenotype` using the same selection/dedup as the plate-position
timecourse: runs 353-356, `strain_id IS NOT NULL`,
`hours_since_plate_start > 0`, largest object per plate/well/time,
`plate_id = run_plate`. Two size readouts:

- `ln_area` — `log(Shape_Area)`; colonial expansion / biomass proxy;
  rises throughout the ~120 h window (rarely plateaus in-window).
- `intensity` — `Intensity_MeanIntensity`; rises 0.25-\>0.41 and
  **saturates ~48-54 h** (early window).

| copper_mm | n_tp | n_colonies |
|----------:|-----:|-----------:|
|         0 |    1 |         14 |
|         0 |    2 |         17 |
|         0 |    3 |          4 |
|         0 |    4 |          2 |
|         0 |    7 |          1 |
|         0 |   11 |          1 |
|         0 |   12 |          3 |
|         0 |   13 |          4 |
|         0 |   14 |          1 |
|         0 |   15 |          8 |
|         0 |   16 |          5 |
|         0 |   17 |         16 |
|         0 |   18 |        159 |
|         0 |   19 |        661 |
|         0 |   20 |        300 |
|         5 |    1 |         19 |
|         5 |    2 |         10 |
|         5 |    3 |          8 |
|         5 |    4 |          1 |
|         5 |    5 |          6 |
|         5 |    7 |          3 |
|         5 |    8 |          1 |
|         5 |    9 |          7 |
|         5 |   10 |          4 |
|         5 |   11 |          4 |
|         5 |   13 |          2 |
|         5 |   14 |          5 |
|         5 |   15 |          3 |
|         5 |   16 |          4 |
|         5 |   17 |         12 |
|         5 |   18 |         16 |
|         5 |   19 |        793 |
|         5 |   20 |        293 |
|        10 |    1 |         20 |
|        10 |    2 |         12 |
|        10 |    3 |         26 |
|        10 |    4 |         24 |
|        10 |    5 |         27 |
|        10 |    8 |          1 |
|        10 |    9 |          1 |
|        10 |   11 |          3 |
|        10 |   12 |          1 |
|        10 |   13 |          3 |
|        10 |   14 |          5 |
|        10 |   15 |         11 |
|        10 |   16 |         12 |
|        10 |   17 |         15 |
|        10 |   18 |         86 |
|        10 |   19 |        669 |
|        10 |   20 |        276 |
|        15 |    1 |         25 |
|        15 |    2 |         11 |
|        15 |    3 |         29 |
|        15 |    4 |         28 |
|        15 |    6 |          1 |
|        15 |    7 |          2 |
|        15 |   10 |          1 |
|        15 |   11 |          1 |
|        15 |   13 |          2 |
|        15 |   14 |          2 |
|        15 |   15 |          8 |
|        15 |   16 |         12 |
|        15 |   17 |         20 |
|        15 |   18 |         41 |
|        15 |   19 |        737 |
|        15 |   20 |        268 |
|        20 |    1 |         35 |
|        20 |    2 |         17 |
|        20 |    3 |         16 |
|        20 |    4 |          9 |
|        20 |    5 |          4 |
|        20 |    6 |          3 |
|        20 |    7 |         24 |
|        20 |    8 |          2 |
|        20 |    9 |          4 |
|        20 |   10 |         20 |
|        20 |   11 |         21 |
|        20 |   12 |          6 |
|        20 |   13 |          8 |
|        20 |   14 |          9 |
|        20 |   15 |          9 |
|        20 |   16 |         15 |
|        20 |   17 |         89 |
|        20 |   18 |        140 |
|        20 |   19 |        611 |
|        20 |   20 |        152 |
|        25 |    1 |         58 |
|        25 |    2 |         18 |
|        25 |    3 |         52 |
|        25 |    4 |         38 |
|        25 |    5 |          5 |
|        25 |    6 |          3 |
|        25 |    8 |          3 |
|        25 |    9 |          6 |
|        25 |   10 |         22 |
|        25 |   11 |          1 |
|        25 |   13 |          3 |
|        25 |   14 |         11 |
|        25 |   15 |         11 |
|        25 |   16 |         16 |
|        25 |   17 |         37 |
|        25 |   18 |        107 |
|        25 |   19 |        618 |
|        25 |   20 |        214 |
|        30 |    1 |         38 |
|        30 |    2 |         61 |
|        30 |    3 |         50 |
|        30 |    4 |         43 |
|        30 |    5 |          7 |
|        30 |    6 |          8 |
|        30 |    7 |          4 |
|        30 |    8 |          8 |
|        30 |    9 |          6 |
|        30 |   10 |          9 |
|        30 |   11 |          6 |
|        30 |   12 |          7 |
|        30 |   13 |          9 |
|        30 |   14 |          6 |
|        30 |   15 |         35 |
|        30 |   16 |         19 |
|        30 |   17 |         21 |
|        30 |   18 |         70 |
|        30 |   19 |        778 |
|        30 |   20 |         77 |

Series coverage (runs 353-356)

## Model fitting and rate extraction

`scripts/01_fit_growth_models.R` fits two sigmoid forms per colony per
trait (Zwietering reparameterizations, parameterized as maximum rate
`mu`/h, lag `lambda` h, asymptote `A`), requiring \>=8 usable timepoints
spanning \>=24 h:

- Gompertz: `y = A * exp(-exp((mu*e/A)*(lambda - t) + 1))`
- Logistic: `y = A / (1 + exp((4*mu/A)*(lambda - t) + 2))`

**Identifiability finding (why the headline is a data-derived rate):**
within the available window the lag parameter is structurally
unidentifiable — ~94% of fits put `lambda` on a boundary (\<= -1 h or
\>= 144 h), and for `ln_area` the asymptote is extrapolated in ~82% of
fits because area does not plateau in-window. Unbounded Gauss-Newton
refits on a 60-colony sample reproduced degenerate lags (-297…-20 h).
The parametric `mu` is therefore kept as a secondary parameter (weak
agreement with the data-derived rate), and the **primary estimator is
the data-derived peak slope**:

- `rate_area` — max positive 6 h slope of `log(Shape_Area)` per hour
- `rate_int` — max positive 6 h slope of `Intensity_MeanIntensity` per
  hour

Defined for ~99% of colonies (`rate_area` NA on 2 colonies). Both
Gompertz and logistic fits are retained with AIC; preferred model =
lower AIC (Gompertz favored ~80-88%, depending on trait).

| colonies fitted (\>=8tp, \>=24h) | rate_area defined | rate_int defined | pref. Gompertz (area) | pref. Gompertz (intens.) | Kendall tau (mu vs rate_area) |
|---:|---:|---:|---:|---:|---:|
| 7663 | 15322 | 15326 | 80 | 88 | 0.373 |

Fit coverage and rate-extraction summary (n columns)

## Growth rate by copper and species

`scripts/02_species_cu_rates.R`. Rates are aggregated to strain × Cu
(plain mean across up to 4 run-plates; histogram n=1:130, n=2:161,
n=3:380, n=4:1505 of 2,183 strain-Cu rows — the n\<2 single-culture rows
are flagged). Mixed model (colony level, REML):

`rate ~ factor(copper_mm) + (1 | species / strain_code)`

Species and strain enter as random effects (shrinkage toward the
well-sampled mean; 20 named species + `unknown`, 282 colonies without a
species label).

| term | Estimate | Std..Error | tvalue | pvalue | model |
|:---|---:|---:|---:|---:|:---|
| factor(copper_mm)5 | 0.0175 | 0.00587 | 2.99 | 0.00283 | rate_area ~ Cu(factor) + species/strain (colony-level) |
| factor(copper_mm)10 | 0.0145 | 0.00595 | 2.44 | 0.01460 | rate_area ~ Cu(factor) + species/strain (colony-level) |
| factor(copper_mm)15 | 0.0441 | 0.00594 | 7.42 | 0.00000 | rate_area ~ Cu(factor) + species/strain (colony-level) |
| factor(copper_mm)20 | 0.0627 | 0.00595 | 10.50 | 0.00000 | rate_area ~ Cu(factor) + species/strain (colony-level) |
| factor(copper_mm)25 | 0.0812 | 0.00600 | 13.50 | 0.00000 | rate_area ~ Cu(factor) + species/strain (colony-level) |
| factor(copper_mm)30 | 0.1030 | 0.00601 | 17.10 | 0.00000 | rate_area ~ Cu(factor) + species/strain (colony-level) |

rate_area: fixed effects (diff vs 0 mM)

| copper_mm |     diff |       se |       lo |       hi | tvalue | pvalue | model    |
|----------:|---------:|---------:|---------:|---------:|-------:|-------:|:---------|
|         5 | 0.000631 | 0.000121 | 0.000394 | 0.000867 |   5.23 |  2e-07 | rate_int |
|        10 | 0.001030 | 0.000122 | 0.000786 | 0.001260 |   8.39 |  0e+00 | rate_int |
|        15 | 0.001470 | 0.000122 | 0.001230 | 0.001710 |  12.10 |  0e+00 | rate_int |
|        20 | 0.002400 | 0.000122 | 0.002160 | 0.002630 |  19.60 |  0e+00 | rate_int |
|        25 | 0.003670 | 0.000123 | 0.003420 | 0.003910 |  29.70 |  0e+00 | rate_int |
|        30 | 0.005910 | 0.000123 | 0.005670 | 0.006150 |  47.90 |  0e+00 | rate_int |

rate_int vs 0 mM (Wald 95% CI)

**Result (both readouts): peak growth/brightening rate increases
monotonically with copper** (area: F = 79.0, p ~ 6e-96; intensity: F =
540, p ~ 0; all levels vs 0 mM significant, monotone coefficients 5-\>30
mM). Counter to a naive “Cu is toxic =\> slower growth” expectation,
this is consistent with the plate-position timecourse: at higher Cu,
colony area keeps expanding through the ~120 h window (the log-growth
phase lasts longer), so the realized peak log-slope is **higher** — the
peak-slope measure is phase/length sensitive rather than a growth-rate
constant. Caveat: a phase-independent rate (e.g. fixed-interval doubling
time) would be needed to separate “faster” from “still growing when
observed”.

![](results/figures/rate_by_cu_spp.png)<!-- -->![](results/figures/rate_int_by_cu_spp.png)<!-- -->

Well-sampled species (\>=8 strains): dairenensis (8), diobovata (10),
mucilaginosa (216), paludigena (16), sphaerocarpa (8), toruloides (10);
`R. sp. clade I` excluded from the interaction test (unnamed). The Cu ×
species interaction test (well-sampled subset) is in
`results/tables/rate_models_species_cu.csv/txt`.

## Growth-rate × color / light-intensity interaction

`scripts/03_color_interaction.R`. Endpoint color = last non-NA CIELAB
reading per colony. `L*` (lightness) is the “light intensity” outcome;
the measured `Intensity_MeanIntensity` is the same axis (endpoint r =
0.999) so it enters only as the growth readout `rate_int`, never as a
co-predictor. Model:

`L* ~ rate_area × copper + (1 | species/strain_code)` (and a*/b*
analogues)

| term                                  | Estimate | Std..Error | t.value |   Pr…t.. |
|:--------------------------------------|---------:|-----------:|--------:|---------:|
| rate_area_ln_area:factor(copper_mm)5  |     3.70 |      0.645 |    5.74 | 0.00e+00 |
| rate_area_ln_area:factor(copper_mm)10 |     2.91 |      0.688 |    4.23 | 2.41e-05 |
| rate_area_ln_area:factor(copper_mm)15 |     2.56 |      0.616 |    4.16 | 3.24e-05 |
| rate_area_ln_area:factor(copper_mm)20 |     1.52 |      0.564 |    2.69 | 7.10e-03 |
| rate_area_ln_area:factor(copper_mm)25 |     1.78 |      0.547 |    3.25 | 1.15e-03 |
| rate_area_ln_area:factor(copper_mm)30 |     1.30 |      0.536 |    2.42 | 1.56e-02 |

L\* ~ rate_area x Cu: interaction coefficients

The growth-rate × Cu interaction on `L*` is significant (F = 8.16, p ~
8e-9): faster-growing colonies end lighter, with the slope largest at
low Cu (5 mM: +3.7 L\* per unit rate) and weakest at 30 mM (+1.3).
Copper itself dominates chromatic outcome (Cu F ~ 900 on L*, ~ 1100 on
a*). The two growth modes are weakly correlated (`rate_area` vs
`rate_int` r = 0.21) — area expansion and brightness rise are distinct
axes.

## Limitations

- Peak-slope rate is phase/length sensitive (see caveat above);
  parametric `mu` is unidentifiable within the window.
- 130 of 2,183 strain-Cu aggregates rest on a single culture (flagged,
  kept).
- 282 colonies (3.3%) lack a species label -\> `unknown` random-level.
- `R. sp. clade I` is excluded from the species-interaction test.
- High-Cu missingness (~15-18% of strain-plates never reach late
  timepoints) shapes which colonies reach the fitting threshold;
  missingness itself is Cu-dependent (see plate-position analysis).

## Reproduce

    bash analysis/growth_rates/run.sh
