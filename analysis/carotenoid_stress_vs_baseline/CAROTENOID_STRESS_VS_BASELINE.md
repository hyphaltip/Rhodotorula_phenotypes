# a* (carotenoid proxy) vs colony size under metal stress

Question: does metal stress raise CIELAB a* (induced carotenoid), or are some strains inherently high-a* producers?

## Reproduce
`sbatch -p short -c 4 --mem=24G -t 60 --wrap="cd $PWD && bash analysis/carotenoid_stress_vs_baseline/run.sh"` (about 2 min). Outputs go to `results/` (gitignored, rebuildable).

## Design (all decisions)
- Source: `heavy_metal_measurement`, column `ColorLab_a*GeoMedian`, size = ln(`Shape_Area`).
- Dropped: NULL `strain_id` (10,921 rows) and `Control-N` (4,399 rows). 502,051 of 517,371 rows kept.
- Replicate unit = well (Metal, run, plate_position, Grid_RowNum, Grid_ColNum). Each plate has one dose. Cr, Cu, Fe, Pb: median 3-4 replicate wells per strain x dose, on different plates of one run. Zinc: 1 well per strain x dose.
- Several objects can occupy one well in one image (up to 28). The largest object per well-image is kept. A `--all-objects` sensitivity table is built (`wells_allobj.csv`) but NOT yet modelled.
- Time: hours since the plate's first image. Window = [T-24, T] h, T = median plate span per metal (Cr 114, Cu 114, Fe 89.9, Pb 107.5, Zn 107.6). A plate-relative window was rejected after review because spans differ between plates. Wells need >= 2 distinct window images: 269 of 27,631 wells dropped; the loss by dose is in `wells_lost_by_dose.csv` (highest at the top doses; Cr 1.2: 158 of 1,210).
- Models (lme4, REML), per metal, dose_s = conc / max conc: `a ~ dose_s [+ lnA_c] + (1 + dose_s | strain) + (1 | plate)`. M0 omits size (total effect). M1 adds size centred on the control mean (effect at the same size). M2 uses dose as a factor.
- Zinc is descriptive only: 9 plates, and dose is confounded with run (d000388 = doses 0-15, d000390 = 10-30). No Zinc strain has both 0 and top dose.

## Results (numbers from `mixed_model_summary.csv`, `dose_factor_effects_M2.csv`)
Dose effect on a*, 0 to top dose (SE in text of csv):

| Metal | M0 total | M1 at same size | strain repeatability at 0 / top dose (M1) |
|---|---|---|---|
| Chromium | -13.5 | -8.4 | 0.55 / 0.24 |
| Copper | -15.9 | -13.4 | 0.51 / 0.28 |
| Iron | +2.4 | +3.8 | 0.82 / 0.91 |
| Lead | -16.8 (singular fit) | -6.1 | 0.49 / 0.03 |

Dose-as-factor (M2, at fixed size, a* units vs 0 dose): Cr +2.2 (0.2), +3.1 (0.4), +1.2, +0.9, -3.3, -9.1 (1.2); Cu 0.05 (5), then falling to -12.3 (30); Fe +2.5 (5) rising to +4.6 (30); Pb +4.4 (5), +7.5 (10), +7.6 (15), then -0.9 to -1.4 (20-30). SE is 0.2-0.5.

## Reading
- Not a simple "stress raises a*". Cu falls monotonically; Fe rises monotonically; Cr and Pb rise at low-mid doses and fall at the highest doses (the linear-dose models hide this, so use M2).
- Fig 1: at the same colony size, Cu high-dose curves sit below the control curve, so the Cu drop is not only a size effect.
- Fe: inherent differences dominate (repeatability 0.8-0.9 at both ends), and strains with higher baseline a* gain more (corr 0.44).
- Cr, Cu, Pb: strain differences shrink at high dose (repeatability falls to 0.24, 0.28, 0.03): strains converge to a common low a*. Baseline and induction BLUPs correlate at -0.8 to -0.99. That is what convergence to a common value produces, and it is partly built into the model (Pb M1 is near the boundary). Do not read it as strain-specific biology.

## Caveats
- Size is partly an effect of stress, so M1 is a direct effect, not the whole induction.
- M1 uses one size slope for all strains and extrapolates it to the small colonies at high dose (0.0-3.9% of stressed colonies fall below the control 1st percentile size in the earlier pooled fit; not recomputed for wells). A strain whose a*-size relation differs can look like a baseline or induction effect.
- a* of very small colonies may include background pixels; not checked.
- Survivor selection: wells with no growth have no objects and are absent from the data. Top-dose estimates describe wells that grew.
- Cr dose scale (0-1.2) and units for all metals are not given in the source.
- Fe and Pb results include plates with short imaging spans (Fe 36-96 h); wells with < 2 window images were dropped.
- Not viewed or checked: Fig 2 and Fig 4 were generated but not inspected.
