# a* (carotenoid proxy) under metal stress: stratified analyses

Generated 2026-10-08. All numbers come from `report/tables/*.csv`; scripts are in `scripts/` and `scripts/strat/`; the data-quality caveats are in `analysis/heavy_metal_data_problems/HEAVY_METAL_DATA_PROBLEMS.md`.

**Question.** When Rhodotorula strains grow on plates with a heavy metal at increasing doses, does CIELAB a\* (red-green axis, used here as a carotenoid proxy) rise because of stress, or do strains differ inherently in a\*? Colony size confounds this, so every dose effect is shown with and without size.

## Summary

1. **There is no general "stress raises a\*".** The effect depends on the metal. Over the full dose range, at fixed colony size: Cr -8.4, Cu -13.4, Pb -6.1 (a\* units), Fe +3.8. Zn (-3.7, p = 0.44) is not informative (7 plates, one per dose).
2. **The dose-response is not linear.** Dose-as-factor models show Cr and Pb a\* *rise* at low to middle doses (Cr up to +3.1 at dose 0.4; Pb up to about +7.5 at doses 10-15, at fixed size) and fall at the highest doses. Cu a\* falls from the first doses, even where colonies are not smaller than control (area ratio 1.0-1.2 at doses 5-15). Fe rises steadily.
3. **Top-dose a\* in Cr, Cu and Pb looks like a floor, and small colonies are a concern in Pb and Cr.** At the top dose a\* collapses to a narrow band (about 5-6) for every strain, whatever its baseline. In Pb the smallest colonies (below about 1,850 px) have a\* of 6-7 with SD 2-2.7 even at dose 0. Dropping colonies below 1,000-5,000 px (section 10) leaves Cu and Fe unchanged. The Cr top-dose drop persists in larger colonies (-9.1 with no filter, -10.7 at 5,000 px). The Pb rise at doses 5-15 persists (+4.9 to +9.5 at 2,000 px), but almost no Pb colonies at doses 20-30 pass the larger thresholds (242-265 wells at 2,000 px, 8-17 at 5,000 px), so the Pb inhibitory range cannot be tested.
4. **Genetic structure matters more than species labels.** R. mucilaginosa populations differ in baseline a\* in all four metals (omnibus p < 0.001; population 5 is lowest, -2.6 to -4.0 a\* units vs population 1). Species differences are metal-specific. R. toruloides, R. dairenensis and R. diobovata have lower estimates than R. mucilaginosa in all four metals, but each is individually significant in only one to three of them.
5. **Phylogenetic signal is strong.** Pagel's lambda for baseline a\* is 0.90-0.98 in Cr, Cu and Pb, also within R. mucilaginosa alone. Fe is lower (0.51; 0.23 size-adjusted).
6. **Regression to the mean inflates the baseline-vs-induction pattern.** Selecting the high-baseline third of strains on one half of the replicate wells and measuring the change on the other half gives a smaller gap than the naive analysis: Cu -1.7 vs -5.3, Cr -5.5 vs -6.9, Pb -4.2 vs -6.7. Fe is unaffected (+1.9 vs +1.6).
7. **Strain-specific change in a\* is only partly repeatable.** Correlation of the change between replicate halves: Fe 0.70, Cr 0.54, Pb 0.31, Cu 0.13.
8. **Run (batch) is a small share of variance (0-2%).** Only Cr shows heterogeneity in the dose effect between runs (I2 = 0.68, Q p = 0.025): the effect weakens from -12.6 (run d000320) to -5.6 (run d000323). Each run is a different strain subset, so run and strain set cannot be separated.

## 1. Data and method (short)

- Source: DuckDB table `heavy_metal_measurement` (517,371 colony rows; 5 metals). Named strains only (502,051 rows). Unit of analysis is the **well**: the largest object per well and image, median over an absolute late window [T-24 h, T] after plate start, T = median plate imaging span per metal (Cr 114, Cu 114, Fe 89.9, Pb 107.5, Zn 107.6 h); wells need at least 2 images in the window.
- Models (lme4, REML): strain, run and plate are random effects. `dose_s` = concentration / maximum concentration. "At fixed size" adds ln(area) centred on the control mean.
- Dose units are not given in the source. Cr spans 0-1.2 and the other metals 0-30, so metals are never pooled.
- Zinc: the two runs are treated as one experiment. It has 7 usable plates, one per dose, 1 well per strain per dose, and two disjoint strain sets (doses 0-15 and 10-30). Dose and plate cannot be separated, so Zinc is shown for completeness only and is not used in the stratified models.

## 2. Overall dose response

![](figures/fig1_astar_vs_size_by_dose.png)

*Figure 1. a\* against colony size by dose (well medians, binned). Overlapping curves mean a size effect only; separated curves mean a shift at the same size.*


| Metal | model | n_wells | n_strains | dose effect (0 to top) | SE | p | strain share @0 | strain share @top | singular |
|--------|-------------|-------|---------|----------------------|----|------|---------------|-----------------|--------|
| Chromium | total | 8108 | 316 | -13.54 | 1.07 | <0.001 | 0.50 | 0.24 | False |
| Chromium | at fixed size | 8108 | 316 | -8.36 | 0.91 | <0.001 | 0.55 | 0.24 | False |
| Copper | total | 7652 | 319 | -15.86 | 0.39 | <0.001 | 0.47 | 0.36 | False |
| Copper | at fixed size | 7652 | 319 | -13.37 | 0.31 | <0.001 | 0.51 | 0.28 | False |
| Iron | total | 4073 | 220 | 2.39 | 0.36 | <0.001 | 0.80 | 0.89 | False |
| Iron | at fixed size | 4073 | 220 | 3.83 | 0.35 | <0.001 | 0.82 | 0.91 | False |
| Lead | total | 6755 | 312 | -16.77 | 1.54 | <0.001 | 0.51 | 0.00 | True |
| Lead | at fixed size | 6755 | 312 | -6.10 | 1.15 | <0.001 | 0.49 | 0.03 | False |
| Zinc | total | 505 | 151 | -16.93 | 4.49 | 0.01 | 0.26 | 0.72 | False |
| Zinc | at fixed size | 505 | 151 | -3.66 | 4.46 | 0.44 | 0.41 | 0.66 | False |

*Table 1. Dose effect on a\* from 0 to the top dose. "total" has no size term; "at fixed size" adds ln(area). "strain share" is strain variance / (strain + plate + residual) at dose 0 and at the top dose. The Pb total model is a singular fit.*

| dose | Chromium | Copper | Iron | Lead | Zinc |
|-----|--------|------|----|-----|-----|
| 0.20 | 2.17 |  |  |  |  |
| 0.40 | 3.10 |  |  |  |  |
| 0.60 | 1.18 |  |  |  |  |
| 0.80 | 0.94 |  |  |  |  |
| 1.00 | -3.31 |  |  |  |  |
| 1.20 | -9.12 |  |  |  |  |
| 5.00 |  | 0.05 | 2.53 | 4.44 | 5.43 |
| 10.00 |  | -3.21 | 3.09 | 7.49 | 6.80 |
| 15.00 |  | -5.56 | 2.87 | 7.56 | 6.91 |
| 20.00 |  | -8.61 | 3.58 | -0.90 | 5.38 |
| 25.00 |  | -10.56 | 4.28 | -0.65 | 1.37 |
| 30.00 |  | -12.30 | 4.64 | -1.41 | -1.76 |

*Table 2. a\* difference from 0 dose at fixed size (dose-as-factor model; SE 0.2-0.5 except Zn, about 8.8).*

![](figures/fig2_baseline_vs_stressed.png)

*Figure 2. Per-strain a\* (top) and ln area (bottom), unstressed against top dose (mean of replicate wells, bars = SE). For Zinc the top dose shown is 15, the highest dose shared with dose 0.*


![](figures/fig3_strain_baseline_vs_induction.png)

*Figure 2b. Model-based strain effects: baseline a\* (x) against strain-specific extra change in a\* at the top dose (y), total (top row) and at fixed size (bottom row). Pb total is a singular fit and is not drawn. Pb at fixed size (r = -0.99) and Zn total (r = +1.00) lie on a straight line: the random-effect correlation is at the model boundary, so these points are a rescaling of one random effect and carry no separate information.*


![](figures/fig4_repeatability.png)

*Figure 3. Share of a\* variance due to strain at 0 dose and at the top dose (size-adjusted model).*


**Reading.** In Cr, Cu and Pb strains converge to a common low a\* at the top dose (strain share falls to 0.24, 0.28 and 0.03). Fe keeps its strain differences (0.82 to 0.91).

## 3. Which traits correlate with a\*

![](figures/fig5_astar_trait_correlations.png)

*Figure 4. Spearman correlation of a\* with the traits most correlated with it (well level, all doses).*


All doses (excluding `ColorLab_a*Medoid`, which is a\* itself):

| trait | Chromium | Copper | Iron | Lead | Zinc | mean |rho| |
|-----------------------------|--------|------|-----|----|----|----------|
| ColorHSV_SaturationRobustMean | 0.86 | 0.50 | 0.90 | 0.94 | 0.95 | 0.83 |
| ColorLab_LabTotalVariance | 0.83 | 0.02 | 0.75 | 0.87 | 0.78 | 0.65 |
| ColorLab_b*Medoid | 0.66 | 0.38 | 0.67 | 0.60 | 0.72 | 0.61 |
| ColorLab_b*GeoMedian | 0.66 | 0.37 | 0.66 | 0.59 | 0.71 | 0.60 |
| ColorHSV_HSVConeVariance | 0.74 | 0.05 | 0.60 | 0.77 | 0.67 | 0.57 |
| Intensity_ConvexDensity | 0.58 | 0.66 | -0.31 | 0.71 | 0.57 | 0.57 |
| Intensity_IntegratedIntensity | 0.59 | 0.65 | -0.24 | 0.71 | 0.57 | 0.55 |
| ColorHSV_ValueRobustMean | 0.84 | 0.83 | 0.16 | 0.08 | 0.81 | 0.55 |
| Shape_MinFeretDiameter | 0.62 | 0.63 | -0.12 | 0.72 | 0.59 | 0.53 |
| Shape_MinorAxisLength | 0.62 | 0.62 | -0.12 | 0.72 | 0.59 | 0.53 |

Control dose only:

| trait | Chromium | Copper | Iron | Lead | Zinc | mean |rho| |
|------------------------------------|--------|------|-----|-----|-----|----------|
| ColorHSV_SaturationRobustMean | 0.85 | 0.63 | 0.90 | 0.90 | 0.91 | 0.84 |
| ColorLab_LabTotalVariance | 0.79 | 0.74 | 0.85 | 0.81 | 0.93 | 0.83 |
| ColorHSV_HSVConeVariance | 0.68 | 0.62 | 0.79 | 0.74 | 0.84 | 0.73 |
| ColorHSV_ValueRobustMean | 0.69 | 0.55 | 0.57 | 0.70 | 0.47 | 0.60 |
| ColorLab_b*Medoid | 0.40 | 0.31 | 0.41 | 0.48 | 0.44 | 0.41 |
| Intensity_StandardDeviationIntensity | -0.29 | -0.56 | -0.38 | -0.42 | -0.37 | 0.40 |
| ColorLab_b*GeoMedian | 0.39 | 0.30 | 0.40 | 0.47 | 0.42 | 0.40 |
| Intensity_MedianIntensity | -0.36 | -0.42 | -0.38 | -0.20 | -0.55 | 0.38 |
| ColorLab_DeltaE2000MeanFromMedoid | 0.31 | 0.25 | 0.51 | 0.26 | 0.55 | 0.38 |
| Texture_Correlation-deg090-scale05 | -0.24 | -0.46 | -0.36 | -0.30 | -0.43 | 0.36 |

Within strain and dose (replicate wells only; Zinc has no replicates):

| trait | Chromium | Copper | Iron | Lead | mean |rho| |
|-----------------------------|--------|------|----|----|----------|
| ColorHSV_SaturationRobustMean | 0.83 | 0.63 | 0.88 | 0.84 | 0.79 |
| ColorHSV_ValueRobustMean | 0.75 | 0.50 | 0.63 | 0.68 | 0.64 |
| ColorLab_b*Medoid | 0.67 | 0.45 | 0.73 | 0.60 | 0.61 |
| ColorLab_b*GeoMedian | 0.67 | 0.44 | 0.74 | 0.60 | 0.61 |
| ColorLab_LabTotalVariance | 0.67 | 0.27 | 0.55 | 0.69 | 0.55 |
| ColorHSV_HSVConeVariance | 0.58 | 0.21 | 0.44 | 0.63 | 0.47 |
| Shape_MinorAxisLength | 0.35 | 0.44 | 0.54 | 0.50 | 0.46 |
| Shape_Solidity | 0.35 | 0.44 | 0.54 | 0.50 | 0.46 |
| Shape_MinFeretDiameter | 0.34 | 0.44 | 0.54 | 0.50 | 0.46 |
| Shape_MedianRadius | 0.32 | 0.44 | 0.53 | 0.50 | 0.45 |

**Reading.** Saturation, b\* and the colour-variance traits correlate most strongly with a\*; these are partly the same colour information. Size traits (area, radius, Feret diameter, integrated intensity) correlate at 0.5-0.7 in Cr, Cu, Pb and Zn but with the opposite sign in Fe. Many size traits are near-duplicates of each other, so treat them as one family. Replicate wells within a strain and dose still show a size-a\* correlation of about 0.3-0.54, so the link is not only between strains.

## 4. Species and population (stratification item 1)

Species is a fixed effect (reference R. mucilaginosa); strain, run and plate are random. Only species with at least 5 strains are included. Population labels exist for 201 R. mucilaginosa strains (file `analysis/gwas/data/prior_run_state/pop_assignment_at_run.csv`).

| Metal | model | F | p | p_BH |
|--------|--------------------------------------|-----|------|------|
| Chromium | baseline | 7.74 | <0.001 | <0.001 |
| Chromium | baseline_size_adjusted | 6.40 | <0.001 | <0.001 |
| Chromium | dose_response_size_adjusted | 9.14 | <0.001 | <0.001 |
| Chromium | population_baseline_size_adjusted | 8.46 | <0.001 | <0.001 |
| Chromium | population_dose_response_size_adjusted | 14.18 | <0.001 | <0.001 |
| Copper | baseline | 1.66 | 0.12 | 0.13 |
| Copper | baseline_size_adjusted | 1.48 | 0.17 | 0.18 |
| Copper | dose_response_size_adjusted | 0.64 | 0.72 | 0.72 |
| Copper | population_baseline_size_adjusted | 5.46 | <0.001 | <0.001 |
| Copper | population_dose_response_size_adjusted | 6.56 | <0.001 | <0.001 |
| Iron | baseline | 4.40 | <0.001 | <0.001 |
| Iron | baseline_size_adjusted | 2.97 | 0.01 | 0.01 |
| Iron | dose_response_size_adjusted | 1.96 | 0.07 | 0.09 |
| Iron | population_baseline_size_adjusted | 18.96 | <0.001 | <0.001 |
| Iron | population_dose_response_size_adjusted | 7.68 | <0.001 | <0.001 |
| Lead | baseline | 15.76 | <0.001 | <0.001 |
| Lead | baseline_size_adjusted | 4.33 | <0.001 | <0.001 |
| Lead | dose_response_size_adjusted | 8.75 | <0.001 | <0.001 |
| Lead | population_baseline_size_adjusted | 12.00 | <0.001 | <0.001 |
| Lead | population_dose_response_size_adjusted | 13.73 | <0.001 | <0.001 |

*Table 3. Omnibus tests (F test; Benjamini-Hochberg across all rows).*

![](figures/s1_species_baseline_forest.png)

*Figure 5. Species effect on baseline a\* relative to R. mucilaginosa (size-adjusted).*


| Metal | species | n_strains | estimate | se | p | p_BH |
|--------|---------------|---------|--------|----|------|------|
| Chromium | R. dairenensis | 9 | -2.48 | 1.23 | 0.05 | 0.11 |
| Chromium | R. diobovata | 8 | -2.22 | 1.23 | 0.07 | 0.14 |
| Chromium | R. graminis | 5 | -4.76 | 1.99 | 0.02 | 0.05 |
| Chromium | R. paludigena | 15 | -2.45 | 1.03 | 0.02 | 0.05 |
| Chromium | R. sp_clade_I | 9 | 3.56 | 1.26 | 0.01 | 0.03 |
| Chromium | R. sphaerocarpa | 6 | -0.04 | 1.41 | 0.98 | 0.98 |
| Chromium | R. taiwanensis | 6 | 3.42 | 1.41 | 0.02 | 0.05 |
| Chromium | R. toruloides | 10 | -4.63 | 1.10 | <0.001 | <0.001 |
| Copper | R. dairenensis | 10 | -1.09 | 1.02 | 0.28 | 0.40 |
| Copper | R. diobovata | 8 | -1.50 | 1.12 | 0.18 | 0.30 |
| Copper | R. paludigena | 18 | -1.04 | 0.82 | 0.21 | 0.31 |
| Copper | R. sp_clade_I | 9 | -0.60 | 1.10 | 0.59 | 0.71 |
| Copper | R. sphaerocarpa | 7 | -1.22 | 1.23 | 0.32 | 0.41 |
| Copper | R. taiwanensis | 6 | 0.10 | 1.33 | 0.94 | 0.98 |
| Copper | R. toruloides | 10 | -2.60 | 1.05 | 0.01 | 0.05 |
| Iron | R. dairenensis | 8 | -2.83 | 1.16 | 0.02 | 0.05 |
| Iron | R. diobovata | 5 | -0.73 | 1.41 | 0.61 | 0.71 |
| Iron | R. paludigena | 5 | 1.77 | 1.42 | 0.21 | 0.31 |
| Iron | R. sp_clade_I | 5 | -2.85 | 1.44 | 0.05 | 0.11 |
| Iron | R. taiwanensis | 6 | 2.33 | 1.30 | 0.07 | 0.14 |
| Iron | R. toruloides | 6 | -1.94 | 1.30 | 0.14 | 0.24 |
| Lead | R. dairenensis | 10 | -3.67 | 1.05 | <0.001 | 0.01 |
| Lead | R. diobovata | 8 | -3.53 | 1.13 | 0.00 | 0.02 |
| Lead | R. paludigena | 12 | 0.04 | 1.17 | 0.97 | 0.98 |
| Lead | R. sp_clade_I | 9 | -2.39 | 1.26 | 0.06 | 0.13 |
| Lead | R. sphaerocarpa | 7 | -3.87 | 1.38 | 0.01 | 0.03 |
| Lead | R. taiwanensis | 6 | -0.14 | 1.28 | 0.91 | 0.98 |
| Lead | R. toruloides | 10 | -1.06 | 1.06 | 0.32 | 0.41 |

*Table 4. Species baseline contrasts (size-adjusted).*

![](figures/s1_species_slope_forest.png)

*Figure 6. Species-specific change in a\* over the full dose range at fixed size.*


| Metal | species | n_strains | slope_full_range | se | p vs mucilaginosa |
|--------|---------------|---------|----------------|----|-----------------|
| Chromium | R. mucilaginosa | 216 | -8.49 | 0.90 |  |
| Chromium | R. dairenensis | 9 | -6.74 | 1.37 | 0.10 |
| Chromium | R. diobovata | 8 | -5.19 | 1.38 | 0.00 |
| Chromium | R. graminis | 5 | -4.64 | 1.91 | 0.03 |
| Chromium | R. paludigena | 15 | -5.34 | 1.24 | <0.001 |
| Chromium | R. sp_clade_I | 9 | -13.77 | 1.39 | <0.001 |
| Chromium | R. sphaerocarpa | 6 | -7.70 | 1.51 | 0.53 |
| Chromium | R. taiwanensis | 6 | -11.02 | 1.52 | 0.04 |
| Chromium | R. toruloides | 10 | -4.83 | 1.29 | <0.001 |
| Copper | R. mucilaginosa | 215 | -13.77 | 0.30 |  |
| Copper | R. dairenensis | 10 | -13.99 | 0.91 | 0.80 |
| Copper | R. diobovata | 8 | -12.83 | 0.96 | 0.32 |
| Copper | R. paludigena | 18 | -14.10 | 0.72 | 0.63 |
| Copper | R. sp_clade_I | 9 | -13.82 | 0.96 | 0.95 |
| Copper | R. sphaerocarpa | 7 | -12.38 | 1.12 | 0.21 |
| Copper | R. taiwanensis | 6 | -13.10 | 1.13 | 0.55 |
| Copper | R. toruloides | 10 | -12.73 | 0.89 | 0.24 |
| Iron | R. mucilaginosa | 174 | 3.77 | 0.37 |  |
| Iron | R. dairenensis | 8 | 5.30 | 1.36 | 0.25 |
| Iron | R. diobovata | 5 | 1.52 | 1.35 | 0.10 |
| Iron | R. paludigena | 5 | 6.77 | 1.46 | 0.04 |
| Iron | R. sp_clade_I | 5 | 4.49 | 1.39 | 0.61 |
| Iron | R. taiwanensis | 6 | 5.59 | 1.31 | 0.17 |
| Iron | R. toruloides | 6 | 5.25 | 1.24 | 0.23 |
| Lead | R. mucilaginosa | 214 | -6.83 | 1.14 |  |
| Lead | R. dairenensis | 10 | -2.26 | 1.75 | 0.00 |
| Lead | R. diobovata | 8 | -0.11 | 1.76 | <0.001 |
| Lead | R. paludigena | 12 | -4.63 | 1.72 | 0.11 |
| Lead | R. sp_clade_I | 9 | -2.90 | 1.87 | 0.01 |
| Lead | R. sphaerocarpa | 7 | 0.93 | 1.96 | <0.001 |
| Lead | R. taiwanensis | 6 | -5.71 | 1.96 | 0.49 |
| Lead | R. toruloides | 10 | -3.02 | 1.69 | 0.00 |

*Table 5. Species-specific dose slopes (a\* change over the full dose range).*

| Metal | n_other_strains | baseline_diff | baseline_se | baseline_p | slope_diff | slope_se | slope_p |
|--------|---------------|-------------|-----------|----------|----------|--------|-------|
| Chromium | 68 | -1.38 | 0.53 | 0.01 | 1.29 | 0.48 | 0.01 |
| Copper | 68 | -1.19 | 0.45 | 0.01 | 0.34 | 0.39 | 0.38 |
| Iron | 35 | -0.84 | 0.60 | 0.17 | 1.02 | 0.60 | 0.09 |
| Lead | 62 | -2.09 | 0.51 | 0.00 | 4.27 | 0.63 | 0.00 |

*Table 6. R. mucilaginosa against all other species pooled (contrast only).*

![](figures/s1_population.png)

*Figure 7. Populations within R. mucilaginosa: baseline a\* (top) and dose response (bottom).*


| Metal | population | n_strains | estimate | se | p |
|--------|----------|---------|--------|----|------|
| Chromium | pop2 | 28 | 0.42 | 0.67 | 0.53 |
| Chromium | pop3 | 25 | -2.05 | 0.67 | 0.00 |
| Chromium | pop4 | 37 | 0.45 | 0.55 | 0.42 |
| Chromium | pop5 | 16 | -3.39 | 0.77 | <0.001 |
| Chromium | pop6 | 20 | -2.17 | 0.68 | 0.00 |
| Copper | pop2 | 28 | 0.58 | 0.62 | 0.35 |
| Copper | pop3 | 25 | -1.03 | 0.63 | 0.10 |
| Copper | pop4 | 36 | 1.11 | 0.52 | 0.04 |
| Copper | pop5 | 16 | -2.58 | 0.75 | <0.001 |
| Copper | pop6 | 20 | -0.69 | 0.66 | 0.29 |
| Iron | pop2 | 15 | 0.16 | 0.65 | 0.81 |
| Iron | pop3 | 16 | -1.97 | 0.63 | 0.00 |
| Iron | pop4 | 32 | 1.52 | 0.46 | 0.00 |
| Iron | pop5 | 15 | -3.96 | 0.63 | <0.001 |
| Iron | pop6 | 20 | -2.57 | 0.54 | <0.001 |
| Lead | pop2 | 28 | 0.79 | 0.58 | 0.18 |
| Lead | pop3 | 24 | -1.16 | 0.62 | 0.06 |
| Lead | pop4 | 36 | 0.86 | 0.52 | 0.10 |
| Lead | pop5 | 16 | -3.83 | 0.72 | <0.001 |
| Lead | pop6 | 20 | -2.51 | 0.64 | <0.001 |

*Table 7. Population baseline contrasts against pop1 (size-adjusted).*

| Metal | population | n_strains | slope_full_range | se | z_p | vs_pop1_p |
|--------|----------|---------|----------------|----|----|---------|
| Chromium | pop1 | 75 | -9.18 | 0.96 | 0.00 |  |
| Chromium | pop2 | 28 | -9.61 | 1.03 | 0.00 | 0.48 |
| Chromium | pop3 | 25 | -6.39 | 1.03 | 0.00 | 0.00 |
| Chromium | pop4 | 37 | -10.07 | 1.00 | 0.00 | 0.07 |
| Chromium | pop5 | 16 | -6.06 | 1.09 | 0.00 | 0.00 |
| Chromium | pop6 | 20 | -6.82 | 1.06 | 0.00 | 0.00 |
| Copper | pop1 | 75 | -13.74 | 0.38 | 0.00 |  |
| Copper | pop2 | 28 | -15.11 | 0.53 | 0.00 | 0.02 |
| Copper | pop3 | 25 | -12.88 | 0.56 | 0.00 | 0.16 |
| Copper | pop4 | 36 | -15.29 | 0.48 | 0.00 | 0.00 |
| Copper | pop5 | 16 | -12.08 | 0.70 | 0.00 | 0.02 |
| Copper | pop6 | 20 | -13.14 | 0.62 | 0.00 | 0.34 |
| Iron | pop1 | 69 | 2.82 | 0.42 | 0.00 |  |
| Iron | pop2 | 15 | 2.41 | 0.92 | 0.01 | 0.65 |
| Iron | pop3 | 16 | 3.34 | 0.86 | 0.00 | 0.54 |
| Iron | pop4 | 32 | 6.01 | 0.53 | 0.00 | 0.00 |
| Iron | pop5 | 15 | 4.46 | 0.83 | 0.00 | 0.05 |
| Iron | pop6 | 20 | 3.58 | 0.62 | 0.00 | 0.23 |
| Lead | pop1 | 75 | -10.50 | 1.34 | 0.00 |  |
| Lead | pop2 | 28 | -10.56 | 1.43 | 0.00 | 0.94 |
| Lead | pop3 | 24 | -9.18 | 1.46 | 0.00 | 0.12 |
| Lead | pop4 | 36 | -12.24 | 1.39 | 0.00 | 0.01 |
| Lead | pop5 | 16 | -4.62 | 1.53 | 0.00 | 0.00 |
| Lead | pop6 | 20 | -8.36 | 1.48 | 0.00 | 0.01 |

*Table 8. Population dose slopes.*

![](figures/fig6_baseline_astar_by_species.png)

*Figure 8. Baseline a\* by species (strain means; species with at least 5 strains).*


**Reading.**
- The population effect is the most consistent result in this section. Population 5 is lowest in all four metals (Cr -3.4, Cu -2.6, Fe -4.0, Pb -3.8). Population 6 is also low in Cr, Fe and Pb, and population 3 in Cr and Fe. Population 4 has the highest or tied-highest estimate in every metal, but differs from pop1 significantly only in Cu and Fe. The same pattern in four independent metal screens suggests a real genetic effect.
- Species effects are inconsistent across metals (for example sp_clade_I is +3.6 in Cr and -2.4 in Pb), and the non-mucilaginosa species have 5-18 strains each, so single-species contrasts are noisy. Cu shows no species effect overall (F p = 0.17).
- Pooled "other species" are lower than R. mucilaginosa in baseline a\* in all four metals (-0.8 to -2.1), and respond differently in Pb (slope difference +4.3).

## 5. Phylogenetic signal (item 2)

![](figures/s2_phylo_signal.png)

*Figure 9. Pagel's lambda of baseline a\* and of the change in a\* at the top dose.*


| Metal | scope | trait | n_strains | pagel_lambda | p (lambda > 0) |
|--------|--------------------|---------------------------|---------|------------|--------------|
| Chromium | all_strains_with_tip | baseline a* | 256 | 0.97 | <0.001 |
| Chromium | all_strains_with_tip | baseline a* (size-adjusted) | 256 | 0.97 | <0.001 |
| Chromium | all_strains_with_tip | change in a* at dose 1.2 | 254 | 0.96 | <0.001 |
| Chromium | R_mucilaginosa_only | baseline a* | 201 | 0.97 | <0.001 |
| Chromium | R_mucilaginosa_only | baseline a* (size-adjusted) | 201 | 0.97 | <0.001 |
| Chromium | R_mucilaginosa_only | change in a* at dose 1.2 | 199 | 0.96 | <0.001 |
| Copper | all_strains_with_tip | baseline a* | 264 | 0.93 | <0.001 |
| Copper | all_strains_with_tip | baseline a* (size-adjusted) | 264 | 0.95 | <0.001 |
| Copper | all_strains_with_tip | change in a* at dose 30 | 259 | 0.68 | <0.001 |
| Copper | R_mucilaginosa_only | baseline a* | 201 | 0.93 | <0.001 |
| Copper | R_mucilaginosa_only | baseline a* (size-adjusted) | 201 | 0.95 | <0.001 |
| Copper | R_mucilaginosa_only | change in a* at dose 30 | 201 | 0.67 | <0.001 |
| Iron | all_strains_with_tip | baseline a* | 197 | 0.51 | <0.001 |
| Iron | all_strains_with_tip | baseline a* (size-adjusted) | 197 | 0.23 | 0.00 |
| Iron | all_strains_with_tip | change in a* at dose 30 | 133 | 0.93 | <0.001 |
| Iron | R_mucilaginosa_only | baseline a* | 167 | 0.52 | <0.001 |
| Iron | R_mucilaginosa_only | baseline a* (size-adjusted) | 167 | 0.26 | 0.01 |
| Iron | R_mucilaginosa_only | change in a* at dose 30 | 113 | 0.93 | <0.001 |
| Lead | all_strains_with_tip | baseline a* | 254 | 0.90 | <0.001 |
| Lead | all_strains_with_tip | baseline a* (size-adjusted) | 254 | 0.90 | <0.001 |
| Lead | all_strains_with_tip | change in a* at dose 30 | 252 | 0.86 | <0.001 |
| Lead | R_mucilaginosa_only | baseline a* | 199 | 0.90 | <0.001 |
| Lead | R_mucilaginosa_only | baseline a* (size-adjusted) | 199 | 0.90 | <0.001 |
| Lead | R_mucilaginosa_only | change in a* at dose 30 | 198 | 0.86 | <0.001 |
| Zinc | all_strains_with_tip | baseline a* | 75 | 0.02 | 0.85 |
| Zinc | all_strains_with_tip | baseline a* (size-adjusted) | 75 | 0.51 | 0.04 |
| Zinc | all_strains_with_tip | change in a* at dose 10 | 75 | 0.95 | <0.001 |
| Zinc | R_mucilaginosa_only | baseline a* | 65 | 0.02 | 0.89 |
| Zinc | R_mucilaginosa_only | baseline a* (size-adjusted) | 65 | 0.60 | 0.02 |
| Zinc | R_mucilaginosa_only | change in a* at dose 10 | 65 | 0.95 | <0.001 |

*Table 9. Pagel's lambda (maximum likelihood under Brownian motion on the PHYling FastTree). 254-264 strains with a tree tip in Cr, Cu and Pb; 197 in Fe; 75 in Zn.*

![](figures/s2_distance_decay.png)

*Figure 10. Mean absolute difference in size-adjusted baseline a\* between pairs of strains, by patristic distance.*


![](figures/s2_tree_with_astar.png)

*Figure 11. The tree with baseline a\* per metal.*


**Reading and limits.**
- Closely related strains share baseline a\*. Lambda is 0.90-0.98 for baseline a\* in Cr, Cu and Pb, also inside R. mucilaginosa alone (0.90-0.98), and the pairwise difference in a\* increases with patristic distance (Spearman 0.35 Cr, 0.17 Cu, 0.21 Pb). Fe is weaker (0.51 raw, 0.23 size-adjusted). Baseline colony size has lambda near 0 in Cr, Cu and Pb but 0.83 in Fe and 0.99 in Zn.
- **Blomberg's K was dropped.** The tree has 22 zero-length tips and 74 tips with a neighbour closer than 1e-5, so the covariance matrix is nearly singular (condition number about 1e10). K varied from 1e-7 to 4e-2 with the diagonal jitter (Table 10), so it is not interpretable. Lambda was stable across jitter 1e-8 to 1e-2 (0.89 to 0.99 for Cr, Cu and Pb).
- Lambda near 1 partly reflects near-identical (clonal) strains having similar a\*. It does not show that a\* evolves under Brownian motion.
- 266 of 321 strains have a unique tree tip; strains without a tip are excluded. These are mostly outside the sequenced set.

| Metal | jitter | cond | lam | K | pK |
|--------|------|--------|-----|--------|-----|
| Copper | 1e-08 | 1.17e+10 | 0.95 | 4.88e-07 | 0.145 |
| Copper | 1e-06 | 1.17e+08 | 0.955 | 1.89e-05 | 0.19 |
| Copper | 0.0001 | 1.17e+06 | 0.999 | 0.000772 | 0.06 |
| Copper | 0.001 | 1.17e+05 | 0.999 | 0.00481 | 0.005 |
| Copper | 0.01 | 1.17e+04 | 0.999 | 0.0309 | 0.005 |
| Lead | 1e-08 | 9.93e+09 | 0.894 | 9.91e-07 | 0.03 |
| Lead | 1e-06 | 9.93e+07 | 0.901 | 3.64e-05 | 0.06 |
| Lead | 0.0001 | 9.93e+05 | 0.998 | 0.00121 | 0.01 |
| Lead | 0.001 | 9.93e+04 | 0.998 | 0.0066 | 0.005 |
| Lead | 0.01 | 9.93e+03 | 0.991 | 0.042 | 0.005 |
| Chromium | 1e-08 | 1.03e+10 | 0.965 | 4.1e-07 | 0.23 |
| Chromium | 1e-06 | 1.03e+08 | 0.969 | 2.06e-05 | 0.18 |
| Chromium | 0.0001 | 1.03e+06 | 0.991 | 0.000855 | 0.045 |
| Chromium | 0.001 | 1.03e+05 | 0.975 | 0.00554 | 0.005 |
| Chromium | 0.01 | 1.03e+04 | 0.985 | 0.0404 | 0.005 |

*Table 10. Sensitivity of lambda and K to the diagonal jitter added to the covariance matrix.*

## 6. Dose regimes (item 3)

Regimes are defined from colony size: a dose is **sub-inhibitory** when the median (across strains) of area at that dose / area at dose 0 is at least 0.5, and **inhibitory** below that. The 0.5 cutoff is a choice, not a measured threshold.

![](figures/s3_regimes.png)

*Figure 12. Top: area ratio by dose with the 0.5 cutoff. Bottom: a\* difference from 0 dose, with and without size.*


| Metal | conc | median_ratio | q25 | q75 | n_strains | regime |
|--------|-----|------------|----|----|---------|--------------|
| Chromium | 0.20 | 1.01 | 0.98 | 1.04 | 306 | sub-inhibitory |
| Chromium | 0.40 | 0.93 | 0.89 | 0.96 | 306 | sub-inhibitory |
| Chromium | 0.60 | 0.40 | 0.33 | 0.49 | 305 | inhibitory |
| Chromium | 0.80 | 0.37 | 0.29 | 0.46 | 305 | inhibitory |
| Chromium | 1.00 | 0.09 | 0.06 | 0.13 | 301 | inhibitory |
| Chromium | 1.20 | 0.04 | 0.02 | 0.06 | 305 | inhibitory |
| Copper | 5.00 | 1.24 | 1.14 | 1.35 | 319 | sub-inhibitory |
| Copper | 10.00 | 1.13 | 1.03 | 1.23 | 318 | sub-inhibitory |
| Copper | 15.00 | 1.04 | 0.92 | 1.16 | 318 | sub-inhibitory |
| Copper | 20.00 | 0.89 | 0.77 | 1.04 | 317 | sub-inhibitory |
| Copper | 25.00 | 0.76 | 0.58 | 0.97 | 315 | sub-inhibitory |
| Copper | 30.00 | 0.46 | 0.29 | 0.71 | 314 | inhibitory |
| Iron | 5.00 | 1.04 | 0.97 | 1.11 | 219 | sub-inhibitory |
| Iron | 10.00 | 0.87 | 0.80 | 0.93 | 218 | sub-inhibitory |
| Iron | 15.00 | 0.82 | 0.75 | 0.89 | 218 | sub-inhibitory |
| Iron | 20.00 | 0.68 | 0.61 | 0.75 | 218 | sub-inhibitory |
| Iron | 25.00 | 0.70 | 0.62 | 0.77 | 147 | sub-inhibitory |
| Iron | 30.00 | 0.66 | 0.60 | 0.73 | 146 | sub-inhibitory |
| Lead | 5.00 | 0.84 | 0.76 | 0.91 | 303 | sub-inhibitory |
| Lead | 10.00 | 0.83 | 0.77 | 0.90 | 302 | sub-inhibitory |
| Lead | 15.00 | 0.36 | 0.30 | 0.46 | 304 | inhibitory |
| Lead | 20.00 | 0.04 | 0.03 | 0.06 | 301 | inhibitory |
| Lead | 25.00 | 0.04 | 0.03 | 0.06 | 302 | inhibitory |
| Lead | 30.00 | 0.05 | 0.03 | 0.06 | 302 | inhibitory |

*Table 11. Colony area relative to 0 dose.*

![](figures/s3_regime_slopes.png)

*Figure 13. Slope of a\* per 10% of the maximum dose, by regime.*


| Metal | regime | doses | n_wells | size_adjusted | slope_per_0.1_of_max_dose | se | p |
|--------|--------------|------------------|-------|-------------|-------------------------|----|------|
| Chromium | sub-inhibitory | 0,0.2,0.4 | 3577 | False | 0.89 | 0.11 | <0.001 |
| Chromium | sub-inhibitory | 0,0.2,0.4 | 3577 | True | 0.92 | 0.11 | <0.001 |
| Chromium | inhibitory | 0,0.6,0.8,1,1.2 | 5728 | False | -1.23 | 0.13 | <0.001 |
| Chromium | inhibitory | 0,0.6,0.8,1,1.2 | 5728 | True | -0.74 | 0.11 | <0.001 |
| Copper | sub-inhibitory | 0,5,10,15,20,25 | 6614 | False | -1.53 | 0.05 | <0.001 |
| Copper | sub-inhibitory | 0,5,10,15,20,25 | 6614 | True | -1.37 | 0.03 | <0.001 |
| Iron | sub-inhibitory | 0,5,10,15,20,25,30 | 4073 | False | 0.25 | 0.03 | <0.001 |
| Iron | sub-inhibitory | 0,5,10,15,20,25,30 | 4073 | True | 0.39 | 0.03 | <0.001 |
| Lead | sub-inhibitory | 0,5,10 | 3061 | False | 2.07 | 0.22 | <0.001 |
| Lead | sub-inhibitory | 0,5,10 | 3061 | True | 2.34 | 0.09 | <0.001 |
| Lead | inhibitory | 0,15,20,25,30 | 4717 | False | -1.42 | 0.16 | <0.001 |
| Lead | inhibitory | 0,15,20,25,30 | 4717 | True | -0.59 | 0.13 | <0.001 |

*Table 12. Regime-specific slopes (a\* change per 10% of the metal's maximum dose; random effects for strain, run and plate).*

**Reading.**
- Cr: sub-inhibitory doses (0-0.4) raise a\* by +0.9 per 10% of the range; inhibitory doses lower it (-0.7 at fixed size, -1.2 without size).
- Pb: sub-inhibitory doses (0-10) raise a\* by +2.1 to +2.3 per 10%; inhibitory doses (15-30) lower it (-0.6 at fixed size, -1.4 without size). Pb colonies at doses 20-30 are about 4% of control area, so the inhibitory-regime estimate rests on very small colonies.
- Cu: no inhibitory regime except the top dose, yet a\* falls steadily (-1.4 per 10% at fixed size): Cu lowers a\* without reducing colony size.
- Fe: area never falls below 0.66 of control, and a\* rises slightly (+0.4 per 10% at fixed size).

## 7. Size-matched comparison (item 4)

Colonies are placed in five size bins (quantiles of ln area within each metal) and the dose effect is estimated inside each bin.

![](figures/s4_size_matched.png)

*Figure 14. Wells per size bin and dose (top), a\* against dose within bins (middle), and a\* spread by size bin (bottom).*


| Metal | size_bin | n_wells | n_doses | mean_a | sd_a | slope_full_range | se | p |
|--------|--------|-------|-------|------|----|----------------|----|------|
| Chromium | S1 | 1622 | 7 | 6.66 | 4.45 | -5.07 | 2.03 | 0.01 |
| Chromium | S2 | 1621 | 7 | 16.36 | 5.16 | -10.43 | 1.54 | <0.001 |
| Chromium | S3 | 1622 | 7 | 19.38 | 4.89 | -5.35 | 1.29 | <0.001 |
| Chromium | S4 | 1621 | 7 | 20.92 | 4.85 | -8.40 | 1.46 | <0.001 |
| Chromium | S5 | 1622 | 7 | 20.50 | 4.72 | -12.07 | 2.01 | <0.001 |
| Copper | S1 | 1531 | 7 | 5.10 | 6.15 | -10.16 | 0.62 | <0.001 |
| Copper | S2 | 1530 | 7 | 10.81 | 6.09 | -13.52 | 0.35 | <0.001 |
| Copper | S3 | 1530 | 7 | 13.90 | 5.74 | -14.22 | 0.35 | <0.001 |
| Copper | S4 | 1530 | 7 | 15.81 | 5.39 | -14.10 | 0.44 | <0.001 |
| Copper | S5 | 1531 | 7 | 17.78 | 4.81 | -13.62 | 0.58 | <0.001 |
| Iron | S1 | 815 | 7 | 20.71 | 5.34 | 3.53 | 0.68 | <0.001 |
| Iron | S2 | 814 | 7 | 21.24 | 4.17 | 4.52 | 0.37 | <0.001 |
| Iron | S3 | 815 | 7 | 20.62 | 3.89 | 4.20 | 0.39 | <0.001 |
| Iron | S4 | 814 | 7 | 20.17 | 3.68 | 4.03 | 0.40 | <0.001 |
| Iron | S5 | 815 | 7 | 19.91 | 3.90 | 3.14 | 0.49 | <0.001 |
| Lead | S1 | 1351 | 7 | 7.06 | 2.61 | 4.71 | 0.53 | <0.001 |
| Lead | S2 | 1351 | 7 | 6.25 | 2.58 | 4.56 | 0.90 | <0.001 |
| Lead | S3 | 1351 | 7 | 19.15 | 8.37 | -14.52 | 1.68 | <0.001 |
| Lead | S4 | 1351 | 4 | 23.36 | 6.01 | 17.02 | 1.35 | <0.001 |
| Lead | S5 | 1351 | 4 | 23.63 | 5.73 | 19.75 | 1.23 | <0.001 |

*Table 13. Dose effect within size bins (a\* change over the full dose range; random effects for strain, run, plate).*

**Reading.**
- Cu: a\* falls with dose in every size bin (-10 to -14), and all seven doses are present in every bin. This is the cleanest size-matched result: Cu lowers a\* at the same size.
- Fe: a\* rises in every bin (+3.1 to +4.5).
- Cr: negative in all bins (-5 to -12).
- Pb: size and dose are strongly confounded. The two smallest bins (below about 1,850 px) have a\* of 6.3-7.1 and SD 2.6 at all doses, and the two largest bins contain only doses 0-15. Slopes in the Pb bins differ in sign (+4.7 in the small bins, -14.5 in S3, +17 to +20 in the large bins), so there is no single size-matched Pb estimate.
- **Coverage is uneven.** High doses populate the small bins; in Cr and Pb some size-by-dose cells contain fewer than 10 wells.

## 8. Selection on baseline, split-half (item 5)

Strains with at least 2 replicate wells at dose 0 and at the top dose. Replicate wells are split at random into halves A and B (1,000 splits). **Naive:** select the top and bottom baseline thirds from all wells and measure the change in the same wells. **Split-half:** select on half A, measure the change in half B.

![](figures/s5_split_half_summary.png)

*Figure 15. Gap in the change in a\* between the high- and low-baseline thirds (left) and reliability of the strain-specific change (right).*


![](figures/s5_split_half_scatter.png)

*Figure 16. Baseline against change in a\*: same wells (top) and independent halves (bottom).*


| Metal | top_dose | n_strains | naive_gap | naive_ci_lo | naive_ci_hi | split_gap_mean | split_gap_ci_lo | split_gap_ci_hi | reliability_change_A_vs_B | reliability_baseline_A_vs_B |
|--------|--------|---------|---------|-----------|-----------|--------------|---------------|---------------|-------------------------|---------------------------|
| Chromium | 1.20 | 289 | -6.91 | -7.96 | -5.78 | -5.48 | -6.81 | -4.17 | 0.54 | 0.72 |
| Copper | 30.00 | 303 | -5.29 | -6.16 | -4.45 | -1.68 | -3.08 | -0.27 | 0.13 | 0.25 |
| Iron | 30.00 | 141 | 1.64 | 0.38 | 2.74 | 1.89 | 0.70 | 3.07 | 0.70 | 0.86 |
| Lead | 30.00 | 273 | -6.65 | -7.42 | -5.79 | -4.21 | -5.28 | -2.87 | 0.31 | 0.57 |

*Table 14. Size-adjusted a\*. `naive_gap` and `split_gap_mean` are the change in the top-baseline third minus the bottom third (95% CI from a bootstrap over strains). `reliability_change_A_vs_B` is the Spearman correlation of the strain-specific change between the two halves.*

**Reading.**
- Part of the strong negative baseline-vs-change relationship in Cr, Cu and Pb is regression to the mean: the gap shrinks from -5.3 to -1.7 in Cu, from -6.9 to -5.5 in Cr and from -6.7 to -4.2 in Pb. Fe is unchanged (+1.6 naive, +1.9 split).
- In Cu the strain-specific change is barely repeatable (0.13), and baseline a\* in Cu replicate wells is also noisy (0.25 between halves).
- The reliabilities are for halves of 1-2 wells and understate full-data reliability.

## 9. Batch (item 6)

Each metal is analysed separately, with run and plate as random effects (all models above). Each run is a different set of about 70-84 strains across all doses, so run and strain set are confounded. Zinc is excluded.

![](figures/s6_batch.png)

*Figure 17. Dose effect at fixed size by run (left four panels) and variance components (right).*


| Metal | k_runs | pooled_effect | Q | Q_p | I2 |
|--------|------|-------------|----|----|----|
| Chromium | 4 | -8.07 | 9.35 | 0.03 | 0.68 |
| Copper | 4 | -13.63 | 5.77 | 0.12 | 0.48 |
| Iron | 3 | 3.83 | 2.76 | 0.25 | 0.28 |
| Lead | 4 | -4.55 | 1.61 | 0.66 | 0.00 |

*Table 15. Heterogeneity of the dose effect across runs (I2, Cochran Q).*

| Metal | run | n_wells | n_strains | n_doses | dose_effect | se | p |
|--------|-------|-------|---------|-------|-----------|----|------|
| Chromium | d000320 | 2254 | 84 | 7 | -12.62 | 2.01 | <0.001 |
| Chromium | d000321 | 2086 | 81 | 7 | -9.52 | 1.69 | <0.001 |
| Chromium | d000322 | 2104 | 81 | 7 | -6.24 | 1.51 | <0.001 |
| Chromium | d000323 | 1664 | 70 | 7 | -5.62 | 1.72 | 0.00 |
| Copper | d000353 | 2127 | 84 | 7 | -13.77 | 0.51 | <0.001 |
| Copper | d000354 | 1976 | 84 | 7 | -13.08 | 0.36 | <0.001 |
| Copper | d000355 | 1920 | 83 | 7 | -13.87 | 0.60 | <0.001 |
| Copper | d000356 | 1629 | 68 | 7 | -14.72 | 0.61 | <0.001 |
| Iron | d000406 | 1950 | 76 | 7 | 4.20 | 0.40 | <0.001 |
| Iron | d000407 | 1768 | 73 | 7 | 3.54 | 0.42 | <0.001 |
| Iron | d000408 | 355 | 71 | 5 | 1.58 | 1.89 | 0.46 |
| Lead | d000399 | 2017 | 82 | 7 | -6.69 | 2.32 | 0.01 |
| Lead | d000400 | 1663 | 78 | 7 | -4.76 | 2.17 | 0.04 |
| Lead | d000401 | 1605 | 83 | 7 | -4.87 | 2.44 | 0.06 |
| Lead | d000402 | 1470 | 69 | 7 | -3.07 | 1.73 | 0.09 |

*Table 16. Dose effect (a\* change over the full dose range at fixed size) per run.*

| Metal | n_runs | share_strain | share_run | share_plate | share_resid | singular |
|--------|------|------------|---------|-----------|-----------|--------|
| Chromium | 4 | 0.35 | 0.02 | 0.34 | 0.28 | False |
| Copper | 4 | 0.33 | 0.02 | 0.04 | 0.61 | False |
| Iron | 3 | 0.82 | 0.02 | 0.03 | 0.14 | False |
| Lead | 4 | 0.21 | 0.00 | 0.38 | 0.42 | True |

*Table 17. Variance components (shares of random variance).*

**Reading.** Run explains 0-2% of the random variance. Plate explains 34% (Cr) and 38% (Pb) but only 3-4% in Cu and Fe. The Cr dose effect shrinks monotonically from run d000320 to d000323 (-12.6, -9.5, -6.2, -5.6); this could be a run-order effect or a difference between strain subsets.

## 10. Minimum colony area (small-colony check)

Very small colonies may have a\* dominated by background pixels. Evidence in the data: in Pb, wells below about 1,850 px have a\* of 6-7 with SD 2-2.7 at every dose (a hinge fit of a\* on ln area finds a flat segment below 1,846 px, with 30% of Pb wells below it); in Cr, wells below about 1,100-1,800 px have a\* of 4-6 with SD under 3.6. Cu, Fe and Zn colonies are rarely that small (5th percentile of area: Cu 6,759 px, Fe 13,829 px).

To test how much this matters, well-images whose largest object is smaller than a minimum area (1,000, 2,000, 3,000 or 5,000 px) were dropped, and the dose models (strain, run and plate random effects, at fixed size) were refitted. 90%, 78%, 70% and 59% of well-images remain. The 2,000 px value is the main sensitivity setting; the others show how much the result depends on that choice. Zinc is excluded.

![](figures/s8_minarea.png)

*Figure 18. Dose effect against area threshold (top), per-dose effect at fixed size by threshold (middle), and wells retained by dose (bottom).*


| Metal | 0 | 1000 | 2000 | 3000 | 5000 |
|--------|------|------|------|------|------|
| Chromium | -8.36 | -8.41 | -8.54 | -9.14 | -9.33 |
| Copper | -13.37 | -13.38 | -13.39 | -13.31 | -13.28 |
| Iron | 3.83 | 4.35 | 4.35 | 4.36 | 4.34 |
| Lead | -6.10 | -3.69 | -4.31 | -6.03 | 0.15 |

*Table 18. Dose effect on a\* over the full dose range, at fixed size, by minimum colony area (0 = no filter).*

| dose | effect @0 | effect @1000 | effect @2000 | effect @3000 | effect @5000 | wells @0 | wells @1000 | wells @2000 | wells @3000 | wells @5000 |
|----|---------|------------|------------|------------|------------|--------|-----------|-----------|-----------|-----------|
| 0 |  |  |  |  |  | 1197 | 1189 | 1188 | 1188 | 1187 |
| 0.2 | 2.2 | 2.2 | 2.2 | 2.2 | 2.2 | 1192 | 1190 | 1188 | 1187 | 1187 |
| 0.4 | 3.1 | 3.1 | 3.1 | 3.1 | 3.1 | 1188 | 1187 | 1186 | 1185 | 1183 |
| 0.6 | 1.2 | 1.1 | 1.2 | 0.9 | 0.4 | 1171 | 1168 | 1166 | 1158 | 1140 |
| 0.8 | 0.9 | 0.9 | 0.9 | 0.7 | 0.1 | 1179 | 1174 | 1172 | 1164 | 1131 |
| 1 | -3.3 | -3.2 | -2.9 | -3.1 | -3.2 | 1129 | 1011 | 899 | 751 | 491 |
| 1.2 | -9.1 | -9.6 | -9.8 | -10.5 | -10.7 | 1052 | 753 | 488 | 296 | 98 |

*Table 19. Chromium: a\* difference from 0 dose at fixed size, by threshold, and wells retained.*

| dose | effect @0 | effect @1000 | effect @2000 | effect @3000 | effect @5000 | wells @0 | wells @1000 | wells @2000 | wells @3000 | wells @5000 |
|----|---------|------------|------------|------------|------------|--------|-----------|-----------|-----------|-----------|
| 0 |  |  |  |  |  | 1023 | 1015 | 978 | 975 | 974 |
| 5 | 4.4 | 4.8 | 4.9 | 4.9 | 4.9 | 1025 | 997 | 972 | 962 | 958 |
| 10 | 7.5 | 7.8 | 8 | 7.9 | 8 | 1013 | 994 | 979 | 975 | 971 |
| 15 | 7.6 | 8.9 | 9.5 | 9.2 | 9.3 | 1004 | 968 | 925 | 906 | 864 |
| 20 | -0.9 | 0.6 | -0 | -1.8 | -1.5 | 891 | 572 | 242 | 93 | 13 |
| 25 | -0.6 | 1.2 | -0.3 | -2.2 | -2.6 | 924 | 657 | 245 | 90 | 17 |
| 30 | -1.4 | 0.4 | -0.2 | -2.4 | -1.9 | 875 | 634 | 265 | 90 | 8 |

*Table 20. Lead: a\* difference from 0 dose at fixed size, by threshold, and wells retained.*

| Metal | min_area | top_dose | n_strains_top | top_dose_mean_a | top_dose_sd_across_strains | median_area_top_px |
|--------|--------|--------|-------------|---------------|--------------------------|------------------|
| Chromium | 0 | 1.20 | 311 | 4.43 | 2.52 | 1728.99 |
| Copper | 0 | 30.00 | 314 | 4.06 | 2.64 | 16541.21 |
| Iron | 0 | 30.00 | 146 | 21.64 | 4.38 | 19508.39 |
| Lead | 0 | 30.00 | 307 | 6.48 | 1.09 | 1429.57 |
| Chromium | 1000 | 1.20 | 284 | 4.95 | 3.04 | 2405.72 |
| Copper | 1000 | 30.00 | 309 | 4.34 | 2.38 | 16974.50 |
| Iron | 1000 | 30.00 | 146 | 21.64 | 4.38 | 19508.39 |
| Lead | 1000 | 30.00 | 291 | 6.32 | 1.17 | 1794.31 |
| Chromium | 2000 | 1.20 | 241 | 5.12 | 3.57 | 3281.38 |
| Copper | 2000 | 30.00 | 308 | 4.45 | 2.35 | 17287.13 |
| Iron | 2000 | 30.00 | 146 | 21.64 | 4.38 | 19508.39 |
| Lead | 2000 | 30.00 | 189 | 5.87 | 0.77 | 2712.98 |
| Chromium | 3000 | 1.20 | 177 | 5.53 | 4.11 | 4101.66 |
| Copper | 3000 | 30.00 | 307 | 4.52 | 2.34 | 17960.48 |
| Iron | 3000 | 30.00 | 146 | 21.64 | 4.38 | 19508.39 |
| Lead | 3000 | 30.00 | 81 | 5.67 | 0.56 | 3476.00 |
| Chromium | 5000 | 1.20 | 74 | 7.13 | 6.15 | 6516.56 |
| Copper | 5000 | 30.00 | 306 | 4.66 | 2.37 | 19181.04 |
| Iron | 5000 | 30.00 | 145 | 21.81 | 3.87 | 19522.91 |
| Lead | 5000 | 30.00 | 8 | 5.28 | 0.98 | 6610.00 |

*Table 21. Strain-mean a\* at the top dose by threshold.*

**Reading.**
- **Cu and Fe are unaffected** (Cu -13.4 to -13.3; Fe +3.8 to +4.4). Their colonies are large at every dose.
- **Cr: the top-dose drop is not a small-colony artifact.** The effect at dose 1.2 is -9.1 with no filter and -10.7 with colonies of at least 5,000 px (98 wells, 74 strains remain). The effects at doses 0.2-0.8 do not change. The overall dose effect stays between -8.4 and -9.3.
- **Pb: the rise at doses 5-15 holds or strengthens** with a threshold (+4.9, +8.0 and +9.5 at doses 5, 10 and 15 at 2,000 px, against +4.4, +7.5 and +7.6 with no filter). At doses 20-30 only 242-265 wells remain at 2,000 px, 90 at 3,000 px and 8-17 at 5,000 px (572-657 at 1,000 px), so the inhibitory range of Pb cannot be tested. In those few colonies the fixed-size difference from control is between -2.6 and +1.2. The Pb overall dose effect is unstable across thresholds (-6.1, -3.7, -4.3 [singular fit], -6.0, +0.2) because it depends on how many top-dose wells are left, and should not be quoted.
- **A floor remains in Pb.** Among Pb colonies of at least 3,000 px at dose 30 (81 strains), strain-mean a\* is still 5.7 with SD 0.56 across strains. Area thresholds alone do not tell us whether this is background or a real loss of pigment.
- Fixed-size effects extrapolate a size-a\* slope fitted mostly on larger colonies; thresholds keep the comparison within the range where colonies of similar size exist, which is why the Pb inhibitory doses lose almost all their wells.


## Limits

- Dose units are not given in the source, and metals are not pooled.
- Zinc is not modelled in the stratified analyses (one plate per dose).
- Size is partly an effect of stress; models "at fixed size" estimate a direct effect, not the whole induction. One size slope is used for all strains.
- Survivor selection: wells that never produced objects are absent, and wells with fewer than 2 window images were dropped (see `report/tables/wells_lost_by_dose.csv`).
- 16 strains still have no species and are excluded from species tests (see the data-problems report). Species tests use species with at least 5 strains; R. mucilaginosa has 174-216 strains and the others have 5-18.
- Population labels come from an earlier GWAS run (201 strains).
- The all-objects sensitivity table is built (`results/wells_allobj.csv`) but not modelled.
- Old-vs-new strain assignment for Copper is unresolved (see the data-problems report). This report uses the new data's own strain labels.
