# a* (carotenoid proxy) under metal stress: stratified analyses

Generated 2026-10-08. All numbers come from `report/tables/*.csv`. Scripts are in `scripts/` and `scripts/strat/`. Data-quality caveats are in `analysis/heavy_metal_data_problems/HEAVY_METAL_DATA_PROBLEMS.md`.

**Question.** Strains grow on plates with a heavy metal at increasing doses. Does CIELAB a\* (red-green axis, used here as a carotenoid proxy) rise because of stress? Or do strains differ inherently in a\*? Colony size confounds the answer. Every dose effect is therefore shown with and without size.

## Main dataset (used everywhere unless stated)

- **Wells.** One row per replicate well. The largest object in each well and image is kept.
- **Area cutoff: colony area of at least 2,000 px.** Well-images whose largest object is smaller are dropped (section 10 gives the evidence).
- **Time window.** Median over the last 24 h up to a metal-specific end point T, counted from the first image of each plate. T is the median plate imaging span: Cr 114 h, Cu 114 h, Fe 89.9 h, Pb 107.5 h. **Zinc uses T = 80 h**, because two Zinc plates stop imaging at 83.7 and 89.7 h (section 13).
- **Wells kept.** Wells need at least 2 images in the window. Cr 7,287 wells (307 strains), Cu 7,484 (319), Fe 4,064 (218), Pb 4,606 (307), Zn 502 (150). {{T_RETAINED_NOTE}}
- **Original unfiltered results** (no area cutoff; median-span window for Zinc) are in section 14, beside the main numbers.

## Summary

1. **There is no general "stress raises a\*".** The effect depends on the metal. Over the full dose range at fixed colony size: Cr -8.5, Cu -13.4, Pb -4.3 (singular fit), Fe +4.4 a\* units. Zn cannot be estimated reliably (one plate per dose, so dose and plate effects are confounded; section 13).
2. **The dose-response is not linear.** Dose-as-factor models at fixed size show Cr a\* rising at low doses (+2.2 and +3.1 at doses 0.2 and 0.4) and falling at the top (-9.8 at 1.2). Pb rises at doses 5-15 (+4.9, +8.0, +9.5) and is near 0 at doses 20-30 (-0.0, -0.3, -0.2). Cu falls from the first doses (-3.2 at dose 10, -12.4 at dose 30), even though colonies are not smaller than control there. Fe rises steadily (+2.5 to +4.9). Zn rises to dose 15 (+7.7) and returns to 0 at dose 30; its standard errors are not valid (section 13).
3. **a\* depends on colony area even without stress, so the main dataset drops colonies below 2,000 px.** In unstressed colonies a\* is low and rises slowly below about 1,500-2,000 px (slope under 1.7 a\* per ln unit) and then rises three to four times faster (slope 5-7) to about 20 at 33,000 px (section 10). The cutoff was chosen from the unstressed Cr and Cu curves, where stressed and unstressed colonies share the same a\* below about 2,000 px. It is assumed for Pb, Fe and Zn. It removes mostly top-dose wells, so **top-dose estimates describe larger survivors** (Cr top-dose strains 311 to 241; median top-dose area 1,729 to 3,281 px). The cutoff leaves Cu unchanged (-13.4). Cr is nearly unchanged (-8.4 to -8.5). Fe moves from +3.8 to +4.4 (and its size slope from 2.9 to 3.9) when only 9 wells are removed, so treat Fe as stable only within about 0.5 units. Pb is unstable (-6.1 with no cutoff, -4.3 at 2,000 px, +0.2 at 5,000 px), because few Pb colonies at doses 20-30 pass any cutoff.
4. **Genetic structure matters more than species labels.** R. mucilaginosa populations differ in baseline a\* in all four metals (omnibus p < 0.001, corrected across the 20 omnibus tests). Population 5 is lowest everywhere (Cr -3.4, Cu -2.7, Fe -4.2, Pb -3.6 against population 1). Species effects are metal-specific. R. toruloides, R. dairenensis and R. diobovata have lower estimates than R. mucilaginosa in all four metals, but after correction each is individually significant in only one of them.
5. **Phylogenetic signal is strong.** Pagel's lambda for baseline a\* is 0.97 (Cr), 0.95 (Cu), 0.90 (Pb) for all strains with a tree tip, and the same inside R. mucilaginosa alone (0.98, 0.95, 0.90). Fe is lower (0.51; 0.33 size-adjusted).
6. **Regression to the mean inflates the baseline-vs-change pattern.** Selecting the high-baseline third of strains on half of the replicate wells and measuring the change on the other half gives a smaller gap than the naive analysis: Cu -1.6 against -5.0, Cr -3.5 against -5.1, Pb -4.3 against -5.5. Fe is unaffected (+1.8 against +1.6).
7. **Strain-specific change in a\* is only partly repeatable.** Correlation of the change between replicate halves: Fe 0.70, Pb 0.56, Cr 0.49, Cu 0.11.
8. **Run (batch) is a small share of variance (0-3%). Only Cr shows clear heterogeneity of the dose effect between runs** (I2 = 0.77, Q p = 0.004; Cu 0.49, p = 0.12; Fe 0.41, p = 0.18; Pb 0.53, p = 0.095). The Cr effect weakens from -13.8 (run d000320) to -4.3 (run d000323). Each run is a different strain subset, so run and strain set cannot be separated.
9. **The dose-response shape is shared across species (with one common size slope).** At the top dose, Cu is -11.6 to -12.7 in every species and Fe is +2.0 to +8.9 in every species. Cr is biphasic in every species (+1.7 to +3.8 at doses 0.2-0.4; -4.5 to -14.3 at dose 1.2). Pb is biphasic in every species (+2.8 to +10.8 at doses 5-15; -5.7 to +4.5 at doses 20-30). Species and populations differ mainly in baseline level. The population-by-dose test is also significant in all four metals (F 6.9-12.7), so population dose responses differ somewhat as well (section 11).
10. **Colour and morphology (section 12).** b\* moves much less than a\*. The hue angle shifts from red toward yellow in all four metals (+1.0 to +3.8 SD, total effect). Solidity falls and eccentricity rises under stress, mainly in Cr and Pb. Cu b\* jumps between dose 0 and dose 5 while a\* and size barely change. This looks like a plate or media difference and is not interpreted.
11. **Zinc is usable for strain-level questions (section 13).** With a window ending at 80 h, 8 of 9 plates give data (152 strains; 150 after the area cutoff). Strain a\* ranks agree between plates (Spearman 0.40-0.85). The two Zinc sets show no detectable difference at the shared dose 15. The population-mean dose effect cannot be separated from the plate effect.

## 1. Models (short)

- Source: DuckDB table `heavy_metal_measurement` (517,371 colony rows; 5 metals). Named strains only (502,051 rows). a\* is `ColorLab_a*GeoMedian` (the geometric-median CIELAB a\* of the colony pixels); b\* and L\* are the matching `GeoMedian` columns.
- Models (lme4, REML). The overall models (Tables 1 and 2, Figures 1-3) use strain and plate random effects. Plate is unique within a run, so the overall models have no separate run term. Tables 1 and 2 differ: Table 1 has a strain intercept and a correlated strain dose slope; Table 2 (dose as a factor) has strain and plate intercepts only. The species, population, regime, size-bin and per-run models have a strain intercept, an uncorrelated strain dose slope, and run and plate intercepts (per-run models have no run term). The per-species strata in section 11 (Table 26) have strain and plate intercepts only. Models without a strain slope give standard errors that are too small when strains respond differently, so Table 2 and Table 26 intervals are optimistic.
- `dose_s` is concentration divided by the maximum concentration. "At fixed size" adds ln(area) centred on the control mean. The size term is one common slope per model.
- Dose units are not given in the source. Cr spans 0-1.2 and the other metals 0-30. Metals are never pooled.
- Zinc is treated as one experiment (two runs, two strain sets, one plate per dose). It is shown for completeness in sections 2 and 13 and is left out of the stratified models.

## 2. Overall dose response

{{IMG:fig1_astar_vs_size_by_dose|Figure 1. a* against colony size by dose (well medians, binned). Overlapping curves mean a size effect only. Separated curves mean a shift at the same size.}}

{{T_MM}}

*Table 1. Dose effect on a\* from 0 to the top dose. "total" has no size term. "at fixed size" adds ln(area). "strain share" is strain variance / (strain + plate + residual) at dose 0 and at the top dose. Both Pb models are singular fits.*

{{T_M2}}

*Table 2. a\* difference from 0 dose at fixed size (dose-as-factor model). Standard errors are 0.2-0.6 for Cr, Cu and Fe, 0.4-0.5 for Pb. The Zn standard errors (0.7-1.4) are not valid: with one plate per dose the plate variance collapses to zero, so plate-to-plate noise is ignored.*

{{IMG:fig2_baseline_vs_stressed|Figure 2. Per-strain a* (top) and ln area (bottom), unstressed against top dose (mean of replicate wells, bars = SE). For Zinc the top dose shown is 15, the highest dose shared with dose 0.}}

{{IMG:fig3_strain_baseline_vs_induction|Figure 2b. Model-based strain effects: baseline a* (x) against the strain-specific extra change in a* at the top dose (y), total (top row) and at fixed size (bottom row). Singular fits are not drawn. Points that lie on a straight line (a correlation near -1 or +1) mean the random-effect correlation is at the model boundary. Those points are a rescaling of one random effect and carry no separate information.}}

{{IMG:fig4_repeatability|Figure 3. Share of a* variance due to strain at 0 dose and at the top dose (size-adjusted model).}}

**Reading.** In Cr, Cu and Pb, strains converge to a common low a\* at the top dose (strain share falls to 0.31, 0.31 and 0.10 against 0.52, 0.51 and 0.43 at dose 0). Fe keeps its strain differences (0.83 to 0.92).

## 3. Which traits correlate with a\*

{{IMG:fig5_astar_trait_correlations|Figure 4. Spearman correlation of a* with the traits most correlated with it (well level, all doses).}}

All doses (without `ColorLab_a*Medoid`, which is a\* itself):

{{T_COR_ALL}}

Control dose only:

{{T_COR_0}}

Within strain and dose (replicate wells only; Zinc has no replicates):

{{T_COR_WITHIN}}

**Reading.** Saturation, b\* and the colour-variance traits correlate most strongly with a\*. These share colour information with a\*. Size traits (area, radius, Feret diameter, integrated intensity) correlate at 0.4-0.65 in Cr, Cu, Pb and Zn, but with the opposite sign in Fe. Many size traits are near-duplicates, so treat them as one family. Replicate wells within a strain and dose still show a size-a\* correlation of about 0.3-0.5, so the link is not only between strains.

## 4. Species and population (stratification item 1)

Species is a fixed effect (reference R. mucilaginosa). Strain, run and plate are random. Species with at least 5 strains are included. Population labels exist for 201 R. mucilaginosa strains (`analysis/gwas/data/prior_run_state/pop_assignment_at_run.csv`).

{{T_OM}}

*Table 3. Omnibus tests (F test). `p_BH` is the Benjamini-Hochberg adjustment across all 20 rows, which mixes species and population tests, baseline and dose-response tests.*

{{IMG:s1_species_baseline_forest|Figure 5. Species effect on baseline a* relative to R. mucilaginosa (size-adjusted).}}

{{T_SB}}

*Table 4. Species baseline contrasts (size-adjusted).*

{{IMG:s1_species_slope_forest|Figure 6. Species-specific change in a* over the full dose range at fixed size.}}

{{T_SS}}

*Table 5. Species-specific dose slopes (a\* change over the full dose range).*

{{T_PO}}

*Table 6. R. mucilaginosa against all other species pooled (contrast only).*

{{IMG:s1_population|Figure 7. Populations within R. mucilaginosa: baseline a* (top) and dose response (bottom).}}

{{T_PB}}

*Table 7. Population baseline contrasts against pop1 (size-adjusted).*

{{T_PP}}

*Table 8. Population dose slopes.*

{{IMG:fig6_baseline_astar_by_species|Figure 8. Baseline a* by species (strain means; species with at least 5 strains).}}

**Reading.**
- The population effect is the most consistent result here. It holds in four independent metal screens. Population 5 is lowest in all four metals. Populations 6 (Cr, Fe, Pb) and 3 (Cr, Fe) are also low. Population 4 has the highest or tied-highest estimate in every metal, but differs from pop1 significantly only in Fe (and borderline in Cu and Pb).
- Species effects are inconsistent across metals (for example sp_clade_I is +3.0 in Cr and -2.1 in Pb). The non-mucilaginosa species have 5-18 strains each, so single-species contrasts are noisy. Cu shows no species effect (F p = 0.22 after size adjustment; 0.24 after correction).
- The species effect on baseline a\* shrinks in Pb once colonies below 2,000 px are removed (F 15.8 with no cutoff, 6.0 with the cutoff). The species effect on the Pb dose response is no longer significant (p = 0.059 after correction). Part of the earlier Pb species signal came from tiny colonies.
- Pooled other species are lower than R. mucilaginosa in baseline a\* in all four metals (-0.5 to -1.8). They respond differently in Cr (+1.1) and Pb (+1.5).

## 5. Phylogenetic signal (item 2)

{{IMG:s2_phylo_signal|Figure 9. Pagel's lambda of baseline a* and of the change in a* at the top dose.}}

{{T_PH}}

*Table 9. Pagel's lambda (maximum likelihood under Brownian motion on the PHYling FastTree). 249-264 strains with a tree tip in Cr, Cu and Pb; 196 in Fe; 75 in Zn.*

{{IMG:s2_distance_decay|Figure 10. Mean absolute difference in size-adjusted baseline a* between pairs of strains, by patristic distance.}}

{{IMG:s2_tree_with_astar|Figure 11. The tree with baseline a* per metal.}}

**Reading and limits.**
- Closely related strains share baseline a\*. Lambda is 0.90-0.97 for baseline a\* in Cr, Cu and Pb, and also inside R. mucilaginosa alone (0.90-0.98). The pairwise difference in a\* increases with patristic distance. Fe is weaker (0.51 raw, 0.33 size-adjusted).
- **Blomberg's K was dropped.** The tree has 22 zero-length tips and 74 tips with a neighbour closer than 1e-5, so the covariance matrix is nearly singular (condition number about 1e10). K varied from 1e-7 to 4e-2 with the diagonal jitter (Table 10), so it is not interpretable. Lambda was stable across jitter 1e-8 to 1e-2.
- Lambda near 1 partly reflects near-identical (clonal) strains with similar a\*. It does not show that a\* evolves under Brownian motion.
- Lambda of the change in a\* at the top dose (0.73-1.00) mostly restates the baseline. The change is top-dose a\* minus baseline a\*, and top-dose a\* sits on a floor with little variation (Pb SD 0.77 across strains), so the change is mostly minus the baseline.
- Strains within species and populations are not independent (lambda 0.9-0.97). The species and population tests in section 4 treat strain as the only grouping, so their p-values are too small.
- 266 of 321 strains have a unique tree tip. Strains without a tip are excluded.

{{T_SJ}}

*Table 10. Sensitivity of lambda and K to the diagonal jitter added to the covariance matrix (run on the unfiltered data).*

## 6. Dose regimes (item 3)

A dose is **sub-inhibitory** when the median (across strains) of area at that dose / area at dose 0 is at least 0.5. Otherwise it is **inhibitory**. The 0.5 cutoff is a choice, not a measured threshold. The regime is defined from colony size, and size also relates to a\*, so the split is descriptive and partly circular.

{{IMG:s3_regimes|Figure 12. Top: area ratio by dose with the 0.5 cutoff. Bottom: a* difference from 0 dose, with and without size.}}

{{T_SR}}

*Table 11. Colony area relative to 0 dose.*

{{IMG:s3_regime_slopes|Figure 13. Slope of a* per 10% of the maximum dose, by regime.}}

{{T_RS}}

*Table 12. Regime-specific slopes (a\* change per 10% of the metal's maximum dose; random effects for strain, run and plate).*

**Reading.**
- Cr: sub-inhibitory doses (0-0.4) raise a\* by +0.9 per 10% of the range. Inhibitory doses lower it (-0.8 at fixed size, -1.1 without size).
- Pb: sub-inhibitory doses (0-10) raise a\* by +2.0 (without size) to +2.4 (at fixed size) per 10%. Inhibitory doses (15-30) lower it (-0.6 at fixed size, -1.6 without size). The inhibitory estimate rests on few wells after the area cutoff (section 10).
- Cu: no inhibitory regime except the top dose, yet a\* falls steadily (-1.4 per 10% at fixed size). Cu lowers a\* without reducing colony size.
- Fe: area never falls below 0.66 of control, and a\* rises slightly (+0.4 per 10% at fixed size).

## 7. Size-matched comparison (item 4)

Colonies are placed in five size bins (quantiles of ln area within each metal, among colonies of at least 2,000 px). The dose effect is estimated inside each bin.

{{IMG:s4_size_matched|Figure 14. Wells per size bin and dose (top), a* against dose within bins (middle), and a* spread by size bin (bottom).}}

{{T_S4}}

*Table 13. Dose effect within size bins. `slope_full_range` is the a\* change over the full dose range; `change_over_observed_range` multiplies it by the highest dose in the bin (random effects: strain intercept and slope, run, plate).*

**Reading.**
- Cu: a\* falls with dose in every size bin (-12.4 to -14.2 over the full range), and all seven doses are present in every bin. This is the cleanest size-matched result: Cu lowers a\* at the same size.
- Fe: a\* rises in every bin (+2.2 to +4.4).
- Cr: negative in every bin (-14.1 to -2.0). The effect in bin S4 (-2.0) is not significant (p = 0.16), and S4 lacks dose 1.2.
- Pb: size and dose stay confounded. The three largest bins contain only doses 0-15. Their slopes (+19 to +23 per full dose range) extrapolate twice beyond the observed range. The change over the observed range is +9.7, +11.4 and +11.2 (last two columns of Table 13). The smallest bin is negative (-5.0) and the second smallest positive (+8.7). There is no single size-matched Pb estimate.
- **Coverage is uneven.** High doses populate the small bins, and some size-by-dose cells in Cr and Pb hold few wells. Bins are quantiles of pooled ln area, so dose and size remain confounded inside a bin.

## 8. Selection on baseline, split-half (item 5)

Strains need at least 2 replicate wells at dose 0 and at the top dose. Replicate wells are split at random into halves A and B (1,000 splits). **Naive:** select the top and bottom baseline thirds from all wells and measure the change in the same wells. **Split-half:** select on half A, measure the change in half B.

{{IMG:s5_split_half_summary|Figure 15. Gap in the change in a* between the high- and low-baseline thirds (left) and reliability of the strain-specific change (right).}}

{{IMG:s5_split_half_scatter|Figure 16. Baseline against change in a*: same wells (top) and independent halves (bottom).}}

{{T_S5}}

*Table 14. Size-adjusted a\*. `naive_gap` and `split_gap_mean` are the change in the top-baseline third minus the bottom third (95% CI from a bootstrap over strains). `reliability_change_A_vs_B` is the Spearman correlation of the strain-specific change between the two halves.*

**Reading.**
- Part of the strong negative baseline-vs-change relationship in Cr, Cu and Pb is regression to the mean. The gap shrinks from -5.0 to -1.6 in Cu, from -5.1 to -3.5 in Cr and from -5.5 to -4.3 in Pb. Fe is unchanged (+1.6 naive, +1.8 split).
- In Cu the strain-specific change is barely repeatable (0.11), and baseline a\* in Cu replicate wells is also noisy (0.25 between halves).
- Only 62 Pb strains have enough wells at dose 30 after the area cutoff, so the Pb estimate is uncertain.
- The reliabilities are for halves of 1-2 wells and understate full-data reliability.
- "Size-adjusted" here is the residual of a\* on ln area from a line fitted to the dose-0 wells of each metal. That differs from the ln(area) term in the mixed models, and it extrapolates a linear slope to small top-dose colonies.
- The cutoff shrinks the sample (Cr 289 to 148 strains, Pb 273 to 62; section 14).

## 9. Batch (item 6)

Each metal is analysed separately, with run and plate as random effects (all models above). Each run is a different set of 67-84 strains across all doses, so run and strain set are confounded. Zinc is excluded.

{{IMG:s6_batch|Figure 17. Dose effect at fixed size by run (left four panels) and variance components (right).}}

{{T_H6}}

*Table 15. Heterogeneity of the dose effect across runs (I2, Cochran Q).*

{{T_P6}}

*Table 16. Dose effect (a\* change over the full dose range at fixed size) per run.*

{{T_V6}}

*Table 17. Variance components (shares of random variance).*

**Reading.** Run explains 0-3% of the random variance. Plate explains 35% (Cr) and 39% (Pb), but only 2-4% in Cu and Fe. Only Cr shows clear heterogeneity of the dose effect between runs (I2 = 0.77, Q p = 0.004). The Cr effect shrinks from run d000320 to d000323 (-13.8, -10.5, -6.4, -4.3). This could be a run-order effect or a difference between strain subsets. The Cu runs lie within -12.8 to -14.2 and the heterogeneity is not significant (p = 0.12). These models have a random strain dose slope; models without it gave smaller standard errors and larger I2.

## 10. Colony area cutoff (2,000 px)

Very small colonies may have a\* dominated by background pixels. Every image was used to test this. Colonies are small early at every dose. So unstressed small colonies (early) can be compared with stressed small colonies (late) at the same area.

{{IMG:area_floor_by_metal|Figure 18. Median a* against colony area, unstressed (dose 0) and top dose, all images. The dotted line is 2,000 px.}}

{{T_S12_BINS}}

*Table 18. Median a\* in area bins for unstressed colonies (dose 0). Bins with fewer than 30 colonies are blank.*

{{T_S12_HINGE}}

*Table 19. Change-point of a\* against ln area (piecewise linear fit; area in px where the slope increases).*

**Reading.**
- **a\* depends on area in every metal, even without stress.** At dose 0, median a\* is about 0-4 below 1,500 px, about 4-5 at 2,000 px, about 8 at 4,500 px and about 20 at 33,000 px. The curve is nearly the same in Cr, Cu, Fe, Pb and Zn. The first-48-h curve (dashed) matches the all-times curve, so area and colony age cannot be separated.
- **Cr and Cu: stressed and unstressed colonies have the same a\* below about 2,000 px.** a\* cannot tell them apart there. Above that, stressed colonies are clearly lower.
- **Fe and Pb: stressed colonies of the same small size have higher a\* than unstressed ones.** Fe top-dose a\* is about 10.5 at 200-2,000 px, against 3-4 at dose 0. Pb top-dose a\* is 5-6.6 at 200-1,500 px, against 0-2.6 at dose 0. So the Pb floor near 6 is partly a stress effect and not only a pixel effect.
- **Zinc has no unstressed colonies below 3,000 px.** The cutoff cannot be tested for Zinc directly. The assumption is that it transfers, because all runs use the same imager. Zinc top-dose a\* is flat (about 3-4.6) from 500 to 9,000 px.
- **The change-point estimates differ by metal** (Table 19: 1,700-1,800 px for Cr and Cu at dose 0; 3,200-4,400 px for Pb, Fe and Zn). The slope below the change-point is 0.05-0.9 for Cr and Cu but 1.7 for Pb and Fe, against 5-7 above, so "flat" holds for Cr and Cu and only "low and slowly rising" for Pb and Fe. Estimates for Zn and Fe sit at the edge of their data. The per-bin medians show a step at about 1,500-2,000 px in all four metals with small colonies. The 2,000 px value was therefore chosen from the data and applied to all metals. It is assumed, not tested, for Pb, Fe and Zn.
- **The cutoff does not remove the size confound.** a\* keeps rising with area above 2,000 px, so every model still includes a size term.

**What the cutoff removes** (share of wells in the unfiltered table that are lost):

{{T_RETAINED}}

*Table 20. Wells in the unfiltered table and in the main table, by metal and dose. Most losses are at the top doses of Cr, Pb and Zn.*

{{IMG:s8_minarea|Figure 19. Sensitivity to the area cutoff: dose effect against threshold (top), per-dose effect at fixed size (middle), wells retained by dose (bottom). Threshold 0 is the original unfiltered run.}}

{{T_S8_MODEL}}

*Table 21. Dose effect on a\* (at fixed size) by minimum colony area (0 = no filter). Zinc is not included.*

{{T_S8_CR}}

*Table 22. Chromium: a\* difference from 0 dose at fixed size, by threshold, and wells retained.*

{{T_S8_PB}}

*Table 23. Lead: a\* difference from 0 dose at fixed size, by threshold, and wells retained.*

{{T_S8_TOP}}

*Table 24. Strain-mean a\* at the top dose by threshold.*

**Reading.**
- **Cu is unaffected** (-13.4 to -13.3). **Fe is stable within about 0.5 units** (+3.8 with no cutoff, +4.3 to +4.4 with any cutoff). The step from 3.8 to 4.4 comes from removing 9 wells, and the Fe size slope moves from 2.9 to 3.9 at the same time, so the Fe fixed-size estimate depends on a few small colonies.
- **Survivor selection.** After the cutoff the top-dose estimates describe larger colonies. In Cr, 241 strains remain at the top dose (311 without a cutoff) and the median top-dose area rises from 1,729 to 3,281 px (Table 24). In Pb, 27-30% of the unfiltered wells at doses 20-30 remain. Compare the main and unfiltered numbers with this in mind.
- **Cr: the top-dose drop is not a small-colony artifact.** The effect at dose 1.2 is -9.1 with no filter and -10.7 with colonies of at least 5,000 px. The overall Cr dose effect stays between -8.4 and -9.3.
- **Pb: the rise at doses 5-15 holds or strengthens** (+4.9, +8.0 and +9.5 at 2,000 px, against +4.4, +7.5 and +7.6 with no filter). At doses 20-30 only 242-265 wells remain at 2,000 px and 8-17 at 5,000 px, so the Pb inhibitory range cannot be tested. The overall Pb dose effect is unstable across cutoffs and should not be quoted as one number.
- **A floor remains in Pb.** Among Pb colonies of at least 3,000 px at dose 30 (81 strains), strain-mean a\* is 5.7 with SD 0.56 across strains.

## 11. Stratified by species

R. mucilaginosa (216 strains) is shown alone. Each other species with at least 5 strains in at least two metals is shown alone against R. mucilaginosa as a grey reference. R. graminis (5 strains, Cr only, 2 doses) is in the tables but not the plots. Zinc is left out: only R. mucilaginosa has enough strains there. Lines show strain means with 95% CI across strains.

{{T_S9_COUNTS}}

*Table 25. Strains per species and metal.*

### 11a. R. mucilaginosa alone

{{IMG:sp_muc_a_vs_size|Figure 20. R. mucilaginosa only: a* against colony size by dose.}}

{{IMG:sp_muc_by_population|Figure 21. R. mucilaginosa only: a* (top) and ln colony area (bottom) by dose, one line per population.}}

{{IMG:sp_muc_baseline_vs_top|Figure 22. R. mucilaginosa only: baseline against top-dose a* per strain, coloured by population.}}

**Reading.** The curves of the six populations have similar shapes in Cr, Cu and Pb and differ mainly in level. The population-by-dose test is nevertheless significant in all four metals (F 12.7 Cr, 6.9 Cu, 7.1 Fe, 8.6 Pb; p < 0.001; Table 3), so the shapes also differ. For example the Pb change over the full dose range is -1.9 for pop5 and -5.0 to -7.5 for the other populations (Table 8). Population 4 stands out in Fe (a\* 24-27 at doses 15-30 against 18-22 for the other populations) and in the Pb peak at dose 10 (a\* about 31 against 23-27).

### 11b. Each other species alone

{{IMG:sp_species_dose_response_a|Figure 23. a* by dose for each species (colour) against R. mucilaginosa (grey).}}

{{IMG:sp_species_dose_response_size|Figure 24. Colony size (ln area) by dose for each species against R. mucilaginosa.}}

{{IMG:sp_species_effect_fixed_size|Figure 25. Dose effect on a* at fixed size for each species against R. mucilaginosa (95% CI). One model per metal with a dose-by-species interaction and one common size slope, as in Table 5.}}

{{T_S9_MODEL}}

*Table 26. Total dose effect on a\* over the full dose range (no size term), by species. Random intercepts for strain and plate only (plate is nested in run), so the intervals are optimistic. "Other species" pools all species with at least 5 strains except R. mucilaginosa. A per-species size slope cannot be estimated reliably in small strata, because dose and size are confounded, so fixed-size effects by species come from the shared-slope models (Table 5, Figure 25).*

{{IMG:sp_species_baseline_vs_top|Figure 26. Baseline against top-dose a* per strain, each species (colour) over R. mucilaginosa (grey).}}

{{IMG:sp_chromium_a_vs_size_by_species|Figure 27. Chromium: a* against colony size by species (colour = dose).}}

{{IMG:sp_copper_a_vs_size_by_species|Figure 28. Copper: a* against colony size by species (colour = dose).}}

{{IMG:sp_iron_a_vs_size_by_species|Figure 29. Iron: a* against colony size by species (colour = dose).}}

{{IMG:sp_lead_a_vs_size_by_species|Figure 30. Lead: a* against colony size by species (colour = dose).}}

**Reading.**
- **Cu is the same in every species.** The effect at the top dose and fixed size is -11.6 to -12.7 in all eight groups (Figure 25). Fe is positive in every species (+2.0 to +8.9 at dose 30).
- **Cr and Pb are biphasic in every species** (Figure 25). Cr is +1.7 to +3.8 at doses 0.2-0.4 and -4.5 to -14.3 at dose 1.2. Pb is +2.8 to +10.8 at doses 5-15 and -5.7 to +4.5 at doses 20-30. The Pb curves at doses 20-30 rest on 39-124 wells per species after the area cutoff, and few of those at the top doses.
- **R. sp_clade_I responds most strongly in Cr** (-14.3 at dose 1.2 against -9.4 for R. mucilaginosa). R. sphaerocarpa responds least (-4.5).
- **Total effects** (no size term, Table 26) range from -8.6 to -16.4 in Cr, -14.9 to -16.9 in Cu, -0.1 to +5.4 in Fe and -10.0 to -20.2 in Pb. R. mucilaginosa has the largest Pb total effect (-20.2).
- **Species mostly differ in baseline level, not in shape.** For example R. sphaerocarpa, R. diobovata and R. dairenensis start at 14-16 in Pb against 19 for R. mucilaginosa. All species converge to a\* of about 5-6 at Pb doses 20-30.
- **Single-species curves are noisy.** Most species have 5-10 strains, so the confidence intervals are wide, especially in Fe.

## 12. b\* and morphology

Colour traits: L\*, a\*, b\*, chroma (square root of a\*^2 + b\*^2) and hue angle (atan2(b\*, a\*) in degrees; 0 is red and 90 is yellow). Morphology: circularity, solidity, eccentricity, compactness, extent and aspect ratio (major / minor axis). Wells and window are the main dataset. For the models, each trait is scaled by the mean and SD of the 0-dose wells of its metal, so effects are in SD units of the unstressed wells. Lines show strain means with 95% CI across strains, for R. mucilaginosa and all other species pooled.

{{IMG:bm_colour_response|Figure 31. Colour traits by dose: R. mucilaginosa and other species.}}

{{IMG:bm_ab_plane|Figure 32. Trajectory in the a*-b* plane under stress (colour = dose; labels give the concentration).}}

{{IMG:bm_morphology_response|Figure 33. Morphology by dose: R. mucilaginosa and other species. Shape metrics of very small colonies are noisy.}}

{{IMG:bm_dose_effect_heatmap|Figure 34. Mixed-model dose effect on colour and morphology traits in SD units (stars: unadjusted p < 0.05, 0.01, 0.001). Rows of panels: main dataset and no area filter. Columns: total effect and effect at fixed size.}}

{{T_S10_TOTAL}}

*Table 27. Dose effect over the full dose range, main dataset, total effect (SD of the unstressed wells).*

{{T_S10_FIXED}}

*Table 28. Same, at fixed colony size.*

{{T_CUB}}

*Table 28b. Copper: strain-mean a\*, b\* and ln area by dose (main dataset).*

{{IMG:bm_morph_colour_corr|Figure 35. Strain-level Spearman correlation of morphology (rows) with colour (columns): baseline (top) and change at the top dose (bottom).}}

{{IMG:bm_species_baseline_b_circ|Figure 36. Baseline b* and circularity by species (strain means at 0 dose).}}

**Reading.**
- **b\* moves much less than a\*** (total effect: Cr -1.2, Cu -0.3, Pb -0.4 SD; a\* is -2.8, -3.2, -4.3). In Fe both rise (b\* +2.3, a\* +0.8 SD).
- **The hue angle rises in all four metals** (+1.0 Cr, +3.0 Cu, +1.6 Fe, +3.8 Pb SD, total). Colour shifts from red toward yellow with stress, mainly because a\* falls faster than b\*.
- **Lightness (L\*) falls in Cr, Cu and Fe and rises in Pb** (+7.0 SD). Pb colonies at doses of 20 and above are tiny and pale.
- **Colonies become less compact under stress.** Solidity falls (Cr -7.0 SD, Pb -4.9, Fe -2.4, Cu -2.1) and eccentricity rises (Cr +4.3, Pb +3.6). The effect is largest in the two metals that inhibit growth most. At fixed size most of the solidity effect disappears (Cr -1.5, Pb +0.7), so it follows colony size.
- **Do not read the Cr aspect-ratio and compactness effects as magnitudes** (+31 and +129 SD for aspect ratio, total and fixed size; +29 SD for compactness at fixed size). The unstressed wells have almost no variance in those traits, so the SD scale blows up. Only the direction (more elongated, less compact) is meaningful.
- **Cu b\* jumps from 18.8 at dose 0 to 32.0 at dose 5** (Table 28b) while a\* (18.7 to 19.4) and ln area (10.3 to 10.5) barely change. This looks like a plate or media effect between the dose-0 plate and the treated plates, not biology. Cu b\* and hue results should be read with that in mind. The same kind of step appears in Zinc (b\* 5.2 to 9.1 between dose 0 and 5, set A) and is not interpreted there either.
- **Strain-level links between morphology and colour are weak at baseline in all four metals** (|rho| at most 0.27; Table not shown, Figure 35). The largest baseline links are in Pb (a\* with solidity +0.14, eccentricity -0.20, aspect ratio -0.23) and Cu (L\* with solidity +0.27), and they track colony size (ln area with a\*: +0.15 in Pb). For the change at the top dose, Cu shows the strongest links: the change in b\* correlates with the change in ln area (rho 0.61) and in solidity (0.58). In Fe the change in aspect ratio correlates with the change in a\* (-0.34) and L\* (+0.41).
- Morphology of very small colonies is unreliable (few pixels). The main dataset already drops colonies below 2,000 px.

## 13. Zinc

Zinc has 9 plates in 2 runs, one plate per dose, and 1 well per strain per dose. Run d000388 (set A) holds 80 strains at doses 0, 5, 10 and 15. Run d000390 (set B) holds about 70 other strains at doses 10-30. The two runs are treated as one experiment. The two sets share no strains.

{{IMG:zn_window_coverage|Figure 37. Zinc: usable wells per plate for different window end points, and plate imaging spans.}}

{{T_ZN_COVER}}

*Table 29. Wells per Zinc plate by window end (24 h window, at least 2 images).*

**Rescue.** The plate-relative late window (end 107.6 h) dropped plate d000388/41 (dose 15, stops at 83.7 h) and cut plate d000390/100 (dose 30, stops at 89.7 h). A window ending at 80 h keeps both. The only plate that cannot be rescued is d000390/96 (dose 10, 7 images up to 35.8 h). Result: 8 plates, 588 wells and 152 strains before the area cutoff. **Everything in this section uses the main dataset (colonies of at least 2,000 px)**, which leaves 502 wells and 150 strains. Table 30b shows the same series without the cutoff.

{{IMG:zn_paired_dose_response|Figure 38. Zinc, window [56, 80] h, colonies of at least 2,000 px: strains followed across doses within each set. The mean line also contains the plate effect.}}

{{IMG:zn_composite_series|Figure 39. Zinc as one experiment: set A (doses 0-15) and set B (doses 15-30) bridged at the shared dose 15.}}

{{T_ZN_COMP}}

*Table 30. Zinc composite series (main dataset): strain means per dose, with 95% CI across strains.*

{{T_ZN_COMP_NF}}

*Table 30b. The same series without the area cutoff (588 wells). The dose-25 and dose-30 means are lower because small colonies are kept.*

{{IMG:zn_strain_level_consistency|Figure 40. Zinc: strain rank agreement between dose plates (left) and agreement of the strain-level change with the same strains in other metals (right).}}

{{T_ZN_RANK}}

*Table 31. Strain rank agreement (Spearman) between pairs of Zinc dose plates, with the mean difference in the trait.*

**What can and cannot be said.**
- **Strain comparisons within a plate are largely protected from the plate effect.** A plate effect that adds the same amount to every strain on the plate does not change ranks or between-strain differences. This assumes the plate effect is additive. With one well per strain per plate, position effects and strain-by-plate interaction cannot be separated from the strain, so rank agreement shows that strain differences are reproducible despite them, not that they are absent. Strain a\* ranks agree between different dose plates (Spearman 0.84 between doses 0 and 5, 0.55 between 0 and 10, 0.40 between 0 and 15; set B 0.85 between doses 15 and 20, 0.84 between 15 and 25 with 34 strains, 0.58 between 15 and 30 with 34 strains). Strain ranks in colony size agree between doses 5 and 15 (0.61-0.82), but unstressed size does not predict stressed size (0.33, -0.11 and 0.04).
- **The two sets show no detectable difference at the shared dose 15** (set A against set B: a\* 18.7 against 18.4, p = 0.81; ln area 9.02 against 9.12, p = 0.18; b\* 10.1 against 10.3). With 73 and 67 strains this is a non-rejection and not proof of equivalence. It supports reading the two sets as one composite series. It does not prove it.
- **Composite shape.** Total a\* rises from 17.3 (dose 0) to 20.3 (dose 5), stays near 18-19 to dose 15, then falls to 7.3 at dose 30. b\* rises from 5.2 to 10.2 at dose 15 and falls to 6.6. ln area falls from 10.2 to 8.4. At fixed size the effect at dose 30 is -0.3 (Table 2), so the total fall in a\* at doses 20-30 goes together with smaller colonies. The pattern resembles Cr and Pb, but Zinc cannot show it independently.
- **The dose-5 rise (+3.0 in a\*, +3.9 in b\*) may be plate noise or a media step.** In the other metals the plate explains up to 39% of the variance. The same b\* step appears in Cu and is not interpreted.
- **The population-mean dose effect cannot be separated from the plate effect.** There is one plate per dose. Zinc dose-effect standard errors and p-values are not valid, because the plate variance collapses to zero. The overall Zinc slopes (-14.4 total, -3.5 at fixed size) also compare disjoint strain sets across doses.
- **Set B has no unstressed baseline.** Its strains are not on the dose-0 plate.
- **Survivor selection.** After the cutoff only 34-36 strains remain in set B at doses 25 and 30. Their means describe the larger colonies (Table 30b shows what changes).
- **Cross-metal strain consistency is weak and exploratory.** The strongest links in set A (73 strains) are |rho| 0.51-0.57: Zinc colony size at dose 15 against the Fe change in a\* (-0.57), and the Zinc change in ln area against the Cu change in ln area (+0.56). In set B (27 strains) the Zinc change in a\* correlates +0.53 with the Cr change. Many pairs were compared, with no multiplicity correction.
- **The 2,000 px cutoff for Zinc is assumed, not tested** (section 10). It removes 19% of Zn wells at dose 20 and about half at doses 25-30.

## 14. Sensitivity: results without the area cutoff

The original analysis used no area cutoff and the median-span window for Zinc (107.6 h). Key numbers from that run (tables in `report/sensitivity_unfiltered/`):

{{T_SENS_MM}}

*Table 32. Dose effect (0 to top dose) in the original unfiltered run and in the main run.*

{{T_SENS_OTHER}}

*Table 33. Other key numbers, original and main.*

**Reading.**
- **Cu and Cr are the same in both runs; Fe moves by about 0.5.**
- **Pb changes most.** The Pb dose-effect estimate moves (-6.1 to -4.3). The species effect on baseline a\* without size adjustment falls (F 15.8 to 6.0), but with size adjustment it does not change (4.3 to 4.3). The species-by-dose test falls from F 8.75 (p < 0.001) to F 2.06 (corrected p = 0.059).
- **Population effects and phylogenetic signal keep their conclusions.**
- **The comparison is not strictly like with like.** The split-half samples shrink (Cr 289 to 148 strains, Pb 273 to 62, because few strains keep at least 2 top-dose wells), so those rows also reflect fewer strains. Run heterogeneity is not compared: the original run used models without a strain dose slope, and the main run uses the slope.
- **Zinc changes** because the 80 h window restores plate d000388/41 (dose 15, set A) and the cutoff removes small colonies at the top doses.

## Limits

- Dose units are not given in the source, and metals are not pooled.
- Species and population tests treat strains as independent, although relatives share a\* (lambda 0.9-0.97). Their p-values are too small. Benjamini-Hochberg is applied across all 20 omnibus tests together.
- The dose-factor model (Table 2) and the per-species strata (Table 26) have no strain dose slope, so their intervals are optimistic.
- The regime split (section 6) is defined from colony size, which also relates to a\*.
- Zinc is not modelled in the stratified analyses. Its plate and dose effects are confounded.
- Size is partly an effect of stress. Models "at fixed size" estimate a direct effect, not the whole induction. One size slope is used for all strains.
- The area cutoff is justified for Cr and Cu from unstressed colonies and assumed for Pb, Fe and Zn. The Pb inhibitory range (doses 20-30) has too few colonies above any cutoff to test.
- Survivor selection: wells that never produced objects are absent, wells with fewer than 2 window images were dropped, and the area cutoff removes small colonies at the top doses. Top-dose estimates therefore describe larger colonies (Cr, Pb, Zn).
- 16 strains still have no species and are excluded from species tests (see the data-problems report). Species tests use species with at least 5 strains. R. mucilaginosa has 174-216 strains and the others have 5-18.
- Population labels come from an earlier GWAS run (201 strains).
- The all-objects sensitivity table is built (`results/wells_allobj.csv`) but not modelled.
- The Cu b\* jump between dose 0 and dose 5 is not explained.
- Old-vs-new strain assignment for Copper is unresolved (see the data-problems report). This report uses the new data's own strain labels.
