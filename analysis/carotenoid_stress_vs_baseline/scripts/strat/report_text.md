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

{{IMG:fig1_astar_vs_size_by_dose|Figure 1. a* against colony size by dose (well medians, binned). Overlapping curves mean a size effect only; separated curves mean a shift at the same size.}}

{{T_MM}}

*Table 1. Dose effect on a\* from 0 to the top dose. "total" has no size term; "at fixed size" adds ln(area). "strain share" is strain variance / (strain + plate + residual) at dose 0 and at the top dose. The Pb total model is a singular fit.*

{{T_M2}}

*Table 2. a\* difference from 0 dose at fixed size (dose-as-factor model; SE 0.2-0.5 except Zn, about 8.8).*

{{IMG:fig2_baseline_vs_stressed|Figure 2. Per-strain a* (top) and ln area (bottom), unstressed against top dose (mean of replicate wells, bars = SE). For Zinc the top dose shown is 15, the highest dose shared with dose 0.}}

{{IMG:fig3_strain_baseline_vs_induction|Figure 2b. Model-based strain effects: baseline a* (x) against strain-specific extra change in a* at the top dose (y), total (top row) and at fixed size (bottom row). Pb total is a singular fit and is not drawn. Pb at fixed size (r = -0.99) and Zn total (r = +1.00) lie on a straight line: the random-effect correlation is at the model boundary, so these points are a rescaling of one random effect and carry no separate information.}}

{{IMG:fig4_repeatability|Figure 3. Share of a* variance due to strain at 0 dose and at the top dose (size-adjusted model).}}

**Reading.** In Cr, Cu and Pb strains converge to a common low a\* at the top dose (strain share falls to 0.24, 0.28 and 0.03). Fe keeps its strain differences (0.82 to 0.91).

## 3. Which traits correlate with a\*

{{IMG:fig5_astar_trait_correlations|Figure 4. Spearman correlation of a* with the traits most correlated with it (well level, all doses).}}

All doses (excluding `ColorLab_a*Medoid`, which is a\* itself):

{{T_COR_ALL}}

Control dose only:

{{T_COR_0}}

Within strain and dose (replicate wells only; Zinc has no replicates):

{{T_COR_WITHIN}}

**Reading.** Saturation, b\* and the colour-variance traits correlate most strongly with a\*; these are partly the same colour information. Size traits (area, radius, Feret diameter, integrated intensity) correlate at 0.5-0.7 in Cr, Cu, Pb and Zn but with the opposite sign in Fe. Many size traits are near-duplicates of each other, so treat them as one family. Replicate wells within a strain and dose still show a size-a\* correlation of about 0.3-0.54, so the link is not only between strains.

## 4. Species and population (stratification item 1)

Species is a fixed effect (reference R. mucilaginosa); strain, run and plate are random. Only species with at least 5 strains are included. Population labels exist for 201 R. mucilaginosa strains (file `analysis/gwas/data/prior_run_state/pop_assignment_at_run.csv`).

{{T_OM}}

*Table 3. Omnibus tests (F test; Benjamini-Hochberg across all rows).*

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
- The population effect is the most consistent result in this section. Population 5 is lowest in all four metals (Cr -3.4, Cu -2.6, Fe -4.0, Pb -3.8). Population 6 is also low in Cr, Fe and Pb, and population 3 in Cr and Fe. Population 4 has the highest or tied-highest estimate in every metal, but differs from pop1 significantly only in Cu and Fe. The same pattern in four independent metal screens suggests a real genetic effect.
- Species effects are inconsistent across metals (for example sp_clade_I is +3.6 in Cr and -2.4 in Pb), and the non-mucilaginosa species have 5-18 strains each, so single-species contrasts are noisy. Cu shows no species effect overall (F p = 0.17).
- Pooled "other species" are lower than R. mucilaginosa in baseline a\* in all four metals (-0.8 to -2.1), and respond differently in Pb (slope difference +4.3).

## 5. Phylogenetic signal (item 2)

{{IMG:s2_phylo_signal|Figure 9. Pagel's lambda of baseline a* and of the change in a* at the top dose.}}

{{T_PH}}

*Table 9. Pagel's lambda (maximum likelihood under Brownian motion on the PHYling FastTree). 254-264 strains with a tree tip in Cr, Cu and Pb; 197 in Fe; 75 in Zn.*

{{IMG:s2_distance_decay|Figure 10. Mean absolute difference in size-adjusted baseline a* between pairs of strains, by patristic distance.}}

{{IMG:s2_tree_with_astar|Figure 11. The tree with baseline a* per metal.}}

**Reading and limits.**
- Closely related strains share baseline a\*. Lambda is 0.90-0.98 for baseline a\* in Cr, Cu and Pb, also inside R. mucilaginosa alone (0.90-0.98), and the pairwise difference in a\* increases with patristic distance (Spearman 0.35 Cr, 0.17 Cu, 0.21 Pb). Fe is weaker (0.51 raw, 0.23 size-adjusted). Baseline colony size has lambda near 0 in Cr, Cu and Pb but 0.83 in Fe and 0.99 in Zn.
- **Blomberg's K was dropped.** The tree has 22 zero-length tips and 74 tips with a neighbour closer than 1e-5, so the covariance matrix is nearly singular (condition number about 1e10). K varied from 1e-7 to 4e-2 with the diagonal jitter (Table 10), so it is not interpretable. Lambda was stable across jitter 1e-8 to 1e-2 (0.89 to 0.99 for Cr, Cu and Pb).
- Lambda near 1 partly reflects near-identical (clonal) strains having similar a\*. It does not show that a\* evolves under Brownian motion.
- 266 of 321 strains have a unique tree tip; strains without a tip are excluded. These are mostly outside the sequenced set.

{{T_SJ}}

*Table 10. Sensitivity of lambda and K to the diagonal jitter added to the covariance matrix.*

## 6. Dose regimes (item 3)

Regimes are defined from colony size: a dose is **sub-inhibitory** when the median (across strains) of area at that dose / area at dose 0 is at least 0.5, and **inhibitory** below that. The 0.5 cutoff is a choice, not a measured threshold.

{{IMG:s3_regimes|Figure 12. Top: area ratio by dose with the 0.5 cutoff. Bottom: a* difference from 0 dose, with and without size.}}

{{T_SR}}

*Table 11. Colony area relative to 0 dose.*

{{IMG:s3_regime_slopes|Figure 13. Slope of a* per 10% of the maximum dose, by regime.}}

{{T_RS}}

*Table 12. Regime-specific slopes (a\* change per 10% of the metal's maximum dose; random effects for strain, run and plate).*

**Reading.**
- Cr: sub-inhibitory doses (0-0.4) raise a\* by +0.9 per 10% of the range; inhibitory doses lower it (-0.7 at fixed size, -1.2 without size).
- Pb: sub-inhibitory doses (0-10) raise a\* by +2.1 to +2.3 per 10%; inhibitory doses (15-30) lower it (-0.6 at fixed size, -1.4 without size). Pb colonies at doses 20-30 are about 4% of control area, so the inhibitory-regime estimate rests on very small colonies.
- Cu: no inhibitory regime except the top dose, yet a\* falls steadily (-1.4 per 10% at fixed size): Cu lowers a\* without reducing colony size.
- Fe: area never falls below 0.66 of control, and a\* rises slightly (+0.4 per 10% at fixed size).

## 7. Size-matched comparison (item 4)

Colonies are placed in five size bins (quantiles of ln area within each metal) and the dose effect is estimated inside each bin.

{{IMG:s4_size_matched|Figure 14. Wells per size bin and dose (top), a* against dose within bins (middle), and a* spread by size bin (bottom).}}

{{T_S4}}

*Table 13. Dose effect within size bins (a\* change over the full dose range; random effects for strain, run, plate).*

**Reading.**
- Cu: a\* falls with dose in every size bin (-10 to -14), and all seven doses are present in every bin. This is the cleanest size-matched result: Cu lowers a\* at the same size.
- Fe: a\* rises in every bin (+3.1 to +4.5).
- Cr: negative in all bins (-5 to -12).
- Pb: size and dose are strongly confounded. The two smallest bins (below about 1,850 px) have a\* of 6.3-7.1 and SD 2.6 at all doses, and the two largest bins contain only doses 0-15. Slopes in the Pb bins differ in sign (+4.7 in the small bins, -14.5 in S3, +17 to +20 in the large bins), so there is no single size-matched Pb estimate.
- **Coverage is uneven.** High doses populate the small bins; in Cr and Pb some size-by-dose cells contain fewer than 10 wells.

## 8. Selection on baseline, split-half (item 5)

Strains with at least 2 replicate wells at dose 0 and at the top dose. Replicate wells are split at random into halves A and B (1,000 splits). **Naive:** select the top and bottom baseline thirds from all wells and measure the change in the same wells. **Split-half:** select on half A, measure the change in half B.

{{IMG:s5_split_half_summary|Figure 15. Gap in the change in a* between the high- and low-baseline thirds (left) and reliability of the strain-specific change (right).}}

{{IMG:s5_split_half_scatter|Figure 16. Baseline against change in a*: same wells (top) and independent halves (bottom).}}

{{T_S5}}

*Table 14. Size-adjusted a\*. `naive_gap` and `split_gap_mean` are the change in the top-baseline third minus the bottom third (95% CI from a bootstrap over strains). `reliability_change_A_vs_B` is the Spearman correlation of the strain-specific change between the two halves.*

**Reading.**
- Part of the strong negative baseline-vs-change relationship in Cr, Cu and Pb is regression to the mean: the gap shrinks from -5.3 to -1.7 in Cu, from -6.9 to -5.5 in Cr and from -6.7 to -4.2 in Pb. Fe is unchanged (+1.6 naive, +1.9 split).
- In Cu the strain-specific change is barely repeatable (0.13), and baseline a\* in Cu replicate wells is also noisy (0.25 between halves).
- The reliabilities are for halves of 1-2 wells and understate full-data reliability.

## 9. Batch (item 6)

Each metal is analysed separately, with run and plate as random effects (all models above). Each run is a different set of about 70-84 strains across all doses, so run and strain set are confounded. Zinc is excluded.

{{IMG:s6_batch|Figure 17. Dose effect at fixed size by run (left four panels) and variance components (right).}}

{{T_H6}}

*Table 15. Heterogeneity of the dose effect across runs (I2, Cochran Q).*

{{T_P6}}

*Table 16. Dose effect (a\* change over the full dose range at fixed size) per run.*

{{T_V6}}

*Table 17. Variance components (shares of random variance).*

**Reading.** Run explains 0-2% of the random variance. Plate explains 34% (Cr) and 38% (Pb) but only 3-4% in Cu and Fe. The Cr dose effect shrinks monotonically from run d000320 to d000323 (-12.6, -9.5, -6.2, -5.6); this could be a run-order effect or a difference between strain subsets.

## 10. Minimum colony area (small-colony check)

{{MINAREA}}

## Limits

- Dose units are not given in the source, and metals are not pooled.
- Zinc is not modelled in the stratified analyses (one plate per dose).
- Size is partly an effect of stress; models "at fixed size" estimate a direct effect, not the whole induction. One size slope is used for all strains.
- Survivor selection: wells that never produced objects are absent, and wells with fewer than 2 window images were dropped (see `report/tables/wells_lost_by_dose.csv`).
- 16 strains still have no species and are excluded from species tests (see the data-problems report). Species tests use species with at least 5 strains; R. mucilaginosa has 174-216 strains and the others have 5-18.
- Population labels come from an earlier GWAS run (201 strains).
- The all-objects sensitivity table is built (`results/wells_allobj.csv`) but not modelled.
- Old-vs-new strain assignment for Copper is unresolved (see the data-problems report). This report uses the new data's own strain labels.
