# Finding: Copper effect on Rhodotorula growth is extent-limited, not rate-limited

- **Date**: 2026-08-15
- **Source**: `analysis/growth_rates/scripts/04_doubling_time.R`, `05_species_cu_sensitivity.R`
- **Topic**: copper-tolerance, growth-phenotyping

## Result

1. **Doubling time is copper-neutral.** Exponential-region doubling time
   (best-R2 4-point window on log(area), per colony) is flat ~37 h across
   0–30 mM Cu (fold-change at 30 mM vs 0 mM = 0.99x; mixed model F = 10.85,
   p ~ 5e-12 but negligible effect size). Kendall tau between peak-slope
   `rate_area` and doubling time = 0.002 — the two estimators are orthogonal.
2. **The previously reported "rate increases with Cu" is a phase artifact.**
   Peak-slope rate (6 h log-slope) rises with Cu because at low Cu most
   colonies reach a plateau before the window (saturated 95.4% at 0 mM vs
   73.9% at 30 mM), so their realized log-slope is lower; at high Cu colonies
   keep expanding through the ~120 h window. The real copper phenotype is
   extent-limited: lower final max area, lower plateau-reachability (t50 drops
   72.6 -> 64.1 h).
3. **Species/strains differ in copper sensitivity (extent-based).**
   `log2 extent ratio = log2(median max-area 25-30 mM / 0-5 mM)`:
   - Most tolerant: `R. glutinis` (+1.01), `R. pacifica`, `R. evergladensis`,
     and among well-sampled taxa **`R. mucilaginosa` is the most tolerant**
     (-0.72, n = 217 strains).
   - Most sensitive: `R. kratochvilovae` (-4.13), `R. paludigenum` (-3.94),
     `R. araucariae` (-3.68), `R. taiwanensis` (-3.06); among well-sampled,
     **`R. diobovata`** (-2.00).
   - Cu x species interaction on log(max_area) significant (F = 4.35,
     p ~ 5.9e-4).

## Caveats

- Species labels contain apparent near-duplicates (paludigenum/paludigena,
  evergladensis/evergladiensis) left unmerged; single-strain, low-n taxa
  (kratochvilovae, glutinis) need confirmation against per-strain tables.
- High-Cu missingness (~15-18%) is itself part of the sensitivity signal.
