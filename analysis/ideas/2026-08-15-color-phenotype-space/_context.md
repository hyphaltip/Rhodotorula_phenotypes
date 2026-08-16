# Project Context — Rhodotorula Copper Phenotype Screen (2026-08-15)

This context is given to every persona-subagent in this ideation session. Read it
carefully and ground every idea in the concrete data/analyses described here.

## What the project is

A high-throughput colony time-lapse imaging screen of **~320 yeast strains
(3 genera: Rhodotorula, Cystobasidium, Pseudomicrostroma)** grown on
**YPDN agar at 30 °C** in a 0–30 mM copper (Cu) dose gradient (7 dose levels:
0,5,10,15,20,25,30 mM; 16 replicate plates per dose; 112 experimental + 8
control plates across 5 imager runs). Colonies are segmented per image and a
~200-column feature vector is extracted per colony (CellProfiler-style).

## Data structure (DuckDB `v_phenotype` view — the common ABUNDANT layer)

One row = **one colony object at one image (timepoint)**. Per row:
- **Keying**: strain_id, unit strain_code, species (also genus derivable),
  run_number, plate_number, well_position, object_label, image_name,
  hours_since_plate_start (0–117 h), copper_mm (0–30), imaged_at.
- **Morphology (~20 cols)**: Shape_Area (px), Perimeter, Circularity,
  ConvexArea, Median/Mean/MaxRadius, Min/MaxFeretDiameter, Eccentricity,
  Solidity, Extent, BboxArea, Major/MinorAxisLength, Compactness, Orientation.
- **Intensity (~10 cols)**: IntegratedIntensity, Min, Max, Mean, Median, StdDev,
  Lower/UpperQuartile, InterquartileRange, Density, ConvexDensity.
- **ColorLab (~24 cols)**: L*, a*, b* — each with Min, Q1, Mean, Median, Q3,
  Max, StdDev, CoeffVar (INTRACOLONY pixel statistics, not across-colony).
- **ColorHSV (~24 cols)**: Hue, Saturation, Brightness — same 8 stats each.
- **Colorxy (~16 cols)**: chromaticity x, y — same 8 stats each.
- **ColorLab chroma**: ChromaEstimatedMean, ChromaEstimatedMedian.
- **TextureGray (~52 cols)**: Haralick features (AngularSecondMoment, Contrast,
  Correlation, HaralickVariance, InverseDifferenceMoment, SumAverage,
  SumVariance, SumEntropy, Entropy, DiffVariance, DiffEntropy, InfoCorrelation1/2)
  at 4 angles (000/045/090/135) + angle-averaged, scale-05.

Also available: per-strain metadata (`strain`): origin (China 180, Italy 119,
America 9) and environment (marsh_tidalflat 162, soil 45, plant 24, food 24,
unknown 12, marine 4, cave 3, air_cloud 2, insect 2, sand 2, snow_ice 2, rock 2,
built env 1, water 1).

## Existing analyses and key findings (DO NOT re-propose these)

1. **explore-plate-position** — within-Cu plate & grid-position variance is
   small (L* ~0.1–4.7%) but plate identity is nested in Cu; a naive pooled
   model over-estimated plate variance. Neighbor adjacency signal only for L*.
2. **growth-rates** — data-derived peak 6 h slope of log(area) & of consumer-
   light intensity as "rate". Peak-slope RISES with Cu (apparent toxicity
   paradox). Doubling-time test showed exponential doubling time is flat
   (~37 h across 0–30 mM): the rising peak slope is a PHASE ARTIFACT. Cu's
   real effect is extent-limited: saturation fraction 95.4% → 73.9%, t50
   72.6 → 64.1 h. L* is ~perfectly collinear with Intensity (>0.99).
   rate_area vs rate_int only r = 0.21 (two distinct growth axes).
3. **species Cu sensitivity** — extent-based log2 sensitivity index
   (25–30 vs 0–5 mM max area): R. mucilaginosa most tolerant (−0.72),
   R. kratochvilovae (−4.13) / araucariae (−3.68) / taiwanensis (−3.06)
   most sensitive; significant Cu × species interaction on log(max_area)
   (F = 4.34, p ~ 6e-4).
4. **control-late-timepoint phenotype** — per-strain tables (3 windows:
   70–80, 80–90, 90–110 h) of colony size + per-colony L*/a*/b* Median →
   median/mean/var/sd across replicate colonies on Cu = 0 only.

Essentially ALL existing analyses use only: Shape_Area, per-colony
ColorLab {L*,a*,b*}.Median (endpoint or time-course), Intensity_MeanIntensity.
**The other ~170 feature columns and their within-colony heterogeneity
(StdDev/CoeffVar), the HSV/xy/chroma descriptors, the full TextureGray set,
most morphology shape descriptors, and most of the temporal dynamics are
UNUSED so far.**

## Open directions the user cares about

"Phenotypic data space" and "color space" of this dataset:
- Which features/axes best describe strain diversity (beyond L*a*b* medians)?
- Intra-colony heterogeneity as a phenotype (StdDev/CoeffVar), surface/edge vs
  core color (Min/Q1 vs Q3/Max of the pixel histogram).
- Temporal trajectories of color (not just endpoint) — darkening, pigment
  appearance time, hue shifts.
- Structure of the 3-D+ color space: hue vs chroma vs lightness; whether
  strains cluster into pigment "morphs"; relationship to species, origin,
  environment, or Cu tolerance.
- Robustness/redundancy: how much of the ~200-col feature space is redundant
  (Info theory), and what the few independent axes are.
- Possible production/BRET relevance: Rhodotorula spp. are industrially
  important carotenoid (beta-carotene/torulene/astaxanthin) producers.

Note the caveat: L* collinear with intensity; a* (green↔red) and b*
(blue↔yellow) carry chromatic signal; hue is circular (careful with arithmetic
means).
