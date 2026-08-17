# Next-Generation GWAS Design — Color, Growth & Copper-Response, Applying Expert Review 08b

**Analysis**: `2026-08-15-color-phenotype-space`
**Design owner**: (this document); Expert priors: `08b-quantitative-geneticist-review.md` (agent QG review)
**Scope**: GWAS of R. mucilaginosa (201 genotyped strains) for (i) colony **color** at control
(YPD, Cu=0), (ii) **colony size growth rate**, and (iii) **copper-response** traits —
AUC / resistance / sensitivity / resilience to increasing copper. **Reuses the existing
pruned SNP panel + GEMMA platform** already built for the 4 original traits.

---

## 0. Reusable platform (already in hand — do NOT rebuild)

| Asset | Location | Notes |
|-------|----------|-------|
| Full scan matrix | `/scratch/jstajich/27384933/gwas/work/gwas.{bed,bim}` | **404,706** biallelic SNP × **201** Rmuc strains, integers 1–23 chr |
| LD-pruned kinship variants | `gwas.pruned.{bed,bim,fam}` | **20,769** variants (`--indep-pairwise 50 5 0.2`) |
| Kinship matrix | `output/kins.cXX.txt` | 201×201, already computed |
| GEMMA valid model | `gemma -lmm 4 -k kins` (no covariates) | produces `output/gwas_<t>.assoc.txt` |
| fam-phenotype convention | per-trait `.fam` col-6 | trait baked into col-6, strain order must match `gwas.fam` |

**Corollary (from 08b decision D-9)**: kinship-only LMM is the valid primary. Explicit PC /
population covariates are structurally degenerate (singular GRM, 22 near-clones) — do NOT add
them as fixed covariates.

---

## 1. Phenotype engineering (highest-leverage work, in order)

### 1.1 Clone-mean phenotyping — apply to ALL traits (08b rec #1)
- Currently each trait is one value per strain. The 16 replicate plates/dose × replicate
  colonies give **between-genotype vs within-genotype** variance.
- For each new trait: build per-`(strain, replicate_plate)` values at the reference window,
  then feed **clone-mean** (`mean` over replicates) to GEMMA; weight residual by `1/n_rep`.
- Expected power gain: single-isolate rel≈0.5 → clone-mean-of-3 rel≈0.75 → detectable per-SNP
  R² ~0.16 → ~0.11 (08b §1).
- Where replicates are unbalanced, retain strains with ≥2 replicates at the reference window.

### 1.2 (i) Colony color at control (YPD, Cu=0) — the new color-GWAS arm
- **Reference window**: fixed late time (100–110 h) at **`Copper concentration = 0`** only.
- **Candidate color axes** (from `v_phenotype` / `colony_measurement`, all present):
  - `ColorLab_ChromaEstimatedMean/Median`
  - `ColorHSV_SaturationMean/Median`, `ColorHSV_BrightnessMean/Median`
  - (intra-colony spread as secondary, 08b-derived: `*_StdDev/CoeffVar` — treat separately)
- Derivation per arm: per-colony → mean over colonies in the window → **clone-mean** over
  replicate plates → 1 value per strain, matched to `gwas.fam` order.
- Run an **ICC audit** first (08b design 04, Part A) to confirm each color axis has
  strain-level repeatability before GWAS; report ICC alongside.
- **Primary color trait**: `ColorLab_ChromaEstimatedMean` (control chroma = "color range/vividness").
  Secondary: `SaturationMean`, `BrightnessMean`; multi-trait-Fisher across the 3 (rec #5).

### 1.3 (ii) Colony size growth rate
- Already have `l10med_fixed` (colony size) — reuse as-is (clone-mean it).
- Optional firmer rate phenotype: linear slope of `log10(Shape_Area)` over the exponential
  phase from the `v_growth_timeseries` per strain at Cu=0; correlates with the existing trait.

### 1.4 (iii) Copper-response phenotypes (AUC / resistance / sensitivity / resilience)
All from `v_growth_timeseries` (`factors['Copper concentration']`) × 306 strains × 0–117 h.
Restrict to the 201 genotyped strains; clone-mean where replicates allow.

Per strain per dose:
- **AUC** = integral of growth curve over time: `∫ Shape_Area dh`, approx by trapezoid over the
  0–117 h timepoints (mirrors `idea_08` onset machinery). This is a *robust total-growth* proxy
  (robust to peak-slope phase, per 04 design).
- **Endpoint biomass** `max_area` at reference time.

Per strain dose-response (across the 7 doses 0–30 mM):
- **Sensitivity** = slope of `log(max_area)` (or `log(AUC)`) vs mM Cu (a negative slope =
  inhibition with dose). This is the growth analogue of the existing `slope_logchroma_per_mM`.
- **Resistance / tolerance** = **IC50_est** = interpolated Cu mM at which growth drops to 50% of
  control (safe linear-interp bracketing, robust-analysis: require monotone flanking points; where
  non-monotone, flag not-NA).
  - Alternatively a distribution-free index: **AUC_ratio** = `AUC(Cu)/AUC(0)` at a reference dose
    (e.g. 10 or 20 mM) — measures residual growth / tolerance.
- **Resilience** = ability to grow at the highest dose: `AUC(30)/AUC(0)` (retained growth under max stress).
- **Overall response to increasing Cu** = slope of the AUC-dose-response curve (negative →
  sensitive; ~0 → resilient). This is the direct "map resistance/sensitivity to condition" trait.
- Record Ceiling/floor clipping: at high Cu colonies may be below detection (04 warn, detection
  curve degrades >20–25 mM); restrict robust response inference to **Cu ≤ 20 mM** for the slope/IC50
  arm, and report AUC at max dose separately (measurement caveat).

**pH note (your caveat)**: Cu²⁺ speciation/activity varies with pH, so "response to increasing
Cu" may be partly pH-coupled. We have no per-plate pH measurement in `v_phenotype`; state this as
a covariate-caveat in the report and treat copper-response as a phenotype under the buffered plate
condition, not pure Cu²⁺ activity, unless pH data is added.

---

## 2. GWAS analysis tiers (apply per phenotype)

### Tier A — single-SNP marginal (baseline, reuse existing GEMMA)
```
gemma -lmm 4 -k output/kins.cXX.txt -bfile gwas_<t>  # <t>.fam col-6 = clone-mean phenotype
```
- **LOCO kinship** (08b rec #2, crucial for 23 scaffolds): build per-scaffold GRM
  (leave-one-chromosome-out) so a SNP's own-scaffold kinship doesn't absorb its signal.
  GEMMA: recompute kinship excluding each scaffold (gcta-style per-scaffold `--make-grm`,
  or GEMMA `-lmm 4` with per-scaffold GRM). This fixes proximal contamination.

### Tier B — set/gene-based tests (PRIMARY per 08b rec #3)
- From the `-lmm 4` per-SNP p-values, run **SKAT-O / burden / VEGAS**.
- **Set universe = top-dxy/π windows from pixy** (rec: focus on divergence, ~5× fewer tests).
- A locus of 5 same-direction variants reaches genome-wide signal at per-SNP R² 0.02–0.04 — this
  is where power actually is at n≈178–201.

### Tier C — architecture (rec #7)
- **BSLMM** (`gemma -bslmm/-ebv`): PVE/PGE, sparse-vs-polygenic → interpret whether color /
  copper-response is polygenic (predicts no genome-wide single-hit) or sparse.

### Threshold (rec #4) — honest, not nominal
- Compute **Meff** (Li & Ji) on the 20,769-pruned or on the LD structure of the scan set, OR
  **`--pheno-permute`** GEMMA permutation null. Report **nominal + empirical / FDR q-values**.

> **Meff proviso (2026-08-16)**: naive PCA-based Meff on n=201 is **bounded by sample size
> (~197) and is NOT the independent-test count** — it reflects sample relatedness, not SNP LD.
> Meff must be computed from the SNP–SNP LD correlation (M×M), which is intractable at
> n≪M by SVD. Practical bound per expert review (08b): **~5k–50k effective tests** on this
> near-clonal panel; nominal **5e-8 is already conservative**. Use `--pheno-permute` (or
> FDR on top-panged) as the defensible empirical threshold rather than a faulty scalar Meff.
- Expect effective independent tests 5k–50k → top p ~1e-5 may be significant → rescues negatives.

### Multi-trait combination (rec #5)
- Fisher / min-P / TATES across the color axes (chroma, sat, bright) and across the
  copper-response axes (AUC, IC50, AUC_ratio, dose-slope).
- Optional bivariate via REML (`gcta --reml-multivariate`) for the correlated color+growth pair.

### Lambda diagnostic (rec #9)
- For each trait, confirm kinship-only lambda (expect <1 from structure, cf. 0.36–0.64 observed).
- **Permuted-phenotype null** on a subset: if permuted≈observed (<1) → model-driven over-
  correction (structure); if permuted≈1 while observed<1 → partly real structure. This defends the
  negative.

---

## 3. Strain/resampling robustness (validity, rec #6)
- **Cull to ~178 max-informative**: greedy / `gcta --grm-subset` drop IBD≥0.95 duplicates to ONE;
  keep near-clones as §1.1 replicate phenotypes. Run both all-201 and culled-178, compare top-hit
  rank stability to catch clone-driven artifacts.

---

## 4. Negative-result write-up plan (if no genome-wide SNPs)
- Report a **power curve** at n_eff=178 (80% power at per-SNP R² 0.25–0.30 @5e-8; ~0.11 at
  set/SKAT level), nominal + empirical FDR, top-hit effect sizes + 95% CI + rank stability.
- Use **pixy dxy/fst co-localization**: ask if trait divergence (color / copper-response strain
  means) co-localizes with high-dxy windows — a region-level architecture claim fully supported
  without genome-wide SNPs. Check replication of the prior growth-rate locus `chr13:13_30149`.

---

## 5. Proposed execution order (concrete)

| # | Task | Tail | Uses |
|---|------|------|------|
| 1 | ICC audit of control color axes (chroma/sat/brightness/stdev) + clone-mean derivations | 04/08b | `v_phenotype` Cu=0 |
| 2 | Build copper-response phenotypes: AUC, dose-slope, IC50_est, AUC(20)/AUC(0), AUC(30)/AUC(0) for 201 strains | this doc §1.4 | `v_growth_timeseries` |
| 3 | Cull to ~178 informative strains; build clone-mean `gwas_<t>.fam` (col-6, matching order) | §3 | `gwas.fam` |
| 4 | LOCO kinship per scaffold | §2A | `gwas.pruned.*` |
| 5 | Tier A single-SNP GEMMA (kinship only, LOCO) for color+copper traits | §2A | GEMMA |
| 6 | Tier B SKAT-O/burden on pixy high-dxy windows | §2B | pixy genome results |
| 7 | Tier C BSLMM architecture + Meff/permute threshold | §2C/D | GEMMA bslmm |
| 8 | Multi-trait Fisher/TATES + pixy co-localization + replication | §2E/§4 | combined |
| 9 | Write `GWAS_REPORT.md` Stage 6 + power curve + `.living` updates | — | — |

---

## 6. Confirmed decisions (2026-08-16)

1. **Copper-response**: derive **ALL** copper traits (AUC, dose-slope, IC50_est, AUC_ratio,
   AUC(30)/AUC(0)) and run **multi-trait Fisher/TATES** across the correlated copper-response block.
2. **Color arm**: primary `ColorLab_ChromaEstimatedMean` + secondary `HSV SaturationMean` and
   `HSV BrightnessMean`, combined multi-trait. (Intra-colony spread StdDev/CoeffVar NOT in the
   color GWAS arm for now.)
3. **Strains**: run **BOTH all-201 AND culled-178**, compare top-hit rank stability.
4. **pH**: **ignore for now** — no per-plate pH exists; do not hard-code pH coupling into the
   report at this stage (deferred).

## 6b. (superseded by §6) Decisions to confirm with you
1. **Copper-response primary trait** — do you want (a) **dose-slope** (one number/trait, like existing
   color-slope), (b) **IC50_est**, (c) **AUC/AUC_ratio** at a reference dose, or (d) all, reported as a
   correlated block? (§1.4)
2. **Color arm scope** — primary `ChromaEstimatedMean` + secondary saturation/brightness, or also the
   intra-colony spread (StdDev/CoeffVar) axes? (§1.2)
3. **Strain culling** — run both all-201 AND culled-178 (recommended), or commit to culled-178 only 
   for speed? (§3)
4. **pH caveat** — confirm no per-plate pH is available before I hard-code the buffer-condition caveat
   into the report. (§1.4)