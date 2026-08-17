# Quantitative Geneticist Review — GWAS Power & Population-Structure Strategy

**Consultation**: 2026-08-16. Reviewed by expert quantitative geneticist (agent).
**Subject**: R. mucilaginosa GWAS (201 strains, 404,706 variants, LMM + 10 PC covariates) and
pixy diversity runs, given near-clone continuum + small clonal populations → ~178 effective
haplotypes, observed LMM lambda 0.36–0.64.

## Executive summary

Design is mechanistically correct but **bleeds power in three places**:
1. near-clone redundancy raises effective multiple-testing burden without adding info;
2. single-variant scans cannot reach the required R² at n≈178–201;
3. LMM with **10 PCs over-corrects** → conservative (deflated) lambda hiding real signal.

Lambda ≤0.64 is a **diagnostic of over-fitting structure, not a biology verdict**. The
highest-leverage change: **phenotype on clone-mean** to shrink environmental noise, then run
**set-based (SKAT/gene) tests** as primary instead of a single-SNP marginal scan.

## Key recommendations (ranked by power gain, feasibility, FP risk)

### 1. Phenotype on clone-mean (HIGHEST leverage, unused)
- Convert the 22 near-clone strains from a structure liability into **replicate phenotype
  reliability**: split total variance into between-genotype vs within-genotype (measurement+plasticity).
- reliability = σ²_g / (σ²_g + σ²_res/m). If σ²_res≈σ²_g: single isolate rel=0.5, clone-mean of
  3 replicates rel=0.75 → detectable per-SNP R² drops ~0.16→~0.11, allele shift 0.71→~0.58 SD.
- Feed clone-mean (~178 rows, not 201) to GEMMA; weight by 1/n_rep in residual covariance.

### 2. Cut PCs to ≤3 (or 0) + LOCO kinship
- 10 PCs are **redundant with kinship** (top PCs ≈ leading eigenvectors of G) → absorb df and
  variance, drive lambda to 0. Defensible compromise: **2–3 PCs**, or none (rely on kinship).
- **LOCO** (leave-one-chromosome-out) kinship crucial here — a SNP's own scaffold kinship can
  over-absorb its own signal (proximal contamination, 23 scaffolds). GEMMA needs per-scaffold
  GRM (gcta `--loocv-loco` or per-scaffold make-grm).
- DAPC rejected: discriminant method for *assigning* demes, not GWAS covariate (bias w/ small
  clonal groups, circularity). Use only as sanity check of 6-group assignment.

### 3. Set/gene-based tests as primary (4–8× locus power)
- SKAT-O / burden / VEGAS (from existing `-lmm 4` p-values). A locus with 5 variants same-direction
  needs only per-SNP R² 0.02–0.04 to reach same joint signal.
- Restrict set universe to **top-dxy/π windows** (pixy) → ~5× fewer tests, focus on divergence.

### 4. Honest threshold via Meff/permutation (not nominal 5e-8)
- Near-clone LD booms → effective independent tests far below 404,706 (plausibly 5k–50k).
- Meff (Li&Ji) or `--pheno-permute`: top p ~1e-5 could be significant under ~20k effective tests →
  most likely way to rescue "negative" → positive.

### 5. Bivariate + multi-trait combination
- chroma + cu_slope bivariate (gcta `--reml-multivariate`), + Fisher/min-P/TATES across the 4
  single-trait scans. Correlated color/growth biology → borrow strength, ~10–40% power gain.

### 6. Cull to ~178 maximally-informative strains (validity, not just power)
- Greedy subset / gcta `--grm-subset`, drop IBD≥0.95 pairs to one; keep clones as §1 replicates.
- Both-vs-one test: top rank stability all-201 vs culled-178 catches clone-driven artifacts.

### 7. BSLMM for architecture (interpretation power)
- GEMMA `-bslmm/-ebv` → PVE/PGE, sparse-vs-polygenic. If PGE small + PVE moderate → polygenic
  (predicts no genome-wide hits at this power curve). If PVE≈0 → environmentally/plasticity-dominated.

### 8. Frequency-aware rare-variant handling
- Keep MAC 5–20 in scan; SKAT weights freq-adaptively; frequency-stratified hit table.

### 9. Lambda diagnostics
- Rerun with 0/2/5/10 PCs + **permuted-phenotype null**. Permuted≈observed (both <1) →
  model-driven over-correction; permuted≈1 while observed<1 → partly real structure.

## Negative-result write-up (if few/no genome-wide hits)

Report a **power curve** (n_eff=178: 80% power at per-SNP R²≈0.25–0.30 @5e-8; ≈0.11 at set/SKAT
level). Use Meff/permutation threshold, report nominal + empirical-FDR q-values, top-hit effect
sizes + 95% CI + rank stability. Use pixy dxy/fst to ask if **trait divergence co-localizes with
high-dxy windows** (region-level architecture claim, fully supported without a genome-wide SNP).
Check replication of prior growth-rate locus `chr13 13_30149` at the set level.

## Bottom line

Not under-powered to answer "what is the genetic architecture"; under-powered to name the single
causal SNP (no reanalysis at 20 Mb will force that). Three concrete fixes: clone-mean phenotype,
≤3 PCs + LOCO, set-based + permutation threshold. The negative result is reportable as a
quantitative *architecture* finding with pixy + BSLMM evidence.