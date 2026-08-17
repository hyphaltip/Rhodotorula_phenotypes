# GWAS + Pixy Live Progress Log

**Project**: Rhodotorula phenotypes (analysis `2026-08-15-color-phenotype-space`)
**Started**: 2026-08-16 ~08:30
**Working dir**: `/scratch/jstajich/27384933/gwas/work` (results copied back to `results/gwas/`)
**Status**: pixy pilot done; **full-genome (23-scaffold) pixy COMPLETE** (job 27502734, 6h35m); next-gen Tier-A GEMMA (24 scans) + **Tier-C BSLMM (5 traits) DONE**; manhattan/QQ figures (png+pdf) generated; **Post-Tier-A follow-up §13 COMPLETE** (LOCO sensitivity + Tier B set tests + dxy/Fst co-localization); remaining: dataviz consult on new figures

> Update this log after every step. Record: steps run, steps tested (incl. failed paths),
> results found, wall-clock runtimes, and next actions. Companion docs: `GWAS_REPORT.md`
> (methods), `04-quantitative-geneticist.md` (design), `FINDINGS.md` (science).

---

## 1. GWAS pipeline — SNP matrix (IDE11 reusing NRRL-Y2510 panel)

| Step | Action | Result | Runtime |
|------|--------|--------|---------|
| S0 | `bcftools view -S our200.txt` → `sub.our201.vcf.gz` | 201 samples × 728,581 sites | ~3–4 min (detached) |
| S1 | exclude `exclude.bed` (scaffold_21 whole + scaffold_6/11/18 high-depth windows) | 728,581 → **628,841** sites | ~1 min |
| S2 | `bcftools view -m2 -M2` biallelic | 628,841 (already biallelic); haploid verified (0 het, 0 missing) | <1 min |
| S3 | `plink2 --geno 0.1 --mind 0.1 --mac 5 --make-bed` | **201 × 404,706** variants; 0 samples dropped | ~1 min |
| S3b | LD-prune kinship `--indep-pairwise 50 5 0.2` | **20,769** pruned variants (`gwas.pruned.*`) | <1 min |

**GEMMA gotchas (FIXED, learning-worthy):**
1. Needs integer chromosome codes — renamed `.bim` per `genome/Chrom_Mapping.tab` (scaffold_1..23 → 1..23). Without this → "number of analyzed individuals = 0".
2. Needs `-p` even for `-gk 1` kinship (0 analyzed individuals otherwise).
3. **`-bfile` + `-p` phenotype IGNORED — GEMMA reads fam col-6** (was all `-9` → beta=se=0 for all SNPs). FIXED by baking trait into `.fam` col-6 and dropping `-p`. Verified via linear run (404,681 nonzero betas, top scaffold_2 p=1.4e-10).

## 2. GWAS — GEMMA LMM (uncorrected results)

| Step | Action | Result | Runtime |
|------|--------|--------|---------|
| S4a | `gemma -gk 1` kinship (pruned 20,769) | `output/kins.cXX.txt` (201×201) | ~1 min |
| S4b | `gemma -lmm 4 -k kins` per trait (size, chroma, cu_slope, pace) | `output/gwas_<t>.assoc.txt` | ~3 min/trait |

**Uncorrected summary** (via `summarize.py` + PLINK clump `--clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.5 --clump-kb 250`):
- Lambda: **0.357–0.638** (all <1 → conservative, over-corrected by near-clone relatedness).
- Raw→independent loci: **size** 44→13, **cu_slope** 147→16, **pace** 14→13.
- Top size hit: `scaffold_16:121473` p_wald=**6.2e-11**; size hits cluster on scaffold_16 (32/44).

## 3. GWAS robustness — population-stratification correction

| Step | Action | Result | Runtime |
|------|--------|--------|---------|
| T1 | Discrete 6-pop dummy covariates (`cov_pops.txt`) | **GSL "matrix is singular" error** — near-clones related AND same-pop → redundant with kinship | <1 min (FAILED — abandoned) |
| T2 | DAPC considered for correction | **Rejected** — supervised, needs labels, circularity risk, poor for near-clone continuum | decision |
| T3 | PCA of LD-pruned genotypes `plink2 --pca 10` | `pca.eigenvec/eigenval`; top PCs 32.4/13.7/11.6/8.6…% ≈ ~100% variance (near-clone-dominated); `cov_pcs.txt` (IID + PC1–10) | <1 min |
| T4 | `gemma -lmm 4 -k kins -c cov_pcs.txt` (10 PCs) | size/chroma DONE; **LUInvert + PVE=0.99997, SE=NaN** (degenerate); ABORTED per expert | ~55 min/trait (too slow + broken) |
| T5 | `gemma -lmm 4 -k kins -c cov_pcs3.txt` (3 PCs) | size done; SAME degeneration (PVE=0.99997, NaN SE, singular). 3-PC does NOT fix it → parked | ~8 min/trait (faster but still broken) |

**Critical finding (T4/T5)**: BOTH the 10-PC AND 3-PC covariate LMMs are numerically degenerate
(pve=0.99997, se(pve)=NaN, `GSL ERROR: matrix is singular in lu.c`). This is **structural, not
PC-count-dependent**: the kinship matrix is singular/ill-conditioned from the 22 near-clone strains,
so adding ANY fixed covariate collapses GEMMA. Confirms (empirically) the earlier GSL singular error
with discrete population dummies, and the expert review's over-correction thesis.
**Conclusion**: the **kinship-only LMM (`gwas_<t>`, no PCs) is the valid analysis** — it ran cleanly
(lambda 0.357–0.638, top size hit scaffold_16:121473 p=6.2e-11; **null PVE = 0.180 ± 0.077**, a sane
heritability estimate). PC-covariate variants are NOT reportable results; they document that explicit
structure covariates are unsupported by this near-clonal relatedness (singular design). 3-PC run
**stopped** after size (was degenerate); only `output/gwasc3_size.assoc.txt` kept for the record.
**User decision 2026-08-16**: abort 10-PC run; retried 3-PC → also degenerate → parked PC-covariate
approach. Forward path per expert: clone-mean phenotyping, set/gene-based tests, LOCO kinship,
Meff/permutation threshold.

## 4. Pixy pilot (scaffold_1 only) — pi/theta/dxy/fst

Needs VCF **with invariant sites** for unbiased π/θ/dxy/fst; variant-only SNP VCF insufficient →
joint-genotype per-sample haploid GVCFs (from `.../NRRLY2510/gvcf/`, 424 strains; all 201 of ours covered).

| Step | Action | Result | Runtime |
|------|--------|--------|---------|
| P0 | env: `pixi add pixy` | pixy 2.2.3 | - |
| P1 | `gvcf_201.list`, `pixy_pops.txt` (6 pops) | 201 samples | - |
| P2 | job 27498382 (stajichlab) | FAILED signal 53 — **submit cwd was node-local `$SCRATCH` scratch** → WorkDir missing on compute node | 0s |
| P3 | job 27498523 (epyc, 24G) | FAILED **OUT_OF_MEMORY** at 4% of scaffold_1 (201-sample GenomicsDBImport) | 12 min |
| P4 | job 27498706 (epyc, **96G**, `--batch-size 50`) | **GenotypeGVCFs SUCCEEDED** (2,175,785 sites, 26.1 min) but job FAILED at bgzip: `bgzip: command not found` (batch node, no htslib loaded); scratch wiped | 35 min (failed at step 3) |
| P5 | job 27499759 (epyc, 96G, `--batch-size 50`, **+`module load bcftools`**) | bgzip/tabix **FIXED** (passed); but **pixy step died**: `--stats ... theta` invalid → pixy 2.2.3 uses `watterson_theta`; also self-copy `cp` bug (cd-before-cp) | ~36 min (failed at pixy) |
| P6 | job 27502295 | **SUCCEEDED** — cohort VCF persisted (`cohort.scaffold_1.all.vcf.gz`, 403 MB); produced pi, dxy, watterson_theta. **fst missing**: pixy default Weir-Cockerham NOT supported for haploid contigs | ~45 min |
| P7 | job 27502392 (`run_pixy_fst.sh`, `--fst_type hudson`, pixy-only, reads persisted cohort) | **SUCCEEDED** — pi, dxy, **fst (Hudson)**, watterson_theta all produced, 49 s | 49 s |

**Pixy pilot RESULTS (scaffold_1, 6 pops, 22×100 kb windows)**:
- **pi highly heterogeneous across populations**: pop5 mean π=0.044, pop6=0.016, pop3=0.0088 vs pop1/pop4≈5–7e-5 (near-monomorphic). Strong architecture signal.
- **dxy** range 2.9e-5 – 0.053 (max ≈60× min).
- **fst (Hudson)**: overall mean **0.494**; least divergent pop3-pop6=0.251, pop5-pop6=0.300; most divergent pop4-pop1=0.773, pop4-pop2=0.876.
- Interpretation: some populations are near-clonal monomorphic (low π), others highly internally diverse; populations strongly differentiated. Matches near-clone continuum. Full genome will be run next (all 23 scaffolds).

**Also fixed**:
- `$SLURM_JOB_SCRATCH` doesn't exist on UCR HPCC → use `$SCRATCH`.
- input files (`our201.list`, `pixy_pops.txt`) copied to shared pixy results dir so compute node can read.
- **bgzip/tabix NOT on batch-node PATH** — login `module` env not inherited. Fix: `source /etc/profile.d/modules.sh; module load bcftools` (htslib) in script. Verified via env-probe job (27499755): bgzip/tabix/pixi all available on compute nodes after `module load bcftools`; `pixi` at `/rhome/jstajich/.pixi/bin/pixi` ON PATH.
- **pixy 2.2.3 stats keyword**: `theta` → `watterson_theta` (pi/dxy/fst/watterson_theta/tajima_d).
- **persist cohort VCF** to `$OUT` immediately after bgzip/tabix (step 3b) so GATK need not be rebuilt if pixy step fails. Subsequent pixy-only reruns read the persisted `$OUT/cohort.scaffold_1.all.vcf.gz`.

**Note (effort)**: each full pixy pass rebuilds the scaffold_1 cohort VCF from scratch (~30 min GATK + import) because node-local `$SCRATCH` is wiped when a job ends. Now that the VCF is persisted to `$OUT`, future runs are pixy-only (fast) unless the GATK output is deleted.

**Pipeline**: GATK 4.6.2.0 `GenomicsDBImport --consolidate --merge-input-intervals --batch-size 50`
→ `GenotypeGVCFs -all-sites` (keeps invariant sites) → `bgzip`+`tabix` → `pixi run pixy
--stats pi theta dxy fst --window_size 100000 --n_cores 8`.

## 5. Next actions

- [ ] Poll + record PC-LMM (`gemma_pc.done`) results; re-summarize + re-clump; compare vs uncorrected loci.
- [x] Full-genome pixy (see §8): genome_{pi,dxy,fst,watterson_theta}.txt in `results/gwas/pixy/`. REMAINING: dxy/fst trait co-localization.
- [x] Next-gen Tier-C BSLMM architecture (PVE/PGE) + PVE estim (see §7). Remaining: `--pheno-permute` permutation threshold for top traits.
- [ ] Next-gen Tier B: set/SKAT + burden on pixy high-dxy windows (from Tier-A p-values).
- [ ] **LOCO sensitivity follow-up** (D-10, `10-LOCO-RUNBOOK.md`): rerun top trait(s) with per-scaffold GRM to confirm Tier-A sensitivity; compare lambda/top-hit/rank stability.
- [ ] Copy final GWAS outputs (plots, CSVs) into repo `results/gwas/tierA_summary/` (done for Tier-A).
- [ ] Finalize `GWAS_REPORT.md`; add `.living` entries (GEMMA fam-phenotype gotcha, pixy invariant-VCF caveat, SLURM scratch WorkDir gotcha, OOM/mem tuning).

## 6. Next-gen Tier-A GEMMA (D-10; 12 phenotypes × all-201 + culled-173)

Phenotypes engineered (clone-mean over ~4 replicate plates, `build_gwas_phenotypes.py`):
color block `chroma/sat/bright`, size `clone_mean_area`, copper block
`AUC_0/10/20/30`, `AUC_ratio_10`, `resilience_30`, `cu_dose_slope`, `IC50_est`.

| Step | Action | Result | Runtime |
|------|--------|--------|---------|
| N1 | `build_gwas_phenotypes.py` (clone-mean over plates; fam-order aligned) | 12 traits × 201 strains, 201/201 aligned | ~5 min |
| N2 | IBS near-clone cull (IBS0<0.005 greedy) | **173 informative strains**, 28 dropped (`culled_keep.list`) | ~1 min |
| N3 | `gwasc` subset + pruned kinship `kinsc.cXX.txt` (GEMMA -gk) | 173×173 | ~1 min |
| N4 | stage per-trait bfiles (ngwas_<t>, ngwasc_<t>; trait in .fam col-6) | 24 bfiles | ~1 min |
| N5 | `gemma -lmm 4 -k kins` per trait (both sets) | **24 scans DONE** | ~3–4 min/scan (~90 min) |
| N6 | `summarize_tiera.py` → lambda/PVE/top hits + multi-trait Fisher (color/copper blocks) | `tierA_summary/` | ~5–10 min |

**Tier-A headline (all-201):**
- Color multi-trait Fisher **p=7.5e-16 @ scaffold_10:384905** (chroma/sat/bright co-localize).
- Copper block **p=5.9e-12 @ scaffold_16:421208**; AUC_30 `scaffold_4:508257` repeats across both sets.
- Most lambdas <1 (conservative, structure); clone_mean_area & AUC_0 λ≈2.1–2.7 (inflated, flagged).
- **Culled-173 shifts top color hit to scaffold_3:570085** → rank instability flagged for LOCO.
- IC50_est culled 673 SNPs<5e-8 at λ=0.19 → treated as structure inflation, not signal.

**Decision D-10 records**: Tier A → Tier B/C first, LOCO as committed sensitivity follow-up.

## 7. Next-gen Tier-C BSLMM (all-201, 5 representative traits)

GEMMA `-bslmm 1` (Bayesian sparse linear mixed model), 100k MCMC, 20% burn-in discard.
`run_ngwas_tierc.sh`; outputs in `results/gwas/tierC_summary/`.

| Trait | PVE [95%] | PGE | n_gamma | top-PIP loci | architecture |
|-------|-----------|-----|---------|--------------|--------------|
| chroma | 0.244 [0.12,0.42] | 0.39 | 138 | scaffold_3:26959/137161/132903/7440 (PIP 1.0) | moderate h², partly sparse (scaffold_3) |
| AUC_10 | 0.403 [0.27,0.53] | 0.96 | 3 | scaffold_5:533603 + scaffold_2:939277 (PIP 1.0) | **near-oligogenic (~3 variants)** |
| AUC_30 | 0.387 [0.22,0.61] | 0.61 | 31 | scaffold_16 block | sparse+polygenic |
| resilience_30 | 0.252 [0.10,0.45] | 0.62 | 4 | scaffold_16:563722 (0.95), scaffold_2:71223 (0.79) | sparse (scaffold_16) |
| clone_mean_area | 0.071 [0.01,0.23] | 0.52 | 17 | weak | low h² |

**Takeaways**: BSLMM converges with Tier A (resilience scaffold_16 ≈ multi-copper; chroma
scaffold_3 ≈ multi-color). **AUC_10 near-oligogenic** = strongest, most interpretable copper
signal (3 loci carry ~all h²). clone_mean_area weakly heritable (matches its Tier-A inflated λ).
Meff proviso: PCA-Meff on n=201 is bounded by sample size, NOT test count → use `--pheno-permute`
permutation null / FDR; nominal 5e-8 already conservative (~5k–50k effective tests).

## 8. Full-genome pixy COMPLETE (job 27502734, 6h35m, COMPLETED)

All 23 scaffolds joint-genotyped (invariant sites retained) → 1290×100 kb windows.

| Metric | n | mean | median | max |
|--------|---|---|--------|-----|
| π (avg_pi) | 1290 | 0.01162 | - | 0.0501 |
| Watterson θ | 1290 | 0.01039 | - | 0.0443 |
| dxy (avg_dxy) | 3225 | - | 0.01486 | 0.500 |
| Fst (Hudson) | 3223 | 0.4464 | 0.4386 | 0.9970 |

- Results in `results/gwas/pixy/genome_{pi,dxy,fst,watterson_theta}.txt` (all Aug 16 21:30).
- Top dxy windows on scaffold_20 (pop2-3 window_pos_1) and scaffold_19 (pos 200001) → candidate
   divergence hotspots for trait co-localization (Tier B set tests).
- Consistent with scaffold_1 pilot: strong genome-wide population structure (mean Fst 0.45).

## 9. Tier-C completion — significance threshold (FDR; permutation tried, rejected)

**Attempted**: GEMMA `--pheno-permute` does not exist in 0.98.3 → manual maxT null (permute
fam col-6, rerun `-lmm 4`, min-p per perm; 1000 × chroma/AUC_10/resilience_30 on LD-pruned
set). **REJECTED**: null min-p reached 1e-17..1e-18 (expected ~2e-6) — fixed-kinship +
permuted-phenotype variance outliers collapse the tail in this near-clonal panel. A maxT
threshold from it would reject every real hit (L-21, D-11).

**Adopted**: **Benjamini-Hochberg FDR (q=0.05)** on existing full-scan Tier-A p-values
(expert rec #4, design line 119). Output `output/fdr/tiera_fdr_summary.csv`:

| Trait | λ | n FDR05 | top |
|-------|-----|---------|-----|
| chroma | 0.32 | 345 | scaffold_10:384905 (2.4e-8) |
| sat | 0.39 | 40 | scaffold_1:218972 (1.8e-7) |
| bright | 0.34 | 0 | - |
| clone_mean_area | 2.12 | **0** | - (inflated λ → untrustworthy) |
| AUC_0 | 2.74 | **7340** | scaffold_8:1190725 (spurious; inflated λ) |
| AUC_10 | 0.73 | 3989 | scaffold_10:396172 (1.4e-8) |
| AUC_20 | 0.66 | 56 | scaffold_16:421208 |
| AUC_30 | 0.50 | 31 | scaffold_4:508257 |
| AUC_ratio_10 | 0.64 | 0 | - |
| resilience_30 | 0.45 | 524 | scaffold_13:810026 (6.4e-9) |
| cu_dose_slope | 0.35 | 22 | scaffold_3:546068 |
| IC50_est | 0.25 | 1 | scaffold_16:122361 (8.7e-12) |

FDR doubles as an inflation detector: the two inflated-λ traits behave pathologically under
it (AUC_0→7340 "sig", clone_mean_area→0). Validated Tier-C loci: chroma scaffold_10,
AUC_10 scaffold_10 scaffold, resilience scaffold_13/scaffold_16, plus moderate sat/
AUC_20/AUC_30/cu_dose_slope/IC50. Unchanged from before: BSLMM convergence (chroma
scaffold_3, resilience scaffold_16).

## Expert review (quantitative geneticist) — status

Launched 2026-08-16. Reviewing: power constraints (178 effective haplotypes, min per-SNP R²≈0.16),
near-clone genetic structure, PC-vs-DAPC correction choice, and strategies to retain power.
See companion output document when complete.

## 13. Post-Tier-A follow-up (D-13) — COMPLETE 2026-08-16/17

Three follow-ups from expert review, all executed and merged into `GWAS_REPORT.md` §8:

1. **LOCO sensitivity** (6 traits × 20 chrs, GEMMA `-lmm 4 -k` per-chromosome kinship).
   Array 27511xxx, 120/120 outputs, 0 failures. **All Tier-A anchors reproduce at unchanged p**
   (chroma 2.46e-8, AUC_10 1.44e-8, resilience_30 6.04e-9 all-chr; gwasc anchors 2.3e-5–5e-8).
   LOCO λ (medians 0.34–0.73) tracks Tier-A λ → the λ<1 deflation is stable near-clonal
   structure, NOT per-chromosome rescue by LOCO. Plots: `results/gwas/figures/loco_sensitivity.*`.
2. **Tier B set tests** (burden/SKAT/min-p on 178 high-dxy windows, gwas + culled gwasc sets).
   No set-level signal survives FDR(q<0.05). burden over-conservative (λ≈0.5, sign-cancellation),
   SKAT mildly deflated (λ≈0.6–0.8, top p~1.4e-3), min_p recovers Tier-A single-SNP hits.
   SKAT three-moment approx verified vs exact MC (gwas r=0.984 n=50; gwasc r=0.997 n=47,
   worst drift 0.871 log10 at p~1e-3 → still non-significant). Job 27511559 (gwas), 27511892 (gwasc, after L-22 eigdec fix).
3. **dxy/Fst co-localization**: FDR GWAS loci NOT enriched in high-divergence windows
   (dxy OR=0.81 p=0.61; Fst OR=0.41 p=0.065). Phenotypes segregate within the near-clonal focal
   clade (standing variation), decoupled from deep pop splits. `scaffold_20:100001` (dxy 0.070)
   flagged as mis-assembly artifact, excluded.

Outputs: `results/gwas/tierB/{tierb_settests_{gwas,gwasc},tierb_skat_mcver_{gwas,gwasc}}.csv`,
`results/gwas/tierB/coloc/*`, `results/gwas/figures/{tierB_settests,coloc_dxy_fst}.*`,
`results/gwas/loco/loco_merged_{gwas,gwasc}.csv`.