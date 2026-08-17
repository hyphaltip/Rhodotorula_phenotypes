# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### [2026-08-15] L-1 — mycelium scripts require Python >= 3.10
- **Category**: tooling
- **What happened**: `install_convention.py` (and other mycelium core scripts) use PEP 604 union syntax (`Path | None`), which throws `TypeError: unsupported operand type(s) for |` under the system default `python3` (3.9.18 on this host).
- **Why it matters**: Home-page python3 is too old; running the scripts fails loudly but misleadingly during mycelium init.
- **Resolution**: Run with `/usr/bin/python3.12` (available) or under the repo pixi env; also pass `--network-dir` explicitly because the script's auto-detect may miss the plugin marketplace path.
- **mitigation_type**: structural — add a version guard/`sys.version_info` check or a shebang in the scripts.
- **structural_mitigation_candidate**: assert `sys.version_info >= (3,10)` with a clear message at script start.
- **Tags**: tooling, mycelium, python, setup

### [2026-08-15] L-2 — skillpacks/Autonomous-Science repo is no longer public
- **Category**: tooling
- **What happened**: `git clone https://github.com/arjunrajlaboratory/Autonomous-Science.git` fails (`Repository not found`) — the repo is absent from the org's public repos. Only `scientific-agent-skills` and `bioSkills` cloned successfully.
- **Why it matters**: The `skill-bridge` convention's `skill-sources.yaml` lists this as a verified persona source (`personas/library/arjun_raj.json`); the persona routing feature is therefore currently unavailable.
- **Resolution**: Left `skillpacks/` with the two working repos. TODO: re-point `skill-bridge` `skill-sources.yaml` to a live personas source or drop the source entry.
- **mitigation_type**: ambient-awareness — monitor whether the org restores/renames the repo.
- **Tags**: tooling, skill-bridge, skillpacks, third-party

### [2026-08-15] L-3 — dplyr mutate using `named_vec[[.$col[1]]]` broadcasts row 1's value to all rows
- **Category**: R / data-wrangling
- **What happened**: In `05_time_variance_partition.R` I relabeled a bound list's `.id` column with `mutate(trait = traits[[.$trait[1]]])`. Because `.[1]` takes the *first element of the whole tibble*, every row got the first trait's label ("L* (lightness)") instead of its own — 6 distinct model fits collapsed to 6 identical labels, and the `pivot_wider` in validation then threw "values not uniquely identified".
- **Why it matters**: This exact anti-pattern (indexing a lookup vector with `[[...]]` + `[1]` inside a single mutate) silently labels every row with the first group's name; it only surfaces later as duplicate keys. Recovery was easy here because the models had already been fit correctly.
- **Resolution**: Use vectorized name lookup `unname(traits[.$trait])` so each row maps through `traits[name]`. Validate by checking `n_distinct(trait)` or that a `pivot_wider` on the relabeled key doesn't warn.
- **mitigation_type**: structural — prefer `lookup[.$col]` over `lookup[[.$col[1]]]` when relabeling within a mutate.
- **Tags**: R, dplyr, mutate, pivot_wider, label-map, debugging

### [2026-08-15] L-4 — dplyr `summarise(across(...))` with `.names` throws "subscript out of bounds" when a group's lambda returns length-0
- **Category**: R / data-wrangling
- **What happened**: `summarise(across(c(a,b,c), ~ .x[which(!is.na(.x))[...]], .names = "last_{.col}"))` crashed with `Error in names(dots)[[i]] : subscript out of bounds` when a plate/well had an all-NA column (empty index vector). Repro was minimal and isolated to the `across`-with-`.names` path.
- **Why it matters**: The "take last non-NA per group" idiom is common for endpoints; the failure is opaque and can't be debugged from the message.
- **Resolution**: Replace with explicit per-column `summarise(x_late = last_non_na(x), ...)` using a small helper `last_non_na <- function(x){y <- x[!is.na(x)]; if(length(y)==0) NA else y[length(y)]}`.
- **mitigation_type**: structural — avoid variadic `across` when the lambda can return length-0 inside `summarise`.
- **Tags**: R, dplyr, summarise, across, last, endpoint, gotcha

### [2026-08-15] L-5 — lmerTest `difflsmeans()` needs the fitting environment, not a list wrapper; realize such "diff vs ref" tables from fixef/vcov instead
- **Category**: R / mixed-models
- **What happened**: Stashed lmer fits in a list and passed the *list* to `difflsmeans()`, which then failed with `no applicable method` (class "list"); passing the model object itself failed with `object 'copper_mm' not found` because the fitting data.frame was a local that difflsmeans re-evaluates via model.frame.
- **Why it matters**: Post-hoc contrast helpers in lmerTest/lme4 are brittle about env/data retention; the failure wastes a full model rerun cycle.
- **Resolution**: Since `factor(copper_mm)` coefficients are already differences vs the reference level, compute contrasts directly as `fixef(m)` + `sqrt(diag(vcov(m)))` Wald 95% CI (no emmeans needed — the pixi env lacks emmeans/AICcmodavg/minpack.lm).
- **Tags**: R, lme4, lmerTest, difflsmeans, contrasts, mixed-model, gotcha

### [2026-08-15] L-6 — rmarkdown chunk CWD = directory of the .Rmd, not the render() call; pdflatex chokes on Unicode ≥/≤/→/⇒/−/…
- **Category**: R / reporting
- **What happened**: Merged growth-rates report failed reading `"analysis/growth_rates/results/tables/..."` because knitr executes chunks with working directory set to the document's folder (so paths must be `results/tables`). Separately, pdflatex compiled em-dashes/× fine but aborted on U+2265/≤/→/⇒/−/…; xelatex in the shared TeX Live 2022 failed on a stale l3names kernel instead.
- **Why it matters**: Two silent portability traps cost several render cycles. The repo's existing plate-position report already follows the document-relative path convention.
- **Resolution**: In Rmds under this repo use doc-relative `results/tables` and keep Unicode to `—` and `×` (what the plate-position report proves pdflatex accepts); else switch the engine. Detect remaining non-ASCII with a quick python one-liner printing U+ codepoints before rendering.
- **mitigation_type**: structural — a report template/lint for allowed Unicode + doc-relative paths.
- **Tags**: R, rmarkdown, knitr, working-directory, pdflatex, unicode, reporting

### [2026-08-15] L-7 — a metadata-CSV label fix does NOT propagate to DuckDB: `strain` used `ON CONFLICT DO NOTHING` and full re-import is blocked by FK from measurements to imager_run
- **Category**: data / DB import
- **What happened**: Corrected two species-name misspellings in `data/metadata/Copper.Strain_info.csv`, then re-ran `10_import_experiment.py`: (1) the global `strain` upsert was `DO NOTHING`, so existing rows kept the old labels; (2) a full re-import crashed with `Constraint error: key run_number 353 is still referenced by a foreign key` — measurements reference `imager_run.run_number`, so the imager_run upsert cannot mid-flight replace rows already referenced.
- **Why it matters**: Analyses join species from the DB `strain` table, not from the CSV, so the CSV fix silently had no effect on any downstream output until the DB was updated. And the only "official" propagation path (full re-import) is structurally impossible in place once measurements load.
- **Resolution**: Switched the `strain` upsert to `DO UPDATE SET` the metadata columns (CSV authoritative) and added `--strain-only` to `10_import_experiment.py` to sync just the strain table without touching imager_run/well_placement; made the previously-required experiment/plate args optional with a fail-loud guard under `--strain-only`. Running `pixi run python scripts/db/10_import_experiment.py --strain-only --strain-info-csv data/metadata/Copper.Strain_info.csv` fixed strains 84 and 327.
- **mitigation_type**: structural — any metadata edit workflow should use `--strain-only`; the full import is only safe on a fresh DB.
- **Tags**: data, duckdb, import, metadata, strain, FK, copper

### [2026-08-15] L-8 — choosing a "timepoint" in this DB: bin hours_since_plate_start to the imaging pass, do NOT take per-strain max hour
- **Category**: data / phenotype extraction
- **What happened**: For a control-only late-timepoint strain phenotype table, taking each strain's exact maximum `hours_since_plate_start` (<=110 h) to define its "latest timepoint" yielded only 1-2 colonies/strain for 84 strains: within a single imaging pass, colonies of one strain image a few minutes apart, so an exact-max filter split one pass into fragments. Rounding to the nearest integer hour and taking the max *rounded* pass restored the full pass (median 4 colonies/strain).
- **Why it matters**: The imaging rig interleaves two ~3 h cadences; a fixed single hour covers only ~231/320 strains, while pass-of-latest-image-in-window covers 314/320 (286 with >=3 colonies). Any per-strain/plate timepoint pick should use the rounded-hour pass.
- **Resolution**: `with c0 as (select *, round(hours_since_plate_start,0) tp_h ...), latest as (select *, max(tp_h) over (partition by strain_id) tp_hmax from c0) ... where tp_h = tp_hmax` — i.e. round, then window-max the rounded hour.
- **mitigation_type**: structural — a reusable timepoint-binning helper / convention for time-series extractions on this pipeline data.
- **Tags**: data, duckdb, timepoint, imaging-pass, cadence, phenotype, gotcha, copper

### [2026-08-15] L-9 — user-site Python packages shadow the pixi conda env unless PYTHONNOUSERSITE=1
- **Category**: environment / pixi
- **What happened**: `numba` refused to import inside the pixi env with `TypeError: numpy.core.multiarray failed to import` because the conda env's numpy 2.4.x was being shadowed by a newer `~/.local/lib/python3.12/site-packages` numpy (2.5.2) on `sys.path`.
- **Why it matters**: User-site packages are a permanent source of version skew on this host; even after pinning `numpy = ">=2.4,<2.5"` in `pixi.toml`, the shadowing can still occur depending on `PYTHONUSERBASE`/`.pth` behavior.
- **Resolution**: Add to `pixi.toml` `[activation.env] PYTHONNOUSERSITE = "1"` so the env never reads the user-site directory; verify with `python -c "import numpy, sys; print(sys.path)"`.
- **mitigation_type**: structural — always set PYTHONNOUSERSITE in repo pixi.toml activation env.
- **Tags**: pixi, environment, numpy, numba, user-site, activation-env, gotcha

### [2026-08-15] L-10 — statsmodels `anova_lm` cannot compare multiple fitted models; compute the F-test by hand
- **Category**: statistics / pandas / statsmodels
- **What happened**: Passing several model refits to `anova_lm(m1, m2, ...)` for type-II tests of a strain×Cu interaction failed (multi-model comparisons not supported in this api).
- **Why it matters**: Hypothesis tests for GxE / added variables are a routine need; blocking on the statsmodels API was wasting time.
- **Resolution**: Fit the reduced and full model separately and compute `F = ((RSS_reduced - RSS_full) / (df_reduced - df_full)) / (RSS_full / df_full)`, p from `scipy.stats.f.sf`.
- **mitigation_type**: procedural — keep a small helper for hand-rolled nested-model F-tests.
- **Tags**: statsmodels, anova, F-test, nested-model, GxE, stats, gotcha

### [2026-08-15] L-11 — basal chroma noise (~1.5-2) in late Cu=0 colonies forces a minimum onset pigment threshold
- **Category**: analysis / phenotype
- **What happened**: With onset defined as first time chroma crosses a low threshold (2-3), nearly every colony was "onset at t=0" (`t_darkening_h=0` everywhere), and median onset collapsed to 0-12 h; threshold sweep showed median onset jumps 0/12/30/48/84 h as threshold rises 2→3→5→7→10.
- **Why it matters**: Baseline colony color has intrinsic chroma spread (~1.5-2) unrelated to pigment onset; any "time-to-pigment" derived from low thresholds reports imaging noise, not biology.
- **Resolution**: Use chroma threshold >=7 (`ColourLab`-scale chroma) for onset definitions; report the sweep as a robustness figure. Only ~51% of colonies reach chroma >= 7 under the IM scope window, so onset-definition must be a documented parameter.
- **mitigation_type**: analysis, procedural — threshold-sweep before committing to an onset definition.
- **Tags**: onset, chroma, threshold, pigment, timecourse, imaging, noise, gotcha

### [2026-08-15] L-12 — median one-time: naive residual-ON-residual mediation blows up numerically
- **Category**: statistics / modeling
- **What happened**: A rough "regress predictor on mediator residuals, then fit outcome on those residuals" mediation attempt produced wildly unstable / non-interpretable coefficients (collinearity of constructed regressors).
- **Why it matters**: Stats courses illustrate mediation with regressions of residuals, but those models are numerically pathological with correlated phenotype measures; the resulting "effects" are not meaningful.
- **Resolution**: Use a proper structural decomposition: outcome ~ species + area + Cu (total) versus outcome ~ species + Cu (direct effect with the mediator in the model), i.e. compare coefficients of Cu between models with/without the mediator; bootstrap the a*b indirect component.
- **mitigation_type**: analysis, procedural — use product-of-coefficients / model-difference decomposition, never residual-residual hacks.
- **Tags**: mediation, collinearity, regression, stats, bootstrap, gotcha

### [2026-08-15] L-13 — sd(log₁₀ colony area) is a size/floor-coupled metric; pool across strains and Cu and you will fabricate a "dispersion grows with stress" trend
- **Category**: analysis / statistics
- **What happened**: Idea 01's Gibrat check found sd(log₁₀ Asat) across replicate wells rising with Cu (0.30→0.43; pooled Spearman 0.23–0.28, p<1e-27). A within-strain follow-up showed raw paired slopes are indeed positive (+0.0033 per mM, Wilcoxon p≈5e-16) but that sd(log₁₀ Asat) is anti-correlated with colony size (Spearman −0.63, p<1e-200): Cu shrinks colonies, and small/floor-adjacent colonies carry inflated log-dispersion. Once colony size was included as a covariate, the within-strain Cu slope collapsed to ~0 (+0.0002, p=0.46; 53.6% positive).
- **Why it matters**: Any "heterogeneity / disorder / dispersion vs treatment" claim built on log-size (or log-intensity, log-chroma) dispersion across replicates will be confounded by the confound of treatment moving the mean near the detection floor. This subsumes the idea 01b diversity:higher-Cu => more log-variance => uncertainty story — it was an artifact.
- **Resolution**: For dispersion-type statistics, always (1) compute within-unit slopes (paired design), (2) control for the metric's own level (median size) as a covariate, and (3) report the pooled trend only alongside the size-controlled one. Here this revealed genuine lineage exceptions (R. taiwanensis widens, R. paludigena narrows) rather than a universal effect.
- **mitigation_type**: analysis, procedural — dispersion metrics need within-unit + size-controlled design, never naive pooled trends.
- **Tags**: dispersion, gibrat, noise, size-artifact, log-variance, variance, stats, floor-effect, gotcha

### [2026-08-15] L-14 — Blomberg's K collapses to ~0 on "near-comet" trees (giant polytomy of near-duplicate genomes); likelihood lambda becomes numerically unstable there too
- **Category**: phylogenetics / statistics
- **What happened**: The PHYling protein tree for the 278 Rhodotorula strains is a near-comet: a huge polytomy of near-duplicate R. mucilaginosa genomes (167/541 edges ≤1e-7; 22 zero-length pendant edges; patristic-distance min ≈7.6e-9) with 216+/278 tips mucilaginosa. On this tree Blomberg's K returned ≈1e-7 for every trait, and a BM power check reproduced K≈2.25 (should be ~1) on simulated data on the same tree — K is degenerate, not "no signal". Pagel's lambda via geiger `fitContinuous` was also internally inconsistent (white-noise lnL exceeded lambda lnL); phytools `phylosig` lambda agreed directionally (λ≈0.83–0.94) with the true result.
- **Why it matters**: On near-comet/genome-cluster trees, reporting Blomberg's K as "no phylogenetic signal" is a false conclusion driven by topology. Likelihood-based lambda is fragile. This tree shape is common for routine strain panels dominated by one clade.
- **Resolution**: Use a rank-based permutation test on pairwise distances — Mantel (Spearman) between cophenetic/patristic distance and absolute trait distance, 999 permutations (seed 7). This gave a coherent, significant signal (size r=0.43, heterogeneity r=0.31, baseline chroma r=0.19, Cu slope r=0.10; all p≤0.003) where K said ~0. Also report lambda as reference only and disclose the topology caveat.
- **mitigation_type**: analysis, procedural — for phylogenetic signal on near-comet trees use Mantel permutation, never raw K; sanity-check K via a BM power simulation on the actual tree.
- **Tags**: phylogenetics, Blomberg-K, pagel-lambda, mantel, permutation, polytomy, near-comet, apc, ggtree, phylosig, gotcha

### [2026-08-15] L-15 — Existing SNP panels may already cover your strains: check shared Population_Genomics projects before proposing de-novo variant calling
- **Category**: data reuse / variant analysis
- **What happened**: Idea 11 (GWAS power) started from the assumption we'd need to call variants from the 278 scaffold assemblies in `BFD/input/dna/`. That turned out to be completely unnecessary: a finished population-genomics + GWAS project already exists at `/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/` (reference NRRL Y-2510, GATK-pipeline) with `vcf/RmucY2510_v2.All.SNP.combined_selected.vcf.gz` = 728,581 SNP sites x 422 haploid strains. **201 of our 278 phenotyped strains** (every R. mucilaginosa in our panel) carry genotypes, and 200/201 have complete data for all 4 GWAS-ready traits. Their GWAS framework (GEMMA 0.98.3 + plink2 + bcftools in `GWAS/vc2gwas_env/bin`) and even a 218-strain growth-rate GWAS (T4C-T37C/Salt6) are already run; we only need to swap in our color traits.
- **Why it matters**: De-novo SNP calling (bwa→GVCF→joint-call→hard-filter) is days of compute. The single-lab strain collection makes genotype overlap across projects highly likely; locating the existing VCF took minutes. Also: their strain names match ours exactly (DBVPG_/TFCN_/EXF_/CHIFNET_/F6_/F8_), which made the join trivial.
- **Resolution**: Before any variant-calling plan, (1) list `shared/projects/Population_Genomics/`, (2) grep the project's `samples.csv`/VCF headers against our phenotype strain list, (3) verify haploid GT format (R. mucilaginosa is haploid — no 0/1 calls, no HWE filtering needed) and non-missingness, then reuse. Power may be limited by effective n (178 for mucilaginosa) regardless — de-novo calling would only have added sites, not more independent strains.
- **mitigation_type**: workflow, procedural — always survey shared population-genomics projects for existing SNP panels + prior GWAS before building variant-calling pipelines; verify strain-name overlap and ploidy.
- **Tags**: gwas, snp, variant, data-reuse, ploidy, haploid, population-genomics, NRRL-Y2510, power, effective-n, idea11, gotcha

### [2026-08-16] L-16 — GEMMA reads phenotype from fam column 6, silently IGNORES `-p` when `-bfile` is used; also needs integer chromosome codes and `-p` even for `-gk`
- **Category**: tooling, gotcha
- **What happened**: Running `gemma -bfile gwas -k kins -lmm 4 -p pheno.txt` gave beta=se=0 for every SNP with p_wald=NaN — because GEMMA uses the FAM column-6 phenotype exclusively and silently ignores `-p` when `-bfile` is the input. Our fam col-6 was all `-9` (missing). Two more traps: (1) GEMMA requires integer chromosome codes — scaffold names ("scaffold_10") are treated as non-autosomal → "number of analyzed individuals = 0"; fixed by renaming `.bim` per the source project's `genome/Chrom_Mapping.tab` (scaffold_1..23 → 1..23); (2) even the kinship step (`-gk 1`) needs a `-p` phenotype or it also reports 0 individuals.
- **Why it matters**: A GWAS that appears to run and "succeed" can silently produce all-zero results. All three traps manifest as plausible-but-wrong output, not a crash.
- **Resolution**: Bake the trait into FAM column 6 (`awk '{print $1,$2,$3,$4,$5,$VALUE}'`), drop `-p` entirely for the association run; rename `.bim` chromosomes to integers; always pass `-p` (any file) to `-gk`. Verify with a tiny linear run (`-lm 1`) that betas are non-zero before trusting the LMM.
- **mitigation_type**: procedural — audit GEMMA input files (fam col-6, .bim chrom codes) before running; sanity-check that a trivial run gives nonzero betas.
- **Tags**: gwas, gemma, LMM, phenotype, fam, chromosome-code, ploidy, gotcha, idea11

### [2026-08-16] L-17 — Adding structure covariates (PCs or population dummies) to a GEMMA LMM collapses when the kinship is singular from near-clone redundancy
- **Category**: statistical, gotcha
- **What happened**: GEMMA `-lmm 4 -k kins -c <10 genotype PCs>` ran ~55 min/trait and produced `pve estimate = 0.99997`, `se(pve) = -nan`, plus `GSL ERROR: matrix is singular in lu.c`. Retrying with 3 PCs gave the identical degeneration. A kinship-only model (no `-c`) is fine: PVE = 0.180 ± 0.077, sane lambdas 0.36–0.64, real hits. The cause is structural: with 22 near-clone strains the kinship (GRM) is singular/ill-conditioned, so ANY fixed covariate (PCs or population dummies — the dummy run threw GSL singular too) pushes the joint model over the edge; top genotype PCs are also nearly collinear with kinship eigenvectors.
- **Why it matters**: Population-structure correction is standard practice, but in near-clonal/haploid panels it can break the model rather than help. PC count was not the fix (3 vs 10 both failed).
- **Resolution**: Rely on kinship-only correction for the LMM scan in this panel; do NOT add explicit PC/population covariates. If covariate correction is mandatory, use Lochin (LOCO kinship) or collapse near-clone strains to a culled set first. Diagnostics: null PVE stuck at ≈1 with NaN se is the signature.
- **mitigation_type**: analytical — test a zero-/kinship-only model for a sane PVE before adding covariates; watch for PVE→1 + NaN se + GSL unitary/singular errors.
- **Tags**: gwas, gemma, LMM, kinship, PCA, population-structure, near-clone, singular, PVE, gotcha, idea11

### [2026-08-16] L-18 — UCR HPCC batch jobs don't inherit login shell: bgzip/tabix/pixi unavailable by default; and submitting from a node-local `$SCRATCH` cwd kills the job (signal 53)
- **Category**: HPC, gotcha
- **What happened**: (1) A GATK+pixy SLURM job died at `bgzip: command not found` — the login shell's `module load bcftools` environment is NOT inherited by batch jobs, which start with a bare PATH. (2) A trivial `sbatch --wrap="hostname"` failed instantly with `RaisedSignal:53` whenever submitted from a cwd inside my interactive session's node-local `/scratch/jstajich/<jobid>/...` — SLURM sets the batch job's WorkDir to that cwd, which does not exist on the target compute node. (3) `$SLURM_JOB_SCRATCH` does not exist on UCR HPCC; the node-local scratch env var is `$SCRATCH`.
- **Why it matters**: Batch-only tools silently vanish, and a missing-WorkDir submit looks like a random system failure (signal 53) rather than a path bug — very costly to chase.
- **Resolution**: In SLURM scripts, `source /etc/profile.d/modules.sh; module load <pkg>` at the top; `sbatch -D /path/on/shared/fs` from a shared (not node-local scratch) workdir; use `$SCRATCH` (not `$SLURM_JOB_SCRATCH`) for node-local scratch. Verify availability with a tiny env-probe job before a long run.
- **mitigation_type**: procedural/structural — standardize a job header snippet (module loading, -D, $SCRATCH) for all cluster jobs.
- **Tags**: HPC, slurm, scratch, modules, bgzip, tabix, pixy, gotcha, upload, workflow

### [2026-08-16] L-19 — GEMMA 0.98.3 has no single `-loco` flag; LOCO = per-scaffold GRM + repeated runs
- **What happened**: Designing the LOCO sensitivity follow-up for next-gen GWAS (D-10). GEMMA 0.98.3 (from the reused NRRL-Y2510 `vc2gwas_env`) does not expose a built-in leave-one-chromosome-out option via a flag; `-help`/usage shows no `loco`/`chromosome` switch. The only supported relatedness input is a precomputed `-k matrix`, so LOCO must be constructed manually.
- **Why it matters**: LOCO is the key sensitivity safeguard against a SNP's own scaffold being over-absorbed by the kinship matrix (proximal contamination). Assuming a magic flag exists would silently produce the wrong (full-kinship) analysis.
- **Resolution**: Implement LOCO by (a) partitioning the LD-pruned variant set by chromosome (integer codes already in `gwas.bim`), (b) running `gemma -gk 1` on each scaffold-excluded subset to make a per-scaffold GRM, and (c) running `gemma -lmm 4` on scaffold `s`'s SNPs with the GRM built from all OTHER scaffolds. This is 23 GRM builds + per-scaffold scans per trait — cheap to restrict to top trait(s) first.
- **mitigation_type**: technical/procedural — verify tool capability (help/usage) before planning a flag-based approach; build LOCO as a wrapper loop.
- **Tags**: gemma, gwas, LOCO, kinship, chromosome, gotcha, nextgen

### [2026-08-16] L-20 — GEMMA 0.98.3 BSLMM invocation + output format; and PCA-Meff is bounded by n, not a test count
- **What happened**: Running Tier-C BSLMM (`gemma -h 10` → `-bslmm` option 1 = BSLMM). Outputs are `output/bslmm_<o>.hyp.txt` (columns: `h pve rho pge pi n_gamma`, one row per MCMC sample) and `bslmm_<o>.gamma.txt` (one row per MCMC sample, 300 columns = GEMMA's internal LD clusters; cell value = representative variant index (1-based into the scan `.bim`), 0 = cluster inactive). Naively averaging the gamma cells as "inclusion probability" gives nonsense (values are SNP indices, not 0/1); the correct cluster PIP = fraction of post-burn-in samples with nonzero cell. Also, computing "Meff" via PCA/SVD of genotypes on n=201 yields ~197 — bounded by the number of samples, NOT the number of independent tests; Li&Ji Meff needs the SNP×SNP LD correlation (M×M), intractable here (n≪M).
- **Why it matters**: BSLMM is the clean way to call architecture (sparse vs polygenic) on a near-clonal panel: it groups the 404k SNPs into ~300 LD clusters, so n_gamma is interpretable (AUC_10: n_gamma=3, PGE=0.96 → near-oligogenic). Misreading gamma as 0/1 inclusion or misquoting Meff would corrupt the significance story.
- **Resolution**: PIP = mean(active) per cluster after dropping burn-in; map rep index via the trait `.bim` (same SNP order across `ngwas_*` bfiles). For significance use a GEMMA `--pheno-permute` permutation null or FDR — not a scalar Meff on n=201 (which is a sample-size bound, not a test-count bound). Nominal 5e-8 is already conservative on ~5k–50k effective tests.
- **mitigation_type**: analytical — always validate the shape/semantics of a new tool's output (check header/dims) before computing derived statistics; prefer permutation/FDR to a mis-scaled Meff.
- **Tags**: gemma, bslmm, BSLMM, PVE, PGE, PIP, gamma, LD-cluster, meff, permutation, gotcha, nextgen

### [2026-08-16] L-21 — GEMMA `-lmm 4` maxT permutation is unusable on near-clonal panels; use FDR instead
- **What happened**: Attempted Tier-C completion via an empirical genome-wide threshold from phenotype-permutation (shuffle trait in fam col-6, rerun `-lmm 4 -k fixedk kinship`, record min-p per perm; 1000 perms × chroma/AUC_10/resilience_30 on the LD-pruned set). The null min-p distribution was catastrophically uncalibrated: min-p reached **1e-17..1e-18 under permutation** (expected ~1/(0.05*n)=~2e-6 for well-calibrated). Median per-SNP p stayed ~0.5 (bulk correct) but the extreme tail collapsed — a few SNPs got absurd p because, in this near-clonal panel, permuting phenotypes while holding a fixed kinship lets variance-component outliers deflate residual variance. GEMMA 0.98.3 has no `--pheno-permute` (only added in newer builds).
- **Why it matters**: A maxT/empirical-threshold from such a null would reject every real hit (chroma's real 2.4e-8 looks weaker than the broken null's 1e-17) — a false-negative disaster.
- **Resolution**: Use **Benjamini-Hochberg FDR (q=0.05) on the existing full-scan Tier-A p-values** — robust, no new compute, explicitly endorsed by the expert review (rec #4) and `09-NEXT-GWAS-DESIGN.md` line 119. The permutation effort still served its intended diagnostic role: permuted-λ vs observed-λ comparison (design line 129). Also: **FDR doubles as an inflation detector** — the inflated-λ traits (clone_mean_area λ=2.12, AUC_0 λ=2.74) behave pathologically under FDR (AUC_0 → 7340 "significant", clone_mean_area → 0), marking them untrustworthy.
- **mitigation_type**: analytical — verify a permutation null's tail calibration (median should stay uniform, expected min-p ~ 0.05/M) before trusting a maxT threshold; prefer FDR/multi-trait on highly structured/near-clonal panels.
- **Tags**: gwas, gemma, permutation, maxT, FDR, benjamini-hochberg, inflation, lambda, near-clone, threshold, gotcha, nextgen

### [2026-08-17] L-22 — `np.linalg.eigvalsh` raises LinAlgError on near-singular LD matrices; guard with small-shift retry + scipy fallback
- **What happened**: Tier-B MC-verify eigendecomposition of residual window LD matrices (subset after dropping SNPs with missing z) hit `numpy.linalg.LinAlgError: Eigenvalues did not converge` at runtime (gwasc run, job 27511892) — the matrix can be numerically near-singular/poorly conditioned at the intersection subset.
- **Why it matters**: The three-moment SKAT approximation AND the exact MC path both need the eigenvalues of the LD matrix; a crash on one window kills the whole verification step after hours of compute.
- **Resolution**: `_safe_eigvals()` in `scripts/tierb_set_tests.py` — try plain `eigvalsh`, then `eigvalsh` with tiny ridge shifts (1e-8, 1e-6 · I), then `scipy.linalg.eigvalsh(driver='evr')`; drop eigenvalues ≤ 1e-12. Also: SKAT moment-approx p-values were verified against exact Monte Carlo (50k draws) on the top-50 windows — r=0.984 in log10 space, worst drift 0.47 log10 (~3×) on a p≈4e-3 window; the fast approximation is fine for ranking/no-call decisions.
- **mitigation_type**: computational — wrap eigendecompositions of real-data covariance/LD matrices in a shift-retry+fallback helper; validate closed-form moments vs MC on a subset.
- **Tags**: gwas, tierb, skat, ld, eigendecomposition, eigen, numpy, scipy, near-singular, monte-carlo, gotcha

### [2026-08-17] L-23 — figure gotchas: same-SNP p for identity plots; never log-scale a -log10 axis; show CIs on null-result bars
- **What happened**: Dataviz review of the three GWAS follow-up figures found two correctness bugs and two presentation gaps. (1) The LOCO-vs-Tier-A identity plot plotted x = -log10(trait-level max Tier-A p) but y = the anchor SNP's LOCO p — for gwasc-chroma the trait max is a *different SNP* (scaffold_3_570085, the culled-set shift) so the point sat far off the identity line making LOCO look like it lost the signal it actually retained. (2) The Tier-B window-rank panel used `ax.set_yscale("log")` on an axis already holding -log10 p (a log-of-a-log), putting the p=0.05 reference at an unreadable position. (3) The co-localization panel B (observed-vs-genome-wide rate) had NO error bars — inappropriate for a headline null result — and the 80th-pct high-dxy/high-Fst thresholds were in the text but not on the scatter.
- **Why it matters**: Both "bugs" would mislead a reader clinically (apparent LOCO failure; unreadable threshold) and panel B without CIs understates the no-enrichment finding the figure exists to prove.
- **Resolution**: identity plots must compare p at the **same SNP on both axes** (use the merged per-chromosome `top_rs_TierA`/`top_p_TierA`); keep -log10 axes linear and draw threshold bands; put 80th-pct threshold lines on the scatter; add Wilson 95% CIs + Fisher OR/p in-figure (recomputed p matched `coloc_enrichment.txt` exactly). Also: overlay the exact-MC points on the approx-SKAT panel to make the validation claim visible, and de-collide anchor labels that share a window (chroma/AUC_10 both in scaffold_10).
- **mitigation_type**: visualization — for any p-p identity/scatter, verify both axes refer to the same locus; never transform an already-logged axis; every bar of a null-result figure gets a CI; pre-commit a figure review against its own headline claim.
- **Tags**: gwas, visualization, matplotlib, dataviz, identity-plot, log10, wilson-ci, fisher, figure, gotcha

### [2026-08-17] L-24 — Tier D/E gotchas: FDR mult-test columns; positional clump hides 2nd block/chr; BSLMM CSV unquoted commas; ABF must be in z-space; nearest-proxy replication
- **What happened**: Executing Tier D/G/E surfaced four concrete gotchas. (1) GEMMA `-lma` output carries multiple test p-value columns (Wald p, score test p, LRT p); the FDR table and the assoc CSV differ in which p is "the" p (AUC_10 chr13:30134 is p_wald=4.03e-6 vs a second-test p=7.90e-6 — both tiny, but pick one column per file and be explicit). (2) The Tier-D "independent loci" clump selects the single most-significant SNP per chromosome, so a chromosome (or small scaffold) with a second, distant block collapses — chr13's 217 FDR-sig AUC_10 SNPs (spanning 11,701–800,664, chained by <250 kb steps) become ONE lead (13_791853); the prior locus 13_30134 (AF 0.015) and everything else in the block vanish from `tierD_independent_loci.csv` even though FDR-sig. (3) The Tier-C BSLMM summary CSV is malformed (unquoted commas inside the top_PIP field) → `pd.read_csv` silently mis-parses; use `csv.reader` or parse by column count. (4) Wakefield ABF posterior over trait-scale betas is meaningless when traits differ in units/scale by >1e5 (AUC_10 β~8e5 vs chroma β~1): compute in scale-free z-space with the prior on the non-centrality parameter (SD=0.2); and without a candidate filter (p<1e-3) credible sets span thousands of null SNPs.
- **Why it matters**: Each silently changed the headline: last-column p differences, vanished replication loci, mis-parsed PIPs, and nonsense credible sets would each have produced a wrong "gene/replication" conclusion.
- **Resolution**: Pick ONE explicit p column per file (p_wald) and document it; for independent loci use window-based per-trait top-loci (per-chromosome second block allowed) when scaffolds hold multiple blocks; parse BSLMM CSV defensively; ABF in z-space + candidate filter; test replication via the nearest-genotyped proxy (exact pos MAF-filtered) and state the offset (15 bp).
- **mitigation_type**: statistical/parsing — cite the exact p column used; never rely on a per-chromosome single-lead clump to enumerate all blocks; defensively parse pandas-incompatible CSVs; fine-map in z-space; replication = nearest-proxy with stated genomic offset.
- **Tags**: gwas, fdr, pvalue, columns, clump, independent-loci, bslmm, csv, parsing, abf, finemapping, z-space, credible-sets, replication, proxy, gotcha
