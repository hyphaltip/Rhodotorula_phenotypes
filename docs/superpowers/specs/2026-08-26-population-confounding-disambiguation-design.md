# Design: distinguishing real causal loci from population-structure artifacts

**Date**: 2026-08-26
**Status**: revised after expert review (fable statistical geneticist/breeder consult),
ready to implement

## Problem

`check_population_confounding.py` flags a Tier A top hit as "population-confound risk"
using a crude heuristic: max-minus-min alt-allele frequency across the 6 subpopulations
> 0.8, while the overall AF isn't itself near-fixed. By this heuristic, ALL of chroma,
lab_L, lab_a, lab_b, AUC_20, cu_dose_slope, and sat's top hits are flagged in at least
one panel (`analysis/gwas/results/gwas/tierA_summary/population_confounding_{gwas,gwasc}.csv`).
The kinship-only LMM (GEMMA `-lmm 4 -k`) corrects for genome-wide relatedness, but does
not guarantee a SINGLE locus's marginal association isn't primarily an echo of
between-population phenotype differences that happen to correlate with that locus's
population-differentiated allele frequency — population membership and phenotype could
covary for reasons that have nothing to do with that specific SNP being causal (shared
ancestry, environment, or drift at many loci simultaneously).

The AF-swing heuristic alone cannot distinguish:
(a) a locus that is *causal* and *also* happens to differ in frequency across
    populations (real biology, e.g. local adaptation — this is a legitimate finding,
    not an artifact), from
(b) a locus that is *not* causal but is *correlated with* population membership, which
    itself correlates with the phenotype for unrelated reasons (a true artifact).

## Proposed disambiguation battery

Six tests, each attacking the question from a different angle, applied to every flagged
(trait, panel, top SNP) combination:

1. **Within-population re-test** (the breeder's/QTL-mapper's gold standard). For each
   population with enough within-population allelic variance at the locus (require
   both allele classes present with n≥3 each), fit `phenotype ~ genotype_dosage`
   *separately within that population only* — this removes ALL between-population
   variance by construction. Report per-population beta, SE, p, n. A locus whose effect
   direction/magnitude is *consistent across ≥2 independent populations* (even if each
   is individually underpowered) is strong evidence of a real, population-independent
   effect. A locus that is null or sign-flips in every population, with the marginal
   signal only appearing when populations are pooled, is the artifact signature.

2. **Fixed-effect meta-analysis** of the per-population betas from (1) (inverse-variance
   weighted). Compare the meta-analytic p to the marginal kinship-LMM p. Also report
   Cochran's Q / I² for effect heterogeneity across populations — high heterogeneity
   with a population-specific direction flip is itself informative (could mean
   population-specific genetic background modifies the effect, or that only one
   population is driving a spurious pooled signal).

3. **Single-locus population-covariate regression.** `phenotype ~ genotype + population
   (6-level factor)`, ordinary least squares at just this one SNP (not genome-wide —
   avoids the GRM/PC collinearity that made D-9 reject genome-wide PC-covariate GEMMA
   runs; that collapse was a many-SNP, many-covariate problem, not a per-locus one).
   Report genotype's partial F-test / t p-value after population is already in the
   model, and partial R² attributable to genotype beyond population (vs. population's
   own R²).

4. **Continuous-ancestry version of (3)**: same regression but with GRM PC1/PC2 (already
   computed in §5 of GWAS.md) as continuous covariates instead of discrete population
   bins — softer, handles admixed/intermediate strains the 6-way split may misclassify.

5. **Fst-matched empirical null**: for each flagged SNP, compute its per-pair or global
   Fst (or reuse the AF-range statistic already computed). Sample a large (~5,000) set
   of random genome-wide SNPs matched on that same Fst/AF-differentiation band but
   otherwise unselected, and look up each one's GWAS p-value for the same trait (already
   computed, no new GEMMA runs needed — this reuses the existing Tier A assoc files).
   If the flagged SNP's p-value is not more extreme than this matched-null distribution
   (e.g. not in the matched-null's top 1%), it is statistically indistinguishable from
   "any similarly population-differentiated SNP would show an association this strong
   just from tracking population" — i.e., not special.

6. **Chromosome-level LOCO** (reuses `run_loco.sh` infrastructure) for the specific
   chromosome each flagged locus sits on, for traits not yet covered by the existing
   3-trait LOCO run (chroma, lab_L, lab_a, lab_b, AUC_20, cu_dose_slope, sat all need
   this — only chroma is a repeat, at a different locus per panel each time per §7).
   A signal that survives its own chromosome's kinship contribution being zeroed out is
   less likely to be a kinship-absorption artifact (same logic already validated for
   resilience_30/scaffold_13:810026 in §8).

## Deliverable

One script (`scripts/check_population_vs_locus.py` or similar) producing a table, one
row per flagged (trait, panel) combination, columns: marginal kinship-LMM p, within-pop
meta-analysis p + per-pop betas (list), partial-R²-beyond-population, PC-adjusted p,
Fst-matched-null percentile, LOCO p (chromosome-level). A verdict column categorizing
each as `likely_real` (survives within-pop test in ≥2 populations AND beats the
Fst-matched null) / `likely_population_artifact` (fails within-pop test AND/OR
indistinguishable from Fst-matched null) / `ambiguous_underpowered` (too few strains per
population per allele class for a within-pop test to be informative).

## Expert review (fable, statistical geneticist/breeder persona) — verdict: sound in
spirit, over-engineered as drafted; trim, fix the null, and add explicit relatedness
handling. Full response archived below; changes adopted:

1. **Drop test 4** (continuous PC-covariate regression) as redundant with test 3 — with
   only 6 well-differentiated, near-discrete lineages (Fst~0.45, not continuous
   admixture), the discrete population-factor regression is both sufficient and more
   interpretable. Formal admixture mapping / local ancestry inference is overkill for
   this population structure; not added. Phylogenetic GLS/TreeWAS is redundant with
   tests 1+3 given the GRM already encodes the same tree structure; not added.

2. **Test 5 (Fst-matched null) redone as a 2D MAF×Fst match**, not Fst alone — power to
   detect association depends on MAF independent of Fst, so an Fst-only match would pull
   systematically different-power comparator SNPs and bias the percentile. Also exclude
   an LD window around any OTHER significant hit for the same trait before sampling
   comparator SNPs, so the null isn't contaminated by tagging a real signal elsewhere.

3. **Test 1 (within-population re-test) is exactly as vulnerable to near-clone
   pseudoreplication as the genome-wide scan** — this was the sharpest flaw identified.
   Fixes adopted: (a) always run on the **182-strain near-clone-culled panel**, never the
   full 213; (b) as a required secondary layer, additionally collapse any REMAINING
   near-identical genotype clusters within each population (reusing the IBS0
   `cull_near_clones.py` pairwise-distance machinery, applied per-population) and refit —
   an association that survives clone-collapse is much stronger evidence than a naive
   per-population OLS; (c) report both `n_raw` and `n_effective` (post-clone-collapse)
   per population in the output table, so a deceptively small effective n is visible, not
   hidden. Full within-population GEMMA LMM was explicitly rejected as infeasible at
   n=18-77 per population (rank-deficient GRM, too few df) — plain OLS on clone-collapsed
   data is the adopted approach instead.

4. **Verdict rule replaced**: no unanimity requirement, no single combined statistic
   (both were rejected as unworkable given 4 partially non-independent tests and small
   per-population n). Adopted rule-based verdict:
   - `likely_real` = within-population directional consistency (same sign, comparable
     magnitude) in **≥2 populations** (clone-collapse-aware) **AND** at least 2 of the
     remaining 3 quantitative tests (meta-analysis significance, partial-R² of genotype
     beyond population, Fst×MAF-matched-null percentile) support it.
   - `likely_population_artifact` = fails the within-population directional-consistency
     test, **OR** partial-R² of genotype beyond population is ≈0.
   - everything else → `ambiguous_underpowered`.
   - Benjamini-Hochberg FDR applied specifically to the meta-analysis p-values (the one
     test in the battery with a clean null distribution suitable for FDR) — not across
     the whole battery.

5. **MAS-usability gates adopted for anything reaching `likely_real`** (beyond this
   disambiguation battery, before any marker-assisted-selection claim): a minimum
   effect-size threshold independent of p-value (locus must move the phenotype by some
   multiple of within-population SD, not just be statistically distinguishable from
   zero); a held-out replication check if any genuinely disjoint strain subset exists;
   an LD-decay check confirming the flagged SNP isn't a distant tag of a stronger nearby
   signal; and — specific to the color traits (L*/a*/b*/chroma) — a plate/batch-confound
   sanity check, since color measurement is exactly the kind of trait vulnerable to
   batch effects that could themselves correlate with which population was measured on
   which plate/day. These gates are OUT OF SCOPE for the initial disambiguation pass
   (which only answers "is this population structure or not") and are deferred to a
   follow-up once any locus actually reaches `likely_real`.

## Final test battery (post-review)

A. **Within-population re-test**, 182-strain culled panel + per-population clone-collapse
   (secondary layer), OLS `phenotype ~ genotype_dosage` per population with ≥3 strains
   per allele class post-collapse. Reports beta/SE/p/n_raw/n_effective per population.
B. **Fixed-effect inverse-variance meta-analysis** of (A)'s per-population betas, plus
   Cochran's Q/I² for heterogeneity. BH-FDR applied here across all flagged loci.
C. **Single-locus population-covariate regression**: `phenotype ~ genotype + population`
   (6-level factor, both panels), partial F-test / partial R² of genotype beyond
   population.
D. **Fst×MAF-matched empirical null**: ~5,000 genome-wide comparator SNPs matched on a
   2D (Fst quintile × MAF quintile) bin, LD-window-excluded around other same-trait hits,
   percentile rank of the flagged SNP's actual GWAS p-value within that matched null.

Chromosome-level LOCO (originally test 6) is retained as a **supplementary, non-gating**
check reported alongside the verdict (reusing existing `run_loco.sh` infrastructure for
the newly-flagged trait/chromosome combinations) — informative but not part of the formal
decision rule, consistent with the reviewer's verdict rule only naming A-D.
