# GWAS: Color and Copper-Response Phenotypes in *Rhodotorula mucilaginosa*

**Status**: in progress (port from `analysis/ideas/2026-08-15-color-phenotype-space/`) — Tier A (corrected, full marker set, both panels), Tier C (BSLMM), and LOCO complete; Tier B running; pixy reuse flagged with a caveat (not recomputed)
**Spec**: `docs/superpowers/specs/2026-08-25-gwas-port-design.md`
**Plan**: `docs/superpowers/plans/2026-08-25-gwas-port-implementation.md`

This analysis ports the Tier A-G GWAS pipeline (single-SNP kinship-only LMM, SKAT/burden
set tests, BSLMM architecture, LOCO sensitivity, gene annotation, fine-mapping, prior-locus
replication) originally developed in `analysis/ideas/2026-08-15-color-phenotype-space/`
(see that folder's `GWAS_REPORT.md` and `PROGRESS.md` for the full original narrative and
`.living/decisions.md` D-9 through D-14 for the modeling decisions carried forward
unchanged here). All new code for this port lives under `analysis/gwas/` — the ideas
folder is read from (its `db_extract.parquet`, `common.py`, `build_gwas_phenotypes.py`)
but never modified.

New in this port (not done rigorously in the original run):
1. Audited strain-name reconciliation (`scripts/reconcile_strains.py`).
2. Per-strain ploidy/heterozygosity validation (`scripts/check_ploidy.py`).

## 1. Strain reconciliation

`scripts/reconcile_strains.py` matched `data/metadata/Copper.Strain_info.csv`'s `Strain`
column (308 phenotyped strains) against the 422 VCF sample IDs in
`data/raw/genotypes/RmucY2510_v2/`, replacing the prior run's informal `our200.txt`/
`our201.txt` lists:

| tier | count |
|---|---|
| exact | 213 |
| normalized | 1 |
| fuzzy | 62 |
| unmatched | 32 |

All 62 fuzzy matches and the 1 normalized-tier duplicate were reviewed by hand
(`results/strain_reconciliation/NEEDS_REVIEW.md`): every fuzzy match was coincidental
string similarity between distinct, sequentially-numbered strain codes (e.g.
`DBVPG_3768` vs `DBVPG_3776`), not the same physical strain — all rejected. The
normalized-tier duplicate (`TFCN_17_332M_1`, a separator-only respelling of exact-match
`TFCN_17-332M-1`) was also rejected to avoid double-mapping one genotype to two
phenotype records. **213/308 phenotype strains accepted** — 12 more than the prior run's
201, 0 removed.

The reconciliation's collision guard (spec review fix #4) found 47 VCF sample IDs each
claimed by 2+ phenotype strains before rejection resolved them; see `NEEDS_REVIEW.md`
for the full list.

## 2. Ploidy validation

`scripts/check_ploidy.py` computed per-strain heterozygosity (via one
`bcftools stats -s -` pass over all 422 VCF samples, ~56s vs. ~3h for 213 separate
per-sample calls) and cross-referenced mosdepth mean depth for the 213 accepted strains.

**Finding: the VCF encodes every genotype haploid** (`GT` like `0`/`1`, not `0/1`) — every
strain's `nHets` is 0, confirming the prior run's aggregate "0 het" check (D-9/L-16) at
per-strain resolution. With no heterozygosity signal available, the flag collapses to the
panel-relative depth ratio alone (median panel depth 51.5x):

| flag | count | criterion |
|---|---|---|
| haploid_ok | 153 | depth ratio < 1.5x, no het signal |
| watch | 14 | depth ratio in [1.5x, 2.5x] |
| watch_contamination | 46 | depth ratio > 2.5x (up to ~13x) |
| likely_diploid_or_mixed | 0 | requires het outlier AND 1.5-2.5x band — never triggers here |

**User decision (2026-08-25): exclude none for now.** All 60 flagged strains are kept;
the depth flags are informational pending downstream GRM diagnostics, not an exclusion.

## 3. Strain-state diff vs. prior run

`scripts/diff_strain_state.py` compared the reconciled+ploidy-cleared strain set against
the prior run's recorded state (`data/prior_run_state/`): strain list, population
assignment, and near-clone/culled-173 partition (spec §4, review fix #1/#5).

- Strain list: **12 added** (`DBVPG_3239, 3445, 3446, 3538, 4379, 4952, 6094, 6649,
  TFCN_223A-8, TFCN_25-332D-2, TFCN_25-332M-2, TFCN_270H-1`), 0 removed.
- Population assignment: **0 changes** among the 201 strains common to both runs.
- `near_clone_partition_touched`: **false** by the diff script's check — but this check
  only tests whether an added/removed strain intersects the *prior* strain/culled sets,
  which is a no-op for pure additions (an added strain can never have been in the prior
  set by definition). It does **not** verify that the culled-173 partition, if
  recomputed on the new 213-strain set, would come out the same — that recomputation was
  explicitly deferred (see §5). This is a known blind spot in the diff logic, not a
  verified "unaffected" result.
- **Decision: `overall_action = rebuild`** (strain list changed → mandatory GRM
  rebuild + conditioning diagnostic + Tier A/B rescan, per spec §4).

## 4. GRM conditioning diagnostic

`scripts/rebuild_tiers_abc.sh` rebuilt the kinship matrix on the 213-strain panel:
`bcftools view -S` subset → biallelic → `plink2` QC (`--geno 0.1 --mind 0.1 --mac 5`,
**`--set-all-var-ids '@:#:$r:$a'`** — the source VCF has `.` for every variant ID, which
PLINK2's `--indep-pairwise` rejects as non-unique; this fix was needed but undocumented
in the original ad hoc PROGRESS.md run) → LD-pruned to **29,453 variants** (up from the
prior run's 20,769, reflecting both the larger panel and the different MAF spectrum) →
`gemma -gk 1` kinship.

`scripts/check_grm_conditioning.py` computed the eigenvalue spectrum:

```
n=213  condition_number=2.146e+11  near_zero_eigvals=3 (1.4%)
Verdict: singular_risk (by the condition-number threshold)
```

**Investigated, not taken at face value.** Two interpretations were distinguished by
inspecting the near-zero eigenvectors' loadings and the pairwise kinship matrix:
- One eigenvalue (2.07e-11) aligns exactly with the all-ones vector — the expected,
  harmless null space of any mean-centered GRM.
- The other two (machine-zero, ~-5e-17) concentrate almost entirely on **DBVPG_5757,
  DBVPG_5758, DBVPG_5759** — a genuine near-clone triplet (same collection batch: Italy,
  soil, same year code; already grouped in the same population, Pop 3, in
  `Rmuc_PopAssigned.csv`).

This is the same phenomenon D-9 already characterized (near-clone strains degenerate the
GRM: there, 22/201 strains; here, a much smaller 3/213 triplet) — not a new pathology.
Per user decision, this was treated as consistent with D-9's precedent (kinship-only LMM
tolerates it) and the rebuild proceeded.

## 5. Population structure

Before extending to Tier B/C/LOCO, the user asked whether population structure in this
panel needs additional correction beyond the existing kinship-only LMM. Checked directly
by eigen-decomposing the rebuilt 213-strain GRM and cross-referencing
`Rmuc_PopAssigned.csv`'s 6 population labels (sizes: pop1=77, pop2=28, pop3=30, pop4=39,
pop5=18, pop6=21):

- **PC1 explains 33.3% of GRM variance, PC2 9.9%** — both cleanly separate the 6
  populations (per-population PC1/PC2 means are well-resolved, most with tiny within-pop
  SD; pop2 and pop3 show more within-pop spread, consistent with pop3 containing the
  DBVPG_5757/5758/5759 near-clone triplet from §4).
- This confirms real, strong population structure — consistent with pixy's mean Fst
  ~0.45 from the original run (PROGRESS.md §8) — but **kinship-only LMM (D-9) is already
  the correct tool for it**, not a gap needing a new fix: D-9 established that adding
  explicit PC or population-dummy covariates on top of the kinship term collapses the
  model (GSL solver singular) precisely because top genotype PCs are collinear with the
  near-clone structure the GRM already encodes. The GRM captures population membership
  at finer, continuous resolution than 6 discrete labels ever could.
- **Data-cleaning conclusion**: no additional population-structure correction is
  warranted beyond what's already in place. The two concrete cleaning steps that *are*
  warranted and now in place are (a) the near-clone culling in §6 (a large fraction of
  "structure" in a panel like this is redundant near-identical strains, not distinct
  ancestry) and (b) the per-hit population-confounding check in §10 below (catches the
  case where a single locus, not genome-wide relatedness, happens to be a population
  marker — kinship correction does not protect against that).

## 6. Near-clone culling (reconstructed)

The 173-strain "informative subset" culling algorithm from the original run was never
saved as code — only described in prose (PROGRESS.md §6 N2: "IBS0<0.005 greedy").
Reconstructed as `scripts/cull_near_clones.py`: pairwise IBS0 (allele-mismatch rate; this
panel is haploid-encoded genome-wide, so IBS0 reduces to a direct mismatch rate) via
`plink2 --make-king-table` (this cluster's plink2 build has no `--distance` flag — L-25),
then greedy removal of one member of the closest sub-threshold pair, repeated until no
pair remains below 0.005.

**Validated before trusting it on new data**: restricted to the prior run's known
201-strain set, the reconstruction kept 171/201 strains (30 removed) vs. the original's
173/201 (28 removed) — **162/173 (93.6%) membership overlap**. Not exact (the original's
tie-break rule for *which* clone in a pair to drop was never recorded), but close enough
in scale and membership to trust on the new panel.

Applied to the full 213-strain panel: **182/213 kept, 31 removed (14.6% culled)** —
proportionally consistent with the prior run's 13.9%. This "gwasc" (culled) panel got its
own rebuilt kinship + GRM diagnostic: condition_number=2.0e11 (flagged `singular_risk` by
the same threshold as §4), but only **1/182 (0.5%) near-zero eigenvalues** — smaller than
the full panel's 3/213, consistent with the same benign pattern (not re-investigated
eigenvector-by-eigenvector since the fraction is even lower here).

## 7. Tier A — corrected (full unpruned SNP set, both panels)

**A real bug, caught before Tier B**: the first Tier A rebuild (this doc's earlier
version) used the LD-pruned SNP set for the association scan itself, not just for
kinship. The original run used the FULL QC'd unpruned set (~404,706 SNPs) for
association, pruning only for kinship (PROGRESS.md's S3/S3b rows). An independent
quant-genetics consult (D-16) confirmed unpruned-for-association is correct — pruning the
scan silently drops >90% of genotyped positions and can exclude the true signal or its
best proxy; this panel's population-structure LD (Fst~0.45, PC1=33%) is long-range and
ancestry-driven, so pruning-by-window would gut resolution for the wrong reason. The
consult also recommended reporting a Meff proxy + Bonferroni-at-Meff alongside BH-FDR,
and a per-hit population-confounding check (§10) — both added.

`scripts/rebuild_full_genotypes_and_tiera.sh` regenerated + persisted the full unpruned
QC'd bfile for both panels (gwas=213: 498,484 variants; gwasc=182-culled: 498,484
variants pre-filter, same source) and reran Tier A (12 traits × 2 panels = 24 scans,
~13 min/scan). `scripts/summarize_tiera.py` (now with `--n-pruned-snps` for the Meff
proxy):

**gwas (213-strain) panel** (Meff proxy = 29,453 pruned SNPs, Bonferroni-at-Meff =
1.70e-6):

| trait | n SNPs | λ | n FDR05 | n sig (Bonf@Meff) | top SNP | top p_wald |
|---|---|---|---|---|---|---|
| AUC_0 | 496,358 | 1.844 | 9 | 12 | scaffold_2:1248484 | 1.0e-7 |
| AUC_10 | 496,358 | 0.619 | 2,973 | 54 | scaffold_6:228721 | 2.1e-8 |
| AUC_20 | 498,459 | 1.129 | 27 | 19 | scaffold_16:417619 | 1.9e-9 |
| AUC_30 | 498,459 | 0.572 | 28 | 25 | scaffold_6:113610 | 1.2e-8 |
| AUC_ratio_10 | 496,358 | 0.841 | 0 | 0 | scaffold_4:61737 | 2.0e-5 |
| IC50_est | 496,364 | 0.414 | 1 | 1 | scaffold_16:122361 | 3.7e-11 |
| bright | 496,358 | 0.544 | 0 | 0 | scaffold_7:784232 | 1.3e-5 |
| chroma | 496,358 | 0.473 | 0 | 9 | scaffold_8:38068 | 2.5e-7 |
| clone_mean_area | 496,358 | 1.642 | 0 | 0 | scaffold_3:685221 | 1.4e-5 |
| cu_dose_slope | 496,347 | 0.374 | 24 | 20 | scaffold_3:546065 | 1.7e-8 |
| resilience_30 | 496,358 | 0.547 | 665 | 372 | **scaffold_13:810026** | 2.6e-9 |
| sat | 496,358 | 0.411 | 41 | 30 | scaffold_1:208569 | 1.4e-7 |

**gwasc (182-strain culled) panel** (Meff proxy = 28,707, Bonferroni-at-Meff = 1.74e-6):

| trait | n SNPs | λ | n FDR05 | n sig (Bonf@Meff) | top SNP | top p_wald |
|---|---|---|---|---|---|---|
| AUC_0 | 496,068 | 1.764 | 11 | 11 | scaffold_14:153892 | 9.6e-8 |
| AUC_10 | 496,068 | 0.641 | 3,181 | 54 | scaffold_7:68713 | 2.8e-9 |
| AUC_20 | 496,605 | 1.012 | 20 | 20 | scaffold_16:417619 | 1.5e-7 |
| AUC_30 | 496,605 | 0.562 | 158 | 32 | scaffold_6:113610 | 4.5e-9 |
| AUC_ratio_10 | 496,068 | 1.273 | 0 | 0 | scaffold_1:1436119 | 2.9e-5 |
| IC50_est | 493,088 | 0.403 | 847 | 616 | scaffold_10:458650 | 1.0e-9 |
| bright | 496,068 | 0.683 | 0 | 0 | scaffold_2:61063 | 2.6e-5 |
| chroma | 496,068 | 0.412 | 0 | 0 | scaffold_2:1406833 | 2.7e-6 |
| clone_mean_area | 496,068 | 1.654 | 0 | 0 | scaffold_4:123170 | 3.5e-6 |
| cu_dose_slope | 496,068 | 0.356 | 1 | 8 | scaffold_7:897138 | 5.9e-10 |
| resilience_30 | 496,068 | 0.548 | 2,533 | 17 | **scaffold_13:810026** | 8.7e-10 |
| sat | 496,068 | 0.484 | 0 | 0 | scaffold_3:368161 | 3.3e-6 |

**resilience_30/AUC_30's scaffold_13:810026 anchor is the standout robust finding**,
replicating across *every* panel/method variant tested this session: prior 201-strain run
(p=6.4e-9), this session's first (buggy, pruned-only) rebuild (p=2.6e-9), the corrected
full-213 rebuild (p=2.6e-9), and the full-182-culled rebuild (p=8.7e-10). See §8 (LOCO)
and §10 (population confounding) for two more independent lines of evidence this is real.

**chroma is unstable across panel/SNP-set choices** — its top hit moves every time: prior
run `scaffold_10:384905` (λ=0.32, 345 FDR-sig) → this session's pruned-only rebuild
`scaffold_8:831789` (λ=1.01, 1 FDR-sig) → corrected full-213 `scaffold_8:38068` (λ=0.47,
0 FDR-sig but 9 Bonf-sig) → full-182-culled `scaffold_2:1406833` (λ=0.41, 0 FDR-sig).
**Do not cite any chroma locus as a confirmed/replicated finding without further
investigation** (todo: `analysis/gwas/results/gwas/tierA_summary/population_confounding_gwas.csv`
does not flag chroma's hit as population-confounded, so the instability isn't simply that
— possibly a genuinely weak/diffuse signal that different marker sets and strain
compositions pick up different marginal SNPs for).

## 8. LOCO (leave-one-chromosome-out sensitivity)

`scripts/run_loco.sh` (reusing the original's `run_loco_shared.sh`/`plink_arch.sh`,
verified to have been actually saved as code, unlike the culling algorithm) ran LOCO for
the 3 traits the original run actually covered (chroma, AUC_10, resilience_30 — PROGRESS.md's
"6 traits" prose turned out to mean 6 (trait,panel) combinations, not 6 distinct traits,
confirmed by inspecting the real `loco/output/` directory contents) × both panels ×
~21 chromosomes.

**resilience_30's scaffold_13:810026 anchor reproduces almost exactly when chr13 itself
is excluded from the kinship**:

| panel | Tier A p (chr13 in kinship) | LOCO p (chr13 excluded) |
|---|---|---|
| gwas (213) | 2.62e-9 | 2.51e-9 |
| gwasc (182-culled) | 8.68e-10 | 8.34e-10 |

This rules out the signal being a kinship-absorption artifact (a spurious association
that only appears because the causal chromosome's own relatedness structure is baked
into the correction) — the association survives essentially unchanged even when its own
chromosome cannot contribute to the GRM.

Genome-wide LOCO λ (median across chr, gwas panel): AUC_10 0.608, chroma 0.458,
resilience_30 0.538 — tracks the Tier-A-per-chromosome λ (median 0.520) closely, meaning
the λ<1 deflation is stable near-clonal-structure correction, not a LOCO-specific
rescue effect (same conclusion as the original run's D-13 finding, now confirmed on the
213-strain panel too).

## 9. Tier C — BSLMM (5 traits, gwas-213 panel only, matching original scope)

`scripts/run_tierc_bslmm.sh`, `-bslmm 1` on chroma/AUC_10/AUC_30/clone_mean_area/
resilience_30 (the original run's 5 representative traits — confirmed by inspecting
`tierC_summary/` contents, no culled-panel BSLMM was ever run originally either).

**MCMC chain length caught and fixed mid-session (L-28)**: a first attempt guessed
`-w 20000 -s 100000` from PROGRESS.md's "100k MCMC, 20% burn-in" prose read as ~120k
total iterations — this produced only 10,000 recorded posterior samples (`.hyp.txt` row
count), 10x shorter than the original's 100,000. GEMMA's default `-rpace` (record pace)
is 10, so "100k MCMC" in the prose meant 100k *retained* samples, requiring
`-s 1000000` (not `-s 100000`). Verified directly by diffing `.hyp.txt.gz` row counts
against the original (100,001 both) before trusting any posterior estimate. The
short-chain outputs are kept (not deleted) under `tierC_summary/superseded_short_chain/`
as a convergence sensitivity check, not presented as results.

`scripts/summarize_tierc.py` computes PVE (phenotypic variance explained by all SNPs) and
PGE (proportion of that from the sparse/large-effect component) posteriors from
`.hyp.txt.gz`, and per-SNP posterior inclusion probability (PIP) from `.gamma.txt.gz` —
truncating each MCMC sample's row to its own `n_gamma` before tallying (GEMMA
zero-pads gamma rows, indistinguishable from a genuine SNP index 0 by value alone; caught
via a separate bug where `.hyp.txt.gz`'s trailing empty 7th tab-field silently shifted
`n_gamma` to `NaN` under plain `pd.read_csv` — fixed with `index_col=False`):

| trait | PVE [95% CI] | PGE | n_gamma (med) | top PIP locus |
|---|---|---|---|---|
| chroma | 0.213 [0.093, 0.395] | 0.472 | 14.0 | scaffold_9:667240 (PIP=0.018) |
| AUC_10 | 0.398 [0.254, 0.577] | 0.784 | 3.0 | scaffold_5:458816 (PIP=0.440) |
| AUC_30 | 0.407 [0.224, 0.652] | 0.537 | 5.0 | scaffold_16:492282 (PIP=0.218) |
| clone_mean_area | 0.087 [0.013, 0.250] | 0.338 | 7.0 | scaffold_10:274947 (PIP=0.004) |
| resilience_30 | 0.245 [0.094, 0.462] | 0.506 | 72.0 | scaffold_13:793374 (PIP=0.122) |

**Substantively different architecture from the original run.** The original found
near-oligogenic PIP≈1.0 loci for several traits (e.g. AUC_10: 2 loci at PIP=1.0,
"near-oligogenic"; chroma: 4 loci at PIP=1.0 on scaffold_3). This rebuild finds **no locus
reaching PIP≥0.5 for any trait** — the highest is AUC_10's scaffold_5:458816 at PIP=0.44.
Not yet explained; candidate reasons include the larger unpruned marker set here (498,484
vs. the original's 404,706 SNPs, diluting PIP across more candidate tag SNPs) and/or the
different, smaller near-clone structure (§4/§6) changing which variants the sparse
component favors. Flagged for investigation, not presented as a contradiction of the
original's architecture claims.

One reassuring cross-check: **resilience_30's top BSLMM locus (`scaffold_13:793374`,
PIP=0.122) sits only ~17 kb from the Tier A/LOCO-confirmed anchor `scaffold_13:810026`**
(§7/§8) — independent (if lower-confidence) corroboration that scaffold_13 harbors real
signal for this trait, from a completely different modeling approach (sparse Bayesian
selection vs. single-SNP LMM).

Per-SNP PIP tables: `results/gwas/tierC_summary/pip/{trait}_pip.csv`.

## 10. Population-confounding check (per top hit)

`scripts/check_population_confounding.py`: for each trait's top FDR-significant SNP,
cross-tabs alt-allele frequency across the 6 populations; flags a >0.8 AF swing (while
overall AF isn't itself near-fixed) as a population-private-variant pattern kinship
correction alone would not catch.

- **resilience_30's scaffold_13:810026 is clean in BOTH panels** (AF range 0.000–0.042,
  overall AF 0.014) — a third independent line of evidence (with cross-panel replication
  §7 and LOCO §8) that this signal is real, not a population artifact.
- **AUC_20's scaffold_16:417619 is flagged in BOTH panels** (AF range 0.02–0.995) — this
  locus is essentially a population marker; treat any AUC_20 finding at this locus with
  caution.
- `cu_dose_slope` and `sat` also flagged in the gwas-213 panel only (not in gwasc),
  suggesting the culled panel's removal of near-clones changes which population
  drives the apparent signal for those traits.

Full tables: `results/gwas/tierA_summary/population_confounding_{gwas,gwasc}.csv`.

## 11. Tier B (SKAT/burden set tests) — complete

Ran via `scripts/run_tierb.sh` (ported `tierb_set_tests.py`, unmodified logic) against
both panels' full unpruned bfiles for window LD, using pixy's high-dxy window universe
(211-215 windows per panel; ~196/211 with a usable LD basis after skipping windows whose
LD matrix failed to converge). **Caveat, not resolved**: `--pixy-dir` points at the
original run's pixy output, computed on a `cohort.all.vcf.gz` restricted to the prior
**201-strain** cohort (confirmed via `bcftools query -l`, returns 201 samples) — pixy was
**not** recomputed for the 213-strain panel this session (the original pixy run took
6h35m). The high-dxy window *definitions* are a population-genetic property of the genome
and unlikely to shift dramatically from 12 more strains within existing populations, but
this remains an approximation.

**Result matches the original run's pattern exactly: no set-level (burden or SKAT)
signal survives BH-FDR(q<0.05) in either panel** (0/2352 gwas, 0/2340 gwasc for both
burden and SKAT, all-windows or high-dxy-only universe). `min_p` "significance" merely
recovers Tier A's single-SNP hits at window resolution and is not independent evidence
of set-level signal (2340-2340/2352 rows pass FDR on `min_p` — expected, not a discovery).
The moment-approximation SKAT p-values track exact Monte Carlo verification closely on
the top-25 windows (e.g. resilience_30/scaffold_16 window: 0.000676 approx vs. 0.00282
MC — same order of magnitude, consistent direction).

**resilience_30's scaffold_13:810026 anchor is NOT in a high-dxy window** (its
800,001-900,000 window has `is_highdxy=False`, SKAT p=0.32/0.46 gwas/gwasc, burden
p=0.54/0.67 — unremarkable at the set level) — consistent with the original run's finding
that FDR-significant GWAS loci are not enriched in high-divergence regions (standing
variation within the near-clonal focal clade, decoupled from deep population splits).

Full tables: `results/gwas/tierB/tierb_settests_{gwas,gwasc}.csv`,
`tierb_skat_mcver_{gwas,gwasc}.csv`.

## 12. Raw CIELAB traits (L*, a*, b*) — added 2026-08-26 on user request

Neither the original run nor this port's first pass tested raw CIELAB components
directly — only derived `chroma` (Lab-magnitude of a*/b*), `sat`, `bright` (HSV, not
Lab at all). Added `lab_L`, `lab_a`, `lab_b` (median, clone-mean over plates — same
construction as the existing color block) to `scripts/build_gwas_phenotypes.py`, and
ran Tier A (`scripts/run_tiera_lab_traits.sh`, full unpruned SNP set, both panels).

| trait | panel | n FDR05 | top SNP | top p_wald |
|---|---|---|---|---|
| lab_L | gwas (213) | 1,619 | scaffold_3:229210 | 5.8e-9 |
| lab_a | gwas (213) | 758 | scaffold_8:38068 | 9.9e-13 |
| lab_b | gwas (213) | 1,614 | scaffold_11:608491 | 9.2e-13 |
| lab_L | gwasc (182) | 942 | scaffold_3:265939 | 4.1e-10 |
| lab_a | gwasc (182) | 117 | scaffold_5:882396 | 3.5e-11 |
| lab_b | gwasc (182) | 687 | scaffold_3:503756 | 1.6e-10 |

**Far more FDR-significant SNPs than any derived color trait** (chroma/sat/bright top
out at 0-41). `lab_a`'s top hit in the gwas-213 panel, `scaffold_8:38068`, is the **same
SNP as chroma's top hit in that panel** (p=2.5e-7) — consistent with chroma being
mathematically derived from a*/b*, with the raw component carrying a much stronger
signal at that locus (p=9.9e-13).

**But: do not read this as strong single-locus color biology without more work.**
`check_population_confounding.py` flags **all three raw Lab traits' top hits as
population-confound risk in BOTH panels** (AF swings 0.03-0.98+ across the 6
populations — the same pattern already seen for `AUC_20`) — including `lab_a`'s
striking `scaffold_8:38068` hit. The top SNP is also unstable across panels for all
three traits (different locus in gwas-213 vs. gwasc-182 for every one), mirroring
chroma's instability (§7). The abundance of FDR-significant hits likely reflects
population structure showing through at the top-hit level despite kinship correction,
not necessarily a clean color-specific causal locus — this needs the same scrutiny
chroma already got flagged for, not a headline claim.

Full tables: `results/gwas/tierA_summary/{tiera_summary_gwas,tiera_summary_gwasc}.csv`
(now include lab_L/a/b), `population_confounding_{gwas,gwasc}.csv`.

## 13. Population-vs-locus disambiguation — most flagged hits are likely real

The crude AF-swing population-confounding screen (§10, §12) cannot distinguish a locus
that is causal *and* happens to be population-differentiated (legitimate — e.g. local
adaptation) from one that is merely correlated with population membership for unrelated
reasons (artifact). Designed a disambiguation battery
(`docs/superpowers/specs/2026-08-26-population-confounding-disambiguation-design.md`),
independently reviewed by a statistical geneticist/breeder consult before implementation
(that review caught a real flaw in the draft — the within-population test is exactly as
vulnerable to near-clone pseudoreplication as the genome-wide scan kinship correction
exists for — and the fix, per-population clone-collapse on the 182-strain culled panel,
is built into `scripts/check_population_vs_locus.py`).

**Four tests**: (A) within-population re-test (phenotype ~ genotype, fit separately per
population, clone-collapse-aware), (B) fixed-effect meta-analysis of A's per-population
betas, (C) single-locus population-covariate regression (partial R² of genotype beyond
population), (D) Fst×MAF-matched empirical null. **Test D turned out to be
uninformative in practice** — caught empirically, not anticipated in the design: it
returned an "extreme" percentile for literally every flagged locus, including one with
no nominal within-population significance at all, because these loci are each trait's
single most significant genome-wide hit by construction (a winner's-curse artifact, not
evidence). The verdict rule was tightened to require nominal within-population
significance (raw meta_p<0.05) as a hard gate rather than letting test D's circular
result substitute for it.

**Result: 8 of 9 flagged (trait, panel) combinations verdict `likely_real`**:

| trait | panel | top SNP | meta_p (within-pop) | n pops (independent) | partial R² beyond pop | verdict |
|---|---|---|---|---|---|---|
| cu_dose_slope | gwas | scaffold_3:546065 | 4.8e-13 | 2 | 0.153 | likely_real |
| lab_L | gwas | scaffold_3:229210 | 4.3e-5 | 2 | 0.034 | likely_real |
| lab_a | gwas | scaffold_8:38068 | 8.8e-12 | 2 | 0.139 | likely_real |
| lab_b | gwas | scaffold_11:608491 | 2.9e-4 | 3 | 0.113 | likely_real |
| sat | gwas | scaffold_1:208569 | 1.5e-8 | 3 | 0.095 | likely_real |
| lab_L | gwasc | scaffold_3:265939 | 2.2e-3 | 2 | 0.034 | likely_real |
| lab_a | gwasc | scaffold_5:882396 | 1.1e-10 | 2 | 0.137 | likely_real |
| lab_b | gwasc | scaffold_3:503756 | 3.8e-4 | 2 | 0.116 | likely_real |
| **AUC_20** | gwas | scaffold_16:417619 | **0.219** | 2 | 0.021 | **ambiguous_underpowered** |

Concretely, for `lab_a`: only populations 3 and 6 had enough within-population allelic
variance to test, and **both independently show the same direction of effect**
(pop3 β=+3.23, p=0.002; pop6 β=+1.65, p=1.1e-5) — this is genuine independent
replication after removing all between-population variance by construction, not an
artifact of pooling correlated samples. Same pattern for `cu_dose_slope` (pop3 β=+0.034
p=0.0016; pop6 β=+0.023 p=6.4e-6).

**Revises the earlier population-confounding conclusion (§10, §12)**: the AF-swing
heuristic was catching loci that are population-differentiated *and* causally real, not
pure artifacts — chroma remains a separate, still-open question (§7) since it wasn't
re-tested here (its top hit changes across panels, unlike the stable loci above, so it
didn't cleanly enter this battery as a single fixed locus). `AUC_20` remains the one
locus that should still be treated with caution — its within-population meta-analysis
isn't even nominally significant (p=0.22).

**Not yet done — explicit MAS-usability gates** (deferred per the design, not part of
this disambiguation pass): minimum effect-size threshold, held-out replication, LD-decay
check confirming the flagged SNP isn't a distant tag of a stronger nearby signal, and —
specific to the color traits — a plate/batch-confound sanity check (color measurement is
exactly the kind of trait vulnerable to batch effects that could correlate with which
population was measured on which plate/day).

Full results: `results/gwas/population_vs_locus/population_vs_locus.csv`,
per-population detail in `within_pop_{trait}_{panel}.csv`.

## 14. MAS-usability gates for the `likely_real` loci

Deferred at §13 (D-17), run here per user request: minimum effect size, LD-decay/
best-candidate-in-block check, plate/batch-confound check, and (color traits only)
timepoint robustness across early/mid/late growth windows. `scripts/check_mas_gates.py`.

**Effect size**: all 8 loci have large effect sizes by conventional standards (|Cohen's
d| 0.79–1.58, i.e. the genotype groups differ by 0.74–1.28 phenotype SDs) — none are
p-value-only, practically negligible hits.

**LD-decay**: this cluster's plink2 build has no windowed `--r2`/`--ld-snp` report either
(L-25's gotcha recurring) — used `--clump` instead, which directly answers the intended
question. **All 8 flagged SNPs are the clump LEAD in their own ±100kb LD block** — no
more-significant candidate was found nearby, so the flagged marker is at least as good as
any alternative in the region. Caveat: each block has 425–1,640 SNPs at r²≥0.2 with the
lead, reflecting this panel's substantial background LD (consistent with the population
structure already characterized) — the causal *region* is well-supported, but this check
alone does not fine-map to a single gene.

**Batch/plate confounding**: every strain in this panel was measured in exactly ONE of 4
runs (353–356) — `run_number` is a fixed per-strain label, not a repeated covariate, so
it's a real confounding risk if genotype clusters by run. No locus shows a
genotype-vs-run association (all p>0.05), and for 7/8 loci adding run as a covariate
barely moves the genotype term's partial significance (e.g. cu_dose_slope 8.1e-8→2.2e-7).
**Exception: `lab_L` (gwasc panel)** — partial p crosses the 0.05 threshold once batch is
controlled (0.015→0.055) — flagging it as the most batch-sensitive/fragile of the eight;
its `gwas`-panel counterpart is more robust (0.014→0.031, stays nominally significant).

**Timepoint robustness** (color traits only — `cu_dose_slope`/`AUC_20` are dose-response
derived, not single-timepoint color reads, out of scope here): recomputed each trait at
early (18–42h), mid (48–72h), and late (85–110h, the original) growth windows and re-ran
the within-population test at each. **All four color traits (lab_L, lab_a, lab_b, sat)
show directionally consistent significance across all three windows, in both panels** —
most *strengthen* monotonically with developmental time, e.g. `lab_a`: p=7.9e-4 (early) →
2.0e-8 (mid) → 8.8e-12 (late), a biologically sensible pattern (pigment/color signal
accumulating as the colony matures), not an artifact of the single late-window choice the
original phenotype construction used.

**Overall**: 7 of 8 loci (`cu_dose_slope`, `lab_a`, `lab_b`, `sat` in both panels tested)
pass every gate cleanly. `lab_L` remains directionally consistent and significant in 5/6
checks but is the weakest and most batch-sensitive of the set — treat its locus with more
caution than the others pending further work.

**Still not done**: held-out replication (no genuinely disjoint strain subset currently
exists to test on) and a formal fine-mapping pass to narrow each broad LD block toward
candidate genes (natural continuation of the deferred Tier D/E/G work, §"Next steps").

Full results: `results/gwas/mas_gates/mas_gates.csv`, per-window detail in
`timepoint_{trait}_{panel}.csv`.

## Session summary (2026-08-25/26)

All of Tier A (corrected), Tier B, Tier C, and LOCO are now complete for this port. The
single most robust finding across every method and panel variant tested — original
201-strain run, this session's (buggy) pruned-only rebuild, the corrected full-213
rebuild, the full-182-culled rebuild, LOCO (chr13 excluded from kinship), and the
population-confounding check — is **`resilience_30`/`AUC_30`'s anchor at
`scaffold_13:810026`** (p ranging 6.4e-9 to 8.7e-10 depending on panel), independently
corroborated at lower confidence by Tier C BSLMM's nearby top locus
(`scaffold_13:793374`, ~17 kb away). `chroma`'s signal is unstable across every SNP-set
and panel choice tried and should not be cited without further investigation.

## Next steps (not done this session)

- Investigate the chroma signal instability (§7) — it is NOT population-confounded per
  §10, so the cause is still open.
- Investigate lab_L/a/b's population-confounding flags (§12) before treating their large
  FDR-significant hit counts as real color-specific signal — same open question as
  chroma, now affecting three more traits.
- Consider whether pixy should be recomputed on the 213-strain panel (currently a
  documented approximation, §11) — Tier B's headline conclusion (no set-level signal;
  resilience_30 anchor not in a high-dxy window) matches the original run closely enough
  that this is lower priority than initially flagged, but still an approximation.
- Tier B/C/LOCO were not rerun for lab_L/a/b (§12 is Tier A only) — natural follow-up if
  those traits turn out, after the confounding investigation above, to carry real signal.
- Tier D/E/G (gene annotation, fine-mapping, prior-locus replication) — not started this
  session; natural next step, focused on the resilience_30/AUC_30 scaffold_13:810026
  anchor given how robustly it's replicated across every check run so far.
