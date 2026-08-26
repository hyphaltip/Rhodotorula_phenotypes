# GWAS: Color and Copper-Response Phenotypes in *Rhodotorula mucilaginosa*

**Status**: in progress (port from `analysis/ideas/2026-08-15-color-phenotype-space/`) — Tier A rebuild complete; Tier B/C/LOCO/pixy pending
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

## 5. Tier A GEMMA rebuild (complete) — Tier B/C/LOCO/pixy (pending)

**Scope decision (2026-08-25, user):** given the near-clone culling algorithm behind the
173-strain "informative" subset was never saved as code in the original run (only
described in prose, PROGRESS.md §6 N2: "IBS0<0.005 greedy"), reconstructing it plus Tier
B (SKAT/burden), Tier C (BSLMM), and LOCO was scoped **out** of this session. Only the
mandatory minimum (spec §4: GRM rebuild + Tier A/B single-SNP rescan) was completed, and
only the Tier A half of that (single-SNP LMM); Tier B set tests, Tier C BSLMM, and LOCO
are left as an explicit follow-up (see `todo/`).

`scripts/build_gwas_phenotypes.py` (ported from the ideas folder, re-pointed at the
rebuilt 213-strain `.fam`) derived the same 12 traits as the original run — chroma/sat/
bright/clone_mean_area (color+size block) and AUC_0/10/20/30, AUC_ratio_10,
resilience_30, cu_dose_slope, IC50_est (copper-response block) — with **213/213 strains
aligned** to the rebuilt `.fam` order.

`scripts/run_tiera_gemma.sh` ran `gemma -lmm 4 -k` per trait against the rebuilt
kinship (12 scans, ~8 min total). `scripts/summarize_tiera.py` computed genomic inflation
(λ) and BH-FDR(q=0.05) per trait:

| trait | n SNPs | λ | n FDR05 | top SNP | top p_wald |
|---|---|---|---|---|---|
| AUC_0 | 28,885 | 1.276 | 202 | scaffold_1:1221610 | 1.2e-6 |
| AUC_10 | 28,885 | 0.891 | 465 | scaffold_6:228721 | 2.1e-8 |
| AUC_20 | 29,441 | 0.865 | 5 | scaffold_16:418561 | 4.0e-7 |
| AUC_30 | 29,441 | 0.755 | 130 | scaffold_13:810026 | 1.6e-7 |
| AUC_ratio_10 | 28,885 | 1.056 | 0 | scaffold_7:1065753 | 3.0e-4 |
| IC50_est | 28,609 | 1.556 | 0 | scaffold_18:12058 | 1.2e-5 |
| bright | 28,885 | 1.071 | 0 | scaffold_2:1233513 | 1.5e-5 |
| chroma | 28,885 | 1.010 | 1 | scaffold_8:831789 | 1.3e-6 |
| clone_mean_area | 28,885 | 1.336 | 0 | scaffold_5:218638 | 1.6e-5 |
| cu_dose_slope | 28,884 | 0.810 | 243 | scaffold_3:546065 | 1.7e-8 |
| resilience_30 | 28,885 | 0.900 | 190 | **scaffold_13:810026** | 2.6e-9 |
| sat | 28,885 | 1.141 | 0 | scaffold_3:368161 | 3.2e-6 |

**Comparison to the prior 201-strain run (PROGRESS.md §9):**
- **Replicates cleanly:** `resilience_30`/`AUC_30` both anchor on **scaffold_13:810026**
  in both runs (this run p=2.6e-9 resilience_30 vs. prior 6.4e-9) — the 12 added strains
  did not disturb this signal.
- **Changed substantially:** `chroma`'s prior top hit was `scaffold_10:384905`
  (p=2.4e-8, 345 FDR-sig SNPs, λ=0.32); this run's top hit is `scaffold_8:831789`
  (p=1.3e-6, only 1 FDR-sig SNP, λ=1.01). The color signal is markedly weaker in the
  213-strain panel — flagged here rather than silently treated as a like-for-like
  replication; worth investigating whether the 12 added strains specifically dilute the
  chroma association (e.g. via population membership or phenotype range) before citing
  the original chroma finding as reconfirmed.
- λ values are generally closer to 1 here than the prior run's 0.357-0.638 range,
  consistent with less near-clone-driven over-correction on the (slightly) less clonal
  213-strain panel — expected given the added strains and the small (not large)
  near-clone triplet found in §4.

## Next steps (not done this session)

- Reconstruct the near-clone IBS0-culling algorithm as reusable code (never saved
  originally) and recompute the culled-173-equivalent set for 213 strains, to properly
  gate whether pixy/BSLMM/LOCO can be reused as-is or need a rebuild (spec §4).
- Tier B (SKAT/burden set tests) and Tier C (BSLMM) rescans on the rebuilt kinship.
- LOCO sensitivity check for the rebuilt panel's top hits (esp. the changed chroma
  signal).
- Investigate the chroma signal shift (§5) before citing it in any writeup.
