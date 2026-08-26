# Reconstruct near-clone IBS0-culling algorithm for the GWAS port

**Priority**: high
**Status**: open
**Category**: gwas
**Date**: 2026-08-25
**Author**: jstajich

## Context

The original run (`analysis/ideas/2026-08-15-color-phenotype-space/`) computed a
"culled-173" informative-strain subset from the 201-strain panel by greedy IBS0<0.005
near-clone removal (PROGRESS.md section 6, N2), used as the basis for the Tier B set
tests, Tier C BSLMM, and LOCO sensitivity scans (`gwasc_culled173.fam`). This algorithm
was run ad hoc on `$SCRATCH` and never saved as a reusable script — only described in
prose.

The `analysis/gwas/` port (D-15) added 12 strains to the panel (201 -> 213, 0 removed,
0 population changes) and rebuilt the kinship + Tier A single-SNP GEMMA scans. Tier B/C
and LOCO were explicitly NOT rebuilt this session because doing so correctly requires
recomputing the culled-173-equivalent subset for the new 213-strain panel first — the
`diff_strain_state.py` check that nominally covers this (`near_clone_partition_touched`)
is a no-op for pure strain additions and does not actually verify whether the culling
outcome would change.

## What's needed

1. Reconstruct the greedy IBS0<0.005 culling algorithm as a script under
   `analysis/gwas/scripts/` (e.g. `cull_near_clones.py`), computing pairwise IBS0 from
   the rebuilt 213-strain pruned genotype set
   (`analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas.pruned.{bed,bim,fam}`).
2. Verify it reproduces the prior run's 173/201 result when restricted to the original
   201 strains (sanity check before trusting it on the new 213).
3. Recompute the culled-N subset for the full 213-strain panel; diff against the prior
   `gwasc_culled173.fam` membership.
4. Feed the (possibly new) culled subset into Tier B (SKAT/burden), Tier C (BSLMM), and
   LOCO reruns on the rebuilt kinship.
5. Document the recomputed culling threshold sensitivity if it differs from the prior
   run's outcome.

## Related

- `analysis/gwas/GWAS.md` section 5 ("Next steps")
- `.living/decisions.md` D-15
- Prior Tier B/C/LOCO methodology: `.living/decisions.md` D-10, D-13
