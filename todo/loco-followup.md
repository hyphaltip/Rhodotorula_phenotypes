# LOCO Kinship Follow-up — Sensitivity Confirmation for Next-Gen GWAS

**Runbook**: `analysis/ideas/2026-08-15-color-phenotype-space/10-LOCO-RUNBOOK.md` (validated wrapper `run_loco.sh`); decision D-10.

**Status**: open (planned; D-10 decision — apply after Tier A/B/C)
**Priority**: high
**Date**: 2026-08-16
**Author**: jstajich
**Source decision**: `.living/decisions.md` D-10 (Option 2 proceed now, LOCO after)

## Why
Tier-A next-gen scans (12 traits × all-201 + culled-173) use the **full kinship** GRM
(`-lmm 4`). In an LMM the kinship matrix can over-absorb a SNP's own chromosomes' signal
("proximal contamination"), biasing true associations toward under-detection. LOCO
(leave-one-chromosome-out) rebuilds the GRM excluding the scaffold being tested, so the test
SNP's chromosome cannot subtract from its own signal.

The user explicitly wants LOCO applied to **confirm we had appropriate sensitivity** before
finalizing Stage-6 conclusions — not as a replacement but as a rigorous sensitivity check.

## What to do
1. Build **per-scaffold (per-chromosome) GRM** from the LD-pruned variant set for both the
   all-201 and culled-173 panels (GEMMA 0.98.3: use `-k` per scaffold, or partition pruned SNPs
   by chromosome and run `gemma -gk 1` per scaffold).
2. For the **most promising trait(s)** from Tier A (rank by top-signal + Tier B SKAT hits),
   rerun `gemma -lmm 4` **excluding** the scaffold under test from the kinship each time
   (i.e., for each scaffold `s`, test SNPs on `s` with a GRM built from all other scaffolds).
3. Compare lambda, top p-values, and hit rank stability vs the full-kinship Tier-A results.
   Flag any trait where a chromosome-proximal top hit disappears/gains substantially.

## Notes / gotchas
- GEMMA 0.98.3 has no single `-loco` flag; LOCO is achieved by per-scaffold GRM + repeated
  23 runs per trait, or a wrapper loop.
- Reuse the integer chromosome codes already in `gwas.bim` (scaffold_1..23 → 1..23).
- Only apply to top trait(s) first (cheap); expand if a meaningful delta is observed.

## Acceptance criteria
- Lambda and genome-wide hits for top trait(s) stable between full-kinship and LOCO, OR a clear
  dispersion where LOCO gains signal is documented in GWAS_REPORT.md Stage 6.
- Conclusion on whether full-kinship Tier A under-detected proximal loci.
