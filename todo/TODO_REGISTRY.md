# TODO Registry

All future work items, ideas, and planned improvements for this project are tracked here. Each item has a dedicated `.md` file in this directory with full details.

## Status Key

| Status | Meaning |
|--------|---------|
| `open` | Not yet started |
| `in-progress` | Actively being worked on |
| `blocked` | Waiting on something external |
| `complete` | Done (kept for reference) |
| `wont-do` | Decided against (kept for rationale) |

## Priority Key

| Priority | Meaning |
|----------|---------|
| `critical` | Must be done — blocks progress or correctness |
| `high` | Important and should be done soon |
| `medium` | Valuable but not urgent |
| `low` | Nice to have |
| `idea` | Speculative — worth capturing but no commitment |

## Registry

| Item | Priority | Status | Category | Date | Author | File |
|------|----------|--------|----------|------|--------|------|
| Compare with public datasets | idea | open | validation | 2026-03-06 | Arjun Raj | [compare-public-data.md](compare-public-data.md) |
| Re-point skill-bridge personas source (Autonomous-Science repo gone) | medium | open | infrastructure | 2026-08-15 | jstajich | [fix-skill-bridge-personas.md](fix-skill-bridge-personas.md) |

| Re-analyze top GWAS trait(s) with LOCO kinship (per-scaffold GRM) to confirm Tier-A sensitivity | high | done | gwas | 2026-08-16 | jstajich | [loco-followup.md](loco-followup.md) |
| Tier C completion: FDR (q=0.05) significance done in `results/gwas/fdr/`; maxT permutation rejected (uncalibrated on near-clonal panel, L-21/D-11) | high | done | gwas | 2026-08-16 | jstajich | |
| Tier B: SKAT/burden set tests on pixy high-dxy windows (from Tier-A p-values) | medium | done | gwas | 2026-08-16 | jstajich | `tierb_settests_{gwas,gwasc}.csv` + MC-verified `tierb_skat_mcver_*.csv`; no set-level signal beyond single SNP (L-22, D-13) |
| Tier B gwasc MC-verify: sanity-check corr between moment-approx and exact MC for the culled set | low | done | gwas | 2026-08-17 | jstajich | r(log10)=0.997 (n=47), mean drift 0.171, worst cu_dose_slope p1.5e-4->1.1e-3 (0.871) still non-significant; verified in-session |
| dxy/Fst co-localization report section (GWAS_REPORT.md §8.3) figures finalized | low | done | gwas | 2026-08-17 | jstajich | `coloc_dxy_fst.{png,pdf}`, `coloc_enrichment.txt` |
| Dataviz consult for final figures (LOCO + Tier B + co-localization) before publication | medium | open | gwas | 2026-08-17 | jstajich | `results/gwas/figures/{loco_sensitivity,tierB_settests,coloc_dxy_fst}.*` ready for review |
<!-- Add new entries above this line -->
