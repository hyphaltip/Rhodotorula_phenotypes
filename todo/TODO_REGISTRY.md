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
| Dataviz consult for Tier D/E figure before publication | medium | open | gwas | 2026-08-17 | jstajich | `results/gwas/figures/tierde_gene_finemap.{png,pdf}` ready for review |
| Resolve rare-EF loci via denser genotyping/imputation — auc10 DBP3 (scaffold_10:396172) CS n=67 and chr13 block (11,701–800,664, OM429_005439) are wide sets; LD (n=201) limits resolution | high | open | gwas | 2026-08-17 | jstajich | see `todo/rare-ef-imputation.md` |
| Functional follow-up of GWAS candidate genes — OM429_005439 (chr13 replicated locus, unannotated), OM429_004640/DBP3 (AUC_10), OM429_004009 telomerase RT (chroma) | medium | open | gwas | 2026-08-17 | jstajich | see `todo/functional-followup.md` |
| Reconstruct near-clone IBS0-culling algorithm as reusable code and recompute culled-173-equivalent subset for the 213-strain gwas port panel | high | done | gwas | 2026-08-26 | jstajich | `cull_near_clones.py`, validated 93.6% overlap vs. prior 173-set; 182/213 kept on new panel |
| Rebuild Tier B (SKAT/burden) + Tier C (BSLMM) + LOCO on the rebuilt 213-strain kinship (`analysis/gwas/`) once the culling algorithm above is reconstructed | high | done | gwas | 2026-08-26 | jstajich | Tier B: no set-level signal (matches original). Tier C: more polygenic architecture than original (no PIP>=0.5 loci). LOCO: confirms scaffold_13:810026 anchor survives chr13 exclusion from kinship |
| Investigate chroma Tier-A signal instability across every panel/SNP-set variant tried (prior scaffold_10:384905 -> pruned-only scaffold_8:831789 -> full-213 scaffold_8:38068 -> full-culled-182 scaffold_2:1406833); NOT population-confounded per `check_population_confounding.py`, cause still open | medium | open | gwas | 2026-08-26 | jstajich | `analysis/gwas/GWAS.md` S7 |
| Investigate lab_L/lab_a/lab_b population-confounding flags before treating their large FDR-sig hit counts as real color signal | medium | done | gwas | 2026-08-26 | jstajich | check_population_vs_locus.py: lab_L/a/b (both panels), cu_dose_slope, sat all verdict likely_real (independent within-population replication in >=2 pops); only AUC_20 remains ambiguous. See `analysis/gwas/GWAS.md` S13 |
| MAS-usability gates for the likely_real loci from the population-vs-locus disambiguation | high | done | gwas | 2026-08-26 | jstajich | check_mas_gates.py: 7/8 loci pass all gates (large effect sizes, clump-lead LD status, no batch confound, timepoint-stable); lab_L flagged as weakest/most batch-sensitive. See `analysis/gwas/GWAS.md` S14 |
| Fine-map the 8 likely_real loci's broad LD blocks (425-1640 SNPs each at r^2>=0.2 within +/-100kb) toward candidate genes -- natural continuation of the deferred Tier D/E/G work, now with a validated starting locus set | high | open | gwas | 2026-08-26 | jstajich | `analysis/gwas/GWAS.md` S14; reuse `annotate_gwas_loci.py`/`finemap_credible_sets.py` from analysis/ideas/.../scripts/ as a starting point |
| Held-out replication check for the likely_real loci once/if a genuinely disjoint strain subset becomes available (none exists currently) | medium | open | gwas | 2026-08-26 | jstajich | `analysis/gwas/GWAS.md` S14 |
| Investigate why lab_L is more batch-sensitive than lab_a/lab_b/sat (partial p crosses 0.05 once run_number is controlled for, gwasc panel) before treating its locus with the same confidence as the other three color traits | medium | open | gwas | 2026-08-26 | jstajich | `analysis/gwas/results/gwas/mas_gates/mas_gates.csv`; `analysis/gwas/GWAS.md` S14 |
| Resolve chroma's identity separately from the population-vs-locus battery -- its top hit is unstable across every panel/SNP-set tried so it never entered the disambiguation battery as one fixed locus; needs its own investigation (e.g. run check_population_vs_locus.py-style within-pop test on each panel's own top hit, or find a stable candidate first) | medium | open | gwas | 2026-08-26 | jstajich | `analysis/gwas/GWAS.md` S7, S13 |
| Tier D/E/G (gene annotation, fine-mapping, prior-locus replication) for the analysis/gwas/ port, focused on the robustly-replicated resilience_30/AUC_30 scaffold_13:810026 anchor | high | open | gwas | 2026-08-26 | jstajich | `analysis/gwas/GWAS.md` "Next steps"; anchor replicated across 5 independent checks (2 panels x Tier A, LOCO, BSLMM proximity, population-confounding) |
| Consider recomputing pixy on the 213-strain gwas-port panel (Tier B currently reuses the prior 201-strain pixy output as a documented approximation) | low | open | gwas | 2026-08-26 | jstajich | Tier B's headline conclusion (no set-level signal; resilience_30 anchor not in high-dxy window) already matches the original closely, lowering priority |
<!-- Add new entries above this line -->
