# Finding: Color/growth GWAS loci do not co-localize with population-divergence (dxy/Fst) outliers

- **Date**: 2026-08-17
- **Source**: `analysis/ideas/2026-08-15-color-phenotype-space/` — LOCO sensitivity
  (`results/gwas/loco/`), Tier B set tests (`results/gwas/tierB/`), dxy/Fst
  co-localization (`results/gwas/tierB/coloc/`), GWAS_REPORT.md §8
- **Topic**: gwas, co-localization, dxy, fst, pixy, tierb, loco, sensitivity

## Result

1. **LOCO sensitivity: every Tier-A anchor reproduces at unchanged p.**
   Leave-one-chromosome-out GEMMA (kinship recomputed without the candidate chr,
   120 scans × {chroma, AUC_10, resilience_30} × {all-201, culled-173} × chr 1–20)
   gives chroma scaffold_10:384905 p=2.46e-8 (Tier-A 2.40e-8), AUC_10
   scaffold_10:396172 p=1.44e-8 (Tier-A 1.43e-8), resilience_30 scaffold_13:810026
   p=6.04e-9 (Tier-A 6.35e-9), plus the culled-set equivalents. The Tier-A signals
   are not kinship-absorption artifacts. LOCO λ (median across 20 scaffolds: 0.34–0.73
   by trait) tracks full-kinship Tier-A λ — the conservative (λ<1) inflation reflects
   near-clonal structure, and LOCO confirms it is stable rather than per-chr rescue.

2. **Tier B set tests (burden/SKAT over pixy high-dxy windows): no multi-SNP signal.**
   Neither burden (sum of z, λ≈0.5, over-conservative from +/− cancel with no
   direction prior) nor SKAT variance-component (λ≈0.6–0.8, top p≈1.4e-3, not
   replicating across all-201/culled) exceeds single-SNP Tier-A within any window.
   min-p across the window recovers the Tier-A hits (383/384 and 408/408 high-dxy
   windows FDR-significant). Causal content in high-dxy windows is concentrated in a
   few SNPs, consistent with Tier-C BSLMM near-oligogenic architecture (AUC_10
   PGE≈0.96, ~3 LD clusters).

3. **FDR-significant GWAS loci are NOT enriched in high-dxy or high-Fst windows.**
   Of 215 pixy 100-kb windows, 189 contain ≥1 FDR(q<0.05) SNP, but only 37 of these
   are top-20% high-dxy (37 expected; Fisher OR=0.81, p=0.61) and 34 top-20% high-Fst
   (37 expected; OR=0.41, p=0.065). Anchor loci sit at moderate divergence (Fst
   0.35–0.53, dxy 0.020–0.025, none high-dxy). Phenotype-contributing alleles segregate
   within the near-clonal focal clade (standing variation), not as fixed differences
   between the deep pop splits that create the dxy/Fst extremes.

4. **scaffold_20:100001 is a mis-assembly artifact, not a hotspot.** Genome-wide dxy
   max (0.070 vs 0.028 next) with ~0 Fst and only 12 genotyped SNPs → exclude from
   set tests and interpretation.

---

# Finding: Prior growth-rate locus chr13:13_30149 replicates in our AUC_10; gene mapping + fine-mapping of Tier-A anchors (Tier D/E/G)

- **Date**: 2026-08-17
- **Source**: `analysis/ideas/2026-08-15-color-phenotype-space/` — Tier D
  (`scripts/annotate_gwas_loci.py`, `results/gwas/tierD/`), Tier E
  (`scripts/finemap_credible_sets.py`, `results/gwas/tierE/`), Tier G replication;
  GWAS_REPORT.md §9
- **Topic**: gwas, tierd, tierte, tierg, annotation, finemapping, credible-sets, replication

## Result

1. **Independent-locus mapping** (Tier D): 12,348 FDR-sig SNPs → 5,286 independent loci
   (250 kb positional clump) across 9 traits. Notable genes: chroma scaffold_8 → telomerase
   RT (OM429_004009); AUC_10 scaffold_10 → DBP3 RNA-dependent ATPase (OM429_004640) +
   RNA-pol-I TF + endodeoxyribonuclease; BSLMM chroma scaffold_3 → methionine aminopeptidase 1
   (OM429_001379). Caveat: clump keeps only the single best SNP per chromosome → small
   scaffolds with a second distant block collapse (chr13's 217 FDR SNPs form ONE 525 kb block,
   single lead 13_791853 p=2.4e-7).

2. **Fine-mapping** (Tier E, Wakefield ABF z-space, prior SD=0.2, logsumexp, p<1e-3 candidate
   filter): chroma scaffold_10 95% CS n=24, lead pp=0.054, β=1.11±0.19, common AF (0.81) —
   the well-bounded common-variant anchor. Rare-EF loci (AF≈0.015) → wide sets (auc10 DBP3
   CS n=67, rare_driven) or fully-resolved singletons (resil scaffold_13 CS n=1, pp=0.615).

3. **Prior-locus replication** (Tier G): prior lab's chr13:13_30149 growth-rate hit
   (p=1.68e-11) replicates in our AUC_10 via proxy 13_30134 (15 bp away): p_wald=4.03e-6,
   FDR-sig, β=804,778, af=0.015 — same chr13 rare-haplotype block (lead 13_791853). Gene at
   locus = OM429_005439, a hypothetical protein with NO functional annotation (flanked 6.7 kb
   downstream by an Ark1-family Ser/Thr kinase). Other traits null there → replication is
   growth-phenotype-specific. Causal gene under a p≈1e-11 locus is functionally unknown —
   top validation target.

## Interpretation

Tier-A anchors split into common, well-mapped signals (chroma_10) and rare-EF signals that
LD (n=201) cannot resolve (auc10 DBP3, chr13 block). The replicated growth-rate locus maps
to an unannotated hypothetical gene — functional follow-up (OM429_005439, OM429_004640 DBP3,
OM429_004009 telomerase RT) is the bottleneck, not GWAS signal.


## Interpretation

The genetic basis of colour/growth/copper tolerance in R. mucilaginosa is a handful of
moderate-effect, intra-clade variants that are decoupled from population differentiation
outliers — a "standing-variation, not between-lineage divergence" architecture.
