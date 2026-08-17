# Rare-Effect-Size Locus Resolution via Denser Genotyping / Imputation

## Status
open (priority: high), added 2026-08-17

## Why
Tier E fine-mapping (Wakefield ABF, `tierE_credible_sets.csv`) shows the interpretable
Tier-A copper signals are **rare-variant driven** and remain wide at n=201:

- **auc10_lead_10** (scaffold_10:396172, DBP3 / OM429_004640): 95% credible set n=67,
  lead af=0.015, `rare_driven=True` — the set does not narrow the causal SNP.
- **chr13 block** (scaffold_13:11,701–800,664; replicated growth-rate locus inside
  OM429_005439): 217 FDR-sig SNPs collapse to one < 250 kb block; the prior locus proxy
  (scaffold_13_30134, af=0.015) is indistinguishable from the rest of the block without
  denser markers.
- BSLMM chroma scaffold_3 clusters (anchor 3a/b/c): credible sets of n=98–131, lead pp≈0.01
  each — the BSLMM chroma signal is diffuse.

The ca. 0.015 allele-frequency floor and n=201 panel set an LD-limited resolution ceiling
for these loci (see GWAS_REPORT §9.2).

## Proposed work
1. Identify a variant-dense source: whole-genome calls for the same 201 strains at higher
   sensitivity (Variant QualityScore recalibration / lower MAC floor, or a second sequencing
   pass), or imputation against a broader Rhodotorula panel if available.
2. Re-run GEMMA `-lmm 4 -k` at the DBP3 and chr13 loci with the denser/ imputed marker set
   and repeat fine-mapping on the same anchors.
3. Report whether credible sets narrow below n≈10–20, and whether the chr13 causal variant
   can be localised under OM429_005439 vs the flanking Ark1-family kinase (OM429_005441).

## Depends on
- Source data (deeper genotyping or an imputation panel) — not yet identified.
- Coordinate names/positions per the chrom-mapping used in Tier A (`genome/Chrom_Mapping.tab`).

## Success criteria
- 95% ABF credible-set size decreases materially for auc10_lead_10 and/or the chr13 block.
- A single unambiguous candidate variant (or small functional set) resolved at chr13.
