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

## Candidate-gene coding variants at the misannotated carotenoid pathway genes (2026-08-26)

Sequence-level pilot of `analysis/candidate_gene_alignment/` on the 2 confirmed carotenoid
pathway genes (`OM429_003333` phytoene synthase/lycopene cyclase, `OM429_003336` phytoene
desaturase — see D-20's misannotation finding) across all 213 panel strains, snpEff-annotated:

- `OM429_003333`: 49 polymorphic CDS sites in the panel; 7 missense + 1
  `splice_donor_variant&intron_variant` (HIGH impact, snpEff) — a splice-donor change in a
  core carotenoid-synthesis gene is a plausible loss/change-of-function candidate worth
  flagging for whoever picks up functional follow-up, though it has not been tested against
  any phenotype here (this pipeline is descriptive, not a new association test).
- `OM429_003336`: 72 polymorphic CDS sites; 23 missense variants, 0 premature stops.
- Neither gene shows a nonsense (premature stop) variant in any of the 213 strains.
- Both genes' CDS+/-2kb regions carry 0 indels in the separate INDEL VCF (screened
  read-only) — the coding-variant picture above is complete at the SNP level for these 2
  genes specifically (not yet checked for other candidate genes).

This is a lead for follow-up (e.g. does the splice-donor genotype at `OM429_003333`
correlate with any color/pigment phenotype), not a validated finding — no statistical
test has been run on it yet.

## Candidate-gene coding variation DOES associate with color — at the GWAS-locus genes, not the carotenoid genes (2026-08-26, corrected)

Sequence pipeline extended to all 10 candidate genes (`analysis/candidate_gene_alignment/`,
see `results/candidate_gene_phenotype_assoc_all10.csv` + `..._all10.csv` and the indel
battery `candidate_indel_phenotype_assoc_all10.csv`), with the SAME population-aware
battery as the disambiguation exercise (within-pop re-test + meta + covariate + FDR).

1. **CORRECTION to the pilot (CRITICAL for any reader of the older finding above)**: the
   pilot's "0 indels" claim was a silent toolchain failure — `screen_indels.sh` reported
   the count of an empty bcftools pipe when bcftools wasn't on PATH (`2>/dev/null` +
   `|| true` masked it). The real INDEL VCF has 433/821 records in the two carotenoid
   genes' CDS±2kb and every one of the 10 target genes has ≥2 segregating indels INSIDE
   its CDS exons (`results/gene_coding_indel_screen.csv`). The SNP-only sequences/tables
   are the SNP layer only; indel genotypes must be tested separately (now done).
2. **The OM429_003333 splice_donor finding above is MONOMORPHIC in the 213-strain panel
   (0/213 alt at a PASS-quality site)** — it cannot be tested against phenotype and
   cannot drive within-panel color variation. The "splice-donor as top validation
   target" framing is answered: no testable variation there.
3. **NEW headline result — sat locus has a coding repeat-copy variant invisible to
   SNP-only analysis**: `OM429_000065` (sat-locus gene) carries c.839C>T p.Ala280Val
   (the lead SNP scaffold_1:208569) AND a +6 bp in-frame GAGCGG-repeat insertion at
   scaffold_1:208398 that is **100% concordant with the lead SNP** across all 213
   strains. Both associate with lab_a/chroma/sat under the population-aware battery
   (meta_p ≈ 8.8e-11, FDR q ≈ 1e-8, partial R² 0.086–0.10, replicated in 3/3
   testable populations). The sat locus's probable molecular driver may be a coding
   repeat-copy-number polymorphism — invisible to any SNP-only pipeline.
4. Other replicated candidate-gene coding hits (FDR-sig + `likely_real`): OM429_001521
   p.Gln245His (bright), OM429_005034 p.Tyr12His (lab_a/chroma/sat), OM429_002663
   p.Asn4Lys (lab_a/chroma/sat), OM429_001533 p.Tyr427Phe (sat/chroma/lab_b). These are
   the *coding layer* of loci already validated in D-17/D-18 — new contribution is the
   specific codon/AA.
5. **Neither carotenoid gene (OM429_003333/003336) shows an FDR-significant AND
   replicated color association** for either SNP or indel surface — consistent with the
   GWAS port's "carotenoid cluster not co-localized with any validated color locus."

Caveat: perfect LD between the sat insertion and lead SNP means insertion-causal vs
SNP-causal cannot be split without denser or functional data; populations 1/2/4 are
near-fixed for these alleles, so only 3/6 are testable for the top loci.
