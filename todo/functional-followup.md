# Functional follow-up of Tier A–E GWAS candidate genes

## Status
open (priority: medium), added 2026-08-17

## Why
Tier D gene mapping (GWAS_REPORT §9.1, `tierD_*_annotated.csv`) assigned several Tier-A
anchors and BSLMM peaks to biologically plausible or novel genes that deserve experimental
follow-up:

| Locus (trait) | Gene(s) | Reasoning |
|---------------|---------|-----------|
| chr13 block (AUC_10; replicates prior growth-rate locus, p_wald=4.0e-6) | OM429_005439 (hypothetical, scaffold_13:30,096–30,670) | Replication locus at the exact prior hit; no functional annotation |
| scaffold_10:396172 (AUC_10) | OM429_004640 / DBP3 (RNA-dependent ATPase) | Rare-variant signal (CS n=67, af=0.015); ATPase plausibly copper/energy related |
| scaffold_8:831789 (chroma) | OM429_004009 telomerase RT | Chromosome-stability / pigment candidate |
| scaffold_3 BSLMM peak (chroma) | OM429_001379 methionine aminopeptidase 1 | Top-PIP cluster; protein-turnover link to pigment maturation |

## Proposed work
1. Expression / RT-qPCR survey of the candidate genes across representative strains
   (low vs high phenotype) under control and copper stress.
2. Targeted sequencing or allele-confirmation of the lead variants in the n=201 panel.
3. If a culturable phenotype assay is feasible (e.g., DBP3 ATPase inhibitors, pigment
   kinetics for telomerase/metAP), test allele-effect predictions from Tier E.
4. Annotate OM429_005439 in the reference annotation (structure/AI function) to break the
   hypothetical-protein barrier.

## Depends on
- Wet-lab capacity; no dependency on code.

## Success criteria
- Expression or phenotypic difference consistent with the GWAS allele direction for at
  least one candidate; a functional hypothesis for OM429_005439.
