# GWAS: Color and Copper-Response Phenotypes in *Rhodotorula mucilaginosa*

**Status**: in progress (port from `analysis/ideas/2026-08-15-color-phenotype-space/`)
**Spec**: `docs/superpowers/specs/2026-08-25-gwas-port-design.md`

This analysis ports the Tier A-G GWAS pipeline (single-SNP kinship-only LMM, SKAT/burden
set tests, BSLMM architecture, LOCO sensitivity, gene annotation, fine-mapping, prior-locus
replication) originally developed in `analysis/ideas/2026-08-15-color-phenotype-space/`
(see that folder's `GWAS_REPORT.md` and `PROGRESS.md` for the full original narrative and
`.living/decisions.md` D-9 through D-13 for the modeling decisions carried forward
unchanged here).

New in this port (not done rigorously in the original run):
1. Audited strain-name reconciliation (`scripts/reconcile_strains.py`) — see
   `results/strain_reconciliation/`.
2. Per-strain ploidy/heterozygosity validation (`scripts/check_ploidy.py`) — see
   `results/ploidy_check/`.

<!-- Sections below filled in by Task 8 -->
## Strain reconciliation
## Ploidy validation
## Strain-state diff vs. prior run
## GWAS pipeline results
