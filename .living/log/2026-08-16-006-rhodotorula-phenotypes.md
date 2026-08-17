# Session 2026-08-16-006 — GWAS run + pixy diversity pilot

**Date**: 2026-08-16
**Work dir**: `/scratch/jstajich/27384933/gwas/work` (results in `results/gwas/`)
**Docs**: `00 analysis/ideas/2026-08-15-color-phenotype-space/{GWAS_REPORT,PROGRESS}.md`,
`08b-quantitative-geneticist-review.md`

## What was done
- **GWAS pipeline** (kinship-only LMM, valid): VCF subset to 201 strains → region
  exclusion (scaffold_21 + aneuploid windows) → biallelic → PLINK2 QC
  (404,706 vars × 201) → LD-pruned kinship (20,769) → GEMMA `-lmm 4` per trait.
  Results: lambda 0.357–0.638, PVE(size)=0.180±0.077, top size hit scaffold_16:121473
  p=6.2e-11 (13 independent loci). Fixed GEMMA gotchas (fam col-6 phenotype, integer
  chrom codes) → learnings L-16.
- **Population-structure robustness runs** (10 PCs, then 3 PCs, then pop dummies) all
  collapsed (PVE=0.99997, se=NaN, GSL singular) from the singular near-clone GRM.
  **Decision: kinship-only LMM is the valid model; PC-covariate runs parked.**
  → learning L-17, decision D-9.
- **Expert quantitative-geneticist review** launched → saved `08b-...review.md` with
  power-retention strategies (clone-mean phenotyping, set/SKAT tests, LOCO kinship,
  Meff/permutation threshold, bivariate, BSLMM architecture).
- **Pixy pilot (scaffold_1)**: built joint-genotyped 201-strain cohort VCF with
  invariant sites (GATK GenotypeGVCFs -all-sites), persisted to `results/gwas/pixy/`;
  computed pi/dxy/fst(Hudson)/watterson_theta. Results strongly population-structured
  (mean Fst 0.494; pi varies 600× across pops) → finding #12.
  Fixed SLURM/HPC gotchas (module loading, $SCRATCH, sbatch -D cwd, pixy keyword,
  pixy fst_type) → learning L-18.

## Decisions
- Park PC-covariate LMM; kinship-only LMM is primary (D-9).

## Next steps
- Full-genome pixy (all 23 scaffolds) to map trait-relevant divergence.
- Optional power-retention analyses from expert review (clone-mean phenotyping,
  SKAT/gene-based, LOCO kinship, permutation threshold).
- Copy final GWAS plots/outputs into repo results; finalize report.

## Tags
gwas, gemma, LMM, kinship, pixy, pi, dxy, fst, population-structure, near-clone,
singular, HPC, slurm, idea11