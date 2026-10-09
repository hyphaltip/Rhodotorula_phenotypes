# gwas_dh4148

GWAS, dose-0 consistency, growth rate and IC50 on the DH4148 reference (Cr, Cu, Pb; 126-strain panel). Report: `report/REPORT.md` and `report/REPORT.pdf`.

Run order (SLURM, from the repo root; `pixi run python` for Python and R):
1. `scripts/01_phenotypes.py` strain traits, growth rate, IC50
2. `scripts/02_genotypes.sh` panel VCF, BIMBAM genotypes, annotation table
3. `scripts/03_dose0_models.R`, `04_dose0_plots.py` dose-0 consistency
4. `scripts/05_make_pheno.py`, `10_lineages_and_adjust.py` GWAS phenotypes, lineages, run adjustment
5. `scripts/06_run_gemma.sh PHENO OUTDIR [COVAR]` (needs `module load gemma`), `07a_prune.sh`
6. `scripts/07_gwas_summary.py OUTDIR PHENO TAG [lineage]`, `11_lineage_marker_check.py`, `08_nonadditive.py`, `09_growth_ic50_plots.py`
7. `scripts/13a_lineage_vcf.sh`, `13b_recombination.py`, `14_cr_tree_signal.py`, `15_candidate_panels.py`, `16_within_lineage_scan.py`
8. `scripts/12_make_report.py` (text in `scripts/report_text.md`, tables in `scripts/report_tables.py`), then pandoc.

Large files (panel VCF, GEMMA output) are not tracked. Main finding: the panel has 14 clonal lineages (D-61).
