#!/usr/bin/env bash
# Reproduce the GWAS port end to end. Run from the repo root: bash analysis/gwas/run.sh
# NOTE: Tasks 3-4's human-review steps (NEEDS_REVIEW.md adjudication, ploidy-flag
# keep/exclude decisions) are NOT automated here -- this script assumes
# strain_match_table.reviewed.csv already reflects those decisions, and that the
# ploidy keep/exclude call has already been made (see GWAS.md S2).
set -euo pipefail
# Run from the repo root (or set GWAS_REPO_ROOT). No BASH_SOURCE: it breaks on SLURM.
cd "${GWAS_REPO_ROOT:-$PWD}"
test -f analysis/gwas/run.sh || { echo "run from the repo root or set GWAS_REPO_ROOT" >&2; exit 1; }
module load bcftools

# Strain table now comes from DuckDB (strain_info view), not data/metadata/Copper.Strain_info.csv
# (removed in D-35). See analysis/gwas/DUCKDB_RETOOL_PLAN.md. The export opens the DB read-only
# (or a $SCRATCH snapshot if write-locked), so run this script inside a SLURM job:
#   sbatch -p short -c 2 --mem=8G -t 30 --wrap="bash analysis/gwas/run.sh"
pixi run python3 analysis/gwas/scripts/duckdb_inputs.py \
  --out analysis/gwas/results/strain_reconciliation_duckdb/strain_info_from_duckdb.csv

pixi run python3 analysis/gwas/scripts/reconcile_strains.py \
  --pheno analysis/gwas/results/strain_reconciliation_duckdb/strain_info_from_duckdb.csv \
  --pheno-col Strain \
  --vcf data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.SNP.combined_selected.vcf.gz \
  --out-dir analysis/gwas/results/strain_reconciliation_duckdb
# NOTE: output goes to strain_reconciliation_duckdb/ so the existing, human-reviewed
# strain_reconciliation/strain_match_table.reviewed.csv is NOT overwritten. Diff the two, then
# carry forward decisions by hand.

echo "STOP: review analysis/gwas/results/strain_reconciliation_duckdb/NEEDS_REVIEW.md and"
echo "produce strain_match_table.reviewed.csv before continuing."
echo ""
echo "Remaining steps (run manually after each human-review gate -- see GWAS.md):"
echo "  1. pixi run python3 analysis/gwas/scripts/check_ploidy.py --strains .../strain_match_table.reviewed.csv --vcf ... --mosdepth-dir ... --cram-dir ... --out analysis/gwas/results/ploidy_check/ploidy_flags.csv"
echo "  2. Review ploidy_flags.csv keep/exclude decision (never auto-excluded)."
echo "  3. pixi run python3 analysis/gwas/scripts/diff_strain_state.py --reviewed ... --exclude <ploidy-excluded ids> --prior-fam-all ... --prior-fam-culled ... --prior-pop ... --current-pop ... --out analysis/gwas/results/strain_state_diff/state_diff_report.json"
echo "  4. If overall_action == rebuild: sbatch --wrap='bash analysis/gwas/scripts/rebuild_tiers_abc.sh' (inside a SLURM job, needs \$SCRATCH)"
echo "  5. Read grm_diagnostic.json's verdict -- if singular_risk, investigate before trusting p-values (see GWAS.md S4)."
echo "  6. phenotypes: build_gwas_phenotypes.py reads the removed db_extract.parquet. Use build_gwas_phenotypes_duckdb.py (writes to results/duckdb_retool/, approximate colour/texture columns; see DUCKDB_RETOOL_PLAN.md)"
echo "  7. sbatch --wrap='bash analysis/gwas/scripts/run_tiera_gemma.sh' (inside a SLURM job)"
echo "  8. pixi run python3 analysis/gwas/scripts/summarize_tiera.py --assoc-dir analysis/gwas/results/gwas/tierA_summary/gemma_output --out-summary analysis/gwas/results/gwas/tierA_summary/tiera_summary.csv --out-fdr-dir analysis/gwas/results/gwas/fdr"
echo "  (Tier B/C/LOCO/pixy rebuild: not yet scripted -- see GWAS.md 'Next steps')"
