#!/usr/bin/env bash
# Reproduce the GWAS port end to end. Run from the repo root: bash analysis/gwas/run.sh
# NOTE: Tasks 3-4's human-review steps (NEEDS_REVIEW.md adjudication, ploidy-flag
# keep/exclude decisions) are NOT automated here -- this script assumes
# strain_match_table.reviewed.csv already reflects those decisions, and that the
# ploidy keep/exclude call has already been made (see GWAS.md S2).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
module load bcftools

pixi run python3 analysis/gwas/scripts/reconcile_strains.py \
  --pheno data/metadata/Copper.Strain_info.csv \
  --pheno-col Strain \
  --vcf data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.SNP.combined_selected.vcf.gz \
  --out-dir analysis/gwas/results/strain_reconciliation

echo "STOP: review analysis/gwas/results/strain_reconciliation/NEEDS_REVIEW.md and"
echo "produce strain_match_table.reviewed.csv before continuing."
echo ""
echo "Remaining steps (run manually after each human-review gate -- see GWAS.md):"
echo "  1. pixi run python3 analysis/gwas/scripts/check_ploidy.py --strains .../strain_match_table.reviewed.csv --vcf ... --mosdepth-dir ... --cram-dir ... --out analysis/gwas/results/ploidy_check/ploidy_flags.csv"
echo "  2. Review ploidy_flags.csv keep/exclude decision (never auto-excluded)."
echo "  3. pixi run python3 analysis/gwas/scripts/diff_strain_state.py --reviewed ... --exclude <ploidy-excluded ids> --prior-fam-all ... --prior-fam-culled ... --prior-pop ... --current-pop ... --out analysis/gwas/results/strain_state_diff/state_diff_report.json"
echo "  4. If overall_action == rebuild: sbatch --wrap='bash analysis/gwas/scripts/rebuild_tiers_abc.sh' (inside a SLURM job, needs \$SCRATCH)"
echo "  5. Read grm_diagnostic.json's verdict -- if singular_risk, investigate before trusting p-values (see GWAS.md S4)."
echo "  6. pixi run python3 analysis/gwas/scripts/build_gwas_phenotypes.py"
echo "  7. sbatch --wrap='bash analysis/gwas/scripts/run_tiera_gemma.sh' (inside a SLURM job)"
echo "  8. pixi run python3 analysis/gwas/scripts/summarize_tiera.py --assoc-dir analysis/gwas/results/gwas/tierA_summary/gemma_output --out-summary analysis/gwas/results/gwas/tierA_summary/tiera_summary.csv --out-fdr-dir analysis/gwas/results/gwas/fdr"
echo "  (Tier B/C/LOCO/pixy rebuild: not yet scripted -- see GWAS.md 'Next steps')"
