#!/usr/bin/env bash
# Population-vs-locus disambiguation (GWAS.md section 13's battery) for the 4
# FDR-significant loci found for cu_doseauc_v0151 (section 18) -- 5 SNPs since
# scaffold_2 has two independent hits ~1Mb apart. check_population_vs_locus.py
# processes one (trait, top_snp) row per invocation cleanly but writes its
# per-population detail file as within_pop_<trait>_<panel>.csv, which would
# collide across loci sharing one trait name -- so this driver runs it once
# per locus with its own --out/--work-dir, then concatenates the 5 one-row
# results.
#
# Run from the repo root, inside a SLURM job (needs $SCRATCH):
#   bash analysis/gwas/scripts/run_population_vs_locus_cu_v0151.sh
set -euo pipefail

REPO_ROOT="$PWD"
[ -f "$REPO_ROOT/pixi.toml" ] || { echo "ERROR: run from repo root." >&2; exit 1; }
: "${SCRATCH:?SCRATCH env var not set -- run inside a SLURM job}"

TOOLCHAIN="/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/GWAS/vc2gwas_env/bin"
VCF="data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.SNP.combined_selected.vcf.gz"
PHENO_CSV="analysis/gwas/results/gwas_next_phenotypes_fam_order.csv"
POP_CSV="analysis/gwas/data/prior_run_state/pop_assignment_at_run.csv"
CULLED_LIST="analysis/gwas/results/gwas/near_clone_culling/culled_keep.txt"
PRUNED_GWAS="analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas.pruned"
PRUNED_GWASC="analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship_culled/gwasc.pruned"
ASSOC_DIR="analysis/gwas/results/gwas/tierA_summary/assoc_csv"
FDR_DIR="analysis/gwas/results/gwas/fdr"
OUTDIR="analysis/gwas/results/gwas/population_vs_locus"
WORK="$SCRATCH/pop_vs_locus_cu_v0151"
mkdir -p "$WORK" "$OUTDIR"

for f in "$VCF" "$PHENO_CSV" "$POP_CSV" "$CULLED_LIST" \
         "$PRUNED_GWAS.bed" "$PRUNED_GWASC.bed" \
         "$ASSOC_DIR/gwas_cu_doseauc_v0151_assoc.csv.gz" "$FDR_DIR/gwas_cu_doseauc_v0151_fdr05.csv"; do
  [ -f "$f" ] || { echo "ERROR: required input missing: $f" >&2; exit 1; }
done

# trait=cu_doseauc_v0151 in every row (must match the actual phenotype column
# name -- check_population_vs_locus.py looks it up as pheno_df[trait]);
# panel=gwas, forced population_confound_risk=True since these are being
# tested by explicit request (D-24), not because they passed the crude
# AF-swing screen.
LOCI="scaffold_9:704260:C:T scaffold_16:455499:C:T scaffold_8:698507:T:C scaffold_2:515984:G:A scaffold_2:1560553:T:A"

RESULTS=()
i=0
for rs in $LOCI; do
  i=$((i + 1))
  LOCUS_WORK="$WORK/locus_$i"
  LOCUS_OUT="$OUTDIR/population_vs_locus_cu_doseauc_v0151_locus${i}.csv"
  mkdir -p "$LOCUS_WORK"
  CONFOUND_CSV="$LOCUS_WORK/confounding_gwas.csv"
  printf 'trait,top_snp,n_fdr05,population_confound_risk\ncu_doseauc_v0151,%s,32,True\n' "$rs" > "$CONFOUND_CSV"

  echo "=== locus $i: $rs ==="
  pixi run python3 analysis/gwas/scripts/check_population_vs_locus.py \
    --confounding-csv "$CONFOUND_CSV" \
    --pheno-csv "$PHENO_CSV" \
    --pop-csv "$POP_CSV" \
    --culled-list "$CULLED_LIST" \
    --vcf "$VCF" \
    --bcftools "$TOOLCHAIN/bcftools" \
    --plink2 "$TOOLCHAIN/plink2" \
    --pruned-bfile-gwas "$PRUNED_GWAS" \
    --pruned-bfile-gwasc "$PRUNED_GWASC" \
    --assoc-csv-dir "$ASSOC_DIR" \
    --fdr-dir "$FDR_DIR" \
    --work-dir "$LOCUS_WORK" \
    --out "$LOCUS_OUT"
  RESULTS+=("$LOCUS_OUT")

  # rescue the per-population detail file (written under a fixed name that
  # would be overwritten by the next locus) into a locus-specific copy
  DETAIL="$OUTDIR/within_pop_cu_doseauc_v0151_gwas.csv"
  if [ -f "$DETAIL" ]; then
    cp "$DETAIL" "$OUTDIR/within_pop_cu_doseauc_v0151_locus${i}_gwas.csv"
  fi
done

echo "Merging ${#RESULTS[@]} per-locus results ..."
pixi run python3 -c "
import pandas as pd
frames = [pd.read_csv(f) for f in [$(printf '\"%s\",' "${RESULTS[@]}")]]
out = pd.concat(frames, ignore_index=True)
out.to_csv('$OUTDIR/population_vs_locus_cu_doseauc_v0151.csv', index=False)
print(out[['trait','panel','top_snp','meta_p','meta_p_fdr','partial_r2','null_percentile','verdict']].to_string(index=False))
"
echo "Wrote $OUTDIR/population_vs_locus_cu_doseauc_v0151.csv"
