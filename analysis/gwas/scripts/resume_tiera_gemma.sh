#!/usr/bin/env bash
# Resume Tier A GEMMA scans (full unpruned SNP set) for whichever (panel,trait)
# combinations are still missing from tierA_summary/gemma_output/ -- companion
# to rebuild_full_genotypes_and_tiera.sh, used if that job's Tier A loop is
# killed by SLURM walltime before finishing all 24 scans. Does NOT regenerate
# genotype sets or kinship -- assumes both panels' persisted bfiles already
# exist under results/gwas/grm_conditioning/{rebuilt_kinship,rebuilt_kinship_culled}/.
#
# Run from the repo root, inside a SLURM job (needs $SCRATCH):
#   sbatch --time=3:00:00 --wrap="bash analysis/gwas/scripts/resume_tiera_gemma.sh"
set -euo pipefail

REPO_ROOT="$PWD"
[ -f "$REPO_ROOT/pixi.toml" ] || { echo "ERROR: run from repo root." >&2; exit 1; }

: "${SCRATCH:?SCRATCH env var not set -- run inside a SLURM job}"
WORK="$SCRATCH/gwas_tiera_resume"
mkdir -p "$WORK"

PHENO_CSV="analysis/gwas/results/gwas_next_phenotypes_fam_order.csv"
TOOLCHAIN="/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/GWAS/vc2gwas_env/bin"
OUTDIR="analysis/gwas/results/gwas/tierA_summary"
mkdir -p "$OUTDIR/gemma_output"

TRAITS="chroma sat bright clone_mean_area AUC_0 AUC_10 AUC_20 AUC_30 AUC_ratio_10 resilience_30 cu_dose_slope IC50_est"

run_tiera_panel() {
  local panel="$1" bfile="$2" kinship="$3"
  for trait in $TRAITS; do
    outfile="$OUTDIR/gemma_output/${panel}_${trait}.assoc.txt"
    if [ -s "$outfile" ]; then
      echo "skip [$panel/$trait]: already present ($outfile)"
      continue
    fi
    echo "=== Tier A [$panel]: $trait ==="
    TWORK="$WORK/${panel}_${trait}"
    mkdir -p "$TWORK"
    cp "$bfile.bed" "$TWORK/g.bed"
    cp "$bfile.bim" "$TWORK/g.bim"
    pixi run python3 -c "
import pandas as pd
fam = pd.read_csv('$bfile.fam', sep=r'\s+', header=None, names=['fid','iid','pid','mid','sex','pheno'])
pheno = pd.read_csv('$PHENO_CSV')
pheno_map = dict(zip(pheno['strain_code'], pheno['$trait']))
vals = fam['iid'].map(pheno_map)
n_missing = vals.isna().sum()
fam['pheno'] = vals.apply(lambda v: 'NA' if pd.isna(v) else v)
fam.to_csv('$TWORK/g.fam', sep=' ', header=False, index=False)
print(f'  [$panel] $trait: {len(fam) - n_missing}/{len(fam)} with phenotype, {n_missing} NA')
"
    "$TOOLCHAIN/gemma" -bfile "$TWORK/g" -k "$kinship" -lmm 4 \
      -o "${panel}_${trait}" -outdir "$OUTDIR/gemma_output" \
      > "$OUTDIR/gemma_output_${panel}_${trait}.log" 2>&1 \
      || { echo "ERROR: GEMMA failed for $panel/$trait" >&2; exit 1; }
    rm -rf "$TWORK"
  done
}

run_tiera_panel gwas \
  "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas" \
  "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/kins.cXX.txt"
run_tiera_panel gwasc \
  "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship_culled/gwasc" \
  "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship_culled/kinsc.cXX.txt"

echo "Tier A resume complete."
