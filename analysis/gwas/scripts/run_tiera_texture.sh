#!/usr/bin/env bash
# Tier A add-on: 13 raw Haralick/GLCM texture metrics (control condition, Cu=0,
# same late window as the color traits, clone-mean over plates) as 13
# additional traits -- added 2026-09-11 to test a PI hypothesis from the
# sibling Rhodotorula_Metabolites project (AHL autoinducer production vs.
# colony morphology; smooth/rough is thought to be a proxy for capsule
# production). Same full unpruned SNP set + per-panel kinship as
# run_tiera_gemma.sh / run_tiera_lab_traits.sh -- skip-if-present so this is
# safe to rerun.
#
# Run from the repo root, inside a SLURM job (needs $SCRATCH):
#   sbatch --partition=stajichlab --time=10:00:00 --wrap="bash analysis/gwas/scripts/run_tiera_texture.sh"
set -euo pipefail

REPO_ROOT="$PWD"
[ -f "$REPO_ROOT/pixi.toml" ] || { echo "ERROR: run from repo root." >&2; exit 1; }

: "${SCRATCH:?SCRATCH env var not set -- run inside a SLURM job}"
WORK="$SCRATCH/gwas_tiera_texture"
mkdir -p "$WORK"

PHENO_CSV="analysis/gwas/results/gwas_next_phenotypes_fam_order.csv"
TOOLCHAIN="/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/GWAS/vc2gwas_env/bin"
OUTDIR="analysis/gwas/results/gwas/tierA_summary"
mkdir -p "$OUTDIR/gemma_output"

TRAITS="tex_AngularSecondMoment tex_Contrast tex_Correlation tex_HaralickVariance tex_InverseDifferenceMoment tex_SumAverage tex_SumVariance tex_SumEntropy tex_Entropy tex_DiffVariance tex_DiffEntropy tex_InfoCorrelation1 tex_InfoCorrelation2"

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

echo "Texture Tier A scan complete: 13 traits x 2 panels."
