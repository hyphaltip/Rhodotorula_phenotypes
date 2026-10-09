#!/usr/bin/env bash
# GEMMA univariate LMM (Wald test) per trait. Usage: 06_run_gemma.sh PHENO_FILE OUT_DIR [COVARIATE_FILE]
# The centred relatedness matrix is built once from all panel SNPs. Run as a SLURM job from the repo root.
set -euo pipefail
module load gemma
G=analysis/gwas_dh4148/results/gwas; B=analysis/gwas_dh4148/results/geno/panel.bimbam.gz
PH=$1; OUT=$G/$2; CV=${3:-}
mkdir -p $OUT
[ -s $G/output/K.cXX.txt ] || { mkdir -p $G/output; gemma -g $B -p $PH -gk 1 -n 1 -o K -outdir $G/output > $G/output/K.stdout 2>&1; }
N=$(wc -l < $G/trait_list.csv); N=$((N-1))
for i in $(seq 1 $N); do
  t=$(sed -n "$((i+1))p" $G/trait_list.csv | cut -d, -f2)
  C=""; [ -n "$CV" ] && C="-c $CV"
  gemma -g $B -p $PH -k $G/output/K.cXX.txt $C -lmm 1 -n $i -maf 0.05 -miss 0.1 -o $t -outdir $OUT > $OUT/$t.stdout 2>&1
  gzip -f $OUT/$t.assoc.txt
done
echo done
