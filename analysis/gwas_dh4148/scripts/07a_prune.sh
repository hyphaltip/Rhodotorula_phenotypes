#!/usr/bin/env bash
# LD-pruned SNP sets of the panel VCF (used for the effective number of tests and for the epistasis scan).
set -euo pipefail
module load bcftools
G=analysis/gwas_dh4148/results/geno; W=${SCRATCH:?}
for r in 0.5 0.2; do
  bcftools +prune -m $r -w 100kb $G/panel.vcf.gz -Oz -o $W/pruned_$r.vcf.gz
  bcftools query -f '%CHROM:%POS:%REF:%ALT\n' $W/pruned_$r.vcf.gz > $G/pruned_r2_$r.snps.txt
  echo "r2 $r: $(wc -l < $G/pruned_r2_$r.snps.txt) SNPs"
done
