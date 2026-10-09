#!/usr/bin/env bash
# Genotypes of each large lineage (L01-L03) from popgen's rmuc_core QC VCF (all SNPs, not only MAF >= 0.05 in the 247 strains).
# Keeps biallelic SNPs with minor allele count >= 2 inside the lineage and <= 10% missing. Output: results/recomb/L0x.gt.tsv.gz (CHROM POS REF ALT, then one column per strain).
set -euo pipefail
: ${SCRATCH:?}
module load bcftools
V=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref/results/RmucDH4148.rmuc_core.qc.annotated.vcf.gz
O=analysis/gwas_dh4148/results/recomb; mkdir -p $O
for L in L01 L02 L03; do
  /usr/bin/python3.12 - $L <<'P'
import sys, pandas as pd
L = sys.argv[1]; d = pd.read_csv("analysis/gwas_dh4148/results/lineages.csv"); d = d[d.lineage == L]
open(f"analysis/gwas_dh4148/results/recomb/{L}.samples.txt", "w").write("\n".join(s + "_" + s for s in d.popgen_strain) + "\n")
P
  N=$(wc -l < $O/$L.samples.txt)
  bcftools view -S $O/$L.samples.txt --force-samples -m2 -M2 -v snps $V -Ou \
   | bcftools +fill-tags -Ou -- -t AC,AN \
   | bcftools view -i "AN>=$((N*9/10)) && AC>=2 && AC<=AN-2" -Oz -o $SCRATCH/$L.vcf.gz
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' $SCRATCH/$L.vcf.gz | sed 's/\t\./\tNA/g' | gzip > $O/$L.gt.tsv.gz
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' $SCRATCH/$L.vcf.gz | gzip > $O/$L.ann.tsv.gz
  echo "$L: $N strains, $(zcat $O/$L.gt.tsv.gz | wc -l) informative SNPs"
done
