#!/usr/bin/env bash
# Panel genotypes for GEMMA (BIMBAM) from popgen's rmuc_core SNP VCF (biallelic SNPs, MAF >= 0.05 in 247 strains).
# Panel = strain_curation.gwas_panel (126 strains, D-49). Run as a SLURM job from the repo root; heavy I/O goes to $SCRATCH.
set -euo pipefail
module load bcftools
OUT=analysis/gwas_dh4148/results/geno; mkdir -p $OUT
V=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref/results/RmucDH4148.rmuc_core.snps.maf.annotated.vcf.gz
W=${SCRATCH:?}
/usr/bin/python3.12 - <<'P'
import pandas as pd
c = pd.read_csv("data/metadata/strain-curation/strain_curation.csv", dtype=str).fillna("")
p = c[c.gwas_panel == "True"]
assert len(p) == 126 and (p.popgen_strain != "").all()
# VCF sample ids are written NAME_NAME
open("analysis/gwas_dh4148/results/geno/panel_vcf_samples.txt", "w").write("\n".join(s + "_" + s for s in p.popgen_strain) + "\n")
p[["strain_id", "sample_name", "popgen_strain"]].to_csv("analysis/gwas_dh4148/results/geno/panel_strains.csv", index=False)
P
bcftools view -S $OUT/panel_vcf_samples.txt --force-samples -m2 -M2 -v snps $V -Ou \
 | bcftools +fill-tags -Ou -- -t AF,AC,AN \
 | bcftools view -i 'AN>=113 && AF>=0.05 && AF<=0.95' -Oz -o $W/panel.vcf.gz
bcftools index -t $W/panel.vcf.gz
cp $W/panel.vcf.gz $W/panel.vcf.gz.tbi $OUT/
echo "panel sites: $(bcftools view -H $W/panel.vcf.gz | wc -l); samples: $(bcftools query -l $W/panel.vcf.gz | wc -l)"
# BIMBAM mean genotype (0/1, haploid; missing = NA) and SNP annotation
bcftools query -l $W/panel.vcf.gz | sed 's/^\(.*\)_\1$/\1/' > $OUT/panel_samples_in_vcf_order.txt
bcftools query -f '%CHROM:%POS:%REF:%ALT,%ALT,%REF[,%GT]\n' $W/panel.vcf.gz | sed 's/,\./,NA/g' > $W/panel.bimbam
cp $W/panel.bimbam $OUT/ && gzip -f $OUT/panel.bimbam
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AN\t%INFO/ANN\n' $W/panel.vcf.gz | gzip > $OUT/panel_sites_ann.tsv.gz
