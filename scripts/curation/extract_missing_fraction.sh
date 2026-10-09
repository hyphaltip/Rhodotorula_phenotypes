#!/usr/bin/env bash
# Per-strain fraction of missing genotypes over all SNP sites of popgen's rmuc_core VCF (the input popgen's declone.py uses).
# Popgen's representative rule needs it ("fewest missing genotypes") but only publishes it for strains in near-identical groups.
# Output is a small derived table, stored with the popgen snapshot. Run as a SLURM job from the repo root; do not rely on BASH_SOURCE.
set -euo pipefail
module load bcftools
VCF=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref/results/RmucDH4148.rmuc_core.qc.annotated.vcf.gz
OUT=data/raw/popgen-callset-metadata/popgen/variant_qc/rmuc_core.missing_by_strain.tsv
{ bcftools query -l "$VCF" | paste -sd '\t'
  bcftools view -v snps -Ou "$VCF" | bcftools query -f '[%GT\t]\n'
} | /usr/bin/python3.12 -c '
import sys
names = sys.stdin.readline().rstrip("\n").split("\t"); n = len(names); miss = [0] * n; sites = 0
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")[:n]; sites += 1
    for i, g in enumerate(f):
        if "." in g: miss[i] += 1
print("strain\tn_snp_sites\tn_missing\tmissing_frac")
def clean(nm):
    h = len(nm) // 2
    return nm[:h] if len(nm) % 2 == 1 and nm[:h] == nm[h + 1:] else nm  # VCF sample ids are written NAME_NAME
for nm, m in zip(names, miss): print(f"{clean(nm)}\t{sites}\t{m}\t{m / sites:.6f}")
' > "$OUT"
wc -l "$OUT"
