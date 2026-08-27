Diamond blastp databases used to produce `blast/*.tsv` were built on the fly from
external shared proteome FASTAs and are NOT retained in this repo (60MB of
regeneratable binary index). To reproduce:

```
module load diamond/2.1.7
DB_SRC=/bigdata/stajichlab/shared/projects/Rhodotorula/MAT_search/db
for f in "$DB_SRC"/*.proteins.fa; do
  name=$(basename "$f" .proteins.fa)
  diamond makedb --in "$f" -d "dbs/$name"
done
# S. cerevisiae (SGD R64)
zcat /bigdata/gen220/shared/data-examples/examples/Saccharomyces_cerevisiae.peps.fa.gz \
  > /tmp/Saccharomyces_cerevisiae.peps.fa
diamond makedb --in /tmp/Saccharomyces_cerevisiae.peps.fa -d dbs/Saccharomyces_cerevisiae
```

Then: `python3 ../../scripts/build_homolog_impact_table.py` (run from `analysis/gwas/`,
after the blastp step below).

```
diamond blastp -q ../sequences/candidate_genes_protein.fa -d dbs/<name>.dmnd \
  --outfmt 6 qseqid sseqid pident length evalue bitscore qstart qend sstart send qlen slen qseq sseq stitle \
  --max-target-seqs 1 -o blast/<name>.tsv
```

Species included: 8 Rhodotorula spp. + 2 Cystobasidium outgroups (the same 11-genome
set as `algorithms/functional_annotation/`'s pigment-HMM scan, minus the R. mucilaginosa
query genome which is instead searched as a same-genome paralog control) + S. cerevisiae.
