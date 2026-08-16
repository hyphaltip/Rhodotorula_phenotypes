<!-- Provenance for the Rhodotorula PHYling protein tree dataset -->
# Provenance: rhodotorula-phyling-protein-tree

## Source

**Type**: computed phylogenetic tree (multi-locus protein concatenation, BUSCO marker set)

**Origin**:
- Pipeline: PHYling `protein_tree` workflow (BUSCO `fungi_odb10` markers, concatenated alignment `protein-Rhodotorula-taxa_278.fungi_odb10.fa`) → FastTree v2.2.0 (LG/CAT, 20 rate categories, SH-like 1000 supports).
- **External source directory (do not delete/relocate — this is the user's shared filesystem tree)**:
  `/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS/BFD/results/phyling_pep/protein/buildtree/fungi_odb10/fasttree/`
- Source files copied into this repo (originals left in place as provenance):
  - `protein-Rhodotorula-taxa_278.fungi_odb10.fasttree.support.treefile` (21,129 B)
  - `protein-Rhodotorula-taxa_278.fungi_odb10.fasttree.nosupport.treefile` (19,809 B)
  - `protein-Rhodotorula-taxa_278.fungi_odb10.fasttree.support.log` (8,697 B)

**Citation / accession**: N/A (internal PHYling run; BFD = "build from database"). BUSCO/FastTree/LG2008 — see log header.

## Acquisition details

**Date acquired**: 2026-08-15 (copied into this repo from the shared source dir above).

**Obtained by**: Jason Stajich lab / PHYling protein_tree pipeline (external, pre-existing).

**Method**: 278 Rhodotorula taxa + outgroups; FastTree 2.2.0 double precision, OpenMP 12 threads; BLOSUM45 joins, balanced; Support SH-like 1000; Search Normal +NNI +SPR (2 rounds range 10) + ML-NNI opt-each=1; ML Model Le-Gascuel 2008, CAT approximation with 20 rate categories. Total branch length 5.637; final LogLk = -1679594.151; runtime ~1116 s; 267/278 unique seqs, 1/264 bad splits (worst delta-LogLk 0.265). Tip labels end in `.proteins` (alignment-derived); see `FINDINGS/idea09` analysis + `scripts/idea_09_phylogeny.R` for how the `.proteins`/`.proteins.fa` suffix is stripped before joining to strain metadata.

**Checksum**: md5 `704de78b2689978279870256b8b416dc` (support treefile), `9b32a896c8fcef6001d361dcd4acabf2` (nosupport treefile).

## Access restrictions

**Restriction level**: none (repo internal; no data-use agreement in place)

**Details**: None.

## Known issues

- **Tip labels end in `.proteins`** (no `.fa` in this naming scheme; user noted some trees may use `.proteins.fa`). Analysis scripts strip both suffixes defensively.
- **Outgroup tips**: 276 of 278 taxa are Rhodotorula; tips `Cystobasidium` and `Pseudomicrostroma` are outgroups, excluded from trait-join.
- **Strain code mismatch**: tip suffix `DH4148` did not match any strain in `strain_metadata.tsv` (277/278 tips matched). See idea 09 analysis.
- Tree reconstructed on protein alignment; branch lengths in substitutions/site, not time.

## Contact

**Primary contact**: Jason Stajich (jstajich)

**Backup contact**: None.

## Version history

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-08-15 | Initial ingestion: support + nosupport treefiles + FastTree log, metadata (provenance, schema, summary) |
| <!-- 1.1 --> | <!-- YYYY-MM-DD --> | <!-- e.g., "Added a time-calibrated tree." --> |
