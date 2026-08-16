# Data Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

### rhodotorula-phyling-protein-tree
```yaml
name: rhodotorula-phyling-protein-tree
type: phylogenetic-tree
source: PHYling protein_tree (BUSCO fungi_odb10) -> FastTree 2.2.0 (LG/CAT); copied from user's shared BFD results dir (see provenance.md)
annotation_sources: analysis/ideas/2026-08-15-color-phenotype-space/data/strain_metadata.tsv (tip-strain join by code)
date_acquired: 2026-08-15
format: Newick (support + nosupport treefiles) + FastTree log
rows: 278 tips / 265 internal nodes
columns: n/a (tree topology + branch lengths)
size: 49.6 KB (3 files)
raw_path: data/raw/rhodotorula-phyling-protein-tree/
metadata_path: data/metadata/rhodotorula-phyling-protein-tree/
status: raw (immutable) + analyzed (idea 09)
known_issues:
  - Tip labels end in .proteins (strip .proteins/.proteins.fa before strain join)
  - 2 outgroup tips (Cystobasidium, Pseudomicrostroma)
  - Tip DH4148 has no matching strain in strain_metadata.tsv
access_restrictions: none
tags: [rhodotorula, phylogeny, phyling, busco, fungi_odb10, fasttree, protein-tree, strains]
```

Maximum-likelihood tree of 278 Rhodotorula-related taxa (276 Rhodotorula + 2 outgroups)
built by the PHYling `protein_tree` workflow and FastTree 2.2.0. Copied from the shared
directory `/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS/BFD/results/phyling_pep/protein/buildtree/fungi_odb10/fasttree/`.
Used by idea 09 (phylogeneticist) to test for phylogenetic signal in strain-level color/
growth phenotypes. See `provenance.md` for full source path and reconstruction settings.

### copper-colony-measurements
```yaml
name: copper-colony-measurements
type: imaging
source: automated time-course colony imaging (imager runs d000353-d000357) + segmentation export
annotation_sources: data/metadata/Copper.Strain_info.csv, data/metadata/Copper.Plate_info.csv
date_acquired: 2026-08-14
format: Parquet (2,398 files) -> DuckDB colony_measurement (181 cols)
rows: 211800
columns: 181
size: 428 MB (preprocessed parquets)
raw_path: data/preprocessed/
metadata_path: data/metadata/copper-colony-measurements/
db_table: colony_measurement
status: processed
known_issues:
  - Variable object count per image (wells not always fully detected)
  - Metadata_Dataset / Metadata_ImageName source columns dropped on import
  - Plate->run/configuration map is an implicit Copper convention applied at import
access_restrictions: none
tags: [copper, colony, rhodotorula, morphology, time-course, segmentation]
```

Per-colony measurements for the Copper exposure phenotyping experiment (runs 353-357,
temperature token 300). One row per segmented colony object per image; 211,800 rows across
2,398 images, keyed by `image_name` + `object_label`. Feature families: Shape, Intensity,
Haralick Texture (gray), Bbox, and color (xy / CIELAB / HSV). Strain/plate annotation lives
in `Copper.Strain_info.csv` and `Copper.Plate_info.csv` (see provenance.md). Fully loaded into
the DuckDB `colony_measurement` table by `scripts/db/`.
