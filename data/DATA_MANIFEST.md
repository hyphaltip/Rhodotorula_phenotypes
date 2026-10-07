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
status: removed
last_present_in: data-v1-pre-metal-replace
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

### RmucY2510_v2-genotypes
```yaml
name: RmucY2510_v2-genotypes
type: variant-calls
source: UCR Population_Genomics Rhodotorula mucilaginosa NRRL Y-2510 GATK hard-filtered VCFs (symlinked, not copied)
annotation_sources: data/metadata/Copper.Strain_info.csv (strain reconciliation), Rmuc_PopAssigned.csv (population groups, referenced by path)
date_acquired: 2026-08-25
format: bgzip VCF + tabix index (SNP + INDEL, biallelic-and-multiallelic mixed, haploid genotype calls)
rows: 422 samples; SNP VCF ~728,581 sites pre-QC (per prior GWAS run's Stage-0/1 filtering)
columns: n/a (VCF)
size: SNP ~662 MB, INDEL ~2.1 GB (symlinked, not duplicated on disk)
raw_path: data/raw/genotypes/RmucY2510_v2/
metadata_path: data/raw/genotypes/RmucY2510_v2/MANIFEST.yaml
status: raw (immutable symlink; source is itself immutable on shared storage)
known_issues:
  - 201 of 422 samples correspond to our phenotyped strains (R. mucilaginosa only); see analysis/gwas/results/strain_reconciliation/
  - scaffold_21 is a collapsed-repeat/aneuploid-like artifact scaffold (mean depth 1411x); excluded downstream, not in this raw copy
access_restrictions: shared lab storage, not redistributable
tags: [rhodotorula, mucilaginosa, gwas, vcf, gemma, snp, indel, genotypes, NRRL-Y2510, versioned]
```

Versioned GATK hard-filtered genotype VCF pair (SNP + INDEL) for 422 R. mucilaginosa strains sequenced
at the Population Genomics Core (UC Riverside). This version (v2, dated 2025-02-04) incorporates updated
sample metadata and filtering criteria compared to v1. Files are symlinked (not copied) from the shared
lab storage at `/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/vcf/`
to preserve storage and maintain a single canonical source. Checksums and resolved paths recorded in
`MANIFEST.yaml`. A subset (201 of 422 samples) overlaps with our phenotyped strains; strain reconciliation
conducted in `analysis/gwas/results/strain_reconciliation/`.

### pigmentation-pathway-hmms
```yaml
name: pigmentation-pathway-hmms
type: other
source: external "pigment bioprotectant discovery" genome/metagenome mining project (hmmbuild, HMMER 3.4); hand-copied from ~/pigment_HMMs/
date_acquired: 2026-08-26
format: HMMER3/f text profile HMM (amino-acid), 133 files
rows: n/a (28 gene profiles + 1 combined library + 1 unverified duplicate + 103 exploratory files in specific/)
columns: n/a
size: ~18 MB
raw_path: data/raw/pigmentation-pathway-hmms/
metadata_path: data/metadata/pigmentation-pathway-hmms/
status: raw (immutable)
known_issues:
  - mysA_test.hmm same size as mysA.hmm but not byte-identical; origin undocumented upstream
  - 9/28 profiles flagged "medium confidence" (broad-domain architecture) by source report
  - No hmmpress index files included; not validated against any Rhodotorula sequence data
access_restrictions: none
tags: [pigmentation, carotenoid, melanin, mysporine, scytonemin, hmm, hmmer, gene-annotation, external, reference]
```

HMMER3 profile HMM library spanning 7 pigment biosynthesis pathways (cyanobacterial MAA/scytonemin/
carotenoid; fungal DHN-melanin/DOPA-melanin/pyomelanin/carotenoid) from an external mining pipeline,
imported as a reference resource for annotating pigmentation genes in Rhodotorula genomes/proteomes.
The fungal and cyanobacterial carotenoid profiles are most directly relevant to this project's
carotenogenic pigmentation phenotypes. Includes a `specific/` subdirectory of exploratory
discriminative-HMM work for t3hnr/t4hnr aimed at reducing false positives against the broad bacterial
SDR superfamily. See `PIGMENTATION_PATHWAY_HMMS.md` for the full gene/pathway table and caveats.

### pigmentation-pathway-reference
```yaml
name: pigmentation-pathway-reference
type: other
source: external "pigment bioprotectant discovery" project report + diagram; hand-copied from ~/report_pigmentation_genes.md and ~/pigmentation_pathways_diagram.png
date_acquired: 2026-08-26
format: Markdown (17 KB) + PNG (1.4 MB)
rows: n/a
columns: n/a
size: 1.4 MB
raw_path: data/raw/pigmentation-pathway-reference/
metadata_path: data/metadata/pigmentation-pathway-reference/
status: raw (immutable)
known_issues:
  - Report references source-project file paths (/mnt/shared-workspace/..., /mnt/results/...) not accessible from this repo
  - No Rhodotorula genomes analyzed in the source report; background/reference material only
access_restrictions: none
tags: [pigmentation, carotenoid, melanin, uv-protection, report, diagram, external, reference]
```

Methods report and pathway-diagram reference material documenting how the `pigmentation-pathway-hmms`
HMM library was built and validated (285 genomes + 3,275 metagenome MAGs, GNPS2 compound linking),
imported as background reading for pigmentation pathway biology. Not primary data generated in this
project — see `PIGMENTATION_PATHWAY_REFERENCE.md` for scope and caveats.

### copper-heavy-metal-screen-v0.15.1
```yaml
name: copper-heavy-metal-screen-v0.15.1
type: other
source: shared ArrayedHeavyMetalScreen project's 0.15.1_Analysis pipeline (R/tidyverse); same imaging runs as copper-colony-measurements, reprocessed
date_acquired: 2026-08-27
format: CSV (16 files copied) + 3 symlinked large sources (large_source/)
rows: 298 strains (copper_auc_mean_by_strain.csv); 1081 strain x configuration rows (copper_auc_all_strains_by_configuration.csv); 306 strains cross-metal (comparative_tolerance_*)
columns: see data/metadata/copper-heavy-metal-screen-v0.15.1/schema.yaml
size: ~925 KB (copied files) + ~1.4 GB symlinked, not duplicated on disk
raw_path: data/raw/copper-heavy-metal-screen-v0.15.1/
metadata_path: data/metadata/copper-heavy-metal-screen-v0.15.1/
status: removed
last_present_in: data-v1-pre-metal-replace

known_issues:
  - AUC here is area-under-(growth-rate-vs-Cu-concentration), NOT area-under-(area-vs-time) like this project's existing AUC_0/AUC_10/AUC_20/AUC_30 GWAS traits -- different phenotype axis, not a drop-in replacement
  - copper_radial_growth_rates_all_concentrations.csv / _30mM.csv cover only strain 254 (prototype), not all strains
  - source dir's copper_measurements_combined.csv vs copper_measurements_combined(1).csv are not byte-identical; only the non-(1) file symlinked
  - strain identity uses shared-project Strain ID/Strain code, not this project's strain_code -- needs reconciliation before merging with GWAS panel
  - 8/306 strains in the cross-metal comparison lack a copper AUC value
access_restrictions: shared lab storage, not redistributable
tags: [copper, heavy-metal, rhodotorula, dose-response, auc, growth-rate, tolerance, comparative, phenotype, update]
```

Updated copper phenotyping outputs (dose-response AUC, linear radial-growth-rate model, cross-metal
comparative tolerance, terminal per-colony QC) from the shared lab's `0.15.1_Analysis` pipeline, covering
the same underlying imaging runs (d000353-d000357) as this project's existing `copper-colony-measurements`
dataset but adding a new growth-model layer (linear radial extension, superseding earlier logistic/power-law
fits per source decisions D-015/D-016/D-018) and a strain-level dose-response AUC not previously computed in
this project. See `COPPER_HEAVY_METAL_SCREEN_V0_15_1.md` for the full file layout, caveats, and suggested
next step (comparing `mean_auc_rate` against the existing `cu_dose_slope`/`AUC_ratio_10`/`resilience_30`
GWAS traits before deciding whether to add or replace).
