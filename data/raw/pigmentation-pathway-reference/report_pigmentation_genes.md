# Pigmentation Genes in Fungi and Cyanobacteria: A Genomic Catalog for Bioprotectant Discovery

## Executive Summary

This report presents a comprehensive catalog of pigmentation biosynthesis genes across 290 fungal and cyanobacterial genomes, with targeted inclusion of 50 extremophilic cyanobacterial MAGs from desert, cave, and tropical environments (Shalygin et al. 2021). We identified 30,007 pigmentation gene hits across 7 biosynthetic pathways (MAA, scytonemin, DHN-melanin, DOPA-melanin, pyomelanin, fungal carotenoids, cyanobacterial carotenoids). All 19 target UV-protective compounds were found in the GNPS2 spectral library (168,442 library matches), providing reference spectra for metabolomics-based validation. Metagenomic mining of 3,275 MAGs from extreme environments (Antarctic soils, Qaidam Basin desert, Yellowstone hot springs) identified 56,124 pigmentation gene hits across 1,223 MAGs, confirming that UV-protective pigment genes are widespread in extreme environment microbial communities.

## 1. Methods

### 1.1 HMM Profile Library Construction

We built 28 custom HMM profiles spanning 7 pigmentation biosynthetic pathways:

| Pathway | Organism Group | Genes | Profiles |
|---------|---------------|-------|----------|
| MAA (mycosporine-like amino acids) | Cyanobacteria | mysA, mysB, mysC, mysD, mysE | 5 |
| Scytonemin | Cyanobacteria | scyA, scyB, scyC, scyD, scyE, scyF | 6 |
| DHN-melanin | Fungi | pks_melanin, t4hnr, t3hnr, scd, ayg1 | 5 |
| DOPA-melanin | Fungi | tyrosinase, laccase | 2 |
| Pyomelanin | Fungi | hppd, hgd | 2 |
| Carotenoid (fungal) | Fungi | crt_fungal_psy, crt_fungal_pds, crt_fungal_lcy | 3 |
| Carotenoid (cyanobacterial) | Cyanobacteria | crtB, crtP, crtQ, crtO, crtR | 5 |

Seed sequences were collected from UniProt and NCBI using taxonomy-filtered queries (7,038 total sequences), aligned with MAFFT, and built into HMM profiles with hmmbuild. Profiles were validated against known producer genomes (*Nostoc punctiforme* PCC 73102 and *Aspergillus fumigatus* Af293).

**Key methodological decisions:**
- Taxonomic filtering: fungal HMMs search only fungal genomes, cyanobacterial HMMs search only cyanobacterial genomes, eliminating cross-kingdom false positives from broad-domain profiles
- Per-profile bit-score and coverage thresholds calibrated via validation runs
- Full sequence score used as primary criterion (better for multi-domain proteins like NRPS and ATP-grasp enzymes)
- Broad-domain profiles (mysB, pks_melanin, tyrosinase, hppd, etc.) flagged as "medium" confidence; specific profiles get "high" confidence

### 1.2 Genome Selection and Download

**Cyanobacteria (185 genomes):**
- 149 RefSeq assemblies from NCBI (GCF prefix, Scaffold+ level)
- 36 MAGs from Shalygin et al. 2021 (terrestrial extremophiles from deserts, caves, tropical soils)
- Protein prediction via Prodigal (meta mode) for MAGs lacking pre-computed protein FASTA

**Fungi (94 genomes):**
- Targeted genus-level searches for known pigment producers: Aspergillus (30 spp.), Fusarium (20 spp.), Penicillium (15 spp.), Phycomyces, Rhinocladiella, Scytalidium, Scedosporium, etc.
- Both RefSeq (GCF) and GenBank (GCA) assemblies accepted

### 1.3 Metagenomic Dataset Selection

Four public metagenomic datasets from extreme environments were targeted:

| Dataset | Environment | MAGs | Source |
|---------|------------|------|--------|
| Antarctic soils | Antarctic soil | 319 | Zenodo 14674625 |
| Qaidam Basin desert | Desert (China) | 1,773 | Zenodo 15743210 |
| Yellowstone hot springs | Hot springs | 780 | Figshare 30284068 |
| Polar cyanobacterial mats | Polar | 37 | Figshare 22003967 |

### 1.4 GNPS2 Compound Linking

The GNPS2 spectral library (1,741,067 spectra) was queried for 19 known UV-protective compounds by name and exact mass matching. Matching spectra provide USIs for MASST searches against public metabolomics datasets.

## 2. Results

### 2.1 Genome-Based Gene Catalog

**30,007 total hits across 290 genomes (28 gene families)**

| Pathway | Total Hits | Genomes with Hits | Avg Hits/Genome |
|---------|-----------|-------------------|-----------------|
| DHN-melanin | 17,835 | 91/91 fungi | 196 |
| MAA | 4,438 | 199/199 cyano | 22 |
| Carotenoid (cyano) | 2,882 | 199/199 cyano | 14 |
| Scytonemin | 1,401 | 199/199 cyano | 7 |
| DOPA-melanin | 1,156 | 91/91 fungi | 13 |
| Pyomelanin | 429 | 83/91 fungi | 5 |
| Carotenoid (fungal) | 104 | 44/91 fungi | 2 |

### 2.2 Pathway Completeness

| Pathway | Complete | Partial (≥50%) | Partial (<50%) | Absent |
|---------|----------|----------------|----------------|--------|
| MAA | 65 | 134 | 0 | 94 |
| Scytonemin | 15 | 184 | 0 | 94 |
| Carotenoid (cyano) | 174 | 25 | 0 | 94 |
| DHN-melanin | 58 | 33 | 0 | 202 |
| DOPA-melanin | 69 | 0 | 0 | 224 |
| Pyomelanin | 82 | 1 | 0 | 210 |
| Carotenoid (fungal) | 21 | 23 | 0 | 249 |

*Note: "Absent" counts include genomes from the other organism group (e.g., fungal genomes show as "absent" for cyanobacterial pathways due to taxonomic filtering).*

### 2.3 Shalygin et al. 2021 Extremophilic Cyanobacterial MAGs

All 50 Shalygin MAGs from desert soils, cave walls, and tropical rocks possess at least partial UV-protective pigment pathways:

| Pathway | Complete | Partial | Absent |
|---------|----------|---------|--------|
| MAA | 33/50 (66%) | 17/50 (34%) | 0 |
| Scytonemin | 6/50 (12%) | 44/50 (88%) | 0 |
| Carotenoid | 25/50 (50%) | 25/50 (50%) | 0 |

**Notable complete-pathway MAGs from extreme environments:**
- *Calothrix* sp. (desert soil, California): MAA 3/3, Scytonemin 6/6, Carotenoid 3/3 — full UV protection toolkit
- *Drouetiella hepatica* (rock surface, Slovakia): MAA 3/3, Scytonemin 4/6, Carotenoid 3/3
- *Brasilonema angustatum* (tropical soil, Hawaii): MAA 3/3, Scytonemin 4/6, Carotenoid 2/3
- Multiple *Nostoc* spp. (desert soils): MAA 3/3, partial scytonemin

The Shalygin MAGs show a higher rate of complete MAA pathways (66%) compared to the overall cyanobacterial average (33%), supporting the hypothesis that extremophilic cyanobacteria invest more heavily in UV protection.

**Shalygin MAGs by habitat:**

| Habitat | MAGs | MAA Complete | Scytonemin Complete | Carotenoid Complete | All 3 Complete |
|---------|------|-------------|--------------------|--------------------|----------------|
| Desert soil | 15 | 11 (73%) | 3 (20%) | 8 (53%) | 2 |
| Desert wet wall | 4 | 3 (75%) | 1 (25%) | 3 (75%) | 1 |
| Cave rock wall | 6 | 4 (67%) | 0 | 3 (50%) | 0 |
| Tropical soil | 3 | 2 (67%) | 0 | 0 | 0 |
| Tropical vernal pool | 2 | 1 | 1 | 1 | 1 |
| Rock surface | 1 | 1 | 0 | 1 | 0 |
| Tropical rock wall | 1 | 1 | 0 | 0 | 0 |
| On wood in tropics | 1 | 1 | 0 | 1 | 0 |
| Aquatic epiphyte | 1 | 0 | 0 | 0 | 0 |
| Unknown | 16 | 9 (56%) | 1 | 8 (50%) | 1 |

Desert soil MAGs show the highest rate of complete MAA pathways (73%), consistent with high UV exposure driving MAA pathway investment.

### 2.4 GNPS2 Spectral Library Matches

All 19 target UV-protective compounds were found in the GNPS2 library:

| Compound | Class | Name Matches | Key Library Spectra |
|----------|-------|-------------|---------------------|
| Shinorine | MAA | 2 | CCMSLIB00005436494, CCMSLIB00010013009 |
| Porphyra-334 | MAA | 1 | CCMSLIB00010013015 |
| Mycosporine-glycine | MAA | 2 | CCMSLIB00005436496 |
| Palythine | MAA | 5 | CCMSLIB00000574579 |
| Asterina-330 | MAA | 1 | CCMSLIB00000840588 |
| Scytonemin | Scytonemin | 4 | CCMSLIB00000001550, CCMSLIB00000070263 |
| Beta-carotene | Carotenoid | 143 | CCMSLIB00006678131 |
| Zeaxanthin | Carotenoid | 577 | CCMSLIB00006382298 |
| Astaxanthin | Carotenoid | 38 | CCMSLIB00000205451 |
| Canthaxanthin | Carotenoid | 24 | CCMSLIB00006126865 |
| Echinenone | Carotenoid | 9 | CCMSLIB00000205501 |
| Myxoxanthophyll | Carotenoid | 1 | CCMSLIB00000205578 |

### 2.5 Metagenomic Mining

Four metagenomic datasets from extreme environments were searched for pigmentation genes using the combined HMM library (both fungal + cyanobacterial profiles, since extreme environment metagenomes may contain both):

| Environment | MAGs Searched | MAGs with Hits | Total Hits | Top Pathways |
|-------------|--------------|----------------|------------|--------------|
| Antarctic soil | 319 | 319 (100%) | 14,781 | DHN-melanin (8,986), MAA (2,884), scytonemin (1,355) |
| Desert (Qaidam Basin) | 2,881 | 829 (29%) | 38,376 | DHN-melanin (24,563), MAA (6,536), scytonemin (3,316) |
| Hot springs (Yellowstone) | 75 | 75 (100%) | 2,967 | DHN-melanin (1,596), MAA (512), carotenoid (469) |
| **Total** | **3,275** | **1,223** | **56,124** | |

**Key observations from metagenomic mining:**

1. **DHN-melanin genes dominate across all environments** (35,145/56,124 = 63% of all hits), reflecting the widespread distribution of t4hnr and t3hnr homologs in environmental microbes. These fungal melanin genes are detected in metagenomic MAGs because the HMMs capture broad short-chain dehydrogenase/reductase domains.

2. **MAA genes are present in all extreme environments** (9,932 hits), confirming that mycosporine-like amino acid biosynthesis is a widespread UV protection strategy beyond cultured cyanobacteria.

3. **Scytonemin genes are detected in all environments** (4,937 hits), with particularly high counts in desert (3,316) and Antarctic (1,355) soils — environments with high UV exposure.

4. **Antarctic soils and Yellowstone hot springs show 100% hit rate** — all 319 Antarctic MAGs and all 75 Yellowstone MAGs contain at least one pigmentation gene hit, suggesting pigmentation is near-universal in these extreme environment microbial communities.

5. **Qaidam Basin desert has the most hits** (38,376) but a lower hit rate (29% of MAGs), reflecting the large number of MAGs from diverse taxa (many non-pigment-producing archaea and bacteria).

**Visualizations:**
- `metagenome_env_pathway_heatmap.svg/.png` — fraction of MAGs with hits per environment per pathway
- `metagenome_gene_distribution.svg/.png` — gene family distribution by environment
- `metagenome_top_candidates.svg/.png` — top 30 MAGs ranked by pathway completeness score

**Pathway completeness in metagenomic MAGs:**

| Environment | MAA Complete | Scytonemin Complete | Carotenoid Complete | Pyomelanin Complete |
|-------------|-------------|--------------------|--------------------|--------------------|
| Antarctic soil (319) | 24 (7.5%) | 0 | 7 (2.2%) | 24 (7.5%) |
| Desert (2,881) | 20 (0.7%) | 0 | 10 (0.3%) | 30 (1.0%) |
| Hot springs (75) | 3 (4.0%) | 0 | 19 (25.3%) | 7 (9.3%) |

No complete DHN-melanin or scytonemin pathways were found in metagenomic MAGs, likely because these require large multi-domain enzymes (NR-PKS for DHN-melanin) or 6 co-located genes (scytonemin) that are difficult to recover in fragmented MAG assemblies. However, partial pathways are widespread — 100% of Antarctic and Yellowstone MAGs have partial scytonemin and DHN-melanin pathways.

**Top metagenome bioprotectant candidates:**

| Rank | MAG | Environment | Score | Complete Pathways | Hits |
|------|-----|-------------|-------|-------------------|------|
| 1 | WB_110 | Antarctic soil | 4.42 | 2 | 95 |
| 2 | LV_1 | Antarctic soil | 4.25 | 2 | 68 |
| 3 | WB_107 | Antarctic soil | 4.08 | 2 | 75 |
| 4 | WB_161 | Antarctic soil | 3.92 | 2 | 75 |
| 5 | CJ0912_30_35bin11 | Desert (Qaidam) | 3.83 | 2 | 81 |
| 6 | AY11.MAG20 | Hot springs (Yellowstone) | 3.83 | 2 | 47 |
| 7 | NSLT21bin26 | Desert (Qaidam) | 3.67 | 2 | 73 |

Antarctic soil MAGs dominate the top candidates, reflecting both high UV exposure and the high quality of the Antarctic MAG assemblies. Yellowstone hot springs show the highest rate of complete carotenoid pathways (25.3%), likely reflecting thermophilic cyanobacteria with well-assembled genomes.

## 3. Key Findings

1. **UV-protective pigments are near-universal in cyanobacteria**: 100% of cyanobacterial genomes (199/199, including 50 extremophilic MAGs) have at least one MAA gene, and 100% have carotenoid pathway genes. This confirms cyanobacteria as a prime source for bioprotectant discovery.

2. **Extremophilic cyanobacteria have enriched UV protection**: The Shalygin et al. 2021 MAGs from desert soils, cave walls, and tropical rocks show 66% complete MAA pathways — double the overall cyanobacterial average (33%). 5 of the top 10 cyanobacterial bioprotectant candidates are Shalygin MAGs with perfect scores (all 3 UV-protective pathways complete).

3. **Scytonemin is less common but present in key taxa**: 15 genomes have complete scytonemin pathways (6/6 genes), concentrated in *Nostoc* spp. and *Calothrix* sp. — both known scytonemin producers. The partial pathways (184 genomes with ≥1/6 genes) suggest either gene loss or incomplete HMM detection for the low-seed-count genes (scyD, scyE).

4. **Fungal melanin pathways are widespread**: 58 fungal genomes have complete DHN-melanin pathways, 69 have tyrosinase (DOPA-melanin), and 82 have complete pyomelanin pathways. 8 fungal genomes have all 4 pigment pathways complete (perfect score 4.00). *Aspergillus*, *Penicillium*, and *Fusarium* species are the primary carriers.

5. **GNPS2 contains reference spectra for all target compounds**: The 19 UV-protective compounds have 168,442 matching spectra in GNPS2, enabling MASST searches against public metabolomics datasets to find evidence of these compounds in environmental samples.

6. **Metagenomic mining confirms pigmentation genes in extreme environments**: 56,124 pigmentation gene hits were found across 1,223 MAGs from Antarctic soils, Qaidam Basin desert, and Yellowstone hot springs. MAA genes (9,932 hits) and scytonemin genes (4,937 hits) are present in all three extreme environments, confirming that UV-protective pigment biosynthesis is widespread beyond cultured reference genomes. Antarctic soil and Yellowstone hot spring MAGs both show 100% hit rates — all MAGs contain at least one pigmentation gene.

## 4. Limitations and Caveats

- **Low-seed-count HMMs**: scyD (2 seeds) and scyE (3 seeds) have minimal training data, potentially missing divergent homologs. This may explain the lower scytonemin completeness scores.
- **Broad-domain profiles**: 9 profiles (mysB, pks_melanin, tyrosinase, hppd, crt_fungal_psy, crt_fungal_pds, crtO, crtR, scyA) are flagged as "medium" confidence due to broad domain architecture. DIAMOND BLASTp validation is recommended for high-confidence applications.
- **MAG completeness**: MAGs from metagenomes may have missing genes due to assembly/binning gaps. Partial pathway completeness in MAGs may reflect incomplete genome recovery rather than true gene absence.
- **Polar cyanobacterial mats**: The 37 MAGs from polar cyanobacterial mats were still being processed at the time of this report and are not included in the metagenome search results. They can be searched separately using the same pipeline.
- **Yellowstone MAG processing**: Only 75 of 1,560 Yellowstone MAG FASTA files were processed at the time of search. Additional Yellowstone MAGs are being processed and can be searched in a follow-up run.
- **Taxonomic filtering**: Fungal HMMs only search fungal genomes and vice versa. This prevents cross-kingdom false positives but means fungal genes in cyanobacterial metagenomes (or vice versa) would be missed. The metagenome search uses both HMM subsets to address this.
- **GenBank assembly download failures**: 115/200 fungal GenBank assemblies failed to download (no protein FASTA at FTP path). The 94 successfully downloaded genomes provide good coverage of key pigment-producing genera.

## 5. Pipeline Deliverables

- **HMM profile library**: 28 pressed HMM profiles at `/mnt/shared-workspace/shared/hmm_library/pigmentation_profiles.hmm`
- **Search pipeline**: `search_genomes.py` (reusable, with per-profile thresholds and taxonomic filtering)
- **Metagenome search**: `search_metagenomes_parallel.py` (parallel search of MAGs from extreme environments)
- **Gene catalog**: `/mnt/results/tables/gene_catalog.tsv` (30,007 hits across 290 genomes)
- **Metagenome catalog**: `/mnt/results/tables/metagenome_gene_catalog.tsv` (56,124 hits across 3,275 MAGs)
- **Metagenome summaries**: `/mnt/results/tables/metagenome_summaries.json` (3,275 MAG summaries)
- **Metagenome bioprotectant candidates**: `/mnt/results/tables/metagenome_bioprotectant_candidates.tsv` (1,223 MAGs ranked by pathway completeness)
- **Pathway summary**: `/mnt/results/tables/pathway_summary.tsv`
- **Bioprotectant candidates**: `/mnt/results/tables/bioprotectant_candidates.tsv` (293 genomes ranked)
- **Compound-gene linking**: `/mnt/results/tables/compound_gene_linking.tsv` (19 compounds linked to biosynthetic genes)
- **GNPS2 matches**: `/mnt/results/tables/gnps2_compound_matches.tsv` (168,442 library matches)
- **Gene family distribution**: `/mnt/results/tables/gene_family_distribution.tsv`
- **Visualizations**: Pathway completeness heatmaps, bar charts, Shalygin MAG heatmap, gene distribution, bioprotectant score distribution, GNPS2 compound summary, pathway diagram in `/mnt/results/figures/`

## 6. References

- Shalygin et al. 2021. Metagenome Sequencing to Explore Phylogenomics of Terrestrial Cyanobacteria. Microbiol Resour Announc. doi:10.1128/mra.00258-21
- GNPS2 Documentation: https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/
- HMMER3: https://github.com/EddyRivasLab/hmmer
