# Pigmentation Pathway Reference

Reference report and pathway diagram for microbial UV-protective pigment biosynthesis (fungal +
cyanobacterial), imported as background/reference material for pigmentation-pathway annotation work
in this project. Documents the methodology and results behind the HMM library in
`data/raw/pigmentation-pathway-hmms/`.

## Source

Copied from the user's home directory (`~/report_pigmentation_genes.md`,
`~/pigmentation_pathways_diagram.png`) on 2026-08-26. Authored by a separate "pigment bioprotectant
discovery" project (external to this repository) — not generated here.

## Contents

- `report_pigmentation_genes.md` — "Pigmentation Genes in Fungi and Cyanobacteria: A Genomic Catalog
  for Bioprotectant Discovery". Covers: HMM library construction (7 pathways, 28 profiles), genome
  selection (185 cyanobacterial + 94 fungal genomes, including 50 Shalygin et al. 2021 extremophilic
  cyanobacterial MAGs), metagenomic mining of 3,275 MAGs from Antarctic soil / Qaidam Basin desert /
  Yellowstone hot springs / polar cyanobacterial mats, GNPS2 spectral-library compound linking, and
  pathway-completeness results/limitations.
- `pigmentation_pathways_diagram.png` — schematic diagram of the 7 pathways (cyanobacterial MAA,
  scytonemin, carotenoid; fungal DHN-melanin, DOPA-melanin, pyomelanin, carotenoid), showing gene
  order and end-product UV-protective compounds per pathway. Useful as a quick visual reference when
  interpreting HMM hits from `data/raw/pigmentation-pathway-hmms/` by pathway/gene order.

## Format

Markdown report (17 KB) + PNG diagram (1.4 MB, ~1536x1024px based on layout).

## Known issues

- This is a **reference document from an external project**, not primary data generated in this
  repository. Genome/metagenome counts, gene catalogs, and file paths referenced in the report
  (e.g. `/mnt/shared-workspace/...`, `/mnt/results/...`) point to the source project's environment
  and are **not accessible from this repository** — only the report text and diagram were copied.
- None of the genomes analyzed in the report are Rhodotorula; relevance to this project is
  restricted to (a) the fungal carotenoid pathway description/diagram as background, and (b) the
  HMM library itself (`pigmentation-pathway-hmms` dataset) which could in principle be run against
  Rhodotorula proteomes.
- Report text uses first person plural ("we built...") referring to the source project's authors,
  not this project's team.

## Related

- `data/raw/pigmentation-pathway-hmms/` — the HMM profile library this report documents.
