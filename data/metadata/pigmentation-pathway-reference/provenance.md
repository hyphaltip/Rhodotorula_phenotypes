# Provenance: pigmentation-pathway-reference

## Source

**Type**: derived (external project report, hand-copied)

**Origin**:
- Derived from: an external "pigment bioprotectant discovery" genome/metagenome mining
  project (not part of this repository). Report and diagram copied from
  `~/report_pigmentation_genes.md` and `~/pigmentation_pathways_diagram.png` in the user's
  home directory.

**Citation / accession**: N/A for the report itself (unpublished internal deliverable). The
report cites Shalygin et al. 2021 (doi:10.1128/mra.00258-21) as the source of its extremophilic
cyanobacterial MAG genomes.

## Acquisition details

**Date acquired**: 2026-08-26

**Obtained by**: Jason Stajich (jasonst@ucr.edu)

**Method**: Manual `cp` from home directory into `data/raw/pigmentation-pathway-reference/`.

**Checksum (SHA256)**:
- `report_pigmentation_genes.md`: `026e39200fe1b63df34018603df7f47e1bc51c444f329098e18f5b98489f45a0`
- `pigmentation_pathways_diagram.png`: `d697dbc800904c8f67734046aad65af001cb1124a859d7b2f6de42fb73969f53`

## Access restrictions

**Restriction level**: none

**Details**: None.

## Known issues

- File paths referenced within the report (`/mnt/shared-workspace/...`, `/mnt/results/...`)
  point to the source project's own environment and are not accessible from this repository —
  only the report text and diagram were copied, not the underlying gene catalogs, figures, or
  pipeline scripts it references.
- No Rhodotorula genomes were analyzed in the source report; relevance here is as background/
  reference material for pigmentation pathway biology and as documentation for the accompanying
  HMM library (`data/raw/pigmentation-pathway-hmms/`).

## Contact

**Primary contact**: Jason Stajich (jasonst@ucr.edu)

**Backup contact**: None.

## Version history

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-08-26 | Initial ingestion |
