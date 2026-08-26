# Provenance: pigmentation-pathway-hmms

## Source

**Type**: derived (external project deliverable, hand-copied)

**Origin**:
- Derived from: an external "pigment bioprotectant discovery" genome/metagenome mining project
  (not part of this repository). HMMs were built with `hmmbuild` (HMMER 3.4, Aug 2023) from
  UniProt/NCBI seed sequences aligned with MAFFT. Copied into this repo from `~/pigment_HMMs/`
  in the user's home directory.

**Citation / accession**: N/A — unpublished internal pipeline output. See
`data/raw/pigmentation-pathway-reference/report_pigmentation_genes.md` for the full methods
write-up and `## 6. References` for the one external citation used by that project
(Shalygin et al. 2021, doi:10.1128/mra.00258-21 — genome source, not the HMM build itself).

## Acquisition details

**Date acquired**: 2026-08-26

**Obtained by**: Jason Stajich (jasonst@ucr.edu)

**Method**: Copied via `cp -r ~/pigment_HMMs/. data/raw/pigmentation-pathway-hmms/` (manual
copy from home directory into the project's `data/raw/`).

**Checksum**: not recorded per-file (133 files); `crtB.hmm` spot-checked byte-identical to
source via `diff -q` at copy time.

## Access restrictions

**Restriction level**: none

**Details**: None. Files were plain user-readable (mode 600) in the source location; permissions
normalized to group-readable on copy.

## Known issues

- `mysA_test.hmm` is same file size as `mysA.hmm` but not byte-identical; origin/purpose
  undocumented upstream (see dataset README).
- 9 profiles flagged "medium confidence" (broad-domain architecture) in the source report —
  see dataset README and `schema.yaml`.
- No pressed (`hmmpress`) index files included.
- Not validated against any Rhodotorula sequence data; thresholds are from the source project's
  validation genomes (*Nostoc punctiforme*, *Aspergillus fumigatus*).

## Contact

**Primary contact**: Jason Stajich (jasonst@ucr.edu)

**Backup contact**: None.

## Version history

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-08-26 | Initial ingestion from `~/pigment_HMMs/` |
