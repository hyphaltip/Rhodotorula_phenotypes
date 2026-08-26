# Summary Statistics: pigmentation-pathway-hmms

<!-- Generated: 2026-08-26 -->
<!-- Script: manual (ls / du / grep at ingestion time) -->

## Overview

| Property | Value |
|----------|-------|
| Files | 133 total (29 top-level .hmm files + 1 top-level .hmm combined library counted above + 103 files in `specific/`) |
| Individual gene profiles | 28 |
| Combined/pressed library file | 1 (`pigmentation_profiles.hmm`, 28 NAME records — a concatenation of the 28 individual profiles, not additional content) |
| Unverified duplicate | 1 (`mysA_test.hmm`) |
| `specific/` (t3hnr/t4hnr discriminative HMM benchmarking) | 103 files |
| File size (apparent) | ~18 MB total |
| Date range | N/A (static reference library; HMMs built 2026-08-15 per file headers) |
| Format | HMMER3/f text profile HMM |

## Column summaries

Not tabular; see `schema.yaml` for per-HMM logical fields. Pathway breakdown of the 28 individual
profiles (from source report, `## 1.1 HMM Profile Library Construction`):

| Pathway | Organism group | Genes | Profiles |
|---------|---------------|-------|----------|
| MAA | Cyanobacteria | mysA, mysB, mysC, mysD, mysE | 5 |
| Scytonemin | Cyanobacteria | scyA, scyB, scyC, scyD, scyE, scyF | 6 |
| Carotenoid (cyanobacterial) | Cyanobacteria | crtB, crtP, crtQ, crtO, crtR | 5 |
| DHN-melanin | Fungi | pks_melanin, t4hnr, t3hnr, scd, ayg1 | 5 |
| DOPA-melanin | Fungi | tyrosinase, laccase | 2 |
| Pyomelanin | Fungi | hppd, hgd | 2 |
| Carotenoid (fungal) | Fungi | crt_fungal_psy, crt_fungal_pds, crt_fungal_lcy | 3 |

## Missing data summary

No pressed index files (`.h3f/.h3i/.h3m/.h3p`) present for any profile — run `hmmpress` before
`hmmscan` if a pressed database is required. No missingness within the text HMM files themselves.

## Quality flags

- `mysA_test.hmm` not byte-identical to `mysA.hmm` despite identical file size — flagged, not
  deleted (raw data immutability).
- 9 of 28 profiles flagged "medium confidence" in source report (see `provenance.md`).
- `specific/` subdirectory is exploratory/intermediate work product (multiple candidate HMM
  variants per gene), not a clean single-profile-per-gene library — consult
  `specific/RECOMMENDED_CONFIG.md` before reusing.

## Notes

Verified with `grep -c '^NAME' pigmentation_profiles.hmm` = 28, consistent with the source
report's "28 pressed HMM profiles" claim (the combined file is a concatenation, not additional
profiles).
