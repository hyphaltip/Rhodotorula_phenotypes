# Pigmentation Pathway HMMs

HMMER3 profile HMM library for 7 microbial pigment biosynthesis pathways (fungal + cyanobacterial),
built for a separate bioprotectant-discovery genome/metagenome mining project and imported here as a
reference resource for annotating pigmentation genes in Rhodotorula genomes/proteomes.

## Source

Copied from `~/pigment_HMMs/` (user's home directory) on 2026-08-26. Built by an external "pigment
bioprotectant discovery" pipeline (see the accompanying report in
`data/raw/pigmentation-pathway-reference/report_pigmentation_genes.md`) — not generated within this
repository. Seed sequences from UniProt/NCBI (taxonomy-filtered), aligned with MAFFT, built with
`hmmbuild` (HMMER 3.4, Aug 2023).

## Contents

28 individual profile HMMs (one per gene), one combined/pressed library (`pigmentation_profiles.hmm`,
28 `NAME` records — concatenation of the 28 individual profiles), one apparent test/duplicate file
(`mysA_test.hmm` — same byte size as `mysA.hmm` but different content; provenance/purpose unclear,
kept as-is since raw data is immutable), and a `specific/` subdirectory of exploratory
discriminative-HMM benchmarking work for `t3hnr`/`t4hnr` (see below).

| Pathway | Organism group | Genes (profiles) |
|---|---|---|
| MAA (mycosporine-like amino acids) | Cyanobacteria | mysA, mysB, mysC, mysD, mysE |
| Scytonemin | Cyanobacteria | scyA, scyB, scyC, scyD, scyE, scyF |
| Carotenoid (cyanobacterial) | Cyanobacteria | crtB, crtP, crtQ, crtO, crtR |
| DHN-melanin | Fungi | pks_melanin, t4hnr, t3hnr, scd, ayg1 |
| DOPA-melanin | Fungi | tyrosinase, laccase |
| Pyomelanin | Fungi | hppd, hgd |
| Carotenoid (fungal) | Fungi | crt_fungal_psy, crt_fungal_pds, crt_fungal_lcy |

Rhodotorula is a carotenogenic basidiomycete yeast, so the **fungal carotenoid** (crt_fungal_psy/pds/lcy)
and, as a cross-check, **cyanobacterial carotenoid** (crtB/P/Q/O/R) profiles are the most directly
relevant to this project's pigmentation phenotypes; the melanin and MAA/scytonemin profiles are
included for completeness but are not expected to have strong hits in Rhodotorula.

### `specific/` subdirectory (103 files)

Exploratory work building **discriminative** HMMs for `t3hnr`/`t4hnr` to reduce false-positive matches
to the broader bacterial SDR (short-chain dehydrogenase/reductase, PF00106) superfamily — includes
multiple candidate HMM/FASTA variants (`t3hnr_disc_*`), validation TSVs, and `RECOMMENDED_CONFIG.md`
documenting the recommended discriminative-column HMM + E-value threshold per gene, with measured
TP/FP rates against 94 fungal genomes and 100 desert metagenome MAGs. Use `RECOMMENDED_CONFIG.md` as
the entry point if reusing these for a t3hnr/t4hnr search.

## Format

HMMER3/f `[3.4 | Aug 2023]` text profile HMM format, amino-acid alphabet. Directly usable with
`hmmscan`/`hmmsearch` (and `hmmpress` if searching with the pressed binary indexes, which are not
included here and would need to be regenerated with `hmmpress`).

## Known issues

- **Broad-domain / lower-confidence profiles**: per the source report, `mysB`, `pks_melanin`,
  `tyrosinase`, `hppd`, `crt_fungal_psy`, `crt_fungal_pds`, `crtO`, `crtR`, `scyA` are flagged
  "medium confidence" due to broad domain architecture — recommend confirmatory BLASTp/DIAMOND
  validation of hits from these profiles.
- **Low-seed-count profiles**: `scyD` (2 seeds) and `scyE` (3 seeds) have minimal training data and
  may miss divergent homologs.
- **`mysA_test.hmm`** is same-size but not byte-identical to `mysA.hmm`; its purpose/origin is not
  documented upstream. Treat as untrusted/experimental until clarified — prefer `mysA.hmm`.
- **No pressed index files** (`.h3f/.h3i/.h3m/.h3p`) included; run `hmmpress` locally before
  `hmmscan` if needed.
- These HMMs were validated against known producer genomes (*Nostoc punctiforme* PCC 73102,
  *Aspergillus fumigatus* Af293) in the source project — not against any Rhodotorula genome. Treat
  bit-score/E-value thresholds from the source report as a starting point, not calibrated for this
  organism.

## Related

- `data/raw/pigmentation-pathway-reference/` — the source project's report and pathway diagram
  documenting how this HMM library was built and validated, and genome/metagenome mining results.
