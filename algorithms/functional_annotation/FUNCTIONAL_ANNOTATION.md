# Functional Annotation (pigment-pathway HMM scan + genome-wide KEGG KO)

## Purpose

A reusable, two-layer functional-annotation practice for genes in the *R. mucilaginosa*
NRRL Y-2510 genome (and, via the same scripts, any funannotate-style proteome + GFF3):

1. **Targeted pathway search** — scan a proteome against a curated HMM library for a
   specific biological question (here: pigmentation biosynthesis), with a data-driven
   calibration step since per-profile score thresholds are frequently not available.
2. **Genome-wide baseline** — KEGG Orthology (KO) assignment via KofamScan for every
   gene, so any future gene-list query (a GWAS locus, a differential-expression hit, …)
   can be checked against KEGG pathway membership directly, not just by keyword-matching
   product names (which this project's first attempt at candidate-gene identification —
   see `analysis/gwas/GWAS.md` §15 — showed misses real pathway genes whose product name
   in the existing annotation is generic, e.g. "hypothetical protein" or even actively
   misleading).

Built for the `analysis/gwas/` fine-mapping follow-up (D-17 through D-19), generalized
into a standing convention per user request (2026-08-26) rather than a one-off script.

## Inputs

| Input | Type | Description |
|-------|------|-------------|
| proteome FASTA | FASTA (amino acid) | e.g. `/bigdata/stajichlab/shared/projects/Rhodotorula/MAT_search/db/Rhodotorula_mucilaginosa_NRRL_Y-2510.proteins.fa` (8,542 proteins; IDs match `gene_index.json`) |
| pigment HMM library | HMMER3 profile HMMs | `data/raw/pigmentation-pathway-hmms/` (28 profiles, 7 pathways; imported resource, see `data/DATA_MANIFEST.md`) |
| gene annotation index | JSON | `analysis/ideas/2026-08-15-color-phenotype-space/results/gwas/tierD/scripts/gene_index.json` (product/GO/InterPro/Pfam per gene, from funannotate) |
| KOfam HMM database | HMMER3 profiles + `ko_list` | `/bigdata/stajichlab/shared/lib/gtotree/kofamscan_data/` (already on shared storage; `eukaryote.hal` subset used for speed) |
| gene list to query (optional) | CSV/list of gene IDs or scaffold:pos | Cross-reference target, e.g. a GWAS locus window |

## Outputs

| Output | Type | Description |
|--------|------|-------------|
| `results/pigment_scan_<genome>/*.domtbl` | HMMER domtblout | Raw per-profile hits, one dir per genome scanned |
| `results/calibrated_thresholds.csv` | CSV | Per-profile empirical bit-score-gap threshold (see Parameters) |
| `results/confident_hits.csv` | CSV | Hits above the calibrated threshold, with existing-annotation concordance for the focal genome |
| `results/kofamscan/<genome>.detail.txt` | KofamScan detail | Every profile hit per protein, `*` marks the KofamScan-recommended best call |
| `results/kofamscan/<genome>.mapper.txt` | TSV | `protein<TAB>KO` — the simple form for joining against a gene list |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--clump-r2` / hmmscan `-E` | 1e-5 (pigment scan) | Permissive HMMER inclusion threshold before gap-calibration narrows it further |
| gap-threshold method | largest score gap | See `calibrate_pigment_hmm_thresholds.py` docstring — an unsupervised "twilight zone" heuristic, NOT equivalent to the source project's labeled true/false-positive calibration (94 fungal genomes + 100 metagenome MAGs) which was not available to reproduce here |
| t3hnr / t4hnr | discriminative HMMs, E<=1e-5 / E<=1.0 | Use `data/raw/pigmentation-pathway-hmms/specific/{t3hnr,t4hnr}_disc_*.hmm` per `RECOMMENDED_CONFIG.md`, NOT the originals (55% false-positive rate against the broad bacterial SDR superfamily) |
| KofamScan `-E` | 1e-4 | Standard KofamScan default, matches the lab's existing usage pattern (`Black_yeasts/KEGG_KOFAM_Scan`) |
| KofamScan profile set | `eukaryote.hal` | ~15,700 of the full 27,233 KO profiles; restricts to eukaryote-relevant KOs for speed without losing recall on a fungal proteome |

## Dependencies

- HMMER 3.4 (`module load hmmer/3.4`) — pigment-pathway HMM scan
- KofamScan 1.3.0 (`module load kofamscan`) — genome-wide KO assignment
- pandas — parsing/joining outputs

## Usage

```bash
# 1. Targeted pathway scan (fast, minutes) -- one genome
bash algorithms/functional_annotation/scripts/run_pigment_hmm_scan.sh \
  <proteome.fa> algorithms/functional_annotation/results/pigment_scan_<name>

# 2. Calibrate thresholds across however many genomes you've scanned (more genomes =
#    a better empirical gap; re-run whenever a new genome is added to results/pigment_scan_*)
pixi run python3 algorithms/functional_annotation/scripts/calibrate_pigment_hmm_thresholds.py \
  --results-glob "algorithms/functional_annotation/results/pigment_scan_*" \
  --focal-genome "pigment_scan_<name>" \
  --gene-index <gene_index.json> \
  --out-thresholds algorithms/functional_annotation/results/calibrated_thresholds.csv \
  --out-hits algorithms/functional_annotation/results/confident_hits.csv

# 3. Genome-wide KEGG KO baseline (slow, hours -- submit as a SLURM job)
sbatch --partition=stajichlab --mem=64G --time=12:00:00 --cpus-per-task=8 \
  --wrap="bash algorithms/functional_annotation/scripts/run_kofamscan.sh <proteome.fa> <out_prefix>"
```

## Implementation Notes

- **Why not just trust product-name keyword search**: the first pass at candidate-gene
  ID for the GWAS loci (`analysis/gwas/scripts/finemap_candidate_genes.py`, GWAS.md §15)
  keyword-searched existing product names for carotenoid-pathway terms and found
  nothing — this pipeline exists because that was the wrong tool: the true fungal
  carotenoid synthase/cyclase gene in this genome (`OM429_003333`) is annotated only
  "hypothetical protein" in the existing GFF3, and the adjacent desaturase
  (`OM429_003336`) is actively mis-annotated as a "transcription factor, PHD finger
  motif" — both were only found by direct HMM homology, not by reading annotation text.
- **Calibration honesty**: the empirical bit-score-gap method and Pfam/InterPro
  concordance check are the best calibration achievable from data already on hand, not
  a claim of equivalent rigor to the source project's actual validated thresholds
  (which used real TP/FP labels from 94 fungal genomes + 100 metagenome MAGs). Treat
  `calibrated_thresholds.csv` as a triage/prioritization aid — profiles flagged "medium
  confidence" in `data/raw/pigmentation-pathway-hmms/PIGMENTATION_PATHWAY_HMMS.md`
  (mysB, pks_melanin, tyrosinase, hppd, crt_fungal_psy, crt_fungal_pds, crtO, crtR, scyA)
  still warrant confirmatory BLASTp/DIAMOND review of any hit before treating it as
  established.
- **Extending the background for gap-finding**: the calibration script pools scores
  across every genome under `results/pigment_scan_*` — scanning more related proteomes
  (any funannotate-style `.proteins.fa`) before re-running the calibration step
  strengthens the empirical gap, especially for pathways with few hits in any single
  genome (e.g. `mysB` had only 1 hit total across 11 genomes and could not be
  meaningfully calibrated — flagged `calibratable=False`).
