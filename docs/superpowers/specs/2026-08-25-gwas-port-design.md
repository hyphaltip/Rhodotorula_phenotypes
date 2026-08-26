# GWAS Port to `analysis/gwas/` — Design

**Date**: 2026-08-25
**Status**: Approved by user, revised after Fable quant-genetics review, pending implementation plan

## Context

A complete GWAS pipeline for *Rhodotorula mucilaginosa* color/copper-response
phenotypes was already built and run to completion (2026-08-16/17) inside
`analysis/ideas/2026-08-15-color-phenotype-space/`, using:

- Genotypes: `RmucY2510_v2.All.{SNP,INDEL}.combined_selected.vcf.gz`
  (422 samples; symlinked into this project at
  `/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rmuc_popgen_NRRLY2510/vcf/`,
  resolving to `/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/vcf/`)
- Population assignment: `.../Population_Genomics/.../Rmuc_PopAssigned.csv` (6 pops)
- Coverage: `.../Population_Genomics/.../coverage/mosdepth/<strain>.{5000,10000,50000}bp.mosdepth.{summary,*.dist}.txt`
  (precomputed from CRAMs in `.../Population_Genomics/.../aln/*.cram`, 424 files)
- Phenotype metadata: `data/metadata/Copper.Strain_info.csv`

Recorded decisions/learnings from that run (D-9 through D-13, L-16, L-17,
L-21/L-22) already establish: kinship-only GEMMA LMM as the valid model (PC/pop
covariates degenerate the GRM — 22 near-clone strains make it singular), the
GEMMA fam-column-6 phenotype gotcha, BH-FDR over permutation for significance,
and the full Tier A→G analysis chain (single-SNP LMM → SKAT/burden → BSLMM →
LOCO sensitivity → gene annotation → fine-mapping → prior-locus replication).

Two things the user flagged as *not* rigorously done in that run:
1. Strain-name reconciliation between phenotype metadata and VCF sample IDs
   was a direct match (`our200.txt`), not an audited fuzzy-match/typo pass.
2. Ploidy validation was only an aggregate check ("0 het across all 201
   combined"), not a per-strain flag.

## Decision

Port the existing methodology into a new top-level `analysis/gwas/` folder
(sibling to `analysis/growth_rates/`, `analysis/explore_plate_position/`,
following this repo's per-analysis convention: UPPER_SNAKE_CASE.md doc,
`scripts/`, `results/`, `run.sh`). Reuse already-computed heavy outputs
(pixy, BSLMM, LOCO) rather than re-running multi-hour SLURM jobs, and add the
two missing audits as new, first-class steps that gate the strain list before
anything downstream consumes it.

## Components

### 1. Versioned data ingestion — `data/raw/genotypes/RmucY2510_v2/`

- Symlinks (not copies) to the two VCFs + `.tbi` + `.QC.tsv` sidecars, at
  their resolved realpath. Rationale: 662MB (SNP) + 2.1GB (INDEL) already
  live immutably on shared storage; copying duplicates without benefit.
- `MANIFEST.yaml` in that directory recording: original path, resolved
  realpath, file size, sha256, version tag (`v2`), ingestion date.
- Future updates land in a sibling `RmucY2510_v3/` etc. — never overwritten
  in place, so any analysis's `data/raw/genotypes/RmucY2510_vN/` pin is
  reproducible after the shared project moves on.
- `data/DATA_MANIFEST.md` gets a new entry pointing at this directory.
- CRAM files (`.../aln/*.cram`, 424 files) are **referenced by absolute path
  at analysis time**, not ingested — they're only consulted if a strain in
  our matched set lacks a precomputed mosdepth summary and one must be
  generated fresh.

### 2. Strain-name reconciliation — `analysis/gwas/scripts/reconcile_strains.py`

- Input: `data/metadata/Copper.Strain_info.csv` (phenotype side) x VCF
  sample list (`bcftools query -l`, 422 IDs).
- Matching tiers, most to least confident:
  1. Exact string match
  2. Normalized match (lowercase, collapse `_`/`-`/whitespace, strip
     leading zeros)
  3. Fuzzy/suffix match — generalizes the existing
     `idea_09_phylogeny.R` "longest strain_code that is a suffix of the
     tip label" trick into a reusable function; also try Levenshtein
     ratio above a threshold as a fallback candidate generator
  4. Unmatched
- Every phenotype strain gets exactly one row in
  `results/strain_reconciliation/strain_match_table.csv` with its tier,
  matched VCF ID (if any), and match score. Fuzzy (tier 3) and unmatched
  (tier 4) rows are called out in
  `results/strain_reconciliation/NEEDS_REVIEW.md` for manual adjudication —
  **the script does not auto-accept fuzzy matches or auto-drop unmatched
  strains**; a human-reviewed `strain_match_table.reviewed.csv` (copy with a
  `decision` column filled in) is the actual input to step 4.
- This replaces the ad hoc `our200.txt`/`our201.txt` strain lists used
  previously.
- **Uniqueness assertion** (added after quant-genetics review): before
  writing `strain_match_table.csv`, assert 1:1 matching — no VCF sample ID
  is claimed by more than one phenotype strain. A collision (two phenotype
  strains fuzzy-matching the same VCF ID) is a hard error surfaced in
  `NEEDS_REVIEW.md`, not silently resolved by taking the higher-scoring
  match.

### 3. Per-strain ploidy/heterozygosity validation — `analysis/gwas/scripts/check_ploidy.py`

- Per-sample het rate on the SNP VCF restricted to the reconciled strain
  set (`bcftools stats -s -`, or equivalent), for every strain individually
  (not aggregated). Het rate defined as `het_sites / callable_sites`
  (callable = non-missing GT calls for that sample), with a minimum
  callable-site floor (e.g. strains with <100k callable sites get
  `flag = insufficient_data` rather than a het-rate verdict — depth-starved
  samples produce noisy het estimates).
- Cross-reference each strain's mean depth from its mosdepth summary file
  (genome-wide mean from the `.mosdepth.summary.txt` total row, not a
  windowed mode/histogram — avoids GC/copy-number-driven window noise);
  if missing, run mosdepth fresh from the strain's CRAM (`.../aln/<id>.cram`).
- **Thresholds are data-driven, not fixed absolutes** (per quant-genetics
  review): compute the panel's own het-rate and depth distributions first.
  - Het outlier: het rate > panel median + 3×MAD(het rate).
  - Depth-ratio band: strain mean depth / panel median depth is treated as
    "elevated" when it falls in **[1.5×, 2.5×]** (consistent with a doubled
    genome; ratios further out — e.g. >3×) are flagged separately as
    likely contamination/mixed-culture rather than diploidy).
- Flag rule: **het outlier AND depth ratio in the 1.5–2.5× band** →
  `ploidy_flag = likely_diploid_or_mixed`. Either signal alone → `watch`.
  Neither → `haploid_ok`. Depth ratio > 2.5× (with or without het signal)
  → `watch_contamination` (distinct from the diploidy call).
- Output: `results/ploidy_check/ploidy_flags.csv` (one row per strain: het
  rate, callable sites, mean depth, depth ratio to panel median, flag,
  the panel median/MAD values used so the thresholds are reproducible)
  plus a short summary in the analysis doc. Flagged strains are surfaced to
  the user for a keep/exclude decision — never auto-excluded.

### 4. GWAS pipeline — ported, reused

- Copy Tier A–G scripts from
  `analysis/ideas/2026-08-15-color-phenotype-space/scripts/` into
  `analysis/gwas/scripts/`, parameterized on the new
  `data/raw/genotypes/RmucY2510_v2/` path and the reviewed strain list
  (step 2) minus any excluded ploidy flags (step 3), instead of the old
  hardcoded `our200.txt`.
- **"Unchanged" is redefined** (per quant-genetics review) as identity of
  the full triple: (1) strain list, (2) each strain's population-group
  assignment (`Rmuc_PopAssigned.csv`), and (3) the near-clone/culled-173
  partition (the prior run scored both an all-201 and a culled-173 set —
  a strain flipping between those groups changes the GRM even if the
  201-name list is untouched). `analysis/gwas/scripts/diff_strain_state.py`
  hashes/diffs all three against the prior run's recorded state and this
  drives the branch below — not a bare strain-name-list comparison.
  - **Unchanged** on all three → copy all existing results (pixy
    π/θ/dxy/Fst, kinship, GEMMA Tier A/B/C, LOCO, fine-mapping, gene
    annotation, replication) into `analysis/gwas/results/`, re-homed and
    re-documented; no SLURM re-run.
  - **Changed** (on any of the three, even a single strain) → mandatory
    steps, in order:
    1. Rebuild the GRM on the new set and re-run the conditioning
       diagnostic (eigenvalue spectrum / condition number) that
       characterized the original singularity (D-9's 22-near-clone
       finding) — confirms whether kinship-only LMM is still well-posed
       before trusting any p-value from it. Record the eigenvalue summary
       in the doc even when it comes back clean.
    2. Re-run kinship + GEMMA Tier A/B/C (~90 min total, cheap) on the new
       set.
    3. Reuse pixy/BSLMM/LOCO from the prior run only if population
       assignments (for pixy) and the near-clone partition (for
       BSLMM/LOCO) are unaffected by the specific strains that changed;
       otherwise flag them explicitly as **stale, computed on the prior
       (N=201) strain set**, pending a rebuild — never silently presented
       as current.
- Carries forward the existing modeling decisions without re-litigating them
  (kinship-only LMM per D-9, BH-FDR per D-11, LOCO as sensitivity check per
  D-10/D-13) — cite them, don't redo the exploration that produced them.

### 5. Documentation & manifests

- `analysis/gwas/GWAS.md`: ported/adapted from `GWAS_REPORT.md`, paths
  updated, plus new sections for strain reconciliation and ploidy check.
- `analysis/gwas/run.sh`: reproduces the port (steps 1–3, then either the
  copy-forward or the re-run branch of step 4) end-to-end.
- `analysis/ANALYSIS_MANIFEST.md`: new `gwas` entry.
- `data/DATA_MANIFEST.md`: new `RmucY2510_v2` genotype entry.
- `.living/decisions.md`: one new entry recording the port and why (reuse
  computed outputs, add the two audits) — references D-9/D-10/D-11/D-13
  rather than duplicating them.

## Out of scope

- Re-litigating the PC-covariate-vs-kinship-only modeling decision (D-9) —
  settled, carried forward.
- Ingesting CRAM files into `data/raw/` — referenced by path only.
- Building GWAS support for species other than *R. mucilaginosa* — this
  panel is R. mucilaginosa-only per the existing 201/278-strain match.
