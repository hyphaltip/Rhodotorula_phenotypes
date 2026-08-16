---
session_id: 2026-08-15-004
project: rhodotorula-phenotypes
branch: main
started: 2026-08-15T20:25:00-0700
ended: 2026-08-15T21:00:00-0700
duration_minutes: 35
files_changed: ~15
---

## Session Log

### 20:25 — Tree ingest
- Copied the 3 PHYling FastTree protein tree files (support + nosupport treefiles + log) from the shared Rhodotorula_Metabolites results dir into `data/raw/rhodotorula-phyling-protein-tree/` (asset name finalized as `rhodotorula-phyling-protein-tree`).
- Wrote `data/metadata/rhodotorula-phyling-protein-tree/{provenance.md,schema.yaml,summary_stats.md}` + registered it in `data/DATA_MANIFEST.md`.

### 20:35 — Idea 09 phylogeny signal analysis
- New `scripts/idea_09_phylogeny.R`: join tree tips to per-strain traits (slope/intercept_logchroma, l10med_fixed, partial_slope_sd_cu, pace_loglog median per strain), prune to 277 matched tips, compute phylo signal in two scopes (all, within R. mucilaginosa ~200 tips).
- Discovered tree is a **near-comet**: 167/541 edges ≤1e-7 (giant near-duplicate mucilaginosa polytomy), 22 zero-length pendant edges (jittered 1e-9). Blomberg's K is degenerate (~1e-7 all traits; BM power check K=2.25 on this tree), lambda numerically unstable (geiger incoherent) -> switched primary statistic to **Mantel permutation** (Spearman, 999 perms).

### 20:50 — Results
- **Phylogenetic signal IS present** (all-strains): size r=0.43 (p=0.001), within-strain heterogeneity r=0.31 (p=0.001), baseline chroma r=0.19 (p=0.001), Cu-sensitivity slope r=0.10 (p=0.003); pigment pace n.s. Within mucilaginosa only chroma (p=0.002) + size (p=0.027) survive -> structure is mostly between-species. Reconciles the idea-05 "continuum" with real relatedness structure.
- Species monophyly (>=3 tips): dairenensis/diobovata/graminis/kratochvilovae/sphaerocarpa/sp. clade I mono; mucilaginosa/paludigena/taiwanensis/toruloides NOT (consistent with the polytomy).

### 21:00 — Documentation
- Wrote `09-phylogeneticist.md`; updated FINDINGS.md (idea 09 section), 00_index (row 09), ANALYSIS_MANIFEST (script/outputs/findings/tags), `.living/findings/phenotype-color-space.md` (F5 signal + F6 K-degeneracy), learning L-14, decision D-7, session log.

### Next
- Commit plus optional push. Possible follow-up: concatenated-loci species tree check; within-clade Mantel on the non-mucilaginosa singleton clades; reconcile with idea 05 atlas via trait-on-tree figures.

### 21:10 — Idea 10: scale of within-species variation (user request)
- `scripts/idea_10_species_variation.py`: builds per-strain trait table (5 key traits incl. median pace), per-species boxplots (n>=3) + exact SS variance decomposition (between vs within species).
- Result: within-species (strain) variance DOMINATES for all traits — Cu slope 87%, chroma 92%, size 62%, heterogeneity 67%, pace 84% within-species fraction; colony size most species-structured (38% between). Reconciles idea 09 Mantel signal (sits in the smaller between-species component) + idea 05 continuum.
- Updated 00_index (row 10), FINDINGS (idea 10 section + synthesis pt 7), ANALYSIS_MANIFEST (script/outputs/findings/tags), `.living/findings/phenotype-color-space.md` (F7).
