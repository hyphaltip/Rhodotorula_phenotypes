# Session Log 2026-08-17-009 — Tier D/E figure + publication methods paragraph

## What was done
- **Tier D/E figure** — wrote `scripts/make_tierde_figure.py` following the repo figure
  conventions (matplotlib Agg, argparse `--out-dir`, PNG+PDF at dpi=150, 2-panel):
  `results/gwas/figures/tierde_gene_finemap.{png,pdf}`.
  - Panel A: gene map of the scaffold_10 anchor region (chroma lead 384,905 & AUC_10/DBP3
    lead 396,172, ~11 kb apart): -log10 p of the FDR-significant SNPs per trait, the lead
    SNP annotated, and a schematic gene track built from the genes overlapping window SNPs
    (`tierD_fdr_snps_annotated.csv`); window 364,905–416,172 → 34 sig SNPs / 10 genes.
  - Panel B: Tier-E 14 anchors' 95% credible-set size (log x) vs lead ABF posterior
    probability, coloured by the rare-driven flag (af<0.05).
- **Publication methods paragraph** — added GWAS_REPORT.md §10 covering the Tier A→G chain
  for direct insertion into a methods section: population (201 haploid isolates, QC 404,706
  variants, 20,769 LD-pruned → centered GRM, clone-mean phenotypes) → Tier A kinship-only
  GEMMA `-lmm 4 -k` (all-201 + culled-173; BH-FDR q<0.05 because GEMMA lacks permutation)
  → Tier B burden/SKAT/min-p on pixy high-dxy windows (MC-verified top-50) → Tier C BSLMM
  (100k MCMC, 20% burn-in) → LOCO sensitivity → Tier D gene annotation (6,799 genes, 250 kb
  clump) → Tier E z-space Wakefield ABF credible sets (NCP SD=0.2, logsumexp, p<1e-3
  candidates, 90/95/99%) → Tier G nearest-proxy replication of chr13:13_30149.
- GWAS_REPORT.md §9.4 cross-linked the new figure + reproduce command.
- **TODO_REGISTRY**: updated the open dataviz item to include the new set of figures; added
  `rare-ef-imputation.md` (high) and `functional-followup.md` (medium) with writeups in
  `todo/`. PROGRESS.md §14 extended with the figure + methods items and cross-links.

## Key numbers
- Locus window (scaffold_10, both anchors): 34 FDR-sig SNPs, 10 overlapping genes.
- Tier-E panel: 14 anchors at 95% CS; log-y for the 1–207 SNP range keeps the n=1
  resolved singleton (resil_lead_13) readable next to the n=131 BSLMM clusters.

## Caveat
- The figure was generated but I cannot visually inspect it in this session (model without
  image input) — pipeline numbers are sane, but the layout (gene-track label overlap in
  Panel A) still warrants a human/dataviz pass before publication (TODO_REGISTRY open item).

## Living-repo updates
- ANALYSIS_MANIFEST: `gwas-tierdeg` entry += `make_tierde_figure.py`, figure outputs,
  §10 methods in GWAS_REPORT, extended reproduce command.
- PROGRESS.md §14 (figure + methods items), `.claude/last-session.md` refreshed.
