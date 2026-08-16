# Idea 09 — Phylogeneticist: is strain phenotype structured by relatedness?

Persona: phylogeneticist. Question: **after building a strain phylogeny, do the
per-strain color/growth traits carry phylogenetic signal** — i.e. are closely
related strains more phenotypically similar than expected by chance?

This post-campaign idea re-examines idea 05's "continuum, no species clusters"
through an explicit tree. The tree is the shared PHYling protein tree
(`branch: protein-Rhodotorula-taxa_278.fungi_odb10.fasttree`, data asset
`rhodotorula-phyling-protein-tree`), which covers 278 Rhodotorula strains, 216+
of them R. mucilaginosa.

## Inputs

- Tree: `data/raw/rhodotorula-phyling-protein-tree/protein-Rhodotorula-taxa_278.fungi_odb10.fasttree.support.treefile`
  (FastTree 2.2.0, LG2008, CAT 20, SH-like 1000; final LogLk −1,679,594; total branch length 5.637).
  Tip labels need normalization: `<species>_<strain>.proteins.fa` then `.proteins` stripped.
- Traits (per strain): `slope_logchroma_per_mM` + `intercept_logchroma`
  (`results/idea04_reaction_norms.csv`, 300 strains); `l10med_fixed` +
  `partial_slope_sd_cu` (`results/idea01b_within_strain.csv`, 295 strains);
  `pace_loglog` = median per strain from `results/idea08_pigment_pace.csv`
  (7,031 rows, strain not unique → median-aggregate, 313 strains).
- Metadata: `data/strain_metadata.tsv` (strain_id → species), with `strain_code`
  dedupe (DBVPG_3853 → strains 192/195; TFCN_48D-10 → 93/252; keep row with
  species, highest strain_id).

## Method

`scripts/idea_09_phylogeny.R`

1. Load tree (ape::read.tree), normalize tip labels.
2. Join traits to tips via `tip_label ⇄ strain_id` through strain_metadata.
   277/278 tips matched; only `Rhodotorula_mucilaginosa_DH4148` unmatched.
3. Prune tree to matched tips; compute the two scopes: **all strains** and
   **within R. mucilaginosa** (200+ tips).
4. **Primary statistic — Mantel permutation test** (Spearman, on pairwise
   patristic distance vs absolute trait distance, 999 permutations, seed 7).
   Rank-based and permutation-nulled, robust to tree-shape degeneracies.
   Positive r ⇒ close relatives more phenotypically similar (signal).
5. Secondary, tree-shape-sensitive references: Blomberg's K and Pagel's lambda
   (phytools), **flagged unreliable here** (see caveats).

Why not trust K/lambda: the tree is a near-"comet" — a giant polytomy of
near-duplicate R. mucilaginosa genomes (167/541 edges ≤1e-7),
22 zero-length pendant edges (jittered to 1e-9 for covariance solves), and
patristic-distance minima ~7.6e-9. On such topologies Blomberg's K collapses to
~1e-7 regardless of true signal, and likelihood-based lambda optimization becomes
numerically unstable (geiger white vs lambda lnL contradicted phylosig). A BM
power check reproduced K≈2.25 (should be ~1) on simulated data on this tree → K is
not architecture-informative here.

## Results

`results/idea09_phylo_signal.csv` (Mantel permutation primary; K/lambda as
reference only):

| scope | trait | n | Mantel r | p | Pagel λ |
|-------|-------|---|---------|----|--------|
| all | slope_logchroma_per_mM (Cu sensitivity) | 272 | 0.101 | **0.003** | 0.87 |
| all | intercept_logchroma (baseline chroma) | 272 | 0.190 | **0.001** | 0.84 |
| all | l10med_fixed (colony size) | 266 | 0.426 | **0.001** | ~0.004 (ns) |
| all | partial_slope_sd_cu (within-strain dispersion widening) | 266 | 0.309 | **0.001** | 0.83 |
| all | pace_loglog (pigment pace) | 270 | 0.058 | 0.071 | 0.94 |
| mucilaginosa | slope_logchroma_per_mM | 200 | 0.018 | n.s. | 0.87 |
| mucilaginosa | intercept_logchroma | 200 | 0.134 | **0.002** | 0.84 |
| mucilaginosa | l10med_fixed | 200 | 0.077 | **0.027** | ~0 |
| mucilaginosa | partial_slope_sd_cu | 200 | 0.022 | n.s. | 0.83 |
| mucilaginosa | pace_loglog | 200 | −0.008 | n.s. | 0.94 |

Species monophyly check (is.monophyletic on >=3-tip species): dairenensis,
diobovata, graminis, kratochvilovae, sphaerocarpa, sp. clade I are monophyletic;
**mucilaginosa, paludigena, taiwanensis, toruloides are NOT** — consistent with a
huge near-zero mucilaginosa polytomy (201 tips) and cross-species strain mixing at
screened loci.

## Headline finding

**Strain phenotypes DO carry phylogenetic signal — the "contrast with idea 05" is
strain-level (relatedness), not discrete clusters.**
Across all strains, 4 of 5 color/growth traits show significant positive Mantel r
— close relatives are more phenotypically similar than chance. Ordering of signal
strength: **colony size (r=0.43) > within-strain heterogeneity broadening (r=0.31)
> baseline chroma (r=0.19) > Cu-sensitivity slope (r=0.10, weak)**; pigment pace
is not structured (r=0.06, p=0.07).

Within R. mucilaginosa, signal survives only for **baseline chroma (r=0.13,
p=0.002)** and weakly colony size (p=0.027); Cu sensitivity and heterogeneity
broadening have *no* within-species signal. So the phylogenetic structure is
mostly **between-species**, concentrated in the near-duplicate genome cluster,
and does not resolve a quantitative gradation within the dominant species.

Caution on interpretation: the near-zero-branch mucilaginosa polytomy means most
of the "signal" is a species-level dichotomy (same-species ⇒ near-identical
genomes ⇒ similar chroma/size), not evidence of fine-scale strain-level heritable
structure. This is consistent with idea 05's continuum within species + idea 04's
high ICC on shape, and refines idea 07's environment result (only weak trait↔env
links).

## Robustness across method choice

Mantel p-values (rank, permutation) are likely stable (no normal assumptions);
lambda agrees directionally (0.83–0.94) on the signal-bearing traits but is
numerically unreliable on this tree. K should NOT be reported as "no signal".
Permutation seed 7; 999 perms ⇒ p resolution ~1e-3.

## Caveats

- K degeneracy ⇒ do not interpret K≈1e-7 as absence of signal on this tree.
- geiger fitContinuous lambda showed internally inconsistent logL (white lnL >
  lambda lnL) on this tree — use phytools phylosig lambda LRT as reference only.
- Skinny scope: 20/278 tips (7.5%) unmatched or missing traits; union of trait
  availability trims n to 266–272 of 277.
- Single gene-family topology (fungi_odb10) — a species tree from concatenated
  loci could differ; not yet checked.
