# Summary Statistics: rhodotorula-phyling-protein-tree

Computed on 2026-08-15 from the copied tree files (see provenance.md).

## Tree statement

| Property | Value |
|----------|-------|
| Tool | FastTree 2.2.0 (double precision, OpenMP 12 threads) |
| Alignment | Protein concatenation (BUSCO fungi_odb10, PHYling protein_tree) |
| Model | Le-Gascuel 2008, CAT approximation, 20 rate categories |
| Support | SH-like 1000 (on internal nodes) |
| Search | Normal +NNI +SPR (2 rounds range 10) + ML-NNI opt-each=1 |
| Total tips | 278 |
| Internal nodes | 265 |
| Unique sequences | 267 / 278 |
| Bad splits | 1 / 264 (worst delta-LogLk 0.265) |
| Final LogLk | -1679594.151 |
| Total branch length | 5.637 (substitutions/site) |
| Tips ending `.proteins` | 278 (all) |
| Non-Rhodotorula tips (outgroups) | 2 (Cystobasidium, Pseudomicrostroma) |

## Files

| File | Size (B) | md5 |
|------|----------|-----|
| `*.fasttree.support.treefile` | 21,129 | 704de78b2689978279870256b8b416dc |
| `*.fasttree.nosupport.treefile` | 19,809 | 9b32a896c8fcef6001d361dcd4acabf2 |
| `*.fasttree.support.log` | 8,697 | - |

## Join to phenotype metadata (computed in idea 09)

| Property | Value |
|----------|-------|
| Tips matching a strain code in `strain_metadata.tsv` | 277 / 278 |
| Unmatched tip | `DH4148` |
| Trait sources | idea04 reaction norms, idea02 heterogeneity, idea01 arrest, idea01b within-strain, idea08 pigment pace, idea05 atlas |

> Note: trait-availability counts (how many of the 277 matched tips carry each trait) are
> reported in `analysis/ideas/2026-08-15-color-phenotype-space/results/idea09_*.csv`.
