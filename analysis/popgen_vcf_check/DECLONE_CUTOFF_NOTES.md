# How to revisit the GWAS de-clone cutoff (decision D-49)

Current choice (2026-10-08): **5 SNPs**, one strain kept per group, representatives re-picked from unflagged strains. **126 strains.** This file says how to change that, for example to a stricter cutoff.

## Panel size by cutoff
Starting set: our R. mucilaginosa strains in popgen's `rmuc_core` after the species fix (173), minus the 22 identity-flagged strains (D-48; 21 flagged in `rmuc_core` plus TFCN_86C-3) = **151 strains**. Groups are single-linkage components of strains within N SNPs; one strain kept per group.

| Cutoff (SNPs) | Panel, 151 unflagged strains | Reference: all 173 | Reference: popgen's whole `rmuc_core` (247) |
|---|---|---|---|
| 0 | 149 | 169 | not published |
| 2 | 137 | 149 | 218 |
| **5 (current)** | **126** | 135 | 203 |
| 10 | 117 | 127 | 190 |
| 20 | 104 | 113 | 174 |
| 50 | 93 | 99 | 160 |

- Popgen notes the data show **no gap** between "clonal" and "unrelated" at any of these cutoffs, so the choice is a judgement, not a measurement.
- A stricter cutoff removes near-identical strains but also removes real genetic diversity. Cutoffs 10-50 lose 9 to 33 strains compared with 5.
- **None of these cutoffs collapses CG001** (120 strains, median about 470 SNPs apart) or the other clone groups at 1e-3 divergence. All our `rmuc_core` strains belong to a clone group. The GWAS has to model kinship regardless of the cutoff.

## How to change the cutoff
1. **Inputs (read only):** popgen `results/variant_qc/declone/rmuc_core.pairwise.tsv.gz` (columns `strain_a, strain_b, diffs, compared`; 30,382 pairs) and `declone/rmuc_core.groups_le5.tsv` (shows the representative rule: `missing_frac`, `mosdepth_mean`). Path: `/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref/`.
2. **Candidate set S:** strains with `species = Rhodotorula mucilaginosa`, `ploidy_status = haploid`, in popgen `rmuc_core`, not `aff_mucilaginosa`, and with no `gwas_exclude_reason` in the curation table (D-47, D-48).
3. **Groups:** connected components of S where `diffs <= N` (code below).
4. **Representative:** best genotype quality (lowest `missing_frac`, then highest `mosdepth_mean`, as in popgen's file). If a group has no unflagged member it contributes no strain.
5. **Record it:** write the panel and the cutoff used to the curation table (`gwas_panel`, `declone_cutoff`) and to the GWAS run notes, then rebuild the DB.
6. **Rerun** every GWAS step that reads the panel. Panels at different cutoffs are nested only approximately, because representatives are re-picked, so do not compare results between cutoffs by deleting strains from a larger panel.

```python
# groups at cutoff N among strain keys `keys` (norm = upper case, letters and digits only, TF_CN -> TFCN)
par = {k: k for k in keys}
def find(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for a, b in zip(sub.a, sub.b):          # sub = pairwise rows with both strains in `keys` and diffs <= N
    par[find(a)] = find(b)
groups = {}
for k in keys: groups.setdefault(find(k), []).append(k)
```

## Make the cutoff easy to change later (proposal, pending the table-design question)
Store group ids at several cutoffs in the curation table, for example `declone_group_le2`, `_le5`, `_le10`, `_le20`, `_le50`, plus a matching `declone_rep_le*` flag. Changing the cutoff is then a filter on one column, not a recomputation. If the unflagged set changes (a flag is resolved), recompute the columns from the pairwise file.

## Revisit triggers
- A reviewer or the GWAS results (inflation, kinship structure) argue for a stricter or looser cutoff.
- An identity flag is resolved (popgen issues #2, #3, #4): the strain returns and groups must be re-picked.
- Popgen rebuilds its groups or pairwise table.
