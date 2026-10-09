# Strain curation, hybrid/ploidy handling and removal of old Y-2510 / old-Copper work

Working record of a requirements interview ("grilling") started 2026-10-08. Nothing below has been applied yet. Decisions live in `.living/decisions.md` (IDs given). Update the status column as items are agreed, applied or revised.

## User instructions (2026-10-08, verbatim intent)
| # | Instruction | Status |
|---|---|---|
| 1 | Fix species for the Species Not Found strains with genomes | agreed, D-45; not applied |
| 2 | Update the database; note the reasons for metadata fixes | agreed in principle; table design pending |
| 3 | Queue a rerun of all figure generation with the corrected species | pending; scope in open questions |
| 4 | Jettison the old GWAS and old Copper data; everything uses the new phenotype data in the DB | pending; scope and method in open questions |
| 5 | No hybrid or diploid strains in the GWAS, haploids only | agreed in principle; coding in D-47 |
| 6 | Add a ploidy status column so hybrids are noted | agreed, D-47; not applied |
| 7 | Throw out all old Y-2510 reference-genome work for now; user will give the "new positioning" | pending; meaning of "new positioning" unclear |
| 8 | Do not use the old population labels; use DH4148 SNP clusters or MASH clusters instead | pending; which source and where |
| 9 | Fix species names and ploidy in the DuckDB strain table, and record them in a new generated table in the data folder that is ingested | agreed in principle; columns pending |
| 10 | For a GWAS panel: de-clone and remove hybrid strains | agreed in principle; method pending |

## Decisions agreed so far
| ID | Decision | Revisit if |
|---|---|---|
| D-45 | Use the sourmash call for the 5 Species Not Found strains with genomes (1A-1-5, 2M-1-3, 86C-3 = R. mucilaginosa; 7-6-3 = R. diobovata; 211C-2 = R. aff. babjevae) | lab or ITS check contradicts a call |
| D-46 | 12 popgen 'other species' strains: 5 become R. frigidialcoholis; 7 stay R. mucilaginosa with an aff_mucilaginosa marker | popgen reclassifies aff. mucilaginosa; user wants them as a separate group |
| D-48 | GWAS panel excludes all 12 identity-flagged strains (133 to 121) | popgen issues #2/#3 resolve a strain; panel too small |
| D-49 | De-clone at 5 SNPs, re-pick representatives from unflagged strains: 126 strains. How to change the cutoff: `analysis/popgen_vcf_check/DECLONE_CUTOFF_NOTES.md` | different cutoff wanted; flags resolved; popgen rebuilds groups |
| D-51 | Population labels: new DH4148 SNP clusters (PCA and tree cut) made by us; old labels dropped | popgen/user supplies labels; MASH disagrees; different k |
| D-52 | Gene-level work waits for DH4148 gene positions/annotation supplied by the user | annotation arrives |
| D-53 | Phenotype figures: 'R. mucilaginosa' = pure haploid strains; hybrids and aff. mucilaginosa as separate groups | too few strains; popgen reclassifies aff. |
| D-54 | DH4148: haploid, in the pure group for phenotype figures, excluded from the GWAS | popgen adds a DH4148 sample |
| D-55 | Popgen files snapshotted with checksums; script generates the curation table | popgen publishes a new callset |
| D-56 | Iron and Zinc held out as 'incomplete source' (partial runs); Cr, Cu, Pb only | owner supplies complete tables, or 0.15.1 tables ingested |
| D-57 | Rerun figures now without populations; SNP clusters built later | user prefers one combined rerun |
| D-47 | `ploidy_status` = haploid / diploid_hybrid / diploid_other / unknown, plus `hybrid_subgroup` | popgen reclassifies; unknown strains get genomes |

## Open questions (queue, asked one at a time)
1. ~~Which identity-flagged strains are excluded from the GWAS panel?~~ Answered: all 12 (D-48).
2. ~~De-clone method~~ Answered: re-pick at 5 SNPs, 126 strains (D-49). CG001 is handled by kinship in the GWAS, not collapsed.
3. Old data removal. Method answered: tag, then `git rm`, with a RETIRED.md (D-50). Still to confirm: the exact directory list.
4. ~~Population labels~~ Answered: new DH4148 SNP clusters made by us (D-51). Still open: what the report shows for population analyses until they exist.
5. ~~What does "new positioning" mean?~~ Answered: DH4148 gene positions / annotation, to be supplied by the user (D-52).
6. ~~Rerun scope for hybrids and other species~~ Answered: they stay in the phenotype figures as separate groups (D-53).
7. ~~Reference strain DH4148~~ Answered (D-54).
8. ~~Snapshot or read by path~~ Answered: snapshot with checksums (D-55).
9. Columns of the curation table (species, species_source, species_previous, ploidy_status, hybrid_subgroup, popgen group, identity flags, notes, date).

## Execution plan (drafted 2026-10-08 for confirmation; nothing applied)
Order of work, each step on a branch and its own PR so it can be reviewed and reverted separately:
1. **Snapshot and curation table (D-55, D-45 to D-49, D-54).** Snapshot the popgen files with checksums into `data/raw/popgen-callset-metadata/`; a script generates `data/metadata/strain-curation/strain_curation.csv` (columns approved in the interview); a note file records every species and ploidy fix and its reason.
2. **Database (D-47, D-56).** The DB build ingests the curation table as `strain_curation`; the `strain_info` view takes species and ploidy from it and carries `ploidy_status`, `hybrid_subgroup`, the GWAS flags and a per-metal source status (Fe and Zn = incomplete_source).
3. **Removal of the old work (D-50).** Tag, `git rm` the retired directories, add `RETIRED.md`, update the manifests.
4. **GWAS panel list (D-48, D-49).** Generate the 126-strain panel (haploid R. mucilaginosa, no identity flags, de-cloned at 5 SNPs) from the snapshot; store the group ids at several cutoffs in the curation table. Preparation only: no GWAS is run.
5. **Figure rerun (D-53, D-56, D-57).** Rerun the a* report pipeline on Cr, Cu, Pb with the corrected groups; rewrite the report (Fe and Zn sections removed or marked pending; populations pending); rebuild the PDF.
6. **Later, separate steps:** DH4148 SNP clusters (D-51); gene-level work when the user supplies the annotation (D-52); Fe and Zn when complete tables exist (D-56).

## Notes and inconsistencies found but not decided (kept for the metadata-fix notes)
- Strain ids 165 and 269 share the sample name TFCN_17-332Y-1 (BY115-H5); both have no genome and are Species Not Found. Both are flagged; unresolved.
- 11 'Species Not Found' strains have no genome; their species stays unresolved. TFCN_223D-8, TFCN_17-332D-2 and TFCN_17-333P-8 are different strains from similar genotyped names (popgen note, user-confirmed) and have no genotype.
- Cr, Cu and Pb have 10,921 colony rows with no strain id, all in the last two runs of each metal; unexplained.
- Cu b* jumps between dose 0 and dose 5 while a* and size barely change (plate or media effect?); not explained.
- TFCN_17-333M-1 (strain 304) is a hybrid diploid in popgen; its species stays R. mucilaginosa with ploidy_status diploid_hybrid.
- Popgen open decisions that may change our table: hybrid renaming (their decision 6), StrainDB species fixes, EXF_1695 (not phenotyped).
