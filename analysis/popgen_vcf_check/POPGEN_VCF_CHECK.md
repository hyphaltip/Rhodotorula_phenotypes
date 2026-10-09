# Phenotype strains against the latest R. mucilaginosa DH4148-reference callset (2026-10-08)

Source: `/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref` (handoff `docs/HANDOFF_2026-10-04.md`, summary `docs/SUMMARY_findings_2026-10-06.md`). Inputs were only read. Script: `scripts/compare_strains.py`. Tables: `results/phenotype_vs_popgen_strains.csv` (one row per phenotype strain, with group membership, QC and identity-issue flags), `results/popgen_strains_without_phenotype.csv`, `results/phenotype_mucilaginosa_not_in_callset.csv`.

Matching: upper case, letters and digits only, `TF_CN` to `TFCN`. VCF sample names are `NAME_NAME` and were halved. 215 of the 216 R. mucilaginosa phenotype strains match by name; no fuzzy matching was needed for them.

## How many R. mucilaginosa strains can be used

216 phenotype strains carry the label R. mucilaginosa in `strain_info`.

| Group (popgen) | Strains | Use |
|---|---|---|
| `rmuc_core` (pure haploid R. mucilaginosa, passes QC) | **170** | main set |
| of these, `rmuc_core_declone` (one strain per group of at most 5 SNPs) | **133** | non-redundant set |
| `hybrid_diploids` | 33 | not pure R. mucilaginosa; analyse separately |
| dropped by QC as "other species" (popgen metadata: 5 R. frigidialcoholis, 7 R. aff. mucilaginosa) | 12 | in `rmuc_core_outgroup`, not in `rmuc_core` |
| reference strain DH4148 itself (no sample in the callset) | 1 | none |

- Choosing phenotyped strains as the representatives of the near-identical groups does not raise the non-redundant count: it is 133 either way.
- No EXF strain is phenotyped, so the EXF mixed-culture issues (EXF_1695 and 7 others) do not touch these data.
- Against the old GWAS VCF (NRRL Y-2510 reference, 213 panel strains): 170 are in both. The other 43 old-panel strains are 31 hybrid diploids and the 12 "other species" strains. `rmuc_core` contains no strain that was not already in the old panel, so the new callset adds no phenotyped strains. 3 phenotyped strains are new to the callset side: 2 hybrid diploids and DH4148.
- 94 callset strains have no phenotype row.

## Lingering issues that touch these strains

1. **Species labels disagree for 12 strains.** We label them R. mucilaginosa. Popgen metadata calls 5 R. frigidialcoholis (TFCN_134A-3, TFCN_3M-1-1, TFCN_1A-14, TFCN_17-332M-1, DBVPG_6660) and 7 R. aff. mucilaginosa (TFCN_25-395P-1, TFCN_33A-4, TFCN_25-333Y-10, DBVPG_4380, TFCN_25-334Y-6, DBVPG_4534, DBVPG_8043). Popgen's StrainDB fix (12 strains to R. frigidialcoholis, with R. aff. mucilaginosa staying R. mucilaginosa "for now") is proposed, not applied. TFCN_3M-1-1 is also a possible mixed culture.
2. **16 strains in `rmuc_core` conflict with the older `ExRhodotorula_Phenotypes/strains.csv`** (another species there; genotype and our labels say R. mucilaginosa). They are enriched in near-identical groups. Whether the phenotyped and sequenced cultures are the same is unresolved (popgen issue #2, ITS check pending).
3. **Species Not Found strains (16 in our table).**
   - Genotype exists for 5: TFCN_86C-3, TFCN_1A-1-5 and TFCN_2M-1-3 are in `rmuc_core` (so R. mucilaginosa by genotype). TFCN_7-6-3 and TFCN_211C-2 were dropped for low coverage and are in no group.
   - This **contradicts the species I proposed from the old table** for TFCN_86C-3 (sphaerocarpa) and TFCN_1A-1-5 (pacifica). The genotype should take precedence. TFCN_86C-3 is also identical (0-1 SNPs) to DBVPG_6742 and DBVPG_4304 (Italy), a likely mix-up; the lab check is pending (popgen issue #4).
   - The other 10 have no genome in the callset.
4. **TFCN_17-333M-1** (our strain 304, set to R. mucilaginosa by your confirmation) is a hybrid diploid in popgen (`affmuc-like` subgroup).
5. **Near-identical and clonal structure.** All 170 phenotyped `rmuc_core` strains belong to a clone group at 1e-3 divergence: CG001 holds 77, CG002 39, CG003 17. CG001 is a group, not one clone (median about 470 SNPs apart; the trees cannot resolve it). 46 phenotyped strains sit in 9 near-identical (at most 5 SNP) groups. Kinship must be handled in any association test.
6. **9003 batch.** DBVPG_3538 and DBVPG_4379 were rebuilt from their original library (fine). DBVPG_6094 is unverified. TFCN_363-1-2 has a wrong species in popgen metadata.
7. **Duplicate sample name in our table:** TFCN_17-332Y-1 (strains 165 and 269). Both are absent from the callset.
8. **Pipeline status.** The callset depends on branch `fix/homref-dp-mask` (nf_genotype_population#7), pushed but not merged. Several decisions are still open (see handoff).

## Integration blockers for the GWAS
- **Different reference.** The callset is on the DH4148 assembly (GCA_058775505.1). The old GWAS, the candidate-gene work and the gene IDs (`OM429_*`) are on NRRL Y-2510. Positions and genes need a mapping, or the GWAS needs a rerun on the new VCF with the DH4148 annotation (the VCFs are SnpEff-annotated).
- **Use the right group VCF.** `RmucDH4148.rmuc_core.snps.maf.annotated.vcf.gz` (or `.rmuc_core_declone.`) for haploid R. mucilaginosa; hybrids separately.
- **Old population labels** (`pop_assignment_at_run.csv`, 6 populations, 201 strains) come from the old callset. The new tree could give new labels. They were not compared.
