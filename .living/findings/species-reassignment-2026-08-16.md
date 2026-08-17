# Species re-assignment note (2026-08-16)

Whole-genome fastANI (screen + targeted) in the metabolomics repo
(Rhodotorula_pheno_MS) confirmed 4 strains are mis-specified. See:
`Rhodotorula_pheno_MS/analysis/integrated_analysis/phase_siderophore/ANI_check/SPECIES_REASSIGNMENT.md`

## Strains re-assigned (all records + this repo's DuckDB strain table updated)

| Strain ID | Strain      | Old species (.fixed) | New species |
|-----------|-------------|----------------------|-------------|
| 28        | TFCN_1A-1-3 | R. pacifica          | R. mucilaginosa |
| 80        | TFCN_1B-1-2 | R. pacifica          | R. mucilaginosa |
| 157       | TFCN_1A-1-2 | R. paludigena        | R. mucilaginosa |
| 105       | TFCN_152C-6 | R. toruloides        | R. taiwanensis |

ANI to new species = 99.98-99.999%; to old species = 77.6-78%.

## Unresolved
- Strain ID 288 `TFCN_1A_1_5` (labeled R. pacifica): no genome available for
  ANI; left as-is pending a genome. The other two R. pacifica strains are
  actually R. mucilaginosa, so this is a strong candidate but unverified.

## Impact on this repo's analyses

Any analysis that groups strains by species is affected. Re-check (and
re-run where species grouping feeds the result):

- `analysis/growth_rates/` results tables (species_sensitivity,
  strain_sensitivity, strain_x_cu_rates) — species labels upstream in this
  repo's DuckDB (already fixed) and `data/metadata/Copper.Strain_info.csv`
  (already fixed).
- `analysis/ideas/2026-08-15-color-phenotype-space/` — species-level color/
  growth summaries.
- `analysis/control_late_timepoint_phenotype/` results.
- The phyling protein tree itself (Rhodotorula_Rodeo) was built from the OLD
  labels; tips are keyed by strain not species, so tip identity is fine, but
  any species-colored render should be rebuilt.