# Popgen callset metadata snapshot

Read-only copy of metadata files from the popgen project. Do not edit. Checksums are in `SHA256SUMS.txt` (verify: `sha256sum -c SHA256SUMS.txt`).

- Source: `/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref/` (files under `popgen/`) and `.../Species_ID_db/` (files under `species_id_db/`).
- Popgen repo HEAD at copy time: `2d5ae40` (2026-10-06). Copied 2026-10-08.
- Callset: 316 strains, GATK joint genotyping on DH4148 (GCA_058775505.1).
- `popgen/variant_qc/rmuc_core.missing_by_strain.tsv` is not a popgen file. `scripts/curation/extract_missing_fraction.sh` made it from the `rmuc_core` VCF (SLURM job 29670578). Check: for the 56 strains in `groups_le5.tsv` the values agree with popgen's `missing_frac` to within 5e-5 (rounding).
- Caveats: popgen's StrainDB fixes are not applied here. The 9003-batch species calls come from one library per strain and are unverified.
- Refresh: copy the files again, rerun `sha256sum`, rerun `scripts/curation/build_strain_curation.py`, and review the diff of `data/metadata/strain-curation/`.
