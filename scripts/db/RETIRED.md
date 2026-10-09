# Retired database scripts (2026-10-08)
Removed because they built or queried the old `colony_measurement` table and its views, which `25_import_heavy_metal.py` replaced (D-34, D-35, D-58).

| Removed | Role |
|---|---|
| `00_init_schema.sql` | dimension tables of the old design |
| `05_generate_metadata.py` | metadata for the old Copper dataset |
| `10_import_experiment.py` | experiment and strain CSV import |
| `20_import_measurements.py` | `colony_measurement` import |
| `30_create_views.sql` | `v_phenotype` and related views |
| `lib/imagename.py` | image-name parser used only by 05 and 20 |
| `query_examples/` | queried `v_phenotype` |

Recover with `git checkout pre-old-db-importers-removal -- scripts/db/<file>`.
Still present: `40_data_dictionary.py` (kept, regenerates SCHEMA.md), `25_import_heavy_metal.py`, `35_create_strain_view.sql`, `check_heavy_metal_db.py`, `lib/db.py`.
Not changed: `DATABASE_DESIGN.md` and `SCHEMA.md` describe the old design and are stale. Analysis folders that read `colony_measurement` are left in place by request.
