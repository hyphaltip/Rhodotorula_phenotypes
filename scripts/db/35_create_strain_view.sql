-- Non-redundant strain table, derived from heavy_metal_measurement.
-- The strain columns repeat on every colony row in the Parquet files. Checked 2026-10-07:
-- no strain_id has conflicting values for any of these columns, across all 5 metals.
-- Some rows carry NULL where other rows of the same strain have a value, so max()
-- (which ignores NULL) picks the one non-NULL value per strain.
-- strain_id is VARCHAR: 321 numeric IDs plus 13 'Control-N' labels that have no metadata.
-- species: strain_species_override (loaded by 25_import_heavy_metal.py from
-- data/metadata/heavy-metal-array-intermediate/strain_species_overrides.tsv) replaces the
-- source value when present; species_source says which. The source text 'Species Not Found'
-- is a placeholder, not a NULL.
-- Rows with NULL strain_id (10,921) belong to no strain and are excluded here.
CREATE OR REPLACE VIEW strain_info AS
SELECT
    m.strain_id,
    try_cast(m.strain_id AS INTEGER) IS NULL        AS is_control,
    max(m.SAMPLE_NAME)                              AS sample_name,
    max(m.STRAIN)                                   AS strain,
    coalesce(any_value(o.species), max(m.SPECIES))  AS species,
    CASE WHEN any_value(o.species) IS NOT NULL THEN any_value(o.source) ELSE 'source data' END AS species_source,
    max(m.MS2_SAMPLE_Cell)                          AS ms2_sample_cell,
    max(m.MS2_SAMPLE_Supernatant)                   AS ms2_sample_supernatant,
    max(m.Location)                                 AS location,
    list_sort(list(DISTINCT m.Metal))               AS metals_tested,
    count(*)                                      AS n_colony_observations
FROM heavy_metal_measurement m
LEFT JOIN strain_species_override o ON o.strain_id = m.strain_id
WHERE m.strain_id IS NOT NULL
GROUP BY m.strain_id;
