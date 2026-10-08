-- Non-redundant strain table, derived from heavy_metal_measurement.
-- The strain columns repeat on every colony row in the Parquet files. Checked 2026-10-07:
-- no strain_id has conflicting values for any of these columns, across all 5 metals.
-- Some rows carry NULL where other rows of the same strain have a value, so max()
-- (which ignores NULL) picks the one non-NULL value per strain.
-- strain_id is VARCHAR: 321 numeric IDs plus 13 'Control-N' labels that have no metadata.
-- Rows with NULL strain_id (10,921) belong to no strain and are excluded here.
CREATE OR REPLACE VIEW strain_info AS
SELECT
    strain_id,
    try_cast(strain_id AS INTEGER) IS NULL        AS is_control,
    max(SAMPLE_NAME)                              AS sample_name,
    max(STRAIN)                                   AS strain,
    max(SPECIES)                                  AS species,
    max(MS2_SAMPLE_Cell)                          AS ms2_sample_cell,
    max(MS2_SAMPLE_Supernatant)                   AS ms2_sample_supernatant,
    max(Location)                                 AS location,
    list_sort(list(DISTINCT Metal))               AS metals_tested,
    count(*)                                      AS n_colony_observations
FROM heavy_metal_measurement
WHERE strain_id IS NOT NULL
GROUP BY strain_id;
