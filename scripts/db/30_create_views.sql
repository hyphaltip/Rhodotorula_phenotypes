-- Analysis-ready views. See DATABASE_DESIGN.md §6.
-- Run after 00_init_schema.sql, 10_import_experiment.py, and
-- 20_import_measurements.py.

-- Per-image time axis, computed with a window function so it's always
-- correct regardless of the order files were imported in.
-- hours_since_plate_start is a DOUBLE (not INTERVAL) on purpose: INTERVAL
-- columns aren't representable in Arrow's fixed interval layout the way
-- Polars/pandas expect it, so R/Python clients reading this view via Arrow
-- would hit a conversion error. epoch-seconds arithmetic keeps it a plain
-- numeric column.
CREATE OR REPLACE VIEW v_image AS
SELECT
    i.*,
    (epoch(i.imaged_at) - MIN(epoch(i.imaged_at)) OVER (PARTITION BY i.run_number, i.plate_number))
        / 3600.0 AS hours_since_plate_start
FROM image i;

-- One row per plate, factors collapsed into a MAP so a plate with more than
-- one manipulated factor (a future temperature x pH experiment) still joins
-- into v_phenotype as exactly one row per plate.
CREATE OR REPLACE VIEW v_condition_factors AS
SELECT
    experiment_id,
    plate_number,
    map(list(factor_name), list(factor_value)) AS factors
FROM condition_plate_factor
GROUP BY experiment_id, plate_number;

-- The main analysis view: one row per colony observation, every trait
-- column, with strain and experimental condition attached. Every join from
-- colony_measurement onward is LEFT so a colony row is never silently
-- dropped by a missing/late-loading dimension -- unmatched dimensions come
-- through as NULL and are easy to filter out downstream
-- (e.g. WHERE strain_id IS NOT NULL).
CREATE OR REPLACE VIEW v_phenotype AS
SELECT
    e.experiment_name,
    cp.replicate_label,
    cp.is_control,
    cf.factors,
    s.strain_id,
    s.strain_code,
    s.strain_name,
    s.species,
    i.image_name,
    i.run_number,
    i.plate_number,
    i.imaged_at,
    i.hours_since_plate_start,
    cm.* EXCLUDE (image_name)
FROM colony_measurement cm
JOIN v_image i                    ON i.image_name = cm.image_name
LEFT JOIN condition_plate cp        ON cp.run_number    = i.run_number
                                   AND cp.plate_number  = i.plate_number
LEFT JOIN experiment e               ON e.experiment_id = cp.experiment_id
LEFT JOIN v_condition_factors cf     ON cf.experiment_id = cp.experiment_id
                                   AND cf.plate_number  = cp.plate_number
LEFT JOIN well_placement wp         ON wp.run_number    = i.run_number
                                   AND wp.configuration = cp.configuration
                                   AND wp.well_position = cm.well_position
LEFT JOIN strain s                   ON s.strain_id = wp.strain_id;

-- Convenience aggregate over a handful of the most-used traits (colony size
-- and intensity). Analysts needing other traits query v_phenotype directly
-- with their own GROUP BY -- this is a starter, not an attempt to
-- pre-aggregate all ~150 traits.
CREATE OR REPLACE VIEW v_strain_experiment_summary AS
SELECT
    experiment_name, factors, replicate_label,
    strain_id, strain_code, strain_name,
    COUNT(*)                             AS n_colonies,
    AVG(Shape_Area)                      AS mean_area,
    STDDEV_SAMP(Shape_Area)              AS sd_area,
    AVG(Intensity_MeanIntensity)         AS mean_intensity,
    STDDEV_SAMP(Intensity_MeanIntensity) AS sd_intensity
FROM v_phenotype
WHERE strain_id IS NOT NULL
GROUP BY ALL;

-- Kinetics view: same grain as v_phenotype, ordered for time-course
-- plotting, restricted to a couple of the most common traits.
CREATE OR REPLACE VIEW v_growth_timeseries AS
SELECT
    experiment_name, factors, replicate_label,
    strain_id, strain_code, image_name, hours_since_plate_start,
    Shape_Area, Intensity_MeanIntensity
FROM v_phenotype
WHERE strain_id IS NOT NULL
ORDER BY strain_id, hours_since_plate_start;
