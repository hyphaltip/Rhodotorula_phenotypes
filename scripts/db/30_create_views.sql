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

-- v_phenotype, v_strain_experiment_summary and v_growth_timeseries were built on the
-- legacy colony_measurement table, dropped 2026-10-07 (D-35). Use heavy_metal_measurement
-- (25_import_heavy_metal.py) and strain_info (35_create_strain_view.sql) instead.
