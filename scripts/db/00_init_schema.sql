-- Dimension tables for the Rhodotorula phenotyping database.
-- See DATABASE_DESIGN.md for the full design rationale.
--
-- colony_measurement (the fact table) is deliberately NOT created here: its
-- trait columns come from whatever the source Parquet files contain, so
-- scripts/db/20_import_measurements.py creates it from the first file's
-- schema and ALTERs it if later files add columns (DATABASE_DESIGN.md §9).
--
-- DuckDB foreign keys are enforced on insert but do not support
-- ON DELETE CASCADE (verified against this DuckDB build) -- the importer
-- deletes colony_measurement rows before deleting their image row on a
-- re-delivered file, rather than relying on the database to cascade it.

CREATE SEQUENCE IF NOT EXISTS experiment_id_seq START 1;

CREATE TABLE IF NOT EXISTS experiment (
    experiment_id     INTEGER PRIMARY KEY DEFAULT nextval('experiment_id_seq'),
    experiment_name   VARCHAR UNIQUE NOT NULL,
    factor_name       VARCHAR NOT NULL,
    factor_unit       VARCHAR,
    plate_info_source VARCHAR,
    strain_info_source VARCHAR,
    notes             VARCHAR
);

CREATE TABLE IF NOT EXISTS imager_run (
    run_number       INTEGER PRIMARY KEY,
    experiment_id    INTEGER NOT NULL REFERENCES experiment(experiment_id),
    library_plate    INTEGER,
    is_control_only  BOOLEAN NOT NULL DEFAULT FALSE
);

CREATE TABLE IF NOT EXISTS strain (
    strain_id     INTEGER PRIMARY KEY,
    strain_code   VARCHAR,
    strain_name   VARCHAR,
    species       VARCHAR,
    origin        VARCHAR,
    environment   VARCHAR
);

CREATE TABLE IF NOT EXISTS well_placement (
    run_number         INTEGER NOT NULL REFERENCES imager_run(run_number),
    configuration      INTEGER NOT NULL,
    well_position       INTEGER NOT NULL,
    strain_id           INTEGER NOT NULL REFERENCES strain(strain_id),
    incubation_temp_c   DOUBLE,
    media               VARCHAR,
    PRIMARY KEY (run_number, configuration, well_position)
);

CREATE TABLE IF NOT EXISTS condition_plate (
    experiment_id    INTEGER NOT NULL REFERENCES experiment(experiment_id),
    plate_number     INTEGER NOT NULL,
    replicate_label  VARCHAR,
    is_control       BOOLEAN NOT NULL DEFAULT FALSE,
    run_number       INTEGER NOT NULL REFERENCES imager_run(run_number),
    configuration    INTEGER,
    PRIMARY KEY (experiment_id, plate_number)
);

CREATE TABLE IF NOT EXISTS condition_plate_factor (
    experiment_id  INTEGER NOT NULL,
    plate_number   INTEGER NOT NULL,
    factor_name    VARCHAR NOT NULL,
    factor_value   DOUBLE,
    factor_unit    VARCHAR,
    PRIMARY KEY (experiment_id, plate_number, factor_name),
    FOREIGN KEY (experiment_id, plate_number) REFERENCES condition_plate(experiment_id, plate_number)
);

CREATE TABLE IF NOT EXISTS image (
    image_name        VARCHAR PRIMARY KEY,
    run_number        INTEGER NOT NULL REFERENCES imager_run(run_number),
    plate_number      INTEGER NOT NULL,
    temperature_token INTEGER,
    imaged_at         TIMESTAMP NOT NULL,
    source_path       VARCHAR,
    file_size_bytes   BIGINT,
    file_mtime        TIMESTAMP
);
