#!/usr/bin/env python3
"""Import one experiment's condition and strain metadata into the database.

Registers/updates one row in `experiment`, loads the plate/condition CSV into
`condition_plate` + `condition_plate_factor`, and loads the strain CSV into
`strain` + `well_placement` + `imager_run`. See DATABASE_DESIGN.md §5.

Default plate-numbering convention (Copper today): plate numbers are laid out
in contiguous blocks of --plates-per-run, each block split into
--configs-per-run configurations x (plates-per-run / configs-per-run) factor
levels. This is the same block arithmetic documented as a fragile,
undocumented convention in README.md quirk #4 -- it is applied here once, at
import time, and validated, never recomputed at query time. A future
experiment whose plate numbering does not follow this convention should pass
--plate-run-map instead of relying on the formula.

Usage (Copper):
    python3 10_import_experiment.py \\
        --experiment-name Copper \\
        --factor-name "Copper concentration" \\
        --factor-unit mM \\
        --plate-info-csv data/metadata/Copper.Plate_info.csv \\
        --strain-info-csv data/metadata/Copper.Strain_info.csv
"""

import argparse
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.db import get_connection  # noqa: E402


def resolve_run_and_configuration(plate_number: int, base_run_number: int,
                                   plates_per_run: int, configs_per_run: int) -> tuple[int, int]:
    levels_per_config = plates_per_run // configs_per_run
    block_index = (plate_number - 1) // plates_per_run
    offset_in_block = (plate_number - 1) % plates_per_run
    run_number = base_run_number + block_index
    configuration = offset_in_block // levels_per_config + 1
    return run_number, configuration


def load_plate_run_map(path: str) -> dict[int, tuple[int, int]]:
    df = pl.read_csv(path)
    return {
        int(row["plate_number"]): (int(row["run_number"]), int(row["configuration"]))
        for row in df.iter_rows(named=True)
    }


def sync_strain(con, strain_df: pl.DataFrame) -> None:
    """Upsert the global `strain` table from the strain CSV (CSV authoritative).

    Reruns are idempotent: on a strain_id conflict the metadata columns are
    refreshed from the CSV, so name/species corrections propagate on re-import.
    """
    strain_rows = (
        strain_df.select(
            pl.col("Strain ID").alias("strain_id"),
            pl.col("Strain").alias("strain_code"),
            pl.col("Strain Name").alias("strain_name"),
            pl.col("Species").alias("species"),
            pl.col("Origin").alias("origin"),
            pl.col("Environment").alias("environment"),
        )
        .unique(subset=["strain_id"])
    )
    con.execute("CREATE OR REPLACE TEMP TABLE _strain_stage AS SELECT * FROM strain_rows")
    con.execute(
        """
        INSERT INTO strain SELECT * FROM _strain_stage
        ON CONFLICT (strain_id) DO UPDATE SET
            strain_code = excluded.strain_code,
            strain_name = excluded.strain_name,
            species = excluded.species,
            origin = excluded.origin,
            environment = excluded.environment
        """
    )
    print(f"strain: {strain_rows.height} rows staged")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=None, help="Path to the DuckDB file (default: db/rhodotorula_phenotypes.duckdb)")
    ap.add_argument("--experiment-name", default=None)
    ap.add_argument("--factor-name", default=None, help="e.g. 'Copper concentration', 'Temperature', 'pH'")
    ap.add_argument("--factor-unit", default=None, help="e.g. 'mM', 'C', 'pH units'")
    ap.add_argument("--plate-info-csv", default=None)
    ap.add_argument("--strain-info-csv", required=True)
    ap.add_argument("--plate-number-column", default="Batch Number")
    ap.add_argument("--factor-value-column", default="Concentration (mM)")
    ap.add_argument("--replicate-column", default="Replicate")
    ap.add_argument("--control-value", default="Control",
                     help="Value in --factor-value-column/--replicate-column that marks a control plate")
    ap.add_argument("--plates-per-run", type=int, default=28,
                     help="Block size for the default run/configuration arithmetic (Copper: 28)")
    ap.add_argument("--configs-per-run", type=int, default=4,
                     help="Configurations per run block (Copper: 4)")
    ap.add_argument("--base-run-number", type=int, default=None,
                     help="Defaults to min(Run Number) in the strain CSV")
    ap.add_argument("--plate-run-map", default=None,
                     help="CSV with columns plate_number,run_number,configuration -- "
                          "use instead of the block-arithmetic default")
    ap.add_argument("--notes", default=None)
    ap.add_argument("--strain-only", action="store_true",
                     help="Sync the global `strain` table from the strain CSV and exit. "
                          "Use after editing Strain_info.csv (e.g. a species-name fix) without "
                          "touching imager_run/well_placement, which full re-runs cannot "
                          "rewrite once measurements reference the runs.")
    args = ap.parse_args()

    con = get_connection(args.db)
    con.execute(open(Path(__file__).resolve().parent / "00_init_schema.sql").read())

    strain_df = pl.read_csv(args.strain_info_csv, infer_schema_length=None)
    strain_df = strain_df.rename({c: c.strip() for c in strain_df.columns})
    n_raw = strain_df.height
    strain_df = strain_df.filter(pl.col("Strain ID").is_not_null())
    n_blank = n_raw - strain_df.height
    if n_blank:
        print(f"Dropped {n_blank} fully-blank padding rows from {args.strain_info_csv} "
              f"(Strain ID null) -- {strain_df.height} real strain rows remain")

    if args.strain_only:
        sync_strain(con, strain_df)
        con.close()
        print("Done (strain-only).")
        return

    missing = [a for a, v in (("--experiment-name", args.experiment_name),
                              ("--factor-name", args.factor_name),
                              ("--plate-info-csv", args.plate_info_csv)) if v is None]
    if missing:
        con.close()
        raise SystemExit(f"error: {', '.join(missing)} is required unless --strain-only is set")

    plate_df = pl.read_csv(args.plate_info_csv, infer_schema_length=None)

    base_run_number = args.base_run_number
    if base_run_number is None:
        base_run_number = int(strain_df["Run Number"].drop_nulls().min())

    known_run_numbers = set(int(x) for x in strain_df["Run Number"].drop_nulls().unique().to_list())

    # --- Register the experiment -------------------------------------------------
    con.execute(
        """
        INSERT INTO experiment (experiment_name, factor_name, factor_unit,
                                 plate_info_source, strain_info_source, notes)
        VALUES (?, ?, ?, ?, ?, ?)
        ON CONFLICT (experiment_name) DO UPDATE SET
            factor_name = excluded.factor_name,
            factor_unit = excluded.factor_unit,
            plate_info_source = excluded.plate_info_source,
            strain_info_source = excluded.strain_info_source,
            notes = excluded.notes
        """,
        [args.experiment_name, args.factor_name, args.factor_unit,
         str(args.plate_info_csv), str(args.strain_info_csv), args.notes],
    )
    experiment_id = con.execute(
        "SELECT experiment_id FROM experiment WHERE experiment_name = ?", [args.experiment_name]
    ).fetchone()[0]
    print(f"experiment_id={experiment_id} for {args.experiment_name!r}")

    # --- Resolve run_number/configuration for every plate ------------------------
    plate_run_map = load_plate_run_map(args.plate_run_map) if args.plate_run_map else None

    resolved = []  # (plate_number, run_number, configuration, replicate_label, is_control, factor_value)
    for row in plate_df.iter_rows(named=True):
        plate_number = int(row[args.plate_number_column])
        raw_factor = row.get(args.factor_value_column)
        raw_replicate = row.get(args.replicate_column)
        is_control = (str(raw_factor) == args.control_value) or (str(raw_replicate) == args.control_value)
        try:
            factor_value = None if is_control else float(raw_factor)
        except (TypeError, ValueError):
            factor_value = None
            is_control = True

        if plate_run_map is not None:
            run_number, configuration = plate_run_map[plate_number]
        else:
            run_number, configuration = resolve_run_and_configuration(
                plate_number, base_run_number, args.plates_per_run, args.configs_per_run
            )

        resolved.append((plate_number, run_number, configuration,
                          str(raw_replicate), is_control, factor_value))

    # --- Validate resolved run numbers against strain metadata --------------------
    resolved_run_numbers = sorted({r[1] for r in resolved})
    unknown_but_present = [rn for rn in resolved_run_numbers if rn not in known_run_numbers]
    print(f"Resolved run numbers: {resolved_run_numbers}")
    if unknown_but_present:
        print(f"  {len(unknown_but_present)} run number(s) have no strain rows "
              f"(expected for control-only runs): {unknown_but_present}")

    # --- imager_run ----------------------------------------------------------------
    library_plate_by_run = {}
    for run_number in resolved_run_numbers:
        sub = strain_df.filter(pl.col("Run Number") == run_number)
        lp = None
        if sub.height > 0:
            lp_vals = sub["Library Plate"].drop_nulls().unique().to_list()
            lp = int(lp_vals[0]) if lp_vals else None
        library_plate_by_run[run_number] = lp
        con.execute(
            """
            INSERT INTO imager_run (run_number, experiment_id, library_plate, is_control_only)
            VALUES (?, ?, ?, ?)
            ON CONFLICT (run_number) DO UPDATE SET
                experiment_id = excluded.experiment_id,
                library_plate = excluded.library_plate,
                is_control_only = excluded.is_control_only
            """,
            [run_number, experiment_id, lp, run_number not in known_run_numbers],
        )

    # --- strain (global, upsert) --------------------------------------------------
    sync_strain(con, strain_df)

    # --- well_placement --------------------------------------------------------------
    well_rows = strain_df.select(
        pl.col("Run Number").alias("run_number"),
        pl.col("Configuration").alias("configuration"),
        pl.col("Well position").alias("well_position"),
        pl.col("Strain ID").alias("strain_id"),
        pl.col("Incubation Temp (°C)").alias("incubation_temp_c"),
        pl.col("Media").alias("media"),
    )
    dupes = (
        well_rows.group_by(["run_number", "configuration", "well_position"])
        .agg(pl.len().alias("n"))
        .filter(pl.col("n") > 1)
    )
    if dupes.height > 0:
        print("Duplicate (run_number, configuration, well_position) rows found:")
        print(dupes)
        raise SystemExit(1)
    con.execute("CREATE OR REPLACE TEMP TABLE _well_stage AS SELECT * FROM well_rows")
    con.execute(
        """
        INSERT INTO well_placement SELECT * FROM _well_stage
        ON CONFLICT (run_number, configuration, well_position) DO UPDATE SET
            strain_id = excluded.strain_id,
            incubation_temp_c = excluded.incubation_temp_c,
            media = excluded.media
        """
    )
    print(f"well_placement: {well_rows.height} rows staged")

    # --- condition_plate + condition_plate_factor -----------------------------------
    for plate_number, run_number, configuration, replicate_label, is_control, factor_value in resolved:
        con.execute(
            """
            INSERT INTO condition_plate (experiment_id, plate_number, replicate_label,
                                          is_control, run_number, configuration)
            VALUES (?, ?, ?, ?, ?, ?)
            ON CONFLICT (experiment_id, plate_number) DO UPDATE SET
                replicate_label = excluded.replicate_label,
                is_control = excluded.is_control,
                run_number = excluded.run_number,
                configuration = excluded.configuration
            """,
            [experiment_id, plate_number, replicate_label, is_control, run_number, configuration],
        )
        if factor_value is not None:
            con.execute(
                """
                INSERT INTO condition_plate_factor (experiment_id, plate_number, factor_name,
                                                      factor_value, factor_unit)
                VALUES (?, ?, ?, ?, ?)
                ON CONFLICT (experiment_id, plate_number, factor_name) DO UPDATE SET
                    factor_value = excluded.factor_value,
                    factor_unit = excluded.factor_unit
                """,
                [experiment_id, plate_number, args.factor_name, factor_value, args.factor_unit],
            )
    print(f"condition_plate: {len(resolved)} rows staged")

    con.close()
    print("Done.")


if __name__ == "__main__":
    main()
