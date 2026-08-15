#!/usr/bin/env python3
"""Import per-colony measurement Parquet files into `colony_measurement`.

Idempotent, file-level, safe to re-run after new files land in
data/preprocessed/. See DATABASE_DESIGN.md §5.

For each file:
  - Skip if already imported with the same size/mtime.
  - If already imported with a *different* size/mtime, treat it as a
    corrected re-delivery: delete its existing rows and reimport.
  - Otherwise import it: parse the filename, compute well_position, and
    insert the image + colony_measurement rows in a single transaction, so a
    crash mid-file never leaves the two tables inconsistent -- the next run
    just retries the whole file.
  - `colony_measurement` is created from the first file's Parquet schema if
    it doesn't exist yet; later files that introduce new columns trigger an
    ALTER TABLE ADD COLUMN (announced, not silent).

Usage:
    python3 20_import_measurements.py --glob "data/preprocessed/*.parquet"
"""

import argparse
import glob as globmod
import sys
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.db import get_connection  # noqa: E402
from lib.imagename import parse_image_name  # noqa: E402

# Metadata_Dataset is a random per-run temp identifier with no analytical
# value (see PhenotypicMeasurements/Copper/Scripts/copper_metadata.py, which
# drops it the same way) -- excluded from colony_measurement.
DROPPED_SOURCE_COLUMNS = {"Metadata_Dataset", "Metadata_ImageName"}


def colony_measurement_exists(con) -> bool:
    return con.execute(
        "SELECT count(*) FROM information_schema.tables WHERE table_name = 'colony_measurement'"
    ).fetchone()[0] > 0


def create_colony_measurement_from(con, parquet_path: Path) -> None:
    cols = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{parquet_path}')").fetchall()
    col_defs = []
    for name, dtype, *_ in cols:
        if name in DROPPED_SOURCE_COLUMNS:
            continue
        col_defs.append(f'"{name}" {dtype}')
    col_defs.append('"image_name" VARCHAR NOT NULL REFERENCES image(image_name)')
    col_defs.append('"object_label" INTEGER NOT NULL')
    col_defs.append('"well_position" INTEGER')
    ddl = f"CREATE TABLE colony_measurement ({', '.join(col_defs)}, PRIMARY KEY (image_name, object_label))"
    con.execute(ddl)
    print("Created colony_measurement from", parquet_path.name)


def sync_columns(con, parquet_path: Path) -> None:
    file_cols = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{parquet_path}')").fetchall()
    existing = {r[0] for r in con.execute("DESCRIBE colony_measurement").fetchall()}
    for name, dtype, *_ in file_cols:
        if name in DROPPED_SOURCE_COLUMNS:
            continue
        if name not in existing:
            con.execute(f'ALTER TABLE colony_measurement ADD COLUMN IF NOT EXISTS "{name}" {dtype}')
            print(f"  schema drift: added column {name!r} ({dtype}) from {parquet_path.name}")


def import_one_file(con, path: Path, source_glob_root: Path) -> str:
    image_name = path.stem
    stat = path.stat()
    size, mtime = stat.st_size, stat.st_mtime

    existing = con.execute(
        "SELECT file_size_bytes, file_mtime FROM image WHERE image_name = ?", [image_name]
    ).fetchone()
    if existing is not None:
        existing_size, existing_mtime = existing
        if existing_size == size and abs((existing_mtime.timestamp() if existing_mtime else 0) - mtime) < 1:
            return "skipped"
        print(f"  {image_name}: size/mtime changed since last import -- treating as corrected re-delivery")
        # Delete the old rows in their own committed transactions, then fall
        # through to the normal insert path below as if this were a new
        # file. The child (colony_measurement) and parent (image) deletes
        # must each be in their *own* committed transaction: deleting both
        # within a single transaction hits a DuckDB limitation where an FK
        # parent row deleted earlier in a transaction is not recognized as
        # clear for a later statement in that same transaction (DuckDB's own
        # "foreign key limitations" docs, referenced directly in the error
        # DuckDB raises for this case). This trades perfect atomicity on the
        # rare re-delivery path for correctness: if the process dies between
        # these transactions, the next run still sees a stale size/mtime and
        # safely retries the same re-delivery cleanup + reimport from
        # scratch (each step is independently idempotent).
        con.execute("BEGIN TRANSACTION")
        con.execute("DELETE FROM colony_measurement WHERE image_name = ?", [image_name])
        con.execute("COMMIT")
        con.execute("BEGIN TRANSACTION")
        con.execute("DELETE FROM image WHERE image_name = ?", [image_name])
        con.execute("COMMIT")

    meta = parse_image_name(image_name)

    if not colony_measurement_exists(con):
        create_colony_measurement_from(con, path)
    else:
        sync_columns(con, path)

    if con.execute("SELECT count(*) FROM imager_run WHERE run_number = ?", [meta["run_number"]]).fetchone()[0] == 0:
        raise RuntimeError(
            f"{image_name}: run_number {meta['run_number']} is not registered in imager_run -- "
            "run 10_import_experiment.py for this experiment first"
        )

    file_cols = [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM read_parquet('{path}')").fetchall()
                 if r[0] not in DROPPED_SOURCE_COLUMNS]
    select_cols = ", ".join(f'"{c}"' for c in file_cols)

    con.execute("BEGIN TRANSACTION")
    try:
        con.execute(
            """
            INSERT INTO image (image_name, run_number, plate_number, temperature_token,
                                imaged_at, source_path, file_size_bytes, file_mtime)
            VALUES (?, ?, ?, ?, ?, ?, ?, to_timestamp(?))
            """,
            [image_name, meta["run_number"], meta["plate_number"], meta["temperature_token"],
             meta["imaged_at"], str(path), size, mtime],
        )
        con.execute(
            f"""
            INSERT INTO colony_measurement
            SELECT {select_cols},
                   ? AS image_name,
                   "ObjectLabel" AS object_label,
                   CAST("Grid_ColMajorIdx" AS INTEGER) + 1 AS well_position
            FROM read_parquet(?)
            """,
            [image_name, str(path)],
        )
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    return "imported"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=None)
    ap.add_argument("--glob", default="data/preprocessed/*.parquet")
    args = ap.parse_args()

    con = get_connection(args.db)
    con.execute(open(Path(__file__).resolve().parent / "00_init_schema.sql").read())

    root = Path(".").resolve()
    files = sorted(Path(p) for p in globmod.glob(args.glob))
    if not files:
        print(f"No files matched glob {args.glob!r}")
        return

    counts = {"imported": 0, "skipped": 0, "error": 0}
    for i, path in enumerate(files, 1):
        try:
            result = import_one_file(con, path, root)
            counts[result] += 1
        except Exception as e:
            counts["error"] += 1
            print(f"  ERROR importing {path.name}: {e}")
        if i % 200 == 0:
            print(f"...{i}/{len(files)} files processed ({counts})")

    print(f"Done. {counts['imported']} imported, {counts['skipped']} already up to date, "
          f"{counts['error']} errors, out of {len(files)} files.")
    con.close()


if __name__ == "__main__":
    main()
