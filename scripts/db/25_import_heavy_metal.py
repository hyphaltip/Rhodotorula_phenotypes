#!/usr/bin/env python3
"""Load the per-metal heavy-metal Parquet files into `heavy_metal_measurement`.

Replaces the legacy `colony_measurement` table (see .living/decisions.md D-34/D-35).
One table, all metals, unioned by column name: 147 columns shared by every metal,
32 present only for Lead/Zinc (NULL elsewhere). `Metal` identifies the source file.
Rebuilds the table from scratch on every run (CREATE OR REPLACE), so it is idempotent.

Usage:
    python3 25_import_heavy_metal.py [--glob "data/preprocessed/heavy_metal_array/*.parquet"] [--drop-legacy]

Also (re)creates the `strain_info` view (35_create_strain_view.sql).

--drop-legacy drops `colony_measurement` and the views built on it
(v_phenotype, v_strain_experiment_summary, v_growth_timeseries). Back up the DB first.
"""
import argparse
import glob as globmod
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.db import get_connection  # noqa: E402

LEGACY_VIEWS = ["v_growth_timeseries", "v_strain_experiment_summary", "v_phenotype"]  # dependents first


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", default="data/preprocessed/heavy_metal_array/*.parquet")
    ap.add_argument("--db", default=None)
    ap.add_argument("--drop-legacy", action="store_true")
    args = ap.parse_args()
    files = sorted(globmod.glob(args.glob))
    if not files:
        sys.exit(f"no files match {args.glob}")
    con = get_connection(args.db)
    expected = sum(con.execute(f"SELECT count(*) FROM read_parquet('{f}')").fetchone()[0] for f in files)
    con.execute("BEGIN")
    con.execute(
        f"CREATE OR REPLACE TABLE heavy_metal_measurement AS "
        f"SELECT * FROM read_parquet('{args.glob}', union_by_name=true)"
    )
    n = con.execute("SELECT count(*) FROM heavy_metal_measurement").fetchone()[0]
    assert n == expected, f"loaded {n} rows, Parquet files hold {expected}"
    if args.drop_legacy:
        for v in LEGACY_VIEWS:
            con.execute(f"DROP VIEW IF EXISTS {v}")
        con.execute("DROP TABLE IF EXISTS colony_measurement")
    con.execute("COMMIT")
    con.execute((Path(__file__).resolve().parent / "35_create_strain_view.sql").read_text())
    ncol = len(con.execute("DESCRIBE heavy_metal_measurement").fetchall())
    print(f"heavy_metal_measurement: {n:,} rows, {ncol} columns from {len(files)} files")
    for metal, k in con.execute("SELECT Metal, count(*) FROM heavy_metal_measurement GROUP BY 1 ORDER BY 1").fetchall():
        print(f"  {metal}: {k:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
