#!/usr/bin/env python3
"""Convert the per-metal ArrayedHeavyMetalScreen intermediate CSVs to Parquet.

One Parquet file per metal. Columns are kept exactly as in the source CSV,
except `Concentration`, which is cast to DOUBLE in every file (Cr is DOUBLE,
the others BIGINT in the source), so the files union cleanly in DuckDB.
`strain_id` stays VARCHAR: it holds numeric IDs, `Control-N` labels and NULLs.

Usage:
    python3 heavy_metal_csv_to_parquet.py --src-dir <raw dir with *.csv> --out-dir <dir>
"""
import argparse
import sys
from pathlib import Path

import duckdb

METALS = {  # CSV stem -> output name
    "chromiummeasurementmeta": "Chromium",
    "coppermeasurementmeta": "Copper",
    "ironmeasurementmeta": "Iron",
    "leadmeasurementmeta": "Lead",
    "zincmeasurementmeta": "Zinc",
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(f"SET threads={args.threads}")
    for stem, metal in METALS.items():
        src = args.src_dir / f"{metal}_{stem}.csv"
        if not src.exists():
            sys.exit(f"missing input: {src}")
        out = args.out_dir / f"{metal}.parquet"
        rd = f"read_csv('{src}', sample_size=-1)"
        n_csv = con.execute(f"SELECT count(*) FROM {rd}").fetchone()[0]
        con.execute(
            f"COPY (SELECT * REPLACE (CAST(Concentration AS DOUBLE) AS Concentration) FROM {rd}) "
            f"TO '{out}' (FORMAT parquet, COMPRESSION zstd)"
        )
        n_pq = con.execute(f"SELECT count(*) FROM read_parquet('{out}')").fetchone()[0]
        ncol = len(con.execute(f"DESCRIBE SELECT * FROM read_parquet('{out}')").fetchall())
        assert n_csv == n_pq, f"{metal}: row count CSV {n_csv} != parquet {n_pq}"
        print(f"{metal}\trows={n_pq}\tcols={ncol}\t{out.stat().st_size/1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
