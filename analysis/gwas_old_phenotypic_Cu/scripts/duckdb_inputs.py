#!/usr/bin/env python3
"""Shared DuckDB input helpers for the GWAS (replaces Copper.Strain_info.csv).

Exports the strain table from the `strain_info` view of
db/rhodotorula_phenotypes.duckdb in the shape reconcile_strains.py expects
(one row per strain_id; `Strain` column = VCF-style sample ID, e.g. TFCN_17-291Y-1).

Mapping (see analysis/gwas/DUCKDB_RETOOL_PLAN.md):
  strain_info.sample_name -> Strain      (the ID matched to VCF sample IDs)
  strain_info.strain      -> Strain_name (lab name, e.g. '17-291Y-1 BY126-B8')
  strain_info.species     -> Species
  strain_info.location    -> Location
Control-N rows (is_control) are excluded: they have no sample_name.

The DB is opened read-only. If another process holds a write lock, the file is
copied to $SCRATCH (SLURM job-local) and the copy is opened instead.

Usage (inside a SLURM job or via sbatch --wrap):
  pixi run python analysis/gwas/scripts/duckdb_inputs.py --out analysis/gwas/results/strain_reconciliation/strain_info_from_duckdb.csv
Do not use BASH_SOURCE in callers; paths here come from this file's location.
"""
from __future__ import annotations

import argparse
import os
import pathlib
import shutil

import duckdb
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[3]
DB = REPO / "db" / "rhodotorula_phenotypes.duckdb"


def connect_readonly(db: pathlib.Path = DB) -> duckdb.DuckDBPyConnection:
    try:
        return duckdb.connect(str(db), read_only=True)
    except duckdb.IOException as e:
        if "lock" not in str(e).lower():
            raise
        scratch = os.environ.get("SCRATCH")
        if not scratch:
            raise RuntimeError("DB is write-locked and $SCRATCH is unset; run inside a SLURM job") from e
        snap = pathlib.Path(scratch) / f"{db.stem}.snapshot.duckdb"
        shutil.copy(db, snap)
        print(f"[duckdb_inputs] DB locked; using snapshot copy {snap}")
        return duckdb.connect(str(snap), read_only=True)


def strain_table(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    df = con.execute("""
        SELECT strain_id,
               sample_name AS "Strain",
               strain      AS "Strain_name",
               species     AS "Species",
               location    AS "Location",
               array_to_string(metals_tested, ';') AS metals_tested,
               n_colony_observations
        FROM strain_info
        WHERE NOT is_control
        ORDER BY try_cast(strain_id AS INTEGER)
    """).df()
    # Fail loudly on shape problems.
    assert len(df) > 0, "strain_info returned 0 non-control strains"
    assert df.strain_id.is_unique, "strain_id not unique in strain_info"
    assert df.Strain.notna().all(), "non-control strain with NULL sample_name"
    dup = df[df.Strain.duplicated(keep=False)]
    df["duplicate_sample_name"] = df.Strain.duplicated(keep=False)
    if len(dup):
        print(f"[duckdb_inputs] WARNING {len(dup)} rows share a sample_name:\n"
              f"{dup[['strain_id','Strain']].to_string(index=False)}")
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--db", default=str(DB))
    args = ap.parse_args()
    con = connect_readonly(pathlib.Path(args.db))
    df = strain_table(con)
    pathlib.Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"[duckdb_inputs] wrote {args.out}: {len(df)} strains, "
          f"{df.Strain.nunique()} distinct sample names")


if __name__ == "__main__":
    main()
