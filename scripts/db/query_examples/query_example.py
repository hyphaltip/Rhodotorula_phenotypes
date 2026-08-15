#!/usr/bin/env python3
"""Minimal example: query the phenotype database from Python.

Run from the repo root:
    pixi run python3 scripts/db/query_examples/query_example.py
"""

import duckdb

DB_PATH = "db/rhodotorula_phenotypes.duckdb"


def main():
    con = duckdb.connect(DB_PATH, read_only=True)

    # 1. Every colony observation for one strain, across every condition and
    #    replicate it appears in.
    strain_code = "TFCN_152C-3"
    df = con.execute(
        """
        SELECT strain_code, factors, replicate_label, hours_since_plate_start,
               Shape_Area, Intensity_MeanIntensity
        FROM v_phenotype
        WHERE strain_code = ?
        ORDER BY hours_since_plate_start
        """,
        [strain_code],
    ).pl()
    print(f"{strain_code}: {len(df)} colony observations")
    print(df.head())

    # 2. Mean colony size per Copper concentration for that strain (pull the
    #    level straight out of the `factors` MAP).
    summary = con.execute(
        """
        SELECT factors['Copper concentration'] AS concentration_mm,
               COUNT(*) AS n, AVG(Shape_Area) AS mean_area
        FROM v_phenotype
        WHERE strain_code = ? AND experiment_name = 'Copper'
        GROUP BY 1
        ORDER BY 1
        """,
        [strain_code],
    ).pl()
    print("\nMean colony area by Copper concentration:")
    print(summary)

    con.close()


if __name__ == "__main__":
    main()
