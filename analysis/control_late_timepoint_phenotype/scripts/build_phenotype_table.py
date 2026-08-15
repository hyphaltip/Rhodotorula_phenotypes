#!/usr/bin/env python3
"""Control-medium (Cu=0) late-timepoint phenotype table.

For each strain, take its latest image within hours_since_plate_start in
[80, 110] on control plates (Copper concentration = 0 mM) and aggregate the
per-colony trait values across all replicate colonies:

  colony size (Shape_Area, px) and CIELAB color (per-colony ColorLab_*Median)
  -> MEDIAN, MEAN, VARIANCE, SD per strain.

Output: results/phenotype_control_late_timepoint.csv
"""
import argparse
from pathlib import Path

import duckdb

DB = Path(__file__).resolve().parents[3] / "db" / "rhodotorula_phenotypes.duckdb"
OUT = Path(__file__).resolve().parents[1] / "results" / "phenotype_control_late_timepoint.csv"

QUERY = """
    with c0 as (
        select p.*, round(p.hours_since_plate_start, 0) as tp_h
        from v_phenotype p
        join imager_run ir using (run_number)
        join condition_plate_factor f
          on f.experiment_id = ir.experiment_id and f.plate_number = p.plate_number
        where f.factor_name = 'Copper concentration'
          and f.factor_value = 0
          and p.hours_since_plate_start between ? and ?
          and p.strain_id is not null
    ),
    latest as (
        select *, max(tp_h) over (partition by strain_id) as tp_hmax
        from c0
    )
    select
        s.strain_id,
        s.strain_code,
        s.strain_name,
        split_part(s.species, ' ', 1) as genus,
        s.species,
        s.origin,
        s.environment,
        count(*) as n_colonies,
        count(distinct run_number || ':' || well_position) as n_replicate_wells,
        max(tp_h) as timepoint_h,
        median(Shape_Area)                      as area_median,
        mean(Shape_Area)                        as area_mean,
        var_samp(Shape_Area)                    as area_var,
        stddev_samp(Shape_Area)                 as area_sd,
        median("ColorLab_L*Median")             as l_median,
        mean("ColorLab_L*Median")               as l_mean,
        var_samp("ColorLab_L*Median")           as l_var,
        stddev_samp("ColorLab_L*Median")        as l_sd,
        median("ColorLab_a*Median")             as a_median,
        mean("ColorLab_a*Median")               as a_mean,
        var_samp("ColorLab_a*Median")           as a_var,
        stddev_samp("ColorLab_a*Median")        as a_sd,
        median("ColorLab_b*Median")             as b_median,
        mean("ColorLab_b*Median")               as b_mean,
        var_samp("ColorLab_b*Median")           as b_var,
        stddev_samp("ColorLab_b*Median")        as b_sd
    from latest l
    join strain s using (strain_id)
    where l.tp_h = l.tp_hmax
    group by all
    order by s.species, s.strain_code
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tmin", type=float, default=80.0, help="window lower bound (h)")
    ap.add_argument("--tmax", type=float, default=110.0, help="window upper bound (h)")
    ap.add_argument("--db", type=Path, default=DB)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    con = duckdb.connect(str(args.db), read_only=True)
    n_total = con.execute("select count(*) from strain where strain_id is not null").fetchone()[0]
    missing = con.execute("""
        select count(*) from (
            select distinct p.strain_id
            from v_phenotype p
            join imager_run ir using (run_number)
            join condition_plate_factor f
              on f.experiment_id = ir.experiment_id and f.plate_number = p.plate_number
            where f.factor_name = 'Copper concentration' and f.factor_value = 0
              and p.hours_since_plate_start between ? and ?
        )
    """, [args.tmin, args.tmax]).fetchone()[0]
    df = con.execute(QUERY, [args.tmin, args.tmax]).df()
    con.close()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    n = len(df)
    n_ge3 = int((df.n_colonies >= 3).sum())
    n_low = n - n_ge3
    print(f"window {args.tmin}-{args.tmax}h | strains with late control image: {n}/{n_total}")
    print(f"strains with any control image in window: {missing} | with >=3 colonies at latest pass: {n_ge3} | flagged (<3): {n_low}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
