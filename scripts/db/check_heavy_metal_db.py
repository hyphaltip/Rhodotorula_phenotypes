#!/usr/bin/env python3
"""Read-only sanity check of heavy_metal_measurement and the strain_info view."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import duckdb  # noqa: E402

DB = Path(__file__).resolve().parents[2] / "db" / "rhodotorula_phenotypes.duckdb"
c = duckdb.connect(str(DB), read_only=True)
q = lambda s: c.execute(s).fetchall()
print("tables/views:", q("select table_name, table_type from information_schema.tables order by 1"))
print("heavy_metal_measurement rows, cols:", q("select count(*) from heavy_metal_measurement")[0][0],
      len(q("describe heavy_metal_measurement")))
print("rows per metal:", q("select Metal, count(*) from heavy_metal_measurement group by 1 order by 1"))
print("strain_info (strains, controls, non-control w/o species, non-control w/o sample_name):",
      q("select count(*), count(*) filter (is_control), count(*) filter (species is null and not is_control), "
        "count(*) filter (sample_name is null and not is_control) from strain_info")[0])
print("strain_info sample:", q("select * from strain_info where not is_control order by try_cast(strain_id as integer) limit 3"))
print("species counts:", q("select species, count(*) from strain_info where not is_control group by 1 order by 2 desc limit 8"))
print("n metals tested per strain:", q("select len(metals_tested) n, count(*) from strain_info group by 1 order by 1"))
print("duplicate strain_id rows:", q("select count(*) from (select strain_id from strain_info group by 1 having count(*)>1)")[0][0])
print("raw MB in the 6 redundant strain text columns of the fact table:",
      q("select sum(strlen(SAMPLE_NAME)+strlen(STRAIN)+strlen(SPECIES)+strlen(MS2_SAMPLE_Cell)"
        "+strlen(MS2_SAMPLE_Supernatant)+strlen(Location))/1e6 from heavy_metal_measurement")[0][0])
