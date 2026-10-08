#!/usr/bin/env python3
"""Write schema.yaml and summary_stats.md for heavy-metal-array-intermediate from the Parquet files."""
import re
import sys
from pathlib import Path

import duckdb

PQ = "data/preprocessed/heavy_metal_array/*.parquet"
OUT = Path("data/metadata/heavy-metal-array-intermediate")
OUT.mkdir(parents=True, exist_ok=True)
con = duckdb.connect()
con.execute(f"CREATE VIEW hm AS SELECT * FROM read_parquet('{PQ}', union_by_name=true)")
cols = con.execute("DESCRIBE hm").fetchall()
metals = ["Chromium", "Copper", "Iron", "Lead", "Zinc"]
present = {m: {r[0] for r in con.execute(
    f"DESCRIBE SELECT * FROM read_parquet('data/preprocessed/heavy_metal_array/{m}.parquet')").fetchall()}
    for m in metals}

GROUPS = [
    (r"^Metadata", "image/dataset metadata (MetadataExperiment_Dataset is a per-run temp ID, absent from Zinc)"),
    (r"^Size_", "object size (Pb/Zn only)"),
    (r"^Shape_", "colony shape (CellProfiler-style morphology)"),
    (r"^Intensity_", "grayscale intensity statistics"),
    (r"^ColorLab_", "CIE Lab color, medoid/geo-median, DeltaE2000 spread"),
    (r"^ColorHSV_", "HSV color, robust mean"),
    (r"^Texture_", "Haralick texture, 4 angles at scale 05"),
    (r"^Bbox_", "bounding box / weighted centres, pixel coordinates"),
    (r"^Grid_|^GridColIdx$", "plate grid position of the colony"),
    (r"^GridLinReg_", "linear-regression grid fit (Pb/Zn only)"),
    (r"^SymZones_", "symmetry zones: core/dense/sparse radii and areas (Pb/Zn only)"),
    (r"^OrientZones_", "outward-rotation orientation metrics by zone (Pb/Zn only)"),
    (r"^(run_number|temperature|plate_position|date|time|capture_datetime|sample_plate|sample_plate_arrangement|Metal|Concentration|strain_id|SAMPLE_NAME|STRAIN|SPECIES|MS2_SAMPLE_Cell|MS2_SAMPLE_Supernatant|Location|Object_Label|index)$",
     "experiment design / identity"),
]
def group(c):
    for pat, d in GROUPS:
        if re.search(pat, c):
            return d
    return "other"

lines = ["# Schema for heavy-metal-array-intermediate",
         "# One row per segmented colony per image. 5 Parquet files (one per metal), unioned by column name.",
         "dataset: heavy-metal-array-intermediate",
         "row_unit: one colony object (MetadataImage_ImageName + Object_Label is unique within each metal)",
         "units_note: Concentration unit is NOT stated in the source. Cu/Fe/Pb/Zn span 0-30; Cr spans 0-1.2. Do not pool Cr with the others without confirming units.",
         "strain_id_note: VARCHAR. Numeric strain IDs, 'Control-N' labels, and NULL (Cr 5923, Cu 3489, Pb 1509 rows).",
         "columns:"]
for name, typ, *_ in cols:
    where = [m for m in metals if name in present[m]]
    lines.append(f'  - name: "{name}"')
    lines.append(f"    type: {typ}")
    lines.append(f"    group: {group(name)}")
    if len(where) < len(metals):
        lines.append(f"    present_in: [{', '.join(where)}]")
(OUT / "schema.yaml").write_text("\n".join(lines) + "\n")

q = lambda s: con.execute(s).fetchall()
md = ["# Summary statistics: heavy-metal-array-intermediate", "",
      f"Total rows: {q('select count(*) from hm')[0][0]:,}; columns in union: {len(cols)}", "",
      "| Metal | rows | cols in file | runs | strains (distinct strain_id) | NULL strain_id | conc levels | conc min-max | first date | last date |",
      "|---|---|---|---|---|---|---|---|---|---|"]
for m in metals:
    r = q(f"""select count(*), count(distinct run_number), count(distinct strain_id),
              count(*) filter (strain_id is null), count(distinct Concentration),
              min(Concentration), max(Concentration), min(date), max(date)
              from read_parquet('data/preprocessed/heavy_metal_array/{m}.parquet')""")[0]
    md.append(f"| {m} | {r[0]:,} | {len(present[m])} | {r[1]} | {r[2]} | {r[3]:,} | {r[4]} | {r[5]:g}-{r[6]:g} | {r[7]} | {r[8]} |")
md += ["", "## Runs per metal", ""]
for m, rn, n in q("select Metal, run_number, count(*) from hm group by 1,2 order by 1,2"):
    md.append(f"- {m} {rn}: {n:,}")
md += ["", "## Missing values in key columns (all metals)", "",
       "| column | NULL count |", "|---|---|"]
for c in ["Shape_Area", "ColorLab_L*GeoMedian", "Intensity_MeanIntensity", "strain_id", "Concentration", "capture_datetime"]:
    md.append(f"| {c} | {q(f'select count(*) filter (\"{c}\" is null) from hm')[0][0]:,} |")
md += ["", "## Columns only in Pb and Zn (32)", ""]
only = sorted(c for c in present['Lead'] if c not in present['Copper'])
md.append(", ".join(f"`{c}`" for c in only))
(OUT / "summary_stats.md").write_text("\n".join(md) + "\n")
print("ok", len(cols), "columns")
