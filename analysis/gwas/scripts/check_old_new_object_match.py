#!/usr/bin/env python3
"""Object-level comparison of old colony_measurement (pre-metal-replace DB backup) vs new
Copper rows: join on (image_name, object_label); compare Shape_Area and the strain
assigned to each object. Read-only. Run in SLURM."""
import sys, pathlib, duckdb, pandas as pd, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import duckdb_inputs as DI
REPO = DI.REPO
old = duckdb.connect(str(REPO/"db/rhodotorula_phenotypes.pre-metal-replace-20261007.duckdb"), read_only=True)
cols = [r[0] for r in old.execute("describe colony_measurement").fetchall()]
print("old colony_measurement cols w/ label/well:", [c for c in cols if "label" in c.lower() or "well" in c.lower()])
o = old.execute("""select cm.image_name, cm.object_label, cm.well_position, cm.Shape_Area as a_old, s.strain_code as old_strain,
                   i.plate_number, i.run_number, cf.factors['Copper concentration']::double as cu_old
                   from colony_measurement cm join image i using (image_name)
                   left join condition_plate cp on cp.run_number=i.run_number and cp.plate_number=i.plate_number
                   left join v_condition_factors cf on cf.experiment_id=cp.experiment_id and cf.plate_number=cp.plate_number
                   left join well_placement wp on wp.run_number=i.run_number and wp.configuration=cp.configuration and wp.well_position=cm.well_position
                   left join strain s on s.strain_id=wp.strain_id""").df()
print("old objects", len(o))
n = pd.read_parquet(REPO/"data/preprocessed/heavy_metal_array/Copper.parquet",
                    columns=["MetadataImage_ImageName","Shape_Area","strain_id","SAMPLE_NAME","Concentration",
                             "Grid_RowMajorIdx","Grid_ColMajorIdx","GridColIdx","Grid_RowNum","Grid_ColNum"])
n = n.rename(columns={"MetadataImage_ImageName":"image_name","Shape_Area":"a_new"})
print("old well_position range", o.well_position.min(), o.well_position.max(), "; new Grid_RowMajorIdx range", n.Grid_RowMajorIdx.min(), n.Grid_RowMajorIdx.max())
# restrict to wells with exactly one object per image in each table
o1 = o[~o.duplicated(["image_name","well_position"], keep=False)]
res = {}
for col in ["Grid_RowMajorIdx","Grid_ColMajorIdx","GridColIdx"]:
    n1 = n[~n.duplicated(["image_name", col], keep=False)]
    for off in (-1, 0, 1):
        j = o1.assign(k=o1.well_position + off).merge(n1, left_on=["image_name","k"], right_on=["image_name", col])
        res[(col, off)] = (len(j), np.isclose(j.a_old, j.a_new, rtol=1e-3).mean() if len(j) else np.nan)
        print(col, off, res[(col, off)])
best = max(res, key=lambda k: (res[k][1] if res[k][1]==res[k][1] else -1))
print("best key", best)
col, off = best
n1 = n[~n.duplicated(["image_name", col], keep=False)]
j = o1.assign(k=o1.well_position + off).merge(n1, left_on=["image_name","k"], right_on=["image_name", col])
print("matched wells", len(j), "area equal rate", np.isclose(j.a_old, j.a_new, rtol=1e-3).mean(), "corr", j.a_old.corr(j.a_new))
b = j[j.old_strain.notna() & j.SAMPLE_NAME.notna()]
print("wells with both strain labels:", len(b), " same strain string:", (b.old_strain==b.SAMPLE_NAME).mean())
bn = b.assign(on=b.old_strain.str.replace("_","-").str.lower(), nn=b.SAMPLE_NAME.str.replace("_","-").str.lower())
print("same after _->- normalisation:", (bn.on==bn.nn).mean())
print("concentration equal:", (j.cu_old==j.Concentration).mean())
mm = b[bn.on!=bn.nn].groupby(["old_strain","SAMPLE_NAME"]).size().sort_values(ascending=False)
print("top mismatched pairs:\n", mm.head(15).to_string()); print("n mismatched pairs:", len(mm))
