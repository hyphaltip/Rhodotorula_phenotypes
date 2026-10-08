import duckdb
c = duckdb.connect('db/rhodotorula_phenotypes.duckdb', read_only=True)
q = lambda s: c.execute(s).fetchall()
base = "from heavy_metal_measurement where strain_id is not null and strain_id not like 'Control%'"
print("grid/well columns:", [r[0] for r in q("describe heavy_metal_measurement") if r[0].startswith(('Grid','sample_plate','plate_position','Object','Location'))])
# one well = (Metal, run, plate_position, Grid_RowNum, Grid_ColNum)
print("\nper (Metal,strain,conc): n distinct runs / plates / wells -- summary of distribution")
for r in q(f"""select Metal, min(nw), quantile_cont(nw,.25), median(nw), quantile_cont(nw,.75), max(nw),
                      median(np), median(nr) from (
   select Metal, strain_id, Concentration,
          count(distinct (run_number, plate_position, Grid_RowNum, Grid_ColNum)) nw,
          count(distinct (run_number, plate_position)) np, count(distinct run_number) nr
   {base} group by 1,2,3) group by 1 order by 1"""): print(r)
print("\nwells per plate-position per strain (is a strain replicated within a plate?):")
print(q(f"""select Metal, median(k), max(k) from (select Metal, strain_id, run_number, plate_position, count(distinct (Grid_RowNum, Grid_ColNum)) k {base} group by 1,2,3,4) group by 1 order by 1"""))
print("\nhow many wells at 0 dose vs max dose per strain, Cu:")
print(q(f"""select Concentration, median(nw), min(nw), max(nw) from (select strain_id, Concentration, count(distinct (run_number, plate_position, Grid_RowNum, Grid_ColNum)) nw {base} and Metal='Copper' group by 1,2) group by 1 order by 1"""))
print("\nobjects per well per image (duplicates per well at a timepoint?):", q(f"select Metal, max(k) from (select Metal, run_number, plate_position, Grid_RowNum, Grid_ColNum, capture_datetime, count(*) k {base} group by all) group by 1 order by 1"))
print("\nrun x strain x dose: do replicates come from different runs? strains seen in >1 run at same dose, Cu:", q(f"select count(*), count(*) filter (nr>1) from (select strain_id, Concentration, count(distinct run_number) nr {base} and Metal='Copper' group by 1,2)"))
