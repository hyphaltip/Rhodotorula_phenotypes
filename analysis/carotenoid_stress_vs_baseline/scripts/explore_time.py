import duckdb, sys
c = duckdb.connect('db/rhodotorula_phenotypes.duckdb', read_only=True)
q = lambda s: c.execute(s).fetchall()
print("cols present:", [r[0] for r in q("describe heavy_metal_measurement") if r[0] in ('sample_plate','plate_position','capture_datetime','ColorLab_a*GeoMedian','ColorLab_a*Medoid','Shape_Area','Grid_RowMajorIdx','temperature')])
# plate = (Metal, run_number, plate_position); time axis = hours since that plate's first image
c.execute("""create temp view t as select Metal, run_number, plate_position, sample_plate, Concentration, capture_datetime,
   epoch(capture_datetime) - min(epoch(capture_datetime)) over (partition by Metal, run_number, plate_position) as s0 from heavy_metal_measurement""")
print("plates (Metal,run,plate_position):", q("select Metal, count(*) from (select distinct Metal, run_number, plate_position from t) group by 1 order by 1"))
print("images per plate (distinct capture_datetime): ", q("select Metal, min(n), median(n), max(n) from (select Metal, run_number, plate_position, count(distinct capture_datetime) n from t group by 1,2,3) group by 1 order by 1"))
print("hours since first image, quantiles:", q("select Metal, round(min(s0)/3600,1), round(quantile_cont(s0,0.5)/3600,1), round(max(s0)/3600,1) from t group by 1 order by 1"))
print("conc levels per plate_position (is concentration a plate property?):", q("select Metal, max(k) from (select Metal, run_number, plate_position, count(distinct Concentration) k from t group by 1,2,3) group by 1 order by 1"))
print("hours of imaging round, histogram (Cu), 12h bins:", q("select floor(s0/3600/12)*12 b, count(*) from t where Metal='Copper' group by 1 order by 1"))
print("a* summary by metal:", q('select Metal, round(min("ColorLab_a*GeoMedian"),1), round(median("ColorLab_a*GeoMedian"),1), round(max("ColorLab_a*GeoMedian"),1) from heavy_metal_measurement group by 1 order by 1'))
print("a* vs area correlation (all rows) by metal:", q('select Metal, round(corr("ColorLab_a*GeoMedian", Shape_Area),3), round(corr("ColorLab_a*GeoMedian", ln(Shape_Area+1)),3) from heavy_metal_measurement group by 1 order by 1'))
