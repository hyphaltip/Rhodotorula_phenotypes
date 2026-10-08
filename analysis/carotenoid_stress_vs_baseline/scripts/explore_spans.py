import duckdb
c = duckdb.connect('db/rhodotorula_phenotypes.duckdb', read_only=True)
q = lambda s: c.execute(s).fetchall()
c.execute("""create temp view p as select Metal, run_number, plate_position, any_value(Concentration) conc, count(distinct capture_datetime) n_img,
  epoch(max(capture_datetime))-epoch(min(capture_datetime)) span_s, count(distinct strain_id) n_strains
  from heavy_metal_measurement group by 1,2,3""")
print("plate span hours (min, p10, median, p90, max) and n plates by metal:")
for r in q("select Metal, count(*), round(min(span_s)/3600,1), round(quantile_cont(span_s,.1)/3600,1), round(median(span_s)/3600,1), round(quantile_cont(span_s,.9)/3600,1), round(max(span_s)/3600,1) from p group by 1 order by 1"): print(r)
print("plates by metal x run x span bucket (hours, 12h):")
for r in q("select Metal, run_number, floor(span_s/3600/12)*12 b, count(*), list_sort(list(distinct conc)) from p group by 1,2,3 order by 1,2,3"): print(r)
