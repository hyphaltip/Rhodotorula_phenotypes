#!/usr/bin/env python3
"""Shared inputs for the stratified analyses: strain table (analysis group from the curated database), tree-tip mapping.
Writes results/strat/strain_table.csv (strain_id, sample_name, species = analysis group, species_db, ploidy_status, tip) and checks the join."""
import re, sys
from pathlib import Path
import pandas as pd, duckdb

R = Path("analysis/carotenoid_stress_vs_baseline/results/strat"); R.mkdir(parents=True, exist_ok=True)
norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper().replace("TF_CN", "TFCN"))
con = duckdb.connect("db/rhodotorula_phenotypes.duckdb", read_only=True)
st = con.execute("""select strain_id, sample_name, species as species_db, ploidy_status, hybrid_subgroup, clade_marker, gwas_panel
                    from strain_info where not is_control""").df()
# D-53: the analysis group `species` splits R. mucilaginosa by ploidy and clade marker. "Rhodotorula mucilaginosa" = pure haploid only.
MUC = "Rhodotorula mucilaginosa"
def grp(r):
    if r.species_db != MUC: return r.species_db
    if r.clade_marker == "aff_mucilaginosa": return "Rhodotorula aff. mucilaginosa"
    if r.ploidy_status == "diploid_hybrid": return "Rhodotorula mucilaginosa hybrid diploid"
    assert r.ploidy_status == "haploid", f"unexpected ploidy for {r.sample_name}: {r.ploidy_status}"
    return MUC
st["species"] = st.apply(grp, axis=1)
st["k"] = st.sample_name.map(norm)
print(f"strains: {len(st)}; analysis groups:\n{st.species.value_counts().to_string()}")
tips = [l.strip().replace(".proteins.fa", "").replace(".proteins", "") for l in open("ignore/tree_tips.txt")]
tipmap = {}
for t in tips:
    m = re.match(r"(Rhodotorula_[a-z]+(?:_sp\._clade_[A-Z]+|_sp_clade_[A-Z]+)?|[A-Z][a-z]+_[a-z]+)_(.*)", t)
    if m: tipmap.setdefault(norm(m.group(2)), []).append(t)
dup = {k: v for k, v in tipmap.items() if len(v) > 1}
st["tip"] = st.k.map(lambda k: tipmap[k][0] if k in tipmap and len(tipmap[k]) == 1 else None)
print(f"tree tips parsed: {len(tips)}; strains with a unique tip: {st.tip.notna().sum()}; ambiguous tip keys: {len(dup)}")
st.drop(columns="k").to_csv(R / "strain_table.csv", index=False)
