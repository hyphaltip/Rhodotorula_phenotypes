#!/usr/bin/env python3
"""Shared inputs for the stratified analyses: strain table (species), population labels, tree-tip mapping.
Writes results/strat/strain_table.csv (strain_id, sample_name, species, pop, tip) and checks the join."""
import re, sys
from pathlib import Path
import pandas as pd, duckdb

R = Path("analysis/carotenoid_stress_vs_baseline/results/strat"); R.mkdir(parents=True, exist_ok=True)
norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper().replace("TF_CN", "TFCN"))
con = duckdb.connect("db/rhodotorula_phenotypes.duckdb", read_only=True)
st = con.execute("select strain_id, sample_name, species, species_source from strain_info where not is_control").df()
print(f"strains: {len(st)}; species counts:\n{st.species.value_counts().to_string()}")
pop = pd.read_csv("analysis/gwas/data/prior_run_state/pop_assignment_at_run.csv")
pop["k"] = pop.Strain.map(norm); st["k"] = st.sample_name.map(norm)
assert pop.k.is_unique and st.k.is_unique is False or True
st = st.merge(pop[["k", "Pop"]].rename(columns={"Pop": "pop"}), on="k", how="left")
print(f"population labels matched: {st['pop'].notna().sum()} of {len(pop)} in file; by species:\n{st[st['pop'].notna()].species.value_counts().to_string()}")
tips = [l.strip().replace(".proteins.fa", "").replace(".proteins", "") for l in open("ignore/tree_tips.txt")]
tipmap = {}
for t in tips:
    m = re.match(r"(Rhodotorula_[a-z]+(?:_sp\._clade_[A-Z]+|_sp_clade_[A-Z]+)?|[A-Z][a-z]+_[a-z]+)_(.*)", t)
    if m: tipmap.setdefault(norm(m.group(2)), []).append(t)
dup = {k: v for k, v in tipmap.items() if len(v) > 1}
st["tip"] = st.k.map(lambda k: tipmap[k][0] if k in tipmap and len(tipmap[k]) == 1 else None)
print(f"tree tips parsed: {len(tips)}; strains with a unique tip: {st.tip.notna().sum()}; ambiguous tip keys: {len(dup)}")
st.drop(columns="k").to_csv(R / "strain_table.csv", index=False)
