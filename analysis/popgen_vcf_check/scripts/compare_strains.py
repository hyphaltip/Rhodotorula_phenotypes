#!/usr/bin/env python3
"""Match the phenotype strains (DuckDB strain_info) to the latest R. mucilaginosa DH4148-reference callset and its groups.
Inputs are read only. Name matching: upper case, keep letters and digits only, TF_CN -> TFCN."""
import re, subprocess, difflib, sys
from pathlib import Path
import pandas as pd

PG = Path("/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref")
OUT = Path("analysis/popgen_vcf_check/results"); OUT.mkdir(parents=True, exist_ok=True)
norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper().replace("TF_CN", "TFCN"))
def half(s):
    h = len(s) // 2
    return s[:h] if len(s) % 2 == 1 and s[:h] == s[h + 1:] else s     # VCF sample names are NAME_NAME
si = pd.read_csv("ignore/strain_info_now.tsv", sep="\t", dtype=str); si = si[si.is_control == "false"].copy(); si["k"] = si.sample_name.map(norm)
print(f"phenotype strains (non-control): {len(si)}; R. mucilaginosa: {(si.species == 'Rhodotorula mucilaginosa').sum()}; duplicate keys: {si.k.duplicated().sum()}")
# ---- callset samples from the VCF header (all strains) and group lists from the YAML
def vcf_samples(f):
    return subprocess.run(["bash", "-lc", f"module load bcftools >/dev/null 2>&1; bcftools query -l {f}"], capture_output=True, text=True).stdout.split()
allv = [half(s) for s in vcf_samples(PG / "results/RmucDH4148.all.qc.annotated.vcf.gz")]
print(f"samples in the all-strain VCF: {len(allv)} (first 3: {allv[:3]})")
groups, cur = {}, None
for line in open(PG / "population_sets.yaml"):
    m = re.match(r"^  ([A-Za-z_]+):\s*$", line); s = re.match(r"^    - (\S+)", line)
    if m: cur = m.group(1); groups[cur] = []
    elif s and cur: groups[cur].append(s.group(1))
groups["all"] = allv
print({g: len(v) for g, v in groups.items()})
for g, v in list(groups.items()):
    ks = [norm(x) for x in v]; dup = len(ks) - len(set(ks))
    if dup: print(f"  WARNING {g}: {dup} duplicate normalized names")
# ---- membership per phenotype strain
for g, v in groups.items(): si["in_" + g] = si.k.isin({norm(x) for x in v})
ex = pd.read_csv(PG / "results/variant_qc/excluded_strains.tsv", sep="\t", dtype=str); ex["k"] = ex.strain.map(norm)
si = si.merge(ex[["k", "decision", "reasons"]].rename(columns={"decision": "qc_decision", "reasons": "qc_reasons"}).drop_duplicates("k"), on="k", how="left")
iss = pd.read_csv(PG / "results/variant_qc/strain_identity_issues_2026-10-06.tsv", sep="\t", dtype=str); iss["k"] = iss.strain.map(norm)
agg = iss.groupby("k").agg(identity_issue_ids=("id", lambda x: ",".join(x)), identity_issue_categories=("category", lambda x: ";".join(sorted(set(x)))), identity_status=("status", lambda x: ";".join(sorted(set(x))))).reset_index()
si = si.merge(agg, on="k", how="left")
pe = pd.read_csv(PG / "results/variant_qc/rmuc_pheno_exclusions.tsv", sep="\t", dtype=str); pe["k"] = pe.strain.map(norm)
si = si.merge(pe[["k", "reason"]].rename(columns={"reason": "popgen_phenotype_note"}).drop_duplicates("k"), on="k", how="left")
cols = ["strain_id", "sample_name", "species", "species_source", "metals", "n_colony_observations"] + ["in_" + g for g in groups] + ["qc_decision", "qc_reasons", "identity_issue_ids", "identity_issue_categories", "identity_status", "popgen_phenotype_note"]
si[cols].to_csv(OUT / "phenotype_vs_popgen_strains.csv", index=False)
mu = si[si.species == "Rhodotorula mucilaginosa"]
print("\n== R. mucilaginosa phenotype strains (n=%d) by callset membership" % len(mu))
for g in groups: print(f"  in {g:20s} {mu['in_' + g].sum():4d}")
print("  in none of the groups:", (~mu[[("in_" + g) for g in groups]].any(axis=1)).sum())
print("\n== all phenotype strains (n=%d)" % len(si)); 
for g in groups: print(f"  in {g:20s} {si['in_' + g].sum():4d}")
# ---- species conflicts: our label vs popgen group
core = si[si.in_rmuc_core]
print("\nphenotype strains in rmuc_core by OUR species label:\n", core.species.value_counts().to_string())
print("\nR. mucilaginosa-labelled phenotype strains NOT in rmuc_core, by where they are:")
nc = mu[~mu.in_rmuc_core]
print("  in hybrid_diploids:", nc.in_hybrid_diploids.sum(), "| in all but no group:", (nc.in_all & ~nc.in_rmuc_core & ~nc.in_hybrid_diploids & ~nc.in_rmuc_core_outgroup & ~nc.in_rmuc_with_hybrids).sum(), "| not in all:", (~nc.in_all).sum())
print("  qc decisions:", nc.qc_decision.value_counts(dropna=False).to_dict())
# ---- issue categories among phenotyped R. mucilaginosa strains
print("\nidentity-issue categories among R. mucilaginosa phenotype strains in rmuc_core:")
x = mu[mu.in_rmuc_core & mu.identity_issue_categories.notna()]; print(x.identity_issue_categories.str.split(";").explode().value_counts().to_string()); print("  strains affected:", len(x))
# ---- unmatched names both ways
pgk = {norm(x): x for x in groups["all"]}; phk = set(si.k)
only_pg = [v for k, v in pgk.items() if k not in phk]
cand = {s: difflib.get_close_matches(norm(s), list(phk), n=1, cutoff=0.88) for s in only_pg}
print(f"\ncallset strains with no phenotype row: {len(only_pg)}; of those with a close name in the phenotype table: {sum(1 for v in cand.values() if v)}")
pd.DataFrame({"popgen_strain": only_pg, "close_phenotype_name_key": [cand[s][0] if cand[s] else "" for s in only_pg]}).to_csv(OUT / "popgen_strains_without_phenotype.csv", index=False)
mu_no = mu[~mu.in_all][["strain_id", "sample_name", "species_source"]]; print(f"R. mucilaginosa phenotype strains not found in the callset by name: {len(mu_no)}")
mu_no.assign(close_popgen_name=[next(iter(difflib.get_close_matches(norm(s), list(pgk), n=1, cutoff=0.85)), "") for s in mu_no.sample_name]).to_csv(OUT / "phenotype_mucilaginosa_not_in_callset.csv", index=False)
# ---- old VCF used for the GWAS (RmucY2510_v2, 422 samples)
old = [l.strip() for l in open("ignore/vcf_samples.txt")]; ok = {norm(x) for x in old}
print(f"\nold GWAS VCF (NRRL Y-2510 reference): {len(old)} samples; phenotype strains in it: {si.k.isin(ok).sum()}; R. mucilaginosa phenotype strains in it: {mu.k.isin(ok).sum()}")
print("R. mucilaginosa phenotype strains in the old VCF AND in new rmuc_core:", (mu.k.isin(ok) & mu.in_rmuc_core).sum(), "| old VCF only:", (mu.k.isin(ok) & ~mu.in_rmuc_core).sum(), "| new rmuc_core only:", (~mu.k.isin(ok) & mu.in_rmuc_core).sum())
print("old VCF samples found in the new all-strain callset:", sum(1 for x in old if norm(x) in pgk), "of", len(old))
