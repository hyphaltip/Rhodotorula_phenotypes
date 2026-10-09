#!/usr/bin/env python3
"""Generate the strain curation table from the checksummed popgen snapshot and the decisions in .living/decisions.md.

Outputs (data/metadata/strain-curation/):
  strain_curation.csv          one row per phenotyped strain (the table the DuckDB build ingests)
  strain_curation_changes.csv  every species / ploidy / marker change with its reason and decision id
Run as a SLURM job from the repo root (reads the Parquet copies of the screen tables):
  pixi run python scripts/curation/build_strain_curation.py

Every rule below carries its decision id. To revisit a decision, edit the rule block and rerun; do not edit the CSV.
"""
import re
import sys
import datetime
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

SNAP = Path("data/raw/popgen-callset-metadata")
POP = SNAP / "popgen"
OUT = Path("data/metadata/strain-curation")
CUTOFFS = [2, 5, 10, 20, 50]          # D-49: stored for every cutoff so the cutoff can be changed by filtering a column
PANEL_CUTOFF = 5                      # D-49: current choice
EXPECTED_PANEL = 126                  # D-49 result; the script stops if the rebuild disagrees
TODAY = datetime.date.today().isoformat()

# ---------------------------------------------------------------- decision rules
MUC = "Rhodotorula mucilaginosa"
# D-45: species for 'Species Not Found' strains that have genomes = the sourmash read call (Species_ID_db calls_reads.tsv).
D45 = {"TFCN_1A-1-5": MUC, "TFCN_2M-1-3": MUC, "TFCN_86C-3": MUC,
       "TFCN_7-6-3": "Rhodotorula diobovata", "TFCN_211C-2": "Rhodotorula aff. babjevae"}
# D-46: popgen species labels that override our 'R. mucilaginosa' label (popgen strain_qc.tsv 'species' column).
D46_RENAME = {"R. frigidialcoholis": "Rhodotorula frigidialcoholis"}          # species changes
D46_MARKER = {"R. aff. mucilaginosa": "aff_mucilaginosa"}                      # species stays; marker set
# D-54: the reference strain has no callset sample.
REFERENCE = "DH4148"

norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper().replace("TF_CN", "TFCN"))
log = lambda m: print(m, flush=True)


def read_yaml_groups(path):
    groups, cur = {}, None
    for line in open(path):
        m = re.match(r"^  ([A-Za-z_]+):\s*$", line)
        s = re.match(r"^    - (\S+)", line)
        if m:
            cur = m.group(1)
            groups[cur] = []
        elif s and cur:
            groups[cur].append(s.group(1))
    return groups


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # ---- 1. screen strains (base table) -----------------------------------
    con = duckdb.connect()
    base = con.execute("""select strain_id, max(SAMPLE_NAME) as sample_name, max(STRAIN) as strain_label, max(SPECIES) as species_screen
                          from read_parquet('data/preprocessed/heavy_metal_array/*.parquet', union_by_name=true)
                          where strain_id is not null and strain_id not like 'Control%' group by strain_id""").df()
    base = base.sort_values("strain_id", key=lambda s: s.astype(int)).reset_index(drop=True)
    base["k"] = base.sample_name.map(norm)
    log(f"screen strains: {len(base)}")
    d = base.copy()
    d["species"] = d.species_screen
    d["species_source"] = "screen data"
    d["clade_marker"] = ""
    notes = {i: [] for i in d.index}
    dec = {i: set() for i in d.index}
    changes = []

    def change(i, field, old, new, reason, decision, evidence=""):
        changes.append(dict(strain_id=d.at[i, "strain_id"], sample_name=d.at[i, "sample_name"], field=field, old_value=old, new_value=new,
                            reason=reason, decision_id=decision, evidence=evidence))
        dec[i].add(decision)

    # ---- 2. earlier user-confirmed species overrides (D-37), folded in ------
    ov = pd.read_csv("data/metadata/heavy-metal-array-intermediate/strain_species_overrides.tsv", sep="\t", dtype=str)
    for _, r in ov.iterrows():
        i = d.index[d.strain_id == str(r.strain_id)][0]
        change(i, "species", d.at[i, "species"], r.species, r.note, "D-37", r.source)
        d.at[i, "species"], d.at[i, "species_source"] = r.species, "user confirmation (D-37)"

    # ---- 3. popgen inputs ----------------------------------------------------
    qc = pd.read_csv(POP / "variant_qc/strain_qc.tsv", sep="\t", dtype=str)
    qc["k"] = qc.strain.map(norm)
    in_callset = set(qc.k)
    log(f"popgen callset strains (strain_qc.tsv): {len(qc)}")
    ploidy = pd.read_csv(POP / "ploidy_overrides.csv", dtype=str)
    ploidy["k"] = ploidy.strain.map(norm)
    hyb = pd.read_csv(POP / "variant_qc/hybrid_diploids.tsv", sep="\t", dtype=str)
    hyb["k"] = hyb.strain.map(norm)
    groups = read_yaml_groups(POP / "population_sets.yaml")
    gk = {g: {norm(x) for x in v} for g, v in groups.items()}
    iss = pd.read_csv(POP / "variant_qc/strain_identity_issues_2026-10-06.tsv", sep="\t", dtype=str)
    iss["k"] = iss.strain.map(norm)
    clone = pd.read_csv(POP / "variant_qc/divergence/clone_groups_1e-3.tsv", sep="\t", dtype=str)
    clone["k"] = clone.strain.map(norm)
    ni = pd.read_csv(POP / "variant_qc/declone/rmuc_core.groups_le5.tsv", sep="\t", dtype=str)
    ni["k"] = ni.strain.map(norm)
    calls = pd.read_csv(SNAP / "species_id_db/calls_reads.tsv", sep="\t", dtype=str)
    calls["name"] = calls["query"].str.split("|").str[-1]
    popname = {k: s for k, s in zip(qc.k, qc.strain)}                    # matched popgen spelling

    d["popgen_strain"] = d.k.map(popname).fillna("")
    d["in_callset"] = d.k.isin(in_callset)
    for g in ("rmuc_core", "rmuc_core_outgroup", "rmuc_with_hybrids", "hybrid_diploids"):
        d["in_" + g] = d.k.isin(gk[g])
    d["clone_group"] = d.k.map(dict(zip(clone.k, clone.clone_group))).fillna("")
    d["near_identical_group"] = d.k.map(dict(zip(ni.near_identical_group.index.map(lambda _: None), []))).fillna("")  # placeholder, set below
    d["near_identical_group"] = d.k.map(dict(zip(ni.k, ni.near_identical_group))).fillna("")

    # ---- 4. species rules (D-45, D-46) ---------------------------------------
    for sname, newsp in D45.items():
        k = norm(sname)
        rows = d.index[d.k == k]
        assert len(rows) == 1, f"D-45 strain {sname}: {len(rows)} rows"
        i = rows[0]
        c = calls[calls.name.map(norm) == k]
        assert len(c) == 1, f"D-45 strain {sname}: expected one sourmash call, found {len(c)}"
        called = c.best_label.iloc[0].replace("aff ", "aff. ")
        assert called == newsp, f"D-45 rule says {newsp} for {sname} but sourmash calls {called}"
        ev = f"sourmash read call {called}, ANI {float(c.best_ani.iloc[0]):.4f}, nearest ref {c.best_ref.iloc[0]}"
        change(i, "species", d.at[i, "species"], newsp, "species from the sourmash call of the single 9003 library (unverified batch)", "D-45", ev)
        d.at[i, "species"], d.at[i, "species_source"] = newsp, "sourmash read call (D-45)"
        notes[i].append(f"D-45: {ev}; single 9003-batch library, identity unverified")
    qcsp = dict(zip(qc.k, qc.species))
    for i in d.index:
        psp = qcsp.get(d.at[i, "k"])
        if psp is None:
            continue
        if psp in D46_RENAME and d.at[i, "species"] == MUC:
            new = D46_RENAME[psp]
            change(i, "species", d.at[i, "species"], new, "popgen genotype/metadata species differs from the screen label", "D-46", f"popgen strain_qc.tsv species = {psp}")
            d.at[i, "species"], d.at[i, "species_source"] = new, "popgen genotype call (D-46)"
        elif psp in D46_MARKER and d.at[i, "species"] == MUC:
            d.at[i, "clade_marker"] = D46_MARKER[psp]
            change(i, "clade_marker", "", D46_MARKER[psp], "popgen genotypes this strain as R. aff. mucilaginosa; species stays R. mucilaginosa for now (popgen decision)", "D-46", f"popgen strain_qc.tsv species = {psp}")
    d["species_changed"] = d.species != d.species_screen

    # ---- 5. ploidy (D-47, D-54) -----------------------------------------------
    pl = dict(zip(ploidy.k, ploidy.ploidy))
    hs = dict(zip(hyb.k, hyb.subgroup))
    status, psrc, hsub = [], [], []
    for i in d.index:
        k, s = d.at[i, "k"], d.at[i, "sample_name"]
        if s == REFERENCE:
            status.append("haploid"); psrc.append("reference assembly strain (D-54)"); hsub.append(""); dec[i].add("D-54")
            notes[i].append("D-54: reference assembly strain; no callset sample")
        elif k in pl:
            if pl[k] == "haploid":
                status.append("haploid"); psrc.append("popgen ploidy_overrides.csv"); hsub.append("")
            elif k in hs:
                status.append("diploid_hybrid"); psrc.append("popgen ploidy_overrides.csv + hybrid_diploids.tsv"); hsub.append(hs[k])
            else:
                status.append("diploid_other"); psrc.append("popgen ploidy_overrides.csv"); hsub.append("")
        else:
            status.append("unknown"); psrc.append("not in popgen callset"); hsub.append("")
        dec[i].add("D-47")
        change(i, "ploidy_status", "", status[-1], "ploidy coding", "D-47", psrc[-1])
    d["ploidy_status"], d["ploidy_source"], d["hybrid_subgroup"] = status, psrc, hsub

    # ---- 6. identity flags (D-48) ----------------------------------------------
    agg = iss.groupby("k").agg(ids=("id", lambda x: ";".join(x)), cats=("category", lambda x: ";".join(sorted(set(x)))))
    d["identity_flag_ids"] = d.k.map(agg.ids).fillna("")
    d["identity_flag_categories"] = d.k.map(agg.cats).fillna("")
    dup = d.k.duplicated(keep=False)
    for i in d.index[dup]:
        others = ",".join(d.loc[(d.k == d.at[i, "k"]) & (d.index != i), "strain_id"])
        d.at[i, "identity_flag_ids"] = (d.at[i, "identity_flag_ids"] + ";" if d.at[i, "identity_flag_ids"] else "") + "build"
        d.at[i, "identity_flag_categories"] = (d.at[i, "identity_flag_categories"] + ";" if d.at[i, "identity_flag_categories"] else "") + "duplicate_sample_name"
        notes[i].append(f"duplicate sample_name with strain_id {others}; no genome; unresolved")
    for i in d.index:
        f = iss[iss.k == d.at[i, "k"]]
        for _, r in f.iterrows():
            notes[i].append(f"popgen {r['id']} ({r.category}): {r.finding}")

    # ---- 7. GWAS eligibility (D-47, D-48, D-54) ---------------------------------
    reasons = []
    for i in d.index:
        r = []
        if d.at[i, "species"] != MUC: r.append(f"species is {d.at[i, 'species']}")
        if d.at[i, "clade_marker"]: r.append(d.at[i, "clade_marker"])
        if d.at[i, "ploidy_status"] != "haploid": r.append(f"ploidy_status is {d.at[i, 'ploidy_status']}")
        if d.at[i, "sample_name"] == REFERENCE: r.append("reference strain; no genotypes")
        elif not d.at[i, "in_rmuc_core"]: r.append("not in popgen rmuc_core")
        if d.at[i, "identity_flag_categories"]: r.append("identity flag: " + d.at[i, "identity_flag_categories"])
        reasons.append("; ".join(r))
    d["gwas_exclude_reason"] = reasons
    elig = d.index[d.gwas_exclude_reason == ""]
    log(f"GWAS-eligible strains before de-cloning: {len(elig)}")

    # ---- 8. de-clone groups at several cutoffs (D-49) ---------------------------
    pw = pd.read_csv(POP / "variant_qc/declone/rmuc_core.pairwise.tsv.gz", sep="\t")
    pw["a"], pw["b"] = pw.strain_a.map(norm), pw.strain_b.map(norm)
    miss = pd.read_csv(POP / "variant_qc/rmuc_core.missing_by_strain.tsv", sep="\t")
    missf = dict(zip(miss.strain.map(norm), miss.missing_frac))
    cov = pd.read_csv(POP / "variant_qc/coverage.tsv", sep="\t")
    depth = dict(zip(cov.strain.map(norm), cov.mean_depth))
    keys = set(d.loc[elig, "k"])
    key2i = {d.at[i, "k"]: i for i in elig}
    for n in CUTOFFS:
        par = {k: k for k in keys}

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        sub = pw[pw.a.isin(keys) & pw.b.isin(keys) & (pw.diffs <= n)]
        for a, b in zip(sub.a, sub.b):
            par[find(a)] = find(b)
        comp = {}
        for k in sorted(keys, key=lambda x: d.at[key2i[x], "sample_name"]):
            comp.setdefault(find(k), []).append(k)
        gid, rep = {}, {}
        for gi, (_, members) in enumerate(sorted(comp.items(), key=lambda kv: d.at[key2i[kv[1][0]], "sample_name"]), start=1):
            best = sorted(members, key=lambda x: (missf.get(x, 1.0), -depth.get(x, 0.0), d.at[key2i[x], "sample_name"]))[0]
            for m in members:
                gid[m] = f"G{n}_{gi:03d}"
                rep[m] = (m == best)
        d[f"declone_group_le{n}"] = [gid.get(d.at[i, "k"], "") for i in d.index]
        d[f"declone_rep_le{n}"] = [bool(rep.get(d.at[i, "k"], False)) for i in d.index]
        log(f"  cutoff <= {n} SNPs: {len(comp)} groups = {sum(rep.values())} representatives")
    d["gwas_panel"] = d[f"declone_rep_le{PANEL_CUTOFF}"]
    log(f"GWAS panel (cutoff {PANEL_CUTOFF}): {int(d.gwas_panel.sum())} strains")
    assert int(d.gwas_panel.sum()) == EXPECTED_PANEL, f"D-49 expects {EXPECTED_PANEL} panel strains, got {int(d.gwas_panel.sum())}"
    # non-representatives of eligible groups get an explicit reason so the column explains every exclusion
    for i in elig:
        if not d.at[i, "gwas_panel"]:
            d.at[i, "gwas_exclude_reason"] = f"not the representative of its <= {PANEL_CUTOFF}-SNP group {d.at[i, f'declone_group_le{PANEL_CUTOFF}']}"

    # ---- 9. population columns (D-51: empty until the clusters are built) --------
    d["population_cluster"], d["population_method"], d["population_k"] = "", "", ""

    # ---- 10. notes and decision ids ----------------------------------------------
    d["curation_notes"] = ["; ".join(notes[i]) for i in d.index]
    d["decision_ids"] = [";".join(sorted(dec[i])) for i in d.index]
    d["curated_on"] = TODAY

    cols = (["strain_id", "sample_name", "popgen_strain", "species_screen", "species", "species_source", "species_changed", "clade_marker",
             "ploidy_status", "ploidy_source", "hybrid_subgroup",
             "in_callset", "in_rmuc_core", "in_rmuc_core_outgroup", "in_rmuc_with_hybrids", "in_hybrid_diploids", "clone_group", "near_identical_group",
             "identity_flag_ids", "identity_flag_categories", "gwas_exclude_reason", "gwas_panel"]
            + [f"declone_group_le{n}" for n in CUTOFFS] + [f"declone_rep_le{n}" for n in CUTOFFS]
            + ["population_cluster", "population_method", "population_k", "curation_notes", "decision_ids", "curated_on"])
    d[cols].to_csv(OUT / "strain_curation.csv", index=False)
    ch = pd.DataFrame(changes)
    ch = ch[~((ch.field == "ploidy_status") & (ch.new_value.isin(["haploid", "unknown"])) & True)]   # keep ploidy rows only where they carry information (diploid / reference)
    ch.to_csv(OUT / "strain_curation_changes.csv", index=False)

    # ---- summary -----------------------------------------------------------------
    log("\nspecies after curation:\n" + d.species.value_counts().to_string())
    log("\nploidy_status:\n" + d.ploidy_status.value_counts().to_string())
    log("\nspecies changes: " + str(len(ch[ch.field == "species"])) + "; marker changes: " + str(len(ch[ch.field == "clade_marker"])))
    log(f"panel: {int(d.gwas_panel.sum())}; wrote {OUT / 'strain_curation.csv'} ({len(d)} rows) and strain_curation_changes.csv ({len(ch)} rows)")


if __name__ == "__main__":
    sys.exit(main())
