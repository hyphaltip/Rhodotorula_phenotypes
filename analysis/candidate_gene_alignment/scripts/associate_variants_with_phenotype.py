#!/usr/bin/env python3
"""Associate candidate-gene coding variants with color/pigment phenotypes across the
panel, using the SAME population-aware battery established by
analysis/gwas/scripts/check_population_vs_locus.py (Test A within-population re-test
with near-clone collapse, Test B inverse-variance meta-analysis, Test C
population-covariate regression).

Why population-aware tests (not a naive pooled regression):
- The panel is near-clonal with real population structure (6 populations); a naive
  genome-wide-style genotype~phenotype regression on all strains is confounded by
  between-population allele-frequency x phenotype differences (this exact failure was
  caught by the whole D-17/D-18 disambiguation exercise, and scaffold_7 -- where both
  pilot genes live -- was the clearest example: lab_L's intergenic scaffold_7:172154
  hit collapsed to a population-4 near-fixation artifact).
- These genes are NOT genome-wide-scan hits (no Tier A p-value to anchor a causal
  claim), so the honest test is: does a coding variant predict the phenotype WITHIN
  populations, replicated across >=2 populations, over and above population structure?

Inputs: per-gene variant_table.csv (gt__ columns), strain context, phenotype CSV,
population map, culled-182 strain list, existing per-population KING collapse maps
(from the population_vs_locus run).

Outputs: assoc_summary.csv (one row per variant x trait test) + per-test details.

This is the phenotype-association step the pilot design explicitly deferred
(CANDIDATE_GENE_ALIGNMENT.md / TODO item "splice_donor genotype vs color phenotypes"):
here it is for ALL coding-impact (HIGH/MODERATE) variants in a gene, not just one.
"""
import argparse
import csv
import glob
import os
import subprocess

import numpy as np
import pandas as pd
from scipy import stats

INDEL_VCF = "data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.INDEL.combined_selected.vcf.gz"


def query_indel_genotypes(vcf: str, scaffold: str, start: int, end: int,
                          accepted: set[str], bcftools: str = "bcftools") -> dict[tuple, dict[str, str]]:
    """Query the INDEL VCF for a CDS region; return {(pos, ref, alt): {strain: '0'|'1'}}.
    Restricts to accepted-panel strains and binarizes to ref/non-ref (haploid panel).
    Skips multiallelic-carrying records (rare) -- they belong in a second pass."""
    r = subprocess.run([bcftools, "query", "-r", f"{scaffold}:{start}-{end}",
                        "-f", "%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n", vcf],
                       capture_output=True, text=True, check=False)
    out = {}
    for line in r.stdout.splitlines():
        parts = line.split("\t")
        pos, ref, alt = parts[0], parts[1], parts[2]
        genos = {}
        any_accepted = False
        for tok in parts[3:]:
            if "=" not in tok:
                continue
            s, gt = tok.split("=")
            if s not in accepted:
                continue
            gt = gt.replace("|", "/").split("/")[0]
            genos[s] = "1" if gt not in ("0", ".") else ("0" if gt == "0" else "missing")
            any_accepted = True
        if any_accepted:
            out[(int(pos), ref, alt)] = genos
    return out


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_variant_table(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    gt_cols = [c for c in df.columns if c.startswith("gt__")]
    df["strain_geno"] = df.apply(
        lambda r: {c[4:]: r[c] for c in gt_cols}, axis=1)
    return df


def load_pop_map(pop_csv: str) -> dict[str, str]:
    pop = {}
    with open(pop_csv) as f:
        for r in csv.DictReader(f):
            pop[r["Strain"]] = r["Pop"]
    return pop


def load_culled_set(culled_txt: str) -> set[str]:
    return set(l.strip() for l in open(culled_txt) if l.strip())


def load_collapse_maps(work_dir: str, threshold: float = 0.005) -> dict[str, dict[str, str]]:
    """Rebuild {population: {strain: representative}} from the precomputed KING
    .kin0 tables of the population_vs_locus run (same greedy IBS0-collapse logic as
    check_population_vs_locus.near_clone_collapse_within_pop, without re-running plink2)."""
    rep_maps: dict[str, dict[str, str]] = {}
    for kin0 in sorted(glob.glob(os.path.join(work_dir, "collapse_pop*_king.kin0"))):
        pop = os.path.basename(kin0).split("_")[1][3:]
        pairs = pd.read_csv(kin0, sep="\t")
        pairs.columns = [c.lstrip("#") for c in pairs.columns]
        if "IBS0" not in pairs.columns:
            continue
        pairs = pairs.sort_values("IBS0").reset_index(drop=True)
        alive = sorted(set(pairs["IID1"]) | set(pairs["IID2"]))
        rep = {s: s for s in alive}
        while True:
            cand = pairs[pairs["IID1"].isin(alive) & pairs["IID2"].isin(alive)]
            if cand.empty:
                break
            row = cand.loc[cand["IBS0"].idxmin()]
            if row["IBS0"] >= threshold:
                break
            a, b = row["IID1"], row["IID2"]
            drop, keep = (b, a) if b > a else (a, b)
            for s in list(rep):
                if rep[s] == drop:
                    rep[s] = keep
            alive.discard(drop)
        rep_maps[pop] = rep
    return rep_maps


# ---------------------------------------------------------------------------
# Test A + B: within-population re-test + meta-analysis (identical to the GWAS
# disambiguation battery)
# ---------------------------------------------------------------------------

def within_pop_test(dosage: dict[str, float], pheno: dict[str, float], pop: dict[str, str],
                     testable: set[str], rep_map: dict[str, str], min_n_per_class: int = 3):
    rows = []
    for p in sorted(set(pop.get(s) for s in testable if pop.get(s))):
        members = [s for s in testable if pop.get(s) == p]
        n_raw = len(members)
        if n_raw < 2 * min_n_per_class:
            rows.append(dict(population=p, n_raw=n_raw, n_effective=None, beta=np.nan,
                              se=np.nan, p=np.nan, n0=0, n1=0, tested=False,
                              reason="n_raw too small"))
            continue
        reps = rep_map or {s: s for s in members}
        collapsed = {}
        for s in members:
            collapsed.setdefault(reps.get(s, s), []).append(s)
        eff_dosage, eff_pheno = [], []
        for r, group in collapsed.items():
            eff_dosage.append(np.mean([dosage[s] for s in group]))
            eff_pheno.append(np.mean([pheno[s] for s in group]))
        eff_dosage = np.array(eff_dosage)
        eff_pheno = np.array(eff_pheno)
        n_eff = len(eff_dosage)
        n0 = int(np.sum(eff_dosage < 0.5))
        n1 = int(np.sum(eff_dosage >= 0.5))
        if n0 < min_n_per_class or n1 < min_n_per_class:
            rows.append(dict(population=p, n_raw=n_raw, n_effective=n_eff, beta=np.nan,
                              se=np.nan, p=np.nan, n0=n0, n1=n1, tested=False,
                              reason=f"allele classes too small ({n0}/{n1})"))
            continue
        lr = stats.linregress(eff_dosage, eff_pheno)
        rows.append(dict(population=p, n_raw=n_raw, n_effective=n_eff, beta=lr.slope,
                          se=lr.stderr, p=lr.pvalue, n0=n0, n1=n1, tested=True, reason=""))
    return pd.DataFrame(rows)


def meta_analysis(within_pop_df: pd.DataFrame):
    valid = within_pop_df[within_pop_df["tested"] & within_pop_df["se"].notna() & (within_pop_df["se"] > 0)]
    if len(valid) < 2:
        return dict(meta_beta=np.nan, meta_se=np.nan, meta_p=np.nan, n_pops=len(valid),
                     cochran_q=np.nan, i_squared=np.nan, n_directionally_consistent=0)
    w = 1.0 / valid["se"].to_numpy() ** 2
    beta = valid["beta"].to_numpy()
    meta_beta = np.sum(w * beta) / np.sum(w)
    meta_se = np.sqrt(1.0 / np.sum(w))
    meta_z = meta_beta / meta_se
    meta_p = 2 * stats.norm.sf(abs(meta_z))
    q = np.sum(w * (beta - meta_beta) ** 2)
    df = len(valid) - 1
    i_sq = max(0.0, (q - df) / q) if q > 0 else 0.0
    n_consistent = int(np.sum(np.sign(beta) == np.sign(meta_beta))) if meta_beta == meta_beta else 0
    return dict(meta_beta=meta_beta, meta_se=meta_se, meta_p=meta_p, n_pops=len(valid),
                cochran_q=q, i_squared=i_sq, n_directionally_consistent=n_consistent)


def population_covariate_test(dosage: dict[str, float], pheno: dict[str, float], pop: dict[str, str],
                               strains: list[str]):
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    df = pd.DataFrame([
        dict(strain=s, genotype=dosage[s], phenotype=pheno[s], population=pop[s])
        for s in strains if s in dosage and s in pheno and pop.get(s)
    ])
    if df["population"].nunique() < 2 or len(df) < 10:
        return dict(partial_r2=np.nan, partial_p=np.nan, n=len(df), n_pops=df["population"].nunique())
    full = smf.ols("phenotype ~ genotype + C(population)", data=df).fit()
    reduced = smf.ols("phenotype ~ C(population)", data=df).fit()
    sse_full = float(np.sum(full.resid ** 2))
    sse_reduced = float(np.sum(reduced.resid ** 2))
    partial_r2 = (sse_reduced - sse_full) / sse_reduced if sse_reduced > 0 else np.nan
    try:
        partial_p = float(full.pvalues.get("genotype", np.nan))
    except Exception:
        partial_p = np.nan
    return dict(partial_r2=partial_r2, partial_p=partial_p, n=len(df), n_pops=df["population"].nunique())


def veridict(meta, cov, min_pops: int = 2) -> str:
    """Apply the check_population_vs_locus verdict logic (test D excluded -- no
    genome-wide null pool exists for these un-scanned genes, and where D was used it
    was uninformative anyway)."""
    meta_p = meta.get("meta_p", np.nan)
    n_pops = meta.get("n_pops", 0)
    n_consistent = meta.get("n_directionally_consistent", 0)
    directionally_consistent = n_pops >= min_pops and n_consistent >= min_pops
    partial_r2 = cov.get("partial_r2", np.nan)
    r2_negligible = (partial_r2 == partial_r2) and partial_r2 < 0.01
    meta_nominally_sig = (meta_p == meta_p) and meta_p < 0.05
    if not directionally_consistent or r2_negligible or not meta_nominally_sig:
        return "likely_population_artifact" if (not directionally_consistent or r2_negligible) else "ambiguous_underpowered"
    corroborated = (partial_r2 == partial_r2) and partial_r2 >= 0.02
    if directionally_consistent and meta_nominally_sig and corroborated:
        return "likely_real"
    return "ambiguous_underpowered"


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant-tables", nargs="+", required=True,
                     help="per-gene variant_table.csv files")
    ap.add_argument("--pop-csv", required=True)
    ap.add_argument("--culled-list", required=True, help="182-strain culled panel list")
    ap.add_argument("--collapse-work-dir", required=True,
                     help="dir with collapse_pop*_king.kin0 from the population_vs_locus run")
    ap.add_argument("--traits", nargs="+", default=["chroma", "sat", "bright", "lab_L", "lab_a", "lab_b"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pop_all = load_pop_map(args.pop_csv)
    culled = load_culled_set(args.culled_list)
    rep_maps = load_collapse_maps(args.collapse_work_dir)
    print(f"Loaded {len(rep_maps)} population collapse maps")

    pheno = pd.read_csv("analysis/gwas/results/gwas_next_phenotypes.csv").set_index("strain_code")

    results = []
    details = []
    for vt in args.variant_tables:
        gene = os.path.basename(os.path.dirname(vt))
        df = load_variant_table(vt)
        coding = df[df["impact"].isin(["HIGH", "MODERATE"])]
        print(f"\n=== {gene}: {len(coding)} coding-impact variants to test ===")
        for _, row in coding.iterrows():
            g = row["strain_geno"]
            variant_label = (f"{gene}:{row['pos']}"
                             f"({row['hgvs_c']}{',' + str(row['hgvs_p']) if pd.notna(row['hgvs_p']) else ''})"
                             f" [{row['consequence']}]")
            # restrict to culled-182 strains with usable genotype + phenotype
            strain_info = {}
            for s, gt in g.items():
                if s not in culled:
                    continue
                if gt not in ("ref", "alt"):
                    continue
                if s not in pop_all:
                    continue
                strain_info[s] = float(0 if gt == "ref" else 1)
            for trait in args.traits:
                phen = {s: float(v) for s, v in pheno[trait].items()
                        if pd.notna(v) and s in strain_info}
                testable = set(strain_info) & set(phen)
                dosage = {s: strain_info[s] for s in testable}
                n_alt = sum(1 for v in dosage.values() if v == 1)
                n_eff_distinct = len({dosage[s] for s in testable})

                a_df = within_pop_test(dosage, phen, pop_all, testable,
                                       {p: rep_maps.get(p, {s: s for s in testable})
                                        for p in sorted(set(pop_all.get(s) for s in testable))})
                meta = meta_analysis(a_df)
                cov = population_covariate_test(dosage, phen, pop_all, list(testable))
                verdict = veridict(meta, cov)

                results.append(dict(
                    gene=gene, pos=row["pos"], variant=variant_label,
                    consequence=row["consequence"], hgvs_c=row["hgvs_c"], hgvs_p=row["hgvs_p"] if pd.notna(row["hgvs_p"]) else "",
                    trait=trait, n_strains=len(testable), n_alt=n_alt,
                    meta_beta=meta["meta_beta"], meta_se=meta["meta_se"], meta_p=meta["meta_p"],
                    n_pops_tested=meta["n_pops"], n_directionally_consistent=meta["n_directionally_consistent"],
                    cochran_q=meta["cochran_q"], i_squared=meta["i_squared"],
                    partial_r2=cov["partial_r2"], partial_p=cov["partial_p"],
                    cov_n=cov["n"], cov_n_pops=cov["n_pops"], verdict=verdict,
                ))
                for _, pr in a_df.iterrows():
                    details.append(dict(gene=gene, pos=row["pos"], variant=variant_label, trait=trait,
                                        population=pr.population, n_raw=pr.n_raw, n_effective=pr.n_effective,
                                        n0=pr.n0, n1=pr.n1, beta=pr.beta, se=pr.se, p=pr.p,
                                        tested=pr.tested, reason=pr.reason))

    out = pd.DataFrame(results)
    # BH-FDR across the full variant x trait test matrix, using meta_p (the one clean
    # null in the battery, per D-17 design).
    valid_meta = out["meta_p"].dropna()
    if len(valid_meta) > 0:
        m = len(valid_meta)
        order = valid_meta.sort_values().index
        p_sorted = valid_meta.loc[order].to_numpy()
        ranks = np.arange(1, m + 1)
        fdr = p_sorted * m / ranks
        fdr = np.minimum.accumulate(fdr[::-1])[::-1]
        out["meta_p_fdr"] = pd.Series(fdr, index=order).reindex(out.index).astype(float)
    else:
        out["meta_p_fdr"] = np.nan
    out["fdr_sig"] = out["meta_p_fdr"] < 0.05

    # Sort: FDR-significant first, then by p
    out = out.sort_values(["fdr_sig", "meta_p"], ascending=[False, True]).reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, index=False)
    details_out = args.out.replace(".csv", "_within_pop_details.csv")
    pd.DataFrame(details).to_csv(details_out, index=False)

    print(f"\nWrote {args.out} ({len(out)} tests) and {details_out}")
    sig = out[out["fdr_sig"]]
    print(f"\nFDR-significant (q<0.05): {len(sig)} / {len(out)} tests")
    if not sig.empty:
        print(sig[["gene", "pos", "hgvs_c", "hgvs_p", "trait", "n_alt", "meta_p", "meta_p_fdr",
                   "partial_r2", "n_pops_tested", "n_directionally_consistent", "verdict"]].to_string(index=False))
    print("\nTop 15 by raw meta_p (any verdict):")
    print(out.head(15)[["gene", "pos", "hgvs_c", "hgvs_p", "trait", "n_alt", "meta_p", "meta_p_fdr",
                        "partial_r2", "n_pops_tested", "n_directionally_consistent", "verdict"]].to_string(index=False))
    print("\nVerdict counts:", out["verdict"].value_counts().to_dict())


if __name__ == "__main__":
    main()
