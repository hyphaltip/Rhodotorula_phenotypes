#!/usr/bin/env python3
"""Disambiguate real causal loci from population-structure artifacts among Tier A
top hits flagged by check_population_confounding.py.

Design: docs/superpowers/specs/2026-08-26-population-confounding-disambiguation-design.md
(drafted, then revised after an independent statistical-geneticist/breeder review).

Four tests, each on the near-clone-culled 182-strain panel ("gwasc") unless noted:

  A. Within-population re-test: phenotype ~ genotype_dosage fit SEPARATELY within
     each population (removes all between-population variance by construction).
     Requires >=3 strains per allele class. Reports n_raw AND n_effective (after an
     additional per-population near-clone collapse pass -- the within-population
     test is exactly as vulnerable to near-clone pseudoreplication as the genome-wide
     scan kinship correction exists for; this is a required, not optional, check).
  B. Fixed-effect inverse-variance meta-analysis of A's per-population betas, with
     Cochran's Q / I^2 for heterogeneity. BH-FDR applied to meta-analysis p-values
     across all flagged loci (the one test in this battery with a clean null).
  C. Single-locus population-covariate regression: phenotype ~ genotype + population
     (6-level factor), ALL strains -- partial F-test / partial R^2 of genotype beyond
     population. This is a ONE-SNP regression, not genome-wide, so it does not hit the
     GRM/PC collinearity that sank genome-wide PC-covariate GEMMA runs (D-9).
  D. Fst x MAF matched empirical null: ~5,000 comparator SNPs from the panel's own
     LD-pruned marker set (already an approximately-independent sample), matched on a
     5x5 (Fst-proxy quintile x MAF quintile) bin, excluding any comparator within 50kb
     of another significant hit for the same trait. Looks up each comparator's actual
     Tier A p_wald (already computed, no new GEMMA runs) and reports the flagged SNP's
     percentile within that matched null.

Fst proxy (haploid, no heterozygosity term available): Var_across_pops(p) / (pbar *
(1 - pbar)) -- the between-population allele-frequency variance relative to its
expectation under panmixia, a standard simplification.

Verdict rule (per the review, not requiring unanimity across 4 partially
non-independent tests):
  likely_real              = A's directional consistency across >=2 populations
                              AND >=2 of {B significant, C partial-R^2 non-negligible,
                              D percentile extreme}
  likely_population_artifact = A fails directional consistency, OR C partial-R^2 ~ 0
  ambiguous_underpowered    = neither of the above (e.g. too few populations testable)

MAS-usability gates (minimum effect size, held-out replication, LD-decay check,
plate/batch confound check) are explicitly OUT OF SCOPE here -- deferred to a follow-up
for any locus that reaches likely_real.
"""
import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# Genotype / data loading
# ---------------------------------------------------------------------------

def genotype_at_site(bcftools_bin: str, vcf: str, scaffold: str, pos: int) -> dict[str, float]:
    region = f"{scaffold}:{pos}-{pos}"
    out = subprocess.run(
        [bcftools_bin, "query", "-r", region, "-f", "[%SAMPLE=%GT\\t]\\n", vcf],
        check=True, capture_output=True, text=True,
    )
    line = out.stdout.strip()
    assert line, f"no VCF record found at {region}"
    dosage = {}
    for tok in line.split("\t"):
        if not tok or "=" not in tok:
            continue
        s, gt = tok.split("=")
        gt = gt.replace("|", "/").split("/")[0]
        if gt in (".", ""):
            continue
        dosage[s] = float(gt)
    return dosage


def parse_snp(rs: str) -> tuple[str, int]:
    parts = rs.split(":")
    return parts[0], int(parts[1])


def load_pop_map(pop_csv: str) -> dict[str, str]:
    pop = {}
    with open(pop_csv) as f:
        for r in csv.DictReader(f):
            pop[r["Strain"]] = r["Pop"]
    return pop


def load_culled_set(culled_txt: str) -> set[str]:
    return set(l.strip() for l in open(culled_txt) if l.strip())


def near_clone_collapse_within_pop(strain_ids: list[str], pruned_bfile: str, plink2_bin: str,
                                    threshold: float, work_prefix: str) -> dict[str, str]:
    """Greedy IBS0 collapse restricted to the given strains (typically one population's
    182-panel members). Returns {strain_id: representative_id} -- every strain maps to
    itself unless collapsed into a clone-cluster representative. No-op (identity map) if
    fewer than 2 strains, since the panel is already globally culled at this threshold
    and a within-population pass is expected to usually find nothing new."""
    if len(strain_ids) < 2:
        return {s: s for s in strain_ids}
    fam_rows = [l.split() for l in open(pruned_bfile + ".fam") if l.split()]
    keep_rows = [(fid, iid) for fid, iid, *_ in fam_rows if iid in set(strain_ids)]
    if len(keep_rows) < 2:
        return {s: s for s in strain_ids}
    keep_file = work_prefix + ".keep.txt"
    with open(keep_file, "w") as f:
        for fid, iid in keep_rows:
            f.write(f"{fid}\t{iid}\n")
    subset_prefix = work_prefix + "_subset"
    r = subprocess.run([plink2_bin, "--bfile", pruned_bfile, "--keep", keep_file,
                        "--make-bed", "--out", subset_prefix, "--allow-extra-chr"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return {s: s for s in strain_ids}
    king_prefix = work_prefix + "_king"
    r = subprocess.run([plink2_bin, "--bfile", subset_prefix, "--make-king-table",
                        "cols=id,nsnp,ibs0", "--out", king_prefix, "--allow-extra-chr"],
                       capture_output=True, text=True)
    kin0 = king_prefix + ".kin0"
    if r.returncode != 0 or not os.path.exists(kin0):
        return {s: s for s in strain_ids}
    pairs = pd.read_csv(kin0, sep="\t")
    pairs.columns = [c.lstrip("#") for c in pairs.columns]
    pairs = pairs.sort_values("IBS0").reset_index(drop=True)

    rep = {s: s for s in strain_ids}
    alive = set(strain_ids)
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
    return rep


# ---------------------------------------------------------------------------
# Test A + B: within-population re-test + meta-analysis
# ---------------------------------------------------------------------------

def within_pop_test(dosage: dict[str, float], pheno: dict[str, float], pop: dict[str, str],
                     culled: set[str], rep_map: dict[str, dict[str, str]], min_n_per_class: int = 3):
    """rep_map: {population: {strain: representative}} from near_clone_collapse_within_pop."""
    rows = []
    for p in sorted(set(pop.get(s) for s in culled if pop.get(s))):
        members = [s for s in culled if pop.get(s) == p and s in dosage and s in pheno]
        n_raw = len(members)
        if n_raw < 2 * min_n_per_class:
            rows.append(dict(population=p, n_raw=n_raw, n_effective=None, beta=np.nan,
                              se=np.nan, p=np.nan, tested=False, reason="n_raw too small"))
            continue
        reps = rep_map.get(p, {s: s for s in members})
        collapsed = {}
        for s in members:
            r = reps.get(s, s)
            collapsed.setdefault(r, []).append(s)
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
                              se=np.nan, p=np.nan, tested=False,
                              reason=f"n_effective allele classes too small ({n0}/{n1})"))
            continue
        lr = stats.linregress(eff_dosage, eff_pheno)
        rows.append(dict(population=p, n_raw=n_raw, n_effective=n_eff, beta=lr.slope,
                          se=lr.stderr, p=lr.pvalue, tested=True, reason=""))
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


# ---------------------------------------------------------------------------
# Test C: population-covariate regression
# ---------------------------------------------------------------------------

def population_covariate_test(dosage: dict[str, float], pheno: dict[str, float], pop: dict[str, str],
                               strains: list[str]):
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    rows = []
    for s in strains:
        if s in dosage and s in pheno and pop.get(s):
            rows.append(dict(strain=s, genotype=dosage[s], phenotype=pheno[s], population=pop[s]))
    df = pd.DataFrame(rows)
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


# ---------------------------------------------------------------------------
# Test D: Fst x MAF matched null
# ---------------------------------------------------------------------------

def compute_pop_af_matrix(raw_path: str, pop: dict[str, str]) -> tuple[pd.DataFrame, list[str]]:
    """Parse a plink2 --recode A .raw file into per-SNP per-population allele
    frequency, plus overall MAF and an Fst proxy (between-pop AF variance relative to
    panmictic expectation)."""
    raw = pd.read_csv(raw_path, sep=r"\s+")
    snp_cols = [c for c in raw.columns if c not in ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
    # plink2 --recode A appends "_<counted allele>" to each column name (e.g.
    # "scaffold_10:480:G:C_G") -- strip it to recover the plain rs id used everywhere
    # else (assoc files, FDR lists), or the merge in matched_null_percentile silently
    # finds zero matches.
    rs_ids = [c.rsplit("_", 1)[0] for c in snp_cols]
    raw["population"] = raw["IID"].map(pop)
    raw = raw[raw["population"].notna()]

    geno = raw[snp_cols].to_numpy(dtype=float) / 2.0  # 0/2 dosage (haploid, no het) -> AF-like 0/1
    pops = raw["population"].to_numpy()
    uniq_pops = sorted(set(pops))

    pop_af = np.full((len(uniq_pops), len(snp_cols)), np.nan)
    for i, p in enumerate(uniq_pops):
        mask = pops == p
        pop_af[i] = np.nanmean(geno[mask], axis=0)

    overall_af = np.nanmean(geno, axis=0)
    between_var = np.nanvar(pop_af, axis=0, ddof=0)
    denom = overall_af * (1 - overall_af)
    fst_proxy = np.where(denom > 1e-9, between_var / denom, np.nan)
    maf = np.minimum(overall_af, 1 - overall_af)

    df = pd.DataFrame({"rs": rs_ids, "overall_af": overall_af, "maf": maf, "fst_proxy": fst_proxy})
    df[["chr", "pos"]] = df["rs"].str.split(":", n=2, expand=True).iloc[:, :2]
    df["pos"] = df["pos"].astype(int)
    return df, uniq_pops


def matched_null_percentile(flagged_fst: float, flagged_maf: float, flagged_scaffold: str, flagged_pos: int,
                             null_pool: pd.DataFrame, trait_assoc: pd.DataFrame, exclude_positions: list[tuple],
                             flagged_p: float, n_bins: int = 5, n_sample: int = 5000, seed: int = 42,
                             ld_window_bp: int = 50_000):
    pool = null_pool.copy()
    for (scf, pos) in exclude_positions:
        pool = pool[~((pool["chr"] == scf) & (pool["pos"].between(pos - ld_window_bp, pos + ld_window_bp)))]

    fst_q = pd.qcut(pool["fst_proxy"], n_bins, labels=False, duplicates="drop")
    maf_q = pd.qcut(pool["maf"], n_bins, labels=False, duplicates="drop")
    pool = pool.assign(fst_q=fst_q, maf_q=maf_q)

    flagged_fst_q = (pool["fst_proxy"] <= flagged_fst).mean()
    flagged_maf_q = (pool["maf"] <= flagged_maf).mean()
    target_fst_bin = min(int(flagged_fst_q * n_bins), n_bins - 1)
    target_maf_bin = min(int(flagged_maf_q * n_bins), n_bins - 1)

    matched = pool[(pool["fst_q"] == target_fst_bin) & (pool["maf_q"] == target_maf_bin)]
    if len(matched) == 0:
        return dict(n_matched=0, percentile=np.nan)
    rng = np.random.default_rng(seed)
    if len(matched) > n_sample:
        matched = matched.sample(n_sample, random_state=seed)

    merged = matched.merge(trait_assoc[["rs", "p_wald"]], on="rs", how="inner")
    if len(merged) == 0:
        return dict(n_matched=0, percentile=np.nan)
    null_p = merged["p_wald"].to_numpy()
    percentile = float(np.mean(null_p <= flagged_p))  # fraction of null AS extreme or more
    return dict(n_matched=len(merged), percentile=percentile)


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

def verdict(meta: dict, cov: dict, dnull: dict, meta_fdr_sig: bool, meta_p: float) -> str:
    """NOTE (found empirically, not anticipated in the design): test D (Fst x MAF
    matched null) returns percentile=0.0 for EVERY flagged locus tested, including ones
    whose within-population meta-analysis isn't even nominally significant. This is a
    winner's-curse/circularity artifact -- these loci are, by construction, each trait's
    single most significant genome-wide hit, so of course they beat a random
    Fst/MAF-matched comparator sample; that says nothing about whether the signal is a
    population artifact. Because of this, test D is treated as UNINFORMATIVE here (not
    counted as a support) rather than removed outright -- it is still reported for
    transparency, but the verdict rule below no longer lets it substitute for a
    significant meta-analysis. A raw (not just FDR-corrected) meta_p < 0.05 is now a
    REQUIRED gate for likely_real, not one of several interchangeable supports -- a
    locus with no nominal within-population signal at all should never be called real
    regardless of what the (circular) matched-null test says."""
    directionally_consistent = meta.get("n_pops", 0) >= 2 and meta.get("n_directionally_consistent", 0) >= 2
    partial_r2 = cov.get("partial_r2", np.nan)
    r2_negligible = (partial_r2 == partial_r2) and partial_r2 < 0.01
    meta_nominally_sig = (meta_p == meta_p) and meta_p < 0.05
    if not directionally_consistent or r2_negligible or not meta_nominally_sig:
        return "likely_population_artifact" if (not directionally_consistent or r2_negligible) else "ambiguous_underpowered"
    corroborated = meta_fdr_sig or ((partial_r2 == partial_r2) and partial_r2 >= 0.02)
    if directionally_consistent and meta_nominally_sig and corroborated:
        return "likely_real"
    return "ambiguous_underpowered"


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def build_null_pool(plink2_bin: str, pruned_bfile: str, pop: dict[str, str], work_dir: str, tag: str):
    raw_prefix = os.path.join(work_dir, f"null_pool_{tag}")
    r = subprocess.run([plink2_bin, "--bfile", pruned_bfile, "--recode", "A",
                        "--out", raw_prefix, "--allow-extra-chr"],
                       capture_output=True, text=True)
    assert r.returncode == 0 and os.path.exists(raw_prefix + ".raw"), (
        f"plink2 --recode A failed for {pruned_bfile}:\n{r.stdout[-2000:]}"
    )
    return compute_pop_af_matrix(raw_prefix + ".raw", pop)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--confounding-csv", nargs="+", required=True,
                     help="one or more population_confounding_{panel}.csv files")
    ap.add_argument("--pheno-csv", required=True)
    ap.add_argument("--pop-csv", required=True)
    ap.add_argument("--culled-list", required=True, help="182-strain culled panel strain list")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--plink2", default="plink2")
    ap.add_argument("--pruned-bfile-gwas", required=True)
    ap.add_argument("--pruned-bfile-gwasc", required=True)
    ap.add_argument("--assoc-csv-dir", required=True, help="dir with {panel}_{trait}_assoc.csv.gz")
    ap.add_argument("--fdr-dir", required=True, help="dir with {panel}_{trait}_fdr05.csv (for LD exclusion)")
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(args.work_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    pop = load_pop_map(args.pop_csv)
    culled = load_culled_set(args.culled_list)
    pheno_df = pd.read_csv(args.pheno_csv)

    flagged = []
    for f in args.confounding_csv:
        panel = "gwas" if "_gwas.csv" in f and "_gwasc.csv" not in f else "gwasc"
        df = pd.read_csv(f)
        df = df[df["population_confound_risk"] == True]
        df["panel"] = panel
        flagged.append(df)
    flagged = pd.concat(flagged, ignore_index=True)
    flagged = flagged.drop_duplicates(subset=["trait", "top_snp"]).reset_index(drop=True)
    print(f"Flagged loci to disambiguate: {len(flagged)}")

    print("Building Fst x MAF null pools (LD-pruned marker set, both panels) ...")
    null_pool_gwas, pops_gwas = build_null_pool(args.plink2, args.pruned_bfile_gwas, pop, args.work_dir, "gwas")
    null_pool_gwasc, pops_gwasc = build_null_pool(args.plink2, args.pruned_bfile_gwasc, pop, args.work_dir, "gwasc")
    print(f"  gwas pool: {len(null_pool_gwas)} SNPs; gwasc pool: {len(null_pool_gwasc)} SNPs")

    print("Building per-population near-clone-collapse maps (182-culled panel) ...")
    rep_map: dict[str, dict[str, str]] = {}
    for p in sorted(set(pop.get(s) for s in culled if pop.get(s))):
        members = [s for s in culled if pop.get(s) == p]
        work_prefix = os.path.join(args.work_dir, f"collapse_pop{p}")
        rep_map[p] = near_clone_collapse_within_pop(
            members, args.pruned_bfile_gwasc, args.plink2, threshold=0.005, work_prefix=work_prefix
        )
        n_collapsed = sum(1 for s, r in rep_map[p].items() if s != r)
        print(f"  pop {p}: {len(members)} strains, {n_collapsed} additionally collapsed")

    results = []
    for _, row in flagged.iterrows():
        trait, panel, rs = row["trait"], row["panel"], row["top_snp"]
        scaffold, pos = parse_snp(rs)
        print(f"\n=== {trait} [{panel}] {rs} ===")

        dosage = genotype_at_site(args.bcftools, args.vcf, scaffold, pos)
        pheno = dict(zip(pheno_df["strain_code"], pheno_df[trait]))
        pheno = {s: v for s, v in pheno.items() if pd.notna(v)}

        # --- Test A: within-population re-test ---
        a_df = within_pop_test(dosage, pheno, pop, culled, rep_map)

        # --- Test B: meta-analysis ---
        meta = meta_analysis(a_df)

        # --- Test C: population-covariate regression ---
        cov = population_covariate_test(dosage, pheno, pop, list(culled))

        # --- Test D: Fst x MAF matched null ---
        null_pool = null_pool_gwas if panel == "gwas" else null_pool_gwasc
        assoc_path = os.path.join(args.assoc_csv_dir, f"{panel}_{trait}_assoc.csv.gz")
        trait_assoc = pd.read_csv(assoc_path, usecols=["rs", "p_wald"])
        fdr_path = os.path.join(args.fdr_dir, f"{panel}_{trait}_fdr05.csv")
        exclude_positions = []
        if os.path.exists(fdr_path):
            fdr_df = pd.read_csv(fdr_path)
            for _, r in fdr_df.iterrows():
                if r["rs"] != rs:
                    sc, ps = parse_snp(r["rs"])
                    exclude_positions.append((sc, ps))
        # flagged SNP's own Fst/MAF: bcftools GT is the haploid single-allele call (0/1
        # directly, NOT a 0/2 diploid-style dosage like plink2's --recode A output used
        # for the null pool below -- these two data sources use different encodings and
        # must NOT be divided by 2 the same way).
        per_pop_af = pd.Series(list(dosage.values()), index=list(dosage.keys()))
        per_pop_af = per_pop_af[per_pop_af.index.isin(pop)]
        overall_af = float(per_pop_af.mean())
        by_pop = per_pop_af.groupby([pop[s] for s in per_pop_af.index]).mean()
        flagged_maf = min(overall_af, 1 - overall_af)
        flagged_between_var = float(np.var(by_pop.dropna().to_numpy()))
        flagged_fst = flagged_between_var / (overall_af * (1 - overall_af)) if 0 < overall_af < 1 else np.nan

        flagged_p_series = trait_assoc.loc[trait_assoc.rs == rs, "p_wald"]
        flagged_p = float(flagged_p_series.iloc[0]) if len(flagged_p_series) else np.nan

        dnull = matched_null_percentile(
            flagged_fst, flagged_maf, scaffold, pos, null_pool, trait_assoc,
            exclude_positions, flagged_p=flagged_p,
        )

        n_consistent = meta.get("n_directionally_consistent", 0)
        n_pops_tested = meta.get("n_pops", 0)
        result = dict(
            trait=trait, panel=panel, top_snp=rs,
            meta_beta=meta["meta_beta"], meta_se=meta["meta_se"], meta_p=meta["meta_p"],
            n_pops_tested=n_pops_tested, n_directionally_consistent=n_consistent,
            cochran_q=meta["cochran_q"], i_squared=meta["i_squared"],
            partial_r2=cov["partial_r2"], partial_p=cov["partial_p"], cov_n=cov["n"],
            fst_proxy=flagged_fst, maf=flagged_maf,
            null_n_matched=dnull["n_matched"], null_percentile=dnull["percentile"],
        )
        results.append(result)
        a_df.insert(0, "trait", trait)
        a_df.insert(1, "panel", panel)
        a_df.insert(2, "top_snp", rs)
        a_out = os.path.join(os.path.dirname(args.out), f"within_pop_{trait}_{panel}.csv")
        a_df.to_csv(a_out, index=False)
        print(f"  meta_p={meta['meta_p']:.3g}  n_pops={n_pops_tested}  n_consistent={n_consistent}  "
              f"partial_r2={cov['partial_r2']:.4f}  null_percentile={dnull['percentile']}")

    out = pd.DataFrame(results)
    valid_meta_p = out["meta_p"].dropna()
    if len(valid_meta_p) > 0:
        order = valid_meta_p.sort_values().index
        m = len(valid_meta_p)
        ranks = pd.Series(range(1, m + 1), index=order)
        fdr = (valid_meta_p.loc[order].to_numpy() * m / ranks.to_numpy())
        fdr = np.minimum.accumulate(fdr[::-1])[::-1]
        fdr_series = pd.Series(fdr, index=order)
        out["meta_p_fdr"] = out.index.map(fdr_series).astype(float)
    else:
        out["meta_p_fdr"] = np.nan
    out["meta_fdr_sig"] = out["meta_p_fdr"] < 0.05

    verdicts = []
    for _, r in out.iterrows():
        meta_d = dict(n_pops=r["n_pops_tested"], n_directionally_consistent=r["n_directionally_consistent"])
        cov_d = dict(partial_r2=r["partial_r2"])
        dnull_d = dict(percentile=r["null_percentile"])
        verdicts.append(verdict(meta_d, cov_d, dnull_d, bool(r["meta_fdr_sig"]), float(r["meta_p"])))
    out["verdict"] = verdicts

    out.to_csv(args.out, index=False)
    print(f"\nWrote {args.out}")
    print(out[["trait", "panel", "top_snp", "meta_p", "meta_p_fdr", "partial_r2", "null_percentile", "verdict"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
