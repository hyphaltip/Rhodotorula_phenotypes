#!/usr/bin/env python3
"""MAS-usability gates for the population-vs-locus disambiguation's `likely_real` loci
(D-17, GWAS.md S13), deferred at that point and run here per user request:

  1. Minimum effect size: genotype effect in phenotype-SD units (not just p-value).
  2. LD-decay check: is the flagged SNP the best candidate in its own LD block, or a
     distant tag of a stronger nearby signal? Uses plink2 pairwise r^2 in a +/-100kb
     window plus the existing Tier A p-values in that window (no new GEMMA runs).
  3. Plate/batch-confound check: every strain in this panel was measured in exactly ONE
     of 4 runs (353-356) -- run_number is a fixed per-strain batch label, not a repeated
     covariate, so it is a real confounding risk if genotype groups cluster by run.
     Tests genotype~run association (Fisher/chi2) and refits the population-covariate
     regression with run_number added as an extra covariate.
  4. Timepoint robustness (color traits only: lab_L/a/b, sat -- cu_dose_slope and AUC_20
     are dose-response-derived, not single-timepoint color reads, so are out of scope
     for this specific check): recomputes each trait at early/mid/late growth windows
     and re-runs the within-population test (same logic as check_population_vs_locus.py
     test A) at each window, to see whether the genotype effect is a stable feature of
     the growth curve or a one-off artifact of the single window originally used.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_population_vs_locus import (  # noqa: E402
    genotype_at_site, load_pop_map, load_culled_set, parse_snp,
    near_clone_collapse_within_pop, within_pop_test, meta_analysis,
)


# ---------------------------------------------------------------------------
# Gate 1: minimum effect size
# ---------------------------------------------------------------------------

def effect_size(dosage: dict, pheno: dict) -> dict:
    common = [s for s in dosage if s in pheno]
    g = np.array([dosage[s] for s in common])
    y = np.array([pheno[s] for s in common])
    g0 = y[g < 0.5]
    g1 = y[g >= 0.5]
    if len(g0) < 2 or len(g1) < 2:
        return dict(n0=len(g0), n1=len(g1), mean0=np.nan, mean1=np.nan,
                    pooled_sd=np.nan, cohens_d=np.nan, pheno_sd=np.nan, frac_sd_shift=np.nan)
    pooled_sd = np.sqrt(((len(g0) - 1) * g0.var(ddof=1) + (len(g1) - 1) * g1.var(ddof=1))
                        / (len(g0) + len(g1) - 2))
    d = (g1.mean() - g0.mean()) / pooled_sd if pooled_sd > 0 else np.nan
    pheno_sd = y.std(ddof=1)
    return dict(n0=len(g0), n1=len(g1), mean0=float(g0.mean()), mean1=float(g1.mean()),
                pooled_sd=float(pooled_sd), cohens_d=float(d), pheno_sd=float(pheno_sd),
                frac_sd_shift=float((g1.mean() - g0.mean()) / pheno_sd) if pheno_sd > 0 else np.nan)


# ---------------------------------------------------------------------------
# Gate 2: LD-decay / best-candidate-in-block check
# ---------------------------------------------------------------------------

def ld_decay_check(plink2_bin: str, bfile: str, rs: str, scaffold: str, pos: int,
                    trait_assoc: pd.DataFrame, work_prefix: str, window_kb: int = 100):
    """This cluster's plink2 build has neither --distance nor a windowed --r2/--ld-snp
    report (only pairwise --ld between two named SNPs) -- L-25's gotcha strikes again.
    Uses --clump instead: restrict the input assoc list to the +/-window_kb region
    around the flagged SNP, then let plink2 group them into LD-based clumps (r^2>=0.2,
    unphased). If the flagged SNP ends up as a SECONDARY member of a clump led by a
    more significant SNP, that lead is a better candidate marker for the same signal.
    """
    window = trait_assoc[(trait_assoc["chr"].astype(str) == str(scaffold).replace("scaffold_", ""))
                          & (trait_assoc["ps"].between(pos - window_kb * 1000, pos + window_kb * 1000))]
    if len(window) < 2:
        return dict(n_snps_in_window=len(window), is_clump_lead=None, n_secondary_in_clump=None,
                     better_candidate_rs=None, better_candidate_p=np.nan,
                     region_kb_p_lt_5eminus2=np.nan, note="too few SNPs in window")

    assoc_file = work_prefix + "_window_assoc.tsv"
    window[["rs", "p_wald"]].rename(columns={"rs": "ID", "p_wald": "P"}).to_csv(
        assoc_file, sep="\t", index=False
    )
    r = subprocess.run(
        [plink2_bin, "--bfile", bfile, "--clump", assoc_file,
         "--clump-id-field", "ID", "--clump-p-field", "P",
         "--clump-p1", "1", "--clump-p2", "0.05", "--clump-r2", "0.2", "--clump-kb", str(window_kb),
         "--clump-unphased", "--out", work_prefix, "--allow-extra-chr"],
        capture_output=True, text=True,
    )
    clumps_path = work_prefix + ".clumps"
    if r.returncode != 0 or not os.path.exists(clumps_path):
        return dict(n_snps_in_window=len(window), is_clump_lead=None, n_secondary_in_clump=None,
                     better_candidate_rs=None, better_candidate_p=np.nan,
                     region_kb_p_lt_5eminus2=np.nan, note=f"plink2 --clump failed: {r.stdout[-500:]}")
    clumps = pd.read_csv(clumps_path, sep="\t")
    clumps.columns = [c.upper() for c in clumps.columns]

    flagged_p = float(trait_assoc.loc[trait_assoc.rs == rs, "p_wald"].iloc[0])
    is_lead = (clumps["ID"] == rs).any()
    if is_lead:
        lead_row = clumps.loc[clumps["ID"] == rs].iloc[0]
        n_secondary = int(lead_row.get("TOTAL", 0))
        better_rs, better_p = None, np.nan
    else:
        lead_row = None
        for _, cr in clumps.iterrows():
            sp2 = str(cr.get("SP2", ""))
            if rs in sp2:
                lead_row = cr
                break
        n_secondary = int(lead_row.get("TOTAL", 0)) if lead_row is not None else 0
        better_rs = lead_row["ID"] if lead_row is not None else None
        better_p = float(lead_row["P"]) if lead_row is not None else np.nan

    region_kb = 2 * window_kb  # by construction of the input window

    return dict(n_snps_in_window=len(window), is_clump_lead=bool(is_lead), n_secondary_in_clump=n_secondary,
                better_candidate_rs=better_rs, better_candidate_p=better_p,
                region_kb_p_lt_5eminus2=float(region_kb), note="")


# ---------------------------------------------------------------------------
# Gate 3: plate/batch confounding
# ---------------------------------------------------------------------------

def batch_confound_check(dosage: dict, pheno: dict, pop: dict, run_of: dict, strains: list):
    import statsmodels.formula.api as smf

    rows = []
    for s in strains:
        if s in dosage and s in pheno and pop.get(s) and s in run_of:
            rows.append(dict(strain=s, genotype=dosage[s], phenotype=pheno[s],
                              population=pop[s], run=str(run_of[s])))
    df = pd.DataFrame(rows)
    if len(df) < 10 or df["run"].nunique() < 2:
        return dict(n=len(df), n_runs=df["run"].nunique() if len(df) else 0,
                     geno_vs_run_p=np.nan, partial_p_with_run=np.nan, partial_p_without_run=np.nan,
                     note="insufficient data")

    # Is genotype itself associated with run (independent of phenotype)?
    ct = pd.crosstab(df["genotype"] >= 0.5, df["run"])
    try:
        chi2, geno_run_p, _, _ = stats.chi2_contingency(ct)
    except Exception:
        geno_run_p = np.nan

    without_run = smf.ols("phenotype ~ genotype + C(population)", data=df).fit()
    with_run = smf.ols("phenotype ~ genotype + C(population) + C(run)", data=df).fit()
    p_without = float(without_run.pvalues.get("genotype", np.nan))
    p_with = float(with_run.pvalues.get("genotype", np.nan))

    return dict(n=len(df), n_runs=df["run"].nunique(), geno_vs_run_p=geno_run_p,
                partial_p_without_run=p_without, partial_p_with_run=p_with, note="")


# ---------------------------------------------------------------------------
# Gate 4: timepoint robustness (color traits only)
# ---------------------------------------------------------------------------

TIMEPOINT_WINDOWS = {"early": (18.0, 42.0), "mid": (48.0, 72.0), "late": (85.0, 110.0)}
LAB_COL = {"lab_L": "ColorLab_L*Median", "lab_a": "ColorLab_a*Median", "lab_b": "ColorLab_b*Median",
           "sat": "ColorHSV_SaturationMedian"}


def build_windowed_phenotype(extract: pd.DataFrame, trait: str, win_lo: float, win_hi: float) -> dict:
    col = LAB_COL[trait]
    d = extract[(extract.copper_mm == 0) & (extract.tp_h >= win_lo) & (extract.tp_h <= win_hi)].copy()
    d[col] = d[col].astype(float)
    plate = (d.groupby(["strain_code", "run_number", "plate_number"])[col].median().reset_index())
    plate_mean = plate.groupby("strain_code")[col].mean()
    return plate_mean.to_dict()


def timepoint_robustness(extract: pd.DataFrame, trait: str, dosage: dict, pop: dict, culled: set,
                          rep_map: dict) -> pd.DataFrame:
    rows = []
    for window_name, (lo, hi) in TIMEPOINT_WINDOWS.items():
        pheno_w = build_windowed_phenotype(extract, trait, lo, hi)
        a_df = within_pop_test(dosage, pheno_w, pop, culled, rep_map)
        meta = meta_analysis(a_df)
        rows.append(dict(window=window_name, win_lo=lo, win_hi=hi, meta_beta=meta["meta_beta"],
                          meta_p=meta["meta_p"], n_pops=meta["n_pops"],
                          n_consistent=meta["n_directionally_consistent"]))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--population-vs-locus-csv", required=True)
    ap.add_argument("--only-verdict", default="likely_real")
    ap.add_argument("--pheno-csv", required=True)
    ap.add_argument("--pop-csv", required=True)
    ap.add_argument("--culled-list", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--plink2", default="plink2")
    ap.add_argument("--gwas-bfile", required=True, help="full unpruned QC'd bfile, gwas panel")
    ap.add_argument("--gwasc-bfile", required=True)
    ap.add_argument("--pruned-bfile-gwasc", required=True)
    ap.add_argument("--assoc-csv-dir", required=True)
    ap.add_argument("--extract-parquet", required=True, help="db_extract.parquet for timepoint robustness")
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(args.work_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    loci = pd.read_csv(args.population_vs_locus_csv)
    loci = loci[loci["verdict"] == args.only_verdict].reset_index(drop=True)
    print(f"Loci passing MAS gates check ({args.only_verdict}): {len(loci)}")

    pop = load_pop_map(args.pop_csv)
    culled = load_culled_set(args.culled_list)
    pheno_df = pd.read_csv(args.pheno_csv)

    print("Loading run_number batch labels ...")
    extract = pd.read_parquet(
        args.extract_parquet,
        columns=["species", "strain_code", "tp_h", "copper_mm", "run_number", "plate_number",
                 "ColorLab_L*Median", "ColorLab_a*Median", "ColorLab_b*Median", "ColorHSV_SaturationMedian"],
    )
    extract = extract[extract.species == "Rhodotorula mucilaginosa"]
    run_of = extract.drop_duplicates("strain_code").set_index("strain_code")["run_number"].to_dict()

    print("Building per-population near-clone-collapse maps (reused from population_vs_locus) ...")
    rep_map = {}
    for p in sorted(set(pop.get(s) for s in culled if pop.get(s))):
        members = [s for s in culled if pop.get(s) == p]
        work_prefix = os.path.join(args.work_dir, f"collapse_pop{p}")
        rep_map[p] = near_clone_collapse_within_pop(
            members, args.pruned_bfile_gwasc, args.plink2, threshold=0.005, work_prefix=work_prefix
        )

    results = []
    for _, row in loci.iterrows():
        trait, panel, rs = row["trait"], row["panel"], row["top_snp"]
        scaffold, pos = parse_snp(rs)
        print(f"\n=== MAS gates: {trait} [{panel}] {rs} ===")
        bfile = args.gwas_bfile if panel == "gwas" else args.gwasc_bfile

        dosage = genotype_at_site(args.bcftools, args.vcf, scaffold, pos)
        pheno = dict(zip(pheno_df["strain_code"], pheno_df[trait]))
        pheno = {s: v for s, v in pheno.items() if pd.notna(v)}

        # --- Gate 1 ---
        es = effect_size(dosage, pheno)
        print(f"  effect size: Cohen's d={es['cohens_d']:.3f}  frac_sd_shift={es['frac_sd_shift']:.3f}  "
              f"(n0={es['n0']}, n1={es['n1']})")

        # --- Gate 2 ---
        assoc_path = os.path.join(args.assoc_csv_dir, f"{panel}_{trait}_assoc.csv.gz")
        trait_assoc = pd.read_csv(assoc_path, usecols=["rs", "chr", "ps", "p_wald"])
        ld_prefix = os.path.join(args.work_dir, f"ld_{trait}_{panel}")
        ld = ld_decay_check(args.plink2, bfile, rs, scaffold, pos, trait_assoc, ld_prefix)
        print(f"  LD block: {ld['n_snps_in_window']} SNPs in +/-100kb window, is_clump_lead={ld['is_clump_lead']}, "
              f"n_secondary={ld['n_secondary_in_clump']}, better candidate={ld['better_candidate_rs']} "
              f"(p={ld['better_candidate_p']})")

        # --- Gate 3 ---
        bc = batch_confound_check(dosage, pheno, pop, run_of, list(culled))
        print(f"  batch confound: geno~run p={bc['geno_vs_run_p']}, "
              f"genotype partial-p without run={bc['partial_p_without_run']:.3g}, "
              f"with run={bc['partial_p_with_run']:.3g}" if bc['partial_p_with_run'] == bc['partial_p_with_run']
              else f"  batch confound: {bc['note']}")

        # --- Gate 4 (color traits only) ---
        tp_rows = None
        if trait in LAB_COL:
            tp_df = timepoint_robustness(extract, trait, dosage, pop, culled, rep_map)
            tp_out = os.path.join(os.path.dirname(args.out), f"timepoint_{trait}_{panel}.csv")
            tp_df.to_csv(tp_out, index=False)
            print("  timepoint robustness:")
            print(tp_df.to_string(index=False))
            tp_rows = tp_df

        result = dict(
            trait=trait, panel=panel, top_snp=rs,
            cohens_d=es["cohens_d"], frac_sd_shift=es["frac_sd_shift"], n0=es["n0"], n1=es["n1"],
            ld_n_snps_window=ld["n_snps_in_window"], ld_is_clump_lead=ld["is_clump_lead"],
            ld_n_secondary=ld["n_secondary_in_clump"],
            ld_better_candidate=ld["better_candidate_rs"], ld_better_candidate_p=ld["better_candidate_p"],
            batch_geno_run_p=bc["geno_vs_run_p"], batch_partial_p_without_run=bc["partial_p_without_run"],
            batch_partial_p_with_run=bc["partial_p_with_run"],
            timepoint_all_significant=(bool((tp_rows["meta_p"] < 0.05).all()) if tp_rows is not None else None),
            timepoint_all_consistent_sign=(bool((np.sign(tp_rows["meta_beta"].dropna()) ==
                                                  np.sign(tp_rows["meta_beta"].dropna().iloc[0])).all())
                                            if tp_rows is not None and tp_rows["meta_beta"].notna().sum() > 1 else None),
        )
        results.append(result)

    out = pd.DataFrame(results)
    out.to_csv(args.out, index=False)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
