#!/usr/bin/env python3
"""dxy/Fst co-localization with GWAS (Tier-A FDR + anchor loci).

Joins:
  - Tier-A FDR-significant SNP loci (any trait)
  - Known anchor loci
onto pixy 100-kb windows (mean dxy & mean Fst across the 15 pop pairs)
and tests whether significant GWAS loci are enriched in high-dxy / high-Fst
windows vs the genome-wide expectation.

Outputs:
  coloc_summary.csv   -- one row per GWAS locus: window, mean dxy, mean Fst,
                         high_dxy (top-20% genome-wide), high_fst (top-20%)
  coloc_enrichment.txt -- Fisher exact enrichment table
"""
import argparse
from pathlib import Path
import pandas as pd
from scipy import stats

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fdr-dir", required=True, help="dir of *_{trait}_fdr05_sig.txt")
    ap.add_argument("--dxy", required=True)
    ap.add_argument("--fst", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--win-size", type=int, default=100000)
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    dxy = pd.read_csv(args.dxy, sep="\t")
    fst = pd.read_csv(args.fst, sep="\t")
    dxy_w = dxy.groupby(["chromosome", "window_pos_1"])["avg_dxy"].mean().reset_index()
    fst_w = fst.groupby(["chromosome", "window_pos_1"])["avg_hudson_fst"].mean().reset_index()
    win = dxy_w.merge(fst_w, on=["chromosome", "window_pos_1"])
    win["high_dxy"] = win["avg_dxy"] >= win["avg_dxy"].quantile(0.80)
    win["high_fst"] = win["avg_hudson_fst"] >= win["avg_hudson_fst"].quantile(0.80)
    print(f"{len(win)} windows; cutoff dxy={win['avg_dxy'].quantile(0.80):.4f} "
          f"fst={win['avg_hudson_fst'].quantile(0.80):.4f}")

    anchors = [  # (trait, chr, rs)
        ("chroma", 10, "scaffold_10_384905"),
        ("AUC_10", 10, "scaffold_10_396172"),
        ("resilience_30", 13, "scaffold_13_810026"),
        ("IC50_est", 16, "scaffold_16_122361"),
        ("chroma(cull)", 3, "scaffold_3_570085"),
    ]

    loci = []
    for f in sorted(Path(args.fdr_dir).glob("*_fdr05_sig.txt")):
        trait = f.name.replace("_fdr05_sig.txt", "")
        t = pd.read_csv(f, sep="\t")
        t["trait"] = trait
        loci.append(t[["trait", "chr", "rs", "ps", "p_wald"]])
    gwas = pd.concat(loci, ignore_index=True)
    gwas["chr_s"] = "scaffold_" + gwas["chr"].astype(str)
    n_loci = len(gwas)
    print(f"{n_loci} FDR-significant SNP loci (any trait)")

    # window each SNP
    def window_of(row):
        m = win[(win["chromosome"] == row["chr_s"]) &
                (win["window_pos_1"] <= row["ps"]) &
                (row["ps"] < win["window_pos_1"] + args.win_size)]
        return m.index[0] if len(m) else pd.NA
    print("... mapping loci to windows")
    gwas["win_idx"] = gwas.apply(window_of, axis=1)
    g = gwas.dropna(subset=["win_idx"])
    print(f"mapped {len(g)} / {n_loci} loci to windows")
    res = g.merge(win.add_prefix("w_"), left_on="win_idx", right_index=True)
    res = res.rename(columns={"w_avg_dxy": "mean_dxy", "w_avg_hudson_fst": "mean_fst",
                              "w_high_dxy": "high_dxy", "w_high_fst": "high_fst"})
    res = res.drop(columns=["w_chromosome", "w_window_pos_1"])
    # Fisher on WINDOWS: of windows containing >=1 FDR SNP, what fraction are
    # high-dxy/high-fst vs windows with no FDR SNP (genome-wide baseline)
    res.to_csv(out / "coloc_final.csv", index=False)
    hit_win = set(g["win_idx"].astype(int))
    win_full = win.copy()
    win_full["has_gwas"] = [i in hit_win for i in range(len(win))]
    n_gwas_hi_dxy = int(win_full[win_full["has_gwas"]]["high_dxy"].sum())
    n_gwas_hi_fst = int(win_full[win_full["has_gwas"]]["high_fst"].sum())
    n_gwas_w = int(win_full["has_gwas"].sum())
    n_win = len(win_full)
    print(f"{n_gwas_w}/{n_win} windows have >=1 FDR SNP")
    print(f"  of which {n_gwas_hi_dxy} high-dxy, {n_gwas_hi_fst} high-fst "
          f"(expect ~{int(win_full['high_dxy'].sum()*n_gwas_w/n_win)}, "
          f"~{int(win_full['high_fst'].sum()*n_gwas_w/n_win)})")
    lines = []
    for col, lab in [("high_dxy", "high-dxy"),
                     ("high_fst", "high-fst")]:
        a = int(win_full[win_full["has_gwas"]][col].sum())
        b = int(win_full["has_gwas"].sum()) - a
        c = int(win_full["high_" + ("dxy" if col == "high_dxy" else "fst")].sum()) - a if col == "high_dxy" else int(win_full[col].sum()) - a
        # simpler: standard 2x2 (has_gwas, topdiff)
        tbl = [[a, b], [int(win_full[col].sum()) - a, n_win - int(win_full[col].sum()) - b]]
        or_, p = stats.fisher_exact(tbl)
        lines.append(f"  {lab}: GWAS-window {a}/{n_gwas_w} vs baseline "
                     f"{int(win_full[col].sum())}/{n_win} -> OR={or_:.3f} p={p:.3g}")
        print(lines[-1])
    with open(out / "coloc_enrichment.txt", "w") as fh:
        fh.write(f"GWAS FDR loci co-localization with pixy high-diff windows\n")
        fh.write(f"windows total={n_win} high_dxy={int(win_full['high_dxy'].sum())} "
                 f"high_fst={int(win_full['high_fst'].sum())}\n")
        fh.write(f"GWAS loci mapped={len(res)} in {n_gwas_w} windows\n")
        for ln in lines:
            fh.write(ln + "\n")

    print("\nAnchor loci windows:")
    av = []
    for trait, chr_, rs in anchors:
        hit = gwas[(gwas["chr"] == chr_) & (gwas["rs"] == rs)]
        if len(hit):
            row = hit.iloc[0]
            m = win[(win["chromosome"] == "scaffold_" + str(chr_)) &
                    (win["window_pos_1"] <= row["ps"]) &
                    (row["ps"] < win["window_pos_1"] + args.win_size)]
            if len(m):
                av.append(dict(trait=trait, rs=rs, ps=row["ps"], mean_dxy=m.iloc[0]["avg_dxy"],
                               mean_fst=m.iloc[0]["avg_hudson_fst"],
                               high_dxy=m.iloc[0]["high_dxy"], high_fst=m.iloc[0]["high_fst"]))
                print(f"  {trait:12s} {rs}: dxy={m.iloc[0]['avg_dxy']:.4f} "
                      f"fst={m.iloc[0]['avg_hudson_fst']:.4f} "
                      f"hi_dxy={m.iloc[0]['high_dxy']} hi_fst={m.iloc[0]['high_fst']}")
    pd.DataFrame(av).to_csv(out / "coloc_anchors.csv", index=False)
    print("\nwrote coloc_final.csv, coloc_anchors.csv, coloc_enrichment.txt")

if __name__ == "__main__":
    main()
