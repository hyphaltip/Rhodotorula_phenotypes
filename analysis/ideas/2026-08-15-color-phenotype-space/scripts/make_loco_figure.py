#!/usr/bin/env python3
"""LOCO sensitivity figure + summary panel.

Two-panel figure per set (gwas/gwasc):
  [A] per-trait boxplot of LOGO lambda (across 20 scaffolds) vs the full-kinship
      Tier-A genome lambda (vertical line)
  [B] Tier-A anchor loci: -log10 p in LOCO scan vs Tier-A p (1:1 identity line)
Saves PNG + PDF.  Also writes loco_sensitivity_summary.csv.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged", required=True, nargs="+",
                    help="loco_merged_{gwas,gwasc}.csv")
    ap.add_argument("--tiera-summary", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    ta = pd.read_csv(args.tiera_summary)

    anchors = [  # (label, chr, rs)
        ("chroma", 10, "scaffold_10_384905"),
        ("AUC_10", 10, "scaffold_10_396172"),
        ("resilience_30", 13, "scaffold_13_810026"),
        ("IC50_est", 16, "scaffold_16_122361"),
    ]
    rows = []
    for f in args.merged:
        pref = Path(f).stem.split("_")[-1]          # gwas | gwasc
        m = pd.read_csv(f)
        # per-trait genome LOCO lambda = median across scaffolds (GC on merged p)
        rows.append(dict(set=pref, metric="loco_lambda_fraction",
                         trait="ALL", value=m["lambda_LOCO"].median()))
        for t, g in m.groupby("trait"):
            rows.append(dict(set=pref, metric="loco_lambda_fraction",
                             trait=t, value=g["lambda_LOCO"].median()))
        # anchors: find LOCO p for rs at its chr
        for lab, chr_, rs in anchors:
            g = m[(m["trait"] == lab if lab in m["trait"].unique() else m["trait"].isin(
                [lab])) & (m["chr"] == chr_)]
            if len(g) == 0:
                continue
            r = g.iloc[0]
            # LOCO top on that chr
            lo_p = r["top_p_LOCO"] if r["top_rs_LOCO"] == rs else np.nan
            rows.append(dict(set=pref, metric="anchor_p", trait=lab,
                             value=-np.log10(max(lo_p, 1e-300)) if lo_p == lo_p else np.nan))
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "loco_sensitivity_summary.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    # --- A: lambda boxplot per set+trait ---
    lb = summary[summary.metric == "loco_lambda_fraction"]
    t_res = lb[lb.trait != "ALL"]
    sets = sorted(lb["set"].unique())
    ax = axes[0]
    xpos = []
    labs = []
    for i, s in enumerate(sets):
        tt = sorted(t_res[t_res.set == s]["trait"].unique())
        sp = [ (i * (len(tt) + 1) + j) for j in range(len(tt)) ]
        vals = [t_res[(t_res.set == s) & (t_res.trait == t)]["value"].iloc[0]
                for t in tt]
        ax.bar(sp, vals, width=0.8, label=f"{s} (LOCO λ)")
        xpos.extend(sp); labs.extend(tt)
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.set_xticks(xpos); ax.set_xticklabels(labs, rotation=45, ha="right")
    ax.set_ylabel("genome-x LOCO λ (median over scaffolds)")
    ax.set_title("(A) LOCO genomic-inflation control vs full-kinship")
    ax.legend()
    # --- B: anchor identity plot ---
    ax = axes[1]
    ta_h = pd.read_csv(args.tiera_summary).set_index("trait")["top_p"].to_dict()
    for i, s in enumerate(sets):
        pts = []
        for lab, chr_, rs in anchors:
            a = summary[(summary.set == s) & (summary.metric == "anchor_p")
                        & (summary.trait == lab) & summary.value.notna()]
            if len(a):
                pts.append((lab, rs, a["value"].iloc[0]))
        xs = [-np.log10(max(ta_h.get(lab, 1), 1e-300)) for lab, rs, _ in pts]
        ys = [v for _, _, v in pts]
        if not pts:
            continue
        ax.scatter(xs, ys, s=80, label=f"{s} (LOCO)", marker="o")
        for (lab, rs, _), x, y in zip(pts, xs, ys):
            ax.annotate(lab, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=8)
    xs_all = []
    ys_all = [v for s in sets for v in
              summary[(summary.set == s) & (summary.metric == "anchor_p") & summary.value.notna()]["value"]]
    for s in sets:
        for lab, chr_, rs in anchors:
            a = summary[(summary.set == s) & (summary.metric == "anchor_p")
                        & (summary.trait == lab) & summary.value.notna()]
            if len(a):
                xs_all.append(-np.log10(max(ta_h.get(lab, 1), 1e-300)))
    if xs_all and ys_all:
        mx = max(max(xs_all), max(ys_all))
        ax.plot([0, mx + .5], [0, mx + .5], "k--", lw=0.8)
    ax.set_xlabel("-log10 p (full-kinship Tier A)")
    ax.set_ylabel("-log10 p (LOCO)")
    ax.set_title("(B) anchor loci: LOCO vs full-kinship")
    ax.legend()
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out / f"loco_sensitivity.{ext}", dpi=150)
    print("wrote loco_sensitivity.{png,pdf} + loco_sensitivity_summary.csv")

if __name__ == "__main__":
    main()
