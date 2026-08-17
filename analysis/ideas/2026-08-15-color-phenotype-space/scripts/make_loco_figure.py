#!/usr/bin/env python3
"""LOCO sensitivity figure + summary panel.

Two-panel figure per set (gwas/gwasc):
  [A] per-trait LOCO lambda (distribution across the 20 per-scaffold scans)
      vs the full-kinship Tier-A genome lambda (x marker) and the null (1.0)
  [B] Tier-A anchor loci: -log10 p in LOCO scan vs the Tier-A p AT THE SAME SNP
      (1:1 identity line)
Saves PNG + PDF.  Also writes loco_sensitivity_summary.csv (lambda medians per trait
+ anchor same-SNP LOCO & Tier-A p).
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# only traits actually re-run under LOCO have merged rows; IC50_est is not in LOCO
anchors = [  # (label, chr, rs)  -- all ran in the LOCO set
    ("chroma", 10, "scaffold_10_384905"),
    ("AUC_10", 10, "scaffold_10_396172"),
    ("resilience_30", 13, "scaffold_13_810026"),
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged", required=True, nargs="+",
                    help="loco_merged_{gwas,gwasc}.csv")
    ap.add_argument("--tiera-summary", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    ta = pd.read_csv(args.tiera_summary)  # set,trait,lambda,...

    merged = {}
    for f in args.merged:
        pref = Path(f).stem.split("_")[-1]          # gwas | gwasc
        merged[pref] = pd.read_csv(f)

    summary_rows = []
    for pref, m in merged.items():
        # per-trait LOCO lambda median across scaffolds + full spread
        for t, g in m.groupby("trait"):
            summary_rows.append(dict(set=pref, metric="loco_lambda_fraction",
                                     trait=t, value=g["lambda_LOCO"].median()))
        summary_rows.append(dict(set=pref, metric="loco_lambda_fraction",
                                 trait="ALL", value=m["lambda_LOCO"].median()))
        # anchors: same-SNP p in LOCO (top_rs_LOCO) vs Tier-A (top_rs_TierA)
        for lab, chr_, rs in anchors:
            g = m[(m["trait"] == lab) & (m["chr"] == chr_)]
            if len(g) == 0:
                continue
            r = g.iloc[0]
            if r["top_rs_LOCO"] != rs or r["top_rs_TierA"] != rs:
                print(f"  WARN {pref} {lab} chr{chr_}: anchor {rs} not the chr top "
                      f"(LOCO={r['top_rs_LOCO']}, TierA={r['top_rs_TierA']}) -- skipped")
                continue
            summary_rows.append(dict(set=pref, metric="anchor_p", trait=lab,
                                     value=-np.log10(max(r["top_p_LOCO"], 1e-300))))
            summary_rows.append(dict(set=pref, metric="anchor_tiera_p", trait=lab,
                                     value=-np.log10(max(r["top_p_TierA"], 1e-300))))
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out / "loco_sensitivity_summary.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    # --- A: LOCO lambda distribution vs Tier-A lambda ---
    ta_l = ta.set_index(["set", "trait"])["lambda"]
    sets = sorted(merged.keys())
    ax = axes[0]
    xpos, labs = [], []
    bw = 0.38                      # bar width
    i = 0
    for s in sets:
        m = merged[s]
        tt = [t for t in ["chroma", "AUC_10", "resilience_30", "ALL"] if t in m["trait"].unique()]
        for j, t in enumerate(tt):
            x = i * 4 + j
            g = m[m["trait"] == t]["lambda_LOCO"] if t != "ALL" else m["lambda_LOCO"]
            vals = g.to_numpy()
            ax.vlines(x, vals.min(), vals.max(), color="0.75", lw=1.2, zorder=1)
            ax.scatter(np.full_like(vals, x), vals, s=10, color="0.6",
                       alpha=.6, edgecolors="none", zorder=2,
                       label=("per-scaffold LOCO λ" if (i == 0 and j == 0) else None))
            med = np.median(vals)
            ax.bar(x, med, width=bw, color="#3182bd", alpha=.85, zorder=3)
            # Tier-A genome lambda as diamond marker
            tkey = t.rstrip("ALL") or "ALL"
            key = (s, t) if (s, t) in ta_l.index else (s, "ALL")
            tl = ta_l.loc[key]
            ax.plot(x + bw / 2 + 0.05, tl, marker="D", ms=7, color="#e31a1c",
                    zorder=4, ls="none",
                    label=("full-kinship Tier-A λ" if (i == 0 and j == 0) else None))
            xpos.append(x); labs.append(f"{s}\n{t}")
        i += 1
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.set_xticks(xpos); ax.set_xticklabels(labs, fontsize=8)
    ax.set_ylabel("λ (genome scan)")
    ax.set_title("(A) LOCO per-scaffold λ vs full-kinship Tier-A λ")
    ax.legend(fontsize=7, loc="upper right")
    ax.set_ylim(0, 2.0)

    # --- B: anchor identity plot (same-SNP p both axes) ---
    ax = axes[1]
    ann_seen = {}
    for s in sets:
        for lab, chr_, rs in anchors:
            a = summary[(summary.set == s) & (summary.metric == "anchor_p")
                        & (summary.trait == lab)]
            b = summary[(summary.set == s) & (summary.metric == "anchor_tiera_p")
                        & (summary.trait == lab)]
            if len(a) == 0 or len(b) == 0:
                continue
            x = float(b.loc[:, "value"].iloc[0]); y = float(a.loc[:, "value"].iloc[0])
            ax.scatter(x, y, s=80, marker="o", edgecolors="k", linewidths=.5,
                       c={"gwas": "#1b9e77", "gwasc": "#d95f02"}[s],
                       label=f"{s} (LOCO)" if s not in ann_seen else None)
            ann_seen.setdefault(s, True)
            tag = lab if (x, y) not in {k[:2] for k in ann_seen} else f"{lab} ({s})"
            ann_seen[(x, y, lab)] = True
            dy = 4 if lab == "AUC_10" else 4
            ax.annotate(tag, (x, y), xytext=(5, dy), textcoords="offset points",
                        fontsize=8)
    if len(summary[summary.metric == "anchor_p"]):
        xs = summary.loc[summary.metric == "anchor_tiera_p", "value"].to_numpy()
        ys = summary.loc[summary.metric == "anchor_p", "value"].to_numpy()
        if len(xs) and len(ys):
            mx = max(xs.max(), ys.max()) + 0.5
            ax.plot([0, mx], [0, mx], "k--", lw=0.8)
            ax.set_xlim(0, mx); ax.set_ylim(0, mx)
    ax.set_xlabel("-log10 p (full-kinship Tier A, same SNP)")
    ax.set_ylabel("-log10 p (LOCO)")
    ax.set_title("(B) anchor loci: LOCO vs Tier-A same-SNP")
    ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out / f"loco_sensitivity.{ext}", dpi=150)
    print("wrote loco_sensitivity.{png,pdf} + loco_sensitivity_summary.csv")

if __name__ == "__main__":
    main()
