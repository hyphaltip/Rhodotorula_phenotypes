#!/usr/bin/env python3
"""Candidate SNP and gene panels for the Cr traits, with window context.
 Panel A: lineage-level partitions behind the Cr Bonferroni loci of scan (b) (relatedness + run adjustment): the SNPs that share the lead SNP's genotype pattern (perfectly linked), with their snpEff impacts and genes.
 Panel B: within-lineage candidates: scan (c) SNPs (run adjustment + lineage covariates) below p = 1e-3 for Cr traits, and the best SNPs of the per-lineage permutation scans (L01, L02, L03).
 Locus-zoom figures: -log10 p of scan (c) and scan (b) in a window, snpEff impact, and the gene models."""
import re, sys, gzip
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
sys.path.insert(0, "analysis/gwas_dh4148/scripts")
from geno_lib import *
G = R / "gwas"; T = Path("analysis/gwas_dh4148/report/tables"); F = Path("analysis/gwas_dh4148/report/figures")
CR = ["Cr_da_top", "Cr_a_auc", "Cr_relarea_top", "Cr_relarea_auc", "Cr_relrgr_auc", "Cr_logIC50"]; WIN = 40000
KEY = re.compile(r"chromat|chromium|sulfate|sulphate|glutathione|thioredoxin|superoxide|catalase|metal|zinc|iron|copper|cadmium|arsen|abc transporter|p-type atpase|flavin|nadph|carotenoid|phytoene|prenyl|isoprenoid|terpen|heat shock|stress|map kinase|vacuol|peroxid|ferric|transporter", re.I)
IMP = {"HIGH": "#D55E00", "MODERATE": "#E69F00", "LOW": "#009E73", "MODIFIER": "#999999", "": "#999999"}
log = lambda m: print(m, flush=True)
ids, X, sid, names = load_geno(); ix = {k: i for i, k in enumerate(ids)}; genes = load_genes()
def read(out, t): a = pd.read_csv(G / out / f"{t}.assoc.txt.gz", sep="\t"); a["chr"], a["pos"] = split_id(a.rs.values); return a[["rs", "chr", "pos", "af", "beta", "se", "p_wald"]]
scb = {t: read("output_adj", t) for t in CR}; scc = {t: read("output_adjlin", t) for t in CR}
# ---- Panel A: lineage partitions
rowsA = []; need = set()
bl = pd.concat([pd.read_csv(T / "gwas_loci_runadj.csv").assign(scan="run-adjusted"), pd.read_csv(T / "gwas_loci_unadjusted.csv").query("bonferroni_sig").assign(scan="unadjusted")])
bl = bl[bl.trait.isin(CR)].sort_values("p_wald")   # run-adjusted scan has no Bonferroni locus for Cr; its suggestive loci and the unadjusted Bonferroni loci are used
seen = set(); parts = []
for _, r in bl.iterrows():
    g = X[ix[r.lead_snp]]; key = tuple(np.where(np.isnan(g), -1, g).astype(int))
    if key in seen: continue
    seen.add(key); parts.append((f'{r.trait} [{r.scan}]', r.lead_snp, g)); 
    if len(parts) >= 8: break
ann_all = {}
with gzip.open(GENO / "panel_sites_ann.tsv.gz", "rt") as f:
    for line in f:
        c = line.rstrip("\n").split("\t"); ann_all[f"{c[0]}:{c[1]}:{c[2]}:{c[3]}"] = (c[3], c[6] if len(c) > 6 else "")
lin = pd.read_csv(R / "lineages.csv").set_index("strain_id").loc[sid, "lineage"].values
for k, (t, snp, g) in enumerate(parts, 1):
    ok = ~np.isnan(g); same = np.where(np.all(((X[:, ok] == g[ok]) | (X[:, ok] == 1 - g[ok]) & False) | (X[:, ok] == g[ok]), axis=1))[0]
    comp = np.where(np.all(X[:, ok] == 1 - g[ok], axis=1))[0]; members = np.concatenate([same, comp]); alt_lin = sorted(set(lin[(g == 1)]))
    ann = [best_ann(ann_all[ids[i]][1], ann_all[ids[i]][0]) for i in members if ids[i] in ann_all]
    hm = pd.DataFrame(ann, columns=["effect", "impact", "gene_id", "hgvs_p"]); hm = hm[hm.impact.isin(["HIGH", "MODERATE"])]
    gp = genes.set_index("gene_id")["product"]; hm["product"] = hm.gene_id.map(gp)
    kg = sorted({f"{r.gene_id} ({r['product'][:35]})" for _, r in hm.iterrows() if isinstance(r["product"], str) and KEY.search(r["product"])})
    rowsA.append(dict(partition=f"P{k}", lead_trait=t, lead_snp=snp, strains_with_alt=int((g == 1).sum()), lineages_with_alt=",".join(alt_lin), n_perfectly_linked_snps=len(members), n_high=int((hm.impact == "HIGH").sum()), n_moderate=int((hm.impact == "MODERATE").sum()),
                      n_genes_high_moderate=hm.gene_id.nunique(), keyword_genes_high_moderate="; ".join(kg[:12])))
A = pd.DataFrame(rowsA); A.to_csv(T / "cr_lineage_partition_panels.csv", index=False); log(A.drop(columns=["keyword_genes_high_moderate"]).to_string(index=False))
# ---- Panel B: within-lineage candidates
rowsB = []
for t in CR:
    a = scc[t]; s = a[a.p_wald < 1e-3]
    for _, r in s.iterrows(): rowsB.append(dict(source="scan c (all lineages, relatedness + lineage covariates)", trait=t, chr=r.chr, pos=int(r.pos), snp=r.rs, af=r.af, beta=r.beta, p=r.p_wald, p_family_wise=np.nan))
wl = pd.read_csv(T / "within_lineage_scan_top_snps.csv")
for _, r in wl[wl.trait.isin(CR)].iterrows():
    if r.p < 1e-2: rowsB.append(dict(source=f"within {r.lineage} (permutation scan)", trait=r.trait, chr=r.chr, pos=r.pos, snp=f"{r.chr}:{r.pos}:{r.ref}:{r.alt}", af=np.nan, beta=r.beta, p=r.p, p_family_wise=r.p_family_wise))
B = pd.DataFrame(rowsB)
def near(c, p, w=5000):
    g = genes[(genes.chr == c) & (genes.end >= p - w) & (genes.start <= p + w)]; return g
ann_cache = ann_for([s for s in B.snp if s in ix])
out = []
for _, r in B.iterrows():
    e = ann_cache.get(r.snp, ("", "", "", "")); g = near(r.chr, r.pos)
    d = dict(r); d.update(effect=e[0], impact=e[1], hgvs_p=e[3], genes_within_5kb="; ".join(g.gene_id), products="; ".join(g["product"].str.slice(0, 40)), keyword_hit=bool(g["product"].map(lambda s: bool(KEY.search(s))).any()) if len(g) else False); out.append(d)
B = pd.DataFrame(out).sort_values("p"); B.to_csv(T / "cr_candidate_snp_panel.csv", index=False); log(f"panel B: {len(B)} SNP rows ({B.snp.nunique()} SNPs, {B.genes_within_5kb.replace('', np.nan).dropna().str.split('; ').explode().nunique()} genes)")
# ---- enrichment of keyword genes among scan-c Cr signal (permutation of gene labels)
gm = []
for _, gr in genes.iterrows():
    best = 0.0
    for t in CR:
        a = scc[t]; s = a[(a.chr == gr.chr) & (a.pos >= gr.start - 2000) & (a.pos <= gr.end + 2000)]
        if len(s): best = max(best, -np.log10(s.p_wald.min()))
    gm.append(best)
genes["score"] = gm; genes["kw"] = genes["product"].map(lambda s: bool(KEY.search(s))); tested = genes[genes.score > 0]; obs = tested[tested.kw].score.mean(); rng = np.random.default_rng(3)
null = np.array([rng.choice(tested.score.values, int(tested.kw.sum()), replace=False).mean() for _ in range(10000)])
E1 = pd.DataFrame([dict(genes_with_snps=len(tested), keyword_genes=int(tested.kw.sum()), mean_neglog10p_keyword=obs, mean_neglog10p_all=tested.score.mean(), perm_p=float((null >= obs).mean()))]); E1.to_csv(T / "cr_keyword_gene_enrichment.csv", index=False); log(E1.round(4).to_string(index=False))
# ---- locus-zoom figures
loci = []
cand = B[B.source.str.startswith("scan c")].sort_values("p")
for _, r in cand.iterrows():
    if all(not (r.chr == c and abs(r.pos - p) < 2 * WIN) for c, p, _ in loci): loci.append((r.chr, r.pos, r.trait))
    if len(loci) >= 6: break
for _, r in B[B.source.str.startswith("within")].sort_values("p").iterrows():
    if all(not (r.chr == c and abs(r.pos - p) < 2 * WIN) for c, p, _ in loci): loci.append((r.chr, r.pos, r.trait + " " + r.source.split(" ")[1]))
    if len(loci) >= 9: break
files = []
for k, (c, p, t) in enumerate(loci, 1):
    tt = t.split(" ")[0]; lo, hi = p - WIN, p + WIN; fig, ax = plt.subplots(3, 1, figsize=(9, 6.2), sharex=True, gridspec_kw=dict(height_ratios=[2, 2, 1.4]))
    for a_, (sc, lab) in zip(ax[:2], ((scc[tt], "scan (c): run-adjusted + lineage covariates"), (scb[tt], "scan (b): run-adjusted, relatedness only"))):
        s = sc[(sc.chr == c) & (sc.pos >= lo) & (sc.pos <= hi)]; a_.scatter(s.pos, -np.log10(s.p_wald), s=8, color="#0072B2"); a_.set_ylabel(f"-log10 p\n{lab}", fontsize=7); a_.axvline(p, color="r", lw=.5, ls=":")
    wi = wl[(wl.chr == c) & (wl.pos >= lo) & (wl.pos <= hi) & (wl.trait == tt)]
    if len(wi): ax[0].scatter(wi.pos, -np.log10(wi.p), s=26, marker="^", color="#D55E00", label="within-lineage scan"); ax[0].legend(fontsize=6)
    gg = genes[(genes.chr == c) & (genes.end >= lo) & (genes.start <= hi)]
    for j, (_, g) in enumerate(gg.iterrows()):
        y = j % 3; s0, e0 = max(g.start, lo), min(g.end, hi); ax[2].add_patch(plt.Rectangle((s0, y), e0 - s0, 0.55, color="#D55E00" if KEY.search(g["product"]) else "#56B4E9")); ax[2].text((s0 + e0) / 2, y + 0.65, (g.symbol or g.gene_id.replace("ACY3AU_", "")) + (" +" if g.strand == "+" else " -"), fontsize=5, ha="center")
    ax[2].set_ylim(-0.2, 3.2); ax[2].set_yticks([]); ax[2].set_xlim(lo, hi); ax[2].set_xlabel(f"{c} position (bp); orange = product matches the metal / stress / pigment keywords"); fig.suptitle(f"Locus {k}: {tt} near {c}:{p:,}", fontsize=9); fig.tight_layout(); fn = f"zoom_cr_locus{k}"; fig.savefig(F / f"{fn}.png", dpi=150); plt.close(fig); files.append((fn, c, p, t, "; ".join(gg.gene_id[:8])))
pd.DataFrame(files, columns=["figure", "chr", "pos", "trait", "genes_in_window"]).to_csv(T / "cr_zoom_loci.csv", index=False); log(f"{len(files)} zoom figures")
