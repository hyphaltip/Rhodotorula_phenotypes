#!/usr/bin/env python3
"""Phylogenetic signal of baseline a*, size-adjusted baseline a*, change in a* at top dose, and baseline colony size.
Pagel's lambda (ML under Brownian motion) and Blomberg's K (999 permutations), on the PHYling FastTree (278 tips; 266 strains have a unique tip).
Scopes: all strains with a tip, and R. mucilaginosa only (is there structure inside the species?)."""
import re, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "analysis/carotenoid_stress_vs_baseline/scripts/strat")
from strat_common import *

rng = np.random.default_rng(7)
TREE = "data/raw/rhodotorula-phyling-protein-tree/protein-Rhodotorula-taxa_278.fungi_odb10.fasttree.nosupport.treefile"

def parse(s):
    s = s.strip().rstrip(";"); par, ln, name, kids = [-1], [0.0], [""], [[]]; cur = 0; i = 0; buf = ""
    def flush(node, buf):
        m = re.match(r"^([^:]*):?([0-9eE.+-]*)$", buf)
        if m: name[node] = m.group(1) if not kids[node] else name[node]; ln[node] = float(m.group(2)) if m.group(2) else 0.0
    for ch in s:
        if ch == "(":
            par.append(cur); ln.append(0.0); name.append(""); kids.append([]); kids[cur].append(len(par) - 1); cur = len(par) - 1
        elif ch == ",":
            flush(cur, buf); buf = ""; p = par[cur]; par.append(p); ln.append(0.0); name.append(""); kids.append([]); kids[p].append(len(par) - 1); cur = len(par) - 1
        elif ch == ")":
            flush(cur, buf); buf = ""; cur = par[cur]
        else: buf += ch
    flush(cur, buf)
    return par, ln, name, kids

par, ln, name, kids = parse(open(TREE).read())

def reroot_on_clade(par, ln, name, kids, prefixes):
    """Reroot on the clade formed by the outgroup tips (names starting with any prefix). The root goes in the middle of the branch that separates
    the outgroup clade from all other tips (ape::root(outgroup, resolve.root = TRUE) does the same), so the outgroup is sister to the ingroup.
    Degree-2 nodes left behind are removed. The outgroup tips must form a clade in the input tree (either side of one branch)."""
    N = len(par); is_tip = [not kids[v] for v in range(N)]
    O = frozenset(v for v in range(N) if is_tip[v] and name[v].startswith(tuple(prefixes))); allt = frozenset(v for v in range(N) if is_tip[v])
    assert len(O) == len(prefixes), f"expected {len(prefixes)} outgroup tips, found {len(O)}"
    pre = [0]
    for v in pre: pre.extend(kids[v])
    clade = {}
    for v in reversed(pre): clade[v] = frozenset([v]) if is_tip[v] else frozenset().union(*[clade[k] for k in kids[v]])
    w = next((v for v in pre[1:] if clade[v] == O), None)
    if w is None: w = next((v for v in pre[1:] if clade[v] == allt - O), None)
    assert w is not None, "outgroup tips do not form a clade in this tree"
    adj = {v: [] for v in range(N)}
    for v in range(1, N): adj[v].append((par[v], ln[v])); adj[par[v]].append((v, ln[v]))
    p, L = par[w], ln[w]
    adj[w] = [(x, l) for x, l in adj[w] if x != p]; adj[p] = [(x, l) for x, l in adj[p] if x != w]
    r = N; adj[r] = [(w, L / 2), (p, L / 2)]; adj[w].append((r, L / 2)); adj[p].append((r, L / 2)); names = list(name) + [""]
    for v in list(adj):                      # suppress degree-2 internal nodes other than the new root
        if v != r and len(adj[v]) == 2 and names[v] == "":
            (a, la), (b, lb) = adj[v]; adj[a] = [(x, l) for x, l in adj[a] if x != v] + [(b, la + lb)]; adj[b] = [(x, l) for x, l in adj[b] if x != v] + [(a, la + lb)]; del adj[v]
    npar, nln, nname, nkids, idx = [], [], [], [], {}
    st = [(r, -1, 0.0)]
    while st:
        v, pv, lv = st.pop(); idx[v] = len(npar); npar.append(idx[pv] if pv >= 0 else -1); nln.append(lv); nname.append(names[v]); nkids.append([])
        if pv >= 0: nkids[idx[pv]].append(idx[v])
        for x, l in adj[v]:
            if x != pv: st.append((x, v, l))
    return npar, nln, nname, nkids

OUTGROUP = ("Cystobasidium", "Pseudomicrostroma")   # the only two non-Rhodotorula tips (same outgroup as analysis/ideas/.../idea_09_phylogeny.R)
par, ln, name, kids = reroot_on_clade(par, ln, name, kids, OUTGROUP)
print(f"rerooted on outgroup clade {OUTGROUP}: root has {len(kids[0])} children; {sum(1 for v in range(len(par)) if not kids[v])} tips")
N = len(par); depth = np.zeros(N)
order = [0]; 
for v in order: order.extend(kids[v])
for v in order[1:]: depth[v] = depth[par[v]] + ln[v]
tips = [v for v in range(N) if not kids[v]]
def leaves_under(v): out = []; st = [v]; \
    [ (out.append(x) if not kids[x] else st.extend(kids[x])) for x in iter(lambda: st.pop() if st else None, None) ] ; return out
tip_ix = {v: i for i, v in enumerate(tips)}; n = len(tips)
C = np.zeros((n, n)); under = {}
for v in reversed(order):
    under[v] = [tip_ix[v]] if not kids[v] else [t for k in kids[v] for t in under[k]]
for v in order:
    ix = under[v]; C[np.ix_(ix, ix)] = depth[v]
tipname = {i: re.sub(r"\.proteins(\.fa)?$", "", name[v]) for v, i in tip_ix.items()}
name2i = {tipname[i]: i for i in range(n)}
print(f"tree: {n} tips, {N-n} internal nodes, max depth {max(depth[tips]):.3f}")

def lam_ll(y, Cs, lam):
    V = lam * Cs + (1 - lam) * np.diag(np.diag(Cs)); k = len(y); Vi = np.linalg.inv(V); one = np.ones(k)
    mu = (one @ Vi @ y) / (one @ Vi @ one); r = y - mu; s2 = (r @ Vi @ r) / k
    return -0.5 * (k * np.log(2 * np.pi * s2) + np.linalg.slogdet(V)[1] + k)
def pagel(y, Cs):
    grid = np.linspace(0, 1, 41); ll = np.array([lam_ll(y, Cs, l) for l in grid]); j = ll.argmax()
    fine = np.linspace(grid[max(j - 1, 0)], grid[min(j + 1, 40)], 41); llf = np.array([lam_ll(y, Cs, l) for l in fine]); jj = llf.argmax()
    from scipy.stats import chi2
    return fine[jj], llf[jj], chi2.sf(2 * (llf[jj] - lam_ll(y, Cs, 0.0)), 1), chi2.sf(2 * (llf[jj] - lam_ll(y, Cs, 1.0)), 1)
def blomberg(y, Vi):
    k = len(y); one = np.ones(k); a = (one @ Vi @ y) / (one @ Vi @ one); r = y - a
    mse0 = r @ r / (k - 1); mse = r @ Vi @ r / (k - 1); return mse0 / mse, (one @ Vi @ one)
def K_perm(y, Cs, B=999):
    Vi = np.linalg.inv(Cs); k = len(y); K0, s1 = blomberg(y, Vi); exp_ = (np.trace(Cs) - k / s1) / (k - 1); K0 /= exp_
    cnt = sum(blomberg(rng.permutation(y), Vi)[0] / exp_ >= K0 for _ in range(B)); return K0, (cnt + 1) / (B + 1)

w = load_wells(); rows = []; traits = {}
for m in METALS:
    g = w[w.Metal == m].copy(); g["a_adj"], sl = size_adjusted(g); dm = top_dose(g)
    b = g[g.conc == 0].groupby("strain_id").agg(a=("a", "mean"), a_adj=("a_adj", "mean"), lnA=("lnA", "mean"), tip=("tip", "first"), species=("species", "first"))
    t = g[g.conc == dm].groupby("strain_id").a.mean(); b["delta"] = t.reindex(b.index) - b.a
    traits[m] = (b, dm)
    for scope in ("all_strains_with_tip", "R_mucilaginosa_only"):
        s = b[b.tip.notna() & b.tip.isin(name2i)]
        if scope.startswith("R_"): s = s[s.species == "Rhodotorula mucilaginosa"]
        ix = np.array([name2i[t_] for t_ in s.tip]); Cs = C[np.ix_(ix, ix)] + 1e-6 * np.mean(np.diag(C)) * np.eye(len(ix))  # jitter: identical sister tips make C singular
        for tr, lab in (("a", "baseline a*"), ("a_adj", "baseline a* (size-adjusted)"), ("delta", f"change in a* at dose {dm:g}"), ("lnA", "baseline ln area")):
            ok = s[tr].notna().values
            if ok.sum() < 30: continue
            y = s[tr].values[ok]; Cc = Cs[np.ix_(ok, ok)]
            lam, ll, p0, p1 = pagel(y, Cc); K, pK = K_perm(y, Cc)
            rows.append(dict(Metal=m, scope=scope, trait=lab, n_strains=int(ok.sum()), pagel_lambda=lam, p_lambda_gt_0=p0, p_lambda_lt_1=p1, blomberg_K=K, p_K=pK))
            print(rows[-1], flush=True)
tab = pd.DataFrame(rows); tab["p_K_BH"] = np.nan
for sc in tab.scope.unique(): i = tab.scope == sc; tab.loc[i, "p_K_BH"] = stats.false_discovery_control(tab.loc[i, "p_K"]) if hasattr(stats, "false_discovery_control") else np.nan
tab.to_csv(REP / "tables/s2_phylo_signal.csv", index=False)
# --- figure 1: lambda and K
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
for k, (col, ttl) in enumerate((("pagel_lambda", "Pagel's lambda (0 = no signal, 1 = Brownian)"), ("blomberg_K", "Blomberg's K (1 = Brownian expectation)"))):
    sub = tab[tab.trait.str.startswith(("baseline a*", "change"))].copy(); sub["lab"] = sub.Metal.str[:2] + " " + sub.trait.str.replace("baseline a* (size-adjusted)", "base adj").str.replace("baseline a*", "base").str.replace(r"change in a\* at dose .*", "delta", regex=True)
    for j, sc in enumerate(("all_strains_with_tip", "R_mucilaginosa_only")):
        s = sub[sub.scope == sc]; x = np.arange(len(s)) + (j - .5) * .38
        ax[k].bar(x, s[col], width=.36, color=("#0072B2", "#D55E00")[j], label=sc.replace("_", " "))
    ax[k].set_xticks(range(len(s))); ax[k].set_xticklabels(s.lab, rotation=75, fontsize=7); ax[k].set_title(ttl, fontsize=9); ax[k].legend(fontsize=7)
fig.tight_layout(); fig.savefig(REP / "figures/s2_phylo_signal.png", dpi=170); plt.close(fig)
# --- figure 2: tree (rooted on Cystobasidium) with a species strip and baseline a* strips for each metal
ypos = {}
dfs_ = []; st_ = [0]   # depth-first order keeps every clade contiguous (a breadth-first order scatters relatives and makes the branches cross)
while st_:
    v_ = st_.pop(); dfs_.append(v_); st_.extend(reversed(kids[v_]))
leaf_order = [v for v in dfs_ if not kids[v]]
for i_, v in enumerate(leaf_order): ypos[v] = i_
for v in reversed(order):
    if kids[v]: ypos[v] = np.mean([ypos[k] for k in kids[v]])
nL = len(leaf_order); xmax = max(depth[tips])
def species_of(tn):
    m_ = re.match(r"(Rhodotorula_sp\._clade_[IVX]+|Rhodotorula_[a-z]+|[A-Z][a-z]+_[a-z]+)", tn); return m_.group(1).replace("_", " ") if m_ else "other"
sp_all = [species_of(tipname[tip_ix[v]]) for v in leaf_order]; sp_top = [s_ for s_, c_ in pd.Series(sp_all).value_counts().items() if c_ >= 5]
PALS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#999999", "#8C564B", "#17BECF"]
spcol = {s_: PALS[i_ % len(PALS)] for i_, s_ in enumerate(sp_top)}
spimg = np.array([[matplotlib.colors.to_rgb(spcol.get(s_, "#ffffff"))] for s_ in sp_all])
MM = np.full((nL, len(METALS)), np.nan)
for mi, m in enumerate(METALS):
    b, dm = traits[m]; vals = b.set_index("tip").a
    for ri, v in enumerate(leaf_order):
        t_ = tipname[tip_ix[v]]
        if t_ in vals.index and not np.isnan(vals[t_]): MM[ri, mi] = vals[t_]
fig = plt.figure(figsize=(12, 17)); gs_ = fig.add_gridspec(1, 4, width_ratios=[6, 0.35, 2.6, 0.22], wspace=0.03)
axt = fig.add_subplot(gs_[0]); axp = fig.add_subplot(gs_[1], sharey=axt); axs = fig.add_subplot(gs_[2], sharey=axt); axc = fig.add_subplot(gs_[3])
for v in order[1:]: axt.plot([depth[par[v]], depth[v]], [ypos[v], ypos[v]], "k-", lw=.45)
for v in order:
    if kids[v]: ys = [ypos[k] for k in kids[v]]; axt.plot([depth[v], depth[v]], [min(ys), max(ys)], "k-", lw=.45)
og_ = [v for v in leaf_order if tipname[tip_ix[v]].startswith(OUTGROUP)]
axt.annotate("outgroup: Cystobasidium sp. DBVPG_10075\n+ Pseudomicrostroma phylloplanum DBVPG_6740", (depth[og_[0]], ypos[og_[0]]), xytext=(0.5 * xmax, nL - 45), fontsize=9, arrowprops=dict(arrowstyle="->", lw=.8))
axt.set_ylim(nL - 0.5, -0.5); axt.set_xlim(-0.01 * xmax, xmax * 1.02); axt.set_yticks([]); axt.set_xlabel("branch length (substitutions/site)", fontsize=10)
axp.imshow(spimg, aspect="auto", interpolation="nearest", extent=(-0.5, 0.5, nL - 0.5, -0.5)); axp.set_xticks([0]); axp.set_xticklabels(["species"], rotation=45, ha="left", fontsize=11); axp.xaxis.tick_top(); axp.set_yticks([])
im = axs.imshow(MM, aspect="auto", cmap="viridis", vmin=0, vmax=30, interpolation="nearest", extent=(-0.5, len(METALS) - 0.5, nL - 0.5, -0.5))
axs.set_xticks(range(len(METALS))); axs.set_xticklabels(METALS, rotation=45, ha="left", fontsize=11); axs.xaxis.tick_top(); axs.set_yticks([]); plt.setp(axs.get_yticklabels(), visible=False)
fig.colorbar(im, cax=axc); axc.set_ylabel("baseline a* (strain mean, dose 0); blank = no data", fontsize=10)
axt.legend(handles=[matplotlib.patches.Patch(color=spcol[s_], label=s_) for s_ in sp_top], loc="center right", fontsize=8, title="species (>= 5 tips)", title_fontsize=8, frameon=True)
fig.suptitle("PHYling tree (278 tips), rooted on the outgroup clade Cystobasidium + Pseudomicrostroma, with baseline a* per metal", fontsize=11, y=0.995)
fig.savefig(REP / "figures/s2_tree_with_astar.png", dpi=130, bbox_inches="tight"); plt.close(fig)
# --- figure 3: trait distance decay
fig, ax = plt.subplots(1, 4, figsize=(18, 4))
for k, m in enumerate(["Chromium", "Copper", "Iron", "Lead"]):
    b, dm = traits[m]; s = b[b.tip.notna() & b.tip.isin(name2i)]; ix = np.array([name2i[t] for t in s.tip]); y = s.a_adj.values
    D = np.add.outer(np.diag(C)[ix], np.diag(C)[ix]) - 2 * C[np.ix_(ix, ix)]; iu = np.triu_indices(len(ix), 1)
    d = D[iu]; dy = np.abs(y[:, None] - y[None, :])[iu]; q = np.quantile(d, np.linspace(0, 1, 13)); bi = np.digitize(d, q[1:-1])
    mean = [dy[bi == j].mean() for j in range(12)]; mid = [np.median(d[bi == j]) for j in range(12)]
    ax[k].plot(mid, mean, "-o", color=COL[m]); ax[k].set_xlabel("patristic distance (quantile bins)"); ax[k].set_ylabel("mean |difference in size-adjusted baseline a*|"); ax[k].set_title(f"{m} ({len(ix)} strains)", fontsize=9)
fig.tight_layout(); fig.savefig(REP / "figures/s2_distance_decay.png", dpi=170); plt.close(fig)
print("done")
