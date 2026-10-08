"""Shared helpers for the stratified a* analyses (no side effects on import)."""
import numpy as np, pandas as pd
from scipy import stats
from pathlib import Path

R = Path("analysis/carotenoid_stress_vs_baseline/results")
REP = Path("analysis/carotenoid_stress_vs_baseline/report"); (REP / "figures").mkdir(parents=True, exist_ok=True); (REP / "tables").mkdir(parents=True, exist_ok=True)
METALS = ["Chromium", "Copper", "Iron", "Lead", "Zinc"]
COL = {"Chromium": "#0072B2", "Copper": "#D55E00", "Iron": "#009E73", "Lead": "#CC79A7", "Zinc": "#E69F00"}

def load_wells(species=True):
    w = pd.read_csv(R / "wells.csv")
    if species:
        st = pd.read_csv(R / "strat/strain_table.csv", dtype={"strain_id": str})
        w["strain_id"] = w.strain_id.astype(str)
        w = w.merge(st[["strain_id", "sample_name", "species", "pop", "tip"]], on="strain_id", how="left")
    return w

def top_dose(g, min_shared=30):
    """Highest dose that shares >= min_shared strains with dose 0 (Zinc: 15; others: their maximum)."""
    n0 = set(g[g.conc == 0].strain_id)
    return max(d for d in g.conc.unique() if d > 0 and len(n0 & set(g[g.conc == d].strain_id)) >= min_shared)

def size_adjusted(g):
    """a* residual from the dose-0 regression of a* on ln area (same metal)."""
    b = g[g.conc == 0]; sl, ic, *_ = stats.linregress(b.lnA, b.a)
    return g.a - (ic + sl * g.lnA), sl
