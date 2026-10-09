#!/usr/bin/env python3
"""Assemble analysis/gwas_dh4148/report/REPORT.md from scripts/report_text.md ({{T_*}} = tables from report/tables, {{IMG:file|caption}} = figures)."""
import re
from pathlib import Path
import numpy as np, pandas as pd
REP = Path("analysis/gwas_dh4148/report"); T = REP / "tables"
def md(df, fmt="{:.2f}", pcols=("p", "p_wald", "p_interaction", "anova_p", "p_epistatic", "min_p", "min_p_lmm")):
    df = df.copy(); df.columns = [str(c) for c in df.columns]; cols = list(df.columns)
    for c in cols:
        if df[c].dtype.kind == "f": df[c] = df[c].map(lambda v: "" if pd.isna(v) else (f"{v:.1e}" if (c in pcols or c.startswith("p_")) and v < 0.001 else (f"{v:.3f}" if c in pcols or c.startswith("p_") else fmt.format(v))))
    esc = lambda x: str(x).replace("|", "\\|").replace("*", "\\*")
    rows = [[esc(x) for x in r] for r in df.itertuples(index=False)]
    wd = [min(40, max(len(c), *(len(r[i]) for r in rows))) if rows else len(c) for i, c in enumerate(cols)]
    out = ["| " + " | ".join(esc(c) for c in cols) + " |", "|" + "|".join("-" * max(3, w) for w in wd) + "|"] + ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)
img = lambda f, cap: "![](figures/%s.png)\n\n*%s*\n" % (f, cap.replace("*", "\\*"))
body = open("analysis/gwas_dh4148/scripts/report_text.md").read()
exec(open("analysis/gwas_dh4148/scripts/report_tables.py").read())   # defines dict TAB
for k, v in TAB.items(): body = body.replace("{{" + k + "}}", v)
for f in sorted(set(re.findall(r"\{\{IMG:([A-Za-z0-9_]+)\|([^}]*)\}\}", body))): body = body.replace("{{IMG:%s|%s}}" % f, img(f[0], f[1]))
assert "{{" not in body, re.findall(r"\{\{[^}]*\}\}", body)[:5]
_l = body.split("\n"); _o = []
for i, l in enumerate(_l):
    if l.startswith("- ") and _o and _o[-1].strip() and not _o[-1].startswith(("- ", "  ", "|")): _o.append("")
    _o.append(l)
(REP / "REPORT.md").write_text("\n".join(_o)); print("REPORT.md written")
