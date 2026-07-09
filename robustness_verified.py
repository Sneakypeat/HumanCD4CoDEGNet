"""Canonical, verified robustness numbers (I recompute rather than trust agent prose).
Consolidates: binning sensitivity, NB-GLM, permutation null, family enrichment, family-exclusion.
Writes artifacts/robust_verified.json.
"""
import json, re
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm

ART = Path.home() / "CoDEG_Tcell" / "artifacts"
df = pd.read_csv(ART / "per_gene_full.csv")
rng = np.random.default_rng(20260709)

def cliffs(x, y):
    nx, ny = len(x), len(y)
    R = stats.rankdata(np.r_[x, y])[:nx].sum()
    return 2 * (R - nx * (nx + 1) / 2) / (nx * ny) - 1

def matched(d, me, mn, metric="nreg_total", nb=10, exclude=None):
    dd = d.loc[me | mn].copy()
    if exclude is not None:
        dd = dd.loc[~exclude.loc[dd.index]]
    dd["ess"] = me.loc[dd.index].values
    dd["b"] = pd.qcut(dd.baseMean, nb, labels=False, duplicates="drop")
    ds, ws, zs = [], [], []
    for _, g in dd.groupby("b"):
        e = g.loc[g.ess, metric].values; n = g.loc[~g.ess, metric].values
        if len(e) < 3 or len(n) < 3: continue
        ds.append(cliffs(e, n)); ws.append(len(e))
        _, p = stats.mannwhitneyu(e, n, alternative="less")
        zs.append(stats.norm.isf(np.clip(p, 1e-12, 1 - 1e-12)))
    z = np.sum(zs) / np.sqrt(len(zs))
    return float(np.average(ds, weights=ws)), float(stats.norm.sf(z)), len(ds)

out = {}

# 1. binning sensitivity
out["binning"] = {}
for name, me, mn in [("Hart", df.ess_hart, ~df.ess_hart), ("CEGv2", df.ess_cegv2, ~df.ess_cegv2)]:
    out["binning"][name] = {str(nb): round(matched(df, me, mn, nb=nb)[0], 3) for nb in (5, 10, 20, 50)}

# 2. NB-GLM, no binning: nreg_total ~ essential + log10(baseMean+1)
out["nb_glm"] = {}
X = sm.add_constant(pd.DataFrame({"log_bm": np.log10(df.baseMean + 1)}))
for name, col in [("Hart", "ess_hart"), ("CEGv2", "ess_cegv2")]:
    Xf = X.copy(); Xf["essential"] = df[col].astype(int)
    y = df.nreg_total.clip(lower=0)  # 14 tiny negatives -> 0 (conservative; all non-essential)
    m = sm.GLM(y, Xf, family=sm.families.NegativeBinomial()).fit()
    out["nb_glm"][name] = {"coef": round(float(m.params["essential"]), 4),
                           "p": float(m.pvalues["essential"]), "note": "negative coef = buffered"}

# 3. permutation null for Hart matched delta (shuffle essential label within baseMean deciles)
obs = matched(df, df.ess_hart, ~df.ess_hart)[0]
dfp = df.copy(); dfp["b"] = pd.qcut(dfp.baseMean, 10, labels=False, duplicates="drop")
n_perm, extreme = 2000, 0
null = np.empty(n_perm)
for i in range(n_perm):
    lab = np.zeros(len(dfp), bool)
    for _, idx in dfp.groupby("b").groups.items():
        idx = np.array(list(idx))
        k = int(df.ess_hart.iloc[idx].sum())
        if k: lab[rng.choice(idx, k, replace=False)] = True
    s = pd.Series(lab, index=dfp.index)
    null[i] = matched(dfp, s, ~s)[0]
p_perm = (np.sum(null <= obs) + 1) / (n_perm + 1)  # one-sided (obs is negative)
out["permutation"] = {"observed_delta": round(obs, 4), "n_perm": n_perm,
                      "null_mean": round(float(null.mean()), 4),
                      "null_ci95": [round(float(np.percentile(null, 2.5)), 3),
                                    round(float(np.percentile(null, 97.5)), 3)],
                      "p_one_sided": float(p_perm)}

# 4. family enrichment in most-buffered decile (resid_total ascending) + exclusion tests
fams = {"OXPHOS": r"^NDUF|^COX\d|^ATP5|^UQCR|^SDH[A-D]", "proteasome": r"^PSM[A-D]",
        "spliceosome": r"^SNRP|^PRPF|^SF3", "aaRS": r"^[A-Z]{1,3}ARS\d?$",
        "ribosomal": r"^RP[LS]\d|^MRP[LS]"}
thr = df.resid_total.quantile(0.10)
buffered = df.resid_total <= thr
out["most_buffered_decile"] = {"n": int(buffered.sum()), "resid_threshold": round(float(thr), 1)}
out["family_enrichment"] = {}
for f, rx in fams.items():
    m = df.gene_name.str.contains(rx, regex=True, na=False)
    if m.sum() == 0:
        out["family_enrichment"][f] = {"n_in_panel": 0, "note": "absent (HVG-excluded)"}
        continue
    tab = [[int((m & buffered).sum()), int((m & ~buffered).sum())],
           [int((~m & buffered).sum()), int((~m & ~buffered).sum())]]
    orr, p = stats.fisher_exact(tab, alternative="greater")
    out["family_enrichment"][f] = {"n_in_panel": int(m.sum()), "OR": round(float(orr), 2), "p": float(p)}

oxphos = df.gene_name.str.contains(fams["OXPHOS"], regex=True, na=False)
prot = df.gene_name.str.contains(fams["proteasome"], regex=True, na=False)
out["exclusion"] = {"none": {}, "OXPHOS": {}, "OXPHOS+proteasome": {}}
for lab, exc in [("none", None), ("OXPHOS", oxphos), ("OXPHOS+proteasome", oxphos | prot)]:
    for name, me, mn in [("Hart", df.ess_hart, ~df.ess_hart), ("CEGv2", df.ess_cegv2, ~df.ess_cegv2)]:
        d2, p2, _ = matched(df, me, mn, exclude=exc)
        out["exclusion"][lab][name] = {"delta": round(d2, 3), "p": float(p2),
                                       "n_excluded": int(exc.sum()) if exc is not None else 0}

(ART / "robust_verified.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
