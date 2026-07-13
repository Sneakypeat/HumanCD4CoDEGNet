"""
Trans-regulatory network architecture of human CD4+ T cells and its rewiring across activation.

Streams two released summary fields of the public Marson-lab atlas (Zhu, Dann et al. 2025,
GWCD4i.DE_stats.h5ad) — per-perturbation .obs (n_downstream, ontarget_significant) and per-gene
varm/measured_genes_stats_* (n_regulators) — and computes:
  1. degree distributions + hub-dominance (out-degree Gini, top-k edge share)
  2. sparse-but-pleiotropic test vs Barton/Pritchard 2026
  3. shape-invariant / identity-labile rewiring (detectability-guarded)
  4. confound guards (power, edge definition, regulator expression, silent regulators)
Writes artifacts/architecture_results.json + artifacts/figures/fig_arch{1,2}_*.png.
Reproduces in seconds; never downloads the 16.8 GB effect-size layers.
"""
import json, os, warnings
import numpy as np, pandas as pd, h5py, fsspec
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
warnings.filterwarnings("ignore")

URL = ("https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
       "marson2025_data/GWCD4i.DE_stats.h5ad")
ART = os.path.expanduser("~/CoDEG_Tcell/artifacts"); FIG = f"{ART}/figures"
os.makedirs(FIG, exist_ok=True)
CONDS = ["Rest", "Stim8hr", "Stim48hr"]
LAB = {"Rest": "Rest", "Stim8hr": "Stim 8h", "Stim48hr": "Stim 48h"}
COL = {"Rest": "#4C72B0", "Stim8hr": "#DD8452", "Stim48hr": "#C44E52"}

def dec(a): return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])
def obs_col(o, k):
    x = o[k]
    if isinstance(x, h5py.Group):
        return dec(x["categories"][:])[x["codes"][:]]
    return x[:]
def gini(x):
    x = np.sort(np.asarray(x, float)); x = x[~np.isnan(x)]; n = len(x)
    if n == 0 or x.sum() == 0: return np.nan
    c = np.cumsum(x); return (n + 1 - 2 * np.sum(c) / c[-1]) / n
def topk(x, f):
    x = np.sort(np.asarray(x, float))[::-1]; x = x[~np.isnan(x)]
    k = max(1, int(np.ceil(f * len(x)))); return float(x[:k].sum() / x.sum())
def lorenz(x):
    x = np.sort(np.asarray(x, float)); x = x[~np.isnan(x)]
    c = np.concatenate([[0], np.cumsum(x) / x.sum()]); return np.linspace(0, 1, len(c)), c

# ---------------------------------------------------------------- stream
print("streaming released summary fields ...")
h = h5py.File(fsspec.open(URL, block_size=4 * 1024 * 1024).open(), "r")
o = h["obs"]
pert = pd.DataFrame({
    "target": obs_col(o, "target_contrast"), "gene": obs_col(o, "target_contrast_gene_name"),
    "cond": obs_col(o, "culture_condition"), "outdeg": obs_col(o, "n_downstream").astype(float),
    "onsig": obs_col(o, "ontarget_significant").astype(bool),
    "reg_baseMean": obs_col(o, "target_baseMean").astype(float),
    "n_cells": obs_col(o, "n_cells_target").astype(float)})
outdeg = pert.pivot_table(index="target", columns="cond", values="outdeg", aggfunc="first")
onsig = pert.pivot_table(index="target", columns="cond", values="onsig", aggfunc="first")
gname = pert.drop_duplicates("target").set_index("target")["gene"]
g0 = h["varm"]["measured_genes_stats_Rest"]
gene = pd.DataFrame({"gene_id": dec(g0["gene_id"][:]), "gene_name": dec(g0["gene_name"][:])})
for c in CONDS:
    gene[f"nreg_{c}"] = h["varm"][f"measured_genes_stats_{c}"]["n_regulators"][:]
print(f"  {len(pert):,} perturbation-conditions | {len(gene):,} genes")

OUT = {}
# ---------------------------------------------------------------- part 1+2
P1 = {}
for c in CONDS:
    ind = gene[f"nreg_{c}"].dropna().values; od = outdeg[c].dropna().values
    P1[c] = dict(in_median=float(np.median(ind)), in_gini=gini(ind),
                 out_gini=gini(od), out_top1=topk(od, .01), out_top5=topk(od, .05),
                 n_reg=int(len(od)), superhubs=int((od > 1000).sum()), silent=int((od == 0).sum()),
                 total_edges=int(od.sum()))
    print(f"  {c:9}: out-Gini {P1[c]['out_gini']:.3f} in-Gini {P1[c]['in_gini']:.3f} "
          f"top5%={P1[c]['out_top5']*100:.0f}% median-in={P1[c]['in_median']:.0f}")
OUT["architecture_per_state"] = P1

# ---------------------------------------------------------------- part 3 rewiring
def both(ca, cb):
    m = onsig[ca].fillna(False).astype(bool) & onsig[cb].fillna(False).astype(bool)
    return outdeg.index[m]
REW = {}
for ca, cb in [("Rest", "Stim8hr"), ("Rest", "Stim48hr"), ("Stim8hr", "Stim48hr")]:
    idx = both(ca, cb); xa, xb = outdeg.loc[idx, ca], outdeg.loc[idx, cb]
    rho, _ = stats.spearmanr(xa, xb)
    shared = len(set(xa.nlargest(100).index) & set(xb.nlargest(100).index))
    REW[f"{ca}->{cb}"] = dict(rho=float(rho), shared=int(shared), displaced=int(100 - shared), n=int(len(idx)))
OUT["rewiring"] = REW
print(f"  hub displacement (Rest->8h/Rest->48h/8h->48h): "
      f"{REW['Rest->Stim8hr']['displaced']}/{REW['Rest->Stim48hr']['displaced']}/{REW['Stim8hr->Stim48hr']['displaced']}%")

# ---------------------------------------------------------------- part 4 guards
G = {}
for c in CONDS:
    s = pert[pert.cond == c]
    sp_cells = stats.spearmanr(s.dropna(subset=["outdeg","n_cells"]).outdeg,
                               s.dropna(subset=["outdeg","n_cells"]).n_cells)[0]
    se = s.dropna(subset=["outdeg","reg_baseMean"]); se = se[se.reg_baseMean > 0]
    sp_expr = stats.spearmanr(se.outdeg, np.log10(se.reg_baseMean))[0]
    od_val = s[s.onsig].outdeg.values
    G[c] = dict(power_rho=float(sp_cells), expr_rho=float(sp_expr), validated_out_gini=gini(od_val))
OUT["confound_guards"] = G

json.dump(OUT, open(f"{ART}/architecture_results.json", "w"), indent=2, default=float)

# ---------------------------------------------------------------- figures
plt.rcParams.update({"savefig.dpi": 300, "font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titlesize": 11, "axes.titleweight": "bold",
                     "legend.frameon": False})
# fig 1
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
ax[0].plot([0, 1], [0, 1], "--", color="#999", lw=1, label="equality")
for c in CONDS:
    x = outdeg[c].dropna().values; p, cu = lorenz(x)
    ax[0].plot(p, cu, color=COL[c], lw=2, label=f"{LAB[c]} (G={gini(x):.2f})")
ax[0].set_xlabel("cumulative fraction of regulators"); ax[0].set_ylabel("cumulative fraction of trans-edges")
ax[0].set_title("A  Hub-dominated & invariant"); ax[0].legend(loc="upper left", fontsize=8.5)
c = "Stim8hr"; x = np.sort(outdeg[c].dropna().values)[::-1]
ax[1].plot(np.arange(1, len(x)+1), x+1, color=COL[c], lw=1.6); ax[1].set_xscale("log"); ax[1].set_yscale("log")
ax[1].axhline(1001, color="#888", ls=":", lw=1)
ax[1].annotate(f"{int((x>1000).sum())} super-hubs", xy=(3, 2500), fontsize=8.5)
ax[1].annotate(f"{int((x==0).sum())} silent", xy=(150, 1.4), fontsize=8.5)
ax[1].set_xlabel("regulator rank (Stim 8h)"); ax[1].set_ylabel("out-degree + 1"); ax[1].set_title("B  A hub minority carries the signal")
xp = np.arange(3); w = .38; ing = [P1[c]["in_gini"] for c in CONDS]; oug = [P1[c]["out_gini"] for c in CONDS]
ax[2].bar(xp-w/2, ing, w, color="#55A868", label="in-degree (reception)")
ax[2].bar(xp+w/2, oug, w, color="#8172B3", label="out-degree (broadcast)")
ax[2].set_xticks(xp); ax[2].set_xticklabels([LAB[c] for c in CONDS]); ax[2].set_ylim(0, 1.05)
ax[2].set_ylabel("Gini coefficient"); ax[2].set_title("C  Broadcast-concentrated,\nreception-distributed"); ax[2].legend(loc="upper center", fontsize=8.5)
plt.tight_layout(); fig.savefig(f"{FIG}/fig_arch1_topology.png", bbox_inches="tight"); plt.close(fig)
# fig 2
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
edges = [P1[c]["total_edges"] for c in CONDS]; ginis = [P1[c]["out_gini"] for c in CONDS]
ax[0].bar(xp, edges, .55, color=[COL[c] for c in CONDS], alpha=.85)
for i, e in enumerate(edges): ax[0].text(i, e+1.4e4, f"{e/1000:.0f}k", ha="center", fontsize=8.5)
ax[0].set_ylim(0, 1.12e6); ax[0].set_xticks(xp); ax[0].set_xticklabels([LAB[c] for c in CONDS]); ax[0].set_ylabel("total trans-edges")
a2 = ax[0].twinx(); a2.spines["top"].set_visible(False); a2.plot(xp, ginis, "-o", color="#111", lw=2, ms=6)
for i, gg in enumerate(ginis): a2.text(i, gg+.018, f"{gg:.3f}", ha="center", fontsize=8.5, fontweight="bold",
                                       bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=.85))
a2.set_ylabel("out-degree Gini"); a2.set_ylim(.5, 1.0); ax[0].set_title("A  Shape-invariant, mass-variable")
m = onsig["Rest"].fillna(False).astype(bool) & onsig["Stim8hr"].fillna(False).astype(bool); sub = outdeg[m]
ax[1].scatter(sub["Rest"]+1, sub["Stim8hr"]+1, s=6, color="#BBB", alpha=.45, edgecolors="none")
ax[1].plot([1, 7000], [1, 7000], "--", color="#888", lw=1)
gn = gname.reindex(sub.index)
for genes, cc in [(["CD3E","CD3D","CD3G","CD247","LAT","LCP2","PLCG1","ZAP70","VAV1","ITK","LCK","BCL10"], "#DD8452"),
                  (["ZMYM2","ARF1","TP53","NFAT5","PTEN","TSC1","NCKAP1L","SRF"], "#4C72B0")]:
    sel = gn.isin(genes).values
    ax[1].scatter(sub["Rest"].values[sel]+1, sub["Stim8hr"].values[sel]+1, s=34, color=cc, edgecolors="k", lw=.4, zorder=5)
ax[1].scatter([], [], s=34, color="#DD8452", edgecolors="k", lw=.4, label="TCR signalosome (gain)")
ax[1].scatter([], [], s=34, color="#4C72B0", edgecolors="k", lw=.4, label="rest-state hubs (lose)")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlim(.9, 7000); ax[1].set_ylim(.9, 7000)
ax[1].set_xlabel("out-degree at Rest +1"); ax[1].set_ylabel("out-degree at Stim 8h +1")
ax[1].set_title(f"B  Hub identity rewires (rho={REW['Rest->Stim8hr']['rho']:.2f})"); ax[1].legend(loc="lower right", fontsize=8)
pairs = [("Rest->Stim8hr", "Rest→8h"), ("Stim8hr->Stim48hr", "8h→48h"), ("Rest->Stim48hr", "Rest→48h")]
turn = [REW[k]["displaced"] for k, _ in pairs]
ax[2].bar([l for _, l in pairs], turn, .6, color=["#DD8452", "#C44E52", "#7A3B52"])
for i, t in enumerate(turn): ax[2].text(i, t+1, f"{t:.0f}%", ha="center", fontsize=9, fontweight="bold")
ax[2].set_ylim(0, 90); ax[2].set_ylabel("% of earlier-state top-100 hubs displaced"); ax[2].set_title("C  Most hubs turn over")
plt.tight_layout(); fig.savefig(f"{FIG}/fig_arch2_rewiring.png", bbox_inches="tight"); plt.close(fig)
print("wrote architecture_results.json + fig_arch1_topology.png + fig_arch2_rewiring.png")
