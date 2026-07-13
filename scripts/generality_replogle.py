"""
GENERALITY TEST — is the hub-dominated / broadcast-concentrated architecture T-cell-specific?
Downloads Replogle 2022 genome-scale Perturb-seq DE gene-sets (K562 genome-wide, K562 essential,
RPE1 essential) from Harmonizome as directed edge lists (perturbation -> DE gene, same concept as
Marson n_downstream) and recomputes the architecture. Also tests, across two cell types (K562 vs RPE1),
whether the SHAPE (out-degree Gini) is conserved while hub IDENTITY turns over — a cross-context analog
of the within-T-cell activation result.

Honest caveats baked in: (1) Harmonizome thresholding != Marson 10% FDR, so ABSOLUTE Gini values are not
directly comparable — only the QUALITATIVE shape and the within-Replogle K562<->RPE1 comparison are.
(2) The essential-only subsets exclude the silent-majority tail that drives hub-dominance, so the
genome-scale K562 network is the fair comparator to Marson; the essential subsets compress the asymmetry.
"""
import os, re, json, numpy as np, pandas as pd, requests
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ART = os.path.expanduser("~/CoDEG_Tcell/artifacts"); FIG = f"{ART}/figures"
CACHE = os.path.expanduser("~/CoDEG_Tcell/data/replogle"); os.makedirs(CACHE, exist_ok=True); os.makedirs(FIG, exist_ok=True)
BASE = "https://maayanlab.cloud/static/hdfs/harmonizome/data"
NET = {"K562\ngenome-wide": "reploglek562genomewide", "K562\nessential": "reploglek562essential",
       "RPE1\nessential": "reploglerpe1essential"}
MARSON_OUT_GINI = 0.929   # Marson CD4+ T, Stim 8h, 10% FDR (see architecture_results.json)
MARSON_IN_GINI = 0.354

def gini(x):
    x = np.sort(np.asarray(x, float)); x = x[~np.isnan(x)]; n = len(x)
    if n == 0 or x.sum() == 0: return np.nan
    c = np.cumsum(x); return (n + 1 - 2 * np.sum(c) / c[-1]) / n
def topk(x, f):
    x = np.sort(np.asarray(x, float))[::-1]; k = max(1, int(np.ceil(f * len(x)))); return float(x[:k].sum() / x.sum())
def load(slug):
    fp = f"{CACHE}/{slug}.edges.txt.gz"
    if not os.path.exists(fp):
        r = requests.get(f"{BASE}/{slug}/gene_attribute_edges.txt.gz", timeout=120); r.raise_for_status()
        open(fp, "wb").write(r.content)
    d = pd.read_csv(fp, sep="\t", compression="gzip")
    d.columns = [c.strip() for c in d.columns]
    return d.rename(columns={"Gene": "gene", "Gene Perturbation": "pert"})[["gene", "pert"]].dropna()
def pgene(p):
    m = re.match(r"^\d+_(.+?)_P\d", p); return m.group(1) if m else p

OUT = {}; data = {}
print(f"{'network':16} {'n_pert':>6} {'edges':>7} | {'out-Gini':>8} {'out-top5%':>9} | {'in-Gini':>7} | {'asym':>6}")
for name, slug in NET.items():
    d = load(slug); data[name] = d
    od = d.groupby("pert").size().values; ind = d.groupby("gene").size().values
    r = dict(n_pert=int(d.pert.nunique()), edges=int(len(d)), out_gini=gini(od), out_top5=topk(od, .05),
             in_gini=gini(ind), out_med=float(np.median(od)), in_med=float(np.median(ind)))
    r["asym"] = r["out_gini"] - r["in_gini"]; OUT[name.replace("\n", "-")] = r
    print(f"{name.replace(chr(10),' '):16} {r['n_pert']:>6} {r['edges']:>7} | {r['out_gini']:>8.3f} {r['out_top5']*100:>8.0f}% | {r['in_gini']:>7.3f} | {r['asym']:>+6.3f}")

# cross-cell-type: K562-ess vs RPE1-ess, matched on perturbed-gene symbol
def by_gene(d):
    s = d.groupby("pert").size(); g = pd.Series({p: pgene(p) for p in s.index})
    return pd.DataFrame({"gene": g.values, "od": s.values}).groupby("gene").od.max()
kk, rr = by_gene(data["K562\nessential"]), by_gene(data["RPE1\nessential"])
shared = sorted(set(kk.index) & set(rr.index)); xa, xb = kk.loc[shared].values, rr.loc[shared].values
rho, _ = stats.spearmanr(xa, xb)
sh = pd.DataFrame({"K562": xa, "RPE1": xb}, index=shared)
ha, hb = set(sh.K562.nlargest(100).index), set(sh.RPE1.nlargest(100).index)
OUT["cross_cell"] = dict(n_matched=len(shared), spearman_rho=float(rho),
                         top100_displaced=int(100 - len(ha & hb)), top100_shared=int(len(ha & hb)))
OUT["marson_reference"] = dict(out_gini=MARSON_OUT_GINI, in_gini=MARSON_IN_GINI)
json.dump(OUT, open(f"{ART}/architecture_replogle_results.json", "w"), indent=2, default=float)
print(f"\ncross-cell K562<->RPE1: {len(shared)} matched | out-Gini {OUT['K562-essential']['out_gini']:.2f} vs "
      f"{OUT['RPE1-essential']['out_gini']:.2f} | Spearman {rho:+.2f} | top-100 displaced {100-len(ha&hb)}%")

# --------------------------------------------------------------- figure
plt.rcParams.update({"savefig.dpi": 300, "font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titlesize": 11, "axes.titleweight": "bold", "legend.frameon": False})
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
labels = ["CD4+ T\n(Marson)"] + list(NET.keys())
og = [MARSON_OUT_GINI] + [OUT[n.replace("\n", "-")]["out_gini"] for n in NET]
ig = [MARSON_IN_GINI] + [OUT[n.replace("\n", "-")]["in_gini"] for n in NET]
cols = ["#2A6F97", "#468FAF", "#9AC0CD", "#C6D8DF"]  # T + K562gw fair (dark), essential subsets (light)
# A: hub-dominance across systems
ax[0].bar(range(4), og, .62, color=cols)
for i, v in enumerate(og): ax[0].text(i, v + .015, f"{v:.2f}", ha="center", fontsize=8.5, fontweight="bold")
ax[0].set_ylim(0, 1.05); ax[0].set_xticks(range(4)); ax[0].set_xticklabels(labels, fontsize=8.5)
ax[0].set_ylabel("out-degree Gini (hub-dominance)")
ax[0].set_title("A  Hub-dominance generalizes\n(genome-scale: T & K562 both >0.85)")
ax[0].axhline(0.85, color="#999", ls=":", lw=1); ax[0].text(3.4, 0.865, "0.85", fontsize=7.5, color="#777", ha="right")
# B: broadcast vs reception asymmetry
xp = np.arange(4); w = .38
ax[1].bar(xp - w/2, ig, w, color="#55A868", label="in-degree (reception)")
ax[1].bar(xp + w/2, og, w, color="#8172B3", label="out-degree (broadcast)")
ax[1].set_xticks(xp); ax[1].set_xticklabels(labels, fontsize=8.5); ax[1].set_ylim(0, 1.0)
ax[1].set_ylabel("Gini coefficient"); ax[1].legend(loc="upper right", fontsize=8)
ax[1].set_title("B  Broadcast > reception\n(genome-scale systems)")
# C: cross-cell-type scatter
ax[2].scatter(xa + 1, xb + 1, s=8, color="#BBB", alpha=.5, edgecolors="none")
lim = [1, max(xa.max(), xb.max()) + 1]
ax[2].plot(lim, lim, "--", color="#888", lw=1)
ax[2].set_xscale("log"); ax[2].set_yscale("log")
ax[2].set_xlabel("out-degree in K562 (essential) +1"); ax[2].set_ylabel("out-degree in RPE1 (essential) +1")
ax[2].set_title(f"C  Cross-cell-type: shape conserved,\nidentity turns over ({100-len(ha&hb)}% of top-100)")
ax[2].annotate(f"ρ={rho:.2f}\n{len(shared)} matched essentials", xy=(.05, .88), xycoords="axes fraction", fontsize=8.5)
plt.tight_layout(); fig.savefig(f"{FIG}/fig_arch3_generality.png", bbox_inches="tight"); plt.close(fig)
print("wrote fig_arch3_generality.png + architecture_replogle_results.json")
