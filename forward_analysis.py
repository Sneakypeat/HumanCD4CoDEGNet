"""Forward / impact demonstrations: (A) buffering is dynamic across T-cell activation;
(B) buffering as a transcriptional-tractability axis for target selection (reconnects to the
track's drug-target prompt). Ships the per-gene buffering-score resource.
"""
import json, io, re, requests
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

ART = Path.home() / "CoDEG_Tcell" / "artifacts"
FIG = ART / "figures"; FIG.mkdir(exist_ok=True)
ESS, NON, ACC, WARM = "#c1121f", "#457b9d", "#2a9d8f", "#e07a5f"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})

df = pd.read_csv(ART / "per_gene_full.csv")
res = pd.read_csv(ART / "results_table.csv")
out = {}

# ---------------- A. Dynamic buffering across activation ----------------
# per-condition matched Cliff's delta (already computed in results_table.csv)
conds = ["Rest", "Stim8hr", "Stim48hr"]
axis_of = {"Hart": "Hart_vs_rest", "CEGv2": "CEGv2_vs_rest"}
dyn = {a: [float(res[(res.axis == axis_of[a]) & (res.scope == c)].matched_delta.values[0]) for c in conds]
       for a in ["Hart", "CEGv2"]}
out["dynamic_buffering"] = {"conditions": conds, **dyn}

# genes most buffered SPECIFICALLY at early activation (resid drop Rest -> Stim8hr)
df["activation_buffering"] = df["resid_Stim8hr"] - df["resid_Rest"]  # more negative = more buffered on activation
top_act = df.nsmallest(15, "activation_buffering")[["gene_name", "resid_Rest", "resid_Stim8hr",
                                                     "activation_buffering", "ess_hart", "ess_cegv2"]]
out["top_activation_buffered_genes"] = top_act.gene_name.tolist()

fig, ax = plt.subplots(figsize=(6.6, 4.2))
x = np.arange(3)
ax.plot(x, dyn["Hart"], "-o", color=ESS, lw=2.4, ms=8, label="Hart essentials")
ax.plot(x, dyn["CEGv2"], "-o", color=WARM, lw=2.4, ms=8, label="CEGv2 essentials")
ax.axhline(0, color="#999", lw=1, ls="--")
ax.set_xticks(x); ax.set_xticklabels(["Rest", "Stim 8 h", "Stim 48 h"])
ax.set_ylabel("essential-gene buffering\n(matched Cliff's $\\delta$; more negative = more buffered)")
ax.set_title("A forward result: essential-gene buffering is DYNAMIC,\npeaking at early T-cell activation (8 h)",
             fontsize=11, loc="left")
ax.legend(frameon=False, loc="lower left")
ax.invert_yaxis()  # so 'more buffered' points up
fig.savefig(FIG / "fig4_dynamic_buffering.png"); plt.close(fig)

# ---------------- B. Buffering as a transcriptional-tractability axis ----------------
base = "https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master/metadata/gene_lists"
cats = ["kinases", "ion_channels", "gpcr_union", "nuclear_receptors", "transporters", "catalytic_receptors"]
def load_symbols(name):
    txt = requests.get(f"{base}/{name}.tsv", timeout=60).text
    syms = set()
    for i, line in enumerate(txt.splitlines()):
        tok = line.split("\t")[0].strip()
        if not tok or (i == 0 and tok.lower() in {"gene", "symbol", "hgnc_symbol", "x", "gene_name"}):
            continue
        if re.match(r"^[A-Z0-9][A-Z0-9orf.-]+$", tok):
            syms.add(tok)
    return syms
drug = {c: load_symbols(c) for c in cats}
druggable = set().union(*drug.values())
df["druggable"] = df.gene_name.isin(druggable)
out["druggable_n_in_panel"] = int(df.druggable.sum())

# buffering distribution: druggable vs non (resid_total; higher = more responsive/tractable)
dd = df.loc[df.druggable, "resid_total"].values
nn = df.loc[~df.druggable, "resid_total"].values
out["druggable_vs_non_resid"] = {"median_druggable": float(np.median(dd)),
                                 "median_non": float(np.median(nn)),
                                 "mwu_p": float(stats.mannwhitneyu(dd, nn).pvalue)}

# TRACTABLE ESSENTIAL TARGETS: essential AND druggable AND least-buffered (most transcriptionally modulatable)
ess = df.ess_hart | df.ess_cegv2
cand = df.loc[ess & df.druggable].copy()
cand["buffering_score"] = -cand["resid_total"]  # higher = more buffered
cand = cand.sort_values("buffering_score")  # least buffered first = most tractable
n_cand = len(cand)
k = min(10, n_cand // 2)  # non-overlapping extremes
out["n_essential_druggable"] = int(n_cand)
out["tractable_essential_druggable"] = cand.head(k).gene_name.tolist()
out["buffered_essential_druggable_hard"] = cand.tail(k).gene_name.tolist()

fig, ax = plt.subplots(figsize=(7, 4.2))
for lab, vals, col in [("all genes", df.resid_total, "#bbb"),
                       ("druggable-genome", df.loc[df.druggable, "resid_total"], ACC),
                       ("essential", df.loc[ess, "resid_total"], ESS)]:
    ax.hist(np.clip(vals, -300, 300), bins=60, histtype="step", lw=2, density=True, label=lab, color=col)
ax.axvline(0, color="#333", lw=1)
ax.set_xlabel("expression-corrected responsiveness   ($\\leftarrow$ buffered   ·   tractable $\\rightarrow$)")
ax.set_ylabel("density")
ax.set_title("Buffering as a target-tractability axis\nessential genes skew buffered; druggable genes skew tractable",
             fontsize=10.5, loc="left")
ax.legend(frameon=False)
fig.savefig(FIG / "fig5_tractability.png"); plt.close(fig)

# ---------------- ship the resource table ----------------
score = df[["gene_id", "gene_name", "nreg_total", "resid_total", "baseMean", "shet",
            "ess_hart", "ess_cegv2", "druggable"]].copy()
score["buffering_score"] = -score["resid_total"]
score["buffering_pctile"] = score["buffering_score"].rank(pct=True).round(3)
score = score.sort_values("buffering_score", ascending=False)
score.to_csv(ART / "buffering_score_resource.csv", index=False)

(ART / "forward_results.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
print("\nwrote fig4_dynamic_buffering.png, fig5_tractability.png, buffering_score_resource.csv")
