"""Figures for the buffering finding. Reads artifacts/*.csv, writes artifacts/figures/*.png."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec

ART = Path.home() / "CoDEG_Tcell" / "artifacts"
FIG = ART / "figures"; FIG.mkdir(exist_ok=True)
ESS, NON, ACC = "#c1121f", "#457b9d", "#2a9d8f"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})

gene = pd.read_csv(ART / "per_gene_full.csv")
res = pd.read_csv(ART / "results_table.csv")

# ================= FIGURE 1 — the sign flip =================
fig = plt.figure(figsize=(11, 4.4))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.25], wspace=0.28)

# Panel A: raw distributions (the naive, confounded answer)
axA = fig.add_subplot(gs[0])
e = gene.loc[gene.ess_hart, "nreg_total"].values
n = gene.loc[~gene.ess_hart, "nreg_total"].values
parts = axA.violinplot([np.log10(n + 1), np.log10(e + 1)], showmedians=True, widths=0.85)
for pc, col in zip(parts["bodies"], [NON, ESS]):
    pc.set_facecolor(col); pc.set_alpha(.6)
for k in ("cbars", "cmins", "cmaxes", "cmedians"):
    parts[k].set_color("#333")
axA.set_xticks([1, 2]); axA.set_xticklabels(["non-essential", "essential\n(Hart)"])
axA.set_ylabel("log10(responsiveness + 1)")
axA.set_title("A  Naive view: essential look MORE responsive\n"
              r"Cliff's $\delta$=+0.26, p=8×10$^{-7}$", fontsize=10.5, loc="left")

# Panel B: within-baseMean-decile medians (the corrected answer)
axB = fig.add_subplot(gs[1])
d = gene.copy()
d["bm_dec"] = pd.qcut(d.baseMean, 10, labels=False, duplicates="drop")
me = d[d.ess_hart].groupby("bm_dec").nreg_total.median()
mn = d[~d.ess_hart].groupby("bm_dec").nreg_total.median()
x = np.arange(10)
axB.plot(x, mn.reindex(x), "-o", color=NON, label="non-essential", lw=2, ms=5)
axB.plot(x, me.reindex(x), "-o", color=ESS, label="essential (Hart)", lw=2, ms=5)
axB.fill_between(x, me.reindex(x), mn.reindex(x),
                 where=(mn.reindex(x).values >= me.reindex(x).values),
                 color=ESS, alpha=.10, interpolate=True)
axB.set_xlabel("baseline expression decile (low → high)")
axB.set_ylabel("median responsiveness")
axB.set_title("B  Matched at equal expression: essential are BUFFERED\n"
              r"within-decile $\delta$=−0.17, p=8×10$^{-6}$ (9/10 deciles)", fontsize=10.5, loc="left")
axB.legend(frameon=False, loc="upper left")
fig.suptitle("Transcriptional buffering of essential genes is masked by an expression confound\n"
             "primary human CD4+ T cells (Marson lab genome-scale Perturb-seq)",
             fontsize=12, y=1.06, x=0.02, ha="left", weight="bold")
fig.savefig(FIG / "fig1_sign_flip.png"); plt.close(fig)

# ================= FIGURE 2 — robustness forest =================
r = res[res.scope == "pooled"].copy()
order = ["Hart_vs_rest", "CEGv2_vs_rest", "shet_continuous"]
labels = {"Hart_vs_rest": "Hart essentials (n=124)", "CEGv2_vs_rest": "CEGv2 essentials (n=377)",
          "shet_continuous": "shet constraint (n=10,166)"}
fig, ax = plt.subplots(figsize=(8.2, 4.6))
rows = []
for a in order:
    raw = res[(res.axis == a) & (res.scope == "pooled")].raw_delta.values[0]
    for scope in ["pooled", "Rest", "Stim8hr", "Stim48hr"]:
        m = res[(res.axis == a) & (res.scope == scope)].matched_delta.values[0]
        rows.append((f"{labels[a]}", scope, raw, m))
lab_y, seen = [], []
y = 0
for a in order:
    ax.text(-0.62, y + 1.5, labels[a], fontsize=10, weight="bold", va="center")
    for scope in ["pooled", "Rest", "Stim8hr", "Stim48hr"]:
        m = res[(res.axis == a) & (res.scope == scope)].matched_delta.values[0]
        raw = res[(res.axis == a) & (res.scope == scope)].raw_delta.values[0]
        ax.plot([raw], [y], "o", color="#bbb", ms=6)
        ax.plot([m], [y], "o", color=ACC, ms=7)
        ax.plot([raw, m], [y, y], "-", color="#ddd", lw=1, zorder=0)
        lab_y.append((y, scope)); y += 1
    y += 1
ax.axvline(0, color="#111", lw=1)
ax.set_yticks([yy for yy, _ in lab_y]); ax.set_yticklabels([s for _, s in lab_y], fontsize=9)
ax.invert_yaxis()
ax.set_xlabel(r"effect size (Cliff's $\delta$ matched / Spearman $\rho$)")
ax.set_xlim(-0.45, 0.45)
ax.plot([], [], "o", color="#bbb", label="raw (confounded)")
ax.plot([], [], "o", color=ACC, label="expression-corrected")
ax.legend(frameon=False, loc="upper right", fontsize=9)
ax.text(0.01, -0.13, "← essential LESS responsive (buffered)", transform=ax.transAxes,
        fontsize=8.5, color=ACC, ha="left")
ax.text(0.99, -0.13, "essential MORE responsive →", transform=ax.transAxes,
        fontsize=8.5, color="#999", ha="right")
ax.set_title("Buffering is robust to essentiality definition and stimulation state\n"
             "raw signal flips sign after expression correction in every stratum",
             fontsize=11, loc="left", pad=12)
fig.savefig(FIG / "fig2_robustness_forest.png"); plt.close(fig)

# ================= FIGURE 3 — shet dose-response =================
fig, ax = plt.subplots(figsize=(6.4, 4.2))
s = gene.dropna(subset=["shet"]).copy()
s["shet_dec"] = pd.qcut(s.shet, 10, labels=False, duplicates="drop")
med = s.groupby("shet_dec").resid_total.median()
ax.axhline(0, color="#999", lw=1, ls="--")
ax.plot(med.index, med.values, "-o", color=ACC, lw=2, ms=6)
ax.set_xlabel("selective constraint (shet) decile, low → high")
ax.set_ylabel("expression-corrected responsiveness\n(summed residual)")
ax.set_title("More-constrained genes are more strongly buffered\n"
             r"Spearman $\rho$=−0.09, p=3×10$^{-20}$ (n=10,166)", fontsize=10.5, loc="left")
fig.savefig(FIG / "fig3_shet_doseresponse.png"); plt.close(fig)

print("wrote:", *[p.name for p in sorted(FIG.glob("*.png"))])
