"""Comprehensive buffering analysis: 3 essentiality axes x (pooled + per-condition),
raw and baseMean-matched. Produces tidy results + per-gene table + figure data.
"""
import json, io, requests
from pathlib import Path
import numpy as np, pandas as pd, h5py, fsspec
from scipy import stats

# Gene lists fetched from public sources (portable; no local paths required).
GWT = "https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master"
BAGEL = "https://raw.githubusercontent.com/hart-lab/bagel/master"
def _fetch(url): return requests.get(url, timeout=60).text
ART = Path.home() / "CoDEG_Tcell" / "artifacts"
URL = ("https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
       "marson2025_data/GWCD4i.DE_stats.h5ad")
CONDS = ["Rest", "Stim8hr", "Stim48hr"]


def dec(a): return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def cliffs_delta(x, y):
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan
    Rx = stats.rankdata(np.concatenate([x, y]))[:nx].sum()
    U = Rx - nx * (nx + 1) / 2.0
    return 2 * U / (nx * ny) - 1


def matched_test(df, mask_e, mask_n, metric, nbins=10):
    """baseMean-decile-stratified Cliff's delta + Stouffer one-sided (essential<non)."""
    d = df.loc[mask_e | mask_n].copy()
    d["ess"] = mask_e.loc[d.index].values
    d["bm_dec"] = pd.qcut(d["baseMean"], nbins, labels=False, duplicates="drop")
    deltas, ns, zs = [], [], []
    for _, g in d.groupby("bm_dec"):
        e = g.loc[g.ess, metric].values
        n = g.loc[~g.ess, metric].values
        if len(e) < 3 or len(n) < 3:
            continue
        deltas.append(cliffs_delta(e, n)); ns.append(len(e))
        _, p = stats.mannwhitneyu(e, n, alternative="less")
        zs.append(stats.norm.isf(np.clip(p, 1e-12, 1 - 1e-12)))
    if not deltas:
        return dict(matched_delta=np.nan, stouffer_p=np.nan, n_strata=0)
    z = np.sum(zs) / np.sqrt(len(zs))
    return dict(matched_delta=float(np.average(deltas, weights=ns)),
                stouffer_p=float(stats.norm.sf(z)), stouffer_z=float(z), n_strata=len(deltas))


# ---------- load per-gene metrics ----------
h = h5py.File(fsspec.open(URL, block_size=8 * 1024 * 1024).open(), "r")
g0 = h["varm"]["measured_genes_stats_Rest"]
df = pd.DataFrame({"gene_id": dec(g0["gene_id"][:]), "gene_name": dec(g0["gene_name"][:])})
for c in CONDS:
    g = h["varm"][f"measured_genes_stats_{c}"]
    df[f"nreg_{c}"] = g["n_regulators"][:]
    df[f"resid_{c}"] = g["expected_n_regulators_residuals"][:]
df["nreg_total"] = df[[f"nreg_{c}" for c in CONDS]].sum(axis=1)
df["resid_total"] = df[[f"resid_{c}" for c in CONDS]].sum(axis=1)
df["baseMean"] = np.nan_to_num(h["layers"]["baseMean"][:1500, :]).mean(0)

# ---------- essentiality axes ----------
hart = set(l.strip() for l in _fetch(f"{GWT}/metadata/gene_lists/core_essentials_hart.tsv").splitlines() if l.strip())
cegv2 = set(pd.read_csv(io.StringIO(_fetch(f"{BAGEL}/CEGv2.txt")), sep="\t")["GENE"])
negv1 = set(pd.read_csv(io.StringIO(_fetch(f"{BAGEL}/NEGv1.txt")), sep="\t")["GENE"])
shet = pd.read_csv(io.StringIO(_fetch(f"{GWT}/src/8_lymphocyte_counts_LoF/input/shet_10bins.txt")), sep="\t")
df["ess_hart"] = df.gene_name.isin(hart)
df["ess_cegv2"] = df.gene_name.isin(cegv2)
df["noness_negv1"] = df.gene_name.isin(negv1)
df["shet"] = df.gene_id.map(dict(zip(shet.ensg, shet.shet)))

print(f"genes={len(df):,} | Hart={df.ess_hart.sum()} | CEGv2={df.ess_cegv2.sum()} "
      f"| NEGv1={df.noness_negv1.sum()} (HVG removes nonessentials!) | shet={df.shet.notna().sum():,}")

# ---------- tests ----------
rows = []
# axis 1&2: Hart-vs-rest, CEGv2-vs-rest. NEGv1 dropped: only 12/928 survive HVG selection
# (nonessential-reference genes are lowly expressed by construction) -> the HVG filter itself
# preferentially retains essential genes, which is part of why the naive comparison is confounded.
axes = [("Hart_vs_rest", df.ess_hart, ~df.ess_hart),
        ("CEGv2_vs_rest", df.ess_cegv2, ~df.ess_cegv2)]
for scope, metric in [("pooled", "nreg_total")] + [(c, f"nreg_{c}") for c in CONDS]:
    for name, me, mn in axes:
        e = df.loc[me, metric].values
        n = df.loc[mn, metric].values
        raw_d = cliffs_delta(e, n)
        _, p_raw = stats.mannwhitneyu(e, n, alternative="two-sided")
        m = matched_test(df, me, mn, metric)
        rows.append(dict(axis=name, scope=scope, n_ess=int(me.sum()), n_non=int(mn.sum()),
                         raw_delta=round(raw_d, 3), raw_p=p_raw,
                         matched_delta=round(m["matched_delta"], 3), matched_p=m["stouffer_p"],
                         n_strata=m["n_strata"]))
# axis 3: shet continuous (residual = power-corrected)
for scope, rescol, nregcol in [("pooled", "resid_total", "nreg_total")] + \
                               [(c, f"resid_{c}", f"nreg_{c}") for c in CONDS]:
    s = df.dropna(subset=["shet"])
    rho_res, p_res = stats.spearmanr(s[rescol], s.shet)
    rho_raw, p_raw = stats.spearmanr(s[nregcol], s.shet)
    rows.append(dict(axis="shet_continuous", scope=scope, n_ess=int(len(s)), n_non=0,
                     raw_delta=round(rho_raw, 3), raw_p=p_raw,
                     matched_delta=round(rho_res, 3), matched_p=p_res, n_strata=np.nan))

res = pd.DataFrame(rows)
pd.set_option("display.width", 200, "display.max_columns", 20)
print("\n", res.to_string(index=False))
res.to_csv(ART / "results_table.csv", index=False)
df.to_csv(ART / "per_gene_full.csv", index=False)

# ---------- figure data: Hart decile money-plot (pooled) ----------
d = df.copy()
d["bm_dec"] = pd.qcut(d.baseMean, 10, labels=False, duplicates="drop")
fig = d.groupby(["bm_dec", "ess_hart"])["nreg_total"].median().unstack()
fig.columns = ["non_essential", "essential"]
fig.to_csv(ART / "fig_decile_hart.csv")
json.dump({"n_genes": len(df)}, open(ART / "analyze_all_meta.json", "w"))
print("\nsaved results_table.csv, per_gene_full.csv, fig_decile_hart.csv")
