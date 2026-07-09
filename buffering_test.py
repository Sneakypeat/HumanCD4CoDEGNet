"""Headline test: are essential / constrained genes transcriptionally buffered
(less responsive) in primary human CD4+ T cells, as in yeast?

Responsiveness = how many perturbations significantly move a gene (its 'n_regulators').
We use the authors' own per-gene metric AND their baseMean power-correction
(expected_n_regulators_residuals, Supp Fig 6), then relate it to essentiality
(Hart core essentials) and constraint (shet) -- the join they did NOT make.

Hypothesis (from YeastCoDEGNet): essential/constrained genes are LESS responsive.
"""
import json, io, requests
from pathlib import Path
import numpy as np, pandas as pd, h5py, fsspec
from scipy import stats

# Gene lists are fetched from their public sources (portable; no local paths required).
GWT = "https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master"
def _fetch(url): return requests.get(url, timeout=60).text
ART = Path.home() / "CoDEG_Tcell" / "artifacts"
URL = ("https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
       "marson2025_data/GWCD4i.DE_stats.h5ad")


def dec(a): return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def cliffs_delta(x, y):
    """P(x>y) - P(x<y) via rank sum. x = essential group."""
    nx, ny = len(x), len(y)
    ranks = stats.rankdata(np.concatenate([x, y]))
    Rx = ranks[:nx].sum()
    U = Rx - nx * (nx + 1) / 2.0          # Mann-Whitney U for x
    return 2 * U / (nx * ny) - 1


# ---- 1. authors' per-gene responsiveness (n_regulators) + power-corrected residual ----
h = h5py.File(fsspec.open(URL, block_size=4 * 1024 * 1024).open(), "r")
conds = ["Rest", "Stim8hr", "Stim48hr"]
g0 = h["varm"]["measured_genes_stats_Rest"]
gene_id = dec(g0["gene_id"][:])
gene_name = dec(g0["gene_name"][:])
nreg = np.zeros(len(gene_id))
resid = np.zeros(len(gene_id))
for c in conds:
    g = h["varm"][f"measured_genes_stats_{c}"]
    nreg += g["n_regulators"][:]
    resid += g["expected_n_regulators_residuals"][:]        # sum of power-corrected residuals
df = pd.DataFrame({"gene_id": gene_id, "gene_name": gene_name,
                   "n_regulators": nreg, "resid_responsiveness": resid})

# my independent responsiveness (raw), as a reimplementation cross-check
Z = np.load(ART / "B_masked.npz", allow_pickle=True)
df["my_gene_response"] = Z["gene_response"].astype(float)

# ---- 2. essentiality: Hart core essentials (their own list) ----
hart = set(l.strip() for l in _fetch(f"{GWT}/metadata/gene_lists/core_essentials_hart.tsv").splitlines() if l.strip())
df["essential_hart"] = df["gene_name"].isin(hart)

# ---- 3. constraint: shet (continuous) ----
shet = pd.read_csv(io.StringIO(_fetch(f"{GWT}/src/8_lymphocyte_counts_LoF/input/shet_10bins.txt")), sep="\t")
shet_map = dict(zip(shet["ensg"], shet["shet"]))
df["shet"] = df["gene_id"].map(shet_map)

print(f"genes: {len(df):,} | Hart essential in set: {df.essential_hart.sum()} "
      f"| shet mapped: {df.shet.notna().sum():,}")

# ---- 4. TESTS ----
out = {"n_genes": int(len(df)), "n_essential": int(df.essential_hart.sum())}

print("\n" + "=" * 70)
print("HYPOTHESIS: essential/constrained genes are LESS responsive (buffered)")
print("=" * 70)

for metric, label in [("resid_responsiveness", "power-corrected residual (authors' Supp Fig 6)"),
                      ("n_regulators", "raw n_regulators (authors)"),
                      ("my_gene_response", "raw responsiveness (my independent build)")]:
    e = df.loc[df.essential_hart, metric].values
    n = df.loc[~df.essential_hart, metric].values
    U, p_two = stats.mannwhitneyu(e, n, alternative="two-sided")
    _, p_less = stats.mannwhitneyu(e, n, alternative="less")   # essential < non
    d = cliffs_delta(e, n)
    print(f"\n[{label}]")
    print(f"  median essential={np.median(e):.2f}  non-essential={np.median(n):.2f}")
    print(f"  Mann-Whitney two-sided p={p_two:.3e}   one-sided(ess<non) p={p_less:.3e}")
    print(f"  Cliff's delta (ess vs non) = {d:+.3f}   "
          f"[{'LESS responsive (supports)' if d < 0 else 'MORE responsive (refutes)'}]")
    out[metric] = {"median_essential": float(np.median(e)), "median_non": float(np.median(n)),
                   "p_two_sided": float(p_two), "p_one_sided_less": float(p_less),
                   "cliffs_delta": float(d)}

# continuous constraint
print("\n" + "-" * 70)
sub = df.dropna(subset=["shet"])
for metric in ["resid_responsiveness", "n_regulators"]:
    rho, p = stats.spearmanr(sub[metric], sub["shet"])
    print(f"Spearman({metric}, shet) = {rho:+.3f}  p={p:.3e}  "
          f"[neg = constrained genes less responsive = supports]")
    out[f"spearman_{metric}_shet"] = {"rho": float(rho), "p": float(p), "n": int(len(sub))}

# shet decile trend (robust, non-parametric)
sub = sub.copy()
sub["shet_dec"] = pd.qcut(sub["shet"], 10, labels=False, duplicates="drop")
trend = sub.groupby("shet_dec")["resid_responsiveness"].median()
out["shet_decile_median_resid"] = {int(k): float(v) for k, v in trend.items()}
print("\nresidual responsiveness by shet decile (0=least,9=most constrained):")
print("  " + "  ".join(f"{v:.1f}" for v in trend.values))

(ART / "buffering_test_results.json").write_text(json.dumps(out, indent=2))
df.to_csv(ART / "responsiveness_essentiality.csv", index=False)
print(f"\nsaved -> {ART/'buffering_test_results.json'}")
