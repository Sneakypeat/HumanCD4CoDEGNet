"""Is the 'essential genes buffered' signal real, or a GLM tail-fitting artifact?

The residual test assumes a Poisson fit of n_regulators ~ log10(baseMean). If essential
genes sit at the extreme high-expression tail where that fit misbehaves, the negative
residual could be artifactual. We replace the parametric residual with a NON-PARAMETRIC
baseMean-matched comparison: within strata of real per-gene baseMean, compare essential
vs non-essential raw responsiveness. If essential stay lower within matched strata, the
buffering is real.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd, h5py, fsspec
from scipy import stats

ART = Path.home() / "CoDEG_Tcell" / "artifacts"
URL = ("https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
       "marson2025_data/GWCD4i.DE_stats.h5ad")


def dec(a): return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


h = h5py.File(fsspec.open(URL, block_size=8 * 1024 * 1024).open(), "r")
base = h["layers"]["baseMean"]
n_obs, n_var = base.shape

# real per-gene baseMean = mean over a CONTIGUOUS block of perturbation rows.
# baseMean is a per-comparison baseline dominated by the gene's own expression level,
# so it is stable across perturbations -> a contiguous 1500-row slice is unbiased and
# reads only ~123 MB (vs re-reading the whole layer for a scattered sample).
blk = base[:1500, :]
gene_baseMean = np.nan_to_num(blk).mean(0)

df = pd.read_csv(ART / "responsiveness_essentiality.csv")
df["baseMean"] = gene_baseMean
print(f"baseMean essential median={df.loc[df.essential_hart,'baseMean'].median():.1f}  "
      f"non={df.loc[~df.essential_hart,'baseMean'].median():.1f}  "
      f"(confirms essential are higher-expressed)")

# ---- baseMean-matched (decile-stratified) Mann-Whitney on RAW responsiveness ----
df["bm_dec"] = pd.qcut(df["baseMean"], 10, labels=False, duplicates="drop")
print("\nWithin baseMean deciles: essential vs non raw n_regulators")
print(f"{'dec':>3} {'n_ess':>5} {'med_ess':>8} {'med_non':>8} {'cliffs_d':>9} {'p_less':>9}")
deltas, ns = [], []
for d, g in df.groupby("bm_dec"):
    e = g.loc[g.essential_hart, "n_regulators"].values
    nn = g.loc[~g.essential_hart, "n_regulators"].values
    if len(e) < 3:
        print(f"{d:>3} {len(e):>5}   (too few essential)"); continue
    U, p = stats.mannwhitneyu(e, nn, alternative="less")
    delta = 2 * (stats.rankdata(np.r_[e, nn])[:len(e)].sum() - len(e) * (len(e) + 1) / 2) / (len(e) * len(nn)) - 1
    deltas.append(delta); ns.append(len(e))
    print(f"{d:>3} {len(e):>5} {np.median(e):>8.1f} {np.median(nn):>8.1f} {delta:>+9.3f} {p:>9.3f}")

wdelta = float(np.average(deltas, weights=ns))
print(f"\nweighted mean Cliff's delta across baseMean strata = {wdelta:+.3f}")
print("  negative => essential LESS responsive at matched expression => buffering is REAL")
print("  positive => raw 'more responsive' was pure power; buffering is a GLM artifact")

# stratified combined test: Stouffer over per-stratum one-sided p (essential < non).
# z_i = Phi^-1(1 - p_i): small p -> large +z (buffering); p>0.5 -> -z (against). No manual flip.
zs = []
for d, g in df.groupby("bm_dec"):
    e = g.loc[g.essential_hart, "n_regulators"].values
    nn = g.loc[~g.essential_hart, "n_regulators"].values
    if len(e) < 3:
        continue
    U, p = stats.mannwhitneyu(e, nn, alternative="less")
    zs.append(stats.norm.isf(np.clip(p, 1e-12, 1 - 1e-12)))
stouffer_z = np.sum(zs) / np.sqrt(len(zs))
stouffer_p = stats.norm.sf(stouffer_z)
print(f"Stouffer combined one-sided p (essential<non within strata) = {stouffer_p:.3e}  "
      f"(z={stouffer_z:.2f}, {len(zs)} strata)")

json.dump({"gene_baseMean_ess_median": float(df.loc[df.essential_hart, "baseMean"].median()),
           "gene_baseMean_non_median": float(df.loc[~df.essential_hart, "baseMean"].median()),
           "weighted_cliffs_delta_matched": wdelta,
           "stouffer_p_matched": float(stouffer_p),
           "per_stratum_delta": [float(x) for x in deltas]},
          open(ART / "robustness_matched.json", "w"), indent=2)
print(f"\nsaved -> {ART/'robustness_matched.json'}")
