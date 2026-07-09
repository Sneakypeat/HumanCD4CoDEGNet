"""Day 1: stream adj_p_value, build the on-target-masked binary matrix, save margins.

Streams ~2.8 GB via byte-range reads. Writes to artifacts/. Idempotent: skips if B exists.
"""
import json, time, sys
from pathlib import Path
import numpy as np, h5py, fsspec, requests

URL = ("https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
       "marson2025_data/GWCD4i.DE_stats.h5ad")
ART = Path.home() / "CoDEG_Tcell" / "artifacts"
ART.mkdir(parents=True, exist_ok=True)
BLOCK, PADJ_THR = 2048, 0.10

if (ART / "B_masked.npz").exists():
    print("B_masked.npz exists; nothing to do."); sys.exit(0)

head = requests.head(URL, timeout=30)
prov = {"url": URL, "etag": head.headers.get("ETag", "").strip('"'),
        "bytes": int(head.headers.get("Content-Length", 0)),
        "accessed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
print(json.dumps(prov))

h = h5py.File(fsspec.open(URL, block_size=8 * 1024 * 1024).open(), "r")
padj = h["layers"]["adj_p_value"]
n_obs, n_var = padj.shape
assert (n_obs, n_var) == (33983, 10282)


def col(name):
    d = h["obs"][name]
    if isinstance(d, h5py.Group):
        cats = d["categories"][:]; codes = d["codes"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
        return cats[codes]
    v = d[:]
    return np.array([x.decode() for x in v]) if v.dtype.kind == "S" else v


keep = (col("ontarget_significant").astype(bool)
        & ~col("distal_offtarget_flag").astype(bool)
        & ~col("neighboring_gene_KD").astype(bool)
        & ~col("single_guide_estimate").astype(bool))
keep_idx = np.flatnonzero(keep)
print(f"kept {keep.sum():,} / {n_obs:,} perturbation-conditions", flush=True)

gene_id = np.array([g.decode() if isinstance(g, bytes) else g for g in h["var"]["gene_ids"][:]])
gene_nm = np.array([g.decode() if isinstance(g, bytes) else g for g in h["var"]["gene_name"][:]])
targets = col("target_contrast")[keep_idx]

B = np.zeros((keep.sum(), n_var), bool)
pos, t0 = 0, time.time()
for s in range(0, n_obs, BLOCK):
    e = min(s + BLOCK, n_obs)
    sel = keep[s:e]
    if not sel.any():
        continue
    p = padj[s:e, :][sel]
    p = np.where(np.isnan(p), 1.0, p)          # NaN padj -> not significant
    B[pos:pos + sel.sum()] = p < PADJ_THR
    pos += sel.sum()
    print(f"  {e/n_obs:6.1%}  {time.time()-t0:6.0f}s", flush=True)

# on-target masking
gpos = {g: i for i, g in enumerate(gene_id)}
masked = 0
for r, t in enumerate(targets):
    j = gpos.get(t)
    if j is not None:
        B[r, j] = False; masked += 1

gene_response = B.sum(0)      # per readout gene: in how many perturbations is it DE
pert_burden   = B.sum(1)      # per perturbation: how many genes it moves

np.savez_compressed(ART / "B_masked.npz", B=np.packbits(B, axis=1), shape=np.array(B.shape),
                    keep_idx=keep_idx, gene_id=gene_id, gene_name=gene_nm,
                    gene_response=gene_response, pert_burden=pert_burden,
                    n_cells_target=col("n_cells_target")[keep_idx])
summary = {**prov, "padj_thr": PADJ_THR, "n_kept": int(keep.sum()),
           "density": float(B.mean()), "nnz": int(B.sum()), "on_target_masked": int(masked),
           "gene_response_median": float(np.median(gene_response)),
           "pert_burden_median": float(np.median(pert_burden)),
           "genes_never_DE": int((gene_response == 0).sum()),
           "elapsed_s": round(time.time() - t0, 1)}
(ART / "day1_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
