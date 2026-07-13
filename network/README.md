# HumanCD4CoDEGNet — causal network resource (NDEx-ready)

The **trans-regulatory hub core** of the primary human CD4⁺ T cell, reconstructed from the
Marson-lab genome-scale CRISPRi Perturb-seq atlas (~22 M cells; Zhu, Dann et al. 2025).

- **Nodes (75):** the top out-degree regulators (validated knockdown) in each activation state,
  plus the druggable disease-gene candidates (ZAP70, ITK, LCK, PTPRC, IL12RB2).
- **Edges (2,222, directed):** `A → B` means knockdown of **A** significantly shifts **B**. Each edge
  is annotated with the activation state(s) it appears in (`Rest` / `Stim8hr` / `Stim48hr`).
- **`hub_status` node attribute** — `stable` (16), `gained` (34), `lost` (14), `peripheral` (11) —
  captures the *shape-invariant but identity-labile* rewiring: hub-dominance holds while the
  identity of the hubs turns over on activation. **858 edges are gained on activation.**

## Files
| file | use |
|---|---|
| `HumanCD4CoDEGNet.cx2` | **upload this to NDEx** (CX2) |
| `HumanCD4CoDEGNet.graphml` | open in Cytoscape |
| `nodes.csv` / `edges.csv` | universal edge/node tables |
| `HumanCD4CoDEGNet_network.png` | preview figure |
| `upload_ndex.py` | one-command NDEx upload (prompts for your login) |

## Publish on NDEx
**Web UI:** sign in at [ndexbio.org](https://www.ndexbio.org) → **Upload** → choose `HumanCD4CoDEGNet.cx2`.
**Or script:** `pip install ndex2 && python upload_ndex.py` (uploads PRIVATE; flip to PUBLIC when ready).

Rebuild from source: `scripts/build_network.py` (streams `.obs` from the public S3 file; no bulk download),
render the preview with `scripts/render_network.py`. Method lineage: **YeastCoDEGNet**
(Nasar, Rehman, Ott & Alam, *NAR* 2026, gkaf1410; NDEx `794dd2cc-bcae-11f0-a218-005056ae3c32`).
