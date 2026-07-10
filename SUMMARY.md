# Summary — A shape-invariant, identity-labile control architecture in human CD4⁺ T cells

**Claude Hackathon: Life Sciences (Researcher Track).** Data: Marson-lab genome-scale CRISPRi
Perturb-seq, ~22M primary human CD4⁺ T cells, Rest / Stim 8 h / Stim 48 h (Zhu, Dann et al. 2025,
bioRxiv 10.64898/2025.12.23.696273).

## One paragraph

We ask a purely descriptive, causal question the atlas makes answerable for the first time: **what is
the *architecture* of the trans-regulatory network in a primary human immune cell, and does that
architecture reorganize when the cell activates?** Treating every validated knockdown → trans-DE gene
as a causal edge, we find the network is **hub-dominated and sparse-but-pleiotropic** in every state —
the top 5 % of regulators drive ~78 % of all trans-effects (out-degree Gini **0.92**), while a typical
gene is controlled by only ~0.5 % of possible regulators. Regulation is **broadcast-concentrated but
reception-distributed** (out-degree Gini 0.92 vs in-degree Gini 0.35). This is the **first causal
confirmation** of the hub-dominated / sparse-but-pleiotropic topology that the Pritchard lab's own 2026
theory paper predicted from twin-study heritability alone (Barton et al., *Cell Genomics* 2026) — a
prediction the atlas paper never tested against its data. The novel twist the theory does **not**
predict: the architecture is **shape-invariant but identity-labile.** Hub-dominance is pinned at Gini
≈ 0.92 across Rest → 8 h → 48 h even as the network densifies +32 % and relaxes, while up to **59 % of
the top-100 hubs are displaced** and the TCR signalosome switches on as the dominant activation-state
broadcaster. *The cell preserves the shape of its control architecture while swapping out who is in
control.*

## The finding, in numbers

| Architecture (per state) | Rest | Stim 8 h | Stim 48 h |
|---|---|---|---|
| Out-degree Gini (hub-dominance) | **0.921** | **0.929** | **0.923** |
| Top 1 % of regulators → share of all edges | 37 % | 40 % | 35 % |
| Top 5 % of regulators → share of all edges | 77 % | 80 % | 78 % |
| In-degree Gini (regulators per gene) | 0.40 | 0.35 | 0.42 |
| Median regulators per gene (of ~11.3 k possible) | 47 (0.4 %) | 66 (0.6 %) | 50 (0.4 %) |
| Super-hubs (>1000 targets) / silent regulators (0) | 139 / 3588 | 195 / 3433 | 185 / 3249 |
| Total trans-edges | 592 k | 780 k | 663 k |

| Rewiring (across states) | value |
|---|---|
| Hub-identity turnover, top-100 (Rest→8h / Rest→48h / 8h→48h) | 44 % / **59 %** / 42 % displaced |
| Hub-rank correlation (Spearman ρ, both-testable) | 0.70 / 0.63 / 0.69 |
| Densification 8h/Rest, relax 48h/Rest | ×1.32 / ×1.12 |
| Activation super-broadcasters (validated KD in both states) | TCR signalosome: CD3E/D/G, CD247, LAT, LCP2, PLCG1, ZAP70, VAV1 (~2→5000 targets) |

## Why it is defensible (every confound guarded)

- **KD-detectability (the confound that sank our earlier attempts).** All rewiring is measured **only
  among regulators whose on-target knockdown was validated in *both* compared states** — so a hub
  "emerging" is genuine state-dependent control, not "the gene became expressible." Only 3/30 top
  gainers were detectability artifacts; the TCR signalosome survives (KD works at Rest, but changes
  ~2 genes; at 8 h the same KD changes ~5000).
- **Power.** Out-degree is *negatively* correlated with cells-per-perturbation (ρ ≈ −0.20) — hubs are
  not a sampling-depth artifact.
- **Edge definition.** Hub-dominance is identical on validated-KD-only edges (Gini 0.91 vs 0.92).
- **Regulator expression.** Out-degree ↔ regulator baseMean ρ ≈ +0.20 only — hubness is not KD
  efficiency. Silent regulators (~30 %) are half validated KDs *with more cells* — genuinely inert, not
  failed.

## What is novel (and what is not)

- **Novel:** (1) first **causal** test of the Barton/Pritchard (2026) hub-dominated / sparse-but-
  pleiotropic prediction — prior evidence was simulation on twin-study heritability, no perturbation
  data; (2) the **shape-invariant / identity-labile** principle — a quantitative architectural constant
  (Gini ≈ 0.92) that holds through massive edge turnover; (3) the **broadcast-vs-reception asymmetry**
  (out-Gini 0.92 vs in-Gini 0.35) in a causal human GRN.
- **Not novel (we do not claim it):** that regulators are context-specific (Zhu/Marson 2025, same data
  — but they report the raw `n_regulators`/`n_downstream` fields without the topology framing, the
  invariance constant, or the guarded turnover); that a "conserved core / rewired periphery" exists in
  activated CD4⁺ T cells (Mihai 2025 — but that is *correlational* multiome, not causal perturbation).

## What this is *not*

A drug or a mechanism. It is a **regulatory-architecture law**: the topology of causal control in a
human immune cell is quantitatively conserved in shape and almost entirely reassigned in identity
across activation. No fragile causal claim — the result is a distribution-level invariant that survives
every confound we could pose.

## Generality (external replication)

The hub-dominated, broadcast-concentrated architecture is **not T-cell-specific**: the identical pipeline
on **Replogle 2022** genome-scale Perturb-seq recovers it in **K562** (out-Gini 0.86, asymmetry +0.40 —
a non-immune line, different lab/platform), and across K562↔RPE1 the shape is roughly conserved while hub
identity turns over ~80 % (a cross-cell-type echo of the within-T-cell result). The *tight dynamic
invariance* remains T-cell-specific — no comparator has a matched activation timecourse.
(`fig_arch3_generality.png`, `generality_replogle.py`.)

## Reproducibility

Runs in **seconds from two released summary fields** — the per-perturbation `.obs`
(`n_downstream`, `ontarget_significant`) and the per-gene `n_regulators` — no 16.8 GB download.
Effect-size sensitivity over the FDR/edge definition included. Code + figures public:
https://github.com/Sneakypeat/HumanCD4CoDEGNet
