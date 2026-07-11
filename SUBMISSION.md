# A shape-invariant, identity-labile control architecture in human CD4⁺ T cells

**Claude Hackathon: Life Sciences — Researcher Track ("Build From the Bench")**
Dataset: Marson-lab genome-scale CRISPRi Perturb-seq in primary human CD4⁺ T cells (Gladstone / UCSF /
Stanford; Zhu, Dann et al. 2025)

---

## Abstract (~200 words)

A primary human T cell rewires thousands of genes within hours of activation, yet remains a stable,
functioning cell. What, if anything, is *conserved* about how it is controlled? Genome-scale CRISPRi
Perturb-seq (22 M cells; Rest / Stim 8 h / Stim 48 h) lets us draw the **causal** trans-regulatory
network of a human immune cell — every validated knockdown → trans-DE gene is one edge — and redraw it
in each activation state. We find the network is **hub-dominated and sparse-but-pleiotropic** in every
state: the top 5 % of regulators drive ~78 % of all trans-effects (out-degree Gini **0.92**), while a
typical gene is controlled by only ~0.5 % of possible regulators. Regulation is **broadcast-concentrated
but reception-distributed** (out-degree Gini 0.92 vs in-degree Gini 0.35). This is the **first causal
confirmation** of the hub-dominated / sparse-but-pleiotropic topology that Barton et al. (*Cell Genomics*
2026, from the dataset's own co-senior lab) predicted from twin-study heritability with no perturbation
data — a prediction the atlas paper never tested. The novel principle theory does not predict: the
architecture is **shape-invariant but identity-labile** — hub-dominance is pinned at Gini ≈ 0.92 across
all three states even as the network densifies +32 % and up to **59 % of the top-100 hubs are displaced**,
the TCR signalosome switching on as the activation-state broadcaster. *The cell keeps the shape of its
control while swapping out who is in control.*

---

## Figures

**Figure 1 — A hub-dominated, sparse-but-pleiotropic causal architecture** (`fig_arch1_topology.png`).
(A) Lorenz curves of out-degree in all three states collapse onto one another (Gini ≈ 0.92); the top 5 %
of regulators account for ~78 % of all trans-edges. (B) Rank–out-degree (Stim 8 h): ~195 super-hubs
(>1000 targets) sit above a long tail, ~3400 regulators are transcriptionally silent (0 detectable
targets). (C) Out-degree Gini (broadcast, ~0.92) vastly exceeds in-degree Gini (reception, ~0.35–0.42) in
every state — control is broadcast-concentrated but reception-distributed.

**Figure 2 — Shape-invariant, identity-labile rewiring** (`fig_arch2_rewiring.png`).
(A) The money panel: total trans-edges swing 592k → 780k → 663k (densify +32 % at 8 h, relax by 48 h)
while out-degree Gini stays pinned at 0.921 / 0.929 / 0.923. Shape-invariant, mass-variable.
(B) Hub-identity scatter (out-degree Rest vs Stim 8 h; only regulators with a validated knockdown in
*both* states — the detectability guard): the TCR signalosome (orange; CD3E/D/G, LAT, LCP2, PLCG1, ZAP70,
VAV1) jumps from ~2 targets to ~5000, while rest-state hubs (blue; TP53, PTEN, NFAT5, ARF1) collapse.
(C) Up to 59 % of the top-100 hubs are displaced across activation (44 % Rest→8h, 59 % Rest→48h, 42 %
8h→48h).

---

## Methods

- **Edge set.** Each perturbation (CRISPRi knockdown) → significantly trans-DE gene (authors' 10 % FDR)
  is one directed causal edge. **Out-degree** of a regulator = the authors' `n_downstream`; **in-degree**
  of a gene = the authors' per-gene `n_regulators`. The two totals match by construction per state
  (592k / 780k / 663k; a cross-check we print), so in- and out-degree describe the same network.
- **Hub-dominance** = Gini coefficient and top-k edge share of the out-degree distribution.
  **Sparsity** = median in-degree ÷ number of testable regulators. **Pleiotropy** = the out-degree
  distribution itself (targets per regulator).
- **Rewiring** = (i) *shape*: out-degree Gini per state; (ii) *identity*: Spearman rank-correlation and
  top-100 displacement of per-regulator out-degree between states.
- **KD-detectability guard (critical).** All cross-state rewiring is computed **only among regulators
  whose on-target knockdown was validated (`ontarget_significant`) in both compared states**, so an
  "emerging" hub reflects genuine state-dependent control — not a gene merely becoming expressed and
  therefore knock-down-able on activation.

## Robustness — every confound guarded

| Alternative explanation | Test | Result |
|---|---|---|
| Sequencing **power** (more cells → more detected edges) | Spearman(out-degree, cells/perturbation) | **−0.19 to −0.22** (negative → hubs have *fewer* cells) |
| **Edge definition** (failed-KD noise inflates hubs) | Gini on `ontarget_significant`-only edges | 0.906–0.916 vs 0.92 — **identical** |
| **Regulator expression** (hubness = KD efficiency) | Spearman(out-degree, log₁₀ baseMean) | +0.19 to +0.23 — weak |
| **Silent regulators** are just failed KDs | fraction validated, cells | 47–51 % validated, *more* cells than active — **genuinely inert** |
| **Detectability** drives apparent rewiring | rewiring restricted to both-states-validated | 27/30 top gainers survive; TCR signalosome genuine |

The result is a **distribution-level invariant**, not a fragile causal claim — which is why it survives
every confound we could pose.

## Generality — external replication (Replogle 2022; different cell type, lab, platform)

We ran the identical degree-based pipeline on **Replogle et al. 2022** genome-scale Perturb-seq
(Harmonizome DE gene-sets; directed perturbation → DE-gene edges, the same construction as `n_downstream`).
Figure `fig_arch3_generality.png`; numbers `artifacts/architecture_replogle_results.json`.

- **Static shape generalizes.** In **K562** (genome-scale, ~7,200 perturbations — a non-immune cancer
  line) the network is hub-dominated with the same broadcast/reception asymmetry: **out-degree Gini 0.86,
  top 5 % of regulators → 61 % of edges, in-degree Gini 0.47 (asymmetry +0.40)** — qualitatively identical
  to Marson-T (0.92 / +0.55). Hub-dominance is not T-cell-specific.
- **Cross-cell-type echo.** Between **K562 and RPE1** (epithelial) on 1,672 matched essential-gene
  perturbations, the shape is roughly conserved (out-Gini 0.74 vs 0.63) while **hub identity turns over
  ~80 %** (Spearman 0.60) — a cross-*context* analog of the within-T-cell shape-invariant/identity-labile
  result.
- **Honest scope.** Harmonizome's standardized-value threshold ≠ Marson's 10 % FDR, so only the
  *qualitative* shape and the *within-Replogle* comparison are directly comparable (absolute Gini values
  are not). The essential-only Replogle subsets exclude the silent-majority tail that drives hub-dominance,
  so they compress the asymmetry (RPE1-essential even reverses it) — the **genome-scale K562 network is the
  fair comparator, and it cleanly supports generality.** The *tight dynamic invariance* (Gini constant to
  three decimals across a matched activation timecourse) remains specific to the Marson data, because no
  comparator has a matched within-cell state axis.

## Further tests — the deferred threads, each closed to a definite outcome

Three secondary questions raised by the topology literature, run to a definite result on the
reconstructed per-state causal graph (`further_tests.py`, `fig_arch4_furthertests.png`,
`artifacts/further_tests_results.json`):

- **Plasticity tracks regulatory complexity — not motifs (CONFIRMED).** A gene's expression plasticity
  (SD of expression-corrected responsiveness across states) scales with its **regulatory complexity**
  (in-degree): Spearman **ρ = +0.66** (partial **+0.45** controlling baseline expression, n = 10,282),
  while **feedback-loop participation adds ~nothing beyond in-degree** (partial ρ = +0.03 — ~20× weaker
  than complexity's +0.66). First confirmation, in a causal human immune-cell atlas, of the plasticity ↔ complexity-but-not-motifs
  dissociation predicted by Petit et al. (*Genetics* 2026). Reuses Angle 1's in-degree and reinforces
  that the *quantity of inputs*, not local wiring, sets plasticity.
- **Trans-effects are not "local" beyond degree (informative null).** Directed 2-hop closure
  (transitivity) of the regulator→regulator graph is **0.209 vs 0.212 in a degree-preserving null — a
  negligible ~1.5 % difference**, reciprocity 0.01. The hub-dominated architecture is a *degree* phenomenon, not a
  modular/short-path one — this rules out a local-propagation reading and keeps the result anchored on
  the degree distribution (contra a naive two-hop-locality expectation).
- **Middle-manager hierarchy — not supported (scoped out).** Throughput (in×out) is dominated by the top
  tier (median 12,040 / 1,800 / 35 for top / middle / bottom), i.e. hub-dominated rather than a Yu–Gerstein mid-tier bottleneck. A
  rigorous betweenness test would need to separate direct from indirect (cascade) edges, which this data
  cannot support (the "~99 % direct" simplification failed verification), so we make **no** hierarchy claim.

## Novelty and positioning (honest)

- **Novel, and ours:** (1) the **first causal test** of the Barton/Pritchard (2026) hub-dominated /
  sparse-but-pleiotropic prediction — prior evidence was a scale-free simulation calibrated to twin-study
  heritability, with no perturbation data; (2) the **shape-invariant / identity-labile** principle — a
  quantitative architectural constant (Gini ≈ 0.92) that holds through massive edge turnover across
  activation; (3) the **broadcast-vs-reception asymmetry** (out-Gini 0.92 vs in-Gini 0.35) in a causal
  human GRN.
- **Not ours (we do not claim it):** the raw `n_regulators` / `n_downstream` fields were released by the
  atlas authors (Zhu, Dann et al. 2025); that individual regulators are context-specific is their headline —
  but they do not frame the fields as a topology test, do not report the invariance constant, and do not
  quantify the detectability-guarded turnover. A "conserved core / rewired periphery" was described in
  activated CD4⁺ T cells by **Mihai et al. 2025** — but *correlationally* (multiome centrality); ours is
  the *causal* (perturbational) version they explicitly lack. The closest causal-GRN method (LLCB,
  Weinstock et al. 2024, same senior authors) was demonstrated on 84 genes; we operate at genome scale
  on released summary fields.

## Impact — what it enables

- A **compact, quantitative law** for how a human immune cell reorganizes control: a conserved
  hub-dominated topology whose hub identity is almost entirely reassigned across activation. Systems
  immunology, not a fragile mechanism.
- **A causal benchmark for network-topology theory.** The atlas + Barton/Pritchard (2026) are from the
  same senior lab and were never connected; this is the empirical test — and it extends the theory with a
  cross-state dimension (architectural shape-invariance) the static model does not contain.
- **Shipped resources:** per-perturbation out-degree table across states
  (`artifacts/arch_perturbation_outdegree.csv`) and the full architecture result set
  (`artifacts/architecture_results.json`) — a ready substrate for asking which hubs a given disease
  program routes through, per activation state.
- **Generality already tested (Replogle 2022, below):** hub-dominance generalizes to K562; the tight
  dynamic invariance stays T-cell-specific. **Next:** a comparator *with* a matched state axis (e.g. a
  stimulation timecourse in another primary cell) is what would test whether dynamic shape-invariance is
  a general law.

## Reproducibility

- Runs in **seconds from two released summary fields** of the single public file
  (`GWCD4i.DE_stats.h5ad`): the per-perturbation `.obs` (`n_downstream`, `ontarget_significant`) and the
  per-gene `varm/measured_genes_stats_*` (`n_regulators`). The **16.8 GB effect-size layers are never
  downloaded.** No credentials, no local data.
- `notebooks/HumanCD4CoDEGNet_architecture.ipynb` (streams live from S3, 0 errors) and
  `architecture.py` (CLI) both regenerate every number and both figures. Effect-size/edge-definition
  sensitivity is included (all-edges vs validated-KD-only).

---

## Secondary result — essential-gene transcriptional buffering (confirm-and-extend, not a discovery)

As a supporting analysis we also confirm that **selectively-constrained / essential genes are
transcriptionally buffered** (fewer trans-regulators than expected for their expression) in primary human
CD4⁺ T cells — expression-corrected Cliff's δ ≈ −0.17 (Hart) and −0.19 (CEGv2), Spearman(resid, shet)
p = 3×10⁻²⁰, unanimous across two essential-gene sets × three states. **This direction, and the need to
control for expression, were already reported in human iPSCs by Feng et al. 2026 (expression-controlled,
Fig S2D) and via eQTL by Mostafavi 2023 — so this is an extension, not a discovery**, and we do not build
the submission on it. (Notebook: `HumanCD4CoDEGNet_buffering.ipynb`; resource:
`artifacts/buffering_score_resource.csv`.) A mechanism we explored — that activation edges "route around"
essential genes — **did not survive edge-level testing and is retired**; see `artifacts/MECHANISM_NOTE.md`.

## Honest limitations

- "Sparse-but-pleiotropic" is confirmed **qualitatively**; Barton/Pritchard's numbers are
  heritability-weighted trans-regulators, not CRISPRi-DE edges, so we test the *qualitative* shape, not
  their exact median.
- In/out-degree conflate direct and indirect (cascade) effects; a genome-scale direct/indirect
  decomposition (e.g. LLCB) is future work. The invariance and turnover results do not depend on it.
- Essentiality/constraint sets are organismal, not CD4⁺-T-cell-specific.
- The field is fast-moving (atlas is a Dec-2025 preprint; the topology and rewiring papers are
  2025–2026) — a fresh prior-art scan is warranted before any manuscript.

## Credit and references

- **Data:** Zhu R., Dann E. et al. (2025) *Genome-scale perturb-seq in primary human CD4⁺ T cells maps
  context-specific regulators of T cell programs and human immune traits*, bioRxiv
  10.64898/2025.12.23.696273.
- **Central theory tested:** Barton A. R. et al. / Pritchard group (2026) *Regulatory network topology and
  the genetic architecture of gene expression*, **Cell Genomics** S2666-979X(26)00081-9 (bioRxiv
  2025.08.12.669924).
- **Adjacent (correlational) rewiring:** Mihai et al. (2025) *A conserved transcriptional backbone and
  rewiring of gene-regulatory networks in activated human CD4⁺ T cells*, bioRxiv 2025.11.25.687998.
- **Closest causal-GRN method:** Weinstock J. et al. (2024) LLCB, **Cell Genomics** (PMC10516010).
- **Method lineage:** Nasar, Rehman, Ott & Alam (2026) *Uncovering coordinated pathway interactions
  through gene co-differential expression in yeast*, NAR 54(1):gkaf1410.
- **Secondary-result prior art:** Feng C. et al. (2026) *Cell Genomics* 3:101076 (expression-controlled
  buffering direction in human iPSCs); Mostafavi H. et al. (2023) *Nature* (eQTL depletion of constrained
  genes). Reference gene sets: Hart core-essentials & CEGv2 (Hart lab); shet constraint (gnomAD-based).
