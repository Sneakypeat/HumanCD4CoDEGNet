# HumanCD4CoDEGNet — a shape-invariant, identity-labile control architecture

**Claude Hackathon: Life Sciences** (Researcher Track), 7–13 July 2026.
📄 **Live poster:** https://sneakypeat.github.io/HumanCD4CoDEGNet/

**Finding.** The *causal* trans-regulatory network of primary human CD4⁺ T cells is **hub-dominated and
sparse-but-pleiotropic** in every activation state — the first causal confirmation of the
Barton/Pritchard (*Cell Genomics* 2026) topology prediction, which was made from twin-study heritability
with no perturbation data. The architecture is **shape-invariant but identity-labile**: hub-dominance is
pinned at **Gini ≈ 0.92** across Rest → 8 h → 48 h even as the network densifies +32 % and up to **59 %
of the top-100 hubs are displaced**, with the TCR signalosome switching on as the activation-state
broadcaster. *The cell keeps the shape of its control while swapping out who is in control.*

Data: Marson-lab genome-scale CRISPRi Perturb-seq, ~22 M primary human CD4⁺ T cells (Zhu, Dann et al.
2025, bioRxiv 10.64898/2025.12.23.696273). Method lineage: **YeastCoDEGNet** (Nasar, **Rehman**, Ott &
Alam, *NAR* 2026, gkaf1410).

## Result in one line

| view (per activation state) | Rest | Stim 8h | Stim 48h |
|---|---|---|---|
| out-degree Gini (hub-dominance) | **0.921** | **0.929** | **0.923** |
| top 5 % of regulators → share of all trans-edges | 77 % | 80 % | 78 % |
| median regulators per gene (of ~11.3k possible) | 47 | 66 | 50 |
| total trans-edges | 592k | 780k | 663k |
| top-100 hubs displaced vs Rest | — | 44 % | 59 % |

Shape (Gini) invariant; identity (which genes are hubs) turns over. Every confound guarded
(power ρ ≈ −0.20; validated-KD-only Gini 0.91; not KD-efficiency).

**Disease lead:** the hubs that switch on *specifically* with activation are ~2× enriched for
monogenic-disease (ClinVar) genes vs stable hubs (OR 2.3, p = 0.01) — a shortlist of druggable
state-specific control points (ZAP70, ITK, LCK, PTPRC, IL12RB2). See `disease_hubs.py`.

## Run it

```bash
cd ~/CoDEG_Tcell
python3.12 -m venv venv && source venv/bin/activate
pip install -r requirements.txt

# either:
jupyter lab notebooks/HumanCD4CoDEGNet_architecture.ipynb   # streams live, 0 errors, figures inline
# or:
python scripts/architecture.py                              # CLI: writes results JSON + both figures
```

Reproduces in **seconds from two released summary fields** of the single public S3 file
(`GWCD4i.DE_stats.h5ad`): per-perturbation `.obs` (`n_downstream`, `ontarget_significant`) and per-gene
`varm/measured_genes_stats_*` (`n_regulators`) — the atlas's authoritative DE statistics, streamed by
byte-range. No credentials, no local data.

## Layout

```
notebooks/
  HumanCD4CoDEGNet_architecture.ipynb  THE submission notebook (streams live, executed, figures inline)
scripts/
  architecture.py           degree distributions, hub-dominance, rewiring, confound guards, figures
  generality_replogle.py    external replication: Replogle 2022 K562/RPE1 (hub-dominance generalizes)
  further_tests.py          plasticity~complexity (not motifs), locality null, hierarchy scoped out
  disease_hubs.py           disease lead: state-specific hubs enriched for monogenic-disease / druggable genes
  analyze_all.py, run_day1.py         build the committed data inputs (per_gene_full.csv, B_masked.npz)
  build_poster.py, build_presenter.py regenerate the poster (docs/index.html) and the presenter
artifacts/
  figures/fig_arch1..5.png             topology · rewiring · generality · further-tests · disease
  architecture_results.json            architecture + rewiring + confound-guard numbers
  architecture_replogle_results.json   generality-test numbers (K562/RPE1)
  further_tests_results.json           plasticity / locality / hierarchy outcomes
  disease_hubs_results.json            disease-enrichment numbers + candidate shortlist
  arch_perturbation_outdegree.csv      per-perturbation out-degree across states (shipped resource)
  per_gene_full.csv, B_masked.npz      committed data inputs (per-gene stats; on-target-masked edge matrix)
  MECHANISM_NOTE.md                    a mechanism we tested and RETIRED (refuted) — kept for honesty
docs/index.html           the poster, served as a GitHub Pages site
video/presenter.html      self-advancing teleprompter for the demo video
SUMMARY.md                one-page pitch          SUBMISSION.md   full writeup
```

## What we verified (not assumed)

- Dataset shapes / obs fields / varm groups read from the file, not the docs.
- **Edge-total cross-check:** sum of out-degree == sum of in-degree per state (592k / 780k / 663k) —
  in- and out-degree describe the same causal network.
- **KD-detectability guard:** all cross-state rewiring restricted to regulators whose knockdown was
  validated (`ontarget_significant`) in *both* compared states — so hub emergence is genuine, not "the
  gene became expressible."
- Hub-dominance is identical on validated-KD-only edges (Gini 0.91 vs 0.92); out-degree is *negatively*
  correlated with cells-per-perturbation (not a power artifact).

## Honesty

- The raw `n_regulators` / `n_downstream` fields were released by the atlas authors; our contribution is the topology
  framing, the causal test of Barton/Pritchard (2026), the invariance constant, and the guarded turnover.
  See `SUBMISSION.md` for full positioning vs Zhu/Marson 2025 and Mihai 2025 (correlational).
- A secondary essential-gene *buffering* result is real but **not a discovery** (Feng et al. 2026 reported
  the direction, expression-controlled, in iPSCs). A mechanism we explored ("activation routes around
  essential genes") was **refuted and retired** — see `artifacts/MECHANISM_NOTE.md`.

## Relevant work

Nasar MI, Rehman SSU, Ott S, Alam MT. *Uncovering coordinated pathway interactions through gene
co-differential expression in yeast.* Nucleic Acids Research 54(1), 2026, gkaf1410.
Yeast network: NDEx `794dd2cc-bcae-11f0-a218-005056ae3c32`.
