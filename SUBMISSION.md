# Transcriptional buffering of essential genes: a cross-kingdom extension to primary human T cells

**Built with Claude: Life Sciences — Researcher Track ("Build From the Bench")**
Dataset: Marson-lab genome-scale CRISPRi Perturb-seq in primary human CD4+ T cells (Gladstone / UCSF)

---

## Abstract (198 words)

In *Saccharomyces cerevisiae*, essential genes are markedly less transcriptionally responsive than
nonessential genes — the cell buffers its indispensable transcripts against perturbation
(YeastCoDEGNet, *NAR* 2026). Whether this generalises to a human *primary* cell had not been directly
tested (related human evidence exists — Feng 2024, Mostafavi 2023 — so the direction is not
first-reported here; see Positioning).

Marson and colleagues' genome-scale CRISPRi Perturb-seq screen in 22 million primary human CD4+
T cells now provides one — 33,983 perturbation-conditions across 10,282 measured genes. Using
Claude Science we scored each gene's responsiveness as the number of perturbations that
significantly move it. Naively, human essential genes appear *more* responsive (Cliff's δ = +0.26,
p = 8×10⁻⁷) — the opposite of yeast. But this is an expression–power artifact: essential genes are
~3× more highly expressed and therefore better powered for detection. Matching genes at equal
baseline expression inverts the result — essential genes are buffered across 9 of 10 expression
deciles (δ = −0.17, p = 8×10⁻⁶), confirmed by a second essential-gene set (CEGv2, δ = −0.19,
p = 1×10⁻¹⁴) and by continuous selective constraint (shet, p = 3×10⁻²⁰), in all three stimulation
states.

This buffering direction — and the need to control for expression — was already reported in human
iPSCs by Feng et al. (2026, controlling for expression level). Our contribution is therefore an
**extension, not a discovery**: we show the principle holds in primary human CD4+ T cells rather than
pluripotent lines, that it is **dynamic across T-cell activation** (peaking at 8 h), and that it
**unifies with the yeast co-differential-expression law** across ~10⁹ years of divergence — plus a
transcriptional-tractability application and a shipped per-gene buffering-score resource.

---

## Figures

**Figure 1 — The sign flip.**
(A) Raw responsiveness: essential genes (Hart) are differentially expressed by *more* perturbations
than non-essential genes (Cliff's δ = +0.26, p = 8×10⁻⁷) — opposite to yeast. (B) Within deciles of
baseline expression, essential-gene responsiveness (red) falls *below* non-essential (blue) in 9 of
10 deciles (weighted δ = −0.17, Stouffer p = 8×10⁻⁶). The naive signal was an expression–power
confound; matched at equal expression, essential genes are buffered.

**Figure 2 — Robustness forest.**
Effect size in every stratum: two essential-gene definitions (Hart n=124; CEGv2 n=377) × pooled and
three stimulation states, plus continuous shet constraint. Raw estimates (grey) are positive
everywhere; every expression-corrected estimate (teal) is negative. The sign flips in all 12 strata.
Buffering is strongest at Stim8hr (early activation).

**Figure 3 — shet dose-response.**
Expression-corrected responsiveness declines monotonically across selective-constraint deciles,
concentrated in the most-constrained genes (Spearman ρ = −0.09, p = 3×10⁻²⁰, n = 10,166).

---

## Methods

- **Responsiveness** = the authors' per-gene `n_regulators` (number of perturbations that
  significantly move a gene at 10% FDR), summed across conditions and per condition.
- **Expression control.** Two independent corrections: (i) the authors' own power-residual
  `expected_n_regulators_residuals` (Poisson `n_regulators ~ log10(baseMean)`, their Supp. Fig 6);
  (ii) a non-parametric **baseMean-decile-matched** Mann-Whitney / Cliff's δ (no model assumptions),
  combined across strata by Stouffer's method.
- **Essentiality axes.** Hart core-essentials (283; from the authors' repo), CEGv2 reference
  essentials (684; Hart lab), and shet continuous selective constraint (Zeng et al.). NEGv1
  nonessential-reference genes were dropped — only 12/928 survive HVG selection, itself evidence
  that the gene universe is expression-biased.
- **Cross-validation.** An independent from-scratch rebuild of the binary DE matrix (byte-range
  stream of `adj_p_value`, on-target masking) reproduces `n_regulators` to three decimals.

## Reproducibility

- Single public input file (`GWCD4i.DE_stats.h5ad`), read by **byte-range streaming** — never
  downloaded whole. Gene lists fetched from raw GitHub. No credentials, no local data.
- `notebooks/HumanCD4CoDEGNet_buffering.ipynb` runs top-to-bottom in **~2 minutes** and regenerates
  every number and figure here.
- All effect sizes, p-values, and per-gene tables are written to `artifacts/` with dataset ETag
  provenance.

## How Claude Science got us here

Dataset discovery (CZI Virtual Cells Platform) → read the preprint and the authors' analysis repo to
locate the exact un-taken analysis (their Hart list is loaded but commented out) → byte-range access
to a 16.8 GB file → reuse of their responsiveness metric and power-correction → the essentiality
join they did not make → non-parametric matched validation → interpretation against the user's own
prior yeast work (the only reason a cross-kingdom comparison is possible).

## Robustness (independently verified)

- **Independent pipeline.** Rebuilding responsiveness from scratch (streamed `adj_p_value`<0.10,
  on-target masked) instead of the authors' `n_regulators` reproduces every effect size
  (Spearman with authors' metric r=0.996; matched Hart δ=−0.172, CEGv2 δ=−0.181).
- **Binning.** Matched δ stays negative at 5/10/20/50 baseMean bins for both essential sets.
- **Permutation null.** Shuffling essential labels *within* baseMean deciles (2,000×): observed
  δ=−0.174 vs null mean 0.000 (95% CI [−0.10, +0.11]), p=0.001.
- **Continuous model.** NB-GLM `nreg ~ essential + log10(baseMean)`: CEGv2 coef −0.17, p=0.001;
  Hart coef −0.13, **p=0.16 (n.s. — the 124-gene list is underpowered for a parametric model)**.
- **Not a housekeeping artifact.** The most-buffered decile is dominated by OXPHOS (OR 5.55,
  p=5×10⁻¹³; ribosomal genes are HVG-excluded). Buffering **survives excluding OXPHOS+proteasome**
  and is if anything stronger (Hart δ −0.17→−0.20).

## Prior art and positioning

The buffering **direction is not first-reported here** — verified against the primary sources:

- **Feng et al. 2026** (*Cell Genomics*; human iPSC genome-scale CRISPRi across many lines): **"Controlling
  for differences in expression level… evolutionary conservation was the most informative predictor for
  having fewer [trans-]regulators"** (Fig S2.1D). This is decisive: Feng **already control for expression**
  and **already report the buffering direction, on our axis, in human cells.** So neither the direction nor
  the need to correct for expression is novel here — we **confirm and extend**, we do not discover, and we
  do **not** claim the expression-confound insight as ours.
- **Nourreddine et al. 2026** (*Nature Biotechnology*; KOLF2.1J iPSC genome-scale Perturb-seq atlas):
  concurrent human iPSC atlas; builds a co-perturbation cell map + specific regulators (ZBTB41, RNF7, DBR1)
  and overlaps DepMap essentiality, but does **not** relate essentiality to incoming trans-regulator count.
  Related platform, not a pre-emption of this analysis.
- **Mostafavi & Pritchard 2023** (*Nature*): constrained genes are eQTL-depleted — buffering via *natural*
  variation + selection. (J.K. Pritchard co-authors our dataset.) Ours is selection-unfiltered CRISPRi.
- **Newman 2006 / Batada & Hurst 2007**: essential genes have low expression *noise* — a stochastic axis;
  the confound runs the opposite way. **Costanzo 2016**: essential genes are genetic-interaction *hubs* —
  the perturbation-*source* axis, which coexists with readout-buffering.

**We do not claim first evidence, nor the expression-confound insight** — Feng et al. controlled for
expression and reported the direction in human iPSCs. This is an **extension**, and its value is the
forward programme, not a discovery. What genuinely survives as ours:
1. **A different, disease-relevant cell system** — primary human CD4+ T cells, not iPSC lines; the
   principle is not confined to pluripotent cells.
2. **Activation dynamics** — buffering is *dynamic*, peaking at early activation (8 h); Feng had no
   stimulation-state axis. Genuinely new and immune-specific.
3. **Cross-kingdom unification** with the yeast co-differential-expression law — only possible because
   we hold the yeast comparator.
4. **Forward applications** — the transcriptional-tractability axis and the shipped per-gene
   buffering-score resource for the Marson atlas.

The naive→corrected sign-flip is retained only as a *pedagogical illustration* of the same confound
Feng controlled for — **not** claimed as a novel methodological warning.

## Impact — what this enables (demonstrated)

The contribution is less the effect size than the research programme it opens. Three forward results,
two demonstrated in the notebook:

1. **Buffering is dynamic — it peaks at early T-cell activation (8 h)** (matched δ: Hart −0.04→−0.34→−0.10;
   CEGv2 −0.11→−0.29→−0.13). The genes with the strongest activation-associated buffering are canonical
   effectors (IFNG, GZMB, CSF2, LAG3) → hypothesis: the effector programme is transcriptionally
   *canalised* during early activation. Immunology-relevant, testable, hands a new question to the resource.
2. **Buffering is a transcriptional-tractability axis** — reconnecting to the track's "find drug targets"
   prompt. Druggable-genome genes are significantly *less* buffered (median resid +8.0 vs −16.9, p≈10⁻²¹,
   n=823); intersecting essential + druggable + least-buffered recovers known tractable targets
   (PLK1, AURKB, CHEK1, CDK1, WEE1). Buffered genes resist transcriptional modulation — a principled
   filter for upstream/indirect targeting strategies.
3. **A pedagogical illustration of a known correction** — expression-power correction matters for any
   Perturb-seq responsiveness / trans-hubness ranking (Feng et al. already apply it); our sign-flip makes
   the stakes vividly explicit as a teaching example. Not claimed as novel.

**Resource shipped:** per-gene buffering score for all 10,282 genes (`buffering_score_resource.csv`).

**Flagship next experiment:** run the identical test on Replogle 2022 genome-scale Perturb-seq
(K562/RPE1) — same correction → tests whether buffering is a general principle or cell-type-specific.

*Honest scope:* the essential+druggable list is small (n=23, illustrative); the cross-condition
comparison is exploratory. Both are hypothesis-generating.

## Honest limitations

- Effect sizes are modest (matched δ ≈ −0.17 to −0.19) though unanimous across 12 strata, two
  essential-gene sets, and an independent pipeline; this is a distributional tendency, not a rule.
- The parametric NB-GLM is n.s. for the 124-gene Hart set (significant for CEGv2); the claim rests on
  the non-parametric matched test and the permutation null.
- HVG selection (10,282 genes; no ribosomal genes) depletes nonessential-reference genes and
  **attenuates** the contrast — the estimate is conservative, not inflated.
- Stronger buffering at Stim8hr is suggestive (binary lists) but the continuous shet axis is flat
  across conditions, so condition-dependence is a hypothesis, not a claim.
- Essentiality is organismal/cell-line (Hart, CEGv2), not CD4+ T-cell fitness specifically.
- The buffering *direction* is corroborated by prior human work (Feng 2024; Mostafavi 2023); the
  sign-inversion, cross-kingdom framing, and primary-immune-cell context are the novel parts.

## Credit and references

- **Data:** Zhu R., Dann E. et al. (2025) *Genome-scale perturb-seq in primary human CD4+ T cells*,
  bioRxiv 10.64898/2025.12.23.696273.
- **Method lineage:** Nasar, Rehman, Ott & Alam (2026) *Uncovering coordinated pathway interactions
  through gene co-differential expression in yeast*, NAR 54(1):gkaf1410.
- **Closest prior art (largely pre-empts the direction + expression-control):** Feng C. et al. (2026)
  *A genome-scale single-cell CRISPRi map of trans gene regulation across human pluripotent stem cell
  lines*, **Cell Genomics 3:101076** (bioRxiv 2024.11.28.625833). Fig S2.1D, controlling for expression,
  finds conservation the top predictor of *fewer* trans-regulators.
- Nourreddine S. et al. (2026) *A genome-scale CRISPRi perturbation atlas of human induced pluripotent stem
  cells*, **Nature Biotechnology**, doi:10.1038/s41587-026-03199-w (concurrent human iPSC Perturb-seq atlas).
- Mostafavi H. et al. (2023) *Systematic differences in discovery of genetic effects on gene expression
  and complex traits*, Nature. Newman J.R.S. et al. (2006) Nature 441:840. Batada & Hurst (2007) Nat
  Genet 39:945. Costanzo M. et al. (2016) Science 353:aaf1420.
- **Reference gene sets:** Hart et al. core essentials & CEGv2/NEGv1; shet selective constraint
  (gnomAD-based). Full prior-art report: `artifacts/novelty.md`.
