# Novelty / prior-art assessment — transcriptional buffering of essential genes

**Claim vetted:** In primary human CD4+ T cells (genome-scale CRISPRi Perturb-seq), essential /
selectively-constrained genes are transcriptionally *buffered* — differentially expressed by FEWER
trans-acting perturbations than expected for their expression level. A raw analysis shows the OPPOSITE
(essential genes look MORE responsive) because they are highly expressed and better-powered for DE
detection; correcting for expression inverts the sign. Generalises a yeast co-DE "buffering law".

**Verdict: PARTIALLY NOVEL.** The *biological direction* (constrained genes have fewer trans-regulators)
is already reported in human cells — most damagingly by Feng et al. in hPSCs — and in humans via eQTL
depletion. What survives as novel and defensible: (i) primary human immune cells across activation states;
(ii) the specific, quantified **sign-inversion** driven by the expression-power confound; (iii) the explicit
yeast→human cross-kingdom framing.

---

## Closest prior works (ranked by threat)

### 1. Feng et al. — genome-scale CRISPRi trans-regulation map in human iPSCs  *(CLOSEST — partial pre-emption)*
- Feng C. et al. (2024/2025) *A genome-scale single-cell CRISPRi map of trans gene regulation across human
  pluripotent stem cell lines.* bioRxiv 2024.11.28.625833 → Cell Genomics 2025.
- URL: https://www.biorxiv.org/content/10.1101/2024.11.28.625833v1 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC12903452/
- **Relevance (high):** Same platform (human, genome-scale, single-cell CRISPRi, *trans* readout). They
  report **evolutionary conservation is the most informative predictor of a gene having FEWER trans-regulators**,
  and they **control for statistical power** (cells per target). That is our buffering direction, power-corrected,
  already in a human system. Separately they find target **essentiality predicts having MORE trans effects** — but
  that is the *perturbation-SOURCE* axis (essential knockdowns broadcast widely), NOT the *readout* axis we claim.
- **What they did NOT do:** primary immune cells / activation-state contrast; and they used a joint regression
  (power as covariate) — they did not exhibit the **marginal naive→corrected SIGN REVERSAL**, nor a cross-kingdom law.

### 2. Nasar, Rehman, Ott & Alam (2026) — the yeast buffering law being generalised  *(self-citation)*
- *Uncovering coordinated pathway interactions through gene co-differential expression in yeast.* NAR 54(1):gkaf1410.
- URL: https://academic.oup.com/nar/article/54/1/gkaf1410/8415822
- **Relevance:** The prior result the human claim generalises — "essential genes exhibit significantly lower
  responsiveness than nonessential genes (P<3.5e−112)." Confirmed by fetch: the yeast paper does **NOT** discuss the
  expression-power confound. So the confound-correction + human generalization is additive, not duplicative. This is
  the user's own method lineage (not an independent pre-emption).

### 3. Mostafavi … Pritchard (2023) — constrained genes are depleted for eQTLs  *(human, natural-variation route)*
- Mostafavi H. et al. (2023) *Systematic differences in discovery of genetic effects on gene expression and complex
  traits.* Nature (bioRxiv 2022.05.07.491045). Related: Glassberg et al. 2019, *Evidence for weak selective constraint
  on human gene expression*, Genetics 211:757.
- URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC12270542/ ; https://web.stanford.edu/group/pritchardlab/publications/pdfs/Glassberg19.pdf
- **Relevance (high, conceptual):** LoF-intolerant / high-pLI genes are depleted for eQTLs and show low allelic
  imbalance — i.e. constrained genes are *buffered against genetic regulatory variation* in humans. **Mechanistic
  distinction that protects us:** this is *natural* variation filtered by **stabilizing selection** ("eQTLs with large
  effects on constrained genes are purged") — an ascertainment/selection story. Our CRISPRi perturbations are
  *experimentally imposed and unfiltered by selection*, so buffering reflects intrinsic regulatory architecture, not
  the removal of variants. **NB: J.K. Pritchard co-authors both this paper and our CD4 dataset** → a reviewer from that
  lab will know this framework; must cite and differentiate.

### 4. Newman 2006 / Batada & Hurst 2007 / Lehner 2008 — essential genes have LOW expression NOISE  *(different axis)*
- Newman J.R.S. et al. (2006) Nature 441:840 (yeast GFP noise); Batada & Hurst (2007) Nat Genet 39:945; Lehner (2008).
- URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC2290932/ ; https://www.nature.com/articles/ng2071
- **Relevance / distinction:** "Low noise" = stochastic cell-to-cell variation under **constant** conditions
  (intrinsic+extrinsic). Our metric = count of **trans-genetic perturbations** that move the gene — a directed,
  causal regulatory-wiring quantity, not stochasticity. Crucially the confound runs the **opposite** way: for noise,
  high expression → *low* noise (aligned, makes essentials look quiet); for DE-power, high expression → *more* detected
  DE (a confound that must be **removed** to reveal buffering). Not a pre-emption, but the nearest classical relative.

### 5. DE-power / dosage-constraint confound + Costanzo hubness  *(the confound is textbook; hubness is source-axis)*
- "Recalibrating differential gene expression by genetic dosage variance prioritizes functionally relevant genes"
  (bioRxiv 2024.04.10.588830 / PMC11030425): expression level & dosage constraint bias DE detection — highly
  expressed genes are easier to call DE; constrained genes are disadvantaged.
  URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC11030425/
- Costanzo M. et al. (2016) Science 353:aaf1420: essential genes are genetic-interaction **hubs** (~5× more
  interactions). URL: https://pubmed.ncbi.nlm.nih.gov/27708008/
- **Relevance:** The expression→DE-power confound is a *known* phenomenon (so the confound's existence is not novel;
  its *severity — reversing the essentiality sign* — is our specific contribution). Costanzo's hubness is the
  **source** axis; it coexists with (does not contradict) readout-buffering, and makes the dissociation more striking.

---

## What is NOVEL vs PRE-EXISTING

**Pre-existing (do NOT claim as new):**
- Constrained/essential genes are transcriptionally buffered / have fewer trans-regulators — shown in yeast
  (Nasar 2026) AND in human cells via eQTL depletion (Mostafavi/Pritchard 2023) AND via genome-scale CRISPRi in
  hPSCs where conservation predicts fewer trans-regulators with power controlled (Feng et al.).
- Essential genes are highly expressed and low-noise (Newman 2006; Batada & Hurst 2007).
- Expression level is a power confound for DE detection; constrained genes lose out in fold-change DE (dosage-recalibration paper).
- Essential genes are network hubs as perturbation *sources* (Costanzo 2016; Feng "essentiality → more trans effects").

**Novel / defensible in this submission:**
1. **Primary human immune cells across activation states** (resting + stimulation) — Feng was hPSC lines; distinct,
   context-specific biology.
2. **The quantified naive→corrected SIGN INVERSION** — raw responsiveness is not merely attenuated but *reversed*
   (δ=+0.26 → δ≈−0.17); no prior source demonstrates that the expression-power confound flips the essentiality sign.
   This reconciles the apparent yeast–human discrepancy and is a genuine methodological cautionary result.
3. **Explicit cross-kingdom unification** of the co-DE buffering law (yeast→human, ~1 Gyr), with the confound as the
   reconciling mechanism, and validated on two essential-gene sets + continuous shet across three states.
4. **Experimental (selection-unfiltered) evidence** rules out the stabilizing-selection/ascertainment explanation that
   underlies the eQTL route (Mostafavi) — a mechanistic complement, not a duplicate.

**Interpretive frame worth stating:** essential genes are buffered as *readouts* yet act as hubs as *sources*
(Costanzo/Feng) — a clean dissociation the field has not articulated together.

## Defensible positioning (one line)
Do NOT claim first evidence that human constrained genes are buffered — cite **Feng et al. (hPSC; conservation→fewer
regulators, power-controlled)** and **Mostafavi/Pritchard (eQTL depletion)** as independent corroboration, and pitch the
contribution as (a) primary human T cells + cross-kingdom generalization and (b) the confound-driven **sign-inversion**
as a methodological warning that naive cross-species responsiveness comparison is actively misleading.

---
### Sources
- https://www.biorxiv.org/content/10.1101/2024.11.28.625833v1
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12903452/
- https://academic.oup.com/nar/article/54/1/gkaf1410/8415822
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12270542/
- https://web.stanford.edu/group/pritchardlab/publications/pdfs/Glassberg19.pdf
- https://pmc.ncbi.nlm.nih.gov/articles/PMC2290932/
- https://www.nature.com/articles/ng2071
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11030425/
- https://pubmed.ncbi.nlm.nih.gov/27708008/
- https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1
