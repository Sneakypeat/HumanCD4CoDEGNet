# Mechanism novelty / prior-art check — activation-induced canalization of the essential program

**Claim vetted (the DYNAMIC mechanism, distinct from the static-buffering claim in `novelty.md`):**
In primary human CD4+ T cells (genome-scale CRISPRi Perturb-seq), T-cell activation *reorganizes* the
trans-regulatory network — total trans edges rise ~32% (peak 8h) — and this new layer preferentially
targets **non-essential immune-response genes** while **sparing core-essential genes**, which instead
switch from activating to repressive trans-inputs and become transcriptionally **buffered / canalized**.
I.e. activation induces a *division of labor* that canalizes the essential program.

**VERDICT: PARTIALLY NOVEL.** Every *ingredient* has a strong precedent (rewiring on T-cell activation;
static buffering of essential/constrained genes; repression-as-canalization; housekeeping-robust vs
responsive-plastic division of labor). No single prior work I found makes the **integrated, activation-
INDUCED** claim — that the essential program's buffering is *created by the state transition* via
essential-sparing densification + an activating→repressive input switch. The defensible novel core is the
**dynamics/induction**, not the existence of buffering. Two same-system preprints (one on the *same
dataset*) already publish the weaker "activation reorganizes regulation" framing, so the contribution must
be positioned as a new structural analysis, not a restatement.

---

## Closest prior works (ranked by threat)

### 1. Zhu, Dann, … Pritchard & Marson (2025) — the SAME dataset  *(CLOSEST; shares the "activation reorganizes regulation" direction)*
- *Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific regulators of T cell
  programs and human immune traits.* bioRxiv 10.64898/2025.12.23.696273 (Rest / Stim-8h / Stim-48h; 22M cells).
- URL: https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1
- **Relevance (highest):** This is the dataset the mechanism is derived from. Their abstract states
  "active regulators and the gene programs they control **change dramatically across stimulation
  conditions**." That establishes context-specific regulatory turnover on activation — the *direction* of
  our claim. **What they do NOT claim (abstract + CZI dataset page):** any quantified edge-count
  densification (~32%/8h), any essential-vs-non-essential division of labor, any activating→repressive
  switch, or the words buffering/canalization. **Residual risk:** confirm no supplementary panel quantifies
  edges-by-condition or stratifies by essentiality (full text was Cloudflare-blocked; grep their PDF/supp
  for "essential", "canaliz", "buffer").

### 2. Mihai, Núñez-Carpintero, Cirillo … Selinger (2025) — "conserved backbone + rewiring" in activated CD4+ T cells  *(near-miss on phrasing)*
- *A conserved transcriptional backbone and rewiring of gene-regulatory networks in activated human CD4+
  T cells.* bioRxiv 10.1101/2025.11.25.687998 (multiome scRNA+scATAC GRN inference; Th1/2/17/Treg).
- URL: https://www.biorxiv.org/content/10.1101/2025.11.25.687998v1
- **Relevance:** Title reads like our "conserved core spared / rest rewired" division of labor. **But on
  inspection it is a different axis:** their "backbone" = **8 highly central TFs (regulators / SOURCE
  hubs)** conserved across differentiation subtypes — *not* essential TARGET genes (the readout axis we
  claim spared). Method is *inferred* GRNs from expression+accessibility, **not CRISPRi trans-perturbation**.
  No essentiality, buffering, canalization, or edge-densification quantification. Distinguishes cleanly.

### 3. Feng et al. (2025/2026, Cell Genomics; bioRxiv 2024.11.28.625833) — static trans-map in hPSCs  *(the STATIC pre-emptor)*
- *A genome-scale single-cell CRISPRi map of trans gene regulation across human pluripotent stem cell lines.*
- URL: https://www.biorxiv.org/content/10.1101/2024.11.28.625833v1 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC12903452/
- **Relevance:** Same platform class; reports **evolutionary conservation → FEWER trans-regulators, power-
  controlled** — i.e. the buffering direction already in a human system. **Key differentiator (given):** it
  is **static** (one iPSC state), *not activation dynamics*. It cannot show that buffering is *induced* by a
  state transition, which is our claim.

### 4. Trapotsi / Rahni / Loeb / Bar-Ziv … Wernig? "Canalizing cell fate by transcriptional repression" (2024, Mol Syst Biol)  *(repression→canalization, wrong target)*
- URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC10912439/
- **Relevance:** Establishes "**repressors canalize** by silencing alternative programs" — matches our
  "essential genes shift to repressive inputs → canalized" language. **But it is cell-FATE canalization**
  (Waddington landscape, alternative-lineage silencing), not buffering of an essential-gene program against
  trans-perturbation. The principle "repression→canalization" is therefore *not* ours to claim as new;
  applying it to activation-induced essential-gene buffering is. (Classic dev precedent: gap-gene mutual-
  repression canalization, Manu/Reinitz PLoS Biol 2009.)

### 5. Static division of labor: robust-housekeeping vs plastic-responsive architecture  *(the concept is old & static)*
- Sigalova & Furlong et al. (2023) *Promoter sequence and architecture determine expression variability and
  confer robustness to genetic variants.* eLife 80943 — housekeeping/robust vs signaling/response-to-stimulus/
  variable promoters. URL: https://elifesciences.org/articles/80943
- Nishiyama et al. (2009) — in ES cells, "**responsive genes** have regulatory functions & become tissue-
  specific; **non-responsive genes have housekeeping functions**" (housekeeping buffered vs TF manipulation).
  URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC3044670/
- **Relevance:** The *static* division of labor (essential/housekeeping = robust/buffered; stimulus-responsive
  = plastic/targeted) is well established at promoter-architecture and TF-responsiveness level. So the
  *existence* of the division is not novel; its **activation-INDUCED, trans-network form** (essential-sparing
  densification) is what remains.

### Adjacent / supporting (not pre-emptions)
- Freimer et al. (2022, Nat Genet, s41588-022-01106-y) — Perturb-seq "architecture of immune networks" in
  human T cells (targeted ~24 regulators; co-functional modules) — establishes Perturb-seq immune-network
  mapping but not essential-sparing/densification/canalization. https://www.nature.com/articles/s41588-022-01106-y
- Schmidt et al. (2022, Science, abj4008) — genome-wide CRISPRa/i **stimulation-response** screens in T cells
  (cytokine readouts, not trans-transcriptome). https://www.science.org/doi/10.1126/science.abj4008
- "Negative feedback buffers effects of regulatory variants" (PMC4332157) — essential genes have stronger
  negative feedback / trans-buffering (static, natural variation). https://pmc.ncbi.nlm.nih.gov/articles/PMC4332157/
- Amit et al. (2009, Science) — dynamic transcriptional network on pathogen stimulation in dendritic cells
  (foundational "activation drives a regulatory network"; no essential-sparing/canalization framing).
- Dynamic regulatory networks of T-cell trajectory (PMC8577129) — rewiring as a hallmark of state transition.

---

## NOVEL vs KNOWN — explicit

**Known / do NOT claim as new:**
- Activation *reorganizes/rewires* T-cell regulation (Zhu/Pritchard/Marson 2025; Mihai 2025; Amit 2009; T-cell
  rewiring literature). The *direction* "activation changes the network" is public, including on this dataset.
- Essential / evolutionarily-constrained / housekeeping genes are transcriptionally **buffered / less
  regulated** — but STATICALLY (Feng hPSC, power-controlled; Nishiyama ES cells; negative-feedback buffering;
  eQTL depletion — see `novelty.md`).
- **Repression canalizes** gene expression / cell fate (Mol Syst Biol 2024; gap-gene mutual repression).
- Static **division of labor**: housekeeping-robust vs stimulus-responsive-plastic promoter architecture
  (Sigalova/Furlong 2023; TFIID-vs-SAGA classic).

**Novel / defensible in this submission:**
1. **Induction by a state transition.** The buffering/canalization of essentials is presented as *created by
   activation*, not an intrinsic static property — no prior work shows a cell-state transition *inducing*
   essential-gene canalization. Feng (the static pre-emptor) explicitly does not do activation dynamics.
2. **Quantified, essential-*sparing* densification.** ~32% trans-edge increase peaking at 8h that spatially
   *avoids* essential genes. "Active regulators change" (dataset paper) is not a net-edge-count increase, and
   nobody stratifies that increase by essentiality.
3. **The activating→repressive input switch on essentials** as the mechanism of the canalization — a specific,
   directional, testable structural claim, distinct from generic "repression canalizes."
4. **Integration into one "division of labor" statement** for the resting→activated transition in primary human
   immune cells, tying densification (on non-essentials) to canalization (of essentials).

---

## Defensible positioning (one line)
Position the claim as the **dynamic INDUCTION** of essential-gene canalization by the activation state-
transition — a quantified, essential-*sparing* densification of the trans-network with an activating→repressive
input switch on essentials — and cite Zhu/Pritchard/Marson 2025 (context-specific regulators, *same data*),
Feng et al. (static buffering), and Mihai 2025 (backbone+rewiring, source-hub axis) as the priors it *extends*,
never as claims it restates. Do NOT claim first evidence of essential-gene buffering or of activation rewiring.

**Action before submission:** obtain full text of bioRxiv 2025.12.23.696273 (Cloudflare-blocked here) and
grep figures/supplement for "essential", "canaliz", "buffer", and any edges-per-condition panel, to be certain
the same-dataset paper did not quantify densification or essentiality stratification in a supplementary figure.

---
### Sources
- https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1
- https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- https://www.biorxiv.org/content/10.1101/2025.11.25.687998v1
- https://www.biorxiv.org/content/10.1101/2024.11.28.625833v1
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12903452/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC10912439/
- https://elifesciences.org/articles/80943
- https://pmc.ncbi.nlm.nih.gov/articles/PMC3044670/
- https://www.nature.com/articles/s41588-022-01106-y
- https://www.science.org/doi/10.1126/science.abj4008
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4332157/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC8577129/
