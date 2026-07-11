# Video transcript — 2.5–3 min (figure-cued, payoff-first)
*(~415 words at ~150 wpm ≈ 2:45. **[SHOW]** cues indicate which figure to display on screen.)*

---

**[0:00 — Hook: the payoff]**

When a human T cell switches on to fight an infection, a specific set of control genes comes online — and
they turn out to be disproportionately *druggable*, and disproportionately genes that cause disease when
they break. Genes like ZAP70 and ITK. We found them by mapping the entire **causal control system** of a
human immune cell and reading off where disease risk concentrates. Here's how.

**[0:28 — The data]**

You need causality, not correlation — you have to perturb genes and watch what moves. The Marson lab just
released exactly that: a genome-scale Perturb-seq atlas. Twenty-two million primary human CD4 T cells,
every expressed gene individually knocked down, at rest and after activation. For the first time we can
draw the causal wiring diagram of a human immune cell — and redraw it in each state.

**[0:52 — The architecture | SHOW Figure 1]**

So we drew it. The network has a very particular shape: it's **hub-dominated**. The top five percent of
regulators drive nearly eighty percent of all causal effects; most genes do almost nothing when you knock
them down. That's exactly the shape a Pritchard-lab theory paper predicted last year from heritability
alone — confirmed causally here for the first time.

**[1:25 — The twist | SHOW Figure 2]**

When the cell activates, does the architecture change? It splits in two. The **shape** doesn't change — the
concentration stays pinned at rest, eight hours, forty-eight hours. But the **identity** of the hubs turns
over almost completely — up to fifty-nine percent replaced. The cell keeps the shape of its control while
swapping out who's in control. And the new hubs are the T-cell receptor signalosome — ZAP70, LAT, ITK —
switching on.

**[2:00 — The payoff delivered | SHOW Figure 5]**

Here's why that matters. Those activation-specific hubs are about **twice as enriched for genes that cause
disease when mutated** as the stable hubs — and they're druggable: ZAP70, ITK, LCK, PTPRC, IL12RB2. The
control layer that activation builds is exactly where druggable disease risk concentrates — candidate
targets that only matter in the activated state.

**[2:30 — Generality and close | SHOW Figure 3]**

It's not a T-cell quirk — the same hub-dominated architecture appears in a completely different cell type.
And it survives every confound we tested: a distribution-level law, not a fragile story, reproducible from
public data, code on GitHub. It's not a drug yet — but it's a map of how a human immune cell controls
itself, and a shortlist of where to intervene. Thanks.
