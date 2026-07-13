# Video script — with self-introduction (first person, < 3 min)
*(~375 words ≈ 2:45 at ~150 wpm. **[SHOW]** cues indicate what's on screen. A self-advancing
teleprompter version of this script is in `video/presenter.html`.)*

---

**[0:00 — Intro]**

Hi, I'm **Syed Sabih ur Rehman, a PhD student at United Arab Emirates University**. My project is
HumanCD4CoDEGNet.

**[0:12 — Hook: the payoff]**

When a human T cell switches on to fight an infection, a specific set of control genes comes online — and
they turn out to be disproportionately *druggable*, and disproportionately genes that cause disease when
they break. Genes like ZAP70 and ITK. I found them by mapping the entire **causal control system** of a
human immune cell. Here's how.

**[0:38 — The data | SHOW poster top]**

You need causality, not correlation — you perturb genes and watch what moves. The Marson lab released
exactly that: a genome-scale Perturb-seq atlas — twenty-two million primary human CD4 T cells, every gene
knocked down, at rest and after activation. For the first time we can draw the causal wiring diagram of a
human immune cell, and redraw it in each state.

**[1:00 — The architecture | SHOW Figure 1]**

So I drew it. The network has a very particular shape: it's **hub-dominated**. The top five percent of
regulators drive nearly eighty percent of all causal effects; most genes do almost nothing when knocked
down. That's exactly the shape a Pritchard-lab theory paper predicted last year from heritability alone —
confirmed causally here for the first time.

**[1:30 — The twist | SHOW Figure 2]**

When the cell activates, does the architecture change? It splits in two. The **shape** doesn't change — the
concentration stays pinned at rest, eight hours, forty-eight hours. But the **identity** of the hubs turns
over almost completely — up to fifty-nine percent replaced. The cell keeps the shape of its control while
swapping out who's in control. The new hubs are the T-cell receptor signalosome — ZAP70, LAT, ITK —
switching on.

**[2:05 — The payoff delivered | SHOW Figure 5]**

Here's why that matters. Those activation-specific hubs are about **twice as enriched for genes that cause
disease when mutated** as the stable hubs — and they're druggable: ZAP70, ITK, LCK, PTPRC, IL12RB2. The
control layer that activation builds is exactly where druggable disease risk concentrates — targets that
only matter in the activated state.

**[2:35 — Close | SHOW Figure 3]**

It's not a T-cell quirk — the same architecture appears in a completely different cell type. And it
survives every confound I tested — reproducible from public data, code on GitHub. It's a map of how a
human immune cell controls itself, and a shortlist of where to intervene. Thanks for watching.
