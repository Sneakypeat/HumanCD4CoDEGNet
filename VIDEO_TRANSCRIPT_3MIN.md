# Video transcript — 2.5–3 min (figure-cued)
*(~440 words at ~150 wpm ≈ 2:55. **[SHOW]** cues indicate which figure to display on screen.)*

---

**[0:00 — Hook | title card]**

A human T cell has to do something almost paradoxical. To fight an infection, it rewires thousands of
genes in a matter of hours — a violent reorganization. And yet it stays a stable, functioning cell the
whole time. So here's the question: when the wiring changes that much, what stays the same? Is there
anything constant about how a T cell is *controlled*?

**[0:22 — The data]**

To answer that you can't just watch genes correlate — you have to perturb them and see what moves. The
Marson lab just released exactly that: a genome-scale Perturb-seq atlas. Twenty-two million primary
human CD4 T cells, every expressed gene individually knocked down, measured at rest and after
activation. For the first time we can draw the *causal* wiring diagram of a human immune cell — and
redraw it in each state.

**[0:48 — The architecture | SHOW Figure 1]**

So we drew it. Every knockdown that moves a gene is one causal edge. And the network has a very
particular shape. It's **hub-dominated**: the top five percent of regulators drive nearly eighty percent
of all the causal effects. Most genes do almost nothing when you knock them down. Control is broadcast
by a tiny minority and received broadly by everyone else. And here's the striking part — this is exactly
the shape a theory paper from the Pritchard lab predicted last year, from heritability data alone, with
no perturbation experiment. This is the first time it's been confirmed causally.

**[1:30 — The twist | SHOW Figure 2]**

Now the new part. When the cell activates, does this architecture change? The answer splits in two. The
**shape** does not change — the concentration stays pinned at the same value, at rest, at eight hours, at
forty-eight hours. Rock steady. But the **identity** of the hubs changes almost completely — up to
fifty-nine percent of the top hubs are replaced. The T-cell receptor signalosome — genes that move a
handful of targets at rest — come online and each command thousands. The cell keeps the *shape* of its
control while swapping out *who* is in control.

**[2:12 — Generality | SHOW Figure 3]**

And it isn't only T cells. Run the same analysis on a completely different cell type — a cancer line,
a different lab, a different platform — and the same hub-dominated shape appears. The architecture is a
general property of how cells are wired.

**[2:32 — Honesty and close]**

Why do I trust this? It isn't a fragile mechanism — it's a distribution-level invariant, and it survives
every confound we could throw at it: sequencing power, knockdown efficiency, detectability. The whole
thing reproduces in seconds from two public summary tables, code on GitHub. It's not a drug. It's a law
about how a human immune cell stays itself while reinventing itself. Thanks for watching.
