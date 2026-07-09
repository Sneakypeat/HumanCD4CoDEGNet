# 3-minute video transcript
*(~450 words, ~150 wpm. Cue timings are guides.)*

---

**[0:00 — Hook]**

A human T cell has to do something almost paradoxical. To fight an infection it rewires thousands of
genes in a matter of hours — a violent reorganization. And yet it stays a stable, functioning cell the
whole time. So here's the question: when the wiring changes that much, what stays the same? Is there
anything constant about how a T cell is *controlled*?

**[0:25 — The data]**

To answer that you need to actually perturb genes and watch what moves — you need causality, not
correlation. The Marson lab just released exactly that: a genome-scale Perturb-seq atlas. Twenty-two
million primary human CD4 T cells, every expressed gene individually knocked down, measured at rest and
after activation. For the first time we can draw the *causal* wiring diagram of a human immune cell —
and redraw it in each activation state.

**[0:55 — The architecture]**

So we drew it. Every validated knockdown that moves a gene is one causal edge. And the network has a
very particular shape. It's **hub-dominated**: the top five percent of regulators drive nearly eighty
percent of all the causal effects. Most genes do almost nothing when you knock them down — a typical
gene is controlled by less than one percent of its possible regulators. Control is broadcast by a tiny
minority and received broadly by everyone else. Strikingly, this is exactly the shape a theory paper
from the Pritchard lab predicted last year — from twin-study heritability alone, with no perturbation
data. This is the first time anyone has confirmed it causally.

**[1:40 — The twist]**

Now the part that's new. We asked: when the cell activates, does this architecture change? And the
answer is beautifully split in two. The **shape** does not change — the concentration, the
hub-dominance, is pinned at the same value, a Gini of zero-point-nine-two, at rest, at eight hours, at
forty-eight hours. Rock steady. But the **identity** of the hubs changes almost completely — up to
fifty-nine percent of the top hubs are replaced. The T-cell receptor signalosome — CD3, ZAP70, LAT —
genes that control just a handful of targets at rest, come online and each command thousands.

**[2:20 — The one-liner]**

So the cell does something elegant: it **preserves the shape of its control architecture while
swapping out who is in control.** A constant grammar, a rewritten vocabulary.

**[2:40 — Honesty and close]**

And I want to be clear about why I trust this. It's not a fragile mechanism — it's a distribution-level
invariant, and it survives every confound we could throw at it: sequencing power, knockdown efficiency,
detectability. The whole thing reproduces in seconds from two public summary tables, code on GitHub.
It's not a drug. It's a law about how a human immune cell stays itself while reinventing itself. Thanks
for watching.
