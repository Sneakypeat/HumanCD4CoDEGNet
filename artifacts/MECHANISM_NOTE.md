# Retired mechanism — "activation edges route around essential genes"

**Status: REFUTED at the edge level. Not claimed anywhere in the submission.**

During exploration we hypothesized that the activation-induced densification of the trans-regulatory
network *spares* essential genes — i.e., that the new regulatory edges preferentially avoid core-essential
targets, producing a "division of labor" (build the response program, canalize the essential one).

We tested this directly (`artifacts/essential_avoidance.json`). It does not hold:

- Genome-wide essential fraction: **0.0367**. Fraction of *newly-added* activation edges landing on
  essential targets: **0.0456** — a ratio of **1.24× baseline** (essential genes get *more* of the new
  edges, not fewer).
- Binomial test for *avoidance* (fewer than baseline): p = 1.0 (i.e., no avoidance).
- Fisher OR = 1.25, p = 1.1×10⁻⁵ in the direction *opposite* to the hypothesis.
- **Verdict: NO.** The network does not route around essential genes.

We therefore **retired this mechanism** and did not build the submission on it. The essential-gene
*buffering* phenomenon itself is real and is reported only as a confirm-and-extend secondary result
(matching Feng et al. 2026); the *mechanistic* story about how activation produces it is withdrawn.

This file is kept for transparency and to prevent the refuted claim from being reintroduced. The flagship
finding (network architecture: shape-invariant, identity-labile) is independent of this and does not rely
on any essentiality claim.
