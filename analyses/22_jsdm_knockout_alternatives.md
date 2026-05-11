# 22. JSDM knockout — alternative-knockout sensitivity

**Question.** If we knock out CSP1-2 (or any keystone) in a fitted joint species distribution model, what is the predicted system response — and how much depends on the choice of focal taxon?

**Method.** Fit Bayesian JSDM (HMSC / sjSDM-style), simulate per-sample abundance with the focal taxon set to zero. Run for CSP1-2 and a battery of alternative knockouts. `cache/jsdm_knockout_*.tsv`.

**Key results.**
- JSDM perturbation predicts modest effects when CSP1-2 is removed.
- **Alternative knockouts** (Nibribacter, Massilia, Aquibacillus, Halomonas, …) produce comparable or larger effects.
- Topological centrality alone does not single out CSP1-2 as uniquely important.

**Interpretation.** This was the early signal that the single-keystone narrative was fragile. Picked up later by the explicit knockout test (#38).

**Status.** solid as a sensitivity check.

**Cross-refs.** 38 (later, decisive), 53.
