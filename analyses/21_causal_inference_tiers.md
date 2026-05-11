# 21. Causal inference — Tier 2 structure learning + Tier 3 interventions

**Question.** Beyond mediation, what causal-graph structure can be inferred from the data, and what would happen under hypothetical interventions?

**Method.**
- **Tier 2** — structure learning via PC, LiNGAM, FCI on per-site climate + diversity + chemistry + dominance variables. `cache/causal_tier2_pc_edges.tsv`, `cache/causal_tier2_lingam_edges.tsv`, `cache/causal_tier2_fci_edges.tsv`, `cache/causal_tier2_fci_stability.tsv`.
- **Tier 3** — counterfactual interventions: vary one variable, predict outcomes from fitted SCM. `scripts/run_cmip6_projections.py`, `cache/causal_tier3_interventions{,_per_compartment,_expanded}.tsv`, `cache/cmip6_interventions.tsv`.

**Key results.**
- PC/LiNGAM/FCI all agree on a small backbone: climate → chemistry → diversity, with weak/uncertain edges into CSP1-2.
- CMIP6 intervention table is the **v1 climate projection** (#43) — uses cross-sectional features and is **artifactual**.

**Status.** Tier 2 structure-learning results are solid as exploratory structure; Tier 3 was superseded by the longitudinal projections in #43–45.

**Cross-refs.** 20, 43, 44, 45, 46.
