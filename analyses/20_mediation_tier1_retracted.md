# 20. Mediation tier-1 — "88% mediated" (RETRACTED, recast)

**Question.** Of CSP1-2's apparent effect on community diversity, what fraction is mediated by intermediate biological pathways?

**Method.** Per-compartment mediation with NDVI / soil-moisture / temperature as mediators. `scripts/run_mediation_per_compartment.py`, `scripts/run_mediation_robustness.py`, `scripts/run_mediation_sensitivity_full.py`. Outputs in `cache/causal_tier1_*.tsv`.

**Original result.** "88% of CSP1-2's effect on community diversity is mediated through NDVI" — headline number used in early drafts.

**Retraction (2026-05-08).** Mediation analysis was overstated:
1. The 88% number was unstable across alternative mediator orderings and sensitivity parameters.
2. The relationship is **non-linear** (Hill-functional, not linear-mediation).
3. Subsequent analysis (#33, #46) shows CSP1-2's apparent diversity effect is partly a **relic-DNA artifact** (#28, #32).

**Recast.** Replaced with **direct-on-guild + non-linear Hill via CSP1-2** framing. Mediation magnitude was dropped from main paper; mechanistic single-keystone claim was replaced by guild-level framing (#38, #52).

**Status.** retracted; replaced by guild/Hill framing.

**Outputs (kept for reproducibility).**
- `cache/causal_tier1_mediation.tsv`
- `cache/causal_tier1_mediation_robustness.tsv`
- `cache/causal_tier1_mediation_sensitivity{,_full}.tsv`
- `cache/causal_tier1_dml_ate.tsv`
- `cache/causal_tier1_panel_fe.tsv`
- `cache/causal_frame_tier1.parquet`

**Cross-refs.** 38, 32, 52.
