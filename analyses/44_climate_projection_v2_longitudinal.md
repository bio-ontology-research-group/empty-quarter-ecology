# 44. Climate projection v2 — longitudinal (precip only)

**Question.** Properly: how does climate change shift A↔B transition rates, using **trip-to-trip transitions** as training data instead of cross-sectional clustering?

**Method.** 453 trip-to-trip transitions as training set. Two logits:
- P(A→B | Δclimate) on A-start transitions (389 transitions, 47 A→B events = 12.1%).
- P(B→A | Δclimate) on B-start transitions (64 transitions, 31 B→A events = 48.4%).

Features: Δd7, Δd30, Δd90, Δd180, Δd365 (precip windows). No temperature in v2. `scripts/two_strategy_climate_projection_v2.py`, `cache/two_strategy_projection_v2/`.

Apply CMIP6 (ΔT, ΔP_pct) → compute 2-state Markov equilibrium π_B = P(A→B) / (P(A→B) + P(B→A)).

**Key results — equilibrium B-fraction trajectory.**

| Scenario      | π_B   | Δ from historical |
|---------------|-------|-------------------|
| Historical    | 0.469 | —                 |
| SSP1-2.6_2100 | 0.455 | −1.4 pp           |
| SSP2-4.5_2100 | 0.445 | −2.4 pp           |
| SSP3-7.0_2100 | 0.420 | **−4.9 pp**       |

**Logistic coefficients for P(A→B).**
- Δd7: +0.747 (acute precip pulse pushes A→B)
- Δd365: +1.874 (chronic wet year pushes A→B)
- Δd180: −2.861 (intermediate window opposite sign)

**Interpretation (per v2 alone).** Drying SUPPRESSES the precipitation pulses that trigger A→B transitions → modest equilibrium drift TOWARD A. B is the **wet-pulse-favored state** (Bacilli/Halomonas exploit ephemeral water); A is the resilient drought default.

**v2's blind spot.** Transition table has no temperature ΔT — model uses **precip only**. Direction is **opposite to v1's artifact**, but the magnitude is small (~5 pp) and temperature pathway is missing.

**Status.** **superseded by v3** (#45). Direction is correct but incomplete.

**Outputs.**
- `cache/two_strategy_projection_v2/scenario_summary_v2.tsv`
- `cache/two_strategy_projection_v2/per_site_AtoB_risk.tsv`

**Cross-refs.** 41, 43, 45.
