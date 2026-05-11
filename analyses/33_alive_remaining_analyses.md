# 33. Remaining analyses on alive subset (iCAMP, mediation, wind-Mantel, distance-decay)

**Question.** Do other macroscopic findings survive the relic filter?

**Method.** Re-run each canonical analysis on the alive feature table. `scripts/relic_alive_only_analyses.py`, `scripts/relic_alive_remaining_analyses.py`, `scripts/relic_alive_only_functions.py`.

**Findings (one line each).**

| Analysis | All-data signal | Alive-only signal | Verdict |
|---|---|---|---|
| iCAMP homogenizing dispersal | ~67% | survives (~62%; mildly weaker) | preserved |
| Wind-dispersal Mantel | r = 0.05–0.34 | r ≈ 0.10–0.15 | weakened — wind operates more on relic |
| Distance-decay (taxonomic) | clear slope | preserved | OK |
| Functional redundancy (Allison-Martiny slope) | 0.21 | **0.42** | **artifact** — doubling reverted (per #28) |
| Mediation (CSP1-2 → diversity) | 88% (#20) | not replicated | retracted |
| Climate-Shannon gradient (#30) | ρ = −0.40 | **ρ ≈ 0** | RELIC ARTIFACT |
| Temporal stability (#31) | moderate | **2× tighter** | alive more stable |

**Interpretation.** Three categories emerge: (1) wind dispersal survives but is **more pronounced in relic** (wind moves dead biomass); (2) climate-diversity is purely relic; (3) functional redundancy doubling was an artifact of how the relic mask was applied.

**Status.** solid.

**Cross-refs.** 6, 7, 12, 30, 31.
