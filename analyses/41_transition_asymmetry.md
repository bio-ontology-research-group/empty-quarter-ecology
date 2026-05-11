# 41. Transition asymmetry — A→A vs B→B vs A→B vs B→A

**Question.** Are the two states equally stable, or does the system have a preferred attractor?

**Method.** For every (site, comp) cell, compute trip-to-trip transitions A→A, A→B, B→A, B→B. Aggregate. `scripts/transition_asymmetry.py`, `cache/transition_asymmetry/`.

**Key results — overall transition matrix (453 trip-pairs).**

| Transition | Count | Fraction of from-state |
|---|---|---|
| A → A | 342 | 88% |
| A → B | 47  | **12%** |
| B → A | 31  | 48% |
| B → B | 33  | 52% |

**Per-trip-pair transition rates.**

| Pair    | P(A→B) | P(B→A) |
|---------|--------|--------|
| T1→T3   | 0.114  | 0.286  |
| T3→T4   | 0.108  | 0.630  |
| T4→T5   | 0.172  | 0.500  |

**Per-cell sequence types.**

| Sequence | Count | Pattern |
|---|---|---|
| stable_A | 108 (66%) | never goes B |
| oscillating | 24 (15%) | 2+ direction changes |
| A_to_B_drift | 19 (12%) | directional drying drift |
| stable_B | 8 (5%) | never goes A |
| B_to_A_drift | 5 (3%) | directional wetting drift |

**Interpretation.**
- **Strategy A is the resilient default** (88% persistence per trip).
- **Strategy B is transient** (52% persistence; ~50% chance of reverting per trip).
- Drift is **asymmetric**: 19 cells drift A→B vs 5 cells drift B→A (≈4×) — net flow is A→B under drying.
- 24 cells truly oscillate (climate-sensitive boundary sites — the ones that will respond most under climate change).

**Status.** solid.

**Outputs.**
- `cache/transition_asymmetry/all_transitions.tsv`
- `cache/transition_asymmetry/per_cell_sequences.tsv`
- `cache/transition_asymmetry/per_cell_trip_dominant.tsv`

**Cross-refs.** 39, 40, 44, 45.
