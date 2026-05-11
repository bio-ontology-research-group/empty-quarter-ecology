# 15. Pulse–reserve precipitation alignment

**Question.** Do EQ communities respond to precipitation pulses in a Noy-Meir "pulse–reserve" pattern (rapid bloom after rain, slow decay)?

**Method.** Per-(site, comp), align trip dominance to recent-precip windows; fit lagged response curves. `scripts/run_pulse_reserve.py`, output `cache/pulse_reserve/`.

**Inputs.** NASA POWER daily precipitation, sample metadata.

**Key results.**
- Precip pulses align with shifts in community structure (especially shifts toward Strategy A — DOM-cyclers).
- Strategy A correlates positively with d7 precip (ρ ≈ +0.21).
- Pulse responses are clearest in rhizosphere (where biomass is highest).

**Interpretation.** Pulse-reserve is real but moderate-effect. Frames the climate projection: precipitation pulses drive A→B transitions (Section I).

**Status.** solid.

**Cross-refs.** 33 (alive-only redo), 41 (transition asymmetry), 44 (climate projection v2).
