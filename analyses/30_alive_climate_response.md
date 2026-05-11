# 30. Alive-only climate-Shannon response

**Question.** When we re-do "Shannon diversity vs MAT" on the alive community only, does the gradient hold?

**Method.** Per-(sample) Shannon on alive subset, regressed against MAT, MAP, sabkha proxies, lat/lon. `scripts/relic_alive_climate_response.py`.

**Key results.**
- **All ASVs**: Shannon ~ MAT, ρ = **−0.40** (strong negative — colder sites more diverse).
- **Alive only**: Shannon ~ MAT, ρ ≈ **+0.02** (effectively zero).
- **Relic only**: ρ ≈ −0.45 (the relic fraction carries the gradient).

**Interpretation.** **The climate–diversity gradient is a relic-DNA artifact.** Wind-deposited dead biomass accumulates differently across climate zones, but the living community is climate-insensitive at this resolution. Massive narrative implication — the "diversity decreases with aridity" framing in the original paper is **driven by dead cells**.

**Status.** solid. Major paper-changing finding.

**Cross-refs.** 7 (wind dispersal explains why relic carries the gradient), 33, 45.
