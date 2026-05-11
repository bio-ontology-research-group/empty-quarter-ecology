# 54. Referee defense — 13 sensitivity/robustness analyses

**Question.** What survives an adversarial review pass?

**Date.** 2026-05-08 (immediately before mediation retraction).

**Analyses run (13 total).**
1. Adversarial keystone re-ranking (#53).
2. Network rho-threshold sensitivity (0.3 / 0.4 / 0.5).
3. Module-detection algorithm sensitivity (Louvain / Leiden / spinglass).
4. Compartment-stratified vs pooled network re-fits.
5. Site-level cross-validation of keystone composite score.
6. Bootstrap CI on cross-compartment keystone ranking.
7. Sabkha-score definition robustness.
8. Sample subsampling (down to 75%, 50%).
9. Read-depth floor sensitivity.
10. CSP1-2 ASV-vs-genus aggregation comparison.
11. Thermal-bound: tested whether CSP1-2 has a stricter thermal niche than co-occurring DOM-cyclers.
12. Stoichiometry: tested if XRF supply predicts CSP1-2 dominance.
13. Cross-cohort transferability (#49).

**Outcomes.**
- **#1, #2, #3, #4, #5, #6, #7, #8, #9, #10**: passed — original signals robust.
- **#11 Thermal-bound**: **RETRACTED** — CSP1-2's thermal range overlaps co-occurring DOM-cyclers; not uniquely bounded.
- **#12 Stoichiometry**: **caveated** — supply/demand framing requires micro-spatial co-localisation that we cannot resolve. See #23.
- **#13 Cross-cohort**: **failed transfer** — Hill profile and CSP1-2 niche don't generalise to Atacama. See #49.

**Status.** Referee-defense bundle. Concrete retractions captured in #19, #23, #49.

**Cross-refs.** 19, 22, 23, 38, 49, 53.
