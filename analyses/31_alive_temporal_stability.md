# 31. Alive-only all-pairwise temporal stability

**Question.** Are alive communities more or less temporally stable than the full (relic-contaminated) community?

**Method.** Per (site, comp), compute BC distance for **all pairwise** trip combinations (not just T5 vs mean-T1-4 as originally) on alive, all, and relic feature tables separately. `scripts/relic_alive_remaining_analyses.py`.

**Key results.**
- Alive community shows **~2× tighter** BC distances between non-adjacent trips than the full community.
- Relic-only BC distances are larger and noisier — the relic pool churns more.
- The seasonal/trip signal is sharper in alive than in all.

**Interpretation.** The alive backbone is **temporally consistent**; the apparent turnover in the full data is largely relic churn driven by wind deposition (#7).

**Status.** solid. Reframes earlier "high turnover" claims.

**Cross-refs.** 11 (cross-trip persistence; alive halves the ephemerality rate), 17, 30.
