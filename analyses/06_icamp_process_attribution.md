# 06. iCAMP — process attribution

**Question.** What fraction of EQ community assembly is deterministic (selection) vs stochastic (dispersal, drift)?

**Method.** iCAMP framework (Ning et al. 2020): bin-level βNTI + RCbray to partition pairwise comparisons among five processes (homogeneous selection, heterogeneous selection, dispersal limitation, homogenizing dispersal, drift). `scripts/run_icamp_rcbray.py`, `scripts/run_bnti.py`. Run per compartment.

**Inputs.**
- `cache/feature_table.parquet`, `cache/trees/*` for phylogeny

**Key results (process_summary_all.tsv).**

| Compartment   | Homogenizing dispersal | Selection (homog + hetero) | Drift |
|---------------|------------------------|----------------------------|-------|
| Surface       | ~67%                   | ~7%                        | ~26%  |
| Deep          | ~67%                   | ~15%                       | ~18%  |
| Rhizosphere   | ~67%                   | ~15%                       | ~18%  |

**Interpretation.** EQ is a **dispersal-driven** system, not a selection-driven one. The relatively uniform 67% homogenizing dispersal across all three compartments is unusual and points to a mechanism that mixes ASVs system-wide regardless of micro-niche (entry-point for wind hypothesis #7).

**Status.** solid. Major candidate reframe of the paper.

**Outputs.**
- `cache/icamp/process_summary_all.{tsv,txt}`
- `cache/icamp/process_attribution_{surface,deep,rhizosphere}.tsv`
- `cache/icamp/RCbray_{surface,deep,rhizosphere}.parquet`

**Cross-refs.** 07, 08, 14 (functional iCAMP near-uniform → an even more striking version of the same finding), 33.
