# 14. Functional iCAMP (Test 1)

**Question.** Same iCAMP framework as #6 — but applied to PICRUSt2 functional pathway distances instead of phylogenetic distances. Is the system functionally as well as taxonomically assembled by dispersal?

**Method.** βNTI/RCbray on pathway-level distances. `scripts/test1_functional_icamp.py`, output `cache/test1_functional_icamp/`.

**Key results.**
- **Functional iCAMP near-uniform (98%+)** dominated by stochastic processes — much more uniform than taxonomic (67%).
- Selection signal at the functional level is **vanishingly small**.

**Interpretation.** Even more striking than the taxonomic finding: at the function level, EQ looks like **noise**. Same potential pathways are everywhere; environment selects negligibly. Fits the high-redundancy framing (#12).

**Status.** solid — but treat as PICRUSt2-prediction-derived, not metagenome-confirmed.

**Cross-refs.** 6, 12, 13.
