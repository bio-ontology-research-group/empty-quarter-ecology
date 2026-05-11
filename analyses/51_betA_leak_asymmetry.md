# 51. betA leak-asymmetry — FAILS (key negative)

**Question.** Are the betA-producer bins purely producers (Black Queen donors), or do they also encode betaine uptake transporters (compromising the donor/recipient asymmetry)?

**Method.** For each of the 152 producer bins (#50), test for presence of betaine uptake transporters (opuD, proU, proV/W/X). `scripts/run_leak_asymmetry_test.py`, `cache/leak_asymmetry_per_bin.tsv`, `cache/leak_asymmetry_test.txt`.

**Key results.**
- **54.8% of producer bins ALSO carry uptake transporters.**
- The clean donor/recipient asymmetry assumed by Black Queen models **fails** at the genome level.
- Producers are not pure donors; they are dual-mode (produce + scavenge).

**Implication for narrative.**
- The original CSP1-2 "Black Queen keystone donor" story does not survive the leak-asymmetry test.
- "Guild indicates stress" (presence of producers correlates with arid conditions) is preserved, but "guild architects diversity by leaking osmolytes" is not.
- Paper must be reframed from causal-architect language to indicator-of-stress language.

**Status.** solid **negative result**. Paper-rewriting finding (alongside #38).

**Outputs.**
- `cache/leak_asymmetry_per_bin.tsv`
- `cache/leak_asymmetry_test.txt`

**Cross-refs.** 13, 50, 52.
