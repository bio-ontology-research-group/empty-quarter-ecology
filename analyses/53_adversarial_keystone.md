# 53. Adversarial keystone (referee defense)

**Question.** If a hostile referee adversarially picks the "best alternative keystone," how does the original CSP1-2 claim stack up?

**Method.** Compute knockout magnitude (loss of connectivity, modularity change) for every candidate genus in the alive co-occurrence network. Rank. `scripts/run_adversarial_keystone.py`, `cache/jsdm_knockout_*` and `cache/keystone_knockout.tsv`.

**Key results.**
- CSP1-2 ranks **5th out of 15** candidate "keystones" by knockout magnitude.
- The 4 above it are: Aquibacillus, Halomonas, Tumebacillus, and Nibribacter.
- CSP1-2's uniqueness is **mechanistic** (it encodes a stress-response gene set) — **not topological** (it's not THE most connected node).

**Interpretation.** Confirms #38 — single-keystone framing is untenable. Paper must talk about CSP1-2 (or whichever taxon) as a **mechanistic indicator**, not a topological hub.

**Status.** solid. Defense-ready.

**Cross-refs.** 22, 32, 38, 54.
