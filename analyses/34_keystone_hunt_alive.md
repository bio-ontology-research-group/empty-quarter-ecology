# 34. Keystone hunt — alive network

**Question.** If CSP1-2 is not THE keystone, who is — on the alive (MAG-prior-corrected) network?

**Method.** Build genus-level alive co-occurrence networks (CLR + Spearman + BH q < 0.05, |ρ| ≥ 0.4) per compartment. Score each genus on degree × betweenness × closeness composite. `scripts/relic_keystone_hunt.py`, output `cache/keystone_hunt/`.

**Key results.**
- Per-compartment top keystones differ; **Nibribacter** is the only **cross-compartment** top hub.
- Cross-compartment composite "validated" keystone score: **Nibribacter = 1.21** vs **CSP1-2 = 0.21** (∼6× larger).
- Nibribacter is a **Bacteroidota DOM-cycler** — the natural alive-network counterpart of what CSP1-2 was claimed to be.

**Caveat (added 2026-05-11).** Nibribacter ranks **6th** in the A-dominant network — behind Tumebacillus, Neobacillus, Rubrobacter, Oceanobacillus, Anseongella (3 of 5 still Bacilli). What makes Nibribacter notable is:
- It's the **highest-ranked Bacteroidota DOM-cycler** in the top tier.
- It's the **most abundant** node in the top tier (~4.2% mean relabund vs ≤0.06% for the higher-ranked Bacilli).

**Status.** solid (with the rank-6 caveat).

**Outputs.**
- `cache/keystone_hunt/cross_compartment_ranking.tsv`
- `cache/keystone_hunt/per_compartment_keystones.tsv`
- `cache/keystone_hunt/per_genus_keystone_*.tsv`
- `cache/keystone_hunt/edges_*.tsv`

**Cross-refs.** 32, 35, 38, 42.
