# 38. Single-keystone vs three-guild — discriminating knockout test

**Question.** Is the alive-network architecture better described by ONE keystone (Nibribacter) or by THREE functional guilds (DOM-cyclers, halotolerants, low-abundance specialists)?

**Method.** Two ablation experiments on the alive co-occurrence network:
1. Single-keystone knockout — set Nibribacter alone to zero; recompute connectivity & module structure.
2. Guild knockouts — set each of three a priori guilds to zero in turn:
   - **Bact_DOM** (Bacteroidota DOM-cyclers): Nibribacter, Flavisolibacter, Solirubrobacter, Telluribacter, Rubellimicrobium, …
   - **Bacilli** (halotolerant Strategy B): Aquibacillus, Oceanobacillus, Halobacillus, Halomonas, …
   - **Low-abundance specialists**: bottom-decile taxa.

Compare connectivity loss. `scripts/keystone_vs_guild_test.py`, `cache/keystone_test/`.

**Key results.**

| Knockout | Connectivity loss |
|---|---|
| Nibribacter alone | −7% to −16% |
| Bact_DOM guild | **−29% to −37%** |
| Bacilli guild | **−46% to −57%** |

**Interpretation.** The system is **guild-organized**, not single-keystone-organized. The Bacilli guild (Strategy B) carries the **largest** connectivity weight despite Nibribacter's individual centrality. Decisive evidence for the three-guild (eventually two-strategy) framing.

**Status.** solid. Paper-rewriting finding — replaces the original single-keystone narrative.

**Cross-refs.** 22 (early JSDM hint), 32, 34, 39, 53.
