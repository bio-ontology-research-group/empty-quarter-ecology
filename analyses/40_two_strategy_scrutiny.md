# 40. Two-strategy scrutiny — compositional + mechanism

**Question.** Is the A-B anti-correlation a compositional artifact? And what are the functional mechanisms of each strategy?

**Method.** `scripts/scrutiny_two_strategy.py`, output `cache/two_strategy_scrutiny/`.

**Compositional artifact check.**

| Method | ρ | p |
|---|---|---|
| Relabund | −0.287 | 1.3e-24 |
| Absolute counts | −0.135 | 2.0e-06 |
| **CLR-mean** | **−0.474** | 1.3e-69 |

Random-pair null (100 random pairs of equivalent sizes): median +0.084, p5..p95 [−0.196, +0.417]. The observed −0.287 sits at the **97th percentile** of the null — moderately extreme, not overwhelming.

**Within-site check (162 site-comp cells, n ≥ 4 samples each).**
- Median ρ = −0.246
- 67% of cells show negative A-B correlation
- 44% show strongly negative (< −0.3)

**Mechanism via metagenomics.**

| Function | Strategy A | Strategy B |
|---|---|---|
| **Sporulation initiation** (spo0A/F/E + spoIIE) | absent | **present** |
| **Ectoine biosynth** (ectA, ectB) | absent | **present** (ectC missing) |
| **Betaine uptake** (opuD) | absent | **present** |
| **Trehalose biosynth** (otsA + otsB) | **present** (Nibribacter, #36) | absent |
| **DNA repair** (specific KOs) | **saturated** (Nibribacter 139% target) | normal |
| **Heat shock** | **saturated** (150%) | normal |

**Interpretation.**
- Strategy A: DOM-cycling, **trehalose-osmoprotected**, DNA-repair-heavy. Relies on persistence + UV/oxidative defenses.
- Strategy B: **sporulation + ectoine biosynthesis + betaine uptake**. Classic halotolerant stress survivor.
- A LACKS osmoprotectant biosynthesis besides trehalose; B LACKS the heavy DNA-repair / heat-shock toolkit.

**Status.** solid. The anti-correlation is statistically robust though moderate (97th pct null). The mechanistic divergence is sharp.

**Major correction included.** The "Nibribacter 96 betaine ORFs" was a regex artifact — actual = 1. See #36.

**Cross-refs.** 36, 38, 39, 41.
