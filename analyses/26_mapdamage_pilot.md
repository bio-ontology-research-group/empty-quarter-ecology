# 26. mapDamage pilot — chemical-degradation signatures (NEGATIVE)

**Question.** Are there ancient-DNA-style chemical damage signatures (5' C→T, 3' G→A) in the EQ metagenomic reads, which would confirm that the relic fraction is chemically degraded?

**Method.**
- **Track A (pilot)** — Run mapDamage on 10 EQ metagenome samples mapped to their own assemblies. `scripts/parse_mapdamage_pilot.py`, output `cache/mapdamage_pilot/`.
- **Robustness check (MAG remap)** — Realign the same reads against curated MAGs to rule out misalignment artifacts. `scripts/parse_mapdamage_mag.py`, output `cache/mapdamage_mag/`.

**Key results.**
- Pilot: 5' C→T at position 1 = **0.31%** (essentially baseline).
- MAG remap: drops further to **0.16%** with cleaner reference.
- **NO aDNA-style damage signature.**

**Interpretation.** EQ relic DNA is **biologically intact** — i.e. recent dead biomass, not chemically aged aDNA. The relic fraction is best modeled as "currently dead cells" not "old preserved DNA."

This **invalidates aDNA-style chemical relic markers** in the indicator (Track A is dead end). Forces reliance on:
- PMA ground truth (#24)
- Amplicon sequence features (#25)
- Bayesian + MAG-presence priors (#27)

**Status.** solid negative result — methodologically important.

**Outputs.** `cache/mapdamage_pilot/`, `cache/mapdamage_mag/`.

**Cross-refs.** 25, 27.
