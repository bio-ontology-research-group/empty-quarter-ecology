# 09. EMP cosmopolitanism (Test 2)

**Question.** Are EQ taxa endemic, or are they globally distributed (consistent with high dispersal)?

**Method.** Match EQ ASVs against Earth Microbiome Project (EMP) global database; classify each as cosmopolitan (≥X global samples) vs endemic. `scripts/test2_emp_cosmopolitanism.py`.

**Inputs.**
- EQ ASVs (`cache/feature_table.parquet`)
- EMP reference (`cache/emp_cosmopolitanism/`)

**Key results.**
- **>50% of EQ taxa** are EMP-cosmopolitan (present in multiple global biomes).
- Endemic fraction is concentrated in low-abundance specialists.

**Interpretation.** EQ doesn't have a unique microbiome in a strict sense — most of its members are everywhere. The "uniqueness" must therefore be in **abundance structure and functional dominance**, not species identity. Reinforces the dispersal-driven framing (#6, #7).

**Status.** solid.

**Cross-refs.** 6, 7, 11 (cross-trip ephemerality).
