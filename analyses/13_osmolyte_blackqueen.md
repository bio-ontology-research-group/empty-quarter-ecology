# 13. Osmolyte uptake vs biosynthesis — Black Queen (Test 5)

**Question.** Does the EQ community look like a Black Queen system, where most members rely on **uptake** of compatible solutes leaked by a small fraction of producers?

**Method.** Aggregate PICRUSt2 KOs into biosynthesis (e.g. ectABC, otsAB, betAB) vs uptake (opu*, proU*) groups; compute ratio across the community. `scripts/test5_osmolyte_blackqueen.py`, output `cache/test5_osmolyte/`.

**Key results.**
- Uptake / biosynthesis **ratio ≈ 230×**.
- Massive bias toward uptake — classic Black Queen signature.

**Interpretation (initial).** A small minority produce compatible solutes; the rest take them up. The original CSP1-2 keystone hypothesis was framed around this Black-Queen "producer" role.

**Critical update (post-metagenomics).** Test 5 used PICRUSt2 (16S → predicted KOs), which inflates redundancy artifacts. Metagenomic guild census (#50) **fails leak-asymmetry** (#51): 54.8% of producer bins **also carry uptake genes**. The neat producer/consumer split doesn't survive metagenomic scrutiny.

**Status.** caveated — Test 5 result stands as amplicon-PICRUSt2 statistic but **mechanistic interpretation overturned** by metagenomic data.

**Cross-refs.** 50, 51, 13 supersedes itself via 51.
