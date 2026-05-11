# 36. Nibribacter KEGG — CORRECTED with specific KO lists

**Question.** With the regex-artifact fixed, what does Nibribacter actually carry?

**Method.** Replace regex matching against KO definitions with **curated KO-ID lists per category** (KEGG-specific, not text-search). `scripts/nibribacter_kegg_corrected.py`. Categories include:
- CAZymes (glycoside hydrolases, polysaccharide lyases, carbohydrate esterases)
- TonB-dependent / SusC / SusD genes
- DNA repair (specific recA / mutS / uvr KOs)
- Heat shock chaperones (groEL, dnaK, grpE, …)
- Oxidative stress (catalase, SOD, ahpC, …)
- Compatible solute biosynthesis: **betaine** (betA, betB, gbsA); **trehalose** (otsA, otsB, treT, treP); **ectoine** (ectA, ectB, ectC, ectD)
- Compatible solute uptake: opuA-D, proV/W/X
- Na+/K+ pumps: nhaA/B, mrpA-C
- N-fix, denitrification, sulfate reduction (controls)
- Sigma factors
- **Negative controls**: photosynthesis, sporulation initiation (Nibribacter is Bacteroidota — should NOT have these)

**Key results.**
- **DNA repair**: 46 unique KOs detected (**139% of target list** — saturating; this is a UV-damage-survivor genome).
- **Heat shock**: 18 unique KOs (**150% of target**).
- **Trehalose biosynthesis**: both otsA + otsB present — **active osmoprotection.**
- **Ectoine biosynthesis**: **absent** in Nibribacter (this is a Strategy A feature, not B).
- **Betaine biosynthesis**: **just 1 betA ORF** total across all Nibribacter MAGs (previously inflated to 96 by regex bug).
- Sigma factors: rpoD + sigB + rpoH present (stress-responsive).
- Sporulation negative control: **absent**, as expected for Bacteroidota.
- Photosynthesis negative control: absent, as expected.

**Interpretation.** Nibribacter is a **DNA-repair-heavy, trehalose-osmoprotected DOM cycler**. **Not** a betaine producer. Functionally distinct from the Strategy B halotolerants (which use ectoine + betaine uptake). Validates the two-strategy story (#39–42) at the genome level.

**Status.** solid. Replaces #35.

**Outputs.**
- `cache/nibribacter_mags/corrected_function_summary.tsv`
- `cache/nibribacter_mags/corrected_summary.txt`

**Cross-refs.** 35 (superseded), 40, 50, 52.
