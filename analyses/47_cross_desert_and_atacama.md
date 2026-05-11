# 47. Cross-desert comparisons — Gurbantunggut, Namib, McMurdo, Atacama

**Question.** How does EQ compare to other hyperarid systems? Are its dominant taxa and community structure shared or unique?

**Method.** Re-process publicly available 16S datasets from:
- **Gurbantunggut** (Mongolia/China) — `cache/crossdesert/gurbantunggut_*` + reprocessed manifest
- **Namib** — `cache/crossdesert/namib_runs.tsv`
- **McMurdo dry valleys** (Antarctica) — `cache/crossdesert/mcmurdo_runs.tsv`
- **Atacama** (Chile) — separate workflow (#48)

Re-run through same DADA2 + SILVA pipeline; compare genus-level relative abundance, Shannon, β-diversity to EQ. `cache/crossdesert/`, scripts `stage1_download.sh` → `stage5_genus_compare.py`.

**Key results — comparison summary** (`cache/crossdesert/comparison_summary.tsv`).
- Gurbantunggut shares ~50% of EQ-dominant genera but with different rank ordering.
- Namib salt-flats overlap heavily on Halomonas / halotolerant Bacilli.
- McMurdo Dry Valleys diverge: psychrophile-dominated, much lower α-diversity.
- The "halotolerant core" (Strategy B in EQ) is **not unique** — present in all salt-influenced deserts.
- The "DOM-cycler" backbone (Strategy A — Bacteroidota + Massilia) is **EQ-distinctive** at the dominance level.

**Salinity calibration.** `cache/cross_desert_salinity_calibration.tsv` — calibrates EQ XRF salinity proxies against measured salinity in other systems.

**Status.** solid. Bookends with #48 (Atacama detail), #19 (thermal calibration), #54 (referee defense).

**Outputs.**
- `cache/crossdesert/comparison_summary.tsv`
- `cache/crossdesert/desert_vs_eq_genus_spearman.tsv`
- `cache/crossdesert/genus_median_rel_abund.tsv`
- `cache/crossdesert/eq_shannon_reference.tsv`

**Cross-refs.** 19, 48, 49, 54.
