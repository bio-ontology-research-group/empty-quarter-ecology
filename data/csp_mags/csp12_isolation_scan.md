# CSP1-2 isolation / enrichment scan — /data/emptyquarter/sequencing-results

**Purpose.**  Task #6 of the pre-submission strengthening pass: check whether any
assembly produced in parallel isolation, enrichment, or culture attempts on
Empty-Quarter samples contains a CSP1-2 genome.

**Date:** 2026-04-19

**Data scanned.**  All `.fa` / `.fasta` / `.fna` files recursively under
`/data/emptyquarter/sequencing-results/` in the following sub-trees
(excluding `sites/`):

| Sub-tree | # Assembly files found | Notes |
|----------|-------------------------|-------|
| `cultures/`         | 51  | MEGAHIT single-sample bins from culture dishes |
| `enrichment/`       | 341 | metaflye long-read assemblies + binned MAGs |
| `isolates/`         | 122 | bacterial isolates, mostly from sites 59 & 60 |
| `rh_sequencing/`    | 10  | rhizosphere targeted sequencing |
| `thermophilic/`     | 0   | directory contains only `culture_conditions.txt` + QC tsv files; no assemblies deposited |
| `other/`            | 0   | contains controls and fastqs only, no assemblies |
| **Total**           | **524** | |

**Method.**

1. Recursively collected all `.fa`/`.fasta`/`.fna` files in the six sub-trees.
2. Barrnap v0.9 extracted 12,995 16S rRNA genes from the 524 assemblies.
3. VSEARCH --usearch_global against the four Empty-Quarter CSP1-2 MAG-associated
   ASVs (`csp_ref.fasta`) at identities 97, 90, 85, 75%.
4. skani v0.3.1 whole-genome ANI between all 524 assemblies and the four
   reference CSP1-2 MAGs
   (V27Dr2\_\_SemiBin\_73, V30PRr1\_\_SemiBin\_7, V32PRr1\_\_SemiBin\_26,
   V38PRr3\_\_SemiBin\_38), min-AF 5%, min-ANI 65%.

**Results.**

| Test | Hits |
|------|------|
| 16S  ≥ 97% V4 identity | 0 |
| 16S  ≥ 90% V4 identity | 0 |
| 16S  ≥ 85% V4 identity | 0 |
| 16S  ≥ 75% V4 identity | 3{,}953 (non-specific, includes any conserved-region match across most bacteria) |
| skani whole-genome ≥ 65% ANI | 0 |

**Conclusion.**

**No CSP1-2 has been isolated, enriched, or captured in any of the 524
parallel-project assemblies.**  The lineage remains uncultivated in the
Empty-Quarter workflow as well as in the published literature.  The
radiation-tolerance, osmo-protection, and alternative-nitrogenase
phenotype predicted from the four MAGs is therefore not currently
testable on a bench culture.

This is consistent with the literature status of *Dadabacteria* as
fully uncultured.  The result reinforces the priority of targeted
enrichment (osmoprotectant + vanadium-supplied medium with desiccation
cycling) and SIP / metatranscriptomics on fresh soil, as recommended in
the Discussion.

**Reproducibility.**

Raw outputs at `unimatrix01:/data/emptyquarter/ecology-paper-runs/csp_search/full_2257/`:
- `all_assemblies.txt` — list of scanned assembly paths
- `all_16S.fasta` — extracted 16S rRNA genes
- `hits16S_id{0.97,0.90,0.85,0.75}.b6` — VSEARCH hits per identity threshold
- `skani_vs_all.tsv` — whole-genome ANI output

Scan script: `unimatrix01:/data/emptyquarter/ecology-paper-runs/csp_scan_full.sbatch`
Slurm job ID: 2257 (2026-04-19 21:47--22:05 UTC+3).
