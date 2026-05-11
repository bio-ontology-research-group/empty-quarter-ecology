# 52. betA x-comparator — multi-genus betA carriage

**Question.** Is betA restricted to CSP1-2 / Nibribacter (single-keystone framing), or does it cross multiple genera?

**Method.** Per-assembly + per-genome census of K00108-bearing bins, with phylogenetic placement of the betA gene. `scripts/run_xcomparator_betA_summary.py`, `scripts/render_xcomparator_betA.py`. Outputs in `cache/xcomparator_betA_*.tsv` and `cache/xcomparator_betA_*.faa`.

**Key results.**
- **Rubellimicrobium** also carries functional betA (alongside Nibribacter, Flavisolibacter, Solirubrobacter, Telluribacter, …).
- Several genera, not one — i.e. betA carriage is a **guild trait**, not a single-genus signature.
- Phylogeny of betA: distinct clades per genus → multiple independent origins, not horizontal transfer from a single source.

**Implication for narrative.**
- **Single-keystone framing retracted** (alongside #38).
- New framing: betA carriage is a **community-level guild trait** that EQ communities deploy under stress, with multiple genera contributing.

**Status.** solid.

**Outputs.**
- `cache/xcomparator_betA_summary.tsv`
- `cache/xcomparator_betA_per_genome.tsv`
- `cache/xcomparator_betA_per_assembly.tsv`
- `cache/xcomparator_betA_aln.faa`
- `cache/xcomparator_betA_top_seqs.faa`

**Cross-refs.** 38, 50, 51.
