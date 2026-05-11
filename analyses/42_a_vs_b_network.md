# 42. A-only vs B-only stratified networks

**Question.** Is the network architecture different at A-dominant vs B-dominant samples?

**Method.** Build alive-only co-occurrence networks **separately** on A-dominant samples (n=1,029) and B-dominant samples (n=207). Genus-level, CLR + Spearman + BH q < 0.05, |ρ| ≥ 0.4. `scripts/network_A_vs_B_sites.py`, `cache/network_A_vs_B/`.

**Key results.**

| Class       | n_samples | n_genera | n_edges | n_modules |
|-------------|-----------|----------|---------|-----------|
| A-dominant  | 1,029     | 61       | 87      | 13        |
| B-dominant  | 207       | 60       | **208** | 5         |

**The B-dominant network is 2.4× denser with fewer modules.**

**Top keystones.**
- A-dominant: Tumebacillus, Neobacillus, Rubrobacter, Oceanobacillus, Anseongella, **Nibribacter (rank 6)**.
- B-dominant: **Aquibacillus (rank 1)**, then mostly other Bacilli/halotolerants. Nibribacter rank ~35.

**Top-20 keystone overlap: 9/20 common** (mostly Bacilli core). A-only keystones include Nibribacter (DOM-cycling); B-only keystones include Halomonadaceae + additional Bacilli.

**Interpretation.**
- B-dominant network is more **integrated** (many co-fluctuating taxa). Consistent with B being a **stress/pulse state** where many taxa share a co-ordinated stress response (sporulation + ectoine).
- A-dominant network is more **partitioned** into modules. Consistent with A being a **steady-state community** with guilds occupying distinct functional niches.

**Important nuance (added 2026-05-11).** Nibribacter rank 6 ≠ "central reporter" — 3 of top 5 in A-dominant are still Bacilli (Tumebacillus, Neobacillus, Oceanobacillus) at tiny mean relabund (≤ 0.06%). The proper claim is: Nibribacter is the **highest-ranked Bacteroidota DOM-cycler in the top tier**, and the most-abundant top-tier member (~4.2%).

**Status.** solid.

**Outputs.**
- `cache/network_A_vs_B/keystone_{A,B}_dominant.tsv`
- `cache/network_A_vs_B/edges_{A,B}_dominant.tsv`
- `cache/network_A_vs_B/network_summary.tsv`

**Cross-refs.** 34, 38, 39, 40.
