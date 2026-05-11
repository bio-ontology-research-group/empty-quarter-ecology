# 17. Lagged compartment coupling

**Question.** Do compartments couple temporally — does last-trip surface predict this-trip deep (or vice versa)?

**Method.** Per (site, ASV), regress current compartment abundance on prior trip's other compartment. `cache/lagged_compartment_coupling.tsv`, `cache/paired_surface_deep.tsv`.

**Key results.**
- Significant surface→deep coupling at one-trip lag; deep→surface much weaker.
- Compartments are not independent — supports a top-down (surface-as-source) framing.

**Status.** solid.

**Cross-refs.** 6 (iCAMP), 33.
