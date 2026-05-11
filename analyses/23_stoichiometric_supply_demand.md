# 23. Stoichiometric supply / demand

**Question.** Is there a stoichiometric mismatch between element supply (XRF-measured soil) and element demand (taxonomic biomass requirements)?

**Method.** Build supply (XRF) and demand (genome-imputed) per element per site, ratio them. `scripts/run_stoichiometric_model.py`, `cache/stoichiometric_supply_demand.tsv`.

**Key results.**
- Site-level supply/demand mismatch is heterogeneous: some sites are P-limited, others N-limited.
- Compartment differences in demand profile align with rhizosphere plant-host effects.

**Caveat.** Original paper hypothesis "stoichiometric supply explains keystone dominance" requires **micro-spatial co-localisation** of supply and consumer, which we cannot resolve at our spatial scale (referee-defense memo #54). Claim was downgraded.

**Status.** caveated.

**Cross-refs.** 2, 54.
