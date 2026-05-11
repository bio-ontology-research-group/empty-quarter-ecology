#!/usr/bin/env python3
"""Hunt for the actual keystone genus / genera in the MAG-augmented alive
co-occurrence networks.

For each compartment:
  1. Build genus-level alive-only network
  2. Compute degree, betweenness, closeness, eigenvector centrality
  3. Compute composite keystone score (Banerjee et al. style)
  4. Rank top-30 candidates per compartment
  5. Cross-reference: which appear in >=2 / 3 compartments at top-N?
  6. Layer in MAG presence + abundance + pulse-response signal

Final ranking criteria for "true keystone":
  - High composite score in >=2 compartments
  - MAG presence (at least one ASV in genus matches a MAG)
  - Reasonable abundance (not pure rarity artifact)
  - Bonus: pulse-responder OR known biocrust/N-fixer/water-stress role

Outputs:
  cache/keystone_hunt/per_compartment_keystones.tsv
  cache/keystone_hunt/cross_compartment_ranking.tsv
  cache/keystone_hunt/edges_<compartment>.tsv
  cache/keystone_hunt/summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "src"))
from _sample_parse import parse_samples_to_df
from eq.network import (
    compositional_correlation, build_network, louvain_modules, keystone_score
)

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "keystone_hunt"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    print("Loading inputs ...", flush=True)
    p = pd.read_csv(CACHE / "relic_priors" /
                     "relic_score_with_mag_prior.tsv", sep="\t")
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)
    tax = pd.read_parquet(CACHE / "taxonomy.parquet").reset_index().rename(
        columns={"ASV": "asv_id"})

    # Alive pool (MAG-augmented)
    alive_mag = set(p.loc[p["relic_score_with_mag"] <= 0.3, "asv_id"])
    print(f"  alive (MAG-aug) ASVs: {len(alive_mag)}", flush=True)
    ft_a = ft.loc[ft.index.isin(alive_mag)]

    # Aggregate to genus
    ft2 = ft_a.copy()
    ft2.index = ft2.index.rename("asv_id")
    m = ft2.reset_index().merge(tax[["asv_id", "phylum", "class", "genus"]],
                                  on="asv_id", how="left")
    m = m.dropna(subset=["genus"])
    m = m[m["genus"].astype(str).ne("")]
    m = m[m["genus"].astype(str).ne("NA")]
    sample_cols = [c for c in m.columns if c not in
                       ("asv_id", "phylum", "class", "genus")]
    gen = m.groupby("genus")[sample_cols].sum()
    print(f"  genera: {gen.shape[0]}", flush=True)

    # Per-compartment networks
    all_keystones = []
    for comp in ("rhizosphere", "surface", "deep"):
        samps = list(set(smeta[smeta["compartment"] == comp]["sample"]) &
                       set(gen.columns))
        if len(samps) < 30: continue
        sub = gen[samps]
        # Network
        try:
            rho, pval = compositional_correlation(sub, min_prevalence=0.10,
                                                       presence_ra=0.0001,
                                                       pseudo=0.5)
            G = build_network(rho, pval, rho_threshold=0.4, q_threshold=0.05)
            mods = louvain_modules(G, seed=42)
        except Exception as e:
            print(f"  [{comp}] failed: {e}", flush=True); continue
        if G.number_of_nodes() == 0: continue

        ks = keystone_score(G)
        ks["module"] = ks["node"].map(lambda n: mods.get(n, -1))
        # Add eigenvector centrality
        try:
            ev = nx.eigenvector_centrality_numpy(G, weight="weight")
            ks["eigenvector"] = ks["node"].map(ev)
        except Exception:
            ks["eigenvector"] = np.nan
        # Connector role: ASVs whose neighbors span multiple modules
        connector_score = {}
        for n in G.nodes:
            ne_modules = [mods.get(nn, -1) for nn in G.neighbors(n)]
            uniq = len(set(ne_modules))
            connector_score[n] = uniq
        ks["n_neighbor_modules"] = ks["node"].map(connector_score)

        # Mean abundance per genus in this compartment
        gen_relabund = sub.div(sub.sum(axis=0).replace(0, 1), axis=1)
        mean_relabund = gen_relabund.mean(axis=1)
        ks["mean_relabund"] = ks["node"].map(mean_relabund)

        ks["compartment"] = comp
        ks.to_csv(OUT / f"per_genus_keystone_{comp}.tsv",
                   sep="\t", index=False)

        # Save edges
        if G.number_of_edges() > 0:
            ed = pd.DataFrame([(u, v, d.get("rho"), d.get("weight"),
                                  d.get("sign"), d.get("q"))
                                  for u, v, d in G.edges(data=True)],
                                 columns=["source", "target", "rho", "weight",
                                            "sign", "q"])
            ed.to_csv(OUT / f"edges_{comp}.tsv", sep="\t", index=False)

        # Print top-15 by keystone score
        print(f"\n--- {comp.upper()} (n_genera={G.number_of_nodes()}, "
              f"n_edges={G.number_of_edges()}) ---", flush=True)
        print(f"Top 15 by composite keystone score:")
        cols = ["node", "degree", "betweenness", "closeness", "keystone",
                  "eigenvector", "n_neighbor_modules", "mean_relabund",
                  "module"]
        print(ks.head(15)[cols].round(4).to_string(index=False), flush=True)

        all_keystones.append(ks)

    # ============================================================
    # Cross-compartment ranking
    # ============================================================
    if not all_keystones:
        print("No networks built"); return
    K = pd.concat(all_keystones, ignore_index=True)
    K.to_csv(OUT / "per_compartment_keystones.tsv", sep="\t", index=False)

    # Per-genus aggregate across compartments
    cross = (K.groupby("node")
             .agg(n_compartments=("compartment", "nunique"),
                  mean_keystone=("keystone", "mean"),
                  max_keystone=("keystone", "max"),
                  mean_degree=("degree", "mean"),
                  mean_betweenness=("betweenness", "mean"),
                  mean_eigenvector=("eigenvector", "mean"),
                  mean_relabund=("mean_relabund", "mean"),
                  mean_n_neighbor_modules=("n_neighbor_modules", "mean"),
                  )
             .reset_index().rename(columns={"node": "genus"}))

    # Add MAG presence at genus level
    asvs_with_mag = set(p.loc[p["has_match"] == 1, "asv_id"])
    asv_to_genus = tax.set_index("asv_id")["genus"].to_dict()
    asv_to_phylum = tax.set_index("asv_id")["phylum"].to_dict()
    genera_with_mag = set()
    for asv in asvs_with_mag:
        g = asv_to_genus.get(asv)
        if isinstance(g, str) and g and g != "NA":
            genera_with_mag.add(g)
    cross["has_mag_evidence"] = cross["genus"].isin(genera_with_mag).astype(int)

    # Add phylum
    genus_to_phylum = (tax.dropna(subset=["genus", "phylum"])
                         .drop_duplicates("genus")
                         .set_index("genus")["phylum"].to_dict())
    cross["phylum"] = cross["genus"].map(genus_to_phylum)

    # Composite cross-compartment score: mean_keystone * n_compartments / 3
    # weighted by MAG evidence
    cross["composite_cross"] = (cross["mean_keystone"] *
                                    (cross["n_compartments"] / 3.0))
    cross["composite_validated"] = (cross["composite_cross"] *
                                         (1 + 0.5 * cross["has_mag_evidence"]))

    # Sort
    cross = cross.sort_values("composite_validated", ascending=False)
    cross.to_csv(OUT / "cross_compartment_ranking.tsv",
                  sep="\t", index=False)

    print("\n========================================================", flush=True)
    print("CROSS-COMPARTMENT TOP 30 KEYSTONE CANDIDATES", flush=True)
    print("========================================================", flush=True)
    print(cross.head(30)[["genus", "phylum", "n_compartments",
                              "mean_keystone", "mean_degree",
                              "mean_betweenness", "mean_eigenvector",
                              "mean_relabund", "has_mag_evidence",
                              "composite_validated"]].round(4).to_string(
        index=False), flush=True)

    # CSP1-2 specifically
    csp = cross[cross["genus"].str.contains("CSP1", na=False)]
    print("\n--- CSP1-2 cross-compartment score (for comparison) ---")
    print(csp.round(4).to_string(index=False), flush=True)

    # Write summary
    with open(OUT / "summary.txt", "w") as fh:
        fh.write("Keystone hunt in MAG-augmented alive networks\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Alive ASVs (MAG-augmented): {len(alive_mag)}\n")
        fh.write(f"Genera tested: {gen.shape[0]}\n\n")

        fh.write("Top 30 by cross-compartment composite_validated score "
                  "(adjusts for MAG presence):\n\n")
        fh.write(cross.head(30)[["genus", "phylum", "n_compartments",
                                      "mean_keystone", "mean_degree",
                                      "mean_betweenness", "mean_eigenvector",
                                      "mean_relabund", "has_mag_evidence",
                                      "composite_validated"]].round(4)
                  .to_string(index=False))

        fh.write("\n\nCSP1-2 reference:\n")
        fh.write(csp.round(4).to_string(index=False))

        # Cross-compartment hubs (in >=2 compartments)
        crossh = cross[cross["n_compartments"] >= 2]
        fh.write(f"\n\nGenera in >=2 compartments: {len(crossh)}\n")
        fh.write(f"Genera in all 3 compartments: "
                  f"{(cross['n_compartments'] == 3).sum()}\n")
    print(f"\nSummary -> {OUT}/summary.txt", flush=True)


if __name__ == "__main__":
    main()
