#!/usr/bin/env python3
"""Build separate co-occurrence networks at A-dominant vs B-dominant samples.

Tests whether:
  - Strategy A and B form a single integrated network (Nibribacter still
    central in B-dominant samples) OR
  - Two distinct networks exist (Nibribacter only central at A-dominant sites)

For each subset:
  - genus-level network (CLR + Spearman + BH)
  - keystone scores per genus
  - module structure
  - Nibribacter centrality

Outputs in cache/network_A_vs_B/.
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
OUT = CACHE / "network_A_vs_B"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    print("Loading inputs ...", flush=True)
    pst = pd.read_csv(CACHE / "two_strategy_temporal" /
                       "per_sample_strategy_with_precip.tsv", sep="\t")
    pst["dominant"] = np.where(pst["log2_A_over_B"] > 0, "A", "B")

    p = pd.read_csv(CACHE / "relic_priors" /
                     "relic_score_with_mag_prior.tsv", sep="\t")
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)
    tax = pd.read_parquet(CACHE / "taxonomy.parquet").reset_index().rename(
        columns={"ASV": "asv_id"})

    alive_mag = set(p.loc[p["relic_score_with_mag"] <= 0.3, "asv_id"])
    ft_a = ft.loc[ft.index.isin(alive_mag)]

    # Genus aggregation
    ft2 = ft_a.copy()
    ft2.index = ft2.index.rename("asv_id")
    m = ft2.reset_index().merge(tax[["asv_id", "phylum", "genus"]],
                                  on="asv_id", how="left")
    m = m.dropna(subset=["genus"])
    m = m[~m["genus"].astype(str).isin(["NA", ""])]
    sample_cols = [c for c in m.columns if c not in
                       ("asv_id", "phylum", "genus")]
    gen = m.groupby("genus")[sample_cols].sum()

    # Map sample -> dominant strategy
    s2d = pst.set_index("sample")["dominant"].to_dict()

    # Per dominance class, build network across all compartments
    rec_summary = []
    for cls in ("A", "B"):
        samps = [s for s in gen.columns if s2d.get(s) == cls]
        print(f"\n=== Building network for {cls}-dominant samples (n={len(samps)}) ===",
              flush=True)
        sub = gen[samps]
        # Filter genera with min prevalence >= 10%
        try:
            rho, pval = compositional_correlation(sub, min_prevalence=0.10,
                                                       presence_ra=0.0001,
                                                       pseudo=0.5)
            G = build_network(rho, pval, rho_threshold=0.4,
                                q_threshold=0.05)
            mods = louvain_modules(G, seed=42)
        except Exception as e:
            print(f"  failed: {e}"); continue

        n_nodes = G.number_of_nodes()
        n_edges = G.number_of_edges()
        print(f"  Nodes: {n_nodes}  Edges: {n_edges}", flush=True)

        if n_nodes == 0: continue
        ks = keystone_score(G)
        ks["module"] = ks["node"].map(mods)
        ks["dominance_class"] = cls
        # Mean abundance per genus in this stratum
        gen_ra = sub.div(sub.sum(axis=0).replace(0, 1), axis=1)
        ks["mean_relabund"] = ks["node"].map(gen_ra.mean(axis=1))
        ks.to_csv(OUT / f"keystone_{cls}_dominant.tsv", sep="\t",
                   index=False)

        # Top 10 keystones
        print(f"  Top 10 keystone genera in {cls}-dominant network:")
        print(ks.head(10)[["node", "degree", "betweenness", "closeness",
                              "keystone", "mean_relabund", "module"]]
              .round(4).to_string(index=False), flush=True)

        # Nibribacter status
        nib = "Nibribacter" in G.nodes
        if nib:
            nib_deg = G.degree("Nibribacter")
            nib_keystone = float(ks[ks["node"] == "Nibribacter"]
                                    ["keystone"].iloc[0])
            nib_rank = int((ks["keystone"] >= nib_keystone).sum())
            nib_neighbors = list(G.neighbors("Nibribacter"))
        else:
            nib_deg = 0; nib_keystone = 0.0; nib_rank = -1
            nib_neighbors = []
        print(f"\n  Nibribacter in {cls}-dominant network:")
        print(f"    Present: {nib}", flush=True)
        if nib:
            print(f"    Degree: {nib_deg}", flush=True)
            print(f"    Keystone score: {nib_keystone:.3f}", flush=True)
            print(f"    Rank: {nib_rank} of {n_nodes}", flush=True)
            print(f"    Neighbors[:10]: {nib_neighbors[:10]}", flush=True)

        # Save edges
        if n_edges > 0:
            ed = pd.DataFrame([(u, v, d.get("rho"), d.get("weight"),
                                  d.get("sign"), d.get("q"))
                                  for u, v, d in G.edges(data=True)],
                                 columns=["source", "target", "rho", "weight",
                                            "sign", "q"])
            ed.to_csv(OUT / f"edges_{cls}_dominant.tsv", sep="\t",
                       index=False)

        rec_summary.append({"dominance_class": cls,
                              "n_samples": len(samps),
                              "n_genera": n_nodes,
                              "n_edges": n_edges,
                              "n_modules": len(set(mods.values())),
                              "nibribacter_present": nib,
                              "nibribacter_degree": nib_deg,
                              "nibribacter_keystone_score": nib_keystone,
                              "nibribacter_rank": nib_rank})

    summary = pd.DataFrame(rec_summary)
    summary.to_csv(OUT / "network_summary.tsv", sep="\t", index=False)
    print("\n=== Summary ===")
    print(summary.to_string(index=False))

    # ----- Compare top keystones across A and B networks -----
    if Path(OUT / "keystone_A_dominant.tsv").exists() and \
       Path(OUT / "keystone_B_dominant.tsv").exists():
        ksA = pd.read_csv(OUT / "keystone_A_dominant.tsv", sep="\t")
        ksB = pd.read_csv(OUT / "keystone_B_dominant.tsv", sep="\t")
        print("\n=== Top-15 genera in A-dominant network ===")
        print(ksA.head(15)[["node", "degree", "keystone",
                              "mean_relabund"]].round(4).to_string(index=False))
        print("\n=== Top-15 genera in B-dominant network ===")
        print(ksB.head(15)[["node", "degree", "keystone",
                              "mean_relabund"]].round(4).to_string(index=False))
        # Overlap
        topA = set(ksA.head(20)["node"])
        topB = set(ksB.head(20)["node"])
        overlap = topA & topB
        print(f"\n  Top-20 keystone overlap A vs B: {len(overlap)} / 20")
        print(f"    Common: {sorted(overlap)}")
        print(f"    A-only: {sorted(topA - topB)}")
        print(f"    B-only: {sorted(topB - topA)}")


if __name__ == "__main__":
    main()
