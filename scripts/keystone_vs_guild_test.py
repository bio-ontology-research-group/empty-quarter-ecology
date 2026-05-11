#!/usr/bin/env python3
"""Discriminating test: single-keystone (Nibribacter) vs three-guild
hypothesis for the EQ alive community.

Tests:
  A. Nibribacter-stratified samples: in samples where Nibribacter is
     low/absent, does the alive community structure collapse (single
     keystone) or do other Bacteroidota take over (guild)?
  B. Knockout-robustness: simulate removal of Nibribacter alone vs random
     single-node removal vs guild-batch removal.
  C. Within-guild functional redundancy: pairwise PICRUSt2-KO BC
     similarity among putative guild members.
  D. Substitution test: at sample level, do other Bacteroidota
     (Cytophaga, Daejeonella) anti-correlate with Nibribacter (substitution)
     or correlate (co-living)?
  E. Module-guild correspondence: do Louvain modules align with phyla?
  F. Nibribacter abundance vs XRF + climate + nibribacter literature focus.
  G. PICRUSt2 functional fingerprint of Nibribacter.

Outputs in cache/keystone_test/.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial.distance import braycurtis, pdist, squareform
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "src"))
from _sample_parse import parse_samples_to_df
from eq.network import (
    compositional_correlation, build_network, louvain_modules, keystone_score
)

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "keystone_test"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260510)

PICRUSt2_KO_PATH = Path("/home/leechuck/Public/software/empty-quarter/data/"
                          "processed/functional/picrust2/KO_predicted.tsv")


def relabund(M):
    return M.div(M.sum(axis=0).replace(0, 1), axis=1)


def shannon(arr):
    a = arr[arr > 0]
    if len(a) == 0: return 0.0
    p = a / a.sum()
    return float(-(p * np.log(p)).sum())


def main():
    print("Loading inputs ...", flush=True)
    p = pd.read_csv(CACHE / "relic_priors" /
                     "relic_score_with_mag_prior.tsv", sep="\t")
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)
    tax = pd.read_parquet(CACHE / "taxonomy.parquet").reset_index().rename(
        columns={"ASV": "asv_id"})

    alive_mag = set(p.loc[p["relic_score_with_mag"] <= 0.3, "asv_id"])
    ft_a = ft.loc[ft.index.isin(alive_mag)]

    # Aggregate to genus
    ft2 = ft_a.copy()
    ft2.index = ft2.index.rename("asv_id")
    m = ft2.reset_index().merge(tax[["asv_id", "phylum", "genus"]],
                                  on="asv_id", how="left")
    m = m.dropna(subset=["genus"])
    m = m[~m["genus"].astype(str).isin(["NA", ""])]
    sample_cols = [c for c in m.columns if c not in
                       ("asv_id", "phylum", "genus")]
    gen = m.groupby("genus")[sample_cols].sum()
    gen_phyl = (m[["genus", "phylum"]].drop_duplicates("genus")
                 .set_index("genus")["phylum"])

    # Define guilds (top candidates from earlier analysis)
    BACT_DOM_CYCLERS = ["Nibribacter", "Daejeonella", "Cytophaga", "Cnuella",
                          "Niastella", "Flavitalea", "Ohtaekwangia",
                          "Tellurirhabdus"]
    HALOTOL_BACILLI = ["Aquibacillus", "Sediminibacillus", "Litchfieldia",
                         "Tumebacillus", "Gracilibacillus", "Oceanobacillus",
                         "Polygonibacillus", "Salirhabdus",
                         "Alkalihalobacillus"]
    PSEUDO_GENERALISTS = ["Massilia", "Halomonas", "Acidibacter", "Lysobacter",
                            "Methylobacillus"]
    GUILDS = {"Bact_DOM": BACT_DOM_CYCLERS,
                 "Bacilli_halo": HALOTOL_BACILLI,
                 "Pseudo_gen": PSEUDO_GENERALISTS}

    # ======================================================================
    # A. Nibribacter-stratified samples
    # ======================================================================
    print("\n========================================================")
    print("[A] Nibribacter-stratified samples")
    print("========================================================")
    if "Nibribacter" in gen.index:
        nib_per_sample = relabund(gen).loc["Nibribacter"]
    else:
        print("Nibribacter not in gen index"); return
    # Stratify into 4 quartiles
    q = nib_per_sample.quantile([0.25, 0.5, 0.75]).values
    def stratum(x):
        if x <= q[0]: return "Q1_low"
        if x <= q[1]: return "Q2"
        if x <= q[2]: return "Q3"
        return "Q4_high"
    sample_strata = nib_per_sample.apply(stratum)
    print(f"  Nibribacter per-sample relabund quartiles: "
          f"Q25={q[0]:.4f}, Q50={q[1]:.4f}, Q75={q[2]:.4f}", flush=True)

    # Per stratum: compute alive richness, Shannon, n_modules in alive network
    # of stratum samples
    stratum_metrics = []
    for s in ("Q1_low", "Q2", "Q3", "Q4_high"):
        samps = sample_strata[sample_strata == s].index.tolist()
        if len(samps) < 30: continue
        # alive-only counts at these samples
        sub = ft_a[samps]
        alive_richness = (sub > 0).sum(axis=0).median()
        alive_shannon = float(np.median([shannon(sub[c].values)
                                           for c in sub.columns]))
        # Bacteroidota richness specifically (the guild)
        bact_asvs = m[m["phylum"] == "Bacteroidota"]["asv_id"]
        bact_in_ft = ft_a.index.isin(bact_asvs)
        bact_richness = (ft_a[samps].loc[bact_in_ft] > 0).sum(axis=0).median()
        # Bacilli + Pseudomonadota
        bacill_asvs = m[m["phylum"] == "Bacillota"]["asv_id"]
        bacill_in_ft = ft_a.index.isin(bacill_asvs)
        bacill_rich = (ft_a[samps].loc[bacill_in_ft] > 0).sum(axis=0).median()
        # Per-guild relabund
        sub_gen = gen[samps]
        sub_gen_rel = relabund(sub_gen)
        guild_means = {}
        for gname, glist in GUILDS.items():
            present = [g for g in glist if g in sub_gen_rel.index]
            if present:
                guild_means[gname] = float(sub_gen_rel.loc[present]
                                              .sum(axis=0).median())
            else:
                guild_means[gname] = 0
        # Network coherence: build network on these samples only, get
        # n_edges / max possible
        try:
            rho, pval = compositional_correlation(sub_gen,
                                                       min_prevalence=0.10,
                                                       presence_ra=0.0001,
                                                       pseudo=0.5)
            G = build_network(rho, pval, rho_threshold=0.4,
                                q_threshold=0.05)
            n_genera = G.number_of_nodes()
            n_edges = G.number_of_edges()
            n_modules = (len(set(louvain_modules(G).values()))
                            if G.number_of_edges() else 0)
        except Exception as e:
            n_genera = n_edges = n_modules = -1
        stratum_metrics.append({
            "stratum": s, "n_samples": len(samps),
            "median_alive_richness": int(alive_richness),
            "median_alive_shannon": alive_shannon,
            "median_Bacteroidota_richness": int(bact_richness),
            "median_Bacillota_richness": int(bacill_rich),
            "guild_BactDOM_relabund": guild_means["Bact_DOM"],
            "guild_BacilliHalo_relabund": guild_means["Bacilli_halo"],
            "guild_PseudoGen_relabund": guild_means["Pseudo_gen"],
            "network_n_genera": n_genera,
            "network_n_edges": n_edges,
            "n_modules": n_modules,
        })
    sm = pd.DataFrame(stratum_metrics)
    sm.to_csv(OUT / "nibribacter_stratified_metrics.tsv", sep="\t",
                index=False)
    print(sm.round(4).to_string(index=False), flush=True)

    # KEY TEST: in low-Nib samples, does Bact_DOM guild relabund stay HIGH (sub-
    # stitution by other Bacteroidota = guild) or DROP (single keystone)?
    if len(sm) == 4:
        low = sm[sm["stratum"] == "Q1_low"].iloc[0]
        high = sm[sm["stratum"] == "Q4_high"].iloc[0]
        print(f"\n  Q1_low Bact_DOM guild relabund = "
              f"{low['guild_BactDOM_relabund']:.4f}", flush=True)
        print(f"  Q4_high Bact_DOM guild relabund = "
              f"{high['guild_BactDOM_relabund']:.4f}", flush=True)
        # If guild: low should have substantial guild relabund (other Bact)
        # If keystone: low should have very low guild relabund (no substitute)

    # ======================================================================
    # B. Knockout robustness
    # ======================================================================
    print("\n========================================================")
    print("[B] Knockout robustness")
    print("========================================================")
    rec_b = []
    for comp in ("rhizosphere", "surface", "deep"):
        samps = list(set(smeta[smeta["compartment"] == comp]["sample"]) &
                       set(gen.columns))
        if len(samps) < 30: continue
        sub = gen[samps]
        try:
            rho, pval = compositional_correlation(sub, min_prevalence=0.10,
                                                       presence_ra=0.0001,
                                                       pseudo=0.5)
            G = build_network(rho, pval, rho_threshold=0.4,
                                q_threshold=0.05)
        except Exception as e:
            print(f"  [{comp}] failed: {e}"); continue
        if G.number_of_nodes() == 0: continue

        def network_metrics(H):
            return {
                "n_nodes": H.number_of_nodes(),
                "n_edges": H.number_of_edges(),
                "n_components": nx.number_connected_components(H),
                "largest_component": (max(len(c) for c in
                                            nx.connected_components(H))
                                         if H.number_of_nodes() > 0 else 0),
                "avg_clustering": (float(nx.average_clustering(H))
                                       if H.number_of_edges() else 0),
            }

        # Baseline
        base_m = network_metrics(G)

        scenarios = []
        # Knockout 1: Nibribacter alone (if present)
        if "Nibribacter" in G:
            H = G.copy(); H.remove_node("Nibribacter")
            scenarios.append(("Nibribacter_alone", network_metrics(H)))
        # Knockout 2: Bact_DOM guild members present
        present = [g for g in BACT_DOM_CYCLERS if g in G]
        if present:
            H = G.copy()
            for g in present: H.remove_node(g)
            scenarios.append((f"Bact_DOM_guild_n{len(present)}",
                                network_metrics(H)))
        # Knockout 3: Bacilli guild
        present = [g for g in HALOTOL_BACILLI if g in G]
        if present:
            H = G.copy()
            for g in present: H.remove_node(g)
            scenarios.append((f"Bacilli_guild_n{len(present)}",
                                network_metrics(H)))
        # Knockout 4: Pseudomonadota guild
        present = [g for g in PSEUDO_GENERALISTS if g in G]
        if present:
            H = G.copy()
            for g in present: H.remove_node(g)
            scenarios.append((f"Pseudo_guild_n{len(present)}",
                                network_metrics(H)))
        # Knockout 5: random single - average over 50 trials
        node_list = list(G.nodes)
        if "Nibribacter" in node_list:
            node_list.remove("Nibribacter")
        random_metrics = {"n_nodes": [], "n_edges": [], "n_components": [],
                            "largest_component": [], "avg_clustering": []}
        for _ in range(50):
            H = G.copy()
            H.remove_node(RNG.choice(node_list))
            for k, v in network_metrics(H).items():
                random_metrics[k].append(v)
        random_avg = {k: float(np.mean(v)) for k, v in random_metrics.items()}
        scenarios.append(("random_single_avg50", random_avg))
        # Knockout 6: random N (where N = max guild size = 9)
        if len(node_list) >= 9:
            random_n_metrics = {"n_nodes": [], "n_edges": [],
                                  "n_components": [], "largest_component": [],
                                  "avg_clustering": []}
            for _ in range(50):
                H = G.copy()
                rem = RNG.choice(node_list, 9, replace=False)
                for r in rem: H.remove_node(r)
                for k, v in network_metrics(H).items():
                    random_n_metrics[k].append(v)
            random_n_avg = {k: float(np.mean(v))
                              for k, v in random_n_metrics.items()}
            scenarios.append(("random_n9_avg50", random_n_avg))

        # Print
        print(f"\n  [{comp}] BASELINE: nodes={base_m['n_nodes']} "
              f"edges={base_m['n_edges']} ncomp={base_m['n_components']} "
              f"largest={base_m['largest_component']} "
              f"clustering={base_m['avg_clustering']:.3f}", flush=True)
        for name, met in scenarios:
            edge_loss = (base_m['n_edges'] - met['n_edges']) / max(base_m['n_edges'], 1)
            print(f"    {name:<35}: nodes={met['n_nodes']} "
                  f"edges={met['n_edges']:.0f} "
                  f"(-{edge_loss:.1%}) "
                  f"comp={met['n_components']:.1f} "
                  f"largest={met['largest_component']:.1f}", flush=True)
            rec_b.append({"compartment": comp, "scenario": name,
                            "n_nodes_after": met["n_nodes"],
                            "n_edges_after": met["n_edges"],
                            "n_components_after": met["n_components"],
                            "edge_loss_pct": edge_loss * 100,
                            "n_components_baseline": base_m["n_components"]})
    pd.DataFrame(rec_b).to_csv(OUT / "knockout_robustness.tsv",
                                  sep="\t", index=False)

    # ======================================================================
    # C. Within-guild functional redundancy via PICRUSt2 KO profiles
    # ======================================================================
    print("\n========================================================")
    print("[C] Within-guild functional redundancy")
    print("========================================================")
    print("  Loading per-ASV PICRUSt2 KO predictions ...", flush=True)
    if PICRUSt2_KO_PATH.exists():
        ko_pred = pd.read_csv(PICRUSt2_KO_PATH, sep="\t", index_col=0,
                                dtype={"sequence": str})
        ko_pred = ko_pred.astype(np.float32)
        # Per-genus KO profile = mean over its ASVs
        # Map ASV -> genus
        asv_to_genus = (m[["asv_id", "genus"]].set_index("asv_id")
                          ["genus"].to_dict())
        # Restrict ko_pred to alive ASVs
        common = ko_pred.index.intersection(set(asv_to_genus))
        ko_a = ko_pred.loc[common]
        gen_for_ko = pd.Series([asv_to_genus[a] for a in common],
                                  index=common)
        # Per-genus mean KO copies
        ko_per_genus = ko_a.groupby(gen_for_ko).mean()
        # Restrict to KOs present (>0) in at least one genus
        ko_per_genus = ko_per_genus.loc[:, ko_per_genus.sum(axis=0) > 0]
        # Normalize to relabund per genus
        ko_per_genus_rel = ko_per_genus.div(ko_per_genus.sum(axis=1)
                                                .replace(0, 1), axis=0)
        # For each guild: pairwise BC similarity
        for gname, glist in GUILDS.items():
            present = [g for g in glist if g in ko_per_genus_rel.index]
            if len(present) < 2: continue
            sub = ko_per_genus_rel.loc[present]
            D = squareform(pdist(sub.values, metric="braycurtis"))
            # Mean off-diagonal BC
            iu = np.triu_indices(len(present), k=1)
            mean_bc = float(D[iu].mean())
            print(f"  {gname} ({len(present)} genera): mean within-guild BC "
                  f"= {mean_bc:.3f}  (lower = more redundant)", flush=True)
            for i, gi in enumerate(present):
                for j, gj in enumerate(present):
                    if j <= i: continue
                    print(f"    {gi:<20} vs {gj:<20}: BC={D[i,j]:.3f}",
                          flush=True)

        # Compare to between-guild BC
        print("\n  Cross-guild BC (members of different guilds):", flush=True)
        for ga, la in GUILDS.items():
            for gb, lb in GUILDS.items():
                if ga >= gb: continue
                pa = [g for g in la if g in ko_per_genus_rel.index]
                pb = [g for g in lb if g in ko_per_genus_rel.index]
                if not pa or not pb: continue
                bcs = []
                for x in pa:
                    for y in pb:
                        bcs.append(braycurtis(ko_per_genus_rel.loc[x],
                                                  ko_per_genus_rel.loc[y]))
                print(f"    {ga} vs {gb}: mean BC = {float(np.mean(bcs)):.3f}",
                      flush=True)
    else:
        print(f"  KO file not found: {PICRUSt2_KO_PATH}")

    # ======================================================================
    # D. Substitution test: do Bacteroidota anti-correlate when Nib is low?
    # ======================================================================
    print("\n========================================================")
    print("[D] Substitution test")
    print("========================================================")
    nib = relabund(gen).loc["Nibribacter"]
    rec_d = []
    for partner in BACT_DOM_CYCLERS + HALOTOL_BACILLI[:3] + ["Massilia",
                                                                    "Halomonas"]:
        if partner == "Nibribacter": continue
        if partner not in gen.index: continue
        v = relabund(gen).loc[partner]
        common = nib.index.intersection(v.index)
        r, p_val = spearmanr(nib.loc[common], v.loc[common])
        rec_d.append({"partner": partner,
                        "spearman_with_Nibribacter": float(r),
                        "p": float(p_val), "n": len(common)})
    df_d = pd.DataFrame(rec_d).sort_values("spearman_with_Nibribacter")
    df_d.to_csv(OUT / "substitution_correlations.tsv", sep="\t", index=False)
    print(df_d.round(3).to_string(index=False), flush=True)

    # ======================================================================
    # E. Module-guild correspondence
    # ======================================================================
    print("\n========================================================")
    print("[E] Module-guild (phylum) correspondence")
    print("========================================================")
    rec_e = []
    for comp in ("rhizosphere", "surface", "deep"):
        samps = list(set(smeta[smeta["compartment"] == comp]["sample"]) &
                       set(gen.columns))
        if len(samps) < 30: continue
        sub = gen[samps]
        try:
            rho, pval = compositional_correlation(sub, min_prevalence=0.10,
                                                       presence_ra=0.0001,
                                                       pseudo=0.5)
            G = build_network(rho, pval, rho_threshold=0.4,
                                q_threshold=0.05)
            mods = louvain_modules(G, seed=42)
        except Exception as e:
            continue
        node_phylum = {n: gen_phyl.get(n, "NA") for n in G.nodes}
        # Per module: phylum composition
        mod_phyl_counts = {}
        for n, mid in mods.items():
            mod_phyl_counts.setdefault(mid, {}).setdefault(
                node_phylum[n], 0)
            mod_phyl_counts[mid][node_phylum[n]] += 1
        for mid, pcount in mod_phyl_counts.items():
            total = sum(pcount.values())
            if total < 3: continue
            top_phylum = max(pcount, key=pcount.get)
            top_frac = pcount[top_phylum] / total
            rec_e.append({"compartment": comp, "module": mid,
                            "module_size": total,
                            "top_phylum": top_phylum,
                            "top_phylum_frac": top_frac,
                            "n_phyla": len(pcount)})
    df_e = pd.DataFrame(rec_e)
    df_e.to_csv(OUT / "module_phylum_composition.tsv", sep="\t", index=False)
    print(df_e.round(3).to_string(index=False), flush=True)

    # Compute mean top_phylum_frac (high = phylum-pure modules = guild
    # structure)
    print(f"\n  Mean module top-phylum fraction: "
          f"{df_e['top_phylum_frac'].mean():.3f}", flush=True)
    print(f"  Random expectation (uniform): "
          f"{1/df_e['n_phyla'].mean():.3f}", flush=True)


if __name__ == "__main__":
    main()
