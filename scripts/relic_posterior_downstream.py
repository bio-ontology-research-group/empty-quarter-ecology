#!/usr/bin/env python3
"""Re-run key population-level analyses with the prior-augmented posterior
score, AND compare each finding to a random-subset null.

Outputs:
  cache/relic_priors/posterior_climate_shannon.tsv
  cache/relic_priors/posterior_redundancy.tsv
  cache/relic_priors/posterior_temporal.tsv
  cache/relic_priors/posterior_csp_network.tsv
  cache/relic_priors/posterior_summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, braycurtis
from scipy.stats import spearmanr, linregress

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "src"))
from _sample_parse import parse_samples_to_df
from eq.network import compositional_correlation, build_network

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "relic_priors"
RNG = np.random.default_rng(20260510)


def shannon(arr):
    a = arr[arr > 0]
    if len(a) == 0: return 0.0
    p = a / a.sum()
    return float(-(p * np.log(p)).sum())


def relabund(M):
    return M.div(M.sum(axis=0).replace(0, 1), axis=1)


def main():
    print("Loading inputs ...", flush=True)
    p = pd.read_csv(OUT / "relic_score_with_priors.tsv", sep="\t")
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)

    geo_t1 = pd.read_csv(DATA / "geodata" / "trip1_geodata.tsv", sep="\t")
    geo_t1 = geo_t1.rename(columns={"Site": "site"})
    geo_t1["site"] = pd.to_numeric(geo_t1["site"],
                                       errors="coerce").astype("Int64")

    # Define posterior alive pool
    alive_p = set(p.loc[p["relic_score_posterior"] <= 0.3, "asv_id"])
    relic_p = set(p.loc[p["relic_score_posterior"] >= 0.7, "asv_id"])
    alive_m = set(p.loc[p["relic_score_full_gb"] <= 0.3, "asv_id"])
    relic_m = set(p.loc[p["relic_score_full_gb"] >= 0.7, "asv_id"])
    print(f"  posterior alive: {len(alive_p)}, relic: {len(relic_p)}",
          flush=True)
    print(f"  model-only alive: {len(alive_m)}, relic: {len(relic_m)}",
          flush=True)

    # ----- Climate-Shannon: posterior alive vs random subset of same size -----
    print("\n[A] Climate-Shannon (alive_post) ...", flush=True)
    def climate_rho(asv_subset):
        sub = ft.loc[ft.index.isin(asv_subset)]
        sh = pd.DataFrame({"sample": sub.columns,
                              "shannon": [shannon(sub[c].values)
                                            for c in sub.columns]})
        sh = sh.merge(smeta, on="sample", how="left")
        sh = sh.merge(geo_t1[["site", "AnnualMeanTemp"]],
                        on="site", how="left").dropna(
                            subset=["AnnualMeanTemp"])
        if len(sh) < 30: return np.nan
        r, _ = spearmanr(sh["shannon"], sh["AnnualMeanTemp"])
        return float(r), float(sh["shannon"].median())

    r_p, med_p = climate_rho(alive_p)
    print(f"  posterior alive ({len(alive_p)} ASVs):  rho={r_p:+.3f}",
          flush=True)
    r_a, med_a = climate_rho(alive_m)
    print(f"  model-only alive ({len(alive_m)}):     rho={r_a:+.3f}",
          flush=True)
    r_all, med_all = climate_rho(set(ft.index))
    print(f"  all (75k):                  rho={r_all:+.3f}", flush=True)

    # Random null at posterior-alive size
    n_iter = 100
    rand = []
    for _ in range(n_iter):
        sub = set(RNG.choice(ft.index.tolist(),
                                len(alive_p), replace=False))
        r, _ = climate_rho(sub)
        if not np.isnan(r): rand.append(r)
    rand = np.array(rand)
    print(f"  random-null at n={len(alive_p)}:", flush=True)
    print(f"    median: {np.median(rand):+.3f}", flush=True)
    print(f"    p5..p95: [{np.percentile(rand, 5):+.3f}, "
          f"{np.percentile(rand, 95):+.3f}]", flush=True)
    p_extreme = float((rand > r_p).mean())
    print(f"    posterior alive vs random: p={1-p_extreme:.3f} for "
          f"more-negative", flush=True)

    # ----- Allison-Martiny redundancy on posterior alive vs random -----
    print("\n[B] Allison-Martiny on posterior alive ...", flush=True)
    path_func = pd.read_csv(DATA / "functional" / "picrust2" /
                              "path_abun_unstrat.tsv",
                              sep="\t", index_col=0)
    rsamps = list(set(smeta[(smeta["compartment"] == "rhizosphere") &
                                (smeta["trip"] == 3)]["sample"]) &
                    set(ft.columns) & set(path_func.columns))
    Fr_ref = relabund(path_func[rsamps]).T
    bc_func = pdist(Fr_ref.values, metric="braycurtis")

    def am_slope(asv_subset):
        sub = ft.loc[ft.index.isin(asv_subset), rsamps]
        Tr = relabund(sub).T
        Tr = Tr[Tr.sum(axis=1) > 0]
        s2 = list(Tr.index)
        if len(s2) < 10: return np.nan
        bc_t = pdist(Tr.values, metric="braycurtis")
        Fr = Fr_ref.loc[s2]
        bc_f = pdist(Fr.values, metric="braycurtis")
        m = (~np.isnan(bc_t)) & (~np.isnan(bc_f))
        if m.sum() < 30: return np.nan
        slope, _, _, _, _ = linregress(bc_t[m], bc_f[m])
        return float(slope)

    s_p = am_slope(alive_p)
    s_m = am_slope(alive_m)
    s_all = am_slope(set(ft.index))
    print(f"  posterior alive:    slope={s_p:.3f}", flush=True)
    print(f"  model-only alive:   slope={s_m:.3f}", flush=True)
    print(f"  all:                slope={s_all:.3f}", flush=True)

    rand = []
    for _ in range(n_iter):
        sub = set(RNG.choice(ft.index.tolist(),
                                len(alive_p), replace=False))
        s = am_slope(sub)
        if not np.isnan(s): rand.append(s)
    rand = np.array(rand)
    print(f"  random-null slope at n={len(alive_p)}:", flush=True)
    print(f"    median: {np.median(rand):.3f}, "
          f"p5..p95: [{np.percentile(rand, 5):.3f}, "
          f"{np.percentile(rand, 95):.3f}]", flush=True)
    p_low = float((rand < s_p).mean())
    print(f"    posterior alive vs random: p={p_low:.3f} for less-positive",
          flush=True)

    # ----- All-pairwise temporal BC -----
    print("\n[C] All-pairwise temporal BC on posterior pools ...", flush=True)
    def temporal_bc(asv_subset):
        sub = ft.loc[ft.index.isin(asv_subset)]
        Mr = relabund(sub)
        rec = []
        for site in sorted(smeta["site"].unique()):
            for comp in ("rhizosphere", "surface", "deep"):
                samps_per_trip = {}
                for t in (1, 2, 3, 4, 5):
                    ss = list(set(smeta[(smeta["site"] == site) &
                                              (smeta["compartment"] == comp) &
                                              (smeta["trip"] == t)]["sample"]) &
                                set(Mr.columns))
                    if ss:
                        v = Mr[ss].mean(axis=1).values
                        if v.sum() > 0: samps_per_trip[t] = v
                if len(samps_per_trip) < 2: continue
                bcs = []
                trips_l = list(samps_per_trip)
                for i, ti in enumerate(trips_l):
                    for tj in trips_l[i+1:]:
                        bcs.append(braycurtis(samps_per_trip[ti],
                                                  samps_per_trip[tj]))
                if bcs:
                    rec.append(float(np.mean(bcs)))
        return float(np.median(rec)) if rec else np.nan

    t_p = temporal_bc(alive_p)
    t_m = temporal_bc(alive_m)
    t_all = temporal_bc(set(ft.index))
    print(f"  posterior alive: median mean_pairwise_bc = {t_p:.3f}",
          flush=True)
    print(f"  model-only alive:                          {t_m:.3f}",
          flush=True)
    print(f"  all:                                       {t_all:.3f}",
          flush=True)

    rand = []
    for _ in range(20):  # less iterations - this is expensive
        sub = set(RNG.choice(ft.index.tolist(),
                                len(alive_p), replace=False))
        t = temporal_bc(sub)
        if not np.isnan(t): rand.append(t)
    rand = np.array(rand)
    print(f"  random-null at n={len(alive_p)} (20 iter):", flush=True)
    print(f"    median: {np.median(rand):.3f}, "
          f"p5..p95: [{np.percentile(rand, 5):.3f}, "
          f"{np.percentile(rand, 95):.3f}]", flush=True)
    p_low = float((rand < t_p).mean())
    print(f"    posterior alive vs random: p={p_low:.3f} for less-than",
          flush=True)

    # ----- CSP1-2 in posterior network -----
    print("\n[D] CSP1-2 in posterior alive network ...", flush=True)
    tax = pd.read_parquet(CACHE / "taxonomy.parquet").reset_index().rename(
        columns={"ASV": "asv_id"})
    ft_alive_p = ft.loc[ft.index.isin(alive_p)]
    ft2 = ft_alive_p.copy()
    ft2.index = ft2.index.rename("asv_id")
    m = ft2.reset_index().merge(tax[["asv_id", "genus"]], on="asv_id",
                                  how="left")
    m = m.dropna(subset=["genus"])
    sample_cols = [c for c in m.columns if c not in ("asv_id", "genus")]
    gen = m.groupby("genus")[sample_cols].sum()
    rec = []
    for comp in ("rhizosphere", "surface", "deep"):
        samps = list(set(smeta[smeta["compartment"] == comp]["sample"]) &
                       set(gen.columns))
        if len(samps) < 30: continue
        sub = gen[samps]
        try:
            rho, pval = compositional_correlation(sub, min_prevalence=0.10,
                                                       presence_ra=0.0001,
                                                       pseudo=0.5)
            G = build_network(rho, pval, rho_threshold=0.4, q_threshold=0.05)
        except Exception as e:
            continue
        csp_node = next((n for n in G.nodes
                              if isinstance(n, str) and "CSP1" in n), None)
        rec.append({"compartment": comp,
                      "n_genera": G.number_of_nodes(),
                      "n_edges": G.number_of_edges(),
                      "csp_present": csp_node is not None,
                      "csp_degree": (G.degree(csp_node)
                                       if csp_node else 0)})
        print(f"  [{comp}] genera={G.number_of_nodes()} edges={G.number_of_edges()}"
              f"  CSP1-2: {'present deg='+str(G.degree(csp_node)) if csp_node else 'absent'}",
              flush=True)
    pd.DataFrame(rec).to_csv(OUT / "posterior_csp_network.tsv",
                                sep="\t", index=False)


if __name__ == "__main__":
    main()
