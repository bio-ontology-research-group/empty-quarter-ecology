#!/usr/bin/env python3
"""Re-run key downstream analyses using the MAG-augmented relic indicator.

Tests:
  1. CSP1-2 status in MAG-augmented alive network
  2. Climate-Shannon with MAG-augmented alive pool + random null
  3. All-pairwise temporal BC + null
  4. Top alive genera + functional decomposition

Outputs in cache/relic_priors/mag_downstream/.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, braycurtis
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "src"))
from _sample_parse import parse_samples_to_df
from eq.network import compositional_correlation, build_network

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "relic_priors" / "mag_downstream"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260510)


def shannon(arr):
    a = arr[arr > 0]
    if len(a) == 0: return 0.0
    p = a / a.sum()
    return float(-(p * np.log(p)).sum())


def relabund(M):
    return M.div(M.sum(axis=0).replace(0, 1), axis=1)


def main():
    p = pd.read_csv(CACHE / "relic_priors" /
                     "relic_score_with_mag_prior.tsv", sep="\t")
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)

    alive_mag = set(p.loc[p["relic_score_with_mag"] <= 0.3, "asv_id"])
    relic_mag = set(p.loc[p["relic_score_with_mag"] >= 0.7, "asv_id"])
    alive_post = set(p.loc[p["relic_score_posterior"] <= 0.3, "asv_id"])
    print(f"alive (MAG-aug): {len(alive_mag)}", flush=True)
    print(f"alive (post only): {len(alive_post)}", flush=True)
    print(f"relic (MAG-aug): {len(relic_mag)}", flush=True)

    tax = pd.read_parquet(CACHE / "taxonomy.parquet").reset_index().rename(
        columns={"ASV": "asv_id"})

    # ----- 1. CSP1-2 in MAG-augmented network -----
    print("\n[1] CSP1-2 in MAG-augmented alive network ...", flush=True)
    ft_a = ft.loc[ft.index.isin(alive_mag)]
    ft2 = ft_a.copy()
    ft2.index = ft2.index.rename("asv_id")
    m = ft2.reset_index().merge(tax[["asv_id", "genus"]], on="asv_id",
                                  how="left").dropna(subset=["genus"])
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
            print(f"  [{comp}] failed: {e}"); continue
        csp_node = next((n for n in G.nodes
                              if isinstance(n, str) and "CSP1" in n), None)
        if csp_node:
            deg = G.degree(csp_node)
            ne = list(G.neighbors(csp_node))
            print(f"  [{comp}] genera={G.number_of_nodes()} edges={G.number_of_edges()}"
                  f"  CSP1-2: PRESENT deg={deg}  neighbors[:5]={ne[:5]}",
                  flush=True)
        else:
            deg = 0
            print(f"  [{comp}] genera={G.number_of_nodes()} edges={G.number_of_edges()}"
                  f"  CSP1-2: ABSENT", flush=True)
        rec.append({"compartment": comp,
                      "n_genera": G.number_of_nodes(),
                      "n_edges": G.number_of_edges(),
                      "csp_present": csp_node is not None,
                      "csp_degree": deg})
    pd.DataFrame(rec).to_csv(OUT / "csp_in_mag_alive_network.tsv",
                                sep="\t", index=False)

    # ----- 2. Climate-Shannon -----
    geo_t1 = pd.read_csv(DATA / "geodata" / "trip1_geodata.tsv", sep="\t")
    geo_t1 = geo_t1.rename(columns={"Site": "site"})
    geo_t1["site"] = pd.to_numeric(geo_t1["site"],
                                       errors="coerce").astype("Int64")

    def climate_rho(asv_subset):
        sub = ft.loc[ft.index.isin(asv_subset)]
        sh = pd.DataFrame({"sample": sub.columns,
                              "shannon": [shannon(sub[c].values)
                                            for c in sub.columns]})
        sh = sh.merge(smeta, on="sample", how="left")
        sh = sh.merge(geo_t1[["site", "AnnualMeanTemp"]], on="site",
                        how="left").dropna(subset=["AnnualMeanTemp"])
        if len(sh) < 30: return np.nan
        r, _ = spearmanr(sh["shannon"], sh["AnnualMeanTemp"])
        return float(r)

    print(f"\n[2] Climate-Shannon ...", flush=True)
    r_mag = climate_rho(alive_mag)
    r_post = climate_rho(alive_post)
    r_all = climate_rho(set(ft.index))
    print(f"  alive (MAG-aug, n={len(alive_mag)}): rho={r_mag:+.3f}",
          flush=True)
    print(f"  alive (post only, n={len(alive_post)}): rho={r_post:+.3f}",
          flush=True)
    print(f"  all (n={ft.shape[0]}): rho={r_all:+.3f}", flush=True)

    rand = []
    for _ in range(100):
        sub = set(RNG.choice(ft.index.tolist(),
                                len(alive_mag), replace=False))
        rr = climate_rho(sub)
        if not np.isnan(rr): rand.append(rr)
    rand = np.array(rand)
    print(f"  random-null at n={len(alive_mag)}: median={np.median(rand):+.3f}"
          f"  p5..p95: [{np.percentile(rand, 5):+.3f}, "
          f"{np.percentile(rand, 95):+.3f}]", flush=True)
    p_extreme = float((rand > r_mag).mean())
    print(f"  alive (MAG) vs random: p={1-p_extreme:.3f} for more-negative",
          flush=True)

    # ----- 3. Temporal BC -----
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
                                              (smeta["trip"] == t)]["sample"])
                                & set(Mr.columns))
                    if ss:
                        v = Mr[ss].mean(axis=1).values
                        if v.sum() > 0: samps_per_trip[t] = v
                if len(samps_per_trip) < 2: continue
                bcs = []
                tl = list(samps_per_trip)
                for i, ti in enumerate(tl):
                    for tj in tl[i+1:]:
                        bcs.append(braycurtis(samps_per_trip[ti],
                                                  samps_per_trip[tj]))
                if bcs: rec.append(float(np.mean(bcs)))
        return float(np.median(rec)) if rec else np.nan

    print(f"\n[3] Temporal BC ...", flush=True)
    t_mag = temporal_bc(alive_mag)
    t_post = temporal_bc(alive_post)
    t_all = temporal_bc(set(ft.index))
    print(f"  alive (MAG-aug):    {t_mag:.3f}", flush=True)
    print(f"  alive (post only):  {t_post:.3f}", flush=True)
    print(f"  all:                {t_all:.3f}", flush=True)
    rand = []
    for _ in range(20):
        sub = set(RNG.choice(ft.index.tolist(),
                                len(alive_mag), replace=False))
        tt = temporal_bc(sub)
        if not np.isnan(tt): rand.append(tt)
    rand = np.array(rand)
    print(f"  random-null (n={len(alive_mag)}, 20 iter):", flush=True)
    print(f"    median={np.median(rand):.3f}  p5..p95: "
          f"[{np.percentile(rand, 5):.3f}, "
          f"{np.percentile(rand, 95):.3f}]", flush=True)
    p_low = float((rand < t_mag).mean())
    print(f"    alive (MAG) vs random: p={p_low:.3f} for less-than",
          flush=True)


if __name__ == "__main__":
    main()
