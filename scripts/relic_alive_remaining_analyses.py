#!/usr/bin/env python3
"""Remaining alive-only re-analyses:

  1. ALL-PAIRWISE temporal BC (fixes the T5-vs-rest framing).
     For each (site, comp), compute mean BC across all 10 unique
     trip pairs, plus per-pair-of-trips BC.

  2. Genus-level co-occurrence network (alive ASVs only): CSP1-2
     keystone status, modularity, edge counts. Per compartment.

  3. iCAMP-style RCbray on alive vs all vs relic. Vectorized,
     per compartment, 99 perms (lighter than original 999).

  4. Mediation: salinity → CSP1-2 → alive_Shannon vs all_Shannon.
     Compare ACME/ADE shares.

Outputs in cache/relic_alive_subset/.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform, braycurtis
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "src"))
from _sample_parse import parse_samples_to_df
from eq.network import (
    compositional_correlation, build_network, louvain_modules
)

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "relic_alive_subset"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260510)


def relabund(M):
    return M.div(M.sum(axis=0).replace(0, 1), axis=1)


def shannon(arr):
    a = arr[arr > 0]
    if len(a) == 0: return 0.0
    p = a / a.sum()
    return float(-(p * np.log(p)).sum())


# =============================================================================
# 1. All-pairwise temporal BC
# =============================================================================
def all_pairwise_temporal(ft, label, smeta):
    print(f"\n=== [1] All-pairwise temporal BC ({label}) ===", flush=True)
    od = OUT / f"temporal_pairwise_{label}"; od.mkdir(exist_ok=True)
    Mr = relabund(ft)
    pair_rows = []; mean_rows = []
    trips = (1, 2, 3, 4, 5)
    for site in sorted(smeta["site"].unique()):
        for comp in ("rhizosphere", "surface", "deep"):
            samps_per_trip = {}
            for t in trips:
                ss = list(set(smeta[(smeta["site"] == site) &
                                       (smeta["compartment"] == comp) &
                                       (smeta["trip"] == t)]["sample"]) &
                            set(Mr.columns))
                if ss:
                    samps_per_trip[t] = Mr[ss].mean(axis=1).values
            if len(samps_per_trip) < 2: continue
            pair_bcs = []
            for ti in samps_per_trip:
                for tj in samps_per_trip:
                    if tj <= ti: continue
                    if samps_per_trip[ti].sum() == 0: continue
                    if samps_per_trip[tj].sum() == 0: continue
                    bc = braycurtis(samps_per_trip[ti], samps_per_trip[tj])
                    pair_rows.append({"site": site, "compartment": comp,
                                       "trip_a": ti, "trip_b": tj,
                                       "bc": float(bc)})
                    pair_bcs.append(float(bc))
            if not pair_bcs: continue
            mean_rows.append({"site": site, "compartment": comp,
                                "n_trips": len(samps_per_trip),
                                "n_pairs": len(pair_bcs),
                                "mean_pairwise_bc": float(np.mean(pair_bcs)),
                                "median_pairwise_bc": float(np.median(pair_bcs))})
    pair_df = pd.DataFrame(pair_rows)
    mean_df = pd.DataFrame(mean_rows)
    pair_df.to_csv(od / "per_trip_pair_bc.tsv", sep="\t", index=False)
    mean_df.to_csv(od / "per_site_comp_mean_bc.tsv", sep="\t", index=False)
    print(f"  per-pool median of mean_pairwise_bc: "
          f"{mean_df['mean_pairwise_bc'].median():.3f}", flush=True)
    print(f"  per (trip_a, trip_b) median bc:")
    p = pair_df.groupby(["trip_a", "trip_b"])["bc"].median().unstack()
    print(p.round(3).to_string())
    return mean_df


# =============================================================================
# 2. Genus-level network on alive only
# =============================================================================
def network_per_pool(ft, label, smeta, tax):
    print(f"\n=== [2] Network ({label}) ===", flush=True)
    od = OUT / f"network_{label}"; od.mkdir(exist_ok=True)
    # Aggregate to genus
    ft2 = ft.copy()
    ft2.index = ft2.index.rename("asv_id")
    m = ft2.reset_index().merge(tax[["asv_id", "genus"]], on="asv_id",
                                  how="left")
    m = m.dropna(subset=["genus"])
    sample_cols = [c for c in m.columns if c not in ("asv_id", "genus")]
    gen = m.groupby("genus")[sample_cols].sum()
    print(f"  genera: {gen.shape[0]}", flush=True)
    summary = []
    for comp in ("rhizosphere", "surface", "deep"):
        samps = list(set(smeta[smeta["compartment"] == comp]["sample"]) &
                       set(gen.columns))
        if len(samps) < 30: continue
        sub = gen[samps]
        # Network construction
        try:
            rho, p = compositional_correlation(sub, min_prevalence=0.10,
                                                 presence_ra=0.0001,
                                                 pseudo=0.5)
            G = build_network(rho, p, rho_threshold=0.4, q_threshold=0.05)
            mods = louvain_modules(G, seed=42)
        except Exception as e:
            print(f"  [{comp}] network failed: {e}", flush=True)
            continue
        n_nodes = G.number_of_nodes()
        n_edges = G.number_of_edges()
        n_iso = sum(1 for n in G.nodes if G.degree(n) == 0)
        # CSP1-2 status
        csp_node = next((n for n in G.nodes if "CSP1" in n or "CSP1-2" in n),
                          None)
        csp_status = {}
        if csp_node:
            deg = G.degree(csp_node)
            mod = mods.get(csp_node)
            ne = list(G.neighbors(csp_node))
            csp_status = {"present": True, "degree": int(deg),
                            "module": int(mod) if mod is not None else None,
                            "n_neighbors": len(ne),
                            "neighbors": ",".join(ne[:10])}
        else:
            csp_status = {"present": False}
        # Save edges + nodes
        if n_edges > 0:
            ed = pd.DataFrame([(u, v, d["rho"], d["weight"], d["sign"], d["q"])
                                  for u, v, d in G.edges(data=True)],
                                 columns=["source", "target", "rho", "weight",
                                            "sign", "q"])
            ed.to_csv(od / f"edges_{comp}.tsv", sep="\t", index=False)
        nd = pd.DataFrame([(n, mods.get(n, -1), G.degree(n))
                              for n in G.nodes],
                             columns=["node", "module", "degree"])
        nd = nd.merge(tax[["genus", "phylum"]].drop_duplicates(),
                        left_on="node", right_on="genus", how="left")
        nd.to_csv(od / f"nodes_{comp}.tsv", sep="\t", index=False)
        print(f"  [{comp}] genera={n_nodes}  edges={n_edges}  "
              f"isolated={n_iso}  CSP1-2: {csp_status}", flush=True)
        summary.append({"compartment": comp, "label": label,
                          "n_genera": n_nodes, "n_edges": n_edges,
                          "n_isolated": n_iso, **csp_status})
    return pd.DataFrame(summary)


# =============================================================================
# 3. iCAMP-style RCbray (lightweight: 99 perms)
# =============================================================================
def rcbray_simplified(ft, label, smeta, n_perm=99):
    print(f"\n=== [3] RCbray-style ({label}, {n_perm} perms) ===", flush=True)
    od = OUT / f"icamp_{label}"; od.mkdir(exist_ok=True)
    # Vectorized RC: for each pair of samples, compute observed BC vs null
    # null: shuffle ASV abundances within each sample (keep richness)
    # Threshold |RC| > 0.95 = stochastic (homogenizing/heterogenizing)
    Mr = relabund(ft)
    rec = []
    for comp in ("rhizosphere", "surface", "deep"):
        samps = list(set(smeta[smeta["compartment"] == comp]["sample"]) &
                       set(Mr.columns))
        if len(samps) < 30: continue
        sub = Mr[samps]
        sub = sub.loc[(sub > 0).any(axis=1)]
        # Sub-sample to 200 samples for computational tractability
        if len(samps) > 200:
            samps = list(RNG.choice(samps, 200, replace=False))
            sub = sub[samps]
        n = len(samps)
        # Observed BC matrix
        BC_obs = squareform(pdist(sub.T.values, metric="braycurtis"))
        # Null: per-sample shuffle of ASV labels within sample (keeps total
        # reads, shuffles which ASVs have which abundance)
        # Simpler null: random reassignment of relative abundances
        n_taxa = sub.shape[0]
        less_count = np.zeros((n, n))
        equal_count = np.zeros((n, n))
        for _ in range(n_perm):
            # Shuffle ASV labels for each sample
            null_M = sub.values.copy()
            for k in range(n):
                idx = RNG.permutation(n_taxa)
                null_M[:, k] = null_M[idx, k]
            BC_null = squareform(pdist(null_M.T, metric="braycurtis"))
            less_count += (BC_null < BC_obs).astype(np.int32)
            equal_count += (BC_null == BC_obs).astype(np.int32)
        RC = 2 * (less_count + 0.5 * equal_count) / n_perm - 1
        # Classify: RC > 0.95 = heterogenizing dispersal
        #           RC < -0.95 = homogenizing dispersal
        #           |RC| <= 0.95 = stochastic / undominated
        iu = np.triu_indices(n, k=1)
        rc_vals = RC[iu]
        n_pairs = len(rc_vals)
        het = float((rc_vals > 0.95).mean())
        hom = float((rc_vals < -0.95).mean())
        und = float((np.abs(rc_vals) <= 0.95).mean())
        rec.append({"compartment": comp, "label": label, "n_samples": n,
                      "n_pairs": n_pairs,
                      "homogenizing_dispersal": hom,
                      "heterogenizing_dispersal": het,
                      "undominated_stochastic": und,
                      "median_BC_obs": float(np.median(BC_obs[iu])),
                      "median_BC_null": float(np.median(squareform(pdist(
                          (sub.values * 0 + sub.values).T, "braycurtis"))[iu]))})
        print(f"  [{comp}] n={n}  hom={hom:.3f}  het={het:.3f}  "
              f"und={und:.3f}", flush=True)
    df = pd.DataFrame(rec)
    df.to_csv(od / "rcbray_summary.tsv", sep="\t", index=False)
    return df


# =============================================================================
# 4. Mediation: salinity -> CSP1-2 -> alive Shannon
# =============================================================================
def mediation_subset(ft, label, smeta, tax, xrf):
    print(f"\n=== [4] Mediation salinity->CSP1-2->Shannon ({label}) ===",
          flush=True)
    od = OUT / "mediation"; od.mkdir(exist_ok=True)
    from sklearn.linear_model import LinearRegression
    # Per-sample alive Shannon
    sh = pd.DataFrame({"sample": ft.columns,
                          "shannon": [shannon(ft[c].values) for c in ft.columns]})
    sh = sh.merge(smeta, on="sample", how="left")
    # CSP1-2 read fraction per sample
    csp_asvs = tax[tax["genus"].astype(str).str.contains("CSP1-2",
                                                              case=False,
                                                              na=False)]
    csp_set = set(csp_asvs["asv_id"])
    csp_in_ft = ft.index.isin(csp_set)
    sample_total = ft.sum(axis=0)
    csp_reads = ft.loc[csp_in_ft].sum(axis=0)
    sh["csp_frac"] = (csp_reads / sample_total.replace(0, 1)).values
    # XRF: per sample
    if "SampleID" in xrf.columns:
        xrf_sub = xrf.copy()
        xrf_sub["sample"] = xrf_sub["SampleID"]
    else:
        xrf_sub = xrf.copy()
    if "Na" not in xrf_sub.columns:
        print("  XRF Na column missing, skipping"); return None
    salinity_cols = [c for c in ("Na", "Cl", "SO3") if c in xrf_sub.columns]
    if not salinity_cols:
        print("  no salinity cols"); return None
    xrf_per_sample = (xrf_sub.groupby("sample")[salinity_cols].mean()
                       .reset_index())
    xrf_per_sample["salinity"] = xrf_per_sample[salinity_cols].mean(axis=1)
    sh = sh.merge(xrf_per_sample[["sample", "salinity"]], on="sample",
                    how="inner")
    # OK now we have shannon, csp_frac, salinity
    print(f"  N samples with all 3: {len(sh)}", flush=True)
    if len(sh) < 30:
        print("  insufficient"); return None
    # Standardize
    for c in ("shannon", "csp_frac", "salinity"):
        sh[c] = (sh[c] - sh[c].mean()) / sh[c].std()
    # Total effect: shannon ~ salinity
    lr1 = LinearRegression().fit(sh[["salinity"]], sh["shannon"])
    total_beta = float(lr1.coef_[0])
    # a path: csp_frac ~ salinity
    lr2 = LinearRegression().fit(sh[["salinity"]], sh["csp_frac"])
    a = float(lr2.coef_[0])
    # b path: shannon ~ csp_frac + salinity
    lr3 = LinearRegression().fit(sh[["csp_frac", "salinity"]], sh["shannon"])
    b = float(lr3.coef_[0])
    direct = float(lr3.coef_[1])
    indirect = a * b
    indirect_share = indirect / total_beta if total_beta != 0 else np.nan
    rec = {"label": label, "n": len(sh),
            "total_beta": total_beta,
            "a_path_sal_to_csp": a, "b_path_csp_to_shannon": b,
            "direct_beta": direct, "indirect_beta": indirect,
            "indirect_share_of_total": indirect_share}
    print(f"  total beta(sal->shannon): {total_beta:+.3f}")
    print(f"  a (sal->csp): {a:+.3f}")
    print(f"  b (csp->shannon|sal): {b:+.3f}")
    print(f"  direct (sal->shannon|csp): {direct:+.3f}")
    print(f"  indirect (a*b): {indirect:+.3f}")
    print(f"  indirect / total: {indirect_share:+.3f}")
    return rec


def main():
    print("Loading inputs ...", flush=True)
    ft_all = pd.read_parquet(CACHE / "feature_table.parquet")
    ft_alive = pd.read_parquet(CACHE / "feature_table_alive.parquet")
    ft_relic = pd.read_parquet(CACHE / "feature_table_relic.parquet")
    smeta = parse_samples_to_df(ft_all.columns)
    smeta["site"] = smeta["site"].astype(int)
    tax = pd.read_parquet(CACHE / "taxonomy.parquet").reset_index().rename(
        columns={"ASV": "asv_id"})
    xrf = pd.read_csv(DATA / "geochemistry" / "xrf_lab_table_all_trips.tsv",
                       sep="\t")

    summary = []
    for label, ft in (("all", ft_all), ("alive", ft_alive),
                          ("relic", ft_relic)):
        print(f"\n{'#'*70}\n# POOL = {label.upper()}\n{'#'*70}", flush=True)
        try:
            t = all_pairwise_temporal(ft, label, smeta)
            summary.append({"analysis": "temporal_pairwise", "label": label,
                              "median_mean_bc": float(t["mean_pairwise_bc"]
                                                          .median())})
        except Exception as e:
            print(f"  temporal failed: {e}")
        try:
            n = network_per_pool(ft, label, smeta, tax)
            for _, r in n.iterrows():
                summary.append({"analysis": "network", "label": label,
                                  "compartment": r["compartment"],
                                  "n_genera": r["n_genera"],
                                  "n_edges": r["n_edges"],
                                  "csp_present": r.get("present", False),
                                  "csp_degree": r.get("degree", None)})
        except Exception as e:
            print(f"  network failed: {e}")
        try:
            ic = rcbray_simplified(ft, label, smeta, n_perm=99)
            for _, r in ic.iterrows():
                summary.append({"analysis": "rcbray", "label": label,
                                  "compartment": r["compartment"],
                                  "homogenizing": r["homogenizing_dispersal"],
                                  "heterogenizing": r["heterogenizing_dispersal"],
                                  "undominated": r["undominated_stochastic"]})
        except Exception as e:
            print(f"  rcbray failed: {e}")
        try:
            med = mediation_subset(ft, label, smeta, tax, xrf)
            if med:
                summary.append({"analysis": "mediation", **med})
        except Exception as e:
            print(f"  mediation failed: {e}")

    sm = pd.DataFrame(summary)
    sm.to_csv(OUT / "remaining_analyses_summary.tsv", sep="\t", index=False)
    print(f"\nWrote summary to {OUT}/remaining_analyses_summary.tsv")


if __name__ == "__main__":
    main()
