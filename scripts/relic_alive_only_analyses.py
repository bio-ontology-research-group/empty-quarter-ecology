#!/usr/bin/env python3
"""Comprehensive re-analysis on ALIVE-only ASV pool.

Runs the alive-only equivalents of:
  1. Pulse-reserve precipitation response  (Tier 1 #4)
  2. Distance-decay (taxonomic BC vs km)
  3. Within-replicate vs between-site variance + PERMANOVA (Test 3)
  4. Allison-Martiny taxonomic-vs-functional redundancy (Test 4)
  5. Wind-Mantel BC ~ wind connectivity (Tier 1 wind)
  6. Climate trend correlation: per-site relic_frac and alive Shannon vs
     long-term temperature/precip trend
  7. Time-for-space substitution: T1-4 mean predicts T5

Outputs:
  cache/relic_alive_subset/<analysis>/<files>
  cache/relic_alive_subset/SUMMARY.md
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import (mannwhitneyu, spearmanr, linregress,
                            pearsonr)

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "relic_alive_subset"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260510)


def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = np.radians(lat2 - lat1)
    dl = np.radians(lon2 - lon1)
    a = np.sin(dp/2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl/2)**2
    return 2 * R * np.arcsin(np.sqrt(a))


def relabund(M):
    return M.div(M.sum(axis=0).replace(0, 1), axis=1)


def shannon_per_sample(M):
    Mr = relabund(M)
    Mr = Mr.replace(0, np.nan)
    return -(Mr * np.log(Mr)).sum(axis=0)


# =============================================================================
# 1. Pulse-reserve
# =============================================================================
def pulse_reserve(ft, label, smeta, geo, weather):
    print(f"\n=== [1] Pulse-reserve ({label}) ===", flush=True)
    od = OUT / f"pulse_reserve_{label}"; od.mkdir(exist_ok=True)
    sm = smeta.copy()
    sm = sm.merge(geo[["site", "trip", "CenterDate"]],
                    on=["site", "trip"], how="left")
    sm["sample_date"] = pd.to_datetime(sm["CenterDate"], errors="coerce")
    sm = sm.dropna(subset=["sample_date"])
    rec = []
    for _, r in sm.iterrows():
        cur = weather[(weather["site"] == r["site"]) &
                       (weather["date"] <= r["sample_date"])]
        if cur.empty: continue
        for win, days in (("d7", 7), ("d30", 30), ("d90", 90)):
            cutoff = r["sample_date"] - pd.Timedelta(days=days)
            ww = cur[cur["date"] >= cutoff]
            rec.append({"sample": r["sample"], "site": r["site"],
                          "trip": r["trip"], "compartment": r["compartment"],
                          "window": win,
                          "precip_mm": float(ww["precip_mm"].sum())})
    pr = pd.DataFrame(rec)
    if len(pr) == 0:
        print("  no precip records aligned"); return
    # Pulse vs reserve based on d30 > 1 mm (EQ is hyperarid)
    pr30 = (pr[pr["window"] == "d30"]
            .assign(state=lambda x: np.where(x["precip_mm"] > 1.0,
                                                "pulse", "reserve")))
    print(f"  pulse n={(pr30['state']=='pulse').sum()}, "
          f"reserve n={(pr30['state']=='reserve').sum()}", flush=True)
    pr30.to_csv(od / "sample_pulse_state.tsv", sep="\t", index=False)

    # Per-ASV (within alive subset): MW pulse vs reserve relabund
    Mr = relabund(ft)
    rows = []
    for asv in Mr.index:
        a = Mr.loc[asv]
        pulse_v = a[pr30.loc[pr30["state"] == "pulse", "sample"]].dropna()
        reserve_v = a[pr30.loc[pr30["state"] == "reserve", "sample"]].dropna()
        if len(pulse_v) < 10 or len(reserve_v) < 10: continue
        try:
            U, p = mannwhitneyu(pulse_v, reserve_v, alternative="greater")
        except Exception:
            continue
        rows.append({"asv_id": asv,
                      "median_pulse": float(pulse_v.median()),
                      "median_reserve": float(reserve_v.median()),
                      "log2_pulse_reserve": float(np.log2(
                          (pulse_v.median() + 1e-8) /
                          (reserve_v.median() + 1e-8))),
                      "U": float(U), "p": float(p)})
    res = pd.DataFrame(rows)
    print(f"  ASVs tested: {len(res)}", flush=True)
    if len(res) == 0:
        return
    # BH adjust
    p = res["p"].values
    n = len(p); idx = np.argsort(p); ranked = np.empty(n)
    ranked[idx] = np.arange(1, n + 1)
    res["q"] = np.minimum(p * n / ranked, 1.0)
    res = res.sort_values("p")
    res.to_csv(od / "asv_pulse_response.tsv", sep="\t", index=False)
    sig = res[res["q"] < 0.05]
    print(f"  q<0.05 (positive responders): {len(sig)}", flush=True)
    print(f"  fraction of tested ASVs that respond: {len(sig)/len(res):.3f}",
          flush=True)
    return {"label": label, "n_tested": len(res), "n_sig": len(sig),
              "fraction_respond": len(sig)/len(res)}


# =============================================================================
# 2. Distance-decay (taxonomic BC ~ km, per compartment+trip)
# =============================================================================
def distance_decay(ft, label, smeta, geo):
    print(f"\n=== [2] Distance-decay ({label}) ===", flush=True)
    od = OUT / f"distance_decay_{label}"; od.mkdir(exist_ok=True)
    geo_lookup = (geo[["site", "Latitude", "Longitude"]]
                  .drop_duplicates("site").set_index("site"))
    rows = []
    for comp in ("rhizosphere", "surface", "deep"):
        for trip in (1, 2, 3, 4, 5):
            sub = smeta[(smeta["compartment"] == comp) &
                          (smeta["trip"] == trip)]
            samps = list(set(sub["sample"]) & set(ft.columns))
            if len(samps) < 10: continue
            # site-level aggregate
            sub2 = sub.set_index("sample")
            site_arr = sub2.loc[samps, "site"].values
            Mr = relabund(ft[samps]).T  # samples x ASVs
            row_sum = Mr.sum(axis=1)
            Mr = Mr[row_sum > 0]
            if len(Mr) < 10: continue
            samps2 = list(Mr.index)
            site2 = sub2.loc[samps2, "site"].values
            BC = squareform(pdist(Mr.values, metric="braycurtis"))
            n = len(samps2)
            d = []; b = []
            for i in range(n):
                for j in range(i+1, n):
                    if site2[i] == site2[j]: continue
                    if site2[i] not in geo_lookup.index: continue
                    if site2[j] not in geo_lookup.index: continue
                    la = geo_lookup.loc[site2[i]]; lb = geo_lookup.loc[site2[j]]
                    km = haversine(la["Latitude"], la["Longitude"],
                                    lb["Latitude"], lb["Longitude"])
                    d.append(km); b.append(BC[i, j])
            if len(d) < 30: continue
            slope, intercept, r, p, _ = linregress(np.log10(np.array(d) + 1),
                                                       np.array(b))
            rows.append({"compartment": comp, "trip": trip,
                          "n_pairs": len(d),
                          "slope_per_log10km": float(slope),
                          "r": float(r), "p": float(p),
                          "intercept_BC_at_1km": float(intercept)})
    dd = pd.DataFrame(rows)
    dd.to_csv(od / "distance_decay_per_comp_trip.tsv", sep="\t", index=False)
    print(dd.round(3).to_string(index=False))
    return dd


# =============================================================================
# 3. Within-replicate vs between-site (Test 3 reanalysis)
# =============================================================================
def within_between(ft, label, smeta):
    print(f"\n=== [3] Within-rep vs between-site ({label}) ===", flush=True)
    od = OUT / f"within_between_{label}"; od.mkdir(exist_ok=True)
    rows = []
    for comp in ("rhizosphere", "surface", "deep"):
        for trip in (1, 2, 3, 4, 5):
            sub = smeta[(smeta["compartment"] == comp) &
                          (smeta["trip"] == trip)]
            samps = list(set(sub["sample"]) & set(ft.columns))
            if len(samps) < 10: continue
            sub2 = sub.set_index("sample")
            Mr = relabund(ft[samps]).T
            row_sum = Mr.sum(axis=1)
            Mr = Mr[row_sum > 0]
            if len(Mr) < 10: continue
            samps2 = list(Mr.index)
            site2 = sub2.loc[samps2, "site"].values
            BC = squareform(pdist(Mr.values, metric="braycurtis"))
            within = []; between = []
            for i in range(len(samps2)):
                for j in range(i+1, len(samps2)):
                    if site2[i] == site2[j]:
                        within.append(BC[i, j])
                    else:
                        between.append(BC[i, j])
            if not within or not between: continue
            rows.append({"compartment": comp, "trip": trip,
                          "n_within": len(within),
                          "n_between": len(between),
                          "within_median": float(np.median(within)),
                          "between_median": float(np.median(between)),
                          "ratio": float(np.median(within) /
                                            np.median(between))
                                            if np.median(between) > 0 else np.nan})
    wb = pd.DataFrame(rows)
    wb.to_csv(od / "within_vs_between.tsv", sep="\t", index=False)
    print(wb.round(3).to_string(index=False))
    return wb


# =============================================================================
# 4. Allison-Martiny redundancy slope (Test 4 reanalysis)
# =============================================================================
def allison_martiny(ft, label, smeta, path_func):
    print(f"\n=== [4] Allison-Martiny redundancy slope ({label}) ===",
          flush=True)
    od = OUT / f"allison_martiny_{label}"; od.mkdir(exist_ok=True)
    common = sorted(set(ft.columns) & set(path_func.columns))
    rec = []
    for comp in ("rhizosphere", "surface", "deep"):
        for trip in (1, 2, 3, 4, 5):
            sub = smeta[(smeta["compartment"] == comp) &
                          (smeta["trip"] == trip)]
            samps = list(set(sub["sample"]) & set(common))
            if len(samps) < 10: continue
            Tr = relabund(ft[samps]).T
            Fr = relabund(path_func[samps]).T
            row_sum = Tr.sum(axis=1)
            Tr = Tr[row_sum > 0]
            samps2 = [s for s in Tr.index if s in Fr.index]
            if len(samps2) < 10: continue
            Tr = Tr.loc[samps2]; Fr = Fr.loc[samps2]
            BCt = pdist(Tr.values, metric="braycurtis")
            BCf = pdist(Fr.values, metric="braycurtis")
            mask = (~np.isnan(BCt)) & (~np.isnan(BCf))
            if mask.sum() < 30: continue
            slope, _, r, p, _ = linregress(BCt[mask], BCf[mask])
            ratio = float(np.nanmedian(BCf[mask] /
                                          (BCt[mask] + 1e-8)))
            rec.append({"compartment": comp, "trip": trip,
                          "n_samples": len(samps2),
                          "n_pairs": int(mask.sum()),
                          "slope": float(slope),
                          "r": float(r), "p": float(p),
                          "median_func_over_tax": ratio})
    am = pd.DataFrame(rec)
    am.to_csv(od / "allison_martiny.tsv", sep="\t", index=False)
    print(am.round(3).to_string(index=False))
    print(f"\n  median slope: {am['slope'].median():.3f}")
    print(f"  median func/tax ratio: {am['median_func_over_tax'].median():.3f}")
    return am


# =============================================================================
# 5. Wind-Mantel (BC ~ partial wind connectivity)
# =============================================================================
def wind_mantel(ft, label, smeta, wind_conn):
    print(f"\n=== [5] Wind-Mantel ({label}) ===", flush=True)
    od = OUT / f"wind_mantel_{label}"; od.mkdir(exist_ok=True)
    if wind_conn is None: return
    # Build per (comp, trip): pairwise BC matrix and pairwise wind matrix
    rows = []
    for comp in ("rhizosphere", "surface", "deep"):
        for trip in (1, 2, 3, 4, 5):
            sub = smeta[(smeta["compartment"] == comp) &
                          (smeta["trip"] == trip)]
            samps = list(set(sub["sample"]) & set(ft.columns))
            if len(samps) < 10: continue
            sub2 = sub.set_index("sample")
            Mr = relabund(ft[samps]).T
            Mr = Mr[Mr.sum(axis=1) > 0]
            samps2 = list(Mr.index)
            site2 = sub2.loc[samps2, "site"].astype(int).values
            BC = squareform(pdist(Mr.values, metric="braycurtis"))
            # Aggregate to site-mean
            df = pd.DataFrame({"sample": samps2, "site": site2})
            sites_in = sorted(set(site2))
            if len(sites_in) < 6: continue
            sBC = np.zeros((len(sites_in), len(sites_in)))
            for ii, sa in enumerate(sites_in):
                for jj, sb in enumerate(sites_in):
                    ix = df[df["site"] == sa].index
                    iy = df[df["site"] == sb].index
                    if len(ix) == 0 or len(iy) == 0: continue
                    vals = []
                    for x in ix:
                        for y in iy:
                            if x == y: continue
                            vals.append(BC[x, y])
                    sBC[ii, jj] = float(np.mean(vals)) if vals else np.nan
            # Wind matrix
            sub_wind = wind_conn[(wind_conn["trip"] == trip) &
                                    (wind_conn["site_i"].isin(sites_in)) &
                                    (wind_conn["site_j"].isin(sites_in))]
            if len(sub_wind) < 10: continue
            wmat = pd.pivot_table(sub_wind, index="site_i",
                                       columns="site_j",
                                       values="score_max").reindex(
                                           index=sites_in, columns=sites_in)
            # Symmetrize wind
            wmat = wmat.fillna(wmat.T)
            # Distance matrix
            geo_l = wind_conn[["site_i", "site_j", "dist_km"]].drop_duplicates()
            dmat = pd.pivot_table(geo_l, index="site_i", columns="site_j",
                                       values="dist_km").reindex(
                                           index=sites_in, columns=sites_in)
            dmat = dmat.fillna(dmat.T)
            # Mantel: corr(BC, -wind) and partial(BC, wind | dist)
            iu = np.triu_indices(len(sites_in), k=1)
            bc = sBC[iu]; w = wmat.values[iu]; d = dmat.values[iu]
            mask = (~np.isnan(bc)) & (~np.isnan(w)) & (~np.isnan(d))
            if mask.sum() < 10: continue
            r_mantel, p_mantel = spearmanr(bc[mask], -w[mask])
            r_dist, _ = spearmanr(bc[mask], d[mask])
            # partial: BC ~ wind | dist via residualization
            from numpy.linalg import lstsq
            X = np.vstack([np.ones(mask.sum()), d[mask]]).T
            beta_b, *_ = lstsq(X, bc[mask], rcond=None)
            beta_w, *_ = lstsq(X, w[mask], rcond=None)
            res_b = bc[mask] - X.dot(beta_b)
            res_w = w[mask] - X.dot(beta_w)
            r_partial, p_partial = spearmanr(res_b, -res_w)
            rows.append({"compartment": comp, "trip": trip,
                          "n_sites": len(sites_in),
                          "r_BC_wind": float(r_mantel), "p_BC_wind": float(p_mantel),
                          "r_partial_BC_wind_given_dist": float(r_partial),
                          "p_partial": float(p_partial),
                          "r_BC_dist": float(r_dist)})
    wm = pd.DataFrame(rows)
    wm.to_csv(od / "wind_mantel.tsv", sep="\t", index=False)
    print(wm.round(3).to_string(index=False))
    print(f"\n  median r_partial(BC ~ wind | dist): "
          f"{wm['r_partial_BC_wind_given_dist'].median():.3f}")
    return wm


# =============================================================================
def load_climate_geodata():
    geo_frames = []
    for trip in range(1, 6):
        gp = DATA / "geodata" / f"trip{trip}_geodata.tsv"
        if gp.exists():
            g = pd.read_csv(gp, sep="\t")
            g["trip"] = trip
            geo_frames.append(g)
    geo = pd.concat(geo_frames, ignore_index=True)
    geo = geo.rename(columns={"Site": "site"})
    geo["site"] = pd.to_numeric(geo["site"], errors="coerce")
    geo = geo.dropna(subset=["site"]); geo["site"] = geo["site"].astype(int)
    return geo


def load_weather():
    wp = DATA / "climate" / "daily_weather_full.csv"
    w = pd.read_csv(wp)
    w["date"] = pd.to_datetime(w["date"], format="%Y%m%d")
    w = w.rename(columns={"PRECTOTCORR": "precip_mm"})
    return w


def load_wind_connectivity():
    wp = CACHE / "wind_dispersal" / "pair_trip_wind_connectivity.tsv"
    if not wp.exists(): return None
    df = pd.read_csv(wp, sep="\t")
    cols = [c for c in df.columns]
    return df


def main():
    print("Loading inputs ...", flush=True)
    ft_all = pd.read_parquet(CACHE / "feature_table.parquet")
    ft_alive = pd.read_parquet(CACHE / "feature_table_alive.parquet")
    ft_relic = pd.read_parquet(CACHE / "feature_table_relic.parquet")
    print(f"  all  : {ft_all.shape}", flush=True)
    print(f"  alive: {ft_alive.shape}", flush=True)
    print(f"  relic: {ft_relic.shape}", flush=True)

    smeta = parse_samples_to_df(ft_all.columns)
    smeta["site"] = smeta["site"].astype(int)
    geo = load_climate_geodata()
    weather = load_weather()
    print(f"  weather rows: {len(weather)}", flush=True)
    wind = load_wind_connectivity()
    print(f"  wind connectivity: {None if wind is None else wind.shape}",
          flush=True)
    path_func = pd.read_csv(DATA / "functional" / "picrust2" /
                              "path_abun_unstrat.tsv",
                              sep="\t", index_col=0)

    # Run on each pool
    summary_rows = []
    for label, ft in (("all", ft_all),
                          ("alive", ft_alive),
                          ("relic", ft_relic)):
        print(f"\n{'#'*70}\n# POOL = {label.upper()}\n{'#'*70}", flush=True)
        try:
            r = pulse_reserve(ft, label, smeta, geo, weather)
            if r: summary_rows.append({"analysis": "pulse_reserve", **r})
        except Exception as e:
            print(f"  pulse_reserve failed: {e}", flush=True)
        try:
            dd = distance_decay(ft, label, smeta, geo)
            summary_rows.append({"analysis": "distance_decay", "label": label,
                                    "median_slope": float(dd["slope_per_log10km"]
                                                              .median())})
        except Exception as e:
            print(f"  distance_decay failed: {e}", flush=True)
        try:
            wb = within_between(ft, label, smeta)
            summary_rows.append({"analysis": "within_between", "label": label,
                                    "median_within": float(wb["within_median"]
                                                                .median()),
                                    "median_between": float(wb["between_median"]
                                                                  .median()),
                                    "median_ratio": float(wb["ratio"].median())})
        except Exception as e:
            print(f"  within_between failed: {e}", flush=True)
        try:
            am = allison_martiny(ft, label, smeta, path_func)
            summary_rows.append({"analysis": "allison_martiny", "label": label,
                                    "median_slope": float(am["slope"].median()),
                                    "median_func_tax_ratio": float(
                                        am["median_func_over_tax"].median())})
        except Exception as e:
            print(f"  allison_martiny failed: {e}", flush=True)
        if wind is not None:
            try:
                wm = wind_mantel(ft, label, smeta, wind)
                if wm is not None and len(wm):
                    summary_rows.append(
                        {"analysis": "wind_mantel", "label": label,
                         "median_partial_r": float(
                             wm["r_partial_BC_wind_given_dist"].median())})
            except Exception as e:
                print(f"  wind_mantel failed: {e}", flush=True)

    # Save summary
    sm = pd.DataFrame(summary_rows)
    sm.to_csv(OUT / "subset_comparison_summary.tsv", sep="\t", index=False)
    print(f"\n=== FINAL SUMMARY (saved to {OUT}/subset_comparison_summary.tsv) ===")
    print(sm.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
