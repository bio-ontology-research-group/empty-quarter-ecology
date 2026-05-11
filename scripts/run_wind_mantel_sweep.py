#!/usr/bin/env python3
"""Comprehensive wind-Mantel sweep:
   sensitivity (4 angles × 3 dust thresholds)
 × score type  (4: max, min, sum, asymmetry)
 × distance stratum (4: all, <100, 100-500, >500 km)
 × time window (5: 7, 14, 30, 90, 365 days)
 × per (compartment × trip).

Strategy: build per (window × angle × threshold) connectivity matrices
ONCE (with both score_ij and score_ji per pair-trip). Then per Mantel
test, derive the score type as needed and apply distance stratum mask.

Inputs: same as run_wind_dispersal_mantel.py.

Outputs:
  cache/wind_dispersal/sweep_mantel_full.tsv      every test
  cache/wind_dispersal/sweep_summary.txt          digestible summary
"""
from __future__ import annotations

import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "wind_dispersal"
OUT.mkdir(parents=True, exist_ok=True)

WINDOWS_DAYS = [7, 14, 30, 90, 365]
ANGLE_TOL_DEG = [15, 30, 45, 60]
DUST_THR_MS = [5.0, 7.0, 10.0]
SCORE_TYPES = ["max", "min", "sum", "asymmetry"]   # all symmetric by construction
DIST_STRATA = [("all", 0, 1e9),
               ("lt100km", 0, 100),
               ("100_500km", 100, 500),
               ("gt500km", 500, 1e9)]
N_PERM = 199          # comprehensive sweep — moderate precision
RNG = np.random.default_rng(20260509)


def load_geo() -> pd.DataFrame:
    geo = []
    for trip in range(1, 6):
        gp = DATA / "geodata" / f"trip{trip}_geodata.tsv"
        if gp.exists():
            g = pd.read_csv(gp, sep="\t")
            g["trip"] = trip
            g["CenterDate"] = pd.to_datetime(g["CenterDate"])
            g["Site_int"] = pd.to_numeric(g["Site"], errors="coerce")
            g = g.dropna(subset=["Site_int", "CenterDate"])
            g["site"] = g["Site_int"].astype(int)
            geo.append(g[["site", "trip", "CenterDate"]])
    return pd.concat(geo).drop_duplicates(["site", "trip"])


def angular_diff(a, b):
    return np.abs((a - b + 180) % 360 - 180)


def site_trip_communities(ft: pd.DataFrame, smeta: pd.DataFrame,
                            compartment: str) -> dict:
    sub = smeta.loc[smeta["compartment"] == compartment]
    out = {}
    rel = ft.div(ft.sum(axis=0).replace(0, 1), axis=1)
    for (s, t), grp in sub.groupby(["site", "trip"]):
        cols = grp["sample"].tolist()
        if len(cols) == 0:
            continue
        v = rel[cols].mean(axis=1).values
        out[(int(s), int(t))] = v
    return out


def pairwise_bc(comm: dict, sites: list, trip: int) -> np.ndarray:
    n = len(sites)
    M = np.full((n, n), np.nan)
    for i in range(n):
        ci = comm.get((sites[i], trip))
        if ci is None: continue
        for j in range(i + 1, n):
            cj = comm.get((sites[j], trip))
            if cj is None: continue
            num = np.abs(ci - cj).sum()
            den = ci.sum() + cj.sum()
            d = num / den if den > 0 else np.nan
            M[i, j] = M[j, i] = d
    np.fill_diagonal(M, 0.0)
    return M


def mantel(D1, D2, mask, n_perm=N_PERM):
    """Mantel test on upper-triangle entries selected by mask."""
    n = D1.shape[0]
    iu = np.triu_indices(n, k=1)
    x = D1[iu]; y = D2[iu]; m = mask[iu]
    valid = np.isfinite(x) & np.isfinite(y) & m
    if valid.sum() < 10:
        return (np.nan, np.nan, int(valid.sum()))
    x = x[valid]; y = y[valid]
    r_obs = np.corrcoef(x, y)[0, 1]
    cnt = 0
    for _ in range(n_perm):
        perm = RNG.permutation(n)
        D2p = D2[perm][:, perm]
        yp = D2p[iu][valid]
        if abs(np.corrcoef(x, yp)[0, 1]) >= abs(r_obs):
            cnt += 1
    p = (cnt + 1) / (n_perm + 1)
    return (float(r_obs), float(p), int(valid.sum()))


def partial_mantel(D_y, D_x, D_z, mask, n_perm=N_PERM):
    n = D_y.shape[0]
    iu = np.triu_indices(n, k=1)
    y = D_y[iu]; x = D_x[iu]; z = D_z[iu]; m = mask[iu]
    valid = np.isfinite(y) & np.isfinite(x) & np.isfinite(z) & m
    if valid.sum() < 10:
        return (np.nan, np.nan, int(valid.sum()))
    y = y[valid]; x = x[valid]; z = z[valid]

    def pr(yv, xv, zv):
        ryx = np.corrcoef(yv, xv)[0, 1]
        ryz = np.corrcoef(yv, zv)[0, 1]
        rxz = np.corrcoef(xv, zv)[0, 1]
        denom = np.sqrt((1 - ryz**2) * (1 - rxz**2))
        return (ryx - ryz * rxz) / denom if denom > 0 else np.nan

    r_obs = pr(y, x, z)
    cnt = 0
    for _ in range(n_perm):
        perm = RNG.permutation(n)
        Dxp = D_x[perm][:, perm]
        xp = Dxp[iu][valid]
        rr = pr(y, xp, z)
        if rr == rr and abs(rr) >= abs(r_obs):
            cnt += 1
    p = (cnt + 1) / (n_perm + 1)
    return (float(r_obs), float(p), int(valid.sum()))


def build_pair_scores(weather: pd.DataFrame, geo: pd.DataFrame,
                       geom: pd.DataFrame, window_days: int,
                       angle_tol: float, dust_thr: float) -> pd.DataFrame:
    """For each (pair, trip) compute score_ij and score_ji with given params."""
    w = weather.set_index(["site", "date"]).sort_index()
    win_cache = {}
    for _, r in geo.iterrows():
        site, trip = int(r["site"]), int(r["trip"])
        d = pd.Timestamp(r["CenterDate"])
        try:
            site_w = w.loc[site]
        except KeyError:
            continue
        mask = (site_w.index >= d - pd.Timedelta(days=window_days - 1)) & \
               (site_w.index <= d)
        sub = site_w.loc[mask, ["WD10M", "WS10M_MAX"]].copy()
        sub["dir_to"] = (sub["WD10M"].values + 180) % 360
        win_cache[(site, trip)] = sub

    rows = []
    for _, r in geom.iterrows():
        si, sj = int(r["site_i"]), int(r["site_j"])
        bij, bji = float(r["bearing_ij_deg"]), float(r["bearing_ji_deg"])
        for trip in range(1, 6):
            wi = win_cache.get((si, trip))
            wj = win_cache.get((sj, trip))
            if wi is None or wj is None or len(wi) == 0 or len(wj) == 0:
                continue
            # i -> j
            dust_i = wi["WS10M_MAX"].values > dust_thr
            align_i = angular_diff(wi["dir_to"].values, bij) <= angle_tol
            score_ij = float(np.sum(wi["WS10M_MAX"].values[align_i & dust_i]))
            # j -> i
            dust_j = wj["WS10M_MAX"].values > dust_thr
            align_j = angular_diff(wj["dir_to"].values, bji) <= angle_tol
            score_ji = float(np.sum(wj["WS10M_MAX"].values[align_j & dust_j]))
            rows.append({"site_i": si, "site_j": sj, "trip": trip,
                         "score_ij": score_ij, "score_ji": score_ji,
                         "dist_km": float(r["dist_km"])})
    return pd.DataFrame(rows)


def main():
    print(f"[load] feature_table, geometry, weather, geodata", flush=True)
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    geom = pd.read_csv(CACHE / "pairwise_geometry.tsv", sep="\t")
    weather = pd.read_csv(DATA / "climate" / "daily_weather_full.csv")
    weather["date"] = pd.to_datetime(weather["date"], format="%Y%m%d")
    weather["site"] = weather["site"].astype(int)
    geo = load_geo()
    smeta = parse_samples_to_df(ft.columns)
    sites = sorted(set(int(s) for s in geo["site"].unique()))
    n_s = len(sites)
    site_idx = {s: i for i, s in enumerate(sites)}

    # Distance matrix (constant across all tests)
    D_dist = np.zeros((n_s, n_s))
    for _, r in geom.iterrows():
        i = site_idx[int(r["site_i"])]
        j = site_idx[int(r["site_j"])]
        D_dist[i, j] = D_dist[j, i] = float(r["dist_km"])

    # BC matrices per (compartment, trip), computed once
    print(f"[BC] computing site-aggregated BC per (comp, trip)...", flush=True)
    bc_matrices = {}
    for comp in ["rhizosphere", "surface", "deep"]:
        comm = site_trip_communities(ft, smeta, comp)
        for trip in range(1, 6):
            bc_matrices[(comp, trip)] = pairwise_bc(comm, sites, trip)

    # Distance stratum masks
    dist_masks = {}
    for label, lo, hi in DIST_STRATA:
        m = (D_dist > lo) & (D_dist <= hi)
        dist_masks[label] = m

    # Iterate the parameter grid
    print(f"[sweep] iterating {len(WINDOWS_DAYS)}*{len(ANGLE_TOL_DEG)}"
          f"*{len(DUST_THR_MS)} = "
          f"{len(WINDOWS_DAYS)*len(ANGLE_TOL_DEG)*len(DUST_THR_MS)} parameter combos",
          flush=True)
    rows = []
    t0 = time.time()
    n_combos_done = 0
    for w_days in WINDOWS_DAYS:
        for ang in ANGLE_TOL_DEG:
            for thr in DUST_THR_MS:
                pair = build_pair_scores(weather, geo, geom, w_days, ang, thr)
                # per trip, build wind-connectivity matrix for each score type
                for trip in range(1, 6):
                    sub = pair[pair["trip"] == trip]
                    if len(sub) == 0:
                        continue
                    # build base score matrices
                    M_ij = np.zeros((n_s, n_s))
                    M_ji = np.zeros((n_s, n_s))
                    for _, r in sub.iterrows():
                        i = site_idx[int(r["site_i"])]
                        j = site_idx[int(r["site_j"])]
                        M_ij[i, j] = float(r["score_ij"])
                        M_ji[i, j] = float(r["score_ji"])
                        # symmetric reflection
                        M_ij[j, i] = float(r["score_ji"])
                        M_ji[j, i] = float(r["score_ij"])
                    # 4 score types
                    score_mats = {
                        "max":       np.maximum(M_ij, M_ji),
                        "min":       np.minimum(M_ij, M_ji),
                        "sum":       M_ij + M_ji,
                        "asymmetry": np.abs(M_ij - M_ji),
                    }
                    # convert to "dissimilarity" by negation (smaller = more connected
                    # for max/min/sum; smaller |diff| = more symmetric for asymmetry)
                    dissim_mats = {
                        "max":       -score_mats["max"],
                        "min":       -score_mats["min"],
                        "sum":       -score_mats["sum"],
                        "asymmetry":  score_mats["asymmetry"],
                    }
                    for comp in ["rhizosphere", "surface", "deep"]:
                        bc = bc_matrices[(comp, trip)]
                        iu = np.triu_indices(n_s, k=1)
                        n_valid_bc = np.isfinite(bc[iu]).sum()
                        if n_valid_bc < 30:
                            continue
                        for stratum_name, mask in dist_masks.items():
                            for stype, D_w in dissim_mats.items():
                                np.fill_diagonal(D_w, 0.0)
                                r_dist, p_dist, n_d = mantel(bc, D_dist, mask)
                                r_wind, p_wind, n_w = mantel(bc, D_w, mask)
                                r_part, p_part, n_p = partial_mantel(bc, D_w, D_dist, mask)
                                rows.append({
                                    "window_days": w_days,
                                    "angle_tol": ang,
                                    "dust_thr": thr,
                                    "score_type": stype,
                                    "stratum": stratum_name,
                                    "compartment": comp,
                                    "trip": trip,
                                    "n_pairs": int(n_p),
                                    "r_BC_dist": r_dist, "p_BC_dist": p_dist,
                                    "r_BC_wind": r_wind, "p_BC_wind": p_wind,
                                    "r_part": r_part, "p_part": p_part,
                                })
                n_combos_done += 1
                el = time.time() - t0
                rate = n_combos_done / el if el > 0 else 0
                eta = (len(WINDOWS_DAYS)*len(ANGLE_TOL_DEG)*len(DUST_THR_MS)
                        - n_combos_done) / max(rate, 0.001)
                print(f"  [{n_combos_done:3d}] w={w_days:>3d}d ang={ang:>2d}° "
                      f"thr={thr:>4.1f} m/s  cumul rows={len(rows)}  "
                      f"elapsed={el:.0f}s ETA={eta:.0f}s",
                      flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "sweep_mantel_full.tsv", sep="\t", index=False)
    print(f"\nWrote {OUT}/sweep_mantel_full.tsv: {len(df)} rows", flush=True)

    # Build digestible summaries
    with open(OUT / "sweep_summary.txt", "w") as fh:
        fh.write("Comprehensive wind-Mantel sweep (Tier-1 followup)\n")
        fh.write("=" * 72 + "\n\n")
        fh.write(f"Total Mantel tests: {len(df)}\n")
        fh.write(f"Permutations per test: {N_PERM}\n")
        fh.write(f"Windows (days):      {WINDOWS_DAYS}\n")
        fh.write(f"Angle tolerance (°): {ANGLE_TOL_DEG}\n")
        fh.write(f"Dust threshold m/s:  {DUST_THR_MS}\n")
        fh.write(f"Score types:         {SCORE_TYPES}\n")
        fh.write(f"Distance strata:     {[s[0] for s in DIST_STRATA]}\n\n")

        # Helper for pivot summaries
        def piv(sub, idx, col, val):
            return sub.pivot_table(index=idx, columns=col, values=val).round(3)

        # 1) Sensitivity: window x angle, fixed score=max stratum=all comp=surface
        fh.write("=" * 72 + "\n")
        fh.write("1. SENSITIVITY (score=max, stratum=all, comp=surface):\n")
        fh.write("   Median r_part across trips per (window x angle), thr=7.0\n")
        s = df[(df["score_type"] == "max") & (df["stratum"] == "all") &
               (df["compartment"] == "surface") & (df["dust_thr"] == 7.0)]
        fh.write(piv(s, "window_days", "angle_tol",
                     "r_part").to_string())
        fh.write("\n\n   Same, thr=5.0:\n")
        s = df[(df["score_type"] == "max") & (df["stratum"] == "all") &
               (df["compartment"] == "surface") & (df["dust_thr"] == 5.0)]
        fh.write(piv(s, "window_days", "angle_tol",
                     "r_part").to_string())
        fh.write("\n\n   Same, thr=10.0:\n")
        s = df[(df["score_type"] == "max") & (df["stratum"] == "all") &
               (df["compartment"] == "surface") & (df["dust_thr"] == 10.0)]
        fh.write(piv(s, "window_days", "angle_tol",
                     "r_part").to_string())

        # 2) Score type comparison (fixed thr=7, ang=30, stratum=all, per comp)
        fh.write("\n\n" + "=" * 72 + "\n")
        fh.write("2. SCORE TYPE comparison (thr=7, ang=30, stratum=all):\n")
        fh.write("   Median r_part across trips per (window x score_type), per comp\n")
        for comp in ["rhizosphere", "surface", "deep"]:
            s = df[(df["dust_thr"] == 7.0) & (df["angle_tol"] == 30) &
                   (df["stratum"] == "all") & (df["compartment"] == comp)]
            fh.write(f"\n   compartment = {comp}\n")
            fh.write(piv(s, "window_days", "score_type",
                         "r_part").to_string())

        # 3) Distance stratification (fixed thr=7, ang=30, score=max, per comp)
        fh.write("\n\n" + "=" * 72 + "\n")
        fh.write("3. DISTANCE STRATIFICATION (score=max, thr=7, ang=30):\n")
        fh.write("   Median r_part across trips per (window x stratum), per comp\n")
        for comp in ["rhizosphere", "surface", "deep"]:
            s = df[(df["dust_thr"] == 7.0) & (df["angle_tol"] == 30) &
                   (df["score_type"] == "max") & (df["compartment"] == comp)]
            fh.write(f"\n   compartment = {comp}\n")
            fh.write(piv(s, "window_days", "stratum",
                         "r_part").to_string())

        # 4) Significance summary (best parameter combo per comp)
        fh.write("\n\n" + "=" * 72 + "\n")
        fh.write("4. BEST parameter combination per (compartment, window):\n")
        fh.write("   (highest median |r_part| across trips, p_part<0.05)\n")
        sig = df[df["p_part"] < 0.05]
        for comp in ["rhizosphere", "surface", "deep"]:
            for w in WINDOWS_DAYS:
                sub = sig[(sig["compartment"] == comp) & (sig["window_days"] == w)]
                if len(sub) == 0:
                    continue
                # group by parameter combo, take median |r_part| across trips
                grp = (sub.assign(absrp=sub["r_part"].abs())
                          .groupby(["score_type", "stratum", "angle_tol", "dust_thr"])
                          .agg(median_absrp=("absrp", "median"),
                               n_trips=("trip", "nunique"))
                          .reset_index()
                          .sort_values("median_absrp", ascending=False))
                if len(grp):
                    top = grp.iloc[0]
                    fh.write(f"  {comp:>11s} w={w:>3d}d: "
                             f"score={top.score_type}, stratum={top.stratum}, "
                             f"angle={top.angle_tol}°, thr={top.dust_thr} m/s "
                             f"=> median|r_part|={top.median_absrp:.3f} "
                             f"(n_trips={top.n_trips})\n")

    print(f"Wrote {OUT}/sweep_summary.txt")


if __name__ == "__main__":
    main()
