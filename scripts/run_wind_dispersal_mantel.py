#!/usr/bin/env python3
"""Wind-connectivity Mantel: test whether iCAMP-detected homogenizing
dispersal in the EQ is mechanistically wind-driven.

Pipeline:
  1. Per (site, trip): cumulative "dust-uplift hours" (WS10M_MAX > thr)
                        and per-day wind direction (degrees, met convention).
  2. Per (pair, trip): wind-connectivity score = sum over the window of
       WS10M_MAX_i,d * 1[direction wind blew TOWARD site_j ± tol]
     (max of i->j and j->i, symmetric).
  3. Per compartment x trip: site-aggregated relative-abundance vectors,
     pairwise Bray-Curtis dissimilarity.
  4. Mantel: BC ~ geographic distance.
     Partial Mantel: BC ~ wind-connectivity controlling for distance.
     Per-compartment, per-trip and pooled.

Inputs:
  data/climate/daily_weather_full.csv      87,661 rows; WS10M_MAX, WD10M, ...
  cache/pairwise_geometry.tsv              1,770 pairs distance + bearings
  cache/feature_table.parquet              ASV x sample
  cache/taxonomy.parquet
  data/geodata/trip{1..5}_geodata.tsv      CenterDate per (site, trip)

Outputs:
  cache/wind_dispersal/site_trip_wind_summary.tsv
  cache/wind_dispersal/pair_trip_wind_connectivity.tsv
  cache/wind_dispersal/bc_per_compartment_trip.parquet  (multi-index DF)
  cache/wind_dispersal/mantel_results.tsv
  cache/wind_dispersal/sensitivity.tsv
  cache/wind_dispersal/summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "wind_dispersal"
OUT.mkdir(parents=True, exist_ok=True)

# Parameters
WINDOWS_DAYS = [30, 90, 365]
DUST_THR_MS = 7.0          # WS10M_MAX threshold for dust uplift (Shao 2008)
ANGLE_TOL_DEG = 30.0       # alignment tolerance for "wind blows TOWARD j"
N_PERM = 999               # Mantel permutations
RNG = np.random.default_rng(20260509)


def load_geo() -> pd.DataFrame:
    geo = []
    for trip in range(1, 6):
        gp = DATA / "geodata" / f"trip{trip}_geodata.tsv"
        if not gp.exists():
            continue
        g = pd.read_csv(gp, sep="\t")
        g["trip"] = trip
        g["CenterDate"] = pd.to_datetime(g["CenterDate"])
        g["Site_int"] = pd.to_numeric(g["Site"], errors="coerce")
        g = g.dropna(subset=["Site_int", "CenterDate"])
        g["site"] = g["Site_int"].astype(int)
        geo.append(g[["site", "trip", "CenterDate"]])
    return pd.concat(geo).drop_duplicates(["site", "trip"])


def angular_diff_deg(a, b):
    """Smallest angle between two directions in degrees, in [0, 180]."""
    d = np.abs((a - b + 180) % 360 - 180)
    return d


def wind_summaries(weather: pd.DataFrame, geo: pd.DataFrame,
                    window_days: int) -> pd.DataFrame:
    """Per (site, trip): wind statistics over the trip's CenterDate window."""
    w = weather.set_index(["site", "date"]).sort_index()
    rows = []
    for _, r in geo.iterrows():
        site = int(r["site"])
        d = pd.Timestamp(r["CenterDate"])
        try:
            site_w = w.loc[site]
        except KeyError:
            continue
        mask = (site_w.index >= d - pd.Timedelta(days=window_days - 1)) & \
               (site_w.index <= d)
        win = site_w.loc[mask]
        if len(win) == 0:
            continue
        dust_hours = float((win["WS10M_MAX"] > DUST_THR_MS).sum())
        # Vector mean wind direction (met convention: 0=N, 90=E, ...)
        # The direction wind is blowing TOWARD = (WD + 180) % 360
        wd_to = (win["WD10M"].values + 180) % 360
        ux = np.mean(np.cos(np.deg2rad(wd_to)))
        uy = np.mean(np.sin(np.deg2rad(wd_to)))
        mean_dir_to = (np.degrees(np.arctan2(uy, ux)) + 360) % 360
        resultant_len = float(np.sqrt(ux**2 + uy**2))
        rows.append({"site": site, "trip": int(r["trip"]),
                     "window_days": window_days,
                     "n_days": int(len(win)),
                     "dust_hours": dust_hours,
                     "mean_WS10M_MAX": float(win["WS10M_MAX"].mean()),
                     "max_WS10M_MAX": float(win["WS10M_MAX"].max()),
                     "mean_dir_to_deg": float(mean_dir_to),
                     "resultant_length": resultant_len})
    return pd.DataFrame(rows)


def pair_connectivity(weather: pd.DataFrame, geo: pd.DataFrame,
                       geom: pd.DataFrame, window_days: int) -> pd.DataFrame:
    """Per (site_i, site_j, trip): wind-connectivity score from i->j and j->i."""
    w = weather.set_index(["site", "date"]).sort_index()

    # Pre-fetch per-(site, trip) wind windows
    win_cache: dict = {}
    for _, r in geo.iterrows():
        site = int(r["site"]); trip = int(r["trip"])
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

    # Geometry: bearings_ij and bearings_ji
    geom_lkp = {(int(r.site_i), int(r.site_j)): (r.bearing_ij_deg, r.bearing_ji_deg)
                for r in geom.itertuples()}

    rows = []
    for _, r in geom.iterrows():
        si, sj = int(r["site_i"]), int(r["site_j"])
        bij, bji = float(r["bearing_ij_deg"]), float(r["bearing_ji_deg"])
        for trip in range(1, 6):
            wi = win_cache.get((si, trip))
            wj = win_cache.get((sj, trip))
            if wi is None or wj is None or len(wi) == 0 or len(wj) == 0:
                continue
            # i -> j: at site_i, did wind blow toward j? (within tol)
            ws_i = wi["WS10M_MAX"].values
            di_to = wi["dir_to"].values
            aligned_i_to_j = angular_diff_deg(di_to, bij) <= ANGLE_TOL_DEG
            uplift_i = ws_i > DUST_THR_MS
            score_ij = float(np.sum(ws_i[aligned_i_to_j & uplift_i]))
            # j -> i
            ws_j = wj["WS10M_MAX"].values
            dj_to = wj["dir_to"].values
            aligned_j_to_i = angular_diff_deg(dj_to, bji) <= ANGLE_TOL_DEG
            uplift_j = ws_j > DUST_THR_MS
            score_ji = float(np.sum(ws_j[aligned_j_to_i & uplift_j]))
            rows.append({"site_i": si, "site_j": sj, "trip": trip,
                         "window_days": window_days,
                         "score_ij": score_ij, "score_ji": score_ji,
                         "score_max": max(score_ij, score_ji),
                         "score_sum": score_ij + score_ji,
                         "dist_km": float(r["dist_km"])})
    return pd.DataFrame(rows)


def site_trip_communities(ft: pd.DataFrame, smeta: pd.DataFrame,
                            compartment: str) -> dict:
    """For each (site, trip) within a compartment, return mean relative-abundance
    vector (averaged across replicates). Returns dict (site, trip) -> 1-D vec."""
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
    """Bray-Curtis dissimilarity matrix among the given sites for this trip."""
    n = len(sites)
    M = np.full((n, n), np.nan)
    for i in range(n):
        ci = comm.get((sites[i], trip))
        if ci is None:
            continue
        for j in range(i + 1, n):
            cj = comm.get((sites[j], trip))
            if cj is None:
                continue
            num = np.abs(ci - cj).sum()
            den = ci.sum() + cj.sum()
            d = num / den if den > 0 else np.nan
            M[i, j] = M[j, i] = d
    np.fill_diagonal(M, 0.0)
    return M


def mantel(D1: np.ndarray, D2: np.ndarray, n_perm: int = N_PERM) -> tuple:
    """Standard Mantel test on two distance matrices. Returns (r, p)."""
    n = D1.shape[0]
    iu = np.triu_indices(n, k=1)
    x = D1[iu]; y = D2[iu]
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 5:
        return (np.nan, np.nan, int(valid.sum()))
    x = x[valid]; y = y[valid]
    r_obs = np.corrcoef(x, y)[0, 1]
    # Permutation: shuffle rows/cols of D2
    cnt = 0
    for _ in range(n_perm):
        perm = RNG.permutation(n)
        D2p = D2[perm][:, perm]
        yp = D2p[iu][valid]
        r = np.corrcoef(x, yp)[0, 1]
        if abs(r) >= abs(r_obs):
            cnt += 1
    p = (cnt + 1) / (n_perm + 1)
    return (float(r_obs), float(p), int(valid.sum()))


def partial_mantel(D_y: np.ndarray, D_x: np.ndarray, D_z: np.ndarray,
                    n_perm: int = N_PERM) -> tuple:
    """Partial Mantel: r(D_y, D_x | D_z). Permutes rows/cols of D_x."""
    n = D_y.shape[0]
    iu = np.triu_indices(n, k=1)
    y = D_y[iu]; x = D_x[iu]; z = D_z[iu]
    valid = np.isfinite(y) & np.isfinite(x) & np.isfinite(z)
    if valid.sum() < 5:
        return (np.nan, np.nan, int(valid.sum()))
    y = y[valid]; x = x[valid]; z = z[valid]

    def partial_r(yv, xv, zv):
        ryx = np.corrcoef(yv, xv)[0, 1]
        ryz = np.corrcoef(yv, zv)[0, 1]
        rxz = np.corrcoef(xv, zv)[0, 1]
        denom = np.sqrt((1 - ryz**2) * (1 - rxz**2))
        return (ryx - ryz * rxz) / denom if denom > 0 else np.nan

    r_obs = partial_r(y, x, z)
    cnt = 0
    for _ in range(n_perm):
        perm = RNG.permutation(n)
        Dxp = D_x[perm][:, perm]
        xp = Dxp[iu][valid]
        r = partial_r(y, xp, z)
        if r == r and abs(r) >= abs(r_obs):
            cnt += 1
    p = (cnt + 1) / (n_perm + 1)
    return (float(r_obs), float(p), int(valid.sum()))


def main():
    print("[load] feature_table, taxonomy, geometry, weather, geodata", flush=True)
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    geom = pd.read_csv(CACHE / "pairwise_geometry.tsv", sep="\t")
    weather = pd.read_csv(DATA / "climate" / "daily_weather_full.csv")
    weather["date"] = pd.to_datetime(weather["date"], format="%Y%m%d")
    weather["site"] = weather["site"].astype(int)
    geo = load_geo()
    smeta = parse_samples_to_df(ft.columns)

    # Per-site-trip wind summaries (one per window length)
    wind_summary_all = []
    for w in WINDOWS_DAYS:
        s = wind_summaries(weather, geo, w)
        wind_summary_all.append(s)
        print(f"  wind summary w={w}d: {len(s)} (site, trip) rows", flush=True)
    pd.concat(wind_summary_all).to_csv(OUT / "site_trip_wind_summary.tsv",
                                        sep="\t", index=False)

    # Pair wind-connectivity (use 90-day window as primary, plus 30 and 365)
    pair_conn = []
    for w in WINDOWS_DAYS:
        pc = pair_connectivity(weather, geo, geom, w)
        pair_conn.append(pc)
        print(f"  pair connectivity w={w}d: {len(pc)} (pair, trip) rows", flush=True)
    pd.concat(pair_conn).to_csv(OUT / "pair_trip_wind_connectivity.tsv",
                                 sep="\t", index=False)

    # Build BC matrices per (compartment, trip)
    sites = sorted(set(int(s) for s in geo["site"].unique()))
    print(f"\n[BC] computing site-aggregated BC per (comp, trip)...", flush=True)

    bc_records = []
    bc_matrices = {}  # (comp, trip) -> ndarray
    for comp in ["rhizosphere", "surface", "deep"]:
        comm = site_trip_communities(ft, smeta, comp)
        for trip in range(1, 6):
            bc = pairwise_bc(comm, sites, trip)
            n_valid = np.isfinite(bc[np.triu_indices(len(sites), k=1)]).sum()
            print(f"  {comp:>11s} trip {trip}: {n_valid} valid pairs", flush=True)
            bc_matrices[(comp, trip)] = bc
            for i, si in enumerate(sites):
                for j, sj in enumerate(sites):
                    if j > i:
                        bc_records.append({"compartment": comp, "trip": trip,
                                            "site_i": si, "site_j": sj,
                                            "bc": float(bc[i, j])})
    bc_df = pd.DataFrame(bc_records)
    bc_df.to_parquet(OUT / "bc_per_compartment_trip.parquet")

    # Build distance matrix
    n_s = len(sites)
    site_idx = {s: i for i, s in enumerate(sites)}
    D_dist = np.zeros((n_s, n_s))
    for _, r in geom.iterrows():
        i = site_idx[int(r["site_i"])]
        j = site_idx[int(r["site_j"])]
        D_dist[i, j] = D_dist[j, i] = float(r["dist_km"])

    # Build per-(window, trip) wind-connectivity matrices and run Mantel
    print(f"\n[mantel] running per-(comp, trip, window) Mantel + partial Mantel",
          flush=True)
    mantel_rows = []
    for w in WINDOWS_DAYS:
        pc_w = pd.concat(pair_conn).query("window_days == @w")
        for trip in range(1, 6):
            sub = pc_w.query("trip == @trip")
            if len(sub) == 0:
                continue
            # Build wind-connectivity matrix (use score_max as primary).
            # Convert to a "dissimilarity" by negation: smaller = more connected
            D_wind = np.full((n_s, n_s), np.nan)
            for _, r in sub.iterrows():
                i = site_idx[int(r["site_i"])]
                j = site_idx[int(r["site_j"])]
                D_wind[i, j] = D_wind[j, i] = -float(r["score_max"])
            np.fill_diagonal(D_wind, 0.0)

            for comp in ["rhizosphere", "surface", "deep"]:
                D_bc = bc_matrices[(comp, trip)]
                # Skip if BC mostly NaN (this comp/trip has few sites)
                iu = np.triu_indices(n_s, k=1)
                n_valid_bc = np.isfinite(D_bc[iu]).sum()
                if n_valid_bc < 30:
                    continue
                r_dist, p_dist, n_dist = mantel(D_bc, D_dist)
                r_wind, p_wind, n_wind = mantel(D_bc, D_wind)
                r_part, p_part, n_part = partial_mantel(D_bc, D_wind, D_dist)
                mantel_rows.append({"compartment": comp, "trip": trip,
                                     "window_days": w,
                                     "n_pairs": int(n_valid_bc),
                                     "r_BC_dist": r_dist, "p_BC_dist": p_dist,
                                     "r_BC_wind": r_wind, "p_BC_wind": p_wind,
                                     "r_BC_wind|dist": r_part,
                                     "p_BC_wind|dist": p_part})
    mantel_df = pd.DataFrame(mantel_rows)
    mantel_df.to_csv(OUT / "mantel_results.tsv", sep="\t", index=False)
    print(f"\nMantel results: {len(mantel_df)} rows", flush=True)
    print(mantel_df.to_string(index=False))

    # Brief summary
    with open(OUT / "summary.txt", "w") as fh:
        fh.write("Wind-connectivity Mantel (Tier-1 followup)\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Sites: {n_s}\n")
        fh.write(f"Windows tested (days): {WINDOWS_DAYS}\n")
        fh.write(f"Dust uplift threshold: WS10M_MAX > {DUST_THR_MS} m/s\n")
        fh.write(f"Angular tolerance: {ANGLE_TOL_DEG}°\n")
        fh.write(f"Mantel permutations: {N_PERM}\n\n")

        fh.write("Mantel: BC ~ geographic distance (baseline)\n")
        b = mantel_df[mantel_df["window_days"] == 90].pivot_table(
            index="compartment", columns="trip", values="r_BC_dist")
        fh.write(b.round(3).to_string())
        fh.write("\n\nMantel: BC ~ wind-connectivity (90-day window)\n")
        b2 = mantel_df[mantel_df["window_days"] == 90].pivot_table(
            index="compartment", columns="trip", values="r_BC_wind")
        fh.write(b2.round(3).to_string())
        fh.write("\n\nPartial Mantel: BC ~ wind | distance (90-day window)\n")
        b3 = mantel_df[mantel_df["window_days"] == 90].pivot_table(
            index="compartment", columns="trip", values="r_BC_wind|dist")
        fh.write(b3.round(3).to_string())
        fh.write("\n\nFraction of significant tests (p<0.05) per compartment:\n")
        for comp in ["rhizosphere", "surface", "deep"]:
            sub = mantel_df[mantel_df["compartment"] == comp]
            fh.write(f"  {comp:>11s}: distance "
                     f"{(sub['p_BC_dist']<0.05).mean():.0%}; "
                     f"wind {(sub['p_BC_wind']<0.05).mean():.0%}; "
                     f"partial {(sub['p_BC_wind|dist']<0.05).mean():.0%}\n")
    print(f"\nWrote {OUT}/summary.txt")


if __name__ == "__main__":
    main()
