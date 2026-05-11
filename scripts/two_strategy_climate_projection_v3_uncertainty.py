#!/usr/bin/env python3
"""Climate projection v3 — uncertainty quantification.

Adds three things to v3:
  (1) Bootstrap CI on pi_B per scenario.
       Resample transitions with replacement; refit both logits on each bootstrap;
       project; report 2.5 / 50 / 97.5 percentiles per SSP.
  (2) Leave-one-trip-pair-out CV.
       For each trip-pair (e.g., T3->T4), train on remaining transitions,
       predict held-out; compute held-out AUC for both logits + projected pi_B
       under SSP1-2.6_2100.
  (3) Leave-one-site-out CV.
       Same idea, by site. Tests spatial generalization.

Outputs in cache/two_strategy_projection_v3/uncertainty/.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "two_strategy_projection_v3" / "uncertainty"
OUT.mkdir(parents=True, exist_ok=True)

SCENARIOS = {
    "Historical": (0.0, 0.0),
    "SSP1-2.6_2100": (1.6, -3.0),
    "SSP2-4.5_2100": (2.7, -5.0),
    "SSP3-7.0_2100": (4.0, -10.0),
}

FEAT_COLS = ["to_d7", "to_d365",
              "delta_d30", "delta_d90", "delta_d180",
              "to_T_d30", "to_T_d90", "to_T_d365"]

N_BOOT = 1000
RNG = np.random.default_rng(42)


def build_features() -> pd.DataFrame:
    """Build the merged transitions+T table used everywhere."""
    tr = pd.read_csv(CACHE / "transition_asymmetry" / "all_transitions.tsv", sep="\t")
    np_daily = pd.read_parquet(CACHE / "nasa_power_daily.parquet")
    np_daily["Date"] = pd.to_datetime(np_daily["Date"])

    geo = []
    for t in (1, 2, 3, 4, 5):
        g = pd.read_csv(DATA / "geodata" / f"trip{t}_geodata.tsv", sep="\t")
        g["trip"] = t
        g["Site"] = pd.to_numeric(g["Site"], errors="coerce")
        g = g.dropna(subset=["Site"])
        g["Site"] = g["Site"].astype(int)
        g["CenterDate"] = pd.to_datetime(g["CenterDate"])
        geo.append(g[["Site", "trip", "CenterDate"]])
    td = pd.concat(geo, ignore_index=True).drop_duplicates(["Site", "trip"])

    rows = []
    for _, r in td.iterrows():
        s = r["Site"]; tt = r["trip"]; cd = r["CenterDate"]
        sub_site = np_daily[np_daily["Site"] == s]
        rec = {"site": s, "trip": tt}
        for w_d, w_label in [(30, "T_d30"), (90, "T_d90"), (365, "T_d365")]:
            w = sub_site[(sub_site["Date"] >= cd - pd.Timedelta(days=w_d)) &
                           (sub_site["Date"] < cd)]
            rec[w_label] = w["TS"].mean() if len(w) else np.nan
        rows.append(rec)
    site_T = pd.DataFrame(rows)
    site_T_to = site_T.rename(columns={"site": "site", "trip": "to_trip",
                                          "T_d30": "to_T_d30", "T_d90": "to_T_d90",
                                          "T_d365": "to_T_d365"})
    site_T_from = site_T.rename(columns={"site": "site", "trip": "from_trip",
                                            "T_d30": "from_T_d30", "T_d90": "from_T_d90",
                                            "T_d365": "from_T_d365"})
    tr = tr.merge(site_T_to, on=["site", "to_trip"], how="left")
    tr = tr.merge(site_T_from, on=["site", "from_trip"], how="left")
    return tr


def fit_and_project(tr: pd.DataFrame,
                     scenarios: dict) -> pd.DataFrame:
    """Fit P(A->B) and P(B->A) on tr, project under scenarios, return one
    row per scenario with P_AtoB, P_BtoA, pi_B."""
    A = tr[tr["from_dom"] == "A"].copy()
    A["AtoB"] = (A["transition"] == "A->B").astype(int)
    A = A.dropna(subset=FEAT_COLS)
    B = tr[tr["from_dom"] == "B"].copy()
    B["BtoA"] = (B["transition"] == "B->A").astype(int)
    B = B.dropna(subset=FEAT_COLS)

    if len(A) < 30 or len(B) < 10 or A["AtoB"].sum() < 3 or B["BtoA"].sum() < 3:
        return None  # under-powered

    sA = StandardScaler().fit(A[FEAT_COLS].values)
    XA = sA.transform(A[FEAT_COLS].values)
    lrA = LogisticRegression(C=1.0, max_iter=2000, class_weight="balanced",
                                  random_state=42).fit(XA, A["AtoB"].values)

    sB = StandardScaler().fit(B[FEAT_COLS].values)
    XB = sB.transform(B[FEAT_COLS].values)
    lrB = LogisticRegression(C=1.0, max_iter=2000, class_weight="balanced",
                                  random_state=42).fit(XB, B["BtoA"].values)

    rows = []
    for scen, (dT, dP_pct) in scenarios.items():
        A_proj = A.copy()
        A_proj["to_d7"] = A_proj["to_d7"] * (1 + dP_pct / 100)
        A_proj["to_d365"] = A_proj["to_d365"] * (1 + dP_pct / 100)
        A_proj["to_T_d30"] = A_proj["to_T_d30"] + dT
        A_proj["to_T_d90"] = A_proj["to_T_d90"] + dT
        A_proj["to_T_d365"] = A_proj["to_T_d365"] + dT
        pA = lrA.predict_proba(sA.transform(A_proj[FEAT_COLS].values))[:, 1].mean()

        B_proj = B.copy()
        B_proj["to_d7"] = B_proj["to_d7"] * (1 + dP_pct / 100)
        B_proj["to_d365"] = B_proj["to_d365"] * (1 + dP_pct / 100)
        B_proj["to_T_d30"] = B_proj["to_T_d30"] + dT
        B_proj["to_T_d90"] = B_proj["to_T_d90"] + dT
        B_proj["to_T_d365"] = B_proj["to_T_d365"] + dT
        pB = lrB.predict_proba(sB.transform(B_proj[FEAT_COLS].values))[:, 1].mean()

        pi_B = pA / (pA + pB) if (pA + pB) > 0 else np.nan
        rows.append({"scenario": scen, "P_AtoB": pA, "P_BtoA": pB, "pi_B": pi_B})
    return pd.DataFrame(rows)


def fit_only(tr: pd.DataFrame):
    """Fit logits, return (scaler_A, lrA, scaler_B, lrB) or None if under-powered."""
    A = tr[tr["from_dom"] == "A"].copy()
    A["AtoB"] = (A["transition"] == "A->B").astype(int)
    A = A.dropna(subset=FEAT_COLS)
    B = tr[tr["from_dom"] == "B"].copy()
    B["BtoA"] = (B["transition"] == "B->A").astype(int)
    B = B.dropna(subset=FEAT_COLS)
    if len(A) < 30 or len(B) < 10 or A["AtoB"].sum() < 3 or B["BtoA"].sum() < 3:
        return None
    sA = StandardScaler().fit(A[FEAT_COLS].values)
    lrA = LogisticRegression(C=1.0, max_iter=2000, class_weight="balanced",
                                  random_state=42).fit(
        sA.transform(A[FEAT_COLS].values), A["AtoB"].values)
    sB = StandardScaler().fit(B[FEAT_COLS].values)
    lrB = LogisticRegression(C=1.0, max_iter=2000, class_weight="balanced",
                                  random_state=42).fit(
        sB.transform(B[FEAT_COLS].values), B["BtoA"].values)
    return sA, lrA, sB, lrB


def main():
    print("Building features (this involves NASA POWER summarization)...",
          flush=True)
    tr = build_features()
    tr_clean = tr.dropna(subset=FEAT_COLS).copy()
    print(f"  Transitions usable: {len(tr_clean)}/{len(tr)}", flush=True)
    print(f"  By from_dom: {tr_clean['from_dom'].value_counts().to_dict()}",
          flush=True)
    print(f"  Trip pairs: ", flush=True)
    print(tr_clean.groupby(['from_trip','to_trip']).size().to_string())

    # ============================================================
    # Point estimate (sanity check vs v3)
    # ============================================================
    print("\n=== Point estimate (full data) ===")
    point = fit_and_project(tr_clean, SCENARIOS)
    print(point.round(4).to_string(index=False))
    point.to_csv(OUT / "point_estimate.tsv", sep="\t", index=False)

    # ============================================================
    # (1) BOOTSTRAP CI on pi_B
    # ============================================================
    print(f"\n=== Bootstrap CI (n_boot={N_BOOT}) ===")
    boot_rows = []
    for b in range(N_BOOT):
        idx = RNG.integers(0, len(tr_clean), size=len(tr_clean))
        boot = tr_clean.iloc[idx].copy()
        proj = fit_and_project(boot, SCENARIOS)
        if proj is None:
            continue
        proj["boot"] = b
        boot_rows.append(proj)
        if (b + 1) % 100 == 0:
            print(f"  bootstrap {b+1}/{N_BOOT}", flush=True)
    boot_df = pd.concat(boot_rows, ignore_index=True)
    boot_df.to_csv(OUT / "bootstrap_raw.tsv", sep="\t", index=False)

    boot_ci = (boot_df.groupby("scenario")
                  .agg(pi_B_med=("pi_B", "median"),
                        pi_B_lo=("pi_B", lambda x: np.percentile(x, 2.5)),
                        pi_B_hi=("pi_B", lambda x: np.percentile(x, 97.5)),
                        P_AtoB_med=("P_AtoB", "median"),
                        P_AtoB_lo=("P_AtoB", lambda x: np.percentile(x, 2.5)),
                        P_AtoB_hi=("P_AtoB", lambda x: np.percentile(x, 97.5)),
                        P_BtoA_med=("P_BtoA", "median"),
                        P_BtoA_lo=("P_BtoA", lambda x: np.percentile(x, 2.5)),
                        P_BtoA_hi=("P_BtoA", lambda x: np.percentile(x, 97.5)),
                        n_boot=("pi_B", "count"))
                  .reset_index())
    # Add merged point estimate
    boot_ci = boot_ci.merge(point[["scenario", "pi_B"]]
                                .rename(columns={"pi_B": "pi_B_point"}),
                              on="scenario")
    # Re-order
    order = list(SCENARIOS.keys())
    boot_ci["_o"] = boot_ci["scenario"].map(lambda s: order.index(s))
    boot_ci = boot_ci.sort_values("_o").drop(columns="_o")
    boot_ci.to_csv(OUT / "bootstrap_ci.tsv", sep="\t", index=False)
    print(boot_ci[["scenario", "pi_B_point", "pi_B_med",
                       "pi_B_lo", "pi_B_hi", "n_boot"]]
          .round(3).to_string(index=False))

    # ============================================================
    # (2) Leave-one-trip-pair-out CV
    # ============================================================
    print("\n=== Leave-one-trip-pair-out CV ===")
    trip_pairs = sorted(set(zip(tr_clean["from_trip"], tr_clean["to_trip"])))
    print(f"  trip pairs: {trip_pairs}", flush=True)
    lto_rows = []
    for ft, tt in trip_pairs:
        held = tr_clean[(tr_clean["from_trip"] == ft) &
                          (tr_clean["to_trip"] == tt)]
        train = tr_clean[~((tr_clean["from_trip"] == ft) &
                              (tr_clean["to_trip"] == tt))]
        if len(held) < 5:
            print(f"  T{ft}->T{tt}: n_held={len(held)}, skipping", flush=True)
            continue
        fit = fit_only(train)
        if fit is None:
            print(f"  T{ft}->T{tt}: training under-powered", flush=True)
            continue
        sA, lrA, sB, lrB = fit
        # AUC on held-out
        A_h = held[held["from_dom"] == "A"].dropna(subset=FEAT_COLS)
        B_h = held[held["from_dom"] == "B"].dropna(subset=FEAT_COLS)
        aucA = aucB = np.nan
        if len(A_h) >= 5 and A_h["transition"].nunique() > 1:
            p = lrA.predict_proba(sA.transform(A_h[FEAT_COLS].values))[:, 1]
            try:
                aucA = roc_auc_score((A_h["transition"] == "A->B").astype(int), p)
            except ValueError:
                aucA = np.nan
        if len(B_h) >= 5 and B_h["transition"].nunique() > 1:
            p = lrB.predict_proba(sB.transform(B_h[FEAT_COLS].values))[:, 1]
            try:
                aucB = roc_auc_score((B_h["transition"] == "B->A").astype(int), p)
            except ValueError:
                aucB = np.nan
        # Project under SSP1-2.6_2100 using the training model
        proj = fit_and_project(train, {"SSP1-2.6_2100": SCENARIOS["SSP1-2.6_2100"],
                                              "Historical": SCENARIOS["Historical"]})
        if proj is not None:
            pi_B_hist = proj[proj["scenario"] == "Historical"]["pi_B"].iloc[0]
            pi_B_ssp = proj[proj["scenario"] == "SSP1-2.6_2100"]["pi_B"].iloc[0]
        else:
            pi_B_hist = pi_B_ssp = np.nan
        lto_rows.append({"held_out": f"T{ft}->T{tt}",
                            "n_held": int(len(held)),
                            "n_train": int(len(train)),
                            "auc_AtoB_held": aucA,
                            "auc_BtoA_held": aucB,
                            "pi_B_hist": pi_B_hist,
                            "pi_B_ssp126_2100": pi_B_ssp,
                            "delta_pi_B_ssp126": pi_B_ssp - pi_B_hist})
    lto_df = pd.DataFrame(lto_rows)
    lto_df.to_csv(OUT / "lto_trip_pair.tsv", sep="\t", index=False)
    print(lto_df.round(3).to_string(index=False))

    # ============================================================
    # (3) Leave-one-site-out CV (test spatial generalization)
    # ============================================================
    print("\n=== Leave-one-site-out CV ===")
    sites = sorted(tr_clean["site"].unique())
    print(f"  sites: {len(sites)}", flush=True)
    los_rows = []
    for s in sites:
        held = tr_clean[tr_clean["site"] == s]
        train = tr_clean[tr_clean["site"] != s]
        if len(held) < 3:
            continue
        fit = fit_only(train)
        if fit is None:
            continue
        sA, lrA, sB, lrB = fit
        A_h = held[held["from_dom"] == "A"].dropna(subset=FEAT_COLS)
        B_h = held[held["from_dom"] == "B"].dropna(subset=FEAT_COLS)
        aucA = aucB = np.nan
        if len(A_h) >= 3 and A_h["transition"].nunique() > 1:
            try:
                p = lrA.predict_proba(sA.transform(A_h[FEAT_COLS].values))[:, 1]
                aucA = roc_auc_score((A_h["transition"] == "A->B").astype(int), p)
            except ValueError:
                aucA = np.nan
        if len(B_h) >= 3 and B_h["transition"].nunique() > 1:
            try:
                p = lrB.predict_proba(sB.transform(B_h[FEAT_COLS].values))[:, 1]
                aucB = roc_auc_score((B_h["transition"] == "B->A").astype(int), p)
            except ValueError:
                aucB = np.nan
        proj = fit_and_project(train, {"SSP1-2.6_2100": SCENARIOS["SSP1-2.6_2100"],
                                              "Historical": SCENARIOS["Historical"]})
        if proj is not None:
            pi_B_hist = proj[proj["scenario"] == "Historical"]["pi_B"].iloc[0]
            pi_B_ssp = proj[proj["scenario"] == "SSP1-2.6_2100"]["pi_B"].iloc[0]
        else:
            pi_B_hist = pi_B_ssp = np.nan
        los_rows.append({"held_site": int(s),
                            "n_held": int(len(held)),
                            "n_train": int(len(train)),
                            "auc_AtoB_held": aucA,
                            "auc_BtoA_held": aucB,
                            "pi_B_hist": pi_B_hist,
                            "pi_B_ssp126_2100": pi_B_ssp,
                            "delta_pi_B_ssp126": pi_B_ssp - pi_B_hist})
    los_df = pd.DataFrame(los_rows)
    los_df.to_csv(OUT / "lto_site.tsv", sep="\t", index=False)
    print(f"  N sites with valid LOSO test: {len(los_df)}", flush=True)
    print("  AUC distribution (across held-out sites):")
    print(f"    A->B model held-out AUC: n={los_df['auc_AtoB_held'].notna().sum()}",
          f"med={los_df['auc_AtoB_held'].median():.3f}",
          f"IQR=[{los_df['auc_AtoB_held'].quantile(0.25):.3f}, "
          f"{los_df['auc_AtoB_held'].quantile(0.75):.3f}]")
    print(f"    B->A model held-out AUC: n={los_df['auc_BtoA_held'].notna().sum()}",
          f"med={los_df['auc_BtoA_held'].median():.3f}",
          f"IQR=[{los_df['auc_BtoA_held'].quantile(0.25):.3f}, "
          f"{los_df['auc_BtoA_held'].quantile(0.75):.3f}]")
    print("  pi_B projections across held-out sites (SSP1-2.6_2100):")
    print(f"    delta_pi_B: med={los_df['delta_pi_B_ssp126'].median():.3f}",
          f"IQR=[{los_df['delta_pi_B_ssp126'].quantile(0.25):.3f}, "
          f"{los_df['delta_pi_B_ssp126'].quantile(0.75):.3f}]")
    print(f"    fraction with positive delta (warming -> more B): "
          f"{(los_df['delta_pi_B_ssp126'] > 0).mean():.2%}")

    # ============================================================
    # Summary file
    # ============================================================
    summary_lines = []
    summary_lines.append("Climate projection v3 — uncertainty bundle\n")
    summary_lines.append("=" * 60 + "\n\n")
    summary_lines.append("(1) Bootstrap CI (1000 resamples)\n")
    summary_lines.append(boot_ci[["scenario", "pi_B_point", "pi_B_med",
                                          "pi_B_lo", "pi_B_hi"]]
                              .round(3).to_string(index=False) + "\n\n")
    summary_lines.append("(2) Leave-one-trip-pair-out\n")
    summary_lines.append(lto_df.round(3).to_string(index=False) + "\n\n")
    summary_lines.append(f"(3) Leave-one-site-out (n={len(los_df)} sites)\n")
    summary_lines.append(f"  A->B held-out AUC: med={los_df['auc_AtoB_held'].median():.3f} "
                              f"IQR=[{los_df['auc_AtoB_held'].quantile(0.25):.3f}, "
                              f"{los_df['auc_AtoB_held'].quantile(0.75):.3f}]\n")
    summary_lines.append(f"  B->A held-out AUC: med={los_df['auc_BtoA_held'].median():.3f} "
                              f"IQR=[{los_df['auc_BtoA_held'].quantile(0.25):.3f}, "
                              f"{los_df['auc_BtoA_held'].quantile(0.75):.3f}]\n")
    summary_lines.append(f"  pi_B SSP1-2.6_2100 delta: med={los_df['delta_pi_B_ssp126'].median():.3f} "
                              f"IQR=[{los_df['delta_pi_B_ssp126'].quantile(0.25):.3f}, "
                              f"{los_df['delta_pi_B_ssp126'].quantile(0.75):.3f}]\n")
    summary_lines.append(f"  Fraction of held-out sites with delta>0: "
                              f"{(los_df['delta_pi_B_ssp126'] > 0).mean():.2%}\n")
    with open(OUT / "summary.txt", "w") as fh:
        fh.write("".join(summary_lines))
    print("\nDONE. Outputs in", OUT)


if __name__ == "__main__":
    main()
