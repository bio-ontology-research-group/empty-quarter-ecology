#!/usr/bin/env python3
"""Critical sensitivity tests on the relic-likelihood indicator.

Tests:
  1. SITE-LEVEL CV: train on C1 PMA pairs only, test on C2 PMA pairs only
     (and reverse). Estimates true cross-site generalization.
  2. SEQUENCE-ONLY MODEL: refit using only sequence-based features
     (gc_content, asv_len, length_dev, pyr_dinuc_density, tc_cc_density,
     intra_otu_min_pid, is_singleton_otu). Compare AUC and re-score all
     ASVs.
  3. RANDOM-SUBSET NULL: 100 random subsets of size 863 from the ASV
     pool. Compute climate-Shannon Spearman and Allison-Martiny slope per
     subset. Compare distributions to alive subset.
  4. THRESHOLD SENSITIVITY: relic_score thresholds 0.5/0.6/0.7/0.8.
  5. DIRECT PMA EVIDENCE FOR CSP1-2: of 24 CSP1-2 ASVs, how many were
     in PMA samples? What were their T/UT ratios?

Outputs in cache/relic_sensitivity/.
"""
from __future__ import annotations

import sys
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, linregress
from scipy.spatial.distance import pdist, squareform
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
RELIC = Path("/home/leechuck/Public/software/empty-quarter/relic-dna")
OUT = CACHE / "relic_sensitivity"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260510)


def shannon(arr):
    a = arr[arr > 0]
    if len(a) == 0: return 0.0
    p = a / a.sum()
    return float(-(p * np.log(p)).sum())


def relabund(M):
    return M.div(M.sum(axis=0).replace(0, 1), axis=1)


def parse_pma_sample(s: str):
    m = re.match(r"^(C[12])([RS])(\d+)(T|UT)$", s)
    if not m: return None
    return {"site": m.group(1), "comp": m.group(2),
              "rep": int(m.group(3)), "treat": m.group(4)}


# =============================================================================
# 1. Site-level CV (C1 vs C2)
# =============================================================================
def site_level_cv():
    print("\n=== [1] Site-level CV (train C1, test C2 + reverse) ===",
          flush=True)
    feats = pd.read_parquet(CACHE / "test6_disconfirmation" /
                              "relic_features_with_damage.parquet")

    # Compute per-EQ-ASV viability for C1 only and C2 only separately
    rel = pd.read_csv(RELIC / "ASV_table_rel_abundance.tsv", sep="\t",
                       index_col=0)
    samp = []
    for s in rel.columns:
        m = parse_pma_sample(s)
        if m: m["sample_orig"] = s; samp.append(m)
    sm = pd.DataFrame(samp)
    pairs = (sm.pivot(index=["site", "comp", "rep"],
                        columns="treat", values="sample_orig")
             .dropna(subset=["T", "UT"]).reset_index())

    # Map PMA ASV -> EQ ASV
    map_df = pd.read_csv(CACHE / "test6_disconfirmation" /
                           "pma_to_eq_match.tsv", sep="\t", header=None,
                           names=["pma_asv", "eq_asv", "pid", "alen", "mm",
                                    "gaps", "qs", "qe", "ts", "te", "ev",
                                    "bs"])
    pma_to_eq = (map_df.drop_duplicates("pma_asv")
                 [["pma_asv", "eq_asv", "pid"]])

    eps = 1e-7
    rec = []
    for _, p in pairs.iterrows():
        T = rel[p["T"]]; UT = rel[p["UT"]]
        det = UT > 1e-5
        for asv in rel.index[det]:
            r = (float(T[asv]) + eps) / (float(UT[asv]) + eps)
            rec.append({"pma_asv": asv, "site": p["site"], "ratio":
                          float(min(r, 10.0)), "ut": float(UT[asv])})
    R = pd.DataFrame(rec)
    R = R.merge(pma_to_eq, on="pma_asv", how="left").dropna(subset=["eq_asv"])

    # Per (eq_asv, site): weighted median ratio
    def weighted_med(g):
        if g["ut"].sum() == 0: return np.nan
        return float(np.average(g["ratio"], weights=g["ut"]))
    by_site = (R.groupby(["eq_asv", "site"])
               .apply(weighted_med, include_groups=False)
               .rename("v").reset_index())
    by_site = by_site.dropna()
    print(f"  per-(ASV, site) viability records: {len(by_site)}",
          flush=True)

    fcols = [
        "persistence_max", "persistence_mean", "n_site_records",
        "log_mean_abund", "log_max_abund", "n_detections", "n_sites_detected",
        "frac_deep", "frac_surface", "frac_rhizosphere",
        "emp_cosmo_90", "emp_cosmo_100", "emp_cosmo_150",
        "asv_len", "gc_content", "length_dev", "pyr_dinuc_density",
        "tc_cc_density", "cluster_fanout", "is_singleton_otu",
        "intra_otu_min_pid",
    ]

    feats = feats.dropna(subset=fcols).copy()
    sc_data = feats.set_index("asv_id")[fcols]

    rows = []
    for train_site, test_site in (("C1", "C2"), ("C2", "C1")):
        tr = by_site[by_site["site"] == train_site]
        te = by_site[by_site["site"] == test_site]
        # Outside ambig band labels
        tr = tr.assign(
            y=lambda x: np.where(x["v"] < 0.1, 1,
                                    np.where(x["v"] > 0.5, 0, np.nan)))
        te = te.assign(
            y=lambda x: np.where(x["v"] < 0.1, 1,
                                    np.where(x["v"] > 0.5, 0, np.nan)))
        tr = tr.dropna(subset=["y"])
        te = te.dropna(subset=["y"])
        # Restrict to ASVs with features
        tr = tr[tr["eq_asv"].isin(sc_data.index)]
        te = te[te["eq_asv"].isin(sc_data.index)]
        Xtr = sc_data.loc[tr["eq_asv"]].values
        Xte = sc_data.loc[te["eq_asv"]].values
        Xtr = np.where(np.isfinite(Xtr), Xtr, 0)
        Xte = np.where(np.isfinite(Xte), Xte, 0)
        ytr = tr["y"].values.astype(int)
        yte = te["y"].values.astype(int)
        if len(set(ytr)) < 2 or len(set(yte)) < 2: continue
        sc = StandardScaler().fit(Xtr)
        Xtr_s = sc.transform(Xtr); Xte_s = sc.transform(Xte)
        gb = GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                            learning_rate=0.05,
                                            random_state=42).fit(Xtr, ytr)
        lr = LogisticRegression(C=1.0, max_iter=2000,
                                  class_weight="balanced",
                                  random_state=42).fit(Xtr_s, ytr)
        pgb = gb.predict_proba(Xte)[:, 1]
        plr = lr.predict_proba(Xte_s)[:, 1]
        auc_gb = roc_auc_score(yte, pgb)
        auc_lr = roc_auc_score(yte, plr)
        ap_gb = average_precision_score(yte, pgb)
        # Compare to in-distribution: 5-fold CV on train_site
        if len(set(ytr)) == 2 and len(ytr) >= 30:
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
            p_in = cross_val_predict(GradientBoostingClassifier(
                n_estimators=200, max_depth=3, learning_rate=0.05,
                random_state=42), Xtr, ytr, cv=cv, method="predict_proba")[:, 1]
            auc_in = roc_auc_score(ytr, p_in)
        else:
            auc_in = np.nan
        rows.append({"train_site": train_site, "test_site": test_site,
                      "n_train": len(ytr), "n_test": len(yte),
                      "auc_gb_holdout": auc_gb,
                      "auc_lr_holdout": auc_lr,
                      "ap_gb_holdout": ap_gb,
                      "auc_gb_in_dist_5fold": auc_in})
        print(f"  train={train_site} (n={len(ytr)}) -> test={test_site}"
              f" (n={len(yte)})", flush=True)
        print(f"    AUC GB held-out:    {auc_gb:.3f}", flush=True)
        print(f"    AUC LR held-out:    {auc_lr:.3f}", flush=True)
        print(f"    AUC GB in-dist 5CV: {auc_in:.3f}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "site_level_cv.tsv", sep="\t", index=False)
    return df


# =============================================================================
# 2. Sequence-only model
# =============================================================================
def sequence_only_model():
    print("\n=== [2] Sequence-only model ===", flush=True)
    feats = pd.read_parquet(CACHE / "test6_disconfirmation" /
                              "relic_features_with_damage.parquet")
    seq_cols = ["asv_len", "gc_content", "length_dev",
                  "pyr_dinuc_density", "tc_cc_density",
                  "intra_otu_min_pid", "is_singleton_otu", "cluster_fanout"]
    train_mask = (feats["weighted_median_ratio"].notna() &
                    ((feats["weighted_median_ratio"] < 0.1) |
                     (feats["weighted_median_ratio"] > 0.5)))
    df = feats[train_mask].copy()
    df["y"] = (df["weighted_median_ratio"] < 0.1).astype(int)
    df = df.dropna(subset=seq_cols)
    print(f"  train n: {len(df)}", flush=True)
    X = df[seq_cols].values
    X = np.where(np.isfinite(X), X, 0.0)
    y = df["y"].values
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    gb = GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                       learning_rate=0.05, random_state=42)
    p = cross_val_predict(gb, X, y, cv=cv, method="predict_proba")[:, 1]
    auc = roc_auc_score(y, p)
    ap = average_precision_score(y, p)
    print(f"  GB CV AUC (sequence only): {auc:.3f}  AP: {ap:.3f}", flush=True)

    # Refit + score all ASVs
    gb_final = GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                              learning_rate=0.05,
                                              random_state=42).fit(X, y)
    full = feats.dropna(subset=seq_cols).copy()
    Xall = full[seq_cols].values
    Xall = np.where(np.isfinite(Xall), Xall, 0.0)
    full["relic_score_seq_only"] = gb_final.predict_proba(Xall)[:, 1]
    full[["asv_id", "relic_score_seq_only"]].to_csv(
        OUT / "asv_relic_score_seq_only.tsv", sep="\t", index=False)
    print(f"  Distribution of seq-only relic_score:", flush=True)
    for q in (10, 25, 50, 75, 90):
        print(f"    p{q}: "
              f"{np.percentile(full['relic_score_seq_only'], q):.3f}",
              flush=True)
    return {"seq_only_auc": auc, "seq_only_ap": ap}


# =============================================================================
# 3. Random-subset null comparison
# =============================================================================
def random_subset_null():
    print("\n=== [3] Random-subset null for climate-Shannon + AM-slope ===",
          flush=True)
    ind = pd.read_csv(CACHE / "test6_disconfirmation" /
                       "relic_indicator_with_damage_per_asv.tsv", sep="\t")
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)
    geo_t1 = pd.read_csv(DATA / "geodata" / "trip1_geodata.tsv", sep="\t")
    geo_t1 = geo_t1.rename(columns={"Site": "site"})
    geo_t1["site"] = pd.to_numeric(geo_t1["site"],
                                       errors="coerce").astype("Int64")

    alive_size = 863
    n_iter = 100

    # Restrict to ASVs that exist in the ft
    candidate_asvs = ft.index.tolist()

    # Helper: compute Shannon vs MAT spearman for given subset
    def shannon_climate_for_subset(asv_subset):
        sub_ft = ft.loc[ft.index.isin(asv_subset)]
        sh = pd.DataFrame({"sample": sub_ft.columns,
                              "shannon": [shannon(sub_ft[c].values)
                                            for c in sub_ft.columns]})
        sh = sh.merge(smeta, on="sample", how="left")
        sh = sh.merge(geo_t1[["site", "AnnualMeanTemp"]], on="site",
                        how="left").dropna(subset=["AnnualMeanTemp"])
        if len(sh) < 30: return np.nan
        r, _ = spearmanr(sh["shannon"], sh["AnnualMeanTemp"])
        return float(r)

    # Allison-Martiny slope for a subset (just one compartment-trip for speed)
    path_func = pd.read_csv(DATA / "functional" / "picrust2" /
                              "path_abun_unstrat.tsv", sep="\t",
                              index_col=0)
    rsamps = list(set(smeta[(smeta["compartment"] == "rhizosphere") &
                                (smeta["trip"] == 3)]["sample"]) &
                    set(ft.columns) & set(path_func.columns))
    print(f"  using {len(rsamps)} samples for AM slope test (rhizo trip 3)",
          flush=True)
    Fr_ref = relabund(path_func[rsamps]).T
    BC_func_ref = pdist(Fr_ref.values, metric="braycurtis")

    def am_slope_for_subset(asv_subset):
        sub_ft = ft.loc[ft.index.isin(asv_subset), rsamps]
        Tr = relabund(sub_ft).T
        Tr = Tr[Tr.sum(axis=1) > 0]
        s2 = list(Tr.index)
        if len(s2) < 10: return np.nan
        bc_t = pdist(Tr.values, metric="braycurtis")
        Fr = Fr_ref.loc[s2]
        bc_f = pdist(Fr.values, metric="braycurtis")
        mask = (~np.isnan(bc_t)) & (~np.isnan(bc_f))
        if mask.sum() < 30: return np.nan
        slope, _, _, _, _ = linregress(bc_t[mask], bc_f[mask])
        return float(slope)

    # Alive baseline
    alive_set = set(ind.loc[ind["relic_score_full_gb"] <= 0.3, "asv_id"])
    alive_climate = shannon_climate_for_subset(alive_set)
    alive_am = am_slope_for_subset(alive_set)
    print(f"  ALIVE Shannon vs MAT rho:        {alive_climate:+.3f}",
          flush=True)
    print(f"  ALIVE Allison-Martiny slope:     {alive_am:+.3f}", flush=True)

    # All baseline
    all_climate = shannon_climate_for_subset(set(candidate_asvs))
    all_am = am_slope_for_subset(set(candidate_asvs))
    print(f"  ALL Shannon vs MAT rho:          {all_climate:+.3f}",
          flush=True)
    print(f"  ALL Allison-Martiny slope:       {all_am:+.3f}", flush=True)

    # Random subsets of same size
    print(f"  computing {n_iter} random subsets of size {alive_size} ...",
          flush=True)
    rand_climate = []
    rand_am = []
    for k in range(n_iter):
        sub = set(RNG.choice(candidate_asvs, alive_size, replace=False))
        rand_climate.append(shannon_climate_for_subset(sub))
        rand_am.append(am_slope_for_subset(sub))
    rand_climate = np.array([x for x in rand_climate if not np.isnan(x)])
    rand_am = np.array([x for x in rand_am if not np.isnan(x)])
    print(f"  RAND Shannon vs MAT (median):    {np.median(rand_climate):+.3f}"
          f"  [p5..p95: {np.percentile(rand_climate, 5):+.3f}, "
          f"{np.percentile(rand_climate, 95):+.3f}]", flush=True)
    print(f"  RAND Allison-Martiny (median):   {np.median(rand_am):+.3f}"
          f"  [p5..p95: {np.percentile(rand_am, 5):+.3f}, "
          f"{np.percentile(rand_am, 95):+.3f}]", flush=True)

    # Where does alive sit in random distribution?
    p_climate = float((rand_climate > alive_climate).mean())
    p_am = float((rand_am < alive_am).mean())
    print(f"\n  Alive Shannon-MAT rho ({alive_climate:+.3f}) is in "
          f"{p_climate:.0%} of random rhos (p={1-p_climate:.3f} for "
          f"more-negative)", flush=True)
    print(f"  Alive AM slope ({alive_am:+.3f}) is in {p_am:.0%} of "
          f"random slopes (p={p_am:.3f} for less-positive)", flush=True)

    pd.DataFrame({"k": range(len(rand_climate)),
                    "rand_shannon_mat_rho": rand_climate,
                    }).to_csv(OUT / "null_random_climate.tsv",
                                sep="\t", index=False)
    pd.DataFrame({"k": range(len(rand_am)),
                    "rand_am_slope": rand_am,
                    }).to_csv(OUT / "null_random_am.tsv",
                                sep="\t", index=False)
    return {"alive_climate_rho": alive_climate, "all_climate_rho": all_climate,
              "rand_climate_median": float(np.median(rand_climate)),
              "rand_climate_p5": float(np.percentile(rand_climate, 5)),
              "rand_climate_p95": float(np.percentile(rand_climate, 95)),
              "alive_am_slope": alive_am, "all_am_slope": all_am,
              "rand_am_median": float(np.median(rand_am)),
              "rand_am_p5": float(np.percentile(rand_am, 5)),
              "rand_am_p95": float(np.percentile(rand_am, 95)),
              "p_alive_climate_extreme_negative":
                  float(1 - p_climate),
              "p_alive_am_extreme_low": float(p_am)}


# =============================================================================
# 4. Threshold sensitivity
# =============================================================================
def threshold_sensitivity():
    print("\n=== [4] Threshold sensitivity (climate gradient + AM slope)"
          " ===", flush=True)
    ind = pd.read_csv(CACHE / "test6_disconfirmation" /
                       "relic_indicator_with_damage_per_asv.tsv", sep="\t")
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)
    geo_t1 = pd.read_csv(DATA / "geodata" / "trip1_geodata.tsv", sep="\t")
    geo_t1 = geo_t1.rename(columns={"Site": "site"})
    geo_t1["site"] = pd.to_numeric(geo_t1["site"],
                                       errors="coerce").astype("Int64")

    rows = []
    for thresh in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
        alive = set(ind.loc[ind["relic_score_full_gb"] <= thresh, "asv_id"])
        sub_ft = ft.loc[ft.index.isin(alive)]
        sh = pd.DataFrame({"sample": sub_ft.columns,
                              "shannon": [shannon(sub_ft[c].values)
                                            for c in sub_ft.columns]})
        sh = sh.merge(smeta, on="sample", how="left")
        sh = sh.merge(geo_t1[["site", "AnnualMeanTemp"]], on="site",
                        how="left").dropna(subset=["AnnualMeanTemp"])
        r, _ = spearmanr(sh["shannon"], sh["AnnualMeanTemp"])
        rows.append({"alive_threshold": thresh, "n_alive_asvs": len(alive),
                      "median_alive_shannon": float(sh["shannon"].median()),
                      "shannon_vs_MAT_rho": float(r),
                      "n_samples": len(sh)})
        print(f"  alive_threshold<={thresh}: n_asvs={len(alive)}  "
              f"shannon~MAT rho={r:+.3f}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "threshold_sensitivity.tsv", sep="\t", index=False)
    return df


# =============================================================================
# 5. Direct PMA evidence for CSP1-2
# =============================================================================
def csp_pma_evidence():
    print("\n=== [5] Direct PMA evidence for CSP1-2 ===", flush=True)
    csp_fasta = CACHE / "csp1-2_asvs.fasta"
    csp_ids = set()
    with open(csp_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                csp_ids.add(line[1:].strip().split()[0])
    print(f"  CSP1-2 ASVs in EQ: {len(csp_ids)}", flush=True)

    # Per-EQ-ASV viability table
    viab = pd.read_csv(CACHE / "test6_disconfirmation" /
                        "per_eq_asv_viability.tsv", sep="\t")
    print(f"  EQ ASVs with PMA proxy: {len(viab)}", flush=True)
    csp_with_pma = viab[viab["eq_asv"].isin(csp_ids)]
    print(f"  CSP1-2 ASVs with direct PMA proxy: {len(csp_with_pma)}",
          flush=True)
    if len(csp_with_pma) == 0:
        print("  *** NO DIRECT PMA EVIDENCE for any CSP1-2 ASV ***",
              flush=True)
        print("  This means the model classified CSP1-2 entirely by extrap"
                "olation from features.", flush=True)
    else:
        print("  CSP1-2 PMA viability:")
        print(csp_with_pma[["eq_asv", "weighted_median_ratio",
                              "n_pma_proxies"]].round(4).to_string(index=False),
              flush=True)
        # Compare to model's relic_score
        ind = pd.read_csv(CACHE / "test6_disconfirmation" /
                           "relic_indicator_with_damage_per_asv.tsv",
                           sep="\t")
        m = csp_with_pma.merge(ind[["asv_id", "relic_score_full_gb"]],
                                  left_on="eq_asv", right_on="asv_id")
        print("\n  Model score vs PMA truth for CSP1-2:")
        print(m[["eq_asv", "weighted_median_ratio",
                  "relic_score_full_gb"]].round(3)
              .to_string(index=False), flush=True)
    csp_with_pma.to_csv(OUT / "csp_pma_evidence.tsv", sep="\t",
                         index=False)
    return csp_with_pma


def main():
    rows = []
    print(f"All output -> {OUT}", flush=True)
    try:
        site_cv = site_level_cv()
        rows.append({"test": "site_cv", **site_cv.iloc[0].to_dict()})
    except Exception as e:
        print(f"site_level_cv failed: {e}", flush=True)
    try:
        seq = sequence_only_model()
        rows.append({"test": "sequence_only", **seq})
    except Exception as e:
        print(f"sequence_only_model failed: {e}", flush=True)
    try:
        nl = random_subset_null()
        rows.append({"test": "random_null", **nl})
    except Exception as e:
        print(f"random_subset_null failed: {e}", flush=True)
    try:
        ts = threshold_sensitivity()
    except Exception as e:
        print(f"threshold_sensitivity failed: {e}", flush=True)
    try:
        csp = csp_pma_evidence()
    except Exception as e:
        print(f"csp_pma_evidence failed: {e}", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
