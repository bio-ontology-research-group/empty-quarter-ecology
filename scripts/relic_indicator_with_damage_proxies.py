#!/usr/bin/env python3
"""Track C: extend the composite relic-likelihood indicator with
amplicon-level damage proxies.

Adds these per-ASV features:
  - cluster_fanout:        # of ASVs collapsing to this ASV's 99% OTU
  - is_singleton_otu:      0/1 = ASV alone at OTU centroid
  - gc_content:            G+C / sequence length
  - length_dev:            |asv_len - median_len|
  - pyr_dinuc_density:     fraction of dinucleotides that are TpC, CpC,
                              TpT, or CpT (UV-vulnerable sites)
  - tc_cc_density:         fraction TpC+CpC (the highest UV-deamination
                              acceptors at CPDs)
  - intra_otu_min_pid:     1 - min sim to other ASVs in same OTU (proxy
                              for sequencing-error inflation)

Refits LR + GB models, compares AUC to baseline.

Inputs:
  data/taxonomy/ASV_seqs-trips1-5.fasta
  cache/test6_disconfirmation/asv_to_otu_99.tsv
  cache/test6_disconfirmation/relic_features.parquet
  cache/test6_disconfirmation/per_eq_asv_viability.tsv

Outputs:
  cache/test6_disconfirmation/relic_features_with_damage.parquet
  cache/test6_disconfirmation/relic_indicator_with_damage_per_asv.tsv
  cache/test6_disconfirmation/relic_model_metrics_with_damage.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import (roc_auc_score, average_precision_score,
                                brier_score_loss, log_loss)
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
OUT = CACHE / "test6_disconfirmation"

RELIC_T_LO = 0.1
ALIVE_T_HI = 0.5


def parse_fasta(path: Path) -> pd.DataFrame:
    rows = []
    cur_id, cur_seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if cur_id is not None:
                    rows.append((cur_id, "".join(cur_seq)))
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.upper())
        if cur_id is not None:
            rows.append((cur_id, "".join(cur_seq)))
    return pd.DataFrame(rows, columns=["asv_id", "seq"])


def compute_seq_features(seqs: pd.DataFrame) -> pd.DataFrame:
    print("Computing per-ASV sequence features ...", flush=True)
    s = seqs["seq"]
    L = s.str.len()
    median_len = float(L.median())
    print(f"  ASVs: {len(seqs)}; median len: {median_len:.0f}", flush=True)
    gc = (s.str.count("G") + s.str.count("C")) / L
    length_dev = (L - median_len).abs()

    # Pyrimidine dinucleotide densities -- UV-vulnerable
    # Count TT, TC, CT, CC; divide by L-1 (n dinucleotides)
    pyr_count = (s.str.count("TT") + s.str.count("TC") +
                  s.str.count("CT") + s.str.count("CC"))
    tc_cc_count = s.str.count("TC") + s.str.count("CC")
    pyr_dinuc_density = pyr_count / (L - 1).clip(lower=1)
    tc_cc_density = tc_cc_count / (L - 1).clip(lower=1)

    return pd.DataFrame({
        "asv_id": seqs["asv_id"],
        "asv_len": L,
        "gc_content": gc,
        "length_dev": length_dev,
        "pyr_dinuc_density": pyr_dinuc_density,
        "tc_cc_density": tc_cc_density,
    })


def compute_otu_fanout(otu_path: Path,
                          seq_features: pd.DataFrame,
                          seqs: pd.DataFrame) -> pd.DataFrame:
    """For each ASV: # of ASVs collapsing to its 99% OTU; min within-OTU
    nucleotide identity to another ASV (proxy for cluster heterogeneity)."""
    print("Computing OTU fan-out + intra-OTU identity ...", flush=True)
    otu_map = pd.read_csv(otu_path, sep="\t", header=None,
                           names=["asv_id", "otu_id"])
    fanout = otu_map.groupby("otu_id").size().rename("cluster_fanout")
    otu_map = otu_map.merge(fanout, on="otu_id", how="left")
    otu_map["is_singleton_otu"] = (otu_map["cluster_fanout"] == 1).astype(int)
    print(f"  total ASV->OTU mappings: {len(otu_map)}", flush=True)
    print(f"  fanout p50/p75/p90/max: "
          f"{int(np.percentile(otu_map['cluster_fanout'], 50))}/"
          f"{int(np.percentile(otu_map['cluster_fanout'], 75))}/"
          f"{int(np.percentile(otu_map['cluster_fanout'], 90))}/"
          f"{otu_map['cluster_fanout'].max()}", flush=True)

    # Intra-OTU nearest-neighbor (Hamming approximation):
    # Within each multi-ASV OTU, compute mean fraction-mismatch among ASVs.
    # Using simple pairwise comparison on prefix-aligned (centered to median
    # length); for performance, sample up to 5 ASVs per OTU.
    print("  computing within-OTU identity (sampled) ...", flush=True)
    seq_lookup = dict(zip(seqs["asv_id"], seqs["seq"]))
    rows = []
    multi = otu_map[otu_map["cluster_fanout"] > 1].groupby("otu_id")
    for otu_id, g in multi:
        members = list(g["asv_id"])
        if len(members) > 5:
            members = list(np.random.RandomState(int(otu_id[-6:], 16)
                                                  if otu_id[-6:].isalnum()
                                                  else 42)
                             .choice(members, 5, replace=False))
        seqs_m = [seq_lookup.get(m, "") for m in members]
        # min pairwise identity for each member
        for i, ai in enumerate(g["asv_id"]):
            si = seq_lookup.get(ai, "")
            if not si: continue
            min_pid = 1.0
            for j, sj in enumerate(seqs_m):
                if not sj or sj == si: continue
                # Compare on min length, simple Hamming
                Lm = min(len(si), len(sj))
                if Lm == 0: continue
                mm = sum(1 for k in range(Lm) if si[k] != sj[k])
                pid = 1 - mm / Lm
                min_pid = min(min_pid, pid)
            if min_pid < 1.0:
                rows.append({"asv_id": ai, "intra_otu_min_pid": min_pid})
    intra_pid = pd.DataFrame(rows)
    print(f"  intra-OTU pid records: {len(intra_pid)}", flush=True)

    out = otu_map[["asv_id", "cluster_fanout", "is_singleton_otu"]].copy()
    out = out.merge(intra_pid, on="asv_id", how="left")
    # Singletons: pid = 1.0 (no neighbor)
    out["intra_otu_min_pid"] = out["intra_otu_min_pid"].fillna(1.0)
    return out


def main():
    seqs = parse_fasta(REPO / "data" / "taxonomy" / "ASV_seqs-trips1-5.fasta")
    seq_feats = compute_seq_features(seqs)
    otu_feats = compute_otu_fanout(OUT / "asv_to_otu_99.tsv", seq_feats, seqs)

    print("Loading baseline relic feature table ...", flush=True)
    base = pd.read_parquet(OUT / "relic_features.parquet")
    print(f"  baseline shape: {base.shape}", flush=True)

    feats = (base
             .merge(seq_feats, on="asv_id", how="left")
             .merge(otu_feats, on="asv_id", how="left"))
    print(f"  merged shape: {feats.shape}", flush=True)
    # Defaults if missing
    feats["cluster_fanout"] = feats["cluster_fanout"].fillna(1).astype(int)
    feats["is_singleton_otu"] = feats["is_singleton_otu"].fillna(1).astype(int)
    feats["intra_otu_min_pid"] = feats["intra_otu_min_pid"].fillna(1.0)
    feats.to_parquet(OUT / "relic_features_with_damage.parquet")

    # Train set
    train_mask = (feats["weighted_median_ratio"].notna() &
                   ((feats["weighted_median_ratio"] < RELIC_T_LO) |
                    (feats["weighted_median_ratio"] > ALIVE_T_HI)))
    df = feats[train_mask].copy()
    df["y_relic"] = (df["weighted_median_ratio"] < RELIC_T_LO).astype(int)
    print(f"\nTraining set: {len(df)} ASVs "
          f"(y=1 relic: {int(df['y_relic'].sum())}, "
          f"y=0 alive: {int((1-df['y_relic']).sum())})", flush=True)

    base_cols = [
        "persistence_max", "persistence_mean", "n_site_records",
        "log_mean_abund", "log_max_abund", "n_detections", "n_sites_detected",
        "frac_deep", "frac_surface", "frac_rhizosphere",
        "emp_cosmo_90", "emp_cosmo_100", "emp_cosmo_150",
    ]
    damage_cols = [
        "asv_len", "gc_content", "length_dev",
        "pyr_dinuc_density", "tc_cc_density",
        "cluster_fanout", "is_singleton_otu", "intra_otu_min_pid",
    ]
    all_cols = base_cols + damage_cols

    df = df.dropna(subset=all_cols).copy()
    print(f"  after dropna: {len(df)}", flush=True)
    y = df["y_relic"].values

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    def fit_and_score(cols, label):
        X = df[cols].values
        X = np.where(np.isfinite(X), X, 0.0)
        sc = StandardScaler().fit(X)
        Xs = sc.transform(X)
        lr = LogisticRegression(penalty="l2", C=1.0, max_iter=2000,
                                  class_weight="balanced", random_state=42)
        gb = GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                            learning_rate=0.05,
                                            random_state=42)
        p_lr = cross_val_predict(lr, Xs, y, cv=cv,
                                    method="predict_proba")[:, 1]
        p_gb = cross_val_predict(gb, X, y, cv=cv,
                                    method="predict_proba")[:, 1]
        out = {
            "label": label,
            "n_features": len(cols),
            "auc_lr": roc_auc_score(y, p_lr),
            "ap_lr": average_precision_score(y, p_lr),
            "brier_lr": brier_score_loss(y, p_lr),
            "auc_gb": roc_auc_score(y, p_gb),
            "ap_gb": average_precision_score(y, p_gb),
            "brier_gb": brier_score_loss(y, p_gb),
        }
        # Refit final
        lr.fit(Xs, y); gb.fit(X, y)
        return out, lr, gb, sc

    print("\n=== Comparing models ===", flush=True)
    res_base, _, _, _ = fit_and_score(base_cols, "baseline")
    res_dmg, _, _, _ = fit_and_score(damage_cols, "damage_only")
    res_full, lr_f, gb_f, sc_f = fit_and_score(all_cols, "full")

    for r in (res_base, res_dmg, res_full):
        print(f"\n  {r['label']}  (n_features={r['n_features']})")
        print(f"    LR: AUC={r['auc_lr']:.3f} AP={r['ap_lr']:.3f} "
              f"Brier={r['brier_lr']:.3f}")
        print(f"    GB: AUC={r['auc_gb']:.3f} AP={r['ap_gb']:.3f} "
              f"Brier={r['brier_gb']:.3f}")

    delta_auc_lr = res_full['auc_lr'] - res_base['auc_lr']
    delta_auc_gb = res_full['auc_gb'] - res_base['auc_gb']
    print(f"\nDelta AUC from adding damage proxies:", flush=True)
    print(f"  LR: {delta_auc_lr:+.3f}", flush=True)
    print(f"  GB: {delta_auc_gb:+.3f}", flush=True)

    # Feature importance from full GB
    fi = pd.DataFrame({"feature": all_cols,
                          "lr_coef": lr_f.coef_[0],
                          "gb_importance": gb_f.feature_importances_,
                          }).sort_values("gb_importance", ascending=False)
    print("\nFeature importance (full model):")
    print(fi.to_string(index=False), flush=True)

    # Score every ASV with full model
    full_all = feats.dropna(subset=all_cols).copy()
    Xall = full_all[all_cols].values
    Xall = np.where(np.isfinite(Xall), Xall, 0.0)
    full_all["relic_score_full_gb"] = gb_f.predict_proba(Xall)[:, 1]
    Xall_s = sc_f.transform(Xall)
    full_all["relic_score_full_lr"] = lr_f.predict_proba(Xall_s)[:, 1]

    out_cols = ["asv_id"] + all_cols + [
        "weighted_median_ratio", "n_pma_proxies",
        "relic_score_full_lr", "relic_score_full_gb",
    ]
    full_all[out_cols].to_csv(
        OUT / "relic_indicator_with_damage_per_asv.tsv",
        sep="\t", index=False)

    # Write summary
    with open(OUT / "relic_model_metrics_with_damage.txt", "w") as fh:
        fh.write("Composite relic-likelihood indicator + amplicon damage "
                  "proxies (Track C)\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Training set: {len(df)} ASVs (y=1 relic: "
                  f"{int(df['y_relic'].sum())}, y=0 alive: "
                  f"{int((1-df['y_relic']).sum())})\n")
        fh.write(f"Application:  {len(full_all)} ASVs scored\n\n")

        fh.write("--- Cross-validated comparison ---\n\n")
        for r in (res_base, res_dmg, res_full):
            fh.write(f"{r['label']} (n_features={r['n_features']}):\n")
            fh.write(f"  LR: AUC={r['auc_lr']:.3f} AP={r['ap_lr']:.3f} "
                      f"Brier={r['brier_lr']:.3f}\n")
            fh.write(f"  GB: AUC={r['auc_gb']:.3f} AP={r['ap_gb']:.3f} "
                      f"Brier={r['brier_gb']:.3f}\n\n")
        fh.write(f"Delta AUC (full - baseline):\n")
        fh.write(f"  LR: {delta_auc_lr:+.3f}\n")
        fh.write(f"  GB: {delta_auc_gb:+.3f}\n\n")

        fh.write("--- Feature importance (full model) ---\n")
        fh.write(fi.to_string(index=False))
        fh.write("\n\n--- Damage-only feature importance "
                  "(testing if damage signal exists at all) ---\n")

        # Quick subset just damage features
        Xd = df[damage_cols].values
        Xd = np.where(np.isfinite(Xd), Xd, 0.0)
        gb_d = GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                              learning_rate=0.05,
                                              random_state=42)
        gb_d.fit(Xd, y)
        fid = pd.DataFrame({"feature": damage_cols,
                              "gb_importance": gb_d.feature_importances_,
                              }).sort_values("gb_importance", ascending=False)
        fh.write(fid.to_string(index=False))
        fh.write("\n")

    print(f"\nWrote {OUT}/relic_model_metrics_with_damage.txt", flush=True)


if __name__ == "__main__":
    main()
