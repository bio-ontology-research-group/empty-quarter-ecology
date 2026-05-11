#!/usr/bin/env python3
"""Transition asymmetry analysis: A->B vs B->A in switching cells.

For each switching (site, comp) cell:
  - Order trips temporally (1 -> 5)
  - Identify directional transitions: A->A, A->B, B->A, B->B
  - For each transition, compute environmental delta (delta_precip d7/d30/d365,
    delta_temp, delta_humidity) between consecutive trips
  - Test asymmetry: A->B vs B->A counts; are they symmetric?
  - Test triggers: which delta predicts which direction?
  - Identify hysteresis sites: A in T1, B in T3, A in T5 = oscillation
                              vs A->A->A->B = directional drying

Outputs in cache/transition_asymmetry/.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "transition_asymmetry"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    pst = pd.read_csv(CACHE / "two_strategy_temporal" /
                       "per_sample_strategy_with_precip.tsv", sep="\t")
    pst["dominant"] = np.where(pst["log2_A_over_B"] > 0, "A", "B")
    pst["trip"] = pst["trip"].astype(int)
    pst["site"] = pst["site"].astype(int)

    # Per-(site, comp, trip) median dominant strategy (across replicates)
    cell = (pst.groupby(["site", "compartment", "trip"])
             .agg(median_log2_AoB=("log2_A_over_B", "median"),
                  d7=("d7", "median"),
                  d30=("d30", "median"),
                  d90=("d90", "median"),
                  d180=("d180", "median"),
                  d365=("d365", "median"),
                  n_samples=("sample", "count"))
             .reset_index())
    cell["dominant"] = np.where(cell["median_log2_AoB"] > 0, "A", "B")
    cell.to_csv(OUT / "per_cell_trip_dominant.tsv", sep="\t", index=False)

    # ============================================================
    # 1. Identify ordered trip transitions per (site, comp)
    # ============================================================
    print("=== Transition asymmetry analysis ===\n", flush=True)
    transitions = []
    for (site, comp), g in cell.groupby(["site", "compartment"]):
        if len(g) < 2: continue
        g = g.sort_values("trip").reset_index(drop=True)
        for i in range(len(g) - 1):
            a = g.iloc[i]
            b = g.iloc[i + 1]
            transitions.append({
                "site": site, "compartment": comp,
                "from_trip": int(a["trip"]),
                "to_trip": int(b["trip"]),
                "from_dom": a["dominant"],
                "to_dom": b["dominant"],
                "transition": f"{a['dominant']}->{b['dominant']}",
                "delta_log2_AoB": float(b["median_log2_AoB"] -
                                          a["median_log2_AoB"]),
                "delta_d7": float(b["d7"] - a["d7"]),
                "delta_d30": float(b["d30"] - a["d30"]),
                "delta_d90": float(b["d90"] - a["d90"]),
                "delta_d180": float(b["d180"] - a["d180"]),
                "delta_d365": float(b["d365"] - a["d365"]),
                "to_d7": float(b["d7"]),
                "to_d365": float(b["d365"]),
            })
    tr = pd.DataFrame(transitions)
    tr.to_csv(OUT / "all_transitions.tsv", sep="\t", index=False)
    print(f"Total transitions: {len(tr)}", flush=True)
    print(f"\nTransition counts:")
    print(tr["transition"].value_counts().to_string())

    # ============================================================
    # 2. Asymmetry test: A->B vs B->A
    # ============================================================
    print("\n=== A->B vs B->A asymmetry ===")
    n_AtoB = (tr["transition"] == "A->B").sum()
    n_BtoA = (tr["transition"] == "B->A").sum()
    n_AtoA = (tr["transition"] == "A->A").sum()
    n_BtoB = (tr["transition"] == "B->B").sum()
    print(f"  A->A: {n_AtoA},  A->B: {n_AtoB}", flush=True)
    print(f"  B->A: {n_BtoA},  B->B: {n_BtoB}", flush=True)
    if n_AtoB + n_BtoA > 0:
        asym = (n_AtoB - n_BtoA) / (n_AtoB + n_BtoA)
        # Two-tailed Fisher on contingency:
        #         to_A   to_B
        # from_A  AtoA   AtoB
        # from_B  BtoA   BtoB
        odds, fp = fisher_exact([[n_AtoA, n_AtoB], [n_BtoA, n_BtoB]])
        print(f"  Asymmetry (A->B - B->A) / total switches = {asym:+.3f}",
              flush=True)
        print(f"  Fisher exact (2x2 directional table): "
              f"odds={odds:.3f}, p={fp:.3g}", flush=True)
        # Probability of staying in state
        p_stay_A = n_AtoA / (n_AtoA + n_AtoB) if n_AtoA + n_AtoB > 0 else np.nan
        p_stay_B = n_BtoB / (n_BtoA + n_BtoB) if n_BtoA + n_BtoB > 0 else np.nan
        print(f"  P(stay A | start A) = {p_stay_A:.3f}", flush=True)
        print(f"  P(stay B | start B) = {p_stay_B:.3f}", flush=True)

    # ============================================================
    # 3. Environmental triggers per transition direction
    # ============================================================
    print("\n=== Environmental deltas per transition direction ===")
    for var in ("delta_d7", "delta_d30", "delta_d90", "delta_d180",
                  "delta_d365", "to_d7", "to_d365"):
        AtoB = tr.loc[tr["transition"] == "A->B", var].dropna()
        BtoA = tr.loc[tr["transition"] == "B->A", var].dropna()
        if len(AtoB) < 5 or len(BtoA) < 5: continue
        try:
            U, p_mw = mannwhitneyu(AtoB, BtoA, alternative="two-sided")
        except Exception:
            continue
        print(f"  {var:<15}  A->B median={AtoB.median():+.2f}  "
              f"B->A median={BtoA.median():+.2f}  MW p={p_mw:.3g}",
              flush=True)

    # ============================================================
    # 4. Reversibility: hysteresis vs directional drift
    # ============================================================
    print("\n=== Reversibility patterns ===")
    rec_pat = []
    for (site, comp), g in cell.groupby(["site", "compartment"]):
        if len(g) < 3: continue
        g = g.sort_values("trip")
        seq = "".join(g["dominant"].tolist())
        n_changes = sum(1 for i in range(len(seq) - 1) if seq[i] != seq[i+1])
        rec_pat.append({"site": site, "compartment": comp,
                          "n_trips": len(g), "sequence": seq,
                          "n_changes": n_changes})
    pat = pd.DataFrame(rec_pat)
    pat.to_csv(OUT / "per_cell_sequences.tsv", sep="\t", index=False)
    print(f"  Cells with >=3 trips: {len(pat)}", flush=True)
    print(f"  n_changes distribution:")
    print(pat["n_changes"].value_counts().sort_index().to_string())
    print(f"\n  Top 15 dominance-sequence patterns:")
    print(pat["sequence"].value_counts().head(15).to_string())

    # Categorize sequences
    pat["category"] = "stable"
    pat.loc[pat["sequence"].str.contains("AB|BA"), "category"] = "switching"
    pat.loc[(pat["n_changes"] == 1) &
              (pat["sequence"].str.startswith("A")) &
              (pat["sequence"].str.endswith("B")), "category"] = "A_to_B_drift"
    pat.loc[(pat["n_changes"] == 1) &
              (pat["sequence"].str.startswith("B")) &
              (pat["sequence"].str.endswith("A")), "category"] = "B_to_A_drift"
    pat.loc[pat["n_changes"] >= 2, "category"] = "oscillating"
    pat.loc[~pat["sequence"].str.contains("B"), "category"] = "stable_A"
    pat.loc[~pat["sequence"].str.contains("A"), "category"] = "stable_B"
    print(f"\n  Sequence categories:")
    print(pat["category"].value_counts().to_string())

    # ============================================================
    # 5. Trip-pair-specific transition rates
    # ============================================================
    print("\n=== Per trip-pair transition matrix ===")
    pivot = (tr.groupby(["from_trip", "to_trip", "transition"]).size()
             .unstack(fill_value=0))
    print(pivot.to_string())

    # Per consecutive trip-pair, A->B rate
    trip_pairs = sorted(set(zip(tr["from_trip"], tr["to_trip"])))
    print(f"\n  Per trip-pair A->B / B->A rates:")
    for ft, tt in trip_pairs:
        sub = tr[(tr["from_trip"] == ft) & (tr["to_trip"] == tt)]
        n_AtoB = (sub["transition"] == "A->B").sum()
        n_BtoA = (sub["transition"] == "B->A").sum()
        n_AtoA = (sub["transition"] == "A->A").sum()
        n_BtoB = (sub["transition"] == "B->B").sum()
        if n_AtoA + n_AtoB > 0:
            p_AtoB = n_AtoB / (n_AtoA + n_AtoB)
        else:
            p_AtoB = np.nan
        if n_BtoA + n_BtoB > 0:
            p_BtoA = n_BtoA / (n_BtoA + n_BtoB)
        else:
            p_BtoA = np.nan
        print(f"    T{ft}->T{tt}: P(A->B)={p_AtoB:.3f} (n_A_start="
              f"{n_AtoA + n_AtoB})  P(B->A)={p_BtoA:.3f} (n_B_start="
              f"{n_BtoA + n_BtoB})", flush=True)


if __name__ == "__main__":
    main()
