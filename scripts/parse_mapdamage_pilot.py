#!/usr/bin/env python3
"""Parse mapDamage pilot output (10 EQ metagenomes) into per-position damage
curves. Tests whether C->T (5') and G->A (3') signatures are detectable
above Illumina background.

Inputs (rsynced from unimatrix01 /data/emptyquarter/relic_dna_mapping/
        mapdamage_pilot/<sample>/misincorporation.txt):
  cache/mapdamage_pilot/<sample>/misincorporation.txt
  cache/mapdamage_pilot/<sample>/lgdistribution.txt

Outputs:
  cache/mapdamage_pilot/per_sample_damage_curve.tsv
  cache/mapdamage_pilot/per_sample_summary.tsv
  cache/mapdamage_pilot/pilot_verdict.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
PILOT = REPO / "cache" / "mapdamage_pilot"

SAMPLE_INFO = {
    "1Dr2":   ("D",  "site1"),
    "14Dr3":  ("D",  "site14"),
    "16Dr1":  ("D",  "site16"),
    "19Dr1":  ("D",  "site19"),
    "1Sr3":   ("S",  "site1"),
    "12Sr2":  ("S",  "site12"),
    "15Sr1":  ("S",  "site15"),
    "1PRr2":  ("PR", "site1"),
    "15PRr3": ("PR", "site15"),
    "18PRr1": ("PR", "site18"),
}

MAX_POS = 25  # positions from end


def parse_misincorp(path: Path) -> pd.DataFrame:
    """Return DataFrame: end (5p/3p), std (+/-), pos, A,C,G,T,Total + subs."""
    df = pd.read_csv(path, sep="\t", comment="#")
    return df


def damage_curve(df: pd.DataFrame, sample: str) -> pd.DataFrame:
    """Compute C->T per position from 5' end, G->A from 3' end.
    Sum counts across strand within each (end, pos)."""
    rows = []
    # 5' end C->T
    for pos in range(1, MAX_POS + 1):
        sub = df[(df["End"] == "5p") & (df["Pos"] == pos)]
        if len(sub) == 0: continue
        # collapse strands
        c_total = float(sub["C"].sum())
        ct_count = float(sub["C>T"].sum())
        ct_freq = ct_count / c_total if c_total > 0 else np.nan
        # also G>A as control (should NOT show damage at 5')
        g_total = float(sub["G"].sum())
        ga_count = float(sub["G>A"].sum())
        ga_freq = ga_count / g_total if g_total > 0 else np.nan
        # transversions for noise reference
        a_total = float(sub["A"].sum())
        ag_count = float(sub["A>G"].sum())
        ag_freq = ag_count / a_total if a_total > 0 else np.nan
        rows.append({"sample": sample, "end": "5p", "pos": pos,
                      "ct_freq": ct_freq, "ga_freq": ga_freq,
                      "ag_freq": ag_freq,
                      "n_C": int(c_total), "n_G": int(g_total)})
    # 3' end G->A
    for pos in range(1, MAX_POS + 1):
        sub = df[(df["End"] == "3p") & (df["Pos"] == pos)]
        if len(sub) == 0: continue
        c_total = float(sub["C"].sum())
        ct_count = float(sub["C>T"].sum())
        ct_freq = ct_count / c_total if c_total > 0 else np.nan
        g_total = float(sub["G"].sum())
        ga_count = float(sub["G>A"].sum())
        ga_freq = ga_count / g_total if g_total > 0 else np.nan
        a_total = float(sub["A"].sum())
        ag_count = float(sub["A>G"].sum())
        ag_freq = ag_count / a_total if a_total > 0 else np.nan
        rows.append({"sample": sample, "end": "3p", "pos": pos,
                      "ct_freq": ct_freq, "ga_freq": ga_freq,
                      "ag_freq": ag_freq,
                      "n_C": int(c_total), "n_G": int(g_total)})
    return pd.DataFrame(rows)


def parse_length(path: Path):
    if not path.exists(): return None
    df = pd.read_csv(path, sep="\t", comment="#")
    if len(df) == 0 or "Occurences" not in df.columns: return None
    df = df.groupby("Length")["Occurences"].sum()
    return df


def main():
    all_curves = []
    summary_rows = []
    for sample, (comp, site) in SAMPLE_INFO.items():
        mfile = PILOT / sample / "misincorporation.txt"
        if not mfile.exists():
            print(f"  [skip] {sample}: missing", flush=True)
            continue
        misincorp = parse_misincorp(mfile)
        curve = damage_curve(misincorp, sample)
        curve["compartment"] = comp
        curve["site"] = site
        all_curves.append(curve)

        # Summary stats
        ct_pos1 = curve.loc[(curve["end"] == "5p") & (curve["pos"] == 1),
                              "ct_freq"].values
        ct_pos2 = curve.loc[(curve["end"] == "5p") & (curve["pos"] == 2),
                              "ct_freq"].values
        ct_pos5 = curve.loc[(curve["end"] == "5p") & (curve["pos"] == 5),
                              "ct_freq"].values
        ga_pos1 = curve.loc[(curve["end"] == "3p") & (curve["pos"] == 1),
                              "ga_freq"].values
        # control: A->G transition at internal positions (modern Illumina ~baseline)
        ag_internal = curve.loc[(curve["pos"] >= 10) &
                                 (curve["pos"] <= 20), "ag_freq"]
        bg_ag = float(ag_internal.mean())
        # control: 5' G->A and 3' C->T (should NOT show damage)
        ga_5p_pos1 = curve.loc[(curve["end"] == "5p") & (curve["pos"] == 1),
                                 "ga_freq"].values
        ct_3p_pos1 = curve.loc[(curve["end"] == "3p") & (curve["pos"] == 1),
                                 "ct_freq"].values

        # length distribution
        ldist = parse_length(PILOT / sample / "lgdistribution.txt")
        median_len = float(np.average(ldist.index,
                                         weights=ldist.values)) \
            if ldist is not None else np.nan
        # weighted percentiles
        if ldist is not None:
            cum = ldist.cumsum() / ldist.sum()
            p25_l = float(ldist.index[(cum >= 0.25).idxmax()])
            p75_l = float(ldist.index[(cum >= 0.75).idxmax()])
        else:
            p25_l = p75_l = np.nan

        summary_rows.append({
            "sample": sample, "compartment": comp, "site": site,
            "ct_5p_pos1": float(ct_pos1[0]) if len(ct_pos1) else np.nan,
            "ct_5p_pos2": float(ct_pos2[0]) if len(ct_pos2) else np.nan,
            "ct_5p_pos5": float(ct_pos5[0]) if len(ct_pos5) else np.nan,
            "ga_3p_pos1": float(ga_pos1[0]) if len(ga_pos1) else np.nan,
            "ag_internal_bg": bg_ag,
            "ga_5p_pos1_ctrl": float(ga_5p_pos1[0]) if len(ga_5p_pos1)
                                 else np.nan,
            "ct_3p_pos1_ctrl": float(ct_3p_pos1[0]) if len(ct_3p_pos1)
                                 else np.nan,
            "wt_mean_readlen": median_len,
            "readlen_p25": p25_l, "readlen_p75": p75_l,
        })

    if not all_curves:
        print("No samples found. Run pilot first.", flush=True)
        sys.exit(1)

    curves = pd.concat(all_curves, ignore_index=True)
    curves.to_csv(PILOT / "per_sample_damage_curve.tsv",
                   sep="\t", index=False)
    summ = pd.DataFrame(summary_rows)
    summ.to_csv(PILOT / "per_sample_summary.tsv", sep="\t", index=False)

    print("\n=== Per-sample damage summary ===")
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 240)
    print(summ.round(5).to_string(index=False))

    # Excess C->T at pos1 above background AG
    summ["ct_excess"] = summ["ct_5p_pos1"] - summ["ag_internal_bg"]
    print("\n=== C->T excess (pos1 5p) over A->G internal background ===")
    print(summ[["sample", "compartment", "ct_5p_pos1",
                  "ag_internal_bg", "ct_excess"]].round(5).to_string(index=False))

    # Compartment comparison
    print("\n=== Compartment medians (5' C->T pos 1) ===")
    print(summ.groupby("compartment")["ct_5p_pos1"].agg(["count", "mean",
                                                            "median", "min",
                                                            "max"])
          .round(5).to_string())

    # Verdict
    print("\n=== VERDICT ===")
    median_ct = summ["ct_5p_pos1"].median()
    median_ctrl = summ["ga_5p_pos1_ctrl"].median()
    median_ratio = (median_ct / median_ctrl
                    if median_ctrl > 0 else float("inf"))
    print(f"  median 5' C->T pos1: {median_ct:.4f}")
    print(f"  median 5' G->A pos1 (negative ctrl): {median_ctrl:.4f}")
    print(f"  C->T : G->A ratio at 5' pos 1: {median_ratio:.2f}")
    surf_ct = summ.loc[summ["compartment"] == "S", "ct_5p_pos1"].median()
    deep_ct = summ.loc[summ["compartment"] == "D", "ct_5p_pos1"].median()
    print(f"  surface median 5' C->T: {surf_ct:.4f}")
    print(f"  deep    median 5' C->T: {deep_ct:.4f}")

    with open(PILOT / "pilot_verdict.txt", "w") as fh:
        fh.write("mapDamage pilot verdict (10 EQ metagenomes)\n")
        fh.write("=" * 70 + "\n\n")
        fh.write("Per-sample summary:\n")
        fh.write(summ.round(5).to_string(index=False))
        fh.write("\n\nC->T excess over background:\n")
        fh.write(summ[["sample", "compartment", "ct_5p_pos1",
                        "ag_internal_bg", "ct_excess"]].round(5)
                  .to_string(index=False))
        fh.write("\n\nCompartment medians (5' C->T pos1):\n")
        fh.write(summ.groupby("compartment")["ct_5p_pos1"]
                  .agg(["count", "mean", "median", "min", "max"]).round(5)
                  .to_string())
        fh.write(f"\n\n--- Diagnostic ---\n")
        fh.write(f"  median 5' C->T pos1: {median_ct:.4f}\n")
        fh.write(f"  median 5' G->A pos1 (neg ctrl): {median_ctrl:.4f}\n")
        fh.write(f"  C->T:G->A ratio at 5' pos1: {median_ratio:.2f}\n")
        fh.write(f"  surf median: {surf_ct:.4f}\n")
        fh.write(f"  deep median: {deep_ct:.4f}\n\n")
        fh.write("--- Interpretation key ---\n")
        fh.write("  aDNA threshold: 5' C->T pos1 > 5% (0.05) is classic aDNA\n")
        fh.write("  Modern Illumina: 5' C->T pos1 ~ 0.1-0.5% (0.001-0.005)\n")
        fh.write("  Intermediate (1-5%): degraded / soil-aged DNA\n")
        fh.write("  Asymmetry: C->T > G->A at 5' AND G->A > C->T at 3'\n"
                  "    indicates true damage signal\n")
        fh.write("  Compartment: surface > deep would support UV-driven\n"
                  "    deamination\n")


if __name__ == "__main__":
    main()
