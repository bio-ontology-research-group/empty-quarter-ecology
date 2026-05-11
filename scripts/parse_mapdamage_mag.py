#!/usr/bin/env python3
"""Parse mapDamage MAG remap output (3 samples mapped to clean MAG references).
Compares to the assembly-mapped pilot to determine if the negative damage
result was methodology (noisy contigs) or biology (no damage to detect).
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
PILOT = REPO / "cache" / "mapdamage_pilot"
MAG = REPO / "cache" / "mapdamage_mag"

SAMPLES = {"1Dr2": "D", "1Sr3": "S", "1PRr2": "PR"}
MAX_POS = 25


def parse_misincorp(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#")


def per_position(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for end in ("5p", "3p"):
        for pos in range(1, MAX_POS + 1):
            sub = df[(df["End"] == end) & (df["Pos"] == pos)]
            if len(sub) == 0: continue
            c = float(sub["C"].sum())
            g = float(sub["G"].sum())
            a = float(sub["A"].sum())
            ct = float(sub["C>T"].sum())
            ga = float(sub["G>A"].sum())
            ag = float(sub["A>G"].sum())
            rows.append({"end": end, "pos": pos,
                          "ct_freq": ct / c if c > 0 else np.nan,
                          "ga_freq": ga / g if g > 0 else np.nan,
                          "ag_freq": ag / a if a > 0 else np.nan,
                          "n_C": int(c), "n_G": int(g)})
    return pd.DataFrame(rows)


def main():
    rows = []
    for sample, comp in SAMPLES.items():
        for label, base in (("assembly", PILOT), ("MAG", MAG)):
            path = base / sample / "misincorporation.txt"
            if not path.exists():
                continue
            df = parse_misincorp(path)
            cur = per_position(df)
            cur["sample"] = sample
            cur["compartment"] = comp
            cur["mapping"] = label
            rows.append(cur)
    if not rows:
        sys.exit("nothing to parse")
    all_df = pd.concat(rows, ignore_index=True)
    all_df.to_csv(MAG / "assembly_vs_mag_damage_curve.tsv",
                   sep="\t", index=False)

    # Compare position 1 5' end and 3' end
    print("=== Position 1 (key aDNA position) ===")
    print(f"{'sample':<8} {'comp':<4} {'mapping':<10} "
          f"{'5p_CT':>8} {'5p_GA':>8} {'3p_CT':>8} {'3p_GA':>8}")
    for sample in SAMPLES:
        for label in ("assembly", "MAG"):
            sub = all_df[(all_df["sample"] == sample) &
                          (all_df["mapping"] == label) &
                          (all_df["pos"] == 1)]
            if len(sub) == 0: continue
            ct5 = float(sub.loc[sub["end"] == "5p", "ct_freq"].iloc[0])
            ga5 = float(sub.loc[sub["end"] == "5p", "ga_freq"].iloc[0])
            ct3 = float(sub.loc[sub["end"] == "3p", "ct_freq"].iloc[0])
            ga3 = float(sub.loc[sub["end"] == "3p", "ga_freq"].iloc[0])
            print(f"{sample:<8} {SAMPLES[sample]:<4} {label:<10} "
                  f"{ct5:>8.4f} {ga5:>8.4f} {ct3:>8.4f} {ga3:>8.4f}")

    # Decay curves (5' C->T pos 1-10, averaged across the 3 samples per mapping)
    print("\n=== Mean 5' C->T frequency per position (3 samples) ===")
    for label in ("assembly", "MAG"):
        sub = all_df[(all_df["mapping"] == label) & (all_df["end"] == "5p") &
                      (all_df["pos"] <= 10)]
        if len(sub) == 0: continue
        avg = sub.groupby("pos")["ct_freq"].mean().round(4)
        print(f"  {label:<10}: {avg.to_dict()}")

    print("\n=== Internal A->G background by mapping (avg pos 10-20) ===")
    for label in ("assembly", "MAG"):
        sub = all_df[(all_df["mapping"] == label) &
                      (all_df["pos"].between(10, 20))]
        avg = sub["ag_freq"].mean()
        print(f"  {label:<10}: {avg:.4f}")

    # Verdict
    print("\n=== VERDICT ===")
    p1_ct_mag = all_df[(all_df["mapping"] == "MAG") &
                         (all_df["end"] == "5p") & (all_df["pos"] == 1)]["ct_freq"].mean()
    p1_ct_asm = all_df[(all_df["mapping"] == "assembly") &
                         (all_df["end"] == "5p") & (all_df["pos"] == 1)]["ct_freq"].mean()
    print(f"  5' C->T pos1 (MAG):      {p1_ct_mag:.4f}")
    print(f"  5' C->T pos1 (assembly): {p1_ct_asm:.4f}")
    print(f"  Delta (MAG - assembly):  {(p1_ct_mag - p1_ct_asm):+.4f}")
    if p1_ct_mag >= 0.05:
        print("  MAG remap REVEALS aDNA damage -- methodology issue with assembly")
    elif p1_ct_mag > p1_ct_asm * 2:
        print("  MAG remap shows MODESTLY elevated damage vs assembly")
    else:
        print("  MAG remap CONFIRMS no aDNA damage signal -- biology, "
                "not methodology")

    with open(MAG / "verdict.txt", "w") as fh:
        fh.write("mapDamage MAG remap robustness check\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"  5' C->T pos1 (MAG mapping):      {p1_ct_mag:.4f}\n")
        fh.write(f"  5' C->T pos1 (assembly mapping): {p1_ct_asm:.4f}\n\n")
        fh.write("Per-sample 5' / 3' position 1:\n")
        for sample in SAMPLES:
            for label in ("assembly", "MAG"):
                sub = all_df[(all_df["sample"] == sample) &
                              (all_df["mapping"] == label) &
                              (all_df["pos"] == 1)]
                if len(sub) == 0: continue
                ct5 = float(sub.loc[sub["end"] == "5p", "ct_freq"].iloc[0])
                ga5 = float(sub.loc[sub["end"] == "5p", "ga_freq"].iloc[0])
                fh.write(f"  {sample:<8} {SAMPLES[sample]:<4} {label:<10}  "
                          f"5pCT={ct5:.4f}  5pGA={ga5:.4f}\n")


if __name__ == "__main__":
    main()
