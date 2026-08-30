#!/usr/bin/env python3
"""Compare the rerun outputs with the published screen (analysis/v3/control_audit).

Prints a Markdown report; run_rerun.sh saves it as outputs/comparison.md.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PUB = ROOT / "analysis/v3/control_audit"
OUT = HERE / "outputs"


def tsv(path):
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def genus(tax):
    p = tax.split(";")
    return p[5].strip() if len(p) > 5 and p[5].strip() else "(unassigned genus)"


def stats(fractions):
    a = np.asarray(fractions, dtype=float) * 100
    if not len(a):
        return "n/a"
    q1, med, q3 = np.percentile(a, [25, 50, 75])
    return f"median {med:.2f}% (IQR {q1:.2f}-{q3:.2f}%), max {a.max():.2f}%"


def main():
    pub_calls = tsv(PUB / "trip5_primary_contaminant_calls.tsv")
    pub_prof = tsv(PUB / "trip5_removal_fraction_by_profile.tsv")
    pub_sum = json.loads((PUB / "summary.json").read_text())
    new_calls = tsv(OUT / "primary_contaminant_calls.tsv")
    new_prof = tsv(OUT / "removal_fraction_by_profile.tsv")
    new_sum = json.loads((OUT / "summary.json").read_text())
    pub_set = {r["feature_id"] for r in pub_calls}
    by_screen = {}
    for r in new_calls:
        by_screen.setdefault(r["screen"], set()).add(r["feature_id"])
    print(f"# Contaminant screen rerun vs published (generated {new_sum['generated_at']})\n")
    print(f"Mode: {new_sum['mode']}; Trip 4 screen: {new_sum['trip4_screen']}\n")
    print("## Candidate ASVs\n")
    print(f"- Published (Trip 5 blanks EB1-EB17 vs 217 profiles): {len(pub_set)}")
    for s, ids in by_screen.items():
        print(f"- Rerun screen {s}: {len(ids)} candidates; overlap with published set {len(ids & pub_set)}; "
              f"new {len(ids - pub_set)}; dropped {len(pub_set - ids)}")
    t5 = by_screen.get("Trip5")
    if t5 is not None:
        print("- Trip 5 screen reproduces the published candidate set exactly: "
              f"{'YES' if t5 == pub_set else 'NO'}")
    t4 = by_screen.get("Trip4") or by_screen.get("pooled")
    if t4:
        c = Counter(genus(r["taxonomy"]) for r in new_calls if r["screen"] in ("Trip4", "pooled"))
        print("- Top genera among Trip 4 / pooled candidates: " + ", ".join(f"{g} ({n})" for g, n in c.most_common(12)))
        c2 = Counter(genus(r["taxonomy"]) for r in new_calls if r["screen"] in ("Trip4", "pooled") and r["feature_id"] not in pub_set)
        print("- Top genera among candidates NOT in the published set: " + ", ".join(f"{g} ({n})" for g, n in c2.most_common(12)))
    print("\n## Profiles\n")
    pub_bio = [r for r in pub_prof if r["role"] == "compatible_biological_profile"]
    new_bio = [r for r in new_prof if r["role"] == "compatible_biological_profile"]
    print(f"- Published: {len(pub_bio)} filtered profiles; reads removed "
          f"{stats([float(r['candidate_contaminant_read_fraction']) for r in pub_bio])}; "
          f"pooled {100 * sum(int(r['candidate_contaminant_reads']) for r in pub_bio) / sum(int(r['total_reads']) for r in pub_bio):.2f}%")
    for trip in ("Trip5", "Trip4"):
        rows = [r for r in new_bio if r["trip"] == trip]
        if not rows:
            print(f"- Rerun {trip}: no filtered profiles (screen not run)")
            continue
        fr = [float(r["candidate_contaminant_read_fraction"]) for r in rows]
        changed = sum(1 for r in rows if int(r["candidate_contaminant_reads"]) > 0)
        below = sum(1 for r in rows if int(r["total_reads"]) - int(r["candidate_contaminant_reads"]) < 25000)
        print(f"- Rerun {trip}: {len(rows)} filtered profiles, {changed} lose >= 1 read, reads removed {stats(fr)}; "
              f"pooled {100 * sum(int(r['candidate_contaminant_reads']) for r in rows) / sum(int(r['total_reads']) for r in rows):.2f}%; "
              f"profiles below 25,000 reads after filtering: {below}")
    # Per-profile differences for Trip 5 (should be zero in per-trip mode)
    pub_by = {r["profile_id"]: r for r in pub_bio}
    diffs = [(r["profile_id"], r["candidate_contaminant_reads"], pub_by[r["profile_id"]]["candidate_contaminant_reads"])
             for r in new_bio if r["profile_id"] in pub_by
             and r["candidate_contaminant_reads"] != pub_by[r["profile_id"]]["candidate_contaminant_reads"]]
    print(f"- Trip 5 profiles whose removed-read count differs from the published run: {len(diffs)}"
          + (f" (first: {diffs[:5]})" if diffs else ""))
    print("\n## Headline numbers the manuscript reports (Methods, 'Assay controls')\n")
    print("| quantity | manuscript / published | rerun |")
    print("|---|---|---|")
    print(f"| candidate ASVs | 351 | {new_sum['candidates']} |")
    print(f"| linked profiles filtered | 217 (Trip 5) | Trip 5 {new_sum['trip5_mapped_profiles']}, Trip 4 {new_sum['trip4_mapped_profiles'] if new_sum['trip4_screen'] == 'run' else 'not filtered'} |")
    print(f"| training blanks | 17 | Trip 5 {len(new_sum['trip5_training_blanks'])}, Trip 4 {len(new_sum['trip4_training_blanks'])} |")
    print("\n## Downstream (25 tracked conclusions)\n")
    print("Not rerun by this script. The published downstream stage "
          "(scripts/controls/build_control_sensitivity_inputs.py + run_control_ecology_sensitivity.sh) "
          "hard-codes 217 Trip 5 profiles and V-prefixed IDs; extending it to Trip 4 needs a patched copy. "
          "Compare analysis/v3/control_sensitivity/headline_result_sensitivity.tsv once that is run.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
