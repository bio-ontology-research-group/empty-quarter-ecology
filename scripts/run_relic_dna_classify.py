#!/usr/bin/env python3
"""Tier-1 #3 step 3: Quantify relic DNA in EQ 16S signal.

Strategy:
  - Per ASV, count samples it appears in (16S detection set S_16S).
  - For each sample, the MAG community is the SemiBin bins of that sample.
  - For each ASV with a MAG match (asv_to_mag.tsv), find the samples whose
    MAG community contains the matched bin (MAG-supported set S_MAG_supp).
  - Classify each ASV-sample observation:
        active-consistent: 16S detected AND MAG of matched bin assembled
        intermittent:      16S detected, matched MAG present in OTHER samples
        MAG-orphan:        16S detected, ASV has NO MAG match anywhere
                           (likely relic OR low-abundance true cell)
  - Aggregate to phylum / compartment / season; report fraction of total
    16S read mass classified as MAG-orphan.

Inputs:
  cache/feature_table.parquet              ASV x sample
  cache/taxonomy.parquet                   ASV ranks
  cache/relic/asv_to_mag.tsv               (downloaded after barrnap+vsearch)
  cache/relic/bin_to_sample.tsv            (downloaded with barrnap output)

Output:
  cache/relic/asv_relic_classification.tsv
  cache/relic/relic_summary_per_compartment.tsv
  cache/relic/summary.txt
"""
from __future__ import annotations

from pathlib import Path
import re
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
OUT = CACHE / "relic"
OUT.mkdir(parents=True, exist_ok=True)


def parse_compartment(sid: str) -> str:
    s = sid.split("_")[-1]
    m = re.match(r"^[A-Z]*[0-9]+([A-Z]+)r[0-9]+$", s)
    if not m:
        return "?"
    return {"PR": "rhizosphere", "S": "surface", "D": "deep"}.get(m.group(1), "?")


def main():
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    tax = pd.read_parquet(CACHE / "taxonomy.parquet")
    print(f"feature_table: {ft.shape}", flush=True)

    if not (OUT / "asv_to_mag.tsv").exists():
        print(f"MISSING: {OUT}/asv_to_mag.tsv -- run barrnap+vsearch first", flush=True)
        return
    asv_to_mag = pd.read_csv(OUT / "asv_to_mag.tsv", sep="\t", header=None,
        names=["ASV", "MAG_label", "pid", "alen", "mismatches", "gaps",
               "qstart", "qend", "tstart", "tend", "evalue", "bitscore"])
    print(f"asv_to_mag: {len(asv_to_mag)} matches; "
          f"{asv_to_mag['ASV'].nunique()} ASVs with MAG support", flush=True)

    # Bin-to-sample membership (which sample each MAG_label was assembled from)
    bin_to_sample = pd.read_csv(OUT / "bin_to_sample.tsv", sep="\t")
    # MAG_label format: SAMPLEID__binNN -> sample = first half
    asv_to_mag["mag_sample"] = asv_to_mag["MAG_label"].str.split("__").str[0]

    # Build per-ASV best-hit MAG-sample set
    asv_mag_samples = asv_to_mag.groupby("ASV")["mag_sample"].apply(set).to_dict()
    asvs_with_mag = set(asv_to_mag["ASV"].unique())

    # 16S presence per ASV-sample
    pres = (ft > 0)
    sample_compartment = {s: parse_compartment(s) for s in ft.columns}
    print(f"compartment counts: {pd.Series(sample_compartment).value_counts().to_dict()}",
          flush=True)

    # Per-ASV summaries
    rows = []
    sample_set = set(ft.columns)
    for asv in pres.index:
        n_samples_16s = int(pres.loc[asv].sum())
        if n_samples_16s == 0:
            continue
        sup = asv_mag_samples.get(asv, set())
        sup_present = sup & sample_set
        rows.append({
            "ASV": asv,
            "n_samples_16S": n_samples_16s,
            "has_mag_match": asv in asvs_with_mag,
            "n_mag_supporting_samples": len(sup_present),
            "rel_total_reads": float(ft.loc[asv].sum() / max(ft.values.sum(), 1)),
        })
    res = pd.DataFrame(rows)

    # Classify
    def classify(r):
        if not r["has_mag_match"]:
            return "mag_orphan_likely_relic"
        if r["n_mag_supporting_samples"] >= 0.5 * r["n_samples_16S"]:
            return "active_consistent"
        return "intermittent_mag"

    res["classification"] = res.apply(classify, axis=1)
    res["genus"] = res["ASV"].map(tax["genus"].fillna("Unclassified"))
    res["phylum"] = res["ASV"].map(tax["phylum"].fillna("Unclassified"))
    res.to_csv(OUT / "asv_relic_classification.tsv", sep="\t", index=False)
    print("Per-ASV classification counts:")
    print(res["classification"].value_counts().to_string())

    # Per-compartment: fraction of 16S reads in each class
    rel_reads_by_class_compartment = []
    for comp in ["rhizosphere", "surface", "deep"]:
        comp_samples = [s for s, c in sample_compartment.items() if c == comp]
        if not comp_samples:
            continue
        sub = ft[comp_samples]
        total = sub.values.sum()
        for cls in res["classification"].unique():
            asvs = res.loc[res["classification"] == cls, "ASV"]
            r = float(sub.loc[sub.index.intersection(asvs)].values.sum() / max(total, 1))
            rel_reads_by_class_compartment.append({
                "compartment": comp, "class": cls,
                "fraction_of_reads": r,
                "n_asvs": int(len(asvs)),
            })
    cs = pd.DataFrame(rel_reads_by_class_compartment)
    cs.to_csv(OUT / "relic_summary_per_compartment.tsv", sep="\t", index=False)
    print("\nFraction of 16S reads per class per compartment:")
    print(cs.to_string(index=False))

    # Phylum-level signal
    phy = res.groupby(["phylum", "classification"])["rel_total_reads"].sum().unstack(fill_value=0)
    phy["total"] = phy.sum(axis=1)
    phy_top = phy.sort_values("total", ascending=False).head(15)
    if "mag_orphan_likely_relic" in phy_top.columns and "total" in phy_top.columns:
        phy_top["frac_relic"] = phy_top["mag_orphan_likely_relic"] / phy_top["total"]

    with open(OUT / "summary.txt", "w") as fh:
        fh.write("Relic DNA quantification (Tier 1 #3)\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Total ASVs: {len(res)}\n")
        fh.write(f"  with MAG match: {res['has_mag_match'].sum()}\n")
        fh.write(f"  MAG-orphan (likely relic): "
                 f"{(res['classification']=='mag_orphan_likely_relic').sum()}\n")
        fh.write(f"  active-consistent: "
                 f"{(res['classification']=='active_consistent').sum()}\n")
        fh.write(f"  intermittent-mag: "
                 f"{(res['classification']=='intermittent_mag').sum()}\n")
        fh.write("\nFraction of 16S reads per class per compartment:\n")
        fh.write(cs.to_string(index=False))
        fh.write("\n\nTop 15 phyla (relic fraction):\n")
        fh.write(phy_top.to_string())
    print(f"\nWrote {OUT}/summary.txt")


if __name__ == "__main__":
    main()
