#!/usr/bin/env python3
"""TEST 5: Osmolyte uptake vs biosynthesis ratio (Black Queen hypothesis).

Predicted (= field expectation): salinity rises -> compatible-solute
biosynthesis pathways (proB/proA/proC, betA, ectA/B/C) increase at
the community level.

Surprising direction: as salinity rises, the COMMUNITY shifts from
biosynthesis (genes-cost-energy) to uptake (genes-need-suppliers),
indicating Black-Queen-like collective osmoadaptation where most
members rely on a few producers.

Approach:
  16S+PICRUSt2 side (full 1227 samples): use EC abundances for the
  biosynthesis enzymes we have (no transporters in PICRUSt2 EC table).
    - betA = EC:1.1.99.1   (choline dehydrogenase)
    - proB = EC:2.7.2.11   (glutamate 5-kinase)
    - proA = EC:1.2.1.41   (glutamate-5-semialdehyde DH)
    - proC = EC:1.5.1.2    (pyrroline-5-carboxylate reductase)
    - ectA = EC:2.3.1.178  (DAB acetyltransferase)
    - ectB = EC:2.6.1.76   (DAB aminotransferase)
    - ectC = EC:4.2.1.108  (ectoine synthase)
    - opuD/proU TRANSPORTERS NOT IN EC TABLE
  Metagenomic side (296 samples; from cache/betA_guild_census + uptake):
    - K00108 (betA) for biosynthesis (presence per bin)
    - proU/opuD/opuB ... for uptake (presence per bin)

Tests:
  A. Per-sample biosynthesis abundance vs salinity (PICRUSt2 EC).
  B. Metagenomic biosynthesis-bin count vs uptake-bin count vs salinity.

Inputs:
  data/functional/picrust2/metagenome_pred_metagenome_unstrat.tsv
  data/geochemistry/xrf_lab_table_all_trips.tsv
  cache/betA_per_sample_summary.tsv
  cache/betaine_uptake_census.tsv

Outputs:
  cache/test5_osmolyte/picrust_biosynth_vs_salinity.tsv
  cache/test5_osmolyte/metagenomic_uptake_vs_biosynth.tsv
  cache/test5_osmolyte/summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "test5_osmolyte"
OUT.mkdir(parents=True, exist_ok=True)

BIOSYNTH_EC = {
    "betA":  "EC:1.1.99.1",
    "proB":  "EC:2.7.2.11",
    "proA":  "EC:1.2.1.41",
    "proC":  "EC:1.5.1.2",
    "ectA":  "EC:2.3.1.178",
    "ectB":  "EC:2.6.1.76",
    "ectC":  "EC:4.2.1.108",
    "treS":  "EC:5.4.99.16",   # trehalose synthase
    "otsA":  "EC:2.4.1.15",    # trehalose-6-P synthase
    "glgA":  "EC:2.4.1.21",    # glycogen synthase
}


def main():
    ec = pd.read_csv(DATA / "functional" / "picrust2" /
                     "metagenome_pred_metagenome_unstrat.tsv",
                     sep="\t", index_col=0)
    print(f"PICRUSt2 EC table: {ec.shape}", flush=True)

    available = {k: v for k, v in BIOSYNTH_EC.items() if v in ec.index}
    missing = {k: v for k, v in BIOSYNTH_EC.items() if v not in ec.index}
    print(f"available biosynthesis ECs: {available}", flush=True)
    print(f"missing biosynthesis ECs:   {missing}", flush=True)

    # Per-sample sum of biosynthesis ECs (relative to total)
    total = ec.sum(axis=0)
    biosynth_sum = ec.loc[list(available.values())].sum(axis=0)
    biosynth_rel = biosynth_sum / total.replace(0, 1)

    smeta = parse_samples_to_df(ec.columns)
    smeta["site"] = smeta["site"].astype(int)
    sm = smeta.set_index("sample")
    sm["biosynth_rel"] = biosynth_rel

    # Per-EC relative
    for nm, ec_id in available.items():
        sm[f"{nm}_rel"] = ec.loc[ec_id] / total.replace(0, 1)

    # Salinity from XRF
    xrf = pd.read_csv(DATA / "geochemistry" / "xrf_lab_table_all_trips.tsv",
                       sep="\t")
    # XRF is per (trip, site, compartment) cell; aggregate
    xrf_clean = (xrf[["trip", "site", "compartment", "S"]]
                 .dropna(subset=["S"])
                 .assign(compartment=xrf["compartment"].str.lower())
                 .groupby(["trip", "site", "compartment"])["S"].mean()
                 .reset_index())
    sm = sm.reset_index()
    sm = sm.merge(xrf_clean, on=["trip", "site", "compartment"], how="left")
    n_with_S = sm["S"].notna().sum()
    print(f"\nsamples with XRF S: {n_with_S}", flush=True)

    # Spearman: biosynth_rel vs S, per compartment
    rows = []
    for comp, grp in sm.dropna(subset=["S"]).groupby("compartment"):
        r, p = spearmanr(grp["biosynth_rel"], grp["S"])
        rows.append({"compartment": comp, "n": len(grp),
                     "rho_biosynth_vs_S": r, "p_value": p,
                     "median_biosynth_rel": float(grp["biosynth_rel"].median()),
                     "median_S": float(grp["S"].median())})
    rdf = pd.DataFrame(rows)
    print("\nPICRUSt2 biosynthesis vs salinity per compartment:")
    print(rdf.to_string(index=False))

    per_ec_rows = []
    for nm, ec_id in available.items():
        for comp, grp in sm.dropna(subset=["S"]).groupby("compartment"):
            r, p = spearmanr(grp[f"{nm}_rel"], grp["S"])
            per_ec_rows.append({"compartment": comp, "gene": nm,
                                 "EC": ec_id, "n": len(grp),
                                 "rho_vs_S": r, "p_value": p})
    pec = pd.DataFrame(per_ec_rows)
    pec.to_csv(OUT / "picrust_biosynth_vs_salinity.tsv", sep="\t", index=False)

    print("\nPer-EC Spearman vs salinity:")
    print(pec.pivot_table(index="gene", columns="compartment",
                           values="rho_vs_S").round(3).to_string())

    # Metagenomic uptake-vs-biosynthesis: pull from existing betA + uptake census
    bs = pd.read_csv(CACHE / "betA_per_sample_summary.tsv", sep="\t")
    print(f"\nbetA per_sample_summary: {len(bs)}", flush=True)

    upt_path = CACHE / "betaine_uptake_census.tsv"
    if upt_path.exists():
        upt = pd.read_csv(upt_path, sep="\t")
        print(f"betaine_uptake_census: {len(upt)} bin records, cols: "
              f"{list(upt.columns)[:8]}", flush=True)
        # Aggregate to per-sample
        if "sample_id" in upt.columns:
            ucol = [c for c in upt.columns if "K0" in c or "uptake" in c.lower()]
            if ucol:
                upt_per_sample = upt.groupby("sample_id")[ucol].sum().reset_index()
                upt_per_sample["uptake_total"] = upt_per_sample[ucol].sum(axis=1)
                merged = bs.merge(upt_per_sample[["sample_id", "uptake_total"]],
                                   on="sample_id", how="outer").fillna(0)
                merged["ratio_uptake_over_biosynth"] = (
                    merged["uptake_total"] /
                    (merged["n_K00108_strict"].replace(0, 1)))
                merged.to_csv(OUT / "metagenomic_uptake_vs_biosynth.tsv",
                              sep="\t", index=False)
                print("\nMetagenomic per-sample biosynthesis (K00108) vs uptake total:")
                print(f"  total samples: {len(merged)}")
                print(f"  median K00108 per sample: {merged['n_K00108_strict'].median():.1f}")
                print(f"  median uptake genes per sample: {merged['uptake_total'].median():.1f}")
                print(f"  median ratio uptake/K00108: "
                      f"{merged['ratio_uptake_over_biosynth'].median():.2f}")
            else:
                print("Couldn't identify uptake-gene cols in upt", flush=True)

    # Summary
    with open(OUT / "summary.txt", "w") as fh:
        fh.write("Test 5: Osmolyte uptake vs biosynthesis (Black Queen)\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"PICRUSt2 EC biosynthesis pool:\n")
        for nm, ec_id in BIOSYNTH_EC.items():
            avail = "AVAILABLE" if ec_id in ec.index else "MISSING"
            fh.write(f"  {nm:>6s} ({ec_id}): {avail}\n")
        fh.write(f"\nSamples in PICRUSt2 EC table: {ec.shape[1]}\n")
        fh.write(f"Samples with XRF S: {n_with_S}\n\n")

        fh.write("Total biosynthesis (EC pool) vs salinity per compartment:\n")
        fh.write(rdf.to_string(index=False))

        fh.write("\n\nPer-EC Spearman rho vs salinity (rows=gene, cols=compartment):\n")
        fh.write(pec.pivot_table(index="gene", columns="compartment",
                                  values="rho_vs_S").round(3).to_string())

        fh.write("\n\nINTERPRETATION KEY:\n")
        fh.write("  rho_biosynth_vs_S > 0.2 (and significant) ->\n"
                 "    EXPECTED: biosynthesis genes increase with salinity\n")
        fh.write("  rho_biosynth_vs_S ~ 0 ->\n"
                 "    no community-level biosynthesis response — suggests UPTAKE\n"
                 "    dominates (Black Queen direction)\n")
        fh.write("  rho < 0 -> SURPRISING: biosynthesis decreases with salinity\n")
    print(f"\nWrote {OUT}/summary.txt")


if __name__ == "__main__":
    main()
