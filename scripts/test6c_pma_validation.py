#!/usr/bin/env python3
"""TEST 6C: PMA validation of the ephemeral-OTU finding (revised).

PMA ASVs were derived from a separate DADA2 run with slightly different
primers (RV: GACTACHVGGGTATCTAATCC vs main GGACTACNVGGGTWTCTAAT), so
ASV IDs do not match between datasets. We use vsearch (97% ID) to map
PMA ASVs -> EQ ASVs by sequence similarity.

Then for each EQ ASV with a PMA proxy:
  T:UT relabund ratio per pair (9 pairs, sites C1+C2 trip-5).
  median(T)/median(UT) -> per-EQ-ASV viability score.

Compare viability scores across EQ persistence classes.

Inputs:
  /home/leechuck/Public/software/empty-quarter/relic-dna/ASV_table_rel_abundance.tsv
  cache/test6_disconfirmation/pma_to_eq_match.tsv  (vsearch output)
  cache/test6_persistence/per_OTU_site_persistence.parquet

Output:
  cache/test6_disconfirmation/per_eq_asv_viability.tsv
  cache/test6_disconfirmation/test6c_summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

CACHE = REPO / "cache"
RELIC = Path("/home/leechuck/Public/software/empty-quarter/relic-dna")
OUT = CACHE / "test6_disconfirmation"


def parse_relic_sample(name: str):
    m = re.match(r"^(C[12])([RS])(\d+)(T|UT)$", name)
    if not m: return None
    site, comp, rep, treat = m.groups()
    return {"site": site, "comp_code": comp, "rep": int(rep),
            "treatment": treat}


def main():
    rel = pd.read_csv(RELIC / "ASV_table_rel_abundance.tsv", sep="\t",
                      index_col=0)
    samp = []
    for s in rel.columns:
        m = parse_relic_sample(s)
        if m: m["sample_orig"] = s; samp.append(m)
    sm = pd.DataFrame(samp)
    pairs = sm.pivot(index=["site", "comp_code", "rep"],
                      columns="treatment", values="sample_orig").dropna(
                          subset=["T", "UT"])
    print(f"PMA T/UT pairs: {len(pairs)}", flush=True)

    # PMA ASV -> EQ ASV mapping from vsearch
    map_df = pd.read_csv(OUT / "pma_to_eq_match.tsv", sep="\t", header=None,
                          names=["pma_asv", "eq_asv", "pid", "alen", "mismatches",
                                 "gaps", "qstart", "qend", "tstart", "tend",
                                 "evalue", "bitscore"])
    print(f"PMA->EQ map: {len(map_df)} rows ({map_df['pma_asv'].nunique()} unique PMA, "
          f"{map_df['eq_asv'].nunique()} unique EQ)", flush=True)

    # Per pair, compute T/UT log-ratio per PMA ASV
    eps = 1e-7
    R_recs = []
    for (site, comp_code, rep), row in pairs.iterrows():
        T_col = row["T"]; UT_col = row["UT"]
        T = rel[T_col]; UT = rel[UT_col]
        # ratio
        ratio = (T + eps) / (UT + eps)
        ratio = ratio.clip(upper=10)
        for asv_id, r in ratio.items():
            R_recs.append({"pma_asv": asv_id, "site": site,
                            "comp_code": comp_code, "rep": rep,
                            "T_relabund": float(T[asv_id]),
                            "UT_relabund": float(UT[asv_id]),
                            "ratio": float(r)})
    R = pd.DataFrame(R_recs)
    # Restrict to PMA ASVs with UT_relabund > 0 (true detections in untreated)
    R_present = R[R["UT_relabund"] > 0]
    print(f"PMA ASV-pair detections (UT>0): {len(R_present)}", flush=True)

    # Median ratio per PMA ASV
    pma_viab = (R_present.groupby("pma_asv")
                .agg(median_ratio=("ratio", "median"),
                     mean_ratio=("ratio", "mean"),
                     n_pairs_present=("ratio", "count"),
                     mean_UT=("UT_relabund", "mean"))
                .reset_index())
    print(f"PMA ASVs with viability score: {len(pma_viab)}", flush=True)

    # Map to EQ ASVs (one PMA ASV -> one EQ ASV via vsearch best hit)
    pma_to_eq = map_df.drop_duplicates("pma_asv")[["pma_asv", "eq_asv", "pid"]]
    pma_viab = pma_viab.merge(pma_to_eq, on="pma_asv", how="left")

    # Aggregate EQ ASV viability (if multiple PMA ASVs map to same EQ ASV,
    # average their median ratios weighted by mean_UT)
    eq_viab = (pma_viab.dropna(subset=["eq_asv"])
               .groupby("eq_asv")
               .apply(lambda g: pd.Series({
                   "n_pma_proxies": len(g),
                   "weighted_median_ratio": np.average(g["median_ratio"],
                                                       weights=g["mean_UT"]
                                                                if g["mean_UT"].sum() > 0
                                                                else None),
                   "mean_ratio": float(g["median_ratio"].mean()),
                   "median_ratio": float(g["median_ratio"].median()),
                   "max_pid": float(g["pid"].max()),
               })).reset_index())
    print(f"EQ ASVs with PMA viability proxy: {len(eq_viab)}", flush=True)
    eq_viab.to_csv(OUT / "per_eq_asv_viability.tsv", sep="\t", index=False)

    # Map persistence classes
    persist = pd.read_parquet(CACHE / "test6_persistence" /
                               "per_OTU_site_persistence.parquet")
    per_asv_persist = (persist.groupby("ASV")
                       .agg(max_n_trips=("n_trips_present", "max"),
                            mean_n_trips=("n_trips_present", "mean"),
                            n_site_records=("ASV", "count"))
                       .reset_index().rename(columns={"ASV": "eq_asv"}))
    merged = per_asv_persist.merge(eq_viab, on="eq_asv", how="inner")
    print(f"\nEQ ASVs with both persistence + PMA viability: {len(merged)}",
          flush=True)

    # Distribution of viability ratio per persistence class
    print("\n=== Viability T/UT (median per ASV) by max-persistence class ===")
    agg = (merged.groupby("max_n_trips")["weighted_median_ratio"]
           .agg(["count", "median", "mean",
                 lambda x: float(np.percentile(x, 25)),
                 lambda x: float(np.percentile(x, 75))])
           .rename(columns={"<lambda_0>": "p25", "<lambda_1>": "p75"}))
    print(agg.to_string())

    # Statistical: do core OTUs have higher viability than ephemeral?
    eph = merged[merged["max_n_trips"] == 1]
    core = merged[merged["max_n_trips"] == 5]
    if len(eph) > 5 and len(core) > 5:
        U, p = mannwhitneyu(core["weighted_median_ratio"],
                             eph["weighted_median_ratio"],
                             alternative="greater")
        print(f"\nMann-Whitney (core viability > ephemeral viability):")
        print(f"  core median viability: {core['weighted_median_ratio'].median():.3f} (n={len(core)})")
        print(f"  eph  median viability: {eph['weighted_median_ratio'].median():.3f} (n={len(eph)})")
        print(f"  U={U:.0f}, p={p:.3g}")
    # Kruskal across all classes
    groups = [merged.loc[merged["max_n_trips"] == k,
                         "weighted_median_ratio"].dropna() for k in (1, 2, 3, 4, 5)]
    groups = [g for g in groups if len(g) > 3]
    if len(groups) >= 2:
        H, p_kw = kruskal(*groups)
        print(f"\nKruskal-Wallis across persistence classes: H={H:.2f}, p={p_kw:.3g}")

    with open(OUT / "test6c_summary.txt", "w") as fh:
        fh.write("Test 6C: PMA validation of ephemeral-OTU finding (revised)\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"PMA T/UT pairs: {len(pairs)}\n")
        fh.write(f"PMA ASVs in UT-detection: {pma_viab.shape[0]}\n")
        fh.write(f"PMA->EQ mapped (vsearch 97% ID): {len(pma_to_eq)}\n")
        fh.write(f"EQ ASVs with PMA viability proxy: {len(eq_viab)}\n")
        fh.write(f"EQ ASVs with both persistence + viability: {len(merged)}\n\n")

        fh.write("Per max-persistence class, viability T/UT distribution:\n")
        fh.write(agg.to_string())
        if len(eph) > 5 and len(core) > 5:
            fh.write(f"\n\nCore vs Ephemeral viability:\n")
            fh.write(f"  core (5-trip) median: {core['weighted_median_ratio'].median():.3f}  (n={len(core)})\n")
            fh.write(f"  ephemeral (1-trip) median: {eph['weighted_median_ratio'].median():.3f}  (n={len(eph)})\n")
            fh.write(f"  Mann-Whitney p (core>eph): {p:.3g}\n")
        if len(groups) >= 2:
            fh.write(f"  Kruskal-Wallis across classes: H={H:.2f}, p={p_kw:.3g}\n")

        fh.write("\nDISCONFIRMATION KEY:\n")
        fh.write("  T/UT ~ 1.0 -> ASV is alive (PMA-resistant -> viable cell)\n")
        fh.write("  T/UT << 1.0 -> ASV is mostly relic DNA (PMA crosslinked it)\n")
        fh.write("  If core viability >> ephemeral viability AND p<0.01, the\n"
                 "    'ephemeral 67%' is largely relic DNA. Core IS the true\n"
                 "    living community.\n")
        fh.write("  If both classes show similar viability (~1.0), ephemerality\n"
                 "    is real, alive but episodic.\n")
        fh.write("  CAVEATS:\n")
        fh.write("    - PMA samples: Trip 5 only, sites C1+C2 (not whole biome)\n")
        fh.write("    - Different primers between datasets; vsearch 97% mapping\n")
        fh.write("      may miss some matches; 4277 of 6011 PMA ASVs (71%) mapped.\n")
    print(f"\nWrote {OUT}/test6c_summary.txt")


if __name__ == "__main__":
    main()
