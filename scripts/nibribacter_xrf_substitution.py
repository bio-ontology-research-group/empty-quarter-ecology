#!/usr/bin/env python3
"""Nibribacter abundance vs XRF + climate + Bacilli reverse substitution test."""
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
OUT = CACHE / "keystone_test"
OUT.mkdir(parents=True, exist_ok=True)


def relabund(M):
    return M.div(M.sum(axis=0).replace(0, 1), axis=1)


def main():
    p = pd.read_csv(CACHE / "relic_priors" /
                     "relic_score_with_mag_prior.tsv", sep="\t")
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)
    tax = pd.read_parquet(CACHE / "taxonomy.parquet").reset_index().rename(
        columns={"ASV": "asv_id"})

    alive_mag = set(p.loc[p["relic_score_with_mag"] <= 0.3, "asv_id"])
    ft_a = ft.loc[ft.index.isin(alive_mag)]

    # Per-genus aggregation
    ft2 = ft_a.copy()
    ft2.index = ft2.index.rename("asv_id")
    m = ft2.reset_index().merge(tax[["asv_id", "phylum", "genus"]],
                                  on="asv_id", how="left")
    m = m.dropna(subset=["genus"])
    m = m[~m["genus"].astype(str).isin(["NA", ""])]
    sample_cols = [c for c in m.columns if c not in
                       ("asv_id", "phylum", "genus")]
    gen = m.groupby("genus")[sample_cols].sum()
    rel = relabund(gen)

    # ============================================================
    # 1. Nibribacter vs XRF
    # ============================================================
    print("=== Nibribacter abundance vs XRF chemistry ===\n", flush=True)
    xrf = pd.read_csv(DATA / "geochemistry" /
                       "xrf_lab_table_all_trips.tsv", sep="\t")
    if "site" in xrf.columns:
        xrf["site"] = pd.to_numeric(xrf["site"], errors="coerce")
        xrf = xrf.dropna(subset=["site"])
    if "compartment" in xrf.columns:
        xrf["compartment"] = xrf["compartment"].astype(str).str.lower()

    # Per (site, compartment) median XRF
    xrf_num = xrf.select_dtypes(include="number").drop(columns=["site"],
                                                            errors="ignore")
    xrf_num["site"] = xrf["site"].astype(int).values
    if "compartment" in xrf.columns:
        xrf_num["compartment"] = xrf["compartment"].values
        xrf_per = (xrf_num.groupby(["site", "compartment"]).median()
                    .reset_index())
    else:
        xrf_per = xrf_num.groupby("site").median().reset_index()

    nib = rel.loc["Nibribacter"]
    sm = smeta.merge(nib.rename("nib_relabund"), left_on="sample",
                       right_index=True, how="left")
    sm = sm.merge(xrf_per, on=["site", "compartment"]
                    if "compartment" in xrf_per.columns else ["site"],
                    how="left")

    # Spearman per element
    rec = []
    for col in sm.select_dtypes(include="number").columns:
        if col in ("nib_relabund", "site", "trip"): continue
        sub = sm[["nib_relabund", col]].dropna()
        if len(sub) < 30: continue
        if sub[col].std() == 0: continue
        r, pp = spearmanr(sub["nib_relabund"], sub[col])
        rec.append({"variable": col, "n": len(sub),
                      "rho": float(r), "p": float(pp)})
    df = pd.DataFrame(rec).sort_values("rho", key=lambda x: x.abs(),
                                            ascending=False)
    df.to_csv(OUT / "nibribacter_xrf_corr.tsv", sep="\t", index=False)
    print(df.head(20).round(4).to_string(index=False), flush=True)

    # Climate too
    geo_t1 = pd.read_csv(DATA / "geodata" / "trip1_geodata.tsv", sep="\t")
    geo_t1 = geo_t1.rename(columns={"Site": "site"})
    geo_t1["site"] = pd.to_numeric(geo_t1["site"],
                                       errors="coerce").astype("Int64")
    sm2 = sm.merge(geo_t1[["site", "Latitude", "Longitude",
                                "AnnualMeanTemp", "AnnualTotalPrecip"]],
                     on="site", how="left")
    print("\n  Climate correlations:")
    for v in ("Latitude", "Longitude", "AnnualMeanTemp",
                "AnnualTotalPrecip"):
        if v in sm2.columns:
            sub = sm2[["nib_relabund", v]].dropna()
            r, pp = spearmanr(sub["nib_relabund"], sub[v])
            print(f"    {v}: rho={r:+.3f}  p={pp:.3g}  n={len(sub)}",
                  flush=True)

    # ============================================================
    # 2. Bacilli reverse substitution test
    # ============================================================
    print("\n=== Bacilli reverse substitution test (correlations with Aquibacillus) ===")
    if "Aquibacillus" not in rel.index:
        print("  Aquibacillus not in genus list"); return
    aqui = rel.loc["Aquibacillus"]

    bacilli_partners = ["Sediminibacillus", "Litchfieldia", "Tumebacillus",
                          "Gracilibacillus", "Oceanobacillus",
                          "Polygonibacillus", "Salirhabdus",
                          "Alkalihalobacillus"]
    bact_dom = ["Nibribacter", "Daejeonella", "Cytophaga", "Cnuella",
                  "Niastella", "Flavitalea", "Ohtaekwangia",
                  "Tellurirhabdus"]
    pseudo = ["Massilia", "Halomonas", "Acidibacter", "Lysobacter"]

    rec2 = []
    for partner in bacilli_partners + bact_dom + pseudo:
        if partner == "Aquibacillus": continue
        if partner not in rel.index: continue
        v = rel.loc[partner]
        common = aqui.index.intersection(v.index)
        r, _ = spearmanr(aqui.loc[common], v.loc[common])
        rec2.append({"partner": partner, "rho_with_Aquibacillus": float(r),
                      "guild": ("Bacilli" if partner in bacilli_partners
                                  else "Bact_DOM" if partner in bact_dom
                                  else "Pseudo")})
    df2 = pd.DataFrame(rec2).sort_values("rho_with_Aquibacillus")
    df2.to_csv(OUT / "aquibacillus_substitution.tsv", sep="\t", index=False)
    print(df2.round(3).to_string(index=False), flush=True)

    # ============================================================
    # 3. Per-sample dominant strategy classification
    # ============================================================
    print("\n=== Per-sample dominant strategy ===")
    bact_dom_sum = rel.loc[[g for g in bact_dom if g in rel.index]].sum(axis=0)
    bacilli_sum = rel.loc[[g for g in bacilli_partners
                                if g in rel.index]].sum(axis=0)
    pseudo_sum = rel.loc[[g for g in pseudo if g in rel.index]].sum(axis=0)

    df3 = pd.DataFrame({
        "sample": bact_dom_sum.index,
        "bact_dom_sum": bact_dom_sum.values,
        "bacilli_sum": bacilli_sum.values,
        "pseudo_sum": pseudo_sum.values,
    })
    df3 = df3.merge(smeta, on="sample", how="left")

    # Classify
    def classify(row):
        s = {"BactDOM": row["bact_dom_sum"], "Bacilli": row["bacilli_sum"],
              "Pseudo": row["pseudo_sum"]}
        return max(s, key=s.get)
    df3["dominant_guild"] = df3.apply(classify, axis=1)
    print("\n  Dominant guild per compartment:")
    print(df3.groupby(["compartment", "dominant_guild"]).size()
          .unstack(fill_value=0).to_string())

    # Total alive abundance per strategy
    print("\n  Median guild relabund per dominant_guild:")
    print(df3.groupby("dominant_guild")[["bact_dom_sum", "bacilli_sum",
                                              "pseudo_sum"]].median().round(4)
          .to_string())

    df3.to_csv(OUT / "per_sample_strategy.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
