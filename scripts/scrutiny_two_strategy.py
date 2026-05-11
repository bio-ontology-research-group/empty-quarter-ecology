#!/usr/bin/env python3
"""Scrutiny + mechanism analysis of the two-strategy hypothesis.

Tests:
  1. COMPOSITIONAL ARTIFACT: are A vs B truly anti-correlated, or is it
     just sum-to-1 effect? Use CLR-transformed abundances and
     absolute-abundance proxies.
  2. WITHIN-SITE ANTI-CORRELATION: do A and B anti-correlate WITHIN sites
     (where total reads are similar), not just across sites?
  3. STRATEGY B FUNCTIONAL MECHANISM: pull PICRUSt2 KO content for Bacilli
     + Halomonas members; identify ectoine biosynthesis, sporulation, salt
     tolerance, anaerobic capacity.
  4. TRANSITION MECHANISM: for the 51 switching cells, what climate /
     environmental variable changes between A-dominant and B-dominant trips?
  5. NIBRIBACTER FUNCTIONAL ARTIFACT: list specific betaine biosynthesis KOs
     to verify the 96-ORF count is meaningful or inflated by alcohol-DH
     homologs.

Outputs in cache/two_strategy_scrutiny/.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "two_strategy_scrutiny"
OUT.mkdir(parents=True, exist_ok=True)

PICRUSt2_KO_PATH = Path("/home/leechuck/Public/software/empty-quarter/data/"
                          "processed/functional/picrust2/KO_predicted.tsv")

STRATEGY_A = ["Nibribacter", "Daejeonella", "Cytophaga", "Cnuella",
                  "Niastella", "Flavitalea", "Ohtaekwangia", "Tellurirhabdus",
                  "Massilia", "Lysobacter", "Acidibacter"]
STRATEGY_B = ["Aquibacillus", "Sediminibacillus", "Litchfieldia",
                  "Tumebacillus", "Gracilibacillus", "Oceanobacillus",
                  "Polygonibacillus", "Salirhabdus", "Alkalihalobacillus",
                  "Halomonas"]


def relabund(M):
    return M.div(M.sum(axis=0).replace(0, 1), axis=1)


def clr(M, pseudo=0.5):
    """CLR transform on relabund table (rows=features, cols=samples)."""
    R = M + pseudo
    LR = np.log(R)
    return LR.sub(LR.mean(axis=0), axis=1)


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

    # ==================================================================
    # 1. COMPOSITIONAL ARTIFACT CHECK
    # ==================================================================
    print("\n========================================================")
    print("[1] Compositional artifact check")
    print("========================================================")
    A_present = [g for g in STRATEGY_A if g in rel.index]
    B_present = [g for g in STRATEGY_B if g in rel.index]

    # Method 1: relabund (compositional)
    sA_rel = rel.loc[A_present].sum(axis=0)
    sB_rel = rel.loc[B_present].sum(axis=0)
    r_rel, p_rel = spearmanr(sA_rel, sB_rel)
    print(f"  A vs B Spearman (relabund):         rho={r_rel:+.3f}  "
          f"p={p_rel:.3g}", flush=True)

    # Method 2: ABSOLUTE counts (NOT relabund)
    sA_abs = gen.loc[A_present].sum(axis=0)
    sB_abs = gen.loc[B_present].sum(axis=0)
    r_abs, p_abs = spearmanr(sA_abs, sB_abs)
    print(f"  A vs B Spearman (absolute counts):  rho={r_abs:+.3f}  "
          f"p={p_abs:.3g}", flush=True)

    # Method 3: CLR-transformed
    rel_clr = clr(rel)
    sA_clr = rel_clr.loc[A_present].mean(axis=0)
    sB_clr = rel_clr.loc[B_present].mean(axis=0)
    r_clr, p_clr = spearmanr(sA_clr, sB_clr)
    print(f"  A vs B Spearman (CLR mean):         rho={r_clr:+.3f}  "
          f"p={p_clr:.3g}", flush=True)

    # Method 4: random-genus pair as a null
    print(f"\n  RANDOM-GENUS PAIR NULL (100 random genus subsets):")
    rng = np.random.default_rng(20260510)
    null = []
    all_g = [g for g in rel.index if g not in A_present + B_present]
    for _ in range(100):
        g1 = list(rng.choice(all_g, len(A_present), replace=False))
        g2 = list(rng.choice([g for g in all_g if g not in g1],
                                len(B_present), replace=False))
        s1 = rel.loc[g1].sum(axis=0)
        s2 = rel.loc[g2].sum(axis=0)
        r, _ = spearmanr(s1, s2)
        null.append(r)
    null = np.array(null)
    print(f"    median: {np.median(null):+.3f}  "
          f"p5..p95: [{np.percentile(null, 5):+.3f}, "
          f"{np.percentile(null, 95):+.3f}]", flush=True)
    p_extreme = float((null < r_rel).mean())
    print(f"    A-B rho ({r_rel:+.3f}) is "
          f"{'BELOW' if p_extreme > 0.5 else 'WITHIN'} the null "
          f"distribution (extreme-low p={p_extreme:.3f})", flush=True)

    pd.DataFrame([{
        "method": "relabund_total", "rho_AvB": r_rel, "p": p_rel,
        "rho_null_median": float(np.median(null)),
        "rho_null_p5": float(np.percentile(null, 5)),
        "rho_null_p95": float(np.percentile(null, 95)),
        "p_extreme_negative": float(p_extreme),
    }, {"method": "absolute_counts", "rho_AvB": r_abs, "p": p_abs},
       {"method": "CLR_mean", "rho_AvB": r_clr, "p": p_clr},
    ]).to_csv(OUT / "compositional_artifact.tsv", sep="\t", index=False)

    # ==================================================================
    # 2. WITHIN-SITE ANTI-CORRELATION
    # ==================================================================
    print("\n========================================================")
    print("[2] Within-site anti-correlation")
    print("========================================================")
    df = pd.DataFrame({"sample": sA_rel.index,
                          "A_rel": sA_rel.values,
                          "B_rel": sB_rel.values}).merge(smeta, on="sample")
    rec = []
    for site in sorted(df["site"].unique()):
        for comp in ("rhizosphere", "surface", "deep"):
            sub = df[(df["site"] == site) & (df["compartment"] == comp)]
            if len(sub) < 4: continue
            r, _ = spearmanr(sub["A_rel"], sub["B_rel"])
            rec.append({"site": site, "compartment": comp,
                          "n_samples": len(sub), "rho": float(r)})
    wd = pd.DataFrame(rec).dropna()
    print(f"  Per-(site, comp) rho A vs B (n_samples>=4):")
    print(f"    n cells: {len(wd)}", flush=True)
    print(f"    median rho: {wd['rho'].median():+.3f}", flush=True)
    print(f"    fraction negative: {(wd['rho'] < 0).mean():.3f}", flush=True)
    print(f"    fraction strongly negative (<-0.3): "
          f"{(wd['rho'] < -0.3).mean():.3f}", flush=True)
    wd.to_csv(OUT / "within_site_anti_corr.tsv", sep="\t", index=False)

    # ==================================================================
    # 3. STRATEGY B FUNCTIONAL MECHANISM
    # ==================================================================
    print("\n========================================================")
    print("[3] Strategy B functional mechanism (Bacilli + Halomonas KOs)")
    print("========================================================")
    if PICRUSt2_KO_PATH.exists():
        ko_pred = pd.read_csv(PICRUSt2_KO_PATH, sep="\t",
                                index_col=0, dtype={"sequence": str})
        ko_pred = ko_pred.astype(np.float32)

        asv_to_genus = (m[["asv_id", "genus"]].set_index("asv_id")
                          ["genus"].to_dict())
        common = ko_pred.index.intersection(set(asv_to_genus))
        ko_a = ko_pred.loc[common]
        gen_for_ko = pd.Series([asv_to_genus[a] for a in common],
                                  index=common)
        ko_per_genus = ko_a.groupby(gen_for_ko).mean()
        ko_per_genus = ko_per_genus.loc[:, ko_per_genus.sum(axis=0) > 0]

        # Strategy B mean KO copies (across members)
        B_in = [g for g in B_present if g in ko_per_genus.index]
        A_in = [g for g in A_present if g in ko_per_genus.index]
        B_mean = ko_per_genus.loc[B_in].mean(axis=0)
        A_mean = ko_per_genus.loc[A_in].mean(axis=0)
        # log2 ratio per KO: B vs A
        ratio = np.log2((B_mean + 0.001) / (A_mean + 0.001))

        print("  Strategy-B-enriched KOs (top 30 by log2(B/A)):")
        # Specific marker KOs
        markers = {
            # Sporulation
            "K06375": "spo0A (sporulation initiation regulator)",
            "K06376": "spo0F (sporulation initiation phosphotransfer)",
            "K06378": "spo0E",
            "K06398": "spoIIE",
            "K07669": "spoIIIE",
            "K07343": "spoIVA",
            # Compatible solutes
            "K06718": "ectA (L-2,4-diaminobutyric acid acetyltransferase)",
            "K10674": "ectB (DABA aminotransferase)",
            "K10673": "ectC (ectoine synthase)",
            "K10675": "doeA (ectoine hydroxylase)",
            "K01598": "betA (choline dehydrogenase)",
            "K00130": "betB (betaine aldehyde DH)",
            "K05845": "opuD (glycine betaine transporter)",
            "K05846": "opuC",
            # Trehalose
            "K00697": "otsA (trehalose-6P synthase)",
            "K01087": "otsB (trehalose-6P phosphatase)",
            # Salt-tolerance / Na-pump
            "K03313": "Na+/H+ antiporter NhaA",
            "K07301": "Na+/H+ antiporter MrpA",
            "K05384": "Na+-translocating NADH-quinone reductase",
            # Anaerobic / sulfate respiration
            "K11181": "dsrA",
            "K00958": "sat (sulfate adenylyltransferase)",
            # Stress response
            "K04077": "groEL",
            "K04043": "dnaK",
            "K03657": "uvrD",
            "K03701": "uvrA",
        }
        rec = []
        for ko, desc in markers.items():
            if ko in ko_per_genus.columns:
                B_v = float(ko_per_genus.loc[B_in, ko].mean())
                A_v = float(ko_per_genus.loc[A_in, ko].mean())
                r = float(np.log2((B_v + 0.001) / (A_v + 0.001)))
            else:
                B_v = A_v = r = np.nan
            rec.append({"ko": ko, "function": desc,
                          "B_mean": B_v, "A_mean": A_v,
                          "log2_B_over_A": r})
        kdf = pd.DataFrame(rec)
        print(kdf.round(3).to_string(index=False), flush=True)
        kdf.to_csv(OUT / "strategy_B_marker_KOs.tsv", sep="\t", index=False)

        # Top 20 KOs most enriched in B vs A (any KO)
        ko_diff = pd.DataFrame({"ko": ratio.index,
                                  "B_mean": B_mean.values,
                                  "A_mean": A_mean.values,
                                  "log2_B_over_A": ratio.values})
        ko_diff = ko_diff[ko_diff["B_mean"] > 0.05]  # filter rare
        ko_diff = ko_diff.sort_values("log2_B_over_A", ascending=False)
        print(f"\n  Top 30 B-enriched KOs:")
        print(ko_diff.head(30).round(3).to_string(index=False), flush=True)
        ko_diff.to_csv(OUT / "ko_diff_B_vs_A.tsv", sep="\t", index=False)

    # ==================================================================
    # 4. TRANSITION MECHANISM
    # ==================================================================
    print("\n========================================================")
    print("[4] Transition mechanism (switching cells)")
    print("========================================================")
    pst = pd.read_csv(CACHE / "two_strategy_temporal" /
                       "per_sample_strategy_with_precip.tsv", sep="\t")
    pst["dominant"] = np.where(pst["log2_A_over_B"] > 0, "A", "B")

    # Identify switching cells (have both A and B trips)
    sw = (pst.groupby(["site", "compartment"])
          .agg(unique_dom=("dominant", "nunique"),
               n_trips=("trip", "nunique"))
          .reset_index())
    sw_cells = sw[(sw["unique_dom"] == 2) & (sw["n_trips"] >= 2)]
    print(f"  Switching cells: {len(sw_cells)}", flush=True)

    # For each switching cell, compare environmental conditions in A vs B trips
    merged = pst.merge(sw_cells[["site", "compartment"]],
                          on=["site", "compartment"])
    print(f"  Samples in switching cells: {len(merged)}", flush=True)

    for var in ("d7", "d30", "d90", "d180", "d365"):
        if var not in merged.columns: continue
        A_v = merged.loc[merged["dominant"] == "A", var].dropna()
        B_v = merged.loc[merged["dominant"] == "B", var].dropna()
        if len(A_v) < 5 or len(B_v) < 5: continue
        try:
            U, p_mw = mannwhitneyu(A_v, B_v, alternative="two-sided")
        except Exception:
            continue
        print(f"  {var:<6}  A trips: median={A_v.median():.1f} mm  "
              f"B trips: median={B_v.median():.1f} mm  "
              f"MW p={p_mw:.3g}", flush=True)

    # Also: per-trip switching breakdown — which trips show most B
    print(f"\n  Trip-wise A/B in switching cells:")
    print(merged.groupby(["trip", "dominant"]).size().unstack(fill_value=0)
          .to_string())

    # ==================================================================
    # 5. NIBRIBACTER FUNCTIONAL ARTIFACTS
    # ==================================================================
    print("\n========================================================")
    print("[5] Nibribacter MAG functional artifact check")
    print("========================================================")
    mag_kos = pd.read_csv(CACHE / "nibribacter_mags" /
                            "per_mag_function_full.tsv", sep="\t")
    bet = mag_kos[mag_kos["category"] == "betaine_biosynth"]
    print(f"  betaine_biosynth KOs claimed: {len(bet)}", flush=True)
    print(f"  Unique KOs in betaine_biosynth: {bet['ko'].nunique()}",
          flush=True)
    print(f"  Top 10 KOs assigned to betaine_biosynth:")
    counts = bet.groupby("ko").agg(
        n=("orf", "count"),
        definition=("definition", "first")).reset_index()
    print(counts.sort_values("n", ascending=False).head(10)
          .to_string(index=False), flush=True)
    counts.to_csv(OUT / "nibribacter_betaine_KO_breakdown.tsv",
                   sep="\t", index=False)

    # Check the actual annotations - are they really betaine-specific or
    # broad alcohol-DH?
    print("\n  Are these REAL betaine genes?")
    print(f"    K00108 (betA, choline DH): "
          f"{(bet['ko'] == 'K00108').sum()} ORFs", flush=True)
    print(f"    K00130 (betB, betaine aldehyde DH): "
          f"{(bet['ko'] == 'K00130').sum()} ORFs", flush=True)
    print(f"    K17755 (gbsA): "
          f"{(bet['ko'] == 'K17755').sum()} ORFs", flush=True)
    other = bet[~bet["ko"].isin(["K00108", "K00130", "K17755"])]
    print(f"    OTHER (broad alcohol-DH false positives): "
          f"{len(other)} ORFs", flush=True)
    if len(other):
        print(f"    Top OTHER KOs:")
        print(other.groupby("ko").agg(
            n=("orf", "count"),
            definition=("definition", "first"))
              .sort_values("n", ascending=False).head(5)
              .to_string(), flush=True)


if __name__ == "__main__":
    main()
