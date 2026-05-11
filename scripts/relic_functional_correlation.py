#!/usr/bin/env python3
"""Sample-level: correlate PICRUSt2 KO + pathway abundances with sample
relic_frac. Identifies functions that scale up or down as the relic
fraction grows - i.e., functions enriched in the past vs present community.

Outputs:
  cache/relic_population/func_corr_KO.tsv
  cache/relic_population/func_corr_pathway.tsv
  cache/relic_population/func_corr_summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
OUT = CACHE / "relic_population"
PIC = REPO / "data" / "functional" / "picrust2"


def load_pic_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    print(f"  loaded {path.name}: {df.shape}", flush=True)
    return df


def per_sample_relabund(df: pd.DataFrame) -> pd.DataFrame:
    s = df.sum(axis=0).replace(0, 1)
    return df.div(s, axis=1)


def correlate(func_relabund: pd.DataFrame,
                 relic_frac: pd.Series, label: str) -> pd.DataFrame:
    common = list(set(func_relabund.columns).intersection(relic_frac.index))
    print(f"  {label}: aligning {len(common)} samples", flush=True)
    fr = func_relabund[common]
    rf = relic_frac.loc[common]
    rows = []
    for fn in fr.index:
        v = fr.loc[fn].values
        r, p = spearmanr(v, rf.values, nan_policy="omit")
        rows.append({"feature": fn, "spearman_rho": float(r),
                      "p": float(p),
                      "mean_relabund": float(np.nanmean(v)),
                      "max_relabund": float(np.nanmax(v))})
    return pd.DataFrame(rows).sort_values("spearman_rho")


def main():
    print("Loading per-sample relic fraction ...", flush=True)
    rec = pd.read_csv(OUT / "per_sample_relic_fraction.tsv", sep="\t")
    relic_frac = rec.set_index("sample")["relic_frac"]
    alive_frac = rec.set_index("sample")["alive_frac"]
    print(f"  {len(relic_frac)} samples", flush=True)

    print("\nLoading PICRUSt2 KO predictions ...", flush=True)
    ko = load_pic_table(PIC / "metagenome_pred_metagenome_unstrat.tsv")
    ko_rel = per_sample_relabund(ko)

    print("\nLoading PICRUSt2 pathway predictions ...", flush=True)
    path = load_pic_table(PIC / "path_abun_unstrat.tsv")
    path_rel = per_sample_relabund(path)

    print("\nKO correlations with relic_frac ...", flush=True)
    ko_corr = correlate(ko_rel, relic_frac, "KO")
    ko_corr.to_csv(OUT / "func_corr_KO.tsv", sep="\t", index=False)

    print("\nPathway correlations with relic_frac ...", flush=True)
    path_corr = correlate(path_rel, relic_frac, "pathway")
    path_corr.to_csv(OUT / "func_corr_pathway.tsv", sep="\t", index=False)

    # Adjusted p (BH)
    def bh_adjust(p):
        n = len(p)
        idx = np.argsort(p)
        ranked = np.empty(n)
        ranked[idx] = np.arange(1, n + 1)
        return np.minimum(p * n / ranked, 1.0)

    for tbl, name in ((ko_corr, "KO"), (path_corr, "pathway")):
        tbl["q"] = bh_adjust(tbl["p"].values)
        sig = tbl[tbl["q"] < 0.05]
        print(f"\n  {name}: {len(sig)} of {len(tbl)} features q<0.05", flush=True)

    # Pathway descriptions
    desc_path = PIC / "path_abun_unstrat_descriptions.tsv"
    if desc_path.exists():
        desc = pd.read_csv(desc_path, sep="\t", index_col=0)
        if "description" in desc.columns:
            d = desc["description"]
            path_corr = path_corr.merge(
                d.rename("description").reset_index().rename(
                    columns={d.index.name or "index": "feature"}),
                on="feature", how="left")

    print("\n=== Top 15 pathways POSITIVELY correlated with relic_frac ===")
    print(path_corr.sort_values("spearman_rho", ascending=False)
          .head(15)[["feature", "spearman_rho", "p"] +
                       (["description"] if "description" in path_corr.columns
                        else [])]
          .round(4).to_string(index=False))

    print("\n=== Top 15 pathways NEGATIVELY correlated with relic_frac "
          "(= alive-loaded) ===")
    print(path_corr.sort_values("spearman_rho")
          .head(15)[["feature", "spearman_rho", "p"] +
                       (["description"] if "description" in path_corr.columns
                        else [])]
          .round(4).to_string(index=False))

    # KO highlights — focus on stress response, motility, DNA repair, sporulation
    interesting = {
        "K00525": "ribonucleotide reductase aerobic class III",
        "K01580": "glutamate decarboxylase",
        "K00457": "trehalose synthase",
        "K00958": "sulfate adenylyltransferase",
        "K01590": "fructose-bisphosphatase",
        "K05658": "spore coat",
        "K01092": "phosphate acetyltransferase",
        "K07315": "RecA",
        "K03286": "outer membrane porin",
        "K02014": "iron uptake",
        "K01428": "urease alpha",
        "K07466": "DNA helicase RecQ",
    }
    print("\n=== Selected KO functions ===")
    for kid in interesting:
        sub = ko_corr[ko_corr["feature"] == kid]
        if not sub.empty:
            r = sub["spearman_rho"].iloc[0]
            p = sub["p"].iloc[0]
            print(f"  {kid} ({interesting[kid]}): rho={r:+.3f} p={p:.3g}")

    with open(OUT / "func_corr_summary.txt", "w") as fh:
        fh.write("Functional correlation with sample relic_frac\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Samples used: {len(set(ko_rel.columns).intersection(relic_frac.index))}\n\n")
        fh.write("--- Pathways most POSITIVELY correlated with relic_frac ---\n")
        cols = ["feature", "spearman_rho", "p", "q"]
        if "description" in path_corr.columns:
            cols += ["description"]
        fh.write(path_corr.sort_values("spearman_rho", ascending=False)
                  .head(20)[cols].round(4).to_string(index=False))
        fh.write("\n\n--- Pathways most NEGATIVELY correlated (alive-loaded) "
                  "---\n")
        fh.write(path_corr.sort_values("spearman_rho")
                  .head(20)[cols].round(4).to_string(index=False))
    print(f"\nWrote {OUT}/func_corr_summary.txt", flush=True)


if __name__ == "__main__":
    main()
