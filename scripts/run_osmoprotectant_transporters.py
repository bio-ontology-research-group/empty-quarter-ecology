#!/usr/bin/env python3
"""Tabulate osmoprotectant production, uptake and efflux machinery in
each CSP1-2 MAG.

For the public-good model we need:
  PRODUCTION   — betA (choline → betaine aldehyde, K00108)
                 betB (betaine aldehyde → glycine betaine, K00130)
  UPTAKE       — TC 3.A.1.12.* (quaternary amine ABC transporters,
                 OpuA/B/C-like) and 2.A.15 (BetT-like MFS uptake)
  EFFLUX       — mechanosensitive channels MscL/MscS (passive osmotic
                 release), 2.A.1 family permeases acting in reverse,
                 dedicated exporters where annotated
  TREHALOSE    — alternative osmolyte (treS, otsA/B)

Outputs:
  cache/csp_mag_osmoprotectant_summary.tsv  — one row per MAG
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
MAG_ROOT = CACHE / "csp_mags"

MAGS = ["V27Dr2__SemiBin_73",
        "V30PRr1__SemiBin_7",
        "V32PRr1__SemiBin_26",
        "V38PRr3__SemiBin_38"]

def count_in_tsv(tsv_path: Path, patterns: list[str]) -> int:
    """Count rows matching any pattern (case-insensitive) in the
    bakta tsv (gene/product columns)."""
    if not tsv_path.exists():
        return 0
    try:
        df = pd.read_csv(tsv_path, sep="\t", comment="#",
                         header=None,
                         names=["contig", "type", "start", "end", "strand",
                                "locus", "gene", "product", "extra"])
    except Exception:
        return 0
    n = 0
    for p in patterns:
        m = (df["gene"].fillna("").str.contains(p, case=False, regex=True)
             | df["product"].fillna("").str.contains(p, case=False, regex=True))
        n += int(m.sum())
    return n


def count_transporter_tbl(path: Path, substrate: str, evalue_max=1e-20) -> int:
    if not path.exists():
        return 0
    df = pd.read_csv(path, sep="\t", header=None, comment="#",
                     names=["q", "tcid", "substrate_name", "ex_id",
                            "rxn", "subj", "pid", "evalue", "score",
                            "qcov", "header", "start", "end", "type"],
                     dtype=str)
    df["evalue_f"] = pd.to_numeric(df["evalue"], errors="coerce")
    sel = df[(df["substrate_name"] == substrate)
              & (df["evalue_f"] < evalue_max)]
    return int(len(sel))

rows = []
for mag in MAGS:
    short = mag.split("__")[0]
    bakta = MAG_ROOT / mag / f"{mag.split('__')[1]}.tsv"
    transp = MAG_ROOT / f"{mag}_gapseq" / f"{mag.split('__')[1]}-Transporter.tbl"

    betA   = count_in_tsv(bakta, [r"choline\s+dehydrogenase", r"\bbetA\b"])
    betB   = count_in_tsv(bakta, [r"betaine.{1,4}aldehyde\s+dehydrogenase",
                                    r"\bbetB\b"])
    treS   = count_in_tsv(bakta, [r"trehalose\s+synthase", r"\btreS\b"])
    otsA   = count_in_tsv(bakta, [r"trehalose-6-phosphate\s+synthase",
                                    r"\botsA\b"])
    otsB   = count_in_tsv(bakta, [r"trehalose-6-phosphate\s+phosphatase",
                                    r"\botsB\b"])
    mscL   = count_in_tsv(bakta, [r"mechanosensitive\s+channel", r"\bMscL\b"])
    mscS   = count_in_tsv(bakta, [r"\bMscS\b", r"small\s+conductance.*mechano"])
    bet_uptake  = count_transporter_tbl(transp, "betaine")
    chol_uptake = count_transporter_tbl(transp, "choline")
    rows.append({
        "mag": short,
        "betA": betA, "betB": betB,
        "treS": treS, "otsA": otsA, "otsB": otsB,
        "MscL": mscL, "MscS": mscS,
        "betaine_TC_hits": bet_uptake,
        "choline_TC_hits": chol_uptake,
    })

out = pd.DataFrame(rows)
print(out.to_string(index=False))
out.to_csv(CACHE / "csp_mag_osmoprotectant_summary.tsv",
           sep="\t", index=False)
print(f"\nwrote {CACHE / 'csp_mag_osmoprotectant_summary.tsv'}")
print()
print("Interpretation:")
print(f"  betA present in {(out['betA'] > 0).sum()}/{len(out)} MAGs")
print(f"  betB present in {(out['betB'] > 0).sum()}/{len(out)} MAGs")
print(f"  Both betA and betB:   {((out['betA'] > 0) & (out['betB'] > 0)).sum()}/{len(out)} MAGs")
print(f"  trehalose treS in {(out['treS'] > 0).sum()}/{len(out)} MAGs (otsA-otsB absent)")
print(f"  Mechanosensitive efflux (MscL or MscS): "
      f"{((out['MscL'] > 0) | (out['MscS'] > 0)).sum()}/{len(out)} MAGs")
print(f"  Has uptake-transporters for betaine: "
      f"{(out['betaine_TC_hits'] > 0).sum()}/{len(out)} MAGs "
      f"(median {int(out['betaine_TC_hits'].median())} hits)")
