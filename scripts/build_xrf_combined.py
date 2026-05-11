#!/usr/bin/env python3
"""Build the unified XRF table covering all five trips.

Inputs (under ``data/geochemistry/``):
    xrf_lab_table_filtered.tsv      Trip 5 only (multi-replicate)
    raw/xrf_lab_table_trips1-4.tsv  Trips 1-4 (single replicate per
                                    site x compartment x trip)

Output (under ``data/geochemistry/``):
    xrf_lab_table_all_trips.tsv     all five trips, with parsed
                                    site/compartment/trip/replicate
                                    columns. Element columns are kept
                                    raw (units = % dry mass for
                                    elements; % oxide for *Ox forms).

The Trip-5 file already encodes its trip via the leading "V"; the
Trips-1-4 file uses the same prefix convention as the 16S sample IDs:
""=T1, "T"=T2, "F"=T3, "S"=T4, "V"=T5 (see src/eq/sample_id.py).
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
GEO = REPO / "data" / "geochemistry"
OUT = GEO / "xrf_lab_table_all_trips.tsv"

T5_PATH = GEO / "xrf_lab_table_filtered.tsv"
T14_PATH = GEO / "raw" / "xrf_lab_table_trips1-4.tsv"

TRIP_PREFIX = {"": 1, "T": 2, "F": 3, "S": 4, "V": 5}
COMPARTMENT_FROM_LETTER = {"S": "Surface", "D": "Deep", "PR": "Rhizosphere"}

# match prefix + site + compartment + replicate; tolerate legacy "Best
# Detection" / "Fastscreening" suffixes by matching only the canonical
# sample-ID head.
_RE = re.compile(r"^([TFSV]?)(\d+)(PR|[SD])r(\d+)$")


def parse_id(sid: str) -> dict:
    m = _RE.match(str(sid).strip())
    if m is None:
        return {"trip": pd.NA, "site": pd.NA, "compartment": pd.NA, "replicate": pd.NA}
    prefix, site, comp, rep = m.groups()
    return {
        "trip": TRIP_PREFIX[prefix],
        "site": int(site),
        "compartment": COMPARTMENT_FROM_LETTER[comp],
        "replicate": int(rep),
    }


def main() -> None:
    if not T5_PATH.exists():
        raise SystemExit(f"missing input: {T5_PATH}")
    if not T14_PATH.exists():
        raise SystemExit(f"missing input: {T14_PATH}")

    t5 = pd.read_csv(T5_PATH, sep="\t")
    t14 = pd.read_csv(T14_PATH, sep="\t")

    # Schema check: same 69 columns
    if list(t5.columns) != list(t14.columns):
        only_t5 = set(t5.columns) - set(t14.columns)
        only_t14 = set(t14.columns) - set(t5.columns)
        raise SystemExit(
            f"column mismatch — T5 only: {only_t5}; T1-4 only: {only_t14}"
        )

    combined = pd.concat([t5, t14], ignore_index=True)

    # Drop instrument-metadata columns; keep elements + oxides.
    instrument_cols = [c for c in ("Material", "Method", "Mode", "Diameter")
                       if c in combined.columns]
    combined = combined.drop(columns=instrument_cols)

    # Parse IDs
    parsed = combined["SampleID"].apply(parse_id).apply(pd.Series)
    combined = pd.concat([combined, parsed], axis=1)

    # Move parsed columns to the front, after SampleID/SoilType
    front = ["SampleID", "SoilType", "trip", "site", "compartment", "replicate"]
    others = [c for c in combined.columns if c not in front]
    combined = combined[front + others]

    # Sanity: SoilType (from Excel) should agree with parsed compartment
    mismatch = combined[combined["SoilType"].astype(str).str.strip()
                        != combined["compartment"].astype(str)]
    if len(mismatch):
        print(f"WARNING: {len(mismatch)} rows where SoilType disagrees with parsed "
              f"compartment; first few:\n{mismatch[front].head().to_string()}")

    n_unparsed = combined["trip"].isna().sum()
    if n_unparsed:
        print(f"WARNING: {n_unparsed} rows could not be parsed; first few:")
        print(combined.loc[combined['trip'].isna(), 'SampleID'].head().to_list())

    print(f"input rows: T5={len(t5)}, T1-4={len(t14)}, combined={len(combined)}")
    print("Per-trip x compartment counts:")
    pivot = (combined.dropna(subset=["trip", "compartment"])
             .groupby(["trip", "compartment"]).size().unstack("compartment").fillna(0).astype(int))
    print(pivot.to_string())
    print()
    print(f"writing {OUT} ({combined.shape[0]} rows x {combined.shape[1]} cols)")
    combined.to_csv(OUT, sep="\t", index=False)


if __name__ == "__main__":
    main()
