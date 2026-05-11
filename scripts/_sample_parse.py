"""Sample-ID parsing utilities for the EQ amplicon dataset.

Sample IDs in feature_table.parquet look like 'eNNNN_<token>'.
The token encodes (trip, site, compartment, replicate):

  Trip 1: NO leading letter        e.g. 10Dr2  -> trip 1, site 10, deep,  rep 2
  Trip 2: leading 'T'              e.g. T3Dr2  -> trip 2, site  3, deep,  rep 2
  Trip 3: leading 'F'              e.g. F20Sr3 -> trip 3, site 20, surf,  rep 3
  Trip 4: leading 'S'              e.g. S58Dr1 -> trip 4, site 58, deep,  rep 1
  Trip 5: leading 'V'              e.g. V27Dr2 -> trip 5, site 27, deep,  rep 2

Compartment letters:
  D    -> deep
  S    -> surface
  PR   -> rhizosphere

Trip 5 has additional suffixes (O/T/R/RE) for extraction method, not technical
replicate. Those still parse but the 'replicate' captured here is the rep number
that precedes any suffix.

Source: Public/software/empty-quarter/CLAUDE.md.
"""
from __future__ import annotations

import re

PREFIX_TO_TRIP = {"": 1, "T": 2, "F": 3, "S": 4, "V": 5}
COMP_MAP = {"D": "deep", "S": "surface", "PR": "rhizosphere"}

_TOKEN_RE = re.compile(r"^([A-Z]*)([0-9]+)([A-Z]+)r([0-9]+)([A-Za-z]*)$")


def parse_sample(sample_id: str) -> dict | None:
    """Parse 'eNNNN_<token>' or just '<token>' into trip/site/compartment/rep."""
    tok = sample_id.split("_")[-1] if "_" in sample_id else sample_id
    m = _TOKEN_RE.match(tok)
    if not m:
        return None
    prefix, site, comp_code, rep, suffix = m.groups()
    trip = PREFIX_TO_TRIP.get(prefix)
    if trip is None:
        return None
    return {
        "sample": sample_id,
        "trip": trip,
        "site": int(site),
        "comp_code": comp_code,
        "compartment": COMP_MAP.get(comp_code, "?"),
        "replicate": int(rep),
        "suffix": suffix or None,
        "token": tok,
    }


def parse_samples_to_df(sample_ids):
    import pandas as pd
    rows = [parse_sample(s) for s in sample_ids]
    rows = [r for r in rows if r is not None]
    return pd.DataFrame(rows)
