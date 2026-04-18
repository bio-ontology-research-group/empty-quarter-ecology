"""Parse Empty Quarter sample identifiers.

Canonical form (post-QIIME):
    [eNNNN_]<trip_prefix><site><compartment>r<rep>[<suffix>]

Trip prefixes:
    ""  -> trip 1 (spring 2023)
    "T" -> trip 2 (summer 2023)
    "F" -> trip 3 (winter 2024)
    "S" -> trip 4 (summer 2024)
    "V" -> trip 5 (autumn 2025)

Compartments:
    "S"  -> surface
    "D"  -> deep
    "PR" -> rhizosphere (0--3 cm from roots)
"""

from __future__ import annotations

import re
from dataclasses import dataclass

TRIP_PREFIX_MAP = {"": 1, "T": 2, "F": 3, "S": 4, "V": 5}
TRIP_SEASON_MAP = {1: "spring", 2: "summer", 3: "winter", 4: "summer", 5: "autumn"}
TRIP_YEAR_MAP = {1: 2023, 2: 2023, 3: 2024, 4: 2024, 5: 2025}
COMPARTMENT_MAP = {"S": "surface", "D": "deep", "PR": "rhizosphere", "P": "rhizosphere"}

_SAMPLE_RE = re.compile(
    r"^(?:e\d+_)?([TFSV]?)(\d+)(PR|[SD])r(\d+)(O|T|RE|R)?$"
)


@dataclass(frozen=True, slots=True)
class SampleID:
    raw: str
    trip: int
    site: int
    compartment: str
    replicate: int
    suffix: str
    season: str
    year: int

    @property
    def canonical(self) -> str:
        """Return the form without the Ampliseq eNNNN_ prefix."""
        prefix = {v: k for k, v in TRIP_PREFIX_MAP.items()}[self.trip]
        comp = {"surface": "S", "deep": "D", "rhizosphere": "PR"}[self.compartment]
        return f"{prefix}{self.site}{comp}r{self.replicate}{self.suffix}"


def parse(sid: str) -> SampleID | None:
    """Parse a sample identifier. Returns ``None`` on mismatch."""
    m = _SAMPLE_RE.match(sid)
    if m is None:
        return None
    prefix, site, comp, rep, suffix = m.groups()
    trip = TRIP_PREFIX_MAP.get(prefix or "", 0)
    if trip == 0:
        return None
    return SampleID(
        raw=sid,
        trip=trip,
        site=int(site),
        compartment=COMPARTMENT_MAP[comp],
        replicate=int(rep),
        suffix=suffix or "",
        season=TRIP_SEASON_MAP[trip],
        year=TRIP_YEAR_MAP[trip],
    )


def is_control(sid: str) -> bool:
    """Return True for negative/extraction-blank controls and mock communities."""
    u = sid.upper()
    return (
        u.startswith("EB")
        or u.startswith("NEGATIVE")
        or u.startswith("MOCK")
        or u.startswith("ZYMO")
        or u.startswith("PCR_NEG")
    )
