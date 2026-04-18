"""Pull plant observations from the iNaturalist 'rub-al-khali' project and
match them to rhizosphere sampling sites by coordinate proximity + date window."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from pyinaturalist import get_observations

from . import CACHE_DIR

log = logging.getLogger(__name__)

PROJECT_SLUG = "rub-al-khali"
OUT = CACHE_DIR / "inat_observations.tsv"


def _identification_to_row(obs: dict) -> dict:
    """Flatten a single iNaturalist observation into a row dict."""
    taxon = obs.get("taxon") or {}
    loc = obs.get("geojson", {}).get("coordinates") or [None, None]
    return {
        "inat_id": obs.get("id"),
        "observed_on": obs.get("observed_on_details", {}).get("date"),
        "latitude": loc[1],
        "longitude": loc[0],
        "taxon_id": taxon.get("id"),
        "scientific_name": taxon.get("name"),
        "common_name": taxon.get("preferred_common_name"),
        "rank": taxon.get("rank"),
        "iconic_taxon_name": taxon.get("iconic_taxon_name"),
        "quality_grade": obs.get("quality_grade"),
        "observer": obs.get("user", {}).get("login"),
    }


def fetch_project_observations(force: bool = False) -> pd.DataFrame:
    """Fetch all observations from the rub-al-khali iNaturalist project."""
    if OUT.exists() and not force:
        log.info("loading cached iNat obs from %s", OUT)
        return pd.read_csv(OUT, sep="\t")

    log.info("fetching iNat project %s …", PROJECT_SLUG)
    rows = []
    page = 1
    per_page = 200
    while True:
        resp = get_observations(
            project_id=PROJECT_SLUG,
            per_page=per_page,
            page=page,
            iconic_taxa="Plantae",
        )
        results = resp.get("results", [])
        if not results:
            break
        rows.extend(_identification_to_row(o) for o in results)
        log.info("  page %d -> %d (total %d)", page, len(results), len(rows))
        if len(results) < per_page:
            break
        page += 1

    df = pd.DataFrame(rows)
    df["observed_on"] = pd.to_datetime(df["observed_on"], errors="coerce")
    df.to_csv(OUT, sep="\t", index=False)
    log.info("wrote %d observations to %s", len(df), OUT)
    return df


def haversine_km(lat1: np.ndarray, lon1: np.ndarray,
                 lat2: np.ndarray, lon2: np.ndarray) -> np.ndarray:
    """Haversine great-circle distance in km."""
    R = 6371.0
    lat1_r, lon1_r = np.radians(lat1), np.radians(lon1)
    lat2_r, lon2_r = np.radians(lat2), np.radians(lon2)
    dlat = lat2_r - lat1_r
    dlon = lon2_r - lon1_r
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1_r) * np.cos(lat2_r) * np.sin(dlon / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def match_to_sites(
    obs: pd.DataFrame,
    sites: pd.DataFrame,
    trip_dates: dict[int, tuple[str, str]] | None = None,
    radius_km: float = 2.0,
) -> pd.DataFrame:
    """For every (site, trip) pair, collect plant observations within
    ``radius_km`` of the site whose ``observed_on`` falls within the trip
    date window. Returns a row per (site, trip, plant)."""
    if trip_dates is None:
        trip_dates = {
            1: ("2023-03-15", "2023-03-25"),
            2: ("2023-07-10", "2023-07-20"),
            3: ("2024-02-15", "2024-02-25"),
            4: ("2024-08-15", "2024-08-25"),
            5: ("2025-10-10", "2025-10-20"),
        }

    sites = sites.reset_index().dropna(subset=["Latitude", "Longitude"])
    matches = []
    for (_, site_row) in sites.iterrows():
        site = int(site_row["site"])
        trip = int(site_row["trip"])
        start, end = trip_dates.get(trip, (None, None))
        if start is None:
            continue
        d = haversine_km(
            obs["latitude"].to_numpy(),
            obs["longitude"].to_numpy(),
            np.full(len(obs), site_row["Latitude"]),
            np.full(len(obs), site_row["Longitude"]),
        )
        in_space = d <= radius_km
        in_time = (obs["observed_on"] >= start) & (obs["observed_on"] <= end)
        m = in_space & in_time
        if m.any():
            hit = obs.loc[m, ["scientific_name", "common_name", "rank", "inat_id",
                              "observed_on", "latitude", "longitude"]].copy()
            hit["distance_km"] = d[m]
            hit["site"] = site
            hit["trip"] = trip
            matches.append(hit)

    if not matches:
        return pd.DataFrame()
    out = pd.concat(matches, ignore_index=True)
    return out
