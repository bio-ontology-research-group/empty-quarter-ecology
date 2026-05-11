"""Submit AppEEARS point-extract task for MAIAC AOD at 60 EQ sites.

MCD19A2.061 = MAIAC daily 1km Combined Terra+Aqua aerosol product.
Layers requested:
  Optical_Depth_055   AOD at 550 nm (the standard quantitative variable)
  AOD_Uncertainty     pixel-level uncertainty
  AOD_QA              QA flags (filter by valid bits later)

Earthdata bearer-token auth via ~/.netrc (urs.earthdata.nasa.gov).
Polling and result download done by fetch_appeears_poll.py separately.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import requests

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from eq import CACHE_DIR

GEODATA = REPO / "data" / "geodata"
API = "https://appeears.earthdatacloud.nasa.gov/api"


def login() -> str:
    import netrc
    n = netrc.netrc()
    creds = n.authenticators("urs.earthdata.nasa.gov")
    if not creds:
        raise RuntimeError("No Earthdata credentials in ~/.netrc")
    user, _, pw = creds
    r = requests.post(f"{API}/login", auth=(user, pw), timeout=30)
    r.raise_for_status()
    return r.json()["token"]


def pooled_sites() -> pd.DataFrame:
    frames = []
    for t in (1, 2, 3, 4, 5):
        df = pd.read_csv(GEODATA / f"trip{t}_geodata.tsv", sep="\t")
        df["SiteNum"] = pd.to_numeric(df["Site"], errors="coerce")
        df = df.dropna(subset=["SiteNum"]).drop_duplicates("SiteNum")
        df = df[(df.SiteNum >= 1) & (df.SiteNum <= 60)
                & (df.SiteNum == df.SiteNum.astype(int))]
        frames.append(df[["SiteNum", "Latitude", "Longitude"]])
    return (pd.concat(frames).groupby("SiteNum")[["Latitude", "Longitude"]].mean()
            .reset_index().sort_values("SiteNum").reset_index(drop=True))


def submit_task(token: str, sites: pd.DataFrame,
                start: str = "01-01-2022", end: str = "12-31-2025") -> str:
    coords = [
        {"id": f"site_{int(r.SiteNum):02d}",
         "latitude": float(r.Latitude),
         "longitude": float(r.Longitude),
         "category": "desert"}
        for _, r in sites.iterrows()
    ]
    layers = [
        {"product": "MCD19A2.061", "layer": "Optical_Depth_055"},
        {"product": "MCD19A2.061", "layer": "AOD_Uncertainty"},
        {"product": "MCD19A2.061", "layer": "AOD_QA"},
    ]
    task = {
        "task_type": "point",
        "task_name": "eq_aod_60sites_2022_2025",
        "params": {
            "dates": [{"startDate": start, "endDate": end}],
            "layers": layers,
            "coordinates": coords,
        },
    }
    headers = {"Authorization": f"Bearer {token}"}
    r = requests.post(f"{API}/task", json=task, headers=headers, timeout=60)
    r.raise_for_status()
    tid = r.json()["task_id"]
    print(f"  submitted task_id = {tid}", flush=True)
    return tid


def main():
    out_dir = CACHE_DIR / "appeears_aod"
    out_dir.mkdir(exist_ok=True, parents=True)
    sites = pooled_sites()
    print(f"Submitting AppEEARS MAIAC AOD task for {len(sites)} sites", flush=True)
    token = login()
    tid = submit_task(token, sites)
    (out_dir / "task_id.txt").write_text(tid)
    print(f"\nTask submitted. Poll later with: "
          f"scripts/fetch_appeears_poll.py {tid} appeears_aod", flush=True)


if __name__ == "__main__":
    main()
