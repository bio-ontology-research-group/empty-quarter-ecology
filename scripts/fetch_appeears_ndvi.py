"""Submit AppEEARS point-based task for MODIS Terra/Aqua NDVI at 60 sites.

AppEEARS (https://appeears.earthdatacloud.nasa.gov/) extracts
MODIS time-series at point coordinates — much lighter than full
granule download. Uses Earthdata bearer-token auth via ~/.netrc.

Products:
- MOD13Q1.061 (Terra, 16-day, 250m)   variable: _250m_16_days_NDVI
- MYD13Q1.061 (Aqua,  16-day, 250m)   variable: _250m_16_days_NDVI
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from eq import CACHE_DIR

GEODATA = REPO / "data" / "geodata"
API = "https://appeears.earthdatacloud.nasa.gov/api"


def login() -> str:
    """Get AppEEARS bearer token from ~/.netrc Earthdata credentials."""
    # AppEEARS uses Earthdata login; POST to /login with basic-auth from netrc
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
        df = df[(df.SiteNum >= 1) & (df.SiteNum <= 60) & (df.SiteNum == df.SiteNum.astype(int))]
        frames.append(df[["SiteNum", "Latitude", "Longitude"]])
    return (
        pd.concat(frames).groupby("SiteNum")[["Latitude", "Longitude"]].mean()
        .reset_index().sort_values("SiteNum").reset_index(drop=True)
    )


def submit_task(token: str, sites: pd.DataFrame,
                start: str = "01-01-2022", end: str = "09-30-2025") -> str:
    """Submit an AppEEARS point-extract task. Returns task_id."""
    coords = [
        {"id": f"site_{int(r.SiteNum):02d}",
         "latitude": float(r.Latitude),
         "longitude": float(r.Longitude),
         "category": "desert"}
        for _, r in sites.iterrows()
    ]
    layers = [
        {"product": "MOD13Q1.061", "layer": "_250m_16_days_NDVI"},
        {"product": "MOD13Q1.061", "layer": "_250m_16_days_EVI"},
        {"product": "MOD13Q1.061", "layer": "_250m_16_days_VI_Quality"},
        {"product": "MYD13Q1.061", "layer": "_250m_16_days_NDVI"},
        {"product": "MYD13Q1.061", "layer": "_250m_16_days_EVI"},
    ]
    task = {
        "task_type": "point",
        "task_name": "eq_ndvi_60sites_2022_2025",
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
    print(f"  submitted task_id = {tid}")
    return tid


def wait_for_task(token: str, tid: str, poll_sec: int = 30) -> None:
    headers = {"Authorization": f"Bearer {token}"}
    while True:
        r = requests.get(f"{API}/task/{tid}", headers=headers, timeout=30)
        r.raise_for_status()
        j = r.json()
        status = j.get("status", "?")
        print(f"  [{time.strftime('%H:%M:%S')}] task {tid[:8]}: {status}",
              flush=True)
        if status == "done":
            return
        if status in ("failed", "expired"):
            raise RuntimeError(f"AppEEARS task {status}: {j}")
        time.sleep(poll_sec)


def download_results(token: str, tid: str, out_dir: Path) -> list[Path]:
    headers = {"Authorization": f"Bearer {token}"}
    r = requests.get(f"{API}/bundle/{tid}", headers=headers, timeout=60)
    r.raise_for_status()
    files = r.json()["files"]
    # prefer the MOD13Q1-MYD13Q1 merged CSV
    paths = []
    for f in files:
        if not (f["file_name"].endswith(".csv") or
                f["file_name"].endswith(".json")):
            continue
        url = f"{API}/bundle/{tid}/{f['file_id']}"
        local = out_dir / f["file_name"]
        if not local.exists():
            print(f"  downloading {f['file_name']}", flush=True)
            rr = requests.get(url, headers=headers, timeout=600)
            rr.raise_for_status()
            local.write_bytes(rr.content)
        paths.append(local)
    return paths


def main():
    sites = pooled_sites()
    print("authenticating ...")
    token = login()
    print("  OK")

    print(f"submitting AppEEARS task for {len(sites)} sites ...")
    tid = submit_task(token, sites)
    (CACHE_DIR / "appeears_task_id.txt").write_text(tid)

    print("waiting for task to complete (polling every 30s)...")
    wait_for_task(token, tid)

    out = CACHE_DIR / "appeears_ndvi"
    out.mkdir(exist_ok=True)
    files = download_results(token, tid, out)
    print(f"downloaded {len(files)} files")
    for f in files:
        print(f"  {f.name}  ({f.stat().st_size/1024:.0f} KB)")

    # Read the point CSV; AppEEARS names it <task>-<product>-results.csv
    for f in files:
        if "results" in f.name.lower() and f.name.endswith(".csv"):
            df = pd.read_csv(f)
            print(f"{f.name}: {df.shape}, columns: {list(df.columns)[:8]}")
    print("done.")


if __name__ == "__main__":
    main()
