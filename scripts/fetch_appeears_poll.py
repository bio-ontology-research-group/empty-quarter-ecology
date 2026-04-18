"""Resume polling an already-submitted AppEEARS task and download results.

Task ID is stored in cache/appeears_task_id.txt from a prior submit.
Uses longer timeouts and retries to survive transient network issues.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import pandas as pd
import requests

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from eq import CACHE_DIR

API = "https://appeears.earthdatacloud.nasa.gov/api"


def login() -> str:
    import netrc
    n = netrc.netrc()
    creds = n.authenticators("urs.earthdata.nasa.gov")
    user, _, pw = creds
    r = requests.post(f"{API}/login", auth=(user, pw), timeout=60)
    r.raise_for_status()
    return r.json()["token"]


def safe_get(url, headers, timeout=60, retries=5):
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=headers, timeout=timeout)
            r.raise_for_status()
            return r
        except Exception as e:
            wait = 2 ** attempt
            print(f"  [retry {attempt+1}] {e} — sleeping {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError(f"failed after {retries} retries: {url}")


def main():
    tid_file = CACHE_DIR / "appeears_task_id.txt"
    if not tid_file.exists():
        raise RuntimeError("no cached task_id; run fetch_appeears_ndvi.py first")
    tid = tid_file.read_text().strip()
    print(f"polling AppEEARS task {tid}", flush=True)

    token = login()
    headers = {"Authorization": f"Bearer {token}"}

    while True:
        r = safe_get(f"{API}/task/{tid}", headers)
        status = r.json().get("status", "?")
        print(f"  [{time.strftime('%H:%M:%S')}] {status}", flush=True)
        if status == "done":
            break
        if status in ("failed", "expired"):
            raise RuntimeError(f"task {status}: {r.json()}")
        time.sleep(30)

    # download results
    r = safe_get(f"{API}/bundle/{tid}", headers)
    files = r.json()["files"]
    out = CACHE_DIR / "appeears_ndvi"
    out.mkdir(exist_ok=True)
    for f in files:
        local = out / f["file_name"]
        if local.exists() and local.stat().st_size > 0:
            continue
        print(f"  downloading {f['file_name']}", flush=True)
        rr = safe_get(f"{API}/bundle/{tid}/{f['file_id']}", headers, timeout=600)
        local.write_bytes(rr.content)
    print(f"downloaded {len(files)} files → {out}")

    # Inspect the point CSV
    for f in files:
        name = f["file_name"]
        if name.endswith(".csv") and "MOD13Q1" in name or "MYD13Q1" in name or "results" in name.lower():
            local = out / name
            df = pd.read_csv(local)
            print(f"\n{name}: shape {df.shape}")
            print("  columns:", list(df.columns)[:12])
            if "Value" in df.columns or "_250m_16_days_NDVI" in "".join(df.columns):
                # NDVI files have Value column (scaled by 10000 for MOD13Q1)
                val_col = [c for c in df.columns if "NDVI" in c or c == "Value"][0]
                print(f"  {val_col}: mean={df[val_col].mean():.4f}, "
                      f"range=[{df[val_col].min():.4f}, {df[val_col].max():.4f}]")


if __name__ == "__main__":
    main()
