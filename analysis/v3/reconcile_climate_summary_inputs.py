#!/usr/bin/env python3
"""Record why corrected context runs remain valid after the climate-table rerun.

``environment_associations/climate_site_summary.tsv`` was regenerated after the
Site 52 coordinate correction. Only its ``latitude``, ``longitude`` and
``transect_km`` columns changed. The September 2026 biology and taxon context
runs read the earlier table but consume none of those columns. This program
does not rerun them and does not alter their recorded input digest. It
verifies that every consumed column is identical between the table the run
read (located by its recorded SHA-256) and the current table, then appends an
``input_reconciliation`` record to the run manifest and refreshes the
manifest's line in SHA256SUMS. It fails if any consumed value differs.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import pandas as pd


INPUT = "analysis/v3/environment_associations/climate_site_summary.tsv"
CLIMATE = ["mean_air_temperature_c", "mean_monthly_rain_mm", "mean_relative_humidity_pct"]
ALPHA = ["shannon", "expected_richness_25k", "normalized_evenness"]
CONSUMERS = {
    # Columns read by biology_context.py (landform alpha summaries,
    # dune/pan Mann-Whitney tests, dune-only climate-diversity correlations).
    "analysis/v3/biology_context_corrected_20260909": ["site", *ALPHA, *CLIMATE],
    # Columns read by taxon_context.py (route gradients of the climate means).
    "analysis/v3/taxon_context_corrected_20260909": ["site", *CLIMATE],
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def reconcile(root: Path, previous: Path) -> list[dict[str, object]]:
    current_bytes = (root / INPUT).read_bytes()
    previous_bytes = previous.read_bytes()
    current = pd.read_csv(io.BytesIO(current_bytes), sep="\t")
    earlier = pd.read_csv(io.BytesIO(previous_bytes), sep="\t")
    changed = sorted(
        column
        for column in earlier.columns
        if column not in current or not earlier[column].equals(current[column])
    )
    records = []
    for directory, columns in CONSUMERS.items():
        manifest_path = root / directory / "run_manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        recorded = manifest["inputs"]["climate_site_summary"]["sha256"]
        if recorded != sha256_bytes(previous_bytes):
            raise ValueError(f"{previous} is not the table read by {directory}")
        if not earlier[columns].equals(current[columns]):
            raise ValueError(f"Consumed columns differ for {directory}; rerun it")
        record = {
            "input": INPUT,
            "sha256_read_by_run": recorded,
            "sha256_current": sha256_bytes(current_bytes),
            "reason": (
                "environment_associations rerun on 2026-10-04 with the corrected "
                "Site 52 coordinates (data commit 5a17782); only the listed "
                "changed columns differ"
            ),
            "changed_columns": changed,
            "consumed_columns": columns,
            "consumed_columns_identical": True,
            "verification": (
                "pandas DataFrame.equals on the consumed columns of the table "
                "read by this run and the current table (60 rows); this run "
                "was not repeated"
            ),
            "verified_by": "analysis/v3/reconcile_climate_summary_inputs.py",
        }
        manifest["input_reconciliation"] = {"climate_site_summary": record}
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        sums_path = root / directory / "SHA256SUMS"
        digest = sha256_bytes(manifest_path.read_bytes())
        lines = [
            f"{digest}  run_manifest.json" if line.endswith("  run_manifest.json") else line
            for line in sums_path.read_text(encoding="utf-8").splitlines()
        ]
        sums_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        records.append({"directory": directory, **record})
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument(
        "--previous-table",
        type=Path,
        required=True,
        help="Copy of the climate_site_summary.tsv read by the September runs",
    )
    args = parser.parse_args()
    for record in reconcile(args.project_root.resolve(), args.previous_table.resolve()):
        print(record["directory"], "consumed columns identical; changed:", record["changed_columns"])


if __name__ == "__main__":
    main()
