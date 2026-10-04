#!/usr/bin/env python3
"""Check total catalogue read recruitment independently of the 975-MAG trait subset."""
import argparse
import hashlib
import json
import tarfile
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverm", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    with tarfile.open(args.coverm) as archive:
        for member in archive.getmembers():
            if not member.name.endswith(".tsv"):
                continue
            data = pd.read_csv(archive.extractfile(member), sep="\t")
            data.columns = ["genome", "relative_abundance_pct"]
            mapped = data.loc[data.genome != "unmapped"]
            rows.append({"library": Path(member.name).stem,
                         "recruited_pct": float(mapped.relative_abundance_pct.sum()),
                         "genomes": len(mapped)})
    data = pd.DataFrame(rows).sort_values("library")
    assert len(data) == 150 and set(data.genomes) == {991}
    args.output.mkdir(parents=True, exist_ok=True)
    data.to_csv(args.output / "library_recruitment.tsv", sep="\t", index=False)
    summary = {"source": "metadata/metagenome/coverm_profiles.tar.gz",
               "source_sha256": hashlib.sha256(args.coverm.read_bytes()).hexdigest(),
               "n_libraries": len(data), "n_genomes": 991,
               "mean_recruited_pct": float(data.recruited_pct.mean()),
               "median_recruited_pct": float(data.recruited_pct.median())}
    (args.output / "recruitment_summary.json").write_text(json.dumps(summary, indent=2)+"\n")


if __name__ == "__main__":
    main()
