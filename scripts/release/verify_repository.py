#!/usr/bin/env python3
"""Fail-closed verification of the ecology reproducibility repository."""

from __future__ import annotations

import csv
import hashlib
import re
import subprocess
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_file_manifest(root: Path, failures: list[str]) -> None:
    manifest = root / "FILE_MANIFEST.tsv"
    if not manifest.is_file():
        failures.append("missing FILE_MANIFEST.tsv")
        return
    expected: dict[str, tuple[int, str]] = {}
    with manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["path", "bytes", "sha256"]:
            failures.append("FILE_MANIFEST.tsv has an unexpected header")
            return
        for row in reader:
            expected[row["path"]] = (int(row["bytes"]), row["sha256"])

    for relative, (size, digest) in expected.items():
        path = root / relative
        if path.is_symlink() or not path.is_file():
            failures.append(f"missing or non-regular tracked file: {relative}")
        elif path.stat().st_size != size:
            failures.append(f"wrong byte count: {relative}")
        elif sha256(path) != digest:
            failures.append(f"wrong SHA-256: {relative}")

    output = subprocess.check_output(
        ["git", "-C", str(root), "ls-files", "-z", "--cached"]
    )
    tracked = {
        item.decode("utf-8")
        for item in output.split(b"\0")
        if item and item.decode("utf-8") != "FILE_MANIFEST.tsv"
    }
    for relative in sorted(tracked):
        lower = relative.lower()
        if lower.endswith((".nxml", ".html")):
            failures.append(
                f"copyrighted source text must not be tracked: {relative}"
            )
        if lower.endswith(".pdf") and (
            relative.startswith(("literature/", "review-literature/"))
            or "/repository/resources/pdfs/" in relative
        ):
            failures.append(f"literature PDF must not be tracked: {relative}")
    listed = set(expected)
    for relative in sorted(tracked - listed):
        failures.append(f"tracked file absent from manifest: {relative}")
    for relative in sorted(listed - tracked):
        failures.append(f"manifest path is not tracked: {relative}")


def verify_sha256s(root: Path, failures: list[str]) -> None:
    for checksum_file in root.rglob("SHA256SUMS"):
        if any(part in {".git", "data", "data-paper"} for part in checksum_file.parts):
            continue
        for line in checksum_file.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            parts = line.split(maxsplit=1)
            if len(parts) != 2:
                failures.append(f"malformed checksum line: {checksum_file}")
                continue
            digest, filename = parts
            target = checksum_file.parent / filename.lstrip("*")
            if not target.is_file():
                failures.append(f"checksum target missing: {target.relative_to(root)}")
            elif sha256(target) != digest:
                failures.append(f"stale checksum: {target.relative_to(root)}")


def verify_manuscript(root: Path, failures: list[str]) -> None:
    paper = root / "empty-quarter-amplicon"
    main = (paper / "main.tex").read_text(encoding="utf-8")
    supplement = (paper / "supplement.tex").read_text(encoding="utf-8")
    combined = main + supplement
    for marker in (r"\todo{", r"\paragraph{", r"\subparagraph{"):
        if marker in combined:
            failures.append(f"active manuscript contains {marker}")
    if re.search(r"\\(?:input|include|subfile)\s*\{", main):
        failures.append("main.tex contains an included TeX fragment")
    for required in (
        "main.pdf",
        "supplement.pdf",
        "figures/fig1_landscape.pdf",
        "figures/fig2_soil_position.pdf",
        "figures/fig3_function_controls.pdf",
        "figures/figS_campaign_rainfall.pdf",
        "figures/figure_manifest.tsv",
        "figures/figure_runtime.json",
    ):
        if not (paper / required).is_file():
            failures.append(f"missing manuscript artifact: {required}")

    figure_manifest = paper / "figures/figure_manifest.tsv"
    with figure_manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if row["role"] != "output":
            continue
        path = paper / "figures" / row["file"]
        if not path.is_file():
            failures.append(f"missing figure output: {row['file']}")
        elif path.stat().st_size != int(row["bytes"]):
            failures.append(f"wrong figure byte count: {row['file']}")
        elif sha256(path) != row["sha256"]:
            failures.append(f"wrong figure SHA-256: {row['file']}")


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    failures: list[str] = []
    verify_file_manifest(root, failures)
    verify_sha256s(root, failures)
    verify_manuscript(root, failures)
    if failures:
        print("FAIL: ecology repository verification", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1
    print("PASS: repository, checksums, manuscript, and figures verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
