"""Select the active manuscript explicitly without modifying historical copies.

For a checkout with manuscripts elsewhere, set ECOLOGY_MANUSCRIPT_DIR before
pytest. The packaged default is retained for self-contained release checkouts.
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAPER = Path(os.environ.get("ECOLOGY_MANUSCRIPT_DIR", ROOT / "empty-quarter-amplicon")).resolve()
if not (PAPER / "main.tex").is_file() or not (PAPER / "supplement.tex").is_file():
    raise FileNotFoundError(f"Active ecology manuscript unavailable: {PAPER}; set ECOLOGY_MANUSCRIPT_DIR")

# Companion data repository (pinned by DATA_REPOSITORY.lock). The packaged
# ``data`` link points into it; EQ_DATA_REPO overrides for other checkouts.
_DATA_LINK = ROOT / "data"
DATA_REPO = Path(
    os.environ.get(
        "EQ_DATA_REPO",
        _DATA_LINK.resolve().parent if _DATA_LINK.exists() else ROOT.parent / "data-paper",
    )
).resolve()


def locked_data_commit() -> str:
    for line in (ROOT / "DATA_REPOSITORY.lock").read_text(encoding="utf-8").splitlines():
        field, _, value = line.partition("\t")
        if field == "commit":
            return value.strip()
    raise KeyError("commit missing from DATA_REPOSITORY.lock")


def locked_data_bytes(relative: str) -> bytes:
    """Bytes of a data-repository file at the commit pinned by DATA_REPOSITORY.lock.

    Reading through git keeps the check tied to the locked release even when
    the local data checkout is behind or ahead of it.
    """
    import subprocess

    return subprocess.run(
        ["git", "-C", str(DATA_REPO), "show", f"{locked_data_commit()}:{relative}"],
        check=True,
        capture_output=True,
    ).stdout
