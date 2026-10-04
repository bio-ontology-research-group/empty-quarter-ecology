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
