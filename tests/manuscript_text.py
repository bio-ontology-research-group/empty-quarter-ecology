"""Read the active manuscript as rendered claims rather than LaTeX source.

Manuscript claim tests compare numbers and scientific statements with their
source artifacts. They must not depend on how a value is typeset ($0.28\\,\\%$
versus 0.28%), on line breaks, on comments, or on whether a table is written
inline or read with ``\\input``. This module provides:

* ``source(name)``: the LaTeX source with comments removed and every active
  ``\\input``/``\\include`` expanded recursively relative to the manuscript
  directory;
* ``normalize(text)``: a whitespace-flattened rendering in which numeric and
  typographic markup is reduced to plain text (``$``, ``\\,``, ``{,}``,
  ``\\%``, ``^{*}``, ``\\textit{...}`` and similar);
* ``has_number(text, value)``: whether a number appears as a complete token,
  so ``0.28`` does not match ``10.284``.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

from manuscript_paths import PAPER

_INPUT = re.compile(r"\\(?:input|include)\s*\{([^}]+)\}")


def strip_comments(text: str) -> str:
    """Remove LaTeX comments (unescaped % to end of line)."""
    return "\n".join(re.sub(r"(?<!\\)%.*$", "", line) for line in text.splitlines())


def _expand(text: str, base: Path, depth: int = 0) -> str:
    if depth > 10:
        raise RecursionError("\\input nesting deeper than 10 levels")

    def replace(match: re.Match) -> str:
        name = match.group(1).strip()
        path = base / name
        if not path.suffix:
            path = path.with_suffix(".tex")
        if not path.is_file():
            raise FileNotFoundError(f"active manuscript input missing: {path}")
        return _expand(strip_comments(path.read_text(encoding="utf-8")), base, depth + 1)

    return _INPUT.sub(replace, text)


@lru_cache(maxsize=None)
def source(name: str) -> str:
    """Comment-free LaTeX source of ``main`` or ``supplement`` with inputs expanded."""
    path = PAPER / (name if name.endswith(".tex") else f"{name}.tex")
    return _expand(strip_comments(path.read_text(encoding="utf-8")), PAPER)


def body(name: str) -> str:
    text = source(name)
    return text.split(r"\begin{document}", 1)[-1]


_WRAPPERS = re.compile(r"\\(?:textit|emph|textbf|mathrm|text|textrm|mathit|path|texttt|url)\s*\{([^{}]*)\}")


def normalize(text: str) -> str:
    """Flatten whitespace and reduce numeric/typographic markup to plain text."""
    out = text
    out = out.replace(r"\phantom{-}", "")
    for _ in range(3):
        out = _WRAPPERS.sub(r"\1", out)
    replacements = [
        ("{,}", ","),
        (r"\,", ""),
        (r"\ ", " "),
        ("~", " "),
        (r"\%", "%"),
        ("^{*}", "*"),
        ("$^{*}$", "*"),
        (r"\o{}", "ø"),
        ("$", ""),
        ("---", "—"),
        ("--", "–"),
        ("{", ""),
        ("}", ""),
    ]
    for command, symbol in (("geq", "≥"), ("ge", "≥"), ("leq", "≤"), ("le", "≤"),
                            ("times", "×"), ("rho", "ρ"), ("pm", "±"), ("approx", "≈")):
        out = re.sub(r"\\" + command + r"(?![A-Za-z])", symbol, out)
    for old, new in replacements:
        out = out.replace(old, new)
    # A minus sign typed in math ($-0.23$) and an en dash range both become
    # plain characters; collapse whitespace last.
    return " ".join(out.split())


@lru_cache(maxsize=None)
def text(name: str) -> str:
    """Normalized body text of ``main`` or ``supplement``."""
    return normalize(body(name))


def combined() -> str:
    return text("main") + " " + text("supplement")


def abstract() -> str:
    match = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", source("main"), re.S)
    assert match, "abstract environment missing"
    return normalize(match.group(1))


def section(name: str, start: str, end: str | None = None) -> str:
    """Normalized text from the first occurrence of ``start`` up to ``end``."""
    full = text(name)
    i = full.find(start)
    assert i >= 0, f"{start!r} not found in {name}"
    j = full.find(end, i + len(start)) if end else -1
    return full[i:] if j < 0 else full[i:j]


def has_number(haystack: str, value: str) -> bool:
    """True when ``value`` occurs as a whole numeric token in ``haystack``."""
    pattern = r"(?<![\d.])" + re.escape(value) + r"(?![\d])(?!\.\d)"
    return re.search(pattern, haystack) is not None


def fmt(value: float, decimals: int) -> str:
    """Round half away from zero, as authors round by hand."""
    from decimal import ROUND_HALF_UP, Decimal

    quantum = Decimal(1).scaleb(-decimals)
    return str(Decimal(repr(float(value))).quantize(quantum, rounding=ROUND_HALF_UP))


def pct(value: float, decimals: int = 1) -> str:
    return fmt(100 * value, decimals)


def labels(name: str) -> set[str]:
    return set(re.findall(r"\\label\{([^}]+)\}", source(name)))


def refs(name: str) -> set[str]:
    found: set[str] = set()
    for group in re.findall(r"\\(?:ref|eqref|autoref|pageref)\{([^}]+)\}", source(name)):
        found.update(part.strip() for part in group.split(","))
    return found


_RULES = re.compile(r"\\(?:toprule|midrule|bottomrule|hline|addlinespace|cmidrule\([a-z]*\)[\d-]*|cmidrule)\b")


def table_rows(name: str, label: str) -> list[list[str]]:
    """Cells of every data row in the tabular following ``\\label{label}``.

    Rows are split on ``\\\\`` and cells on ``&`` in the normalized text;
    rule commands and ``\\multicolumn`` group headings are dropped.
    """
    full = text(name)
    marker = r"\label" + label
    start = full.find(marker)
    assert start >= 0, f"{label} not found in {name}"
    begin = full.find(r"\begintabular", start)
    end = full.find(r"\endtabular", begin)
    assert begin >= 0 and end > begin, f"tabular for {label} not found"
    body_text = full[begin:end]
    body_text = body_text.split(" ", 1)[1] if " " in body_text else body_text
    rows = []
    for raw in body_text.split(r"\\"):
        cleaned = _RULES.sub(" ", raw).strip()
        if not cleaned or cleaned.startswith(r"\multicolumn") or "&" not in cleaned:
            continue
        rows.append([cell.strip() for cell in cleaned.split("&")])
    return rows


def table_row(name: str, label: str, first_cell: str) -> list[str]:
    for row in table_rows(name, label):
        if row[0] == first_cell:
            return row
    raise AssertionError(f"row {first_cell!r} missing from {label}")
