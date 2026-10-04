#!/usr/bin/env python3
"""Build review PDFs against an extracted, explicitly identified Overleaf baseline.

Requires latexdiff, pdflatex and bibtex. The caller supplies a clean baseline
export (git archive a5d5623 from the Overleaf repository) and a built current
paper. Included tables are flattened; the bibliography is rendered from the
current source and its exact source changes belong in the accompanying patch.
"""
import argparse
import shutil
import subprocess
from pathlib import Path


def run(command, cwd, log):
    with log.open("a") as handle:
        subprocess.run(command, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    args.output = args.output.resolve()
    old, new = args.output / "baseline", args.output / "current"
    shutil.copytree(args.baseline, old, dirs_exist_ok=True)
    shutil.copytree(args.paper, new, dirs_exist_ok=True, ignore=shutil.ignore_patterns(".git", "archive"))
    for name in ("main", "supplement"):
        for directory in (old, new):
            text = (directory / (name + ".tex")).read_text()
            text = text.replace(r"\bibliography{sample}", "% Bibliography rendered separately for review.")
            (directory / (name + "-comparison.tex")).write_text(text)
        destination = new / (name + "-redline.tex")
        log = args.output / (name + "-redline-build.log")
        with destination.open("w") as output, log.open("w") as errors:
            subprocess.run(["latexdiff", "--flatten", "--encoding=utf8", "--type=UNDERLINE",
                            "--graphics-markup=none", str(old / (name + "-comparison.tex")),
                            str(new / (name + "-comparison.tex"))], stdout=output, stderr=errors, check=True)
        text = destination.read_text().replace(r"\begin{document}",
                "\\providecommand{\\lt}{<}\n\\providecommand{\\gt}{>}\n\\begin{document}", 1)
        # Deleted text includes unresolved references from the baseline itself.
        # Show their source labels explicitly instead of unexplained question marks.
        reference_fallback = r"""
\makeatletter
\let\revieworiginalref\ref
\renewcommand{\ref}[1]{\ifcsname r@#1\endcsname\revieworiginalref{#1}\else\texttt{[old ref]}\fi}
\makeatother
"""
        text = text.replace(r"\begin{document}", r"\begin{document}" + reference_fallback, 1)
        text = text.replace(r"\end{document}", "\\bibliography{sample}\n\\end{document}")
        # The baseline cited a nonexistent key; preserve it as an explicit
        # deleted-source diagnostic rather than inventing a bibliography entry.
        text = text.replace(r"\cite{coverm}", r"\texttt{[old cite]}")
        destination.write_text(text)
        stem = name + "-redline"
        run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", destination.name], new, log)
        run(["bibtex", stem], new, log)
        for _ in range(2):
            run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", destination.name], new, log)
        shutil.copy2(new / (stem + ".pdf"), args.output / (stem + ".pdf"))
        shutil.copy2(destination, args.output / destination.name)
    print("Built main-redline.pdf and supplement-redline.pdf")


if __name__ == "__main__":
    main()
