#!/usr/bin/env python3
"""Render the supplementary tables of the taxon-context module as LaTeX.

The fragment is written to ``generated/taxon_context_tables.tex`` in the
manuscript directory and is included by ``supplement.tex``.  Every number in
the fragment comes from ``analysis/v3/taxon_context/*.tsv``; the regression
test ``tests/test_taxon_context.py`` re-renders the fragment and requires it
to be identical to the committed copy.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "analysis/v3/taxon_context"

PHYLA_ROWS = 10
GENERA_ROWS = 15
REPLACEMENT_ROWS = 10


def tex_escape(text: str) -> str:
    return (
        str(text)
        .replace("&", r"\&")
        .replace("%", r"\%")
        .replace("_", r"\_")
        .replace("#", r"\#")
    )


def pct(value: float, digits: int = 1) -> str:
    return f"{100 * float(value):.{digits}f}"


def italic(name: str) -> str:
    return r"\textit{" + tex_escape(name) + "}"


def render_phyla(frame: pd.DataFrame) -> str:
    rows = frame[frame["phylum"] != "unclassified_phylum"].head(PHYLA_ROWS)
    lines = [
        r"\begin{table}[htbp]",
        r"\centering\footnotesize",
        r"\caption{Leading phyla across the $1{,}227$ core-site profiles. Values are"
        r" mean shares of total reads in percent; prevalence is the fraction of"
        r" profiles with at least one read. Compartment and transect-third columns"
        r" are means within the stratum (thirds of $20$ sites ordered by route"
        r" position). Names follow SILVA 138.2.}",
        r"\label{tab:taxa-phyla}",
        r"\begin{tabular}{lrrrrrrrr}",
        r"\toprule",
        r"Phylum & ASVs & Mean & Prev. & Surface & Shallow & Root-adj. & West & East \\",
        r"\midrule",
    ]
    for row in rows.itertuples():
        lines.append(
            f"{tex_escape(row.phylum)} & {row.n_asvs:,} & {pct(row.mean_relative_abundance)} & "
            f"{pct(row.prevalence, 0)} & {pct(row.mean_surface)} & {pct(row.mean_shallow_subsurface)} & "
            f"{pct(row.mean_root_adjacent)} & {pct(row.mean_west_third)} & {pct(row.mean_east_third)} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def render_genera(frame: pd.DataFrame) -> str:
    rows = frame[frame["genus"] != "unclassified_genus"].head(GENERA_ROWS)
    lines = [
        r"\begin{table}[htbp]",
        r"\centering\footnotesize",
        r"\caption{Leading genera across the $1{,}227$ core-site profiles (mean share"
        r" of total reads in percent, prevalence in percent of profiles, and means"
        r" by compartment and transect third).}",
        r"\label{tab:taxa-genera}",
        r"\begin{tabular}{llrrrrrrrr}",
        r"\toprule",
        r"Genus & Phylum & Mean & Prev. & Surface & Shallow & Root-adj. & West & Central & East \\",
        r"\midrule",
    ]
    for row in rows.itertuples():
        lines.append(
            f"{italic(row.genus)} & {tex_escape(row.phylum)} & {pct(row.mean_relative_abundance, 2)} & "
            f"{pct(row.prevalence, 0)} & {pct(row.mean_surface, 2)} & {pct(row.mean_shallow_subsurface, 2)} & "
            f"{pct(row.mean_root_adjacent, 2)} & {pct(row.mean_west_third, 2)} & "
            f"{pct(row.mean_central_third, 2)} & {pct(row.mean_east_third, 2)} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def render_replacement(frame: pd.DataFrame) -> str:
    supported = frame[frame["supported_q_lt_0_05"]]
    decreasing = supported.sort_values("spearman_rho_route_position", kind="mergesort").head(REPLACEMENT_ROWS)
    increasing = supported.sort_values(
        "spearman_rho_route_position", ascending=False, kind="mergesort"
    ).head(REPLACEMENT_ROWS)
    n_dec = int((supported["direction"] == "decreases_eastward").sum())
    n_inc = int((supported["direction"] == "increases_eastward").sum())
    lines = [
        r"\begin{table}[htbp]",
        r"\centering\footnotesize",
        r"\caption{Genera with the strongest monotonic change along the route among the"
        r" $200$ genera of the route model. $\rho$ is the Spearman correlation of the"
        r" site-level CLR value (campaign-by-compartment means removed) with route"
        r" position; $q$ is Benjamini--Hochberg over $200$ tests"
        f" ({n_dec} genera decreased and {n_inc} increased eastward at $q<0.05$)."
        r" $\Delta$ is the east-minus-west difference in mean CLR between the"
        r" transect thirds with a delete-one-site jackknife $95\,\%$ interval;"
        r" the last two columns give the mean share of total reads in percent in"
        r" the western and eastern thirds.}",
        r"\label{tab:taxa-replacement}",
        r"\begin{tabular}{llrrrrr}",
        r"\toprule",
        r"Genus & Phylum & $\rho$ & $q$ & $\Delta$ CLR [95\,\% CI] & West & East \\",
        r"\midrule",
        r"\multicolumn{7}{l}{\emph{Strongest decreases from west to east}} \\",
    ]

    def row_line(row) -> str:
        return (
            f"{italic(row.genus)} & {tex_escape(row.phylum)} & {row.spearman_rho_route_position:.2f} & "
            f"{row.q_bh_200:.1e} & {row.east_minus_west_mean_clr:.2f} "
            f"[{row.east_minus_west_ci_low:.2f}, {row.east_minus_west_ci_high:.2f}] & "
            f"{pct(row.mean_relative_abundance_west_third, 2)} & "
            f"{pct(row.mean_relative_abundance_east_third, 2)} \\\\"
        )

    lines += [row_line(row) for row in decreasing.itertuples()]
    lines.append(r"\midrule")
    lines.append(r"\multicolumn{7}{l}{\emph{Strongest increases from west to east}} \\")
    lines += [row_line(row) for row in increasing.itertuples()]
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def render_gradients(frame: pd.DataFrame) -> str:
    lines = [
        r"\begin{table}[htbp]",
        r"\centering\footnotesize",
        r"\caption{Direction of the environmental gradients along the route:"
        r" Spearman correlation of site means with route position and the means"
        r" of the western and eastern transect thirds.}",
        r"\label{tab:route-gradients}",
        r"\begin{tabular}{lrrrr}",
        r"\toprule",
        r"Variable & Sites & $\rho$ & West third & East third \\",
        r"\midrule",
    ]
    for row in frame.itertuples():
        lines.append(
            f"{tex_escape(row.description)} & {row.n_sites} & {row.spearman_rho_route_position:.2f} & "
            f"{row.mean_west_third:.2f} & {row.mean_east_third:.2f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def render_overlap(frame: pd.DataFrame) -> str:
    lines = [
        r"\begin{table}[htbp]",
        r"\centering\footnotesize",
        r"\caption{Genus-set overlap after subsampling every profile $20$ times to"
        r" $8{,}000$ reads. A genus counts as detected in a set when it occurs in at"
        r" least $10\,\%$ of subsampled profiles. ``Top 50 A in B'' is the number"
        r" of the $50$ most abundant genera of set A that are detected in set B."
        r" The Atacama pit profiles (PRJEB39249) use V4 primers and intracellular"
        r" DNA; the comparison is presence-based and does not merge feature"
        r" tables.}",
        r"\label{tab:genus-overlap}",
        r"\begin{tabular}{p{3.6cm}p{3.6cm}rrrrrrr}",
        r"\toprule",
        r"Set A & Set B & Profiles A & Profiles B & Genera A & Genera B & Shared & Jaccard & Top 50 A in B \\",
        r"\midrule",
    ]
    for row in frame.itertuples():
        lines.append(
            f"{tex_escape(row.set_a)} & {tex_escape(row.set_b)} & {row.n_profiles_a} & {row.n_profiles_b} & "
            f"{row.n_genera_a} & {row.n_genera_b} & {row.n_shared} & {row.jaccard:.3f} & "
            f"{row.top50_a_detected_in_b} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def render_pathways(
    dominance: pd.DataFrame, classes: pd.DataFrame, supported: pd.DataFrame, autotrophy: pd.DataFrame
) -> str:
    label = {
        "biosynthesis": "Biosynthesis",
        "degradation_utilization": "Degradation or utilization",
        "energy_central_metabolism": "Energy and central metabolism",
        "other": "Other",
    }
    lines = [
        r"\begin{table}[htbp]",
        r"\centering\footnotesize",
        r"\caption{Predicted pathway classes. Share is the class total of mean"
        r" relative pathway abundance across the $1{,}227$ core-site profiles"
        r" (keyword grouping of MetaCyc descriptions; rules archived in"
        r" \texttt{analysis/v3/taxon\_context/run\_manifest.json}). The supported"
        r" counts summarize the route-correlation family ($200$ tests) and the"
        r" compartment-contrast family ($600$ tests) of Section~S8.1 by class,"
        r" with the number of positive (eastward increase or first-named"
        r" compartment higher) and negative estimates.}",
        r"\label{tab:pathway-classes}",
        r"\begin{tabular}{lrrrrrr}",
        r"\toprule",
        r"Class & Pathways & Share (\%) & Route (+/$-$) & Compartment (+/$-$) & & \\",
        r"\midrule",
    ]
    for row in classes.itertuples():
        route = supported[(supported["family"] == "route_correlation_supported") & (supported["pathway_class"] == row.pathway_class)]
        comp = supported[(supported["family"] == "compartment_contrast_supported") & (supported["pathway_class"] == row.pathway_class)]
        r_pos = int(route["n_positive"].sum()); r_neg = int(route["n_negative"].sum())
        c_pos = int(comp["n_positive"].sum()); c_neg = int(comp["n_negative"].sum())
        lines.append(
            f"{label[row.pathway_class]} & {row.n_pathways} & {pct(row.share_of_predicted_pathway_abundance)} & "
            f"{r_pos + r_neg} ({r_pos}/{r_neg}) & {c_pos + c_neg} ({c_pos}/{c_neg}) & & \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]

    top = dominance.head(10)
    lines += [
        r"\begin{table}[htbp]",
        r"\centering\footnotesize",
        r"\caption{The ten most abundant predicted MetaCyc pathways and the"
        r" carbon-fixation pathways present in the PICRUSt2 pathway set (mean"
        r" share of predicted pathway abundance in percent; rank among $462$"
        r" pathways).}",
        r"\label{tab:pathway-dominance}",
        r"\begin{tabular}{rllr}",
        r"\toprule",
        r"Rank & Pathway & Description & Share (\%) \\",
        r"\midrule",
    ]
    for row in top.itertuples():
        lines.append(
            f"{row.rank} & {tex_escape(row.pathway)} & {tex_escape(row.description)} & "
            f"{pct(row.mean_relative_abundance, 2)} \\\\"
        )
    lines.append(r"\midrule")
    lines.append(r"\multicolumn{4}{l}{\emph{Carbon-fixation pathways}} \\")
    for row in autotrophy.itertuples():
        lines.append(
            f"{row.rank} & {tex_escape(row.pathway)} & {tex_escape(row.description)} & "
            f"{pct(row.mean_relative_abundance, 2)} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def render(results: Path) -> str:
    phyla = pd.read_csv(results / "phylum_composition.tsv", sep="\t")
    genera = pd.read_csv(results / "genus_composition.tsv", sep="\t")
    replacement = pd.read_csv(results / "transect_replacement.tsv", sep="\t")
    gradients = pd.read_csv(results / "site_gradients.tsv", sep="\t")
    overlap = pd.read_csv(results / "genus_set_overlap.tsv", sep="\t")
    dominance = pd.read_csv(results / "pathway_dominance.tsv", sep="\t")
    classes = pd.read_csv(results / "pathway_class_share.tsv", sep="\t")
    supported = pd.read_csv(results / "supported_pathways_by_class.tsv", sep="\t")
    autotrophy = pd.read_csv(results / "autotrophy_pathways.tsv", sep="\t")
    header = (
        "% Generated by analysis/v3/render_taxon_context_tex.py from\n"
        "% analysis/v3/taxon_context/*.tsv. Do not edit by hand.\n\n"
    )
    return header + "\n".join(
        [
            render_phyla(phyla),
            render_genera(genera),
            render_replacement(replacement),
            render_gradients(gradients),
            render_overlap(overlap),
            render_pathways(dominance, classes, supported, autotrophy),
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "empty-quarter-amplicon/generated/taxon_context_tables.tex",
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render(args.results.resolve()), encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
