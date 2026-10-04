import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis/v3"))
from render_biology_context_tex import marker_description, scientific, render
from render_taxon_context_tex import render as render_taxon


def test_scientific_notation_preserves_precision_and_math_context():
    assert scientific(2.3e-12) == r"$2.3\times10^{-12}$"
    assert scientific(6.02e-10, digits=2, math_mode=False) == r"6.02\times10^{-10}"


def test_gene_symbols_are_distinct_from_protein_names():
    assert marker_description("UvrABC system protein A (uvrA)") == r"UvrABC system protein A (\textit{uvrA})"
    assert marker_description("RuBisCO large subunit (rbcL/cbbL)") == r"RuBisCO large subunit (\textit{rbcL}/\textit{cbbL})"
    assert marker_description("Dps protein; dps gene") == r"Dps protein; \textit{dps} gene"


def test_rendered_tables_use_scientific_q_values_and_gene_italics():
    base = ROOT / "analysis/v3"
    biology, traits = render(base / "biology_context_corrected_20260909",
                             base / "trait_genes_20260909",
                             base / "ph_context_reanalysis_20260909/ph_genus_correlations.tsv")
    taxon = render_taxon(base / "taxon_context_corrected_20260909")
    assert r"q=6.02\times10^{-10}" in biology
    assert r"$2.3\times10^{-12}$" in taxon
    assert r"\textit{hoxH}" in traits
    assert "UvrABC system protein A" in traits
