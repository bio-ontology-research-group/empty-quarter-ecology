"""Data loading and QC filtering for the Empty Quarter amplicon dataset.

This module reproduces the filtering rules described in the paper methods
(min prevalence 3, min total abundance 50, drop mito/chloroplast, drop
samples with <1000 reads). It reads the raw Ampliseq outputs and writes
analysis-ready parquet caches.

Intended entry points:
    load_feature_table_raw  -> (ASVs x samples) int DataFrame, pre-QC
    load_taxonomy           -> (ASVs x ranks) string DataFrame
    load_site_geodata       -> (site, trip) keyed environment & coordinates
    load_daily_weather      -> (site, date) keyed weather
    build_analysis_dataset  -> filter + join + cache
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import CACHE_DIR, DATA_DIR
from .sample_id import is_control, parse

log = logging.getLogger(__name__)

FT_PATH = DATA_DIR / "taxonomy" / "feature-table-trips1-5.tsv"
TAX_PATH = DATA_DIR / "taxonomy" / "taxonomy-trips1-5.tsv"
FASTA_PATH = DATA_DIR / "taxonomy" / "ASV_seqs-trips1-5.fasta"
XRF_PATH = DATA_DIR / "geochemistry" / "xrf_lab_table_filtered.tsv"
DAILY_WEATHER_PATH = DATA_DIR / "climate" / "daily_weather.tsv"
MONTHLY_WEATHER_PATH = DATA_DIR / "climate" / "monthly_weather_averages.tsv"
GEODATA_DIR = DATA_DIR / "geodata"
PICRUSt2_PATH_ABUN = DATA_DIR / "functional" / "picrust2" / "path_abun_unstrat.tsv"


# ----------------------------------------------------------------------
# Filtering defaults
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class QCConfig:
    min_prevalence: int = 3          # ASV must appear in >= this many samples
    min_total_abundance: int = 50    # ASV total read count across all samples
    min_sample_reads: int = 1000     # Minimum total reads per sample
    drop_mitochondria: bool = True
    drop_chloroplast: bool = True
    drop_archaea: bool = False       # Keep archaea by default
    drop_controls: bool = True
    keep_site_max: int = 60          # Sites 61-64 were special probes (Trip 1 only)
    drop_human_contaminants: bool = True
    # Trip 3 (Winter 2024) had additional expedition participants (tourists);
    # the last-segment sites 47-59 show elevated human gut-associated genera
    # at 5-15% relative abundance. Following Salter et al. 2014 Genome Biol
    # and Eisenhofer et al. 2019 Trends Microbiol, we drop the canonical
    # human-gut / skin contaminant genera before analysis.
    human_contaminant_genera: tuple[str, ...] = (
        "Salmonella", "Escherichia-Shigella", "Escherichia", "Shigella",
        "Streptococcus", "Enterococcus", "Staphylococcus", "Lactobacillus",
        "Corynebacterium", "Cutibacterium", "Propionibacterium",
    )


# ----------------------------------------------------------------------
# Raw loading
# ----------------------------------------------------------------------
def load_feature_table_raw() -> pd.DataFrame:
    """Load the Ampliseq feature table as (ASV x samples) int DataFrame.

    The file starts with a BIOM comment (``# Constructed from biom file``)
    and a header line that begins with ``#OTU ID`` followed by sample
    IDs. We skip the BIOM banner and parse the sample header manually
    so the ``#`` on the real header doesn't collide with pandas'
    ``comment=`` stripping.
    """
    log.info("reading %s", FT_PATH)
    with open(FT_PATH) as fh:
        first = fh.readline()
        if not first.lstrip().startswith("# Constructed from biom file"):
            # Fallback: no BIOM banner; rewind and treat as plain TSV
            fh.seek(0)
            first = ""
        header = fh.readline().rstrip("\n").split("\t")
    # Strip the leading ``#`` from the first column name (``#OTU ID``)
    if header and header[0].lstrip().startswith("#"):
        header[0] = header[0].lstrip("#").strip()
    skiprows = 2 if first else 1
    dtypes = {h: "int64" for h in header[1:]}
    dtypes[header[0]] = "string"
    ft = pd.read_csv(
        FT_PATH,
        sep="\t",
        header=None,
        names=header,
        skiprows=skiprows,
        index_col=0,
        dtype=dtypes,
    )
    if "Taxon" in ft.columns:
        ft = ft.drop(columns=["Taxon"])
    ft.index.name = "ASV"
    log.info("raw feature table: %d ASVs x %d sample columns", *ft.shape)
    return ft


def load_taxonomy() -> pd.DataFrame:
    """Load the taxonomy table and explode the SILVA lineage into ranks."""
    log.info("reading %s", TAX_PATH)
    tax = pd.read_csv(TAX_PATH, sep="\t", index_col=0)
    tax.index.name = "ASV"
    # Ampliseq writes the lineage as ";"-separated, kingdom-first
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    parts = tax["Taxon"].str.split(";", expand=True)
    for i, r in enumerate(ranks):
        if i < parts.shape[1]:
            tax[r] = parts[i].str.strip()
        else:
            tax[r] = ""
    # The final column in Ampliseq SILVA output is a confidence score
    if parts.shape[1] > len(ranks):
        tax["confidence"] = pd.to_numeric(parts[len(ranks)], errors="coerce")
    return tax


def load_site_geodata() -> pd.DataFrame:
    """Concatenate per-trip geodata into one (trip, site)-keyed frame."""
    frames = []
    for trip in range(1, 6):
        p = GEODATA_DIR / f"trip{trip}_geodata.tsv"
        if not p.exists():
            continue
        df = pd.read_csv(p, sep="\t")
        df["trip"] = trip
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out = out.rename(columns={"Site": "site"})
    out["site_raw"] = out["site"].astype(str)
    out["site"] = pd.to_numeric(out["site"], errors="coerce")
    # Drop non-numeric sites (e.g., well names) for the numeric index;
    # they can still be joined via site_raw.
    out = out.dropna(subset=["site"])
    out["site"] = out["site"].astype(int)
    return out.set_index(["trip", "site"])


def load_daily_weather() -> pd.DataFrame:
    w = pd.read_csv(DAILY_WEATHER_PATH, sep="\t", parse_dates=["Date"])
    w = w.rename(columns={"Site": "site", "Date": "date"})
    return w.set_index(["site", "date"])


def load_xrf() -> pd.DataFrame:
    """Load XRF lab table; keep raw sample-level rows (not aggregated)."""
    if not XRF_PATH.exists():
        log.warning("XRF file missing: %s", XRF_PATH)
        return pd.DataFrame()
    xrf = pd.read_csv(XRF_PATH, sep="\t", index_col=0)
    return xrf


def load_picrust_pathways() -> pd.DataFrame:
    """Load PICRUSt2 MetaCyc pathway abundance table."""
    if not PICRUSt2_PATH_ABUN.exists():
        log.warning("PICRUSt2 file missing: %s", PICRUSt2_PATH_ABUN)
        return pd.DataFrame()
    pw = pd.read_csv(PICRUSt2_PATH_ABUN, sep="\t", index_col=0)
    pw.index.name = "pathway"
    return pw


# ----------------------------------------------------------------------
# Filtering + metadata
# ----------------------------------------------------------------------
def build_sample_metadata(sample_ids: list[str]) -> pd.DataFrame:
    """Parse sample IDs and return a DataFrame indexed by the raw IDs."""
    rows = []
    for sid in sample_ids:
        parsed = parse(sid)
        if parsed is None:
            rows.append(
                {"sample": sid, "trip": pd.NA, "site": pd.NA, "compartment": pd.NA,
                 "replicate": pd.NA, "suffix": "", "season": pd.NA, "year": pd.NA,
                 "control": is_control(sid)}
            )
        else:
            rows.append(
                {"sample": sid, "trip": parsed.trip, "site": parsed.site,
                 "compartment": parsed.compartment, "replicate": parsed.replicate,
                 "suffix": parsed.suffix, "season": parsed.season, "year": parsed.year,
                 "control": False}
            )
    df = pd.DataFrame(rows).set_index("sample")
    for c in ("trip", "site", "replicate", "year"):
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    return df


def apply_qc(
    ft: pd.DataFrame, tax: pd.DataFrame, config: QCConfig | None = None
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Apply the QC filters described in the paper methods.

    Returns
    -------
    ft_filt : filtered feature table (ASV x samples)
    tax_filt : filtered taxonomy table (ASV x ranks)
    report : dict of step-by-step counts (for the QC notebook)
    """
    config = config or QCConfig()
    report: dict = {}
    report["start_asvs"] = ft.shape[0]
    report["start_samples"] = ft.shape[1]

    # 1. drop control samples by parsed ID
    if config.drop_controls:
        keep_samples = [c for c in ft.columns if not is_control(c)]
        ft = ft[keep_samples]
    report["after_control_drop_samples"] = ft.shape[1]

    # 2. drop mito/chloroplast/archaea ASVs by taxonomy
    tax_lower = tax["Taxon"].str.lower()
    keep_mask = pd.Series(True, index=tax.index)
    if config.drop_mitochondria:
        keep_mask &= ~tax_lower.str.contains("mitochondria", na=False)
    if config.drop_chloroplast:
        keep_mask &= ~tax_lower.str.contains("chloroplast", na=False)
    if config.drop_archaea:
        keep_mask &= ~tax_lower.str.startswith("archaea", na=False)
    # keep only ASVs present in ft
    keep_mask = keep_mask & keep_mask.index.isin(ft.index)
    asvs_in_tax = set(tax.index[keep_mask])
    ft = ft.loc[ft.index.isin(asvs_in_tax)]
    report["after_taxa_drop_asvs"] = ft.shape[0]

    # 2b. drop human-associated contaminant genera (tourist / lab sources)
    if config.drop_human_contaminants and "Taxon" in tax.columns:
        pat = "|".join(rf"g__{g}(?:;|$)|{g}(?:;|$)" for g in config.human_contaminant_genera)
        contaminant_asvs = tax.index[tax["Taxon"].str.contains(pat, case=False, na=False, regex=True)]
        before = ft.shape[0]
        ft = ft.loc[~ft.index.isin(contaminant_asvs)]
        report["dropped_human_contaminant_asvs"] = before - ft.shape[0]
        report["after_human_contaminant_asvs"] = ft.shape[0]

    # 3. min prevalence and min total abundance
    prev = (ft > 0).sum(axis=1)
    abund = ft.sum(axis=1)
    keep = (prev >= config.min_prevalence) & (abund >= config.min_total_abundance)
    ft = ft.loc[keep]
    report["after_prev_abund_asvs"] = ft.shape[0]

    # 4. drop sites above keep_site_max (e.g., sites 61-64 are Trip 1 only)
    if config.keep_site_max:
        parsed = [(s, parse(s)) for s in ft.columns]
        keep_by_site = [s for s, p in parsed
                        if p is not None and p.site <= config.keep_site_max]
        dropped = [s for s, p in parsed
                   if p is None or p.site > config.keep_site_max]
        ft = ft[keep_by_site]
        report["after_site_cap_samples"] = ft.shape[1]
        report["dropped_for_site_cap"] = len(dropped)

    # 5. min sample reads
    sample_reads = ft.sum(axis=0)
    keep_samples = sample_reads[sample_reads >= config.min_sample_reads].index
    ft = ft[keep_samples]
    report["after_sample_reads_samples"] = ft.shape[1]

    tax_filt = tax.loc[tax.index.isin(ft.index)]

    report["final_asvs"] = ft.shape[0]
    report["final_samples"] = ft.shape[1]
    return ft, tax_filt, report


# ----------------------------------------------------------------------
# Orchestration
# ----------------------------------------------------------------------
def build_analysis_dataset(
    config: QCConfig | None = None, use_cache: bool = True
) -> dict:
    """Run full QC and return analysis-ready tables.

    Writes parquet/tsv caches under ``CACHE_DIR``.
    """
    cache_hit = all(
        (CACHE_DIR / f).exists()
        for f in ("feature_table.parquet", "taxonomy.parquet", "metadata.parquet")
    )
    if use_cache and cache_hit:
        return {
            "feature_table": pd.read_parquet(CACHE_DIR / "feature_table.parquet"),
            "taxonomy": pd.read_parquet(CACHE_DIR / "taxonomy.parquet"),
            "metadata": pd.read_parquet(CACHE_DIR / "metadata.parquet"),
            "qc_report": pd.read_json(CACHE_DIR / "qc_report.json", typ="series").to_dict(),
        }

    ft = load_feature_table_raw()
    tax = load_taxonomy()
    ft_filt, tax_filt, report = apply_qc(ft, tax, config)
    meta = build_sample_metadata(ft_filt.columns.tolist())

    # Persist
    ft_filt.to_parquet(CACHE_DIR / "feature_table.parquet")
    tax_filt.to_parquet(CACHE_DIR / "taxonomy.parquet")
    meta.to_parquet(CACHE_DIR / "metadata.parquet")
    pd.Series(report).to_json(CACHE_DIR / "qc_report.json")

    return {
        "feature_table": ft_filt,
        "taxonomy": tax_filt,
        "metadata": meta,
        "qc_report": report,
    }
