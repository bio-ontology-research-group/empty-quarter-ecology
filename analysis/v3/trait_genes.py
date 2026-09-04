#!/usr/bin/env python3
"""Pre-specified trait-gene screen from two independent gene sources.

The manuscript predicts metabolic potential with PICRUSt2 and validates it
against shotgun-derived KEGG Ortholog (KO) profiles.  This module asks a
narrower, biologically framed question for a fixed list of trait genes that
matter in hyperarid soil: atmospheric trace-gas oxidation, carbon fixation,
nitrogen cycling, osmolyte synthesis, pigmentation, sporulation, DNA repair
and oxidative-stress defence.  Two sources are used and never merged:

1. The genome catalogue of the companion shotgun study: eggNOG-mapper KO
   annotations of every predicted protein in the 990-genome analytical
   subset, and CoverM relative abundance of each genome in 150 shotgun
   libraries.  For each trait we report the fraction of genomes that carry
   at least one member KO, the phylum breakdown of the carriers, and the
   abundance-weighted share of the recruited community that carries the
   trait in each library, by compartment.
2. PICRUSt2 predicted KO abundance in the 1,227 core-site profiles.  For each
   trait we report the mean predicted share of KO abundance by compartment,
   paired site-level compartment contrasts (log10 share; sign-flip tests,
   Benjamini-Hochberg over traits x 3 contrasts), and the Spearman
   correlation with route position.

KO membership is a fixed list declared in this file (TRAITS).  Hydrogenase
KOs cannot resolve the high-affinity group 1h enzyme from other group 1
[NiFe] hydrogenases, so that trait is labelled "[NiFe]-hydrogenase large
subunit (group 1/1h KOs)".  Presence of a gene is potential, not activity.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import re
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import taxon_context as tc  # noqa: E402

SCHEMA_VERSION = "1.0"
POSITIONS = ("Surface", "Deep", "Rhizosphere")
LABEL = {"Surface": "surface", "Deep": "shallow_subsurface", "Rhizosphere": "root_adjacent"}
CONTRASTS = (("Deep", "Surface"), ("Rhizosphere", "Surface"), ("Rhizosphere", "Deep"))

# trait -> (category, description, KO list)
TRAITS: dict[str, tuple[str, str, tuple[str, ...]]] = {
    "coxL_CO_dehydrogenase": ("trace gas", "aerobic CO dehydrogenase large subunit (coxL)", ("K03520",)),
    "NiFe_hydrogenase_large": ("trace gas", "[NiFe]-hydrogenase large subunit (group 1/1h KOs hyaB, hhyL, hoxH)", ("K06281", "K18008", "K00436")),
    "rbcL_RuBisCO": ("carbon fixation", "RuBisCO large subunit (rbcL/cbbL)", ("K01601",)),
    "psbA_photosystem_II": ("carbon fixation", "photosystem II D1 protein (psbA)", ("K02703",)),
    "nifH_nitrogenase": ("nitrogen", "nitrogenase iron protein (nifH)", ("K02588",)),
    "amoA_pmoA_ammonia_monooxygenase": ("nitrogen", "ammonia/methane monooxygenase subunit A (amoA/pmoA)", ("K10944",)),
    "nirK_nirS_nitrite_reductase": ("nitrogen", "dissimilatory nitrite reductase (nirK, nirS)", ("K00368", "K15864")),
    "nosZ_N2O_reductase": ("nitrogen", "nitrous-oxide reductase (nosZ)", ("K00376",)),
    "ureC_urease": ("nitrogen", "urease alpha subunit (ureC)", ("K01428",)),
    "otsA_otsB_trehalose": ("osmolyte", "trehalose synthesis (otsA, otsB)", ("K00697", "K01087")),
    "treY_treZ_trehalose": ("osmolyte", "trehalose from maltooligosaccharides (treY, treZ)", ("K06044", "K01236")),
    "ectABC_ectoine": ("osmolyte", "ectoine synthesis (ectA, ectB, ectC)", ("K06718", "K00836", "K06720")),
    "crtB_phytoene_synthase": ("pigment", "phytoene synthase (crtB), carotenoid pigments", ("K02291",)),
    "spo0A_sporulation": ("dormancy", "sporulation master regulator (spo0A)", ("K07699",)),
    "phrB_photolyase": ("DNA repair", "deoxyribodipyrimidine photolyase (phrB)", ("K01669",)),
    "uvrA_excision_repair": ("DNA repair", "UvrABC system protein A (uvrA)", ("K03701",)),
    "recA_recombination_repair": ("DNA repair", "recombinase A (recA)", ("K03553",)),
    "katE_katG_catalase": ("oxidative stress", "catalase (katE) and catalase-peroxidase (katG)", ("K03781", "K03782")),
    "sodA_superoxide_dismutase": ("oxidative stress", "Fe/Mn superoxide dismutase (sodA)", ("K04564",)),
    "dps_DNA_protection": ("oxidative stress", "DNA protection during starvation protein (dps)", ("K04047",)),
}
GENOME_RE = re.compile(r"^(?P<genome>.+?)@@")
KO_RE = re.compile(r"K\d{5}")
SAMPLE_RE = re.compile(r"^(?P<prefix>[TFSV])?(?P<site>\d+)(?P<comp>PR|P|D|S)r(?P<rep>\d+)$")
PREFIX_CAMPAIGN = {"": 1, "T": 2, "F": 3, "S": 4, "V": 5}
CODE = {"D": "Deep", "S": "Surface", "P": "Rhizosphere", "PR": "Rhizosphere"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, float_format="%.10g", lineterminator="\n")


def parse_library(name: str) -> dict | None:
    m = SAMPLE_RE.match(name)
    if not m:
        return None
    return {"library": name, "campaign": PREFIX_CAMPAIGN[m.group("prefix") or ""], "site": int(m.group("site")), "position": CODE[m.group("comp")]}


def genome_ko_presence(annotations: Path, all_kos: set[str]) -> tuple[dict[str, set[str]], dict[str, int]]:
    """Genome -> set of trait KOs present; genome -> number of annotated proteins."""
    presence: dict[str, set[str]] = {}
    proteins: dict[str, int] = {}
    with gzip.open(annotations, "rt", encoding="utf-8") as handle:
        header = None
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.rstrip("\n").split("\t")
                ko_index = header.index("KEGG_ko")
                continue
            fields = line.rstrip("\n").split("\t")
            match = GENOME_RE.match(fields[0])
            if not match:
                continue
            genome = match.group("genome")
            proteins[genome] = proteins.get(genome, 0) + 1
            kos = fields[ko_index]
            if kos == "-":
                continue
            hits = {k for k in KO_RE.findall(kos) if k in all_kos}
            if hits:
                presence.setdefault(genome, set()).update(hits)
    return presence, proteins


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=HERE.parents[1])
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--permutations", type=int, default=9999)
    args = parser.parse_args()
    root = args.project_root.resolve()
    output = (args.output_dir or root / "analysis/v3/trait_genes").resolve()
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    inputs = {
        "eggnog_annotations": root / "data/metadata/metagenome/eq.emapper.annotations.gz",
        "coverm_profiles": root / "data/metadata/metagenome/coverm_profiles.tar.gz",
        "measured_function_inputs": root / "data/metadata/metagenome/measured_function_inputs.tar.gz",
        "picrust_ko": root / "data/processed/functional/picrust2/merged/ko_metagenome_unstrat.tsv",
        "genus_counts": root / "analysis/v2/review/cache/genus_counts.tsv",
        "site_coordinates": root / "analysis/v3/spatial_turnover_rescue/results/site_coordinates.tsv",
    }
    missing = [str(p) for p in inputs.values() if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing inputs (install the bulk artifacts):\n" + "\n".join(missing))

    all_kos = {k for _, _, kos in TRAITS.values() for k in kos}
    trait_table = pd.DataFrame(
        [{"trait": t, "category": c, "description": d, "kos": ",".join(k)} for t, (c, d, k) in TRAITS.items()]
    )

    # ---- genome catalogue --------------------------------------------------
    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(inputs["measured_function_inputs"]) as tar:
            member = "analysis/v2/review/measured_function/filtered_genomes.tsv"
            tar.extract(member, tmp)
            genomes = pd.read_csv(Path(tmp) / member, sep="\t")
        with tarfile.open(inputs["coverm_profiles"]) as tar:
            tar.extractall(tmp)
            coverm_dir = Path(tmp) / "coverm"
            coverm = {}
            for path in sorted(coverm_dir.glob("*.tsv")):
                frame = pd.read_csv(path, sep="\t")
                frame.columns = ["genome", "relative_abundance_pct"]
                coverm[path.stem] = frame.set_index("genome")["relative_abundance_pct"]
    # The genome manifest appends "sta" to the bin names used by CoverM and
    # eggNOG-mapper (for example F25PRr3_MAGScoT_cleanbin_000002sta versus
    # F25PRr3_MAGScoT_cleanbin_000002); strip it to join the three sources.
    genomes["genome"] = genomes["genome"].str.replace(r"sta$", "", regex=True)
    catalogue = set(genomes["genome"])
    lineage = genomes.set_index("genome")["lineage"].fillna("")
    phylum = lineage.map(lambda s: next((p[3:] for p in s.split(";") if p.startswith("p__")), "unclassified"))

    presence, proteins = genome_ko_presence(inputs["eggnog_annotations"], all_kos)
    annotated_genomes = set(proteins)
    analysed = sorted(catalogue & annotated_genomes)
    genome_rows = []
    for g in analysed:
        row = {"genome": g, "phylum": phylum.get(g, "unclassified"), "n_annotated_proteins": proteins[g]}
        for trait, (_, _, kos) in TRAITS.items():
            row[trait] = bool(presence.get(g, set()) & set(kos))
        genome_rows.append(row)
    genome_matrix = pd.DataFrame(genome_rows)
    write_tsv(genome_matrix, output / "genome_trait_presence.tsv")

    summary_rows = []
    for trait, (category, description, kos) in TRAITS.items():
        carriers = genome_matrix[genome_matrix[trait]]
        by_phylum = carriers["phylum"].value_counts()
        summary_rows.append(
            {
                "trait": trait,
                "category": category,
                "description": description,
                "kos": ",".join(kos),
                "n_genomes_with_trait": int(len(carriers)),
                "fraction_of_genomes": len(carriers) / len(genome_matrix),
                "leading_carrier_phyla": "; ".join(f"{p} {n}" for p, n in by_phylum.head(4).items()),
            }
        )
    genome_summary = pd.DataFrame(summary_rows)

    # Abundance-weighted share of recruited genomes carrying each trait.
    libs = []
    for lib, series in coverm.items():
        meta = parse_library(lib)
        if meta is None:
            continue
        series = series.drop(index="unmapped", errors="ignore")
        series = series.reindex(analysed).fillna(0.0)
        total = float(series.sum())
        if total <= 0:
            continue
        row = {**meta, "recruited_pct": total}
        for trait in TRAITS:
            has = genome_matrix.set_index("genome")[trait].reindex(analysed).to_numpy()
            row[trait] = float(series.to_numpy()[has].sum() / total)
        libs.append(row)
    library_shares = pd.DataFrame(libs)
    unparsed = sorted(lib for lib in coverm if parse_library(lib) is None)
    write_tsv(library_shares, output / "library_trait_shares.tsv")
    genome_by_comp = []
    for trait in TRAITS:
        rec = {"trait": trait}
        for position in POSITIONS:
            block = library_shares[library_shares["position"] == position]
            rec[f"mean_share_{LABEL[position]}"] = float(block[trait].mean()) if len(block) else np.nan
            rec[f"n_libraries_{LABEL[position]}"] = int(len(block))
        rec["mean_share_all_libraries"] = float(library_shares[trait].mean())
        genome_by_comp.append(rec)
    genome_summary = genome_summary.merge(pd.DataFrame(genome_by_comp), on="trait")
    write_tsv(genome_summary, output / "genome_trait_summary.tsv")

    # ---- PICRUSt2 predicted KO shares ----------------------------------------
    genus = pd.read_csv(inputs["genus_counts"], sep="\t", index_col=0, nrows=1)
    profiles = [m for sid in genus.columns if (m := tc.sample_metadata(sid)) is not None]
    core = pd.DataFrame(profiles)
    core = core[core["site"].between(1, 60) & core["position"].isin(POSITIONS)]
    header = pd.read_csv(inputs["picrust_ko"], sep="\t", nrows=0).columns.tolist()
    core_ids = [sid for sid in core["sample_id"] if sid in header]
    totals = None
    trait_sums = {t: None for t in TRAITS}
    for chunk in pd.read_csv(inputs["picrust_ko"], sep="\t", index_col=0, usecols=[header[0], *core_ids], chunksize=2000):
        values = chunk.to_numpy(dtype=float)
        totals = values.sum(axis=0) if totals is None else totals + values.sum(axis=0)
        for trait, (_, _, kos) in TRAITS.items():
            mask = chunk.index.isin(kos)
            part = values[mask].sum(axis=0)
            trait_sums[trait] = part if trait_sums[trait] is None else trait_sums[trait] + part
    shares = pd.DataFrame({t: trait_sums[t] / totals for t in TRAITS}, index=core_ids)
    shares.index.name = "sample_id"
    meta = core.set_index("sample_id").loc[core_ids]
    shares["campaign"] = meta["campaign"].to_numpy()
    shares["site"] = meta["site"].to_numpy()
    shares["position"] = meta["position"].to_numpy()
    grouped = shares.groupby(["campaign", "site", "position"])[list(TRAITS)].mean()
    write_tsv(grouped.reset_index(), output / "picrust_trait_shares_grouped.tsv")

    coordinates = pd.read_csv(inputs["site_coordinates"], sep="\t")
    km = coordinates.set_index("site")["transect_km"]
    pic_rows = []
    contrast_rows = []
    site_means = grouped.groupby("site")[list(TRAITS)].mean()
    for trait, (category, description, kos) in TRAITS.items():
        rec = {"trait": trait, "category": category, "description": description, "kos": ",".join(kos)}
        for position in POSITIONS:
            rec[f"mean_predicted_share_{LABEL[position]}"] = float(shares.loc[shares["position"] == position, trait].mean())
        rec["mean_predicted_share_all"] = float(shares[trait].mean())
        rho, p = stats.spearmanr(km.reindex(site_means.index), site_means[trait])
        rec["spearman_rho_route"] = float(rho)
        rec["p_route"] = float(p)
        pic_rows.append(rec)
        logged = np.log10(grouped[trait].clip(lower=1e-12))
        for first, second in CONTRASTS:
            frame = logged.reset_index()
            wide = frame.pivot_table(index=["campaign", "site"], columns="position", values=trait)
            wide = wide.dropna(subset=[first, second])
            diff = (wide[first] - wide[second]).groupby(level="site").mean()
            d = diff.to_numpy()
            n = d.size
            signs = rng.choice((-1.0, 1.0), size=(args.permutations, n))
            null = signs @ d / n
            p_flip = (1 + int((np.abs(null) >= abs(d.mean())).sum())) / (args.permutations + 1)
            contrast_rows.append(
                {
                    "trait": trait,
                    "category": category,
                    "contrast": f"{LABEL[first]}-{LABEL[second]}",
                    "n_sites": int(n),
                    "mean_log10_ratio": float(d.mean()),
                    "fold_change": float(10 ** d.mean()),
                    "sign_flip_p": p_flip,
                }
            )
    picrust_summary = pd.DataFrame(pic_rows)
    picrust_summary["q_bh_route"] = tc.benjamini_hochberg(picrust_summary["p_route"].to_numpy())
    contrasts = pd.DataFrame(contrast_rows)
    contrasts["q_bh"] = tc.benjamini_hochberg(contrasts["sign_flip_p"].to_numpy())
    contrasts["supported_q_lt_0_05"] = contrasts["q_bh"] < 0.05
    write_tsv(picrust_summary, output / "picrust_trait_summary.tsv")
    write_tsv(contrasts, output / "picrust_trait_compartment_contrasts.tsv")

    # ---- agreement between the two sources in matched libraries ---------
    agree_rows = []
    lib_index = library_shares.set_index("library")
    matched = [lib for lib in lib_index.index if any(sid.endswith("_" + lib) or sid == lib for sid in core_ids)]
    sid_for = {}
    for lib in matched:
        sid_for[lib] = next(sid for sid in core_ids if sid.endswith("_" + lib) or sid == lib)
    for trait in TRAITS:
        x = lib_index.loc[matched, trait].to_numpy()
        y = shares.loc[[sid_for[l] for l in matched], trait].to_numpy()
        if np.std(x) == 0 or np.std(y) == 0:
            rho, p = np.nan, np.nan
        else:
            rho, p = stats.spearmanr(x, y)
        agree_rows.append({"trait": trait, "n_matched_libraries": len(matched), "spearman_rho_genome_vs_picrust": float(rho), "p_value": float(p)})
    write_tsv(pd.DataFrame(agree_rows), output / "source_agreement.tsv")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "status": "trait_gene_screen_complete",
        "scope": "Fixed trait KO list; two gene sources reported separately (genome catalogue with CoverM abundance; PICRUSt2 predictions); presence and predicted share are potential, not activity.",
        "genomes": {"catalogue": len(catalogue), "annotated": len(annotated_genomes), "analysed": len(analysed)},
        "libraries": {"coverm": len(coverm), "parsed": int(len(library_shares)), "unparsed": unparsed, "by_position": library_shares["position"].value_counts().to_dict(), "matched_to_picrust": len(matched)},
        "picrust_profiles": len(core_ids),
        "parameters": {"seed": args.seed, "permutations": args.permutations, "traits": {t: list(k) for t, (_, _, k) in TRAITS.items()}},
        "inputs": {name: {"path": str(p.relative_to(root)), "bytes": p.stat().st_size, "sha256": sha256(p)} for name, p in sorted(inputs.items())},
        "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
    }
    (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    write_tsv(trait_table, output / "trait_definitions.tsv")

    readme = ["# Trait-gene screen", "", f"- Genomes analysed: {len(analysed)} (catalogue {len(catalogue)}); libraries {len(library_shares)}; PICRUSt2 profiles {len(core_ids)}", "", "## Genome catalogue: fraction of genomes carrying the trait, abundance-weighted share by compartment", ""]
    for r in genome_summary.itertuples():
        readme.append(f"- {r.trait}: {100 * r.fraction_of_genomes:.1f} % of genomes; share surface {100 * r.mean_share_surface:.1f} %, shallow {100 * r.mean_share_shallow_subsurface:.1f} %, root {100 * r.mean_share_root_adjacent:.1f} %; carriers: {r.leading_carrier_phyla}")
    readme += ["", "## PICRUSt2 predicted share of KO abundance (x1e4) by compartment; route rho", ""]
    for r in picrust_summary.itertuples():
        readme.append(f"- {r.trait}: surface {1e4 * r.mean_predicted_share_surface:.2f}, shallow {1e4 * r.mean_predicted_share_shallow_subsurface:.2f}, root {1e4 * r.mean_predicted_share_root_adjacent:.2f}; route rho {r.spearman_rho_route:+.2f} (q {r.q_bh_route:.2g})")
    readme += ["", "## Supported PICRUSt2 compartment contrasts (BH q<0.05)", ""]
    for r in contrasts[contrasts["supported_q_lt_0_05"]].itertuples():
        readme.append(f"- {r.trait} {r.contrast}: fold {r.fold_change:.2f} (q {r.q_bh:.2g})")
    readme += ["", "## Agreement genome-derived vs PICRUSt2 (matched libraries)", ""]
    for r in pd.DataFrame(agree_rows).itertuples():
        readme.append(f"- {r.trait}: rho {r.spearman_rho_genome_vs_picrust:+.2f}")
    readme += ["", "## Permitted wording", "", "- Gene presence and predicted share are metabolic potential; no activity, rate or expression is measured.", "- The hydrogenase trait pools group 1 KOs and cannot single out the high-affinity group 1h enzyme.", "- The genome catalogue covers the recruited fraction of 150 shotgun libraries, mostly root-adjacent soil; compartment shares from it are descriptive.", ""]
    (output / "README.md").write_text("\n".join(readme), encoding="utf-8")
    lines = [f"{sha256(p)}  {p.name}" for p in sorted(output.iterdir()) if p.is_file() and p.name != "SHA256SUMS"]
    (output / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(readme))


if __name__ == "__main__":
    main()
