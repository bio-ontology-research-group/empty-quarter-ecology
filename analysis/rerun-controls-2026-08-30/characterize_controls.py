#!/usr/bin/env python3
"""Characterise the 14 control libraries denoised on 30 Aug 2026 (Ibex job in ibex/job.id).

For every library: read counts through DADA2, ASV count, the five most abundant
genera, the read share assigned to ZymoBIOMICS mock-community genera (a
positive-standard signature), and the read share in ASVs that also occur in the
canonical biological profiles (at any prevalence and at >= 5 % prevalence).

Inputs (all relative to this directory unless absolute):
  ibex/results/controls_md5.tsv          ASV x library counts, md5-of-sequence ids
  ibex/results/DADA2_stats.tsv           per-library read tracking
  ibex/results/ASV_seqs.fasta            run ASV id -> sequence
  ibex/results/ASV_tax.silva_138_2.tsv   run ASV id -> SILVA 138.2 lineage
  ../../data/processed/taxonomy/taxon-tables/feature-table-trips1-5.tsv (canonical)

Outputs: outputs/control_characterization.tsv and outputs/control_characterization.md
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RES = HERE / "ibex" / "results"
CANON = ROOT / "data/processed/taxonomy/taxon-tables/feature-table-trips1-5.tsv"
OUT = HERE / "outputs"

MOCK_GENERA = {
    "Pseudomonas", "Escherichia-Shigella", "Salmonella", "Lactobacillus",
    "Limosilactobacillus", "Enterococcus", "Staphylococcus", "Listeria", "Bacillus",
}
LIBRARY_META = {
    "TRIP4_SEB": ("M-25-0684_SEB", "Trip 4 extraction blank (SEB)"),
    "TRIP4_EB1": ("M-25-0770_EB1", "Trip 4 extraction blank"),
    "TRIP4_EB2": ("M-25-0771_EB2", "Trip 4 extraction blank"),
    "TRIP4_EB3": ("M-25-0772_EB3", "Trip 4 extraction blank"),
    "TRIP4_EB4": ("M-25-0773_EB4", "Trip 4 extraction blank"),
    "TRIP4_EB5": ("M-25-0774_EB5", "Trip 4 extraction blank"),
    "CTL_ExtractionCtrlPro_Trip1": ("M-25-0929_Extraction-Ctrl-Pro-Trip1", "Trips 1-3 control, role pending (Rund)"),
    "CTL_PCRCtrl_Trip1": ("M-25-0555_PCR-Ctrl-Trip1", "Trips 1-3 control, role pending (Rund)"),
    "CTL_Ctrl1_Trip1": ("M-25-0323_Ctrl-1-Trip1", "Trips 1-3 control, role pending (Rund)"),
    "CTL_Ctrl2": ("M-25-0553_Ctrl-2", "Trips 1-3 control, role pending (Rund)"),
    "CTL_Ctrl3": ("M-25-0554_Ctrl-3", "Trips 1-3 control, role pending (Rund)"),
    "CTL_NegCtrl1_Trip2": ("M-25-0870_Neg-Ctrl-1-Trip2", "Trips 1-3 control, role pending (Rund)"),
    "CTL_NegCtrl2_Trip2": ("M-25-0871_Neg-Ctrl-2-Trip2", "Trips 1-3 control, role pending (Rund)"),
    "CTL_NTC2": ("M-25-0875_NTC-2", "no-template control on the shared run, trip pending (Rund)"),
}


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def md5(seq: str) -> str:
    return hashlib.md5(seq.upper().encode()).hexdigest()


def load_genus_by_md5(seqs: dict[str, str], tax_path: Path) -> dict[str, str]:
    genus = {}
    with tax_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []
        idcol = cols[0]
        gcol = next((c for c in cols if c.lower() == "genus"), None)
        for row in reader:
            seq = seqs.get(row[idcol])
            if seq is None:
                continue
            g = (row.get(gcol) or "").strip() if gcol else ""
            genus[md5(seq)] = g or "NA"
    return genus


def load_counts(path: Path) -> tuple[list[str], dict[str, dict[str, int]]]:
    with path.open() as fh:
        header = None
        table: dict[str, dict[str, int]] = {}
        for line in fh:
            if line.startswith("#") and header is None:
                if line.startswith("#OTU ID"):
                    header = line.rstrip("\n").split("\t")[1:]
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")[1:]
                continue
            parts = line.rstrip("\n").split("\t")
            table[parts[0]] = {c: int(float(v)) for c, v in zip(header, parts[1:]) if float(v) > 0}
    return header, table


def canonical_prevalence(asv_ids: set[str]) -> tuple[dict[str, int], int]:
    """Number of canonical biological profiles in which each ASV occurs."""
    prevalence: dict[str, int] = {}
    with CANON.open() as fh:
        header = None
        for line in fh:
            if header is None:
                if line.startswith("#OTU ID") or not line.startswith("#"):
                    header = line.rstrip("\n").split("\t")[1:]
                continue
            asv, rest = line.split("\t", 1)
            if asv not in asv_ids:
                continue
            values = rest.rstrip("\n").split("\t")
            prevalence[asv] = sum(1 for c, v in zip(header, values) if c in bio_columns_cache and float(v) > 0)
    return prevalence, len(bio_columns_cache)


bio_columns_cache: set[str] = set()


def canonical_columns() -> list[str]:
    with CANON.open() as fh:
        for line in fh:
            if line.startswith("#OTU ID") or not line.startswith("#"):
                return line.rstrip("\n").split("\t")[1:]
    return []


def main() -> int:
    OUT.mkdir(exist_ok=True)
    seqs = read_fasta(RES / "ASV_seqs.fasta")
    genus_of = load_genus_by_md5(seqs, RES / "ASV_tax.silva_138_2.tsv")
    libraries, counts = load_counts(RES / "controls_md5.tsv")
    stats = {r["sample"]: r for r in csv.DictReader((RES / "DADA2_stats.tsv").open(), delimiter="\t")}

    per_lib: dict[str, dict[str, int]] = defaultdict(dict)
    for asv, row in counts.items():
        for lib, v in row.items():
            per_lib[lib][asv] = v
    all_asvs = set(counts)

    global bio_columns_cache
    cols = canonical_columns()
    bio_columns_cache = {c for c in cols if not re.fullmatch(r"(EB\d+|Negative\d+)", c)}
    prevalence, n_bio = canonical_prevalence(all_asvs)

    rows = []
    for lib in libraries:
        asvs = per_lib.get(lib, {})
        total = sum(asvs.values())
        by_genus: Counter = Counter()
        for asv, v in asvs.items():
            by_genus[genus_of.get(asv, "NA")] += v
        top = by_genus.most_common(5)
        mock = sum(v for g, v in by_genus.items() if g in MOCK_GENERA)
        shared_any = sum(v for asv, v in asvs.items() if prevalence.get(asv, 0) >= 1)
        shared_5pct = sum(v for asv, v in asvs.items() if prevalence.get(asv, 0) >= 0.05 * n_bio)
        st = stats.get(lib, {})
        meta = LIBRARY_META.get(lib, ("", ""))
        rows.append({
            "library": lib,
            "sequencing_library": meta[0],
            "description": meta[1],
            "reads_input": st.get("DADA2_input", ""),
            "reads_nonchimeric": st.get("nonchim", total),
            "asvs": len(asvs),
            "top_genera": "; ".join(f"{g} {v / total:.1%}" for g, v in top) if total else "",
            "mock_genus_read_fraction": f"{mock / total:.4f}" if total else "",
            "reads_in_asvs_shared_with_any_biological_profile": f"{shared_any / total:.4f}" if total else "",
            "reads_in_asvs_at_ge5pct_biological_prevalence": f"{shared_5pct / total:.4f}" if total else "",
            "asvs_absent_from_biological_profiles": sum(1 for asv in asvs if prevalence.get(asv, 0) == 0),
        })

    with (OUT / "control_characterization.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with (OUT / "control_characterization.md").open("w") as fh:
        fh.write("# Control libraries denoised 30 Aug 2026 (14 libraries, one DADA2 run)\n\n")
        fh.write(f"Canonical biological profiles used for prevalence: {n_bio}. Mock genera: {', '.join(sorted(MOCK_GENERA))}.\n\n")
        fh.write("| library | sequencing library | reads in | non-chimeric | ASVs | top genera | mock share | reads in ASVs shared with any biological profile | reads in ASVs at >= 5 % prevalence |\n|---|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            fh.write(f"| {r['library']} | {r['sequencing_library']} | {r['reads_input']} | {r['reads_nonchimeric']} | {r['asvs']} | {r['top_genera']} | {r['mock_genus_read_fraction']} | {r['reads_in_asvs_shared_with_any_biological_profile']} | {r['reads_in_asvs_at_ge5pct_biological_prevalence']} |\n")
    json.dump({"libraries": libraries, "n_biological_profiles": n_bio, "n_asvs_total": len(all_asvs)},
              (OUT / "control_characterization.json").open("w"), indent=2)
    print(f"wrote {OUT / 'control_characterization.tsv'} ({len(rows)} libraries, {len(all_asvs)} ASVs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
