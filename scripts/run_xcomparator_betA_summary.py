#!/usr/bin/env python3
"""Cross-comparator betA assay summary.

Reads:
  cache/{Flavisolibacter,Rubellimicrobium,Telluribacter,Solirubrobacter}_K00108.tbl
  cache/xcomparator_betA_aln.faa

Writes:
  cache/xcomparator_betA_summary.tsv
  cache/xcomparator_betA_per_genome.tsv
"""
from __future__ import annotations
from pathlib import Path
import re
import pandas as pd
from Bio import AlignIO

CACHE = Path(__file__).resolve().parents[1] / "cache"

KOFAM_THRESHOLD = 697.33  # K00108 (Kofam ko_list)

GENERA = ["Flavisolibacter","Rubellimicrobium","Telluribacter","Solirubrobacter"]

def parse_tbl(path):
    rows = []
    if not path.exists(): return pd.DataFrame()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            parts = line.split()
            if len(parts) < 6: continue
            try:
                evalue = float(parts[4])
                bitscore = float(parts[5])
            except ValueError:
                continue
            rows.append({"protein": parts[0], "evalue": evalue, "bitscore": bitscore})
    return pd.DataFrame(rows)

# Load alignment, find His473
aln = AlignIO.read(str(CACHE/"xcomparator_betA_aln.faa"), "fasta")
ref_idx = next(i for i,r in enumerate(aln) if "REF_Ecoli" in r.id)
ref = str(aln[ref_idx].seq)
ungapped = 0
his_col = None
for ap, ch in enumerate(ref):
    if ch != "-":
        ungapped += 1
        if ungapped == 473:
            his_col = ap; break
his_residues = {}
for r in aln:
    if r.id == aln[ref_idx].id: continue
    name = r.id
    his_residues[name] = str(r.seq)[his_col]

# For each genus, classify hits
all_hits = []
for g in GENERA:
    tbl = parse_tbl(CACHE / f"{g}_K00108.tbl")
    if tbl.empty:
        all_hits.append({"genus": g, "n_hits_above_KOfam": 0,
                         "n_hits_GMC_paralog_only": 0,
                         "best_bitscore": None, "His473_conserved_top5": None,
                         "n_genomes_with_K00108": 0})
        continue
    tbl["genus"] = g
    tbl["above_threshold"] = tbl["bitscore"] >= KOFAM_THRESHOLD
    # extract genome contig prefix from protein ID (NCBI accession pattern: NZ_XXXX.1_n)
    tbl["genome"] = tbl["protein"].str.extract(r"^([A-Z]+_[A-Z0-9]+)", expand=False)
    n_above = int(tbl["above_threshold"].sum())
    n_below_but_e30 = int((~tbl["above_threshold"]).sum())
    best = float(tbl["bitscore"].max()) if len(tbl) else None
    n_genomes_above = tbl[tbl["above_threshold"]]["genome"].nunique()

    # His473 of top 5 hits — check conservation
    top5 = tbl.nlargest(5, "bitscore")
    his_check = []
    for p in top5["protein"]:
        candidate = f"{g}__{p}"
        h = his_residues.get(candidate, None)
        his_check.append(h)
    his_str = ",".join([h or "?" for h in his_check])
    his_pct = sum(h == "H" for h in his_check) / max(len(his_check),1) * 100

    all_hits.append({"genus": g,
                     "n_hits_above_KOfam": n_above,
                     "n_hits_GMC_paralog_only": n_below_but_e30,
                     "best_bitscore": best,
                     "His473_conserved_top5": his_str,
                     "His473_pct_top5": his_pct,
                     "n_genomes_with_K00108": n_genomes_above})

# Add CSP1-2 reference row from prior MAG analysis
all_hits.append({"genus": "CSP1-2_EQ_MAGs",
                 "n_hits_above_KOfam": 4,    # 4/4 MAGs
                 "n_hits_GMC_paralog_only": 0,
                 "best_bitscore": None,  # not directly measured here
                 "His473_conserved_top5": "H,H,H,H",
                 "His473_pct_top5": 100.0,
                 "n_genomes_with_K00108": 4})

summary = pd.DataFrame(all_hits)
summary.to_csv(CACHE/"xcomparator_betA_summary.tsv", sep="\t", index=False)
print(summary.to_string(index=False))
print(f"\nWrote {CACHE/'xcomparator_betA_summary.tsv'}")

# Per-genome detail
per_genome_rows = []
for g in GENERA:
    tbl = parse_tbl(CACHE / f"{g}_K00108.tbl")
    if tbl.empty: continue
    tbl["genus"] = g
    tbl["genome"] = tbl["protein"].str.extract(r"^([A-Z]+_[A-Z0-9]+)", expand=False)
    for genome, sub in tbl.groupby("genome"):
        best = sub["bitscore"].max()
        n_above = (sub["bitscore"] >= KOFAM_THRESHOLD).sum()
        per_genome_rows.append({"genus": g, "genome": genome,
                                "best_bitscore": best,
                                "above_KOfam_threshold": int(n_above >= 1),
                                "n_K00108_above": n_above,
                                "n_GMC_paralogs": (~(sub["bitscore"] >= KOFAM_THRESHOLD)).sum()})
pg = pd.DataFrame(per_genome_rows)
pg.to_csv(CACHE/"xcomparator_betA_per_genome.tsv", sep="\t", index=False)
print()
print("Per-genome:")
print(pg.to_string(index=False))
print(f"\nWrote {CACHE/'xcomparator_betA_per_genome.tsv'}")
