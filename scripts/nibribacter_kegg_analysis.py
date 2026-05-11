#!/usr/bin/env python3
"""Functional characterization of Nibribacter MAGs from hmmsearch + ko_list
thresholds.

Identifies:
  - DOM-cycling capacity: CAZymes (KEGG glycoside hydrolases, polysacc lyases),
    TonB-dependent transporters (BUT TonB are PFAM not KEGG; capture via gene
    names), polysaccharide-utilization-loci marker genes (susC, susD)
  - N cycle: nifH, nirK, nirS, narG, etc.
  - Stress: trehalose biosynthesis (otsAB), betaine (proU, betA),
    DNA repair (recA, mutS), heat shock (groEL, dnaK)
  - Membrane / osmotic: porin, water-channel
  - Photosynthesis (negative control - should be absent)
  - Sporulation (negative control - Bacteroidota don't sporulate)

For each MAG: count KOs and KO-set membership.

Outputs:
  cache/nibribacter_mags/per_mag_ko_assignments.tsv
  cache/nibribacter_mags/per_mag_function_summary.tsv
  cache/nibribacter_mags/cazy_kos_per_mag.tsv
  cache/nibribacter_mags/summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import re
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache" / "nibribacter_mags"

MAGS = ["15PRr3_SemiBin_102", "19PRr2_SemiBin_101", "22PRr3_SemiBin_75"]


def parse_hmm_tbl(path: Path) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            parts = line.split()
            if len(parts) < 18: continue
            target_name = parts[0]      # ORF
            target_acc = parts[1]
            query_name = parts[2]       # Kxxxxx
            query_acc = parts[3]
            evalue = float(parts[4])
            score = float(parts[5])
            rows.append({"orf": target_name, "ko": query_name,
                          "evalue": evalue, "score": score})
    return pd.DataFrame(rows)


def main():
    print("Loading ko_list ...", flush=True)
    ko = pd.read_csv(CACHE / "ko_list.txt", sep="\t")
    ko["threshold"] = pd.to_numeric(ko["threshold"], errors="coerce")
    ko_thresh = ko.set_index("knum")["threshold"].to_dict()
    ko_def = ko.set_index("knum")["definition"].to_dict()
    print(f"  KO entries: {len(ko)}, with threshold: "
          f"{ko['threshold'].notna().sum()}", flush=True)

    # Per-MAG KO assignment with KOFAM-style threshold filter
    all_assignments = []
    per_mag_ko = {}
    for mag in MAGS:
        tbl = parse_hmm_tbl(CACHE / f"{mag}.hmm.tbl")
        print(f"  {mag}: {len(tbl)} raw HMM hits", flush=True)
        # Apply threshold filter: keep hits with score >= ko_threshold
        tbl["threshold"] = tbl["ko"].map(ko_thresh)
        tbl["pass"] = tbl["score"] >= tbl["threshold"]
        tbl["mag"] = mag
        # If no threshold available, use evalue<1e-30 as fallback
        no_thresh = tbl["threshold"].isna()
        tbl.loc[no_thresh, "pass"] = tbl.loc[no_thresh, "evalue"] < 1e-30
        kept = tbl[tbl["pass"]]
        # Per ORF: best hit only
        kept = (kept.sort_values("score", ascending=False)
                  .drop_duplicates("orf"))
        print(f"    -> {len(kept)} KO-assigned ORFs (after threshold)",
              flush=True)
        per_mag_ko[mag] = set(kept["ko"])
        kept["definition"] = kept["ko"].map(ko_def)
        all_assignments.append(kept[["mag", "orf", "ko", "score",
                                          "definition"]])
    asg = pd.concat(all_assignments, ignore_index=True)
    asg.to_csv(CACHE / "per_mag_ko_assignments.tsv", sep="\t", index=False)
    print(f"\nTotal KO assignments: {len(asg)}", flush=True)

    # Per-MAG: # unique KOs
    print("\n=== Per-MAG KO breadth ===")
    for mag in MAGS:
        sub = asg[asg["mag"] == mag]
        n_orfs = sub["orf"].nunique()
        n_kos = sub["ko"].nunique()
        print(f"  {mag}:  {n_orfs} ORFs assigned KO, {n_kos} unique KOs",
              flush=True)

    # ==================================================================
    # Function categories
    # ==================================================================
    # CAZymes via KEGG: glycoside hydrolases EC 3.2.1.x, polysaccharide
    # lyases EC 4.2.2.x, carbohydrate esterases. Bacteroidota DOM cyclers
    # are characterized by HUNDREDS of GH and PL genes.
    # We capture by definition-string keyword.
    KEYWORDS = {
        "glycoside_hydrolase": [r"glycoside hydrolase", r"glycosidase",
                                  r"glucosidase", r"galactosidase",
                                  r"mannosidase", r"xylosidase",
                                  r"arabinosidase", r"xylanase", r"mannanase",
                                  r"cellulase", r"chitinase",
                                  r"beta-glucanase", r"alpha-amylase",
                                  r"\[EC:3\.2\.1\."],
        "polysaccharide_lyase": [r"polysaccharide lyase", r"\[EC:4\.2\.2\."],
        "carbohydrate_esterase": [r"acetyl xylan esterase",
                                     r"acetylesterase",
                                     r"\[EC:3\.1\.1\.6", r"\[EC:3\.1\.1\.72"],
        "TonB_dependent_porin": [r"TonB", r"tonB", r"susC", r"susD",
                                    r"sus", r"glycan transporter",
                                    r"polysaccharide-specific"],
        "starch_degradation": [r"alpha-amylase", r"alpha-glucosidase",
                                  r"pullulanase", r"\[EC:3\.2\.1\.1"],
        "cellulose_degradation": [r"cellulase", r"cellobiohydrolase",
                                    r"endoglucanase", r"\[EC:3\.2\.1\.4"],
        "chitin_degradation": [r"chitinase", r"\[EC:3\.2\.1\.14"],
        "xylan_degradation": [r"xylanase", r"xylosidase",
                                 r"\[EC:3\.2\.1\.8", r"\[EC:3\.2\.1\.37"],
        # Stress
        "trehalose_biosynth": [r"trehalose-6-phosphate", r"otsA", r"otsB",
                                  r"trehalose synthase",
                                  r"\[EC:2\.4\.1\.15"],
        "betaine_biosynth": [r"choline dehydrogenase", r"betA",
                                r"betaine aldehyde", r"betB",
                                r"\[EC:1\.1\.99\.1"],
        "betaine_uptake": [r"glycine betaine.*transport",
                              r"proU", r"proV", r"opuA", r"opuB", r"opuC"],
        "DNA_repair": [r"DNA repair", r"recA", r"mutS", r"mutL",
                          r"uvrA", r"uvrB", r"uvrC",
                          r"DNA polymerase IV", r"DNA mismatch"],
        "heat_shock": [r"heat shock", r"groEL", r"groES", r"dnaK", r"dnaJ"],
        "OS_response": [r"alkyl hydroperoxide", r"katE", r"sodA",
                          r"superoxide", r"oxidative stress"],
        # N cycle
        "N_fixation": [r"nifH", r"nifD", r"nifK", r"vnf", r"anf"],
        "denitrification": [r"nirK", r"nirS", r"norB", r"nosZ"],
        "nitrification": [r"amoA", r"hao", r"nxr"],
        "nitrate_reduction": [r"narG", r"narH", r"napA", r"napB"],
        # S cycle
        "S_oxidation": [r"soxA", r"soxB", r"soxC", r"soxX",
                          r"sulfide.*oxido"],
        "S_reduction": [r"dsrA", r"dsrB", r"aprA", r"aprB"],
        # Photosynthesis (NEGATIVE CONTROL)
        "photosynthesis": [r"chlorophyll", r"photosystem", r"psbA", r"psaA",
                              r"reaction center"],
        # Sporulation (Bacteroidota negative control)
        "sporulation": [r"sporulation", r"spoIIIE", r"spoIIE", r"spoIIA",
                          r"germination"],
        # Membrane
        "porin": [r"\bporin\b", r"OmpA", r"outer membrane.*protein"],
        # Carbohydrate uptake
        "ABC_sugar_transport": [r"sugar.*transport", r"ABC.*sugar",
                                   r"saccharide.*transport"],
    }

    rec = []
    for mag in MAGS:
        sub = asg[asg["mag"] == mag]
        for cat, keywords in KEYWORDS.items():
            hits = []
            for _, row in sub.iterrows():
                d = str(row["definition"])
                for kw in keywords:
                    if re.search(kw, d, re.IGNORECASE):
                        hits.append((row["ko"], d))
                        break
            n_unique_ko = len(set(h[0] for h in hits))
            n_orfs = len(hits)
            rec.append({"mag": mag, "category": cat,
                          "n_unique_kos": n_unique_ko,
                          "n_orfs": n_orfs,
                          "example_kos": ",".join(sorted(set(h[0] for h
                                                              in hits))[:5])})
    summary = pd.DataFrame(rec)
    summary.to_csv(CACHE / "per_mag_function_summary.tsv", sep="\t",
                    index=False)

    # Pivot: category x mag
    pivot = summary.pivot_table(index="category", columns="mag",
                                    values="n_orfs", fill_value=0)
    pivot["sum"] = pivot.sum(axis=1)
    pivot = pivot.sort_values("sum", ascending=False)
    print("\n=== Function category counts (n ORFs assigned per MAG) ===")
    print(pivot.to_string())

    # Save with definitions
    rec_full = []
    for mag in MAGS:
        sub = asg[asg["mag"] == mag]
        for cat, keywords in KEYWORDS.items():
            for _, row in sub.iterrows():
                d = str(row["definition"])
                for kw in keywords:
                    if re.search(kw, d, re.IGNORECASE):
                        rec_full.append({"mag": mag, "category": cat,
                                          "orf": row["orf"],
                                          "ko": row["ko"],
                                          "score": row["score"],
                                          "definition": d})
                        break
    full = pd.DataFrame(rec_full)
    full.to_csv(CACHE / "per_mag_function_full.tsv", sep="\t", index=False)

    # CAZyme highlight
    print("\n=== Top 20 CAZyme-related KOs (across MAGs) ===")
    cazy = full[full["category"].isin(["glycoside_hydrolase",
                                            "polysaccharide_lyase",
                                            "carbohydrate_esterase",
                                            "starch_degradation",
                                            "cellulose_degradation",
                                            "chitin_degradation",
                                            "xylan_degradation"])]
    print(f"  Total CAZy-related ORFs across 3 MAGs: {len(cazy)}",
          flush=True)
    print(f"  Unique CAZy KOs: {cazy['ko'].nunique()}", flush=True)
    cazy_summary = (cazy.groupby("ko")
                     .agg(definition=("definition", "first"),
                          n_orfs=("orf", "count"))
                     .sort_values("n_orfs", ascending=False)
                     .reset_index())
    cazy_summary.to_csv(CACHE / "cazy_kos_per_mag.tsv",
                         sep="\t", index=False)
    print(cazy_summary.head(20).to_string(index=False), flush=True)

    # ==================================================================
    # Compare to expected Bacteroidota DOM-cycler signature
    # ==================================================================
    print("\n=== DOM-cycler signature interpretation ===")
    print(f"  Bacteroidota DOM-cyclers typically have:")
    print(f"    - 30-100+ glycoside hydrolases (GH)")
    print(f"    - many TonB-dependent transporters (TBDT) for polysaccharide uptake")
    print(f"    - polysaccharide-utilization-loci (PUL) susC/susD pairs")
    print()
    print(f"  Our 3 Nibribacter MAGs:")
    for mag in MAGS:
        sub = full[full["mag"] == mag]
        n_GH = sub[sub["category"] == "glycoside_hydrolase"]["orf"].nunique()
        n_PL = sub[sub["category"] == "polysaccharide_lyase"]["orf"].nunique()
        n_CE = sub[sub["category"] == "carbohydrate_esterase"]["orf"].nunique()
        n_starch = sub[sub["category"] == "starch_degradation"]["orf"].nunique()
        n_TBDT = sub[sub["category"] == "TonB_dependent_porin"]["orf"].nunique()
        n_porin = sub[sub["category"] == "porin"]["orf"].nunique()
        n_treh = sub[sub["category"] == "trehalose_biosynth"]["orf"].nunique()
        n_betA = sub[sub["category"] == "betaine_biosynth"]["orf"].nunique()
        n_betU = sub[sub["category"] == "betaine_uptake"]["orf"].nunique()
        n_DNA = sub[sub["category"] == "DNA_repair"]["orf"].nunique()
        n_heat = sub[sub["category"] == "heat_shock"]["orf"].nunique()
        n_OS = sub[sub["category"] == "OS_response"]["orf"].nunique()
        n_nif = sub[sub["category"] == "N_fixation"]["orf"].nunique()
        n_phot = sub[sub["category"] == "photosynthesis"]["orf"].nunique()
        n_spo = sub[sub["category"] == "sporulation"]["orf"].nunique()
        print(f"\n  --- {mag} ---")
        print(f"    DOM-cycling:  GH={n_GH:<3}  PL={n_PL:<3}  CE={n_CE:<3}  "
              f"starch={n_starch:<3}  TBDT/SusCD={n_TBDT:<3}  porin={n_porin}")
        print(f"    Stress:       trehalose_biosyn={n_treh}  "
              f"betaine_biosyn={n_betA}  betaine_uptake={n_betU}  "
              f"DNA_repair={n_DNA}  heat_shock={n_heat}  "
              f"OS_response={n_OS}")
        print(f"    Negative ctl: N_fixation={n_nif}  photosynthesis={n_phot}  "
              f"sporulation={n_spo}")

    # Save summary
    with open(CACHE / "summary.txt", "w") as fh:
        fh.write("Nibribacter MAG functional characterization\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"3 partial Nibribacter MAGs from rhizosphere of EQ Trip 1:\n")
        for mag in MAGS:
            n_kos = asg[asg["mag"] == mag]["ko"].nunique()
            n_orfs = asg[asg["mag"] == mag]["orf"].nunique()
            fh.write(f"  {mag}: {n_orfs} ORFs assigned KO, "
                      f"{n_kos} unique KOs\n")
        fh.write("\nFunction category counts (n ORFs):\n")
        fh.write(pivot.to_string())
        fh.write("\n\nTop 20 CAZyme KOs:\n")
        fh.write(cazy_summary.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
