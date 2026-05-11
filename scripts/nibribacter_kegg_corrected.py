#!/usr/bin/env python3
"""Corrected Nibribacter MAG functional analysis using SPECIFIC KO IDs
(no regex artifacts).

Categories defined by KEGG-curated KO lists:
  - CAZymes (KEGG glycoside hydrolases EC 3.2.1.x)
  - Polysaccharide lyases (EC 4.2.2.x)
  - Carbohydrate esterases (EC 3.1.1.x acetyl)
  - TonB-dependent / SusC / SusD genes
  - DNA repair (specific recA/mutS/uvr/rec gene KOs)
  - Heat shock chaperones (specific groEL/dnaK/grpE/hslU/hslV KOs)
  - Oxidative stress (specific catalase/SOD/peroxidase KOs)
  - Compatible solutes:
      * Betaine biosynthesis (K00108 betA, K00130 betB, K17755 gbsA)
      * Trehalose biosynthesis (K00697 otsA, K01087 otsB, K05343)
      * Ectoine biosynthesis (K06718 ectA, K10674 ectB, K10673 ectC, K10675)
  - Compatible solute uptake (specific opu / proU/V KOs)
  - Sigma factors (K03086 rpoD, K03088 sigB, K03089 rpoH, K03092 sigF)

Reads: cache/nibribacter_mags/per_mag_ko_assignments.tsv (already produced
by hmmsearch + threshold filter).

Outputs:
  cache/nibribacter_mags/corrected_function_summary.tsv
  cache/nibribacter_mags/corrected_summary.txt
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache" / "nibribacter_mags"

# ============================================================
# CURATED KO LISTS (KEGG-specific, not regex)
# ============================================================
KO_CATEGORIES = {
    # ---------------- CAZy ----------------
    "GH_glycoside_hydrolase": [
        "K01176", "K01180", "K01181", "K01182", "K01183", "K01184", "K01185",
        "K01186", "K01187", "K01188", "K01191", "K01192", "K01193", "K01194",
        "K01195", "K01196", "K01198", "K01199", "K01200", "K01201", "K01202",
        "K01203", "K01205", "K01206", "K01207", "K01208", "K01209", "K01210",
        "K01211", "K01212", "K01213", "K01214", "K01217", "K01218", "K01220",
        "K01222", "K01225", "K01228", "K01229", "K01230", "K01233", "K01806",
        "K05343", "K05349", "K05350", "K05989", "K05991", "K07406", "K15531",
        "K15532", "K15921", "K15922", "K22934"],
    "PL_polysaccharide_lyase": [
        "K01728", "K01729", "K01730", "K05990", "K06381", "K07082", "K07117",
        "K08305", "K08309", "K10532", "K15531", "K17467"],
    "CE_carbohydrate_esterase": [
        "K01044", "K01051", "K01054", "K01057", "K01067", "K01068", "K01072",
        "K01075", "K03412", "K07106", "K15531"],
    "TonB_SusCD": [
        "K03832", "K03833", "K03834", "K03561", "K03561", "K07028"],

    # ---------------- DNA repair ----------------
    "DNA_repair_specific": [
        "K03553", "K03554", "K03555", "K03556", "K03629", "K03651", "K03654",
        "K03655", "K03571", "K03572", "K03581", "K03546", "K03547", "K03548",
        "K03549", "K03550", "K03551", "K03629", "K03657", "K03659", "K03701",
        "K03702", "K03703", "K03723", "K03722", "K01246", "K02335", "K02340",
        "K02617", "K10773", "K10774", "K10775", "K10776"],

    # ---------------- Heat shock ----------------
    "Heat_shock_specific": [
        "K04043",  # dnaK
        "K04044",  # dnaJ
        "K03687",  # grpE
        "K04077",  # groEL
        "K04078",  # groES
        "K03667",  # hslU
        "K01419",  # hslV
        "K03686",  # clpB
        "K03695",  # clpC
        "K03696",  # clpA
        "K03544",  # clpX
        "K03545"], # clpP

    # ---------------- Oxidative stress ----------------
    "OxStress_specific": [
        "K03781",  # katE catalase
        "K03782",  # katA / katG
        "K04565",  # superoxide dismutase Mn
        "K04564",  # superoxide dismutase Fe
        "K00518",  # SOD CuZn
        "K11065",  # AhpC
        "K03386",  # ahpC peroxiredoxin
        "K00432",  # gpx
        "K03671",  # OxyR
        "K11623"], # CymG

    # ---------------- Compatible solutes - BIOSYNTHESIS ----------------
    "Betaine_biosynth_strict": [
        "K00108",  # betA, choline DH
        "K00130",  # betB, betaine aldehyde DH
        "K17755"], # gbsA

    "Trehalose_biosynth_strict": [
        "K00697",  # otsA, trehalose-6P synthase
        "K01087",  # otsB, trehalose-6P phosphatase
        "K13057",  # treT
        "K00978"], # treP

    "Ectoine_biosynth_strict": [
        "K06718",  # ectA
        "K10674",  # ectB
        "K10673",  # ectC
        "K10675"], # doeA / ectD (hydroxyectoine)

    # ---------------- Compatible solutes - UPTAKE ----------------
    "Betaine_uptake_specific": [
        "K05845",  # opuA-like
        "K05846",  # opuC
        "K05847",  # opuD
        "K05848",  # opuB
        "K02000",  # proV
        "K02001",  # proW
        "K02002"], # proX

    # ---------------- Na+ / K+ pumps ----------------
    "Na_K_pumps": [
        "K03313",  # NhaA
        "K03314",  # NhaB
        "K07301",  # MrpA
        "K05384",  # Na+-NQR
        "K07302",  # MrpB
        "K07303"], # MrpC

    # ---------------- N cycle ----------------
    "N_fixation_strict": [
        "K02588",  # nifH
        "K02586",  # nifD
        "K02591",  # nifK
        "K22896",  # vnfH
        "K22897"], # anfG

    "Denitrification_strict": [
        "K00368",  # nirK
        "K15864",  # nirS
        "K04561",  # norB
        "K00376"], # nosZ

    # ---------------- S cycle ----------------
    "Sulfate_reduction_strict": [
        "K11180",  # dsrA
        "K11181",  # dsrB
        "K00958",  # sat
        "K00394",  # aprA
        "K00395"], # aprB

    # ---------------- Sigma factors ----------------
    "Sigma_factors": [
        "K03086",  # rpoD
        "K03087",  # rpoB
        "K03088",  # sigB / sigE
        "K03089",  # rpoH
        "K03090",  # sigK
        "K03092",  # sigF
        "K03100"], # sigW

    # ---------------- Photosynthesis (NEGATIVE control) ----------------
    "Photosynthesis_negctrl": [
        "K02703",  # psbA
        "K02689",  # psaA
        "K02690",  # psaB
        "K02639",  # bchB
        "K02638"], # bchM

    # ---------------- Sporulation INITIATION (NEGATIVE control for
    # Bacteroidota - Bacilli have, Bacteroidota don't) ----------------
    "Sporulation_initiation_negctrl": [
        "K06375",  # spo0A
        "K06376",  # spo0F
        "K06378",  # spo0E
        "K06398",  # spoIIE
        "K07343",  # spoIVA
        "K06381"], # spoIID
}


def main():
    print("Loading per-MAG KO assignments ...", flush=True)
    asg = pd.read_csv(CACHE / "per_mag_ko_assignments.tsv", sep="\t")
    print(f"  total assignments: {len(asg)}", flush=True)

    # Per category, per MAG: count unique KOs and total ORFs
    rec = []
    for mag in asg["mag"].unique():
        sub = asg[asg["mag"] == mag]
        for cat, ko_list in KO_CATEGORIES.items():
            hits = sub[sub["ko"].isin(ko_list)]
            rec.append({"mag": mag, "category": cat,
                          "n_kos_in_list": len(ko_list),
                          "n_unique_kos_found": hits["ko"].nunique(),
                          "n_orfs": len(hits)})
    summary = pd.DataFrame(rec)
    summary.to_csv(CACHE / "corrected_function_summary.tsv",
                    sep="\t", index=False)

    # Pretty pivot
    print("\n=== Corrected per-MAG function counts (n_orfs) ===")
    pivot_orfs = summary.pivot_table(index="category", columns="mag",
                                          values="n_orfs", fill_value=0)
    pivot_orfs["sum"] = pivot_orfs.sum(axis=1)
    pivot_orfs = pivot_orfs.sort_values("sum", ascending=False)
    print(pivot_orfs.to_string())

    print("\n=== Unique KOs detected per category per MAG ===")
    pivot_uk = summary.pivot_table(index="category", columns="mag",
                                        values="n_unique_kos_found",
                                        fill_value=0)
    pivot_uk["sum"] = pivot_uk.sum(axis=1)
    pivot_uk["of_target"] = (
        pivot_uk["sum"] / summary.drop_duplicates("category")
            .set_index("category")["n_kos_in_list"]).round(2)
    pivot_uk = pivot_uk.sort_values("sum", ascending=False)
    print(pivot_uk.to_string())

    # List actual KOs detected for key categories
    print("\n=== KOs actually detected per category (listed) ===")
    for cat in ["GH_glycoside_hydrolase", "PL_polysaccharide_lyase",
                  "CE_carbohydrate_esterase", "TonB_SusCD",
                  "Betaine_biosynth_strict", "Trehalose_biosynth_strict",
                  "Ectoine_biosynth_strict", "Betaine_uptake_specific",
                  "Na_K_pumps", "Sporulation_initiation_negctrl",
                  "N_fixation_strict", "Photosynthesis_negctrl"]:
        ko_list = KO_CATEGORIES[cat]
        found = asg[asg["ko"].isin(ko_list)]
        if len(found) == 0:
            print(f"\n  {cat}: NONE detected", flush=True)
        else:
            print(f"\n  {cat}:")
            cnt = (found.groupby("ko")
                    .agg(n=("orf", "count"),
                          definition=("definition", "first"))
                    .reset_index())
            print(cnt.to_string(index=False), flush=True)

    # Save summary
    with open(CACHE / "corrected_summary.txt", "w") as fh:
        fh.write("Nibribacter MAG functional characterization (CORRECTED)\n")
        fh.write("=" * 70 + "\n\n")
        fh.write("Using KEGG-curated specific KO IDs (no regex artifacts).\n\n")
        fh.write("=== n_orfs per category per MAG ===\n")
        fh.write(pivot_orfs.to_string())
        fh.write("\n\n=== Unique KOs found per category ===\n")
        fh.write(pivot_uk.to_string())


if __name__ == "__main__":
    main()
