#!/usr/bin/env python3
"""Extraction-blank contaminant screen extended to the Trip 4 blanks.

Same method as scripts/controls/run_assay_aware_control_audit.py (the
published screen): for each ASV, presence in extraction blanks versus presence
in the biological profiles that were extracted in the same batches, one-sided
Fisher exact test (alternative="greater"), candidate when
score < --primary-threshold (0.10), present in >= --minimum-blanks (2) blanks,
and blank prevalence > biological prevalence.  Candidates are removed only from
the mapped profiles, in a sensitivity copy of the table.

Extension: the batch map (inputs/extraction_batch_map_extended.tsv) carries
both the 17 Trip 5 blanks (EB1-EB17, columns of the canonical table) and the
6 Trip 4 blanks (SEB, EB1-EB5 of run novaseq_14_07_25).  The Trip 4 blank
libraries are absent from the canonical table, so their ASV profiles must be
supplied with --trip4-blank-table (TSV '#OTU ID' x library columns, or a BIOM
file, with the same md5 ASV identifiers as the canonical table).  Without it
the Trip 4 screen is reported as not run and the Trip 5 screen reproduces the
published result (validation mode).

--mode per-trip (default): two screens, Trip 5 blanks vs Trip 5 mapped
profiles and Trip 4 blanks vs Trip 4 mapped profiles, each candidate set
applied to its own trip.  --mode pooled: all blanks vs all mapped profiles,
one candidate set applied to both trips.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NEGATIVE_COLUMNS = [f"EB{i}" for i in range(1, 19)] + [
    "Negative1", "Negative2", "Negative4", "Negative5", "Negative6", "Negative7",
]


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            d.update(block)
    return d.hexdigest()


def write_rows(path, rows, fieldnames=None):
    fieldnames = fieldnames or (list(rows[0]) if rows else ["empty"])
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def load_taxonomy(path):
    out = {}
    with path.open(encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        idk = r.fieldnames[0]
        tk = "Taxon" if "Taxon" in r.fieldnames else r.fieldnames[1]
        for row in r:
            out[row[idk]] = row[tk]
    return out


def load_batch_map(path):
    """Return {trip: {blank_label: {"column": blank_col, "samples": [canonical cols]}}}."""
    batches = defaultdict(dict)
    sample_to_blank = {}
    with path.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            trip = row["trip"]
            rec = batches[trip].setdefault(
                row["extraction_blank"], {"column": row["blank_table_column"], "samples": [], "missing": []}
            )
            if row["in_canonical_table"] == "True":
                rec["samples"].append(row["canonical_column"])
                sample_to_blank[row["canonical_column"]] = row["extraction_blank"]
            else:
                rec["missing"].append(row["sample_id"])
    return batches, sample_to_blank


def load_trip4_blank_table(path, wanted_columns):
    """Return (column_names_found, {feature_id: np.array of counts over found columns})."""
    if path.suffix == ".biom":
        from biom import load_table
        t = load_table(str(path))
        ids = list(t.ids(axis="sample"))
        found = [c for c in wanted_columns if c in ids]
        alt = {c: c.replace("TRIP4_", "") for c in wanted_columns}
        found = found or [c for c in wanted_columns if alt[c] in ids]
        if not found:
            raise ValueError(f"none of {wanted_columns} found in {path}; samples: {ids}")
        names = [c if c in ids else alt[c] for c in found]
        data = {}
        obs = list(t.ids(axis="observation"))
        mat = np.column_stack([t.data(n, axis="sample", dense=True) for n in names])
        for i, f in enumerate(obs):
            if mat[i].any():
                data[str(f)] = mat[i].astype(float)
        return found, data
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n")
        while header.startswith("#") and not header.startswith("#OTU"):
            header = fh.readline().rstrip("\n")
        cols = header.split("\t")[1:]
        alt = {c: c.replace("TRIP4_", "") for c in wanted_columns}
        pos = {}
        for c in wanted_columns:
            if c in cols:
                pos[c] = cols.index(c)
            elif alt[c] in cols:
                pos[c] = cols.index(alt[c])
        if not pos:
            raise ValueError(f"none of {wanted_columns} found in {path}; columns: {cols}")
        found = list(pos)
        idx = [pos[c] for c in found]
        data = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            v = np.array([float(f[i + 1]) for i in idx])
            if v.any():
                data[f[0]] = v
        return found, data


def fisher_scores(cp, n_blanks, bp, n_bio):
    scores = np.ones(len(cp), dtype=float)
    for i in np.flatnonzero(cp):
        scores[i] = fisher_exact(
            [[int(cp[i]), n_blanks - int(cp[i])], [int(bp[i]), n_bio - int(bp[i])]],
            alternative="greater",
        ).pvalue
    return scores


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canonical-table", type=Path,
                    default=ROOT / "data/processed/taxonomy/taxon-tables/feature-table-trips1-5.tsv")
    ap.add_argument("--taxonomy", type=Path,
                    default=ROOT / "data/processed/taxonomy/taxon-tables/taxonomy-trips1-5.tsv")
    ap.add_argument("--batch-map", type=Path, default=HERE / "inputs/extraction_batch_map_extended.tsv")
    ap.add_argument("--profile-metadata", type=Path, default=ROOT / "analysis/v2/review/cache/alpha.tsv")
    ap.add_argument("--trip4-blank-table", type=Path, default=None)
    ap.add_argument("--require-trip4", action="store_true")
    ap.add_argument("--mode", choices=("per-trip", "pooled"), default="per-trip")
    ap.add_argument("--primary-threshold", type=float, default=0.10)
    ap.add_argument("--minimum-blanks", type=int, default=2)
    ap.add_argument("--output-dir", type=Path, default=HERE / "outputs")
    args = ap.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    t0 = datetime.now(timezone.utc)
    log = lambda *a: print(f"[{datetime.now(timezone.utc).isoformat(timespec='seconds')}]", *a, flush=True)

    batches, sample_to_blank = load_batch_map(args.batch_map)
    trip5_blanks = sorted(batches["Trip5"], key=lambda b: int(b[2:]))
    trip4_blanks = list(batches["Trip4"])
    trip5_cols = [batches["Trip5"][b]["column"] for b in trip5_blanks]
    trip4_cols = [batches["Trip4"][b]["column"] for b in trip4_blanks]
    trip5_mapped = sorted({s for b in trip5_blanks for s in batches["Trip5"][b]["samples"]})
    trip4_mapped = sorted({s for b in trip4_blanks for s in batches["Trip4"][b]["samples"]})
    log(f"batch map: Trip5 {len(trip5_blanks)} blanks / {len(trip5_mapped)} mapped profiles; "
        f"Trip4 {len(trip4_blanks)} blanks / {len(trip4_mapped)} mapped profiles")

    trip4_found, trip4_data = [], {}
    if args.trip4_blank_table is not None:
        trip4_found, trip4_data = load_trip4_blank_table(args.trip4_blank_table, trip4_cols)
        log(f"Trip 4 blank table: {len(trip4_found)} blank columns {trip4_found}, {len(trip4_data)} non-empty ASVs")
    elif args.require_trip4:
        raise SystemExit("--require-trip4 given but no --trip4-blank-table: the Trip 4 blank "
                         "libraries are not in the canonical table (see README.md)")
    else:
        log("no --trip4-blank-table: Trip 4 screen NOT run (validation of the Trip 5 screen only)")
    trip4_ran = bool(trip4_found)

    # Pass 1: stream the canonical table
    with args.canonical_table.open(encoding="utf-8") as fh:
        provenance = fh.readline().rstrip("\n")
        columns = fh.readline().rstrip("\n").split("\t")[1:]
        index = {c: i for i, c in enumerate(columns)}
        for c in trip5_cols:
            if c not in index:
                raise ValueError(f"Trip 5 blank column {c} absent from canonical table")
        i5b = np.array([index[c] for c in trip5_cols])
        i5m = np.array([index[c] for c in trip5_mapped])
        i4m = np.array([index[c] for c in trip4_mapped])
        feats, cp5, bp5, bp4, cr5, br5, br4 = [], [], [], [], [], [], []
        col_total = np.zeros(len(columns)); col_feats = np.zeros(len(columns), dtype=np.int64)
        n = 0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            v = np.fromiter((float(x) for x in f[1:]), dtype=float, count=len(columns))
            p = v > 0
            feats.append(f[0])
            cp5.append(int(p[i5b].sum())); bp5.append(int(p[i5m].sum())); bp4.append(int(p[i4m].sum()))
            cr5.append(float(v[i5b].sum())); br5.append(float(v[i5m].sum())); br4.append(float(v[i4m].sum()))
            col_total += v; col_feats += p
            n += 1
            if n % 50000 == 0:
                log(f"pass 1: {n} features")
    log(f"pass 1 done: {len(feats)} features, {len(columns)} profiles")
    cp5, bp5, bp4 = map(np.asarray, (cp5, bp5, bp4))
    cr5, br5, br4 = map(np.asarray, (cr5, br5, br4))
    fpos = {f: i for i, f in enumerate(feats)}
    cp4 = np.zeros(len(feats), dtype=np.int64); cr4 = np.zeros(len(feats))
    trip4_blank_only = 0
    for f, v in trip4_data.items():
        i = fpos.get(f)
        if i is None:
            trip4_blank_only += 1
            continue
        cp4[i] = int((v > 0).sum()); cr4[i] = float(v.sum())
    if trip4_ran:
        log(f"Trip 4 blank ASVs: {len(trip4_data)} non-empty, {trip4_blank_only} absent from canonical table")

    # Screens
    screens = {}
    if args.mode == "per-trip":
        screens["Trip5"] = dict(cp=cp5, nb=len(trip5_cols), bp=bp5, nm=len(trip5_mapped), cr=cr5, br=br5,
                                blanks=trip5_blanks, mapped=trip5_mapped)
        if trip4_ran:
            screens["Trip4"] = dict(cp=cp4, nb=len(trip4_found), bp=bp4, nm=len(trip4_mapped), cr=cr4, br=br4,
                                    blanks=[b for b, c in zip(trip4_blanks, trip4_cols) if c in trip4_found],
                                    mapped=trip4_mapped)
    else:
        if not trip4_ran:
            raise SystemExit("pooled mode needs the Trip 4 blank table")
        screens["pooled"] = dict(cp=cp5 + cp4, nb=len(trip5_cols) + len(trip4_found), bp=bp5 + bp4,
                                 nm=len(trip5_mapped) + len(trip4_mapped), cr=cr5 + cr4, br=br5 + br4,
                                 blanks=trip5_blanks + trip4_blanks, mapped=trip5_mapped + trip4_mapped)
    taxonomy = load_taxonomy(args.taxonomy)
    call_rows, sens_rows, called = [], [], {}
    for name, s in screens.items():
        log(f"screen {name}: Fisher on {int((s['cp'] > 0).sum())} blank-present ASVs")
        scores = fisher_scores(s["cp"], s["nb"], s["bp"], s["nm"])
        cprev = s["cp"] / s["nb"]; bprev = s["bp"] / s["nm"]
        primary = None
        for thr in (0.01, 0.05, 0.10):
            for mb in (1, 2, 3):
                c = (scores < thr) & (s["cp"] >= mb) & (cprev > bprev)
                sens_rows.append({"screen": name, "prevalence_score_threshold": thr,
                                  "minimum_blank_prevalence_count": mb, "called_features": int(c.sum()),
                                  "control_reads_in_called_features": int(s["cr"][c].sum()),
                                  "biological_reads_in_called_features": int(s["br"][c].sum())})
                if thr == args.primary_threshold and mb == args.minimum_blanks:
                    primary = c
        if primary is None:
            primary = (scores < args.primary_threshold) & (s["cp"] >= args.minimum_blanks) & (cprev > bprev)
        called[name] = {feats[i] for i in np.flatnonzero(primary)}
        s["called"] = primary
        for i in np.flatnonzero(primary):
            call_rows.append({"feature_id": feats[i], "screen": name, "taxonomy": taxonomy.get(feats[i], ""),
                              "extraction_blanks_present": int(s["cp"][i]), "extraction_blanks_total": s["nb"],
                              "mapped_biological_profiles_present": int(s["bp"][i]),
                              "mapped_biological_profiles_total": s["nm"],
                              "blank_prevalence": f"{cprev[i]:.8f}", "biological_prevalence": f"{bprev[i]:.8f}",
                              "prevalence_score": f"{scores[i]:.10g}", "blank_reads": int(s["cr"][i]),
                              "mapped_biological_reads": int(s["br"][i]), "call": "candidate_contaminant"})
        log(f"screen {name}: {int(primary.sum())} primary candidates")
    call_rows.sort(key=lambda r: (r["screen"], float(r["prevalence_score"]), -int(r["blank_reads"])))
    write_rows(out / "primary_contaminant_calls.tsv", call_rows, [
        "feature_id", "screen", "taxonomy", "extraction_blanks_present", "extraction_blanks_total",
        "mapped_biological_profiles_present", "mapped_biological_profiles_total", "blank_prevalence",
        "biological_prevalence", "prevalence_score", "blank_reads", "mapped_biological_reads", "call"])
    write_rows(out / "filter_sensitivity.tsv", sens_rows)

    # Which candidate set applies to which column
    set_for_col = {}
    if args.mode == "per-trip":
        for c in trip5_mapped + trip5_cols + [c for c in NEGATIVE_COLUMNS if c in index]:
            set_for_col[c] = called["Trip5"]
        if trip4_ran:
            for c in trip4_mapped:
                set_for_col[c] = called["Trip4"]
    else:
        for c in trip5_mapped + trip4_mapped + trip5_cols + [c for c in NEGATIVE_COLUMNS if c in index]:
            set_for_col[c] = called["pooled"]
    all_called = set().union(*called.values())
    filtered_cols = trip5_mapped + (trip4_mapped if trip4_ran else [])
    fpos_cols = [index[c] for c in filtered_cols]

    # Pass 2: removal per column and filtered tables
    rem_reads = np.zeros(len(columns)); rem_feats = np.zeros(len(columns), dtype=np.int64)
    trip5_path = out / "trip5_mapped_feature_table_control_filtered.tsv.gz"
    ext_path = out / "mapped_feature_table_control_filtered.tsv.gz"
    retained5 = retained_ext = 0
    t5pos = [index[c] for c in trip5_mapped]
    with args.canonical_table.open(encoding="utf-8") as fh, \
            gzip.GzipFile(filename="", fileobj=trip5_path.open("wb"), mode="wb", mtime=0) as g5, \
            io.TextIOWrapper(g5, encoding="utf-8", newline="") as h5, \
            gzip.GzipFile(filename="", fileobj=ext_path.open("wb"), mode="wb", mtime=0) as ge, \
            io.TextIOWrapper(ge, encoding="utf-8", newline="") as he:
        fh.readline(); fh.readline()
        h5.write(f"# Control-sensitivity table; source={args.canonical_table}; threshold={args.primary_threshold}; "
                 f"minimum_blanks={args.minimum_blanks}; training_controls=EB1-EB17\n")
        h5.write("#OTU ID\t" + "\t".join(trip5_mapped) + "\n")
        he.write(f"# Control-sensitivity table (rerun {HERE.name}); mode={args.mode}; "
                 f"threshold={args.primary_threshold}; minimum_blanks={args.minimum_blanks}; "
                 f"trip4_screen={'run' if trip4_ran else 'not_run'}\n")
        he.write("#OTU ID\t" + "\t".join(filtered_cols) + "\n")
        n = 0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            fid = f[0]
            n += 1
            if n % 50000 == 0:
                log(f"pass 2: {n} features")
            if fid in all_called:
                v = np.fromiter((float(x) for x in f[1:]), dtype=float, count=len(columns))
                for c, s in set_for_col.items():
                    if fid in s:
                        i = index[c]
                        rem_reads[i] += v[i]; rem_feats[i] += v[i] > 0
            t5set = called.get("Trip5", called.get("pooled"))
            if fid not in t5set:
                sel = [f[p + 1] for p in t5pos]
                if any(float(x) > 0 for x in sel):
                    h5.write(fid + "\t" + "\t".join(sel) + "\n"); retained5 += 1
            sel = []
            for c, p in zip(filtered_cols, fpos_cols):
                sel.append("0" if fid in set_for_col[c] else f[p + 1])
            if any(float(x) > 0 for x in sel):
                he.write(fid + "\t" + "\t".join(sel) + "\n"); retained_ext += 1
    log(f"pass 2 done: trip5 filtered table {retained5} features; extended filtered table {retained_ext} features")

    # Per-profile removal
    with args.profile_metadata.open(newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        meta = {row[r.fieldnames[0]]: row for row in r}
    mapped_set = set(trip5_mapped) | (set(trip4_mapped) if trip4_ran else set())
    prof_rows = []
    for i, c in enumerate(columns):
        in_neg = c in NEGATIVE_COLUMNS
        if c not in mapped_set and not in_neg and c not in trip4_mapped:
            continue
        trip = "Trip5" if c in trip5_mapped or c in trip5_cols else "Trip4" if c in trip4_mapped else ""
        role = ("training_extraction_blank" if c in trip5_cols else
                "characterization_only_extraction_blank" if in_neg else
                "compatible_biological_profile" if c in mapped_set else
                "mapped_biological_profile_screen_not_run")
        tot = col_total[i]
        prof_rows.append({"profile_id": c, "trip": trip or meta.get(c, {}).get("Trip", ""), "role": role,
                          "extraction_blank": sample_to_blank.get(c, ""), "total_reads": int(tot),
                          "total_asvs": int(col_feats[i]), "candidate_contaminant_reads": int(rem_reads[i]),
                          "candidate_contaminant_asvs": int(rem_feats[i]),
                          "candidate_contaminant_read_fraction": f"{rem_reads[i] / tot:.8f}" if tot else "",
                          "application_scope": ("removed_in_sensitivity_table" if c in mapped_set
                                                else "characterized_not_filtered")})
    if trip4_ran:
        for b, c in zip(trip4_blanks, trip4_cols):
            if c not in trip4_found:
                continue
            j = trip4_found.index(c)
            tot = sum(v[j] for v in trip4_data.values()); nf = sum(1 for v in trip4_data.values() if v[j] > 0)
            s = set_for_col.get(trip4_mapped[0], set()) if trip4_mapped else set()
            rr = sum(v[j] for f, v in trip4_data.items() if f in s)
            rf = sum(1 for f, v in trip4_data.items() if f in s and v[j] > 0)
            prof_rows.append({"profile_id": c, "trip": "Trip4", "role": "training_extraction_blank",
                              "extraction_blank": b, "total_reads": int(tot), "total_asvs": int(nf),
                              "candidate_contaminant_reads": int(rr), "candidate_contaminant_asvs": int(rf),
                              "candidate_contaminant_read_fraction": f"{rr / tot:.8f}" if tot else "",
                              "application_scope": "characterized_not_filtered"})
    write_rows(out / "removal_fraction_by_profile.tsv", prof_rows)
    bio = [r for r in prof_rows if r["role"] == "compatible_biological_profile"]
    for gtype, field, fn in (("campaign", "Trip", "removal_fraction_by_campaign.tsv"),
                             ("compartment", "Type", "removal_fraction_by_compartment.tsv")):
        grouped = defaultdict(list)
        for r in bio:
            grouped[meta.get(r["profile_id"], {}).get(field, "unknown")].append(r)
        rows = []
        for g, mem in sorted(grouped.items(), key=lambda kv: str(kv[0])):
            tr = sum(int(r["total_reads"]) for r in mem); rr = sum(int(r["candidate_contaminant_reads"]) for r in mem)
            fr = [float(r["candidate_contaminant_read_fraction"]) for r in mem]
            rows.append({"group_type": gtype, "group_value": g, "biological_profile_count": len(mem),
                         "total_reads": tr, "candidate_contaminant_reads": rr,
                         "pooled_candidate_contaminant_read_fraction": f"{rr / tr:.8f}" if tr else "",
                         "median_candidate_contaminant_read_fraction": f"{np.median(fr):.8f}",
                         "maximum_candidate_contaminant_read_fraction": f"{max(fr):.8f}"})
        write_rows(out / fn, rows)
    by_name = {r["profile_id"]: r for r in prof_rows}
    batch_rows = []
    for trip in ("Trip5", "Trip4"):
        for b, rec in batches[trip].items():
            pres = [s for s in rec["samples"] if s in by_name]
            fr = [float(by_name[s]["candidate_contaminant_read_fraction"]) for s in pres
                  if by_name[s]["candidate_contaminant_read_fraction"] != ""]
            blank = by_name.get(rec["column"], {})
            batch_rows.append({"trip": trip, "extraction_blank": b, "blank_column": rec["column"],
                               "blank_in_analysis": str(bool(blank)),
                               "workbook_sample_count": len(rec["samples"]) + len(rec["missing"]),
                               "canonical_sample_count": len(pres), "canonical_missing_sample_count": len(rec["missing"]),
                               "blank_reads": blank.get("total_reads", ""), "blank_asvs": blank.get("total_asvs", ""),
                               "median_biological_removal_fraction": f"{np.median(fr):.8f}" if fr else "",
                               "maximum_biological_removal_fraction": f"{max(fr):.8f}" if fr else ""})
    write_rows(out / "extraction_batch_summary.tsv", batch_rows)

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "runtime_seconds": (datetime.now(timezone.utc) - t0).total_seconds(),
        "mode": args.mode, "primary_method": "one-sided Fisher prevalence enrichment",
        "primary_prevalence_score_threshold": args.primary_threshold,
        "primary_minimum_blank_prevalence_count": args.minimum_blanks,
        "canonical_table_provenance_line": provenance, "canonical_profiles": len(columns),
        "canonical_features": len(feats),
        "trip5_training_blanks": trip5_blanks, "trip5_mapped_profiles": len(trip5_mapped),
        "trip4_training_blanks": [b for b, c in zip(trip4_blanks, trip4_cols) if c in trip4_found],
        "trip4_blanks_in_map": trip4_blanks, "trip4_mapped_profiles": len(trip4_mapped),
        "trip4_screen": "run" if trip4_ran else "not_run: Trip 4 blank ASV profiles unavailable (libraries never denoised)",
        "trip4_blank_table": str(args.trip4_blank_table) if args.trip4_blank_table else None,
        "trip4_blank_only_features_not_in_canonical_table": trip4_blank_only if trip4_ran else None,
        "candidates": {k: len(v) for k, v in called.items()},
        "candidates_union": len(all_called),
        "trip5_filtered_table_retained_features": retained5,
        "extended_filtered_table_profiles": len(filtered_cols),
        "extended_filtered_table_retained_features": retained_ext,
        "input_sha256": {str(p): sha256(p) for p in
                         [args.canonical_table, args.taxonomy, args.batch_map, args.profile_metadata]
                         + ([args.trip4_blank_table] if args.trip4_blank_table else [])},
        "output_sha256": {p.name: sha256(p) for p in (trip5_path, ext_path)},
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    log("done"); print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
