#!/usr/bin/env python3
"""TEST 2: Cosmopolitanism via Earth Microbiome Project.

Strategy:
  1. Download EMP deblur "cosmopolitan ASV" catalogs (90/100/150 bp,
     present in >=25 EMP samples globally).
     URL: https://ftp.microbio.me/emp/release1/otu_info/deblur/
  2. Trim each EQ ASV (V3-V4 ~ 410-450 bp) to the V4 region (~ first 90/100/150 bp
     after the V4 forward primer 515F = GTGYCAGCMGCCGCGGTAA).
  3. Exact-match (vsearch -id 1.0) EQ V4-trimmed reps vs EMP catalog.
  4. Match-rate = % EQ ASVs that are EMP-cosmopolitan (>=25 EMP samples globally).
  5. Per-rank breakdown via taxonomy.

Outputs:
  cache/emp_cosmopolitanism/emp_*.fa  (downloaded EMP fastas)
  cache/emp_cosmopolitanism/eq_v4_trimmed.fasta
  cache/emp_cosmopolitanism/eq_vs_emp_match.tsv
  cache/emp_cosmopolitanism/cosmopolitanism_summary.tsv
  cache/emp_cosmopolitanism/summary.txt
"""
from __future__ import annotations

import sys
import urllib.request
import subprocess
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "emp_cosmopolitanism"
OUT.mkdir(parents=True, exist_ok=True)

EMP_BASE = "https://ftp.microbio.me/emp/release1/otu_info/deblur/"
EMP_FILES = {
    90:  "emp.90.min25.deblur.seq.fa",
    100: "emp.100.min25.deblur.seq.fa",
    150: "emp.150.min25.deblur.seq.fa",
}

V4_FWD_PRIMER = "GTGYCAGCMGCCGCGGTAA"  # 515F (Earth Microbiome Project)
# Our EQ primers are 341F = CCTACGGGNGGCWGCAG, 805R = GACTACHVGGGTATCTAATCC
# So the V3-V4 EQ ASV (sub-440 bp from forward primer) overlaps with V4 EMP region
# starting at the 515F primer, which is ~170 bp into the V3-V4 product.
# (V3 spans 341-515, V4 starts at 515.)
# The 515F primer sequence is typically present near position 170-200 in our
# V3-V4 reads.

# Define a regex that matches the V4 forward primer region with degeneracy
import re
DEGEN = {"Y": "[CT]", "R": "[AG]", "M": "[AC]", "K": "[GT]", "S": "[CG]",
         "W": "[AT]", "N": "[ACGT]", "V": "[ACG]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]"}
def to_regex(seq):
    return ''.join(DEGEN.get(c, c) for c in seq.upper())
PRIMER_RE = re.compile(to_regex(V4_FWD_PRIMER))


def download_emp():
    paths = {}
    for length, fname in EMP_FILES.items():
        dst = OUT / fname
        if dst.exists() and dst.stat().st_size > 1000:
            print(f"  cached: {fname}", flush=True)
        else:
            url = EMP_BASE + fname
            print(f"  downloading: {url}", flush=True)
            try:
                urllib.request.urlretrieve(url, dst)
                print(f"    done ({dst.stat().st_size:,} bytes)", flush=True)
            except Exception as e:
                print(f"    FAILED: {e}", flush=True)
                continue
        paths[length] = dst
    return paths


def find_eq_asv_fasta() -> Path:
    """Locate the EQ ASV reps fasta (not in repo, look in shared dir)."""
    candidates = [
        Path("/home/leechuck/Public/software/empty-quarter/data/processed/"
             "taxonomy/taxon-tables/ASV_seqs-trips1-5.fasta"),
        REPO / "data" / "taxonomy" / "ASV_seqs-trips1-5.fasta",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError("ASV_seqs-trips1-5.fasta not found")


def trim_eq_to_v4(eq_fasta: Path, length: int, out_fa: Path) -> int:
    """Read EQ ASVs, find 515F primer region in each, output 90/100/150 bp
    starting AFTER the primer match. Returns count of trimmed sequences."""
    n_total = 0; n_trimmed = 0
    with open(eq_fasta) as fin, open(out_fa, "w") as fout:
        hdr = None; seq_lines = []
        def emit():
            nonlocal n_total, n_trimmed
            n_total += 1
            seq = ''.join(seq_lines).upper()
            m = PRIMER_RE.search(seq)
            if m is None:
                # Fallback: take last `length` bp (V4 is at the end of V3-V4 read)
                if len(seq) < length: return
                v4 = seq[-length:]
            else:
                start = m.end()
                if len(seq) < start + length: return
                v4 = seq[start:start + length]
            fout.write(f">{hdr}\n{v4}\n")
            n_trimmed += 1
        for ln in fin:
            ln = ln.rstrip()
            if ln.startswith(">"):
                if hdr is not None and seq_lines:
                    emit()
                hdr = ln[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(ln)
        if hdr is not None and seq_lines:
            emit()
    return n_total, n_trimmed


def vsearch_match(query_fa: Path, db_fa: Path, out_tsv: Path,
                   pid: float = 1.0) -> int:
    """Exact-match query vs db using vsearch."""
    vsearch_bin = "vsearch"  # rely on system or .venv
    # Check if vsearch is available
    try:
        subprocess.run([vsearch_bin, "--version"], capture_output=True, timeout=5)
    except Exception:
        # try other paths
        for cand in ["/usr/bin/vsearch",
                     "/home/leechuck/.local/bin/vsearch"]:
            if Path(cand).exists():
                vsearch_bin = cand
                break
        else:
            raise RuntimeError("vsearch not found")
    cmd = [vsearch_bin, "--usearch_global", str(query_fa),
           "--db", str(db_fa),
           "--id", str(pid),
           "--strand", "both",
           "--top_hits_only",
           "--maxaccepts", "1",
           "--maxrejects", "32",
           "--blast6out", str(out_tsv),
           "--quiet"]
    subprocess.run(cmd, check=True, capture_output=True, timeout=900)
    if out_tsv.exists():
        return sum(1 for _ in open(out_tsv))
    return 0


def main():
    print("[step 1] download EMP deblur catalogs", flush=True)
    emp_paths = download_emp()
    if not emp_paths:
        print("FAILED to download any EMP fasta", flush=True)
        return

    eq_fa = find_eq_asv_fasta()
    print(f"\n[step 2] EQ ASV reps: {eq_fa}", flush=True)

    rows = []
    for length in (90, 100, 150):
        if length not in emp_paths: continue
        eq_v4 = OUT / f"eq_v4_{length}bp.fasta"
        n_total, n_trimmed = trim_eq_to_v4(eq_fa, length, eq_v4)
        print(f"\n[length {length}] trimmed {n_trimmed:,} of "
              f"{n_total:,} EQ ASVs to {length}bp V4", flush=True)

        match_tsv = OUT / f"eq_vs_emp_{length}bp.tsv"
        try:
            n_match = vsearch_match(eq_v4, emp_paths[length], match_tsv, pid=1.0)
        except Exception as e:
            print(f"vsearch failed for {length}bp: {e}", flush=True)
            continue
        # Count unique queries with hits (blast6 might report 1 hit per query)
        if match_tsv.exists() and match_tsv.stat().st_size > 0:
            qm = set()
            with open(match_tsv) as f:
                for ln in f:
                    fields = ln.split("\t")
                    if fields:
                        qm.add(fields[0])
            n_unique_match = len(qm)
        else:
            n_unique_match = 0
        rows.append({"v4_length_bp": length,
                     "n_eq_asvs_trimmed": n_trimmed,
                     "n_emp_seqs": sum(1 for ln in open(emp_paths[length])
                                        if ln.startswith(">")),
                     "n_eq_asvs_matched": n_unique_match,
                     "frac_eq_matched": n_unique_match / n_trimmed
                                          if n_trimmed else 0})
        print(f"  {n_unique_match:,} of {n_trimmed:,} EQ ASVs matched EMP "
              f"(>= 25 sample globally) at {length} bp ID=1.0  "
              f"({n_unique_match/n_trimmed:.1%})", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "cosmopolitanism_summary.tsv", sep="\t", index=False)

    with open(OUT / "summary.txt", "w") as fh:
        fh.write("Test 2: EMP cosmopolitanism (exact-ASV matching)\n")
        fh.write("=" * 70 + "\n\n")
        fh.write("EMP catalogs are 'min25' = present in >= 25 EMP samples globally,\n"
                 "so any EQ ASV match means that ASV is GLOBALLY COSMOPOLITAN.\n\n")
        fh.write(df.to_string(index=False))
        fh.write("\n\nINTERPRETATION KEY:\n")
        fh.write("  frac_eq_matched > 0.5 -> EQ microbiome is dominated by\n"
                 "    cosmopolitan ASVs; hyperarid is not a unique pool.\n")
        fh.write("  frac_eq_matched < 0.1 -> EQ has significant endemism;\n"
                 "    most ASVs are NOT in the global cosmopolitan pool.\n")
        fh.write("\nNote: this is a CONSERVATIVE cosmopolitanism estimate. EQ ASVs\n"
                 "absent from EMP-min25 might still be in EMP at <25 samples.\n"
                 "For full cosmopolitanism analysis, would need EMP BIOM table.\n")
    print(f"\nWrote {OUT}/summary.txt")


if __name__ == "__main__":
    main()
