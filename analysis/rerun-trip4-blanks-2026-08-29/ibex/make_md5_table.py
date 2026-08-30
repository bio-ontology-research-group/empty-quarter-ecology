#!/usr/bin/env python3
"""Rewrite an nf-core/ampliseq dada2/ASV_table.tsv with md5-of-sequence ids.

Usage: make_md5_table.py ASV_table.tsv ASV_seqs.fasta out.tsv
The canonical Trips 1-5 table (a QIIME 2 export) keys ASVs by md5 of the
uppercase sequence; this makes the blank table joinable on '#OTU ID'.
"""
import hashlib
import sys


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
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


def main(table, fasta, out):
    seqs = read_fasta(fasta)
    n = dup = 0
    seen = {}
    with open(table) as fh, open(out, "w") as oh:
        header = fh.readline().rstrip("\n").split("\t")
        oh.write("#OTU ID\t" + "\t".join(header[1:]) + "\n")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            asv = parts[0]
            seq = seqs[asv].upper()
            md5 = hashlib.md5(seq.encode()).hexdigest()
            if md5 in seen:
                dup += 1
            seen[md5] = asv
            oh.write(md5 + "\t" + "\t".join(parts[1:]) + "\n")
            n += 1
    print(f"wrote {out}: {n} ASVs, {len(header) - 1} samples, {dup} duplicate md5 ids", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:4])
