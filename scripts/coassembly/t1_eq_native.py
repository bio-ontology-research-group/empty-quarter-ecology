#!/usr/bin/env python3
"""T1 inventory v2: uses EQ-native Prokka-annotated MAGs for the 7 dependent
genera instead of NCBI PGAP reference genomes, and the new 93%-complete
co-assembly CSP1-2 MAG (SemiBin_370) instead of the older single-sample MAGs.
"""
import re, sys
from pathlib import Path
from collections import Counter
import pandas as pd

# Genes — exact gene= tag match from Prokka GFF (case-insensitive)
BIOSYNTH = {
    'otsA': r'^otsA$',
    'otsB': r'^otsB$',
    'treS': r'^treS$',
    'treY': r'^treY$',
    'treZ': r'^treZ$',
    'treT': r'^treT$',
    'betA': r'^betA$',
    'betB': r'^betB$',
    'gbsA': r'^gbsA$',
    'gbsB': r'^gbsB$',
    'ectA': r'^ectA$',
    'ectB': r'^ectB$',
    'ectC': r'^ectC$',
    'ectD': r'^ectD$',
    'ino1': r'^ino1$',
}
UPTAKE = {
    'opuA':  r'^opuA[ABC]?$',
    'opuB':  r'^opuB[A-D]?$',
    'opuC':  r'^opuC[A-D]?$',
    'opuD':  r'^opuD$',
    'opuE':  r'^opuE$',
    'proV':  r'^proV$',
    'proW':  r'^proW$',
    'proX':  r'^proX$',
    'proP':  r'^proP$',
    'betT':  r'^betT$',
    'osmY':  r'^osmY$',
    'treA':  r'^treA$',
    'treB':  r'^treB$',
    'treC':  r'^treC$',
    'treR':  r'^treR$',
    'treP':  r'^treP$',
    'bccT':  r'^bccT$',
    'betP':  r'^betP$',
}

BIO_C = {k: re.compile(p, re.IGNORECASE) for k, p in BIOSYNTH.items()}
UPT_C = {k: re.compile(p, re.IGNORECASE) for k, p in UPTAKE.items()}

def count_prokka_gff(gff):
    counts_bio = Counter()
    counts_upt = Counter()
    with open(gff) as fh:
        for line in fh:
            if line.startswith('#') or '\t' not in line:
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9 or cols[2] != 'CDS':
                continue
            m = re.search(r'gene=([^;_]+)', cols[8])
            if not m:
                continue
            gene = m.group(1).strip()
            # strip Prokka's underscore-number duplicate suffix like otsA_1
            gene = re.sub(r'_\d+$', '', gene)
            for k, r in BIO_C.items():
                if r.match(gene):
                    counts_bio[k] += 1; break
            for k, r in UPT_C.items():
                if r.match(gene):
                    counts_upt[k] += 1; break
    return counts_bio, counts_upt

PROKKA_DIR = Path('/data/emptyquarter/ecology-paper-runs/csp_mag/prokka_out')
# Known bin classifications (from GTDB-Tk)
GROUPS = {
    'coassembly_CSP12_SemiBin_370': ('CSP1-2 co-assembly (93% complete)', True),
    'Herpetosiphonaceae_BPHZ01_1_SemiBin_120': ('Herpetosiphonaceae EQ-native 1', False),
    'Herpetosiphonaceae_BPHZ01_2_SemiBin_145': ('Herpetosiphonaceae EQ-native 2', False),
    'Rubrobacter_1_SemiBin_66': ('Rubrobacter EQ-native 1', False),
    'Rubrobacter_JALYGY01_SemiBin_234': ('Rubrobacter EQ-native 2', False),
    'Rubrobacter_DATHTK01_SemiBin_9': ('Rubrobacter EQ-native 3', False),
    'Streptomyces_SemiBin_236': ('Streptomyces EQ-native', False),
    'Myxococcus_SemiBin_124': ('Myxococcus EQ-native', False),
    'CSP1-4_Limnocylindria_SemiBin_191': ('CSP1-4 EQ-native (candidate phylum)', None),
    'CSP1-6_Rokubacteriales_SemiBin_69': ('CSP1-6 EQ-native (candidate phylum)', None),
}
rows = []
for subdir, (label, is_csp) in GROUPS.items():
    gff = PROKKA_DIR / subdir / f'{subdir}.gff'
    if not gff.exists():
        print(f'skip {subdir}'); continue
    bio, upt = count_prokka_gff(gff)
    row = {'subdir': subdir, 'label': label, 'is_csp_mag': is_csp,
           'biosynth_distinct': len([k for k,v in bio.items() if v>0]),
           'biosynth_total':    sum(bio.values()),
           'biosynth_genes':    ','.join(sorted([k for k,v in bio.items() if v>0])) or '—',
           'uptake_distinct':   len([k for k,v in upt.items() if v>0]),
           'uptake_total':      sum(upt.values()),
           'uptake_genes':      ','.join(sorted([k for k,v in upt.items() if v>0])) or '—',
          }
    for k in BIOSYNTH: row[f'b_{k}'] = bio.get(k, 0)
    for k in UPTAKE:   row[f'u_{k}'] = upt.get(k, 0)
    rows.append(row)

df = pd.DataFrame(rows)
OUT = Path('/data/emptyquarter/ecology-paper-runs/csp_mag/t1_eq_native_inventory.tsv')
df.to_csv(OUT, sep='\t', index=False)

print()
print('='*100)
print('T1 inventory — EQ-native MAGs (Prokka-annotated, consistent toolchain)')
print('='*100)
print(df[['label','biosynth_distinct','biosynth_genes','uptake_distinct','uptake_genes']].to_string(index=False))
print()
print('='*100)
print('Per-gene detail — BIOSYNTHESIS (key: otsA, otsB, betA, betB)')
print('='*100)
bio_cols = [f'b_{k}' for k in ['otsA','otsB','treS','treY','treZ','betA','betB','ectA','ectB','ectC','ectD']]
print(df[['label'] + bio_cols].to_string(index=False))
print()
print('='*100)
print('Per-gene detail — UPTAKE')
print('='*100)
upt_cols = [f'u_{k}' for k in ['opuA','opuB','opuC','opuD','opuE','proV','proW','proX','proP','betT','treA','treB','treC','treR','bccT','betP']]
print(df[['label'] + upt_cols].to_string(index=False))
