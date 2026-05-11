#!/usr/bin/env python3
"""Project E. coli P17444 (BetA_ECOLI) catalytic residues onto the MSA and
tabulate per-group conservation.

E. coli betA functional residues (literature; UniProt P17444 features):
  - His466   catalytic proton acceptor
  - Ser384   substrate binding (choline hydroxyl)
  - Asn120   FAD binding region
  - Tyr465   part of catalytic His-Tyr pair typical of GMC oxidoreductases

We locate P17444 in the alignment, map its residue positions to alignment
columns, then count per-group residue identity at those columns.
"""
import csv
from pathlib import Path
from collections import defaultdict, Counter
MSA=Path('/data/emptyquarter/ecology-paper-runs/public_metagenomes/msa')
ALN=MSA/'betA_combined.aln.faa'
OUT=MSA/'active_site_conservation.tsv'

if not ALN.exists():
    print(f'[err] missing {ALN}'); raise SystemExit(1)

seqs={}
name=None; buf=[]
for line in open(ALN):
    if line.startswith('>'):
        if name: seqs[name]=''.join(buf)
        name=line[1:].strip(); buf=[]
    else: buf.append(line.strip())
if name: seqs[name]=''.join(buf)
print(f'{len(seqs)} sequences in alignment')

# Find E. coli P17444 sequence
ecoli_key=next((k for k in seqs if 'P17444' in k or 'BETA_ECOLI' in k), None)
if not ecoli_key:
    print('[err] E. coli P17444 not found in alignment. Available refs:')
    for k in seqs:
        if 'REF_' in k: print('  ',k)
    raise SystemExit(1)
print(f'Anchor: {ecoli_key}')

# Map E. coli residue positions -> alignment column
ecoli_aln=seqs[ecoli_key]
resno_to_col={}
resno=0
for col, aa in enumerate(ecoli_aln):
    if aa != '-':
        resno += 1
        resno_to_col[resno] = col

# Targets (E. coli numbering; protein starts at M=1)
targets = {
    'Asn120': (120, 'N'),
    'Ser384': (384, 'S'),
    'Tyr465': (465, 'Y'),
    'His466': (466, 'H'),
}
# Confirm anchor residues
for name_res, (pos, expected) in targets.items():
    col = resno_to_col.get(pos)
    if col is None:
        print(f'[warn] E. coli residue {pos} not found in alignment')
        continue
    obs = ecoli_aln[col]
    print(f'  E. coli {name_res}: align-col {col}, residue={obs} (expected {expected})')

# Extract group assignment from header
def group_of(h):
    if h.startswith('EQ_CSP12'): return 'EQ_CSP12'
    if h.startswith('A_dadabacteria'): return 'A_dadabacteria'
    if h.startswith('B_dependent_family'): return 'B_dependent_family'
    if h.startswith('C_soil_top20'): return 'C_soil_top20'
    if h.startswith('REF_betA'): return 'REF_betA'
    return 'other'

# Count residue at each target position per group
print()
print('===== Per-group residue conservation at catalytic positions =====')
rows=[]
for tname,(pos, expected) in targets.items():
    col = resno_to_col.get(pos)
    if col is None: continue
    # Per group counts
    for grp in ('EQ_CSP12','A_dadabacteria','B_dependent_family','C_soil_top20','REF_betA'):
        cnt=Counter()
        for h,s in seqs.items():
            if group_of(h) != grp: continue
            aa=s[col] if col < len(s) else '-'
            cnt[aa]+=1
        tot=sum(cnt.values())
        if tot==0: continue
        conserved=cnt.get(expected,0)
        rows.append((tname, pos, expected, grp, tot, conserved, round(100*conserved/tot,1), dict(cnt.most_common())))
        print(f'  {tname} (col {col}) {grp:25s} {conserved}/{tot} ({100*conserved/tot:.1f}%) top residues: {dict(cnt.most_common(5))}')
# Save
with open(OUT,'w') as fh:
    w=csv.writer(fh,delimiter='\t')
    w.writerow(['residue','ecoli_pos','expected_aa','group','n','conserved','percent_conserved','residue_distribution'])
    for r in rows: w.writerow(r)
print(f'\nWrote {OUT}')
