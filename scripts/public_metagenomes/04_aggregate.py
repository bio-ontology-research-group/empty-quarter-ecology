#!/usr/bin/env python3
"""Aggregate hmmsearch tbl files -> betA matrix with KOfam trusted thresholds.
Skips KOs with '-' threshold (KOfam convention: no reliable threshold).
"""
import csv
from pathlib import Path
from collections import Counter, defaultdict
ROOT=Path('/data/emptyquarter/ecology-paper-runs/public_metagenomes')
KOTHR=Path('/data/emptyquarter/ecology-paper-runs/t1_mechanism/kofam/ko_thresholds.tsv')
thresholds={}
for line in open(KOTHR):
    parts=line.rstrip().split('\t')
    if len(parts)!=2: continue
    ko,thr=parts
    try: thresholds[ko]=float(thr)
    except ValueError: pass  # skip '-' etc.
print(f'[info] {len(thresholds)} KO thresholds loaded')
tbl_dir=ROOT/'hmm/per_genome'
tbls=list(tbl_dir.glob('*.tbl'))
print(f'[info] {len(tbls)} tbl files')
per_genome={}
for tbl in tbls:
    gid=tbl.stem
    counts=Counter()
    with open(tbl) as fh:
        for line in fh:
            if line.startswith('#'): continue
            parts=line.split()
            if len(parts)<6: continue
            ko=parts[2]
            try: score=float(parts[5])
            except ValueError: continue
            thr=thresholds.get(ko)
            if thr is not None and score>=thr:
                counts[ko]+=1
    per_genome[gid]=counts
# Reload thresholds only for the 29 target KOs in target.hmm
HMM=Path('/data/emptyquarter/ecology-paper-runs/t1_mechanism/kofam/target.hmm')
target_kos=set()
for line in open(HMM):
    if line.startswith('NAME'):
        target_kos.add(line.split()[1])
print(f'[info] {len(target_kos)} target KOs in target.hmm')
sel={r['genome_id']:r for r in csv.DictReader(open(ROOT/'metadata/selected_mags.tsv'),delimiter='\t')}
KOS=sorted(target_kos)
out=ROOT/'out/betA_matrix.tsv'
out.parent.mkdir(exist_ok=True)
with open(out,'w') as fh:
    w=csv.writer(fh,delimiter='\t')
    w.writerow(['genome_id','selection_group','taxonomy','C','X','ecosystem_category','ecosystem_type','habitat']+KOS)
    for gid in sorted(per_genome):
        s=sel.get(gid,{})
        row=[gid, s.get('selection_group',''), s.get('taxonomy',''), s.get('C',''), s.get('X',''),
             s.get('ecosystem_category',''), s.get('ecosystem_type',''), s.get('habitat','')]
        row+=[per_genome[gid].get(k,0) for k in KOS]
        w.writerow(row)
print(f'[done] {out}')
by_grp=defaultdict(lambda:[0,0])
for gid,c in per_genome.items():
    s=sel.get(gid)
    if not s: continue
    g=s['selection_group']
    by_grp[g][1]+=1
    if c.get('K00108',0)>0: by_grp[g][0]+=1
print()
print('betA (K00108) presence by selection group (trusted-threshold):')
for g in ('A_dadabacteria','B_dependent_family','C_soil_control'):
    pos,tot=by_grp[g]
    pct=100*pos/tot if tot else 0
    print(f'  {g:25s} {pos}/{tot} ({pct:.1f}%)')
