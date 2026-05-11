#!/usr/bin/env python3
"""Build the comprehensive betA candidate FASTA for MSA + phylogeny.
Pulls top K00108 hit (any score) per MAG from EQ + GEM, plus curated
high-scoring controls and canonical UniProt references.
"""
import csv, sys, re, gzip
from pathlib import Path
from collections import defaultdict
ROOT=Path('/data/emptyquarter/ecology-paper-runs/public_metagenomes')
MSA=ROOT/'msa'; MSA.mkdir(exist_ok=True)

def parse_fasta(p):
    seqs={}
    name=None; buf=[]
    opener=gzip.open if str(p).endswith('.gz') else open
    with opener(p,'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name: seqs[name]=''.join(buf)
                name=line[1:].split()[0]; buf=[]
            else:
                buf.append(line.strip())
        if name: seqs[name]=''.join(buf)
    return seqs

def top_k00108(tbl):
    best=None
    if not tbl.exists(): return None
    for line in open(tbl):
        if line.startswith('#'): continue
        parts=line.split()
        if len(parts)<6: continue
        if parts[2]!='K00108': continue
        try: score=float(parts[5])
        except: continue
        if best is None or score>best[1]: best=(parts[0], score)
    return best

out_fasta=MSA/'candidate_betA.faa'
out_meta=MSA/'candidate_betA_meta.tsv'
records=[]

# --- EQ CSP1-2 MAGs (5) ---
eq_files={
    'CSP12_V27Dr2__SemiBin_73': 'CSP12_V27Dr2__SemiBin_73.faa',
    'CSP12_V30PRr1__SemiBin_7': 'CSP12_V30PRr1__SemiBin_7.faa',
    'CSP12_V32PRr1__SemiBin_26': 'CSP12_V32PRr1__SemiBin_26.faa',
    'CSP12_V38PRr3__SemiBin_38': 'CSP12_V38PRr3__SemiBin_38.faa',
    'coassembly_CSP12_SemiBin_370': 'coassembly_CSP12_SemiBin_370.faa',
}
for gid,fname in eq_files.items():
    faa=MSA/'eq_faa'/fname
    tbl=MSA/'eq_faa'/f'{gid}.K00108.tbl'
    hit=top_k00108(tbl)
    if not hit: print(f'[skip] {gid}: no K00108 hit'); continue
    seqs=parse_fasta(faa)
    if hit[0] not in seqs: print(f'[skip] {gid}: ORF {hit[0]} not in faa'); continue
    tag=f'EQ_CSP12|{gid}|s={hit[1]:.1f}'
    records.append((tag, seqs[hit[0]], 'EQ_CSP12', gid, hit[1], 'EQ_single_or_coassembly'))

# --- GEM MAGs: top K00108 hit per MAG ---
sel={r['genome_id']:r for r in csv.DictReader(open(ROOT/'metadata/selected_mags.tsv'), delimiter='\t')}
# We include: all A_dadabacteria and B_dependent_family with a hit; for C_soil_control include top-20 by score.
gem_candidates=[]
for gid, s in sel.items():
    tbl=ROOT/'hmm/per_genome'/f'{gid}.tbl'
    hit=top_k00108(tbl)
    if not hit: continue
    gem_candidates.append((gid, s, hit))

# Count
from collections import Counter
hit_groups=Counter([(g[1]['selection_group']) for g in gem_candidates])
print('GEM MAGs with K00108 hit:', hit_groups)

for gid, s, hit in gem_candidates:
    grp=s['selection_group']
    # All A_dadabacteria; all B_dependent_family with score>=100; top-20 C_soil_control by score
    if grp=='A_dadabacteria': include=True
    elif grp=='B_dependent_family' and hit[1]>=100: include=True
    elif grp=='C_soil_control': include='defer'
    else: include=False
    if not include: continue
    faa=ROOT/'orfs'/f'{gid}.faa'
    if not faa.exists(): continue
    seqs=parse_fasta(faa)
    if hit[0] not in seqs: continue
    # Short tax label for header
    phy=re.search(r'p__([^;]+)', s['taxonomy']); fam=re.search(r'f__([^;]+)', s['taxonomy'])
    taxlabel=f"{phy.group(1) if phy else '?'}/{fam.group(1) if fam else '?'}"
    tag=f'{grp}|{gid}|s={hit[1]:.0f}|{taxlabel}'
    if include is True:
        records.append((tag, seqs[hit[0]], grp, gid, hit[1], s['taxonomy']))
# Now top-20 C_soil_control
c_hits=sorted([(gid,s,hit) for gid,s,hit in gem_candidates if s['selection_group']=='C_soil_control'], key=lambda x: -x[2][1])[:20]
for gid,s,hit in c_hits:
    faa=ROOT/'orfs'/f'{gid}.faa'
    seqs=parse_fasta(faa)
    if hit[0] not in seqs: continue
    phy=re.search(r'p__([^;]+)', s['taxonomy']); fam=re.search(r'f__([^;]+)', s['taxonomy'])
    taxlabel=f"{phy.group(1) if phy else '?'}/{fam.group(1) if fam else '?'}"
    tag=f'C_soil_top20|{gid}|s={hit[1]:.0f}|{taxlabel}'
    records.append((tag, seqs[hit[0]], 'C_soil_top20', gid, hit[1], s['taxonomy']))

# Write
print(f'Writing {len(records)} sequences to {out_fasta}')
with open(out_fasta,'w') as fh, open(out_meta,'w') as mh:
    mh.write('tag\tgroup\tgenome_id\tK00108_score\ttaxonomy\tseq_length\n')
    for tag,seq,grp,gid,score,tx in records:
        fh.write(f'>{tag}\n{seq}\n')
        mh.write(f'{tag}\t{grp}\t{gid}\t{score}\t{tx}\t{len(seq)}\n')
# Summary
gc=Counter([r[2] for r in records])
print('Counts by group:', gc)
