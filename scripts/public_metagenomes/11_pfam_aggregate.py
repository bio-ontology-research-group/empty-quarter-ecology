#!/usr/bin/env python3
"""Aggregate Pfam GMC_oxred_N + GMC_oxred_C scans per MAG.
A MAG is considered betA-capable in the Pfam sense if it has ANY protein
hitting both GMC_oxred_N and GMC_oxred_C above the Pfam gathering threshold.
"""
import csv
from pathlib import Path
from collections import Counter, defaultdict
ROOT=Path('/data/emptyquarter/ecology-paper-runs/public_metagenomes')
TBL=ROOT/'pfam/per_genome'
sel={r['genome_id']:r for r in csv.DictReader(open(ROOT/'metadata/selected_mags.tsv'),delimiter='\t')}

per={}
for t in TBL.glob('*.tbl'):
    gid=t.stem
    hits_by_protein=defaultdict(set)
    n_hits=0
    for line in open(t):
        if line.startswith('#'): continue
        parts=line.split()
        if len(parts)<6: continue
        target=parts[0]; query=parts[2]
        hits_by_protein[target].add(query)
        n_hits+=1
    # A protein counts if it has EITHER Pfam hit (hmmsearch --cut_ga filters by GA)
    n_any=sum(1 for p,qs in hits_by_protein.items() if qs)
    # Stricter: both domains (full-length GMC oxidoreductase)
    n_both=sum(1 for p,qs in hits_by_protein.items() if 'GMC_oxred_N' in qs and 'GMC_oxred_C' in qs)
    per[gid]={'n_any_GMC':n_any,'n_both_N_and_C':n_both,'n_hits':n_hits}

rows=[]
for gid, s in sel.items():
    p=per.get(gid,{'n_any_GMC':0,'n_both_N_and_C':0})
    rows.append({'genome_id':gid,'group':s['selection_group'],
                 'any_GMC_protein':p['n_any_GMC']>0,
                 'full_GMC_oxidoreductase_protein':p['n_both_N_and_C']>0,
                 'n_GMC_proteins':p['n_any_GMC']})
out=ROOT/'out/pfam_per_genome.tsv'
out.parent.mkdir(exist_ok=True)
with open(out,'w') as fh:
    w=csv.DictWriter(fh,delimiter='\t',fieldnames=['genome_id','group','any_GMC_protein','full_GMC_oxidoreductase_protein','n_GMC_proteins'])
    w.writeheader(); w.writerows(rows)
print('Wrote', out)

# Summary
summary=[]
for g in ('A_dadabacteria','B_dependent_family','C_soil_control'):
    sub=[r for r in rows if r['group']==g]
    n=len(sub)
    any_g=sum(1 for r in sub if r['any_GMC_protein'])
    both=sum(1 for r in sub if r['full_GMC_oxidoreductase_protein'])
    summary.append({'selection_group':g,'n':n,'any_GMC_protein':any_g,
                    'full_GMC_ox_protein':both,
                    'pct_any':round(100*any_g/n,1),'pct_both':round(100*both/n,1)})
out2=ROOT/'out/pfam_summary.tsv'
with open(out2,'w') as fh:
    w=csv.DictWriter(fh,delimiter='\t',fieldnames=list(summary[0].keys()))
    w.writeheader(); w.writerows(summary)
print('Wrote', out2)
for s in summary: print(s)
