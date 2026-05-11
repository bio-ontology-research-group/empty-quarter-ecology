#!/usr/bin/env python3
"""Select public GEM MAGs for the public-metagenome betA asymmetry test.
Output: selected_mags.tsv with columns genome_id, taxonomy, C, X, habitat, selection_group.
Criteria (C>=70, X<=10):
  A) All Dadabacteria (CSP1-2 / Desulfobacterota_D clade), any habitat
  B) All dependent-family MAGs (Herpetosiphonaceae, Rubrobacteraceae, Streptomycetaceae,
     Paenibacillaceae, Myxococcaceae), any habitat
  C) All other Terrestrial Soil MAGs as phylum-level controls
"""
import csv, sys
from pathlib import Path
ROOT=Path('/data/emptyquarter/ecology-paper-runs/public_metagenomes')
rows=[r for r in csv.DictReader(open(ROOT/'metadata/genome_metadata.tsv'), delimiter='\t')]
def qual(r):
    try: return float(r['completeness'])>=70 and float(r['contamination'])<=10
    except: return False
rows=[r for r in rows if qual(r)]
FAMS={'Herpetosiphonaceae','Rubrobacteraceae','Streptomycetaceae','Paenibacillaceae','Myxococcaceae'}
def sel_group(r):
    tax=r['ecosystem']
    if 'Dadabacteria' in tax: return 'A_dadabacteria'
    if any('f__'+f in tax for f in FAMS): return 'B_dependent_family'
    if r['ecosystem_category']=='Terrestrial' and r['ecosystem_type']=='Soil': return 'C_soil_control'
    return None
out=[]
for r in rows:
    g=sel_group(r)
    if g: 
        out.append({'genome_id':r['genome_id'],'taxonomy':r['ecosystem'],
                    'C':r['completeness'],'X':r['contamination'],
                    'ecosystem_category':r['ecosystem_category'],
                    'ecosystem_type':r['ecosystem_type'],
                    'habitat':r['habitat'],'lat':r['latitude'],'lon':r['longitude'],
                    'selection_group':g})
from collections import Counter
print('Selection counts:', Counter([r['selection_group'] for r in out]))
print('Total unique:', len(out))
out_path=ROOT/'metadata/selected_mags.tsv'
with open(out_path,'w') as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0].keys()),delimiter='\t')
    w.writeheader(); w.writerows(out)
print('Wrote', out_path)
