import pandas as pd,numpy as np,json,hashlib,argparse
from pathlib import Path
parser=argparse.ArgumentParser(description='Reconstruct the archaeal mean read fraction from canonical counts.');parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[3]);root=parser.parse_args().root.resolve()
tax=root/'data/metadata/taxonomy/taxonomy-trips1-5.tsv';alpha=root/'analysis/v2/review/cache/alpha.tsv';counts=root/'data/metadata/taxonomy/feature-table-trips1-5.tsv'
t=pd.read_csv(tax,sep='\t',index_col=0);ids=set(t.index[t.Taxon.str.startswith('Archaea;')]);a=pd.read_csv(alpha,sep='\t',index_col=0);a=a[a.Site.between(1,60)];assert len(a)==1227
sums=np.zeros(len(a));depth=np.zeros(len(a))
for c in pd.read_csv(counts,sep='\t',skiprows=1,index_col=0,chunksize=5000):
 c=c.loc[:,a.index];depth+=c.sum(axis=0).to_numpy();sums+=c.loc[c.index.isin(ids)].sum(axis=0).to_numpy()
assert np.array_equal(depth,a.depth.to_numpy())
pct=float((sums/depth).mean()*100)
out={'method':'Mean archaeal read fraction across 1227 core profiles, reconstructed directly from canonical ASV counts and SILVA domain assignments. All profile depths match cached alpha depths exactly.','mean_archaeal_read_share_pct':pct,'archaeal_asvs':len(ids),'n_profiles':len(a),'inputs':{str(p.relative_to(root)):hashlib.file_digest(p.open('rb'),'sha256').hexdigest() for p in [tax,alpha,counts]}}
(root/'analysis/v3/manuscript_validation_20261004/archaea_share.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
