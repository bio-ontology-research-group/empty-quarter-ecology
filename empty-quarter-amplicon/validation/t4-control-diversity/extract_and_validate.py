from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,sys
from scipy.stats import spearmanr
A=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
r=A/'analysis/rerun-controls-2026-08-30/screen/outputs';ledger=pd.read_csv(r/'removal_fraction_by_profile.tsv',sep='\t');ledger=ledger[(ledger.trip=='Trip4') & (ledger.role=='compatible_biological_profile')].set_index('profile_id');assert len(ledger)==95
calls=pd.read_csv(r/'primary_contaminant_calls.tsv',sep='\t');ids=set(calls[calls.screen=='Trip4'].feature_id);assert len(ids)==7
alpha=pd.read_csv(A/'analysis/v2/review/cache/alpha.tsv',sep='\t',index_col=0).loc[ledger.index]
counts={}
with open(A/'data/processed/taxonomy/taxon-tables/feature-table-trips1-5.tsv') as f:
 next(f);header=next(f).rstrip().split('\t');indices=[header.index(x) for x in ledger.index]
 for line in f:
  key=line.split('\t',1)[0]
  if key in ids:
   row=line.rstrip().split('\t');counts[key]=[float(row[i]) for i in indices]
pd.DataFrame(counts,index=ledger.index).to_csv(out/'removed_asv_counts.tsv',sep='\t',index_label='profile');alpha[['depth','richness_raw','shannon']].to_csv(out/'before.tsv',sep='\t',index_label='profile');C=np.array(list(counts.values()));removed=C.sum(axis=0);present=(C>0).sum(axis=0)
assert np.array_equal(removed,ledger.candidate_contaminant_reads);assert np.array_equal(present,ledger.candidate_contaminant_asvs);assert np.array_equal(alpha.depth,ledger.total_reads);assert np.array_equal(alpha.richness_raw,ledger.total_asvs)
N=alpha.depth.to_numpy();H=alpha.shannon.to_numpy();Nr=N-removed
post=np.log(Nr)-(N*(np.log(N)-H)-np.where(C>0,C*np.log(np.maximum(C,1)),0).sum(axis=0))/Nr
raw=alpha.richness_raw.to_numpy();postraw=raw-present
# Independently reconstruct full-depth entropy and richness directly from counts.
depth_direct=np.zeros(95); raw_direct=np.zeros(95); clog_direct=np.zeros(95)
for chunk in pd.read_csv(A/'data/processed/taxonomy/taxon-tables/feature-table-trips1-5.tsv',sep='\t',skiprows=1,usecols=list(ledger.index),chunksize=20000):
 X=chunk.loc[:,ledger.index].to_numpy(dtype=float)
 depth_direct+=X.sum(axis=0);raw_direct+=(X>0).sum(axis=0);clog_direct+=(X*np.log(np.maximum(X,1))).sum(axis=0)
H_direct=np.log(depth_direct)-clog_direct/depth_direct
assert np.array_equal(depth_direct,N);assert np.array_equal(raw_direct,raw);assert np.allclose(H_direct,H,rtol=0,atol=1e-12)
post_direct=np.log(Nr)-(clog_direct-np.where(C>0,C*np.log(np.maximum(C,1)),0).sum(axis=0))/Nr
assert np.allclose(post_direct,post,rtol=0,atol=1e-12)
res=pd.DataFrame({'profile':ledger.index,'depth_before':N,'removed_reads':removed,'shannon_before':H,'shannon_after':post,'richness_before':raw,'richness_after':postraw});res.to_csv(out/'t4_before_after.tsv',sep='\t',index=False)
summary={'n_profiles':95,'called_asvs':7,'called_asvs_present_in_canonical':len(counts),'shannon_spearman':float(spearmanr(H,post).statistic),'raw_richness_spearman':float(spearmanr(raw,postraw).statistic),'max_abs_shannon_change':float(np.max(abs(post-H))),'direct_count_entropy_max_error':float(np.max(abs(H_direct-H))),'removed_reads_ledger_exact_match':True,'depth_and_richness_ledger_exact_match':True};(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

files=[r/'removal_fraction_by_profile.tsv',r/'primary_contaminant_calls.tsv',A/'analysis/v2/review/cache/alpha.tsv',A/'data/processed/taxonomy/taxon-tables/feature-table-trips1-5.tsv']
manifest={str(f.relative_to(A)):hashlib.file_digest(open(f,'rb'),'sha256').hexdigest() for f in files}
(out/'input-sha256.json').write_text(json.dumps(manifest,indent=2))
