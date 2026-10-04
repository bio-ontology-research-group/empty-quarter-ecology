#!/usr/bin/env python3
"""Self-contained common-denominator pH/geography association partition."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SEED = 20261004
FRACTIONS = ['unique_ph', 'unique_geography', 'shared', 'unexplained']


def verify_inputs():
    manifest = json.loads((HERE / 'input_manifest.json').read_text())
    for item in manifest['files']:
        path = HERE / item['packaged_path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item['packaged_sha256'], path


def load_cohort():
    """Rebuild campaign-site-compartment linkage and the frozen top-200 response."""
    a = pd.read_csv(HERE / 'inputs/alpha.tsv', sep='\t', index_col=0)
    a = a[pd.to_numeric(a.Site, errors='coerce').isin(range(1, 61))].copy()
    a['site'] = a.Site.astype(int)
    a['trip'] = a.Trip.astype(int)
    a['compartment'] = a.Type
    keys = ['trip', 'site', 'compartment']
    ph = pd.read_csv(HERE / 'inputs/accepted_ph.tsv', sep='\t')
    assert ph.disposition.eq('ADMITTED_MEASUREMENT').all()
    assert not ph.sample_id.duplicated().any()
    eligible = a[keys].drop_duplicates()
    linked = ph.merge(eligible, on=keys, validate='many_to_one')
    groups = linked.groupby(keys, sort=True).agg(ph=('ph_value', 'mean'), n_ph_specimens=('ph_value', 'size')).reset_index()
    c = pd.read_csv(HERE / 'inputs/genus_counts.tsv.gz', sep='\t', index_col=0)
    c = c.loc[~c.index.isna()]
    # merge discards original index; select profiles directly using eligible group tuples.
    group_keys = set(map(tuple, groups[keys].to_numpy()))
    a = a[[tuple(row) in group_keys for row in a[keys].to_numpy()]]
    columns = [column for column in c.columns if column in a.index]
    matrix = c[columns].T.copy()
    for key in keys:
        matrix[key] = a.loc[columns, key].to_numpy()
    counts = matrix.groupby(keys, sort=True)[c.index.tolist()].sum()
    counts = counts.loc[counts.sum(axis=1) >= 2000]
    cohort = counts.index.to_frame(index=False).merge(groups, on=keys, validate='one_to_one')
    coords = pd.read_csv(HERE / 'inputs/site_coordinates.tsv', sep='\t')
    cohort = cohort.merge(coords[['site', 'transect_km']], on='site', validate='many_to_one')
    assert list(map(tuple,cohort[keys].to_numpy())) == list(counts.index), 'Count/cohort row alignment changed'
    assert len(cohort) == 560 and cohort.site.nunique() == 60
    assert np.isfinite(cohort[['ph', 'transect_km']]).all().all()
    prevalence = (counts > 0).mean(axis=0)
    relative = counts.div(counts.sum(axis=1), axis=0)
    ranked = relative.loc[:, prevalence >= .2].mean(axis=0).sort_values(ascending=False, kind='mergesort')
    taxa = ranked.index[:200].tolist()
    values = counts[taxa].to_numpy(float)
    logged = np.log(values + .5)
    response = logged - logged.mean(axis=1, keepdims=True)
    ref = pd.read_csv(HERE / 'inputs/reference_groups.tsv', sep='\t')
    check = cohort.merge(ref[keys + ['ph']], on=keys, suffixes=('', '_reference'), validate='one_to_one')
    assert len(check) == len(cohort), 'Reference cohort does not cover every retained group'
    np.testing.assert_allclose(check.ph, check.ph_reference, rtol=0, atol=1e-12)
    return cohort, response, taxa, counts.sum(axis=1).to_numpy(), len(columns)


def design(cohort):
    category = cohort.trip.astype(str) + '|' + cohort.compartment
    nuisance = np.column_stack([np.ones(len(cohort)), pd.get_dummies(category, drop_first=True, dtype=float)])
    route = cohort.transect_km.to_numpy()
    z = (route - route.mean()) / route.std()
    ph = cohort.ph.to_numpy()
    ph = (ph - ph.mean()) / ph.std()
    x = np.column_stack([nuisance, ph, z, z*z])
    k = nuisance.shape[1]
    columns = {'N': np.arange(k), 'P': np.arange(k+1), 'G': np.r_[np.arange(k), k+1, k+2], 'PG': np.arange(k+3)}
    return x, columns


def direct_sse(x, y, weights):
    sx, sy = x*np.sqrt(weights)[:, None], y*np.sqrt(weights)[:, None]
    return float(np.square(sy - sx @ np.linalg.lstsq(sx, sy, rcond=None)[0]).sum())


def partition_from_sse(sse):
    n, p, g, pg = [sse[key] for key in ['N', 'P', 'G', 'PG']]
    assert n > 0
    fractions = dict(zip(FRACTIONS, [(g-pg)/n, (p-pg)/n, (n-p-g+pg)/n, pg/n]))
    assert abs(sum(fractions.values()) - 1) < 1e-10
    assert fractions['unique_ph'] > -1e-9 and fractions['unique_geography'] > -1e-9
    return {**fractions, 'ph_including_shared': (n-p)/n, 'geography_including_shared': (n-g)/n,
            'joint_explained': (n-pg)/n, **{'sse_'+key: value for key, value in sse.items()}}


class SiteStatistics:
    """Sufficient statistics retain full sites and fixed within-site weights."""
    def __init__(self, x, y, sites, weights, columns):
        self.columns = columns
        self.sites = np.unique(sites)
        self.xx, self.xy, self.yy = [], [], []
        for site in self.sites:
            mask = sites == site
            sx, sy = x[mask]*np.sqrt(weights[mask, None]), y[mask]*np.sqrt(weights[mask, None])
            self.xx.append(sx.T @ sx)
            self.xy.append(sx.T @ sy)
            self.yy.append(np.square(sy).sum())
        self.xx, self.xy, self.yy = map(np.asarray, [self.xx, self.xy, self.yy])

    def evaluate(self, multiplicity=None):
        if multiplicity is None:
            multiplicity = np.ones(len(self.sites))
        xx = np.einsum('s,sij->ij', multiplicity, self.xx)
        xy = np.einsum('s,sij->ij', multiplicity, self.xy)
        yy = float(multiplicity @ self.yy)
        sse, ranks = {}, {}
        for name, cols in self.columns.items():
            gram = xx[np.ix_(cols, cols)]
            cross = xy[cols]
            eig, vectors = np.linalg.eigh(gram)
            keep = eig > eig[-1]*1e-12
            ranks[name] = int(keep.sum())
            inverse = (vectors[:, keep]/eig[keep]) @ vectors[:, keep].T
            sse[name] = yy - float(np.sum(cross*(inverse @ cross)))
        return partition_from_sse(sse), ranks


class SiteAverageStatistics:
    """Refit common group nuisance then average Y and predictor residuals by site."""
    def __init__(self, x, y, sites, weights, columns):
        self.k=len(columns['N'])
        self.sites=np.unique(sites)
        nuisance=x[:,:self.k]
        joint=np.column_stack([y,x[:,self.k:]])
        self.dim=y.shape[1]
        self.xx=[]; self.xy=[]; self.mean_n=[]; self.mean_joint=[]
        for site in self.sites:
            mask=sites==site
            n=nuisance[mask]; z=joint[mask]; w=weights[mask]
            self.xx.append(n.T@(n*w[:,None]))
            self.xy.append(n.T@(z*w[:,None]))
            self.mean_n.append(n.mean(axis=0)); self.mean_joint.append(z.mean(axis=0))
        self.xx,self.xy,self.mean_n,self.mean_joint=map(np.asarray,[self.xx,self.xy,self.mean_n,self.mean_joint])

    def evaluate(self,multiplicity=None):
        if multiplicity is None: multiplicity=np.ones(len(self.sites))
        gram=np.einsum('s,sij->ij',multiplicity,self.xx)
        cross=np.einsum('s,sij->ij',multiplicity,self.xy)
        beta=np.linalg.pinv(gram,rcond=1e-12)@cross
        residual=self.mean_joint-self.mean_n@beta
        y=residual[:,:self.dim]
        x=np.column_stack([np.ones(len(y)),residual[:,self.dim:]])
        cols={'N':np.array([0]),'P':np.array([0,1]),'G':np.array([0,2,3]),'PG':np.arange(4)}
        ranks={key:int(np.linalg.matrix_rank(x[:,c]*np.sqrt(multiplicity)[:,None])) for key,c in cols.items()}
        sse={key:direct_sse(x[:,c],y,multiplicity) for key,c in cols.items()}
        return partition_from_sse(sse),ranks


def identifiable(rank):
    return (rank['P']-rank['N'], rank['G']-rank['N'], rank['PG']-rank['N']) == (1, 2, 3)


def draw_sites(rng, order, block_length=1):
    n = len(order)
    if block_length == 1:
        sampled = rng.integers(n, size=n)
    else:
        starts = rng.integers(n-block_length+1, size=int(np.ceil(n/block_length)))
        sampled = np.concatenate([order[start:start+block_length] for start in starts])[:n]
    return np.bincount(sampled, minlength=n)


def run(output, replicates):
    verify_inputs()
    output.mkdir(parents=True, exist_ok=False)
    cohort, y, taxa, reads, profile_count = load_cohort()
    x, columns = design(cohort)
    sites = cohort.site.to_numpy()
    size = cohort.groupby('site').site.transform('size').to_numpy()
    weights = 1/size
    stats = SiteStatistics(x, y, sites, weights, columns)
    primary, ranks = stats.evaluate()
    expected = {key: len(value) for key, value in columns.items()}
    assert ranks == expected, (ranks, expected)
    direct = partition_from_sse({key: direct_sse(x[:, cols], y, weights) for key, cols in columns.items()})
    for key in primary:
        np.testing.assert_allclose(primary[key], direct[key], rtol=1e-9, atol=1e-9)
    between_stats=SiteAverageStatistics(x,y,sites,weights,columns)
    between,between_ranks=between_stats.evaluate()
    assert identifiable(between_ranks)
    rows = [{'analysis': 'primary_equal_site_weight', **primary}, {'analysis':'site_averaged_common_adjustment',**between}]
    unweighted, _ = SiteStatistics(x, y, sites, np.ones(len(sites)), columns).evaluate()
    rows.append({'analysis': 'sensitivity_equal_group_weight', **unweighted})
    pd.DataFrame(rows).to_csv(output / 'partitions.tsv', sep='\t', index=False, float_format='%.12g')
    cohort['group_reads'] = reads
    cohort['primary_weight'] = weights
    cohort.to_csv(output / 'cohort.tsv', sep='\t', index=False, float_format='%.12g')
    pd.Series(taxa, name='genus').to_csv(output / 'selected_taxa.tsv', sep='\t', index=False)
    # Geographic partial collinearity is computed after precisely the same nuisance projection.
    k = len(columns['N'])
    ph_n = direct_sse(x[:, columns['N']], x[:, k:k+1], weights)
    ph_g = direct_sse(x[:, columns['G']], x[:, k:k+1], weights)
    collinearity = {'ph_R2_geography_given_nuisance': 1-ph_g/ph_n, 'ph_VIF_given_nuisance': ph_n/ph_g}
    beta=np.linalg.pinv(between_stats.xx.sum(axis=0),rcond=1e-12)@between_stats.xy.sum(axis=0)
    adjusted=between_stats.mean_joint-between_stats.mean_n@beta
    site_ph=adjusted[:,y.shape[1]:y.shape[1]+1]
    site_g=adjusted[:,y.shape[1]+1:]
    one=np.ones((60,1))
    site_ph_n=direct_sse(one,site_ph,np.ones(60))
    site_ph_g=direct_sse(np.column_stack([one,site_g]),site_ph,np.ones(60))
    collinearity.update({'between_site_ph_R2_geography':1-site_ph_g/site_ph_n,
                         'between_site_ph_VIF':site_ph_n/site_ph_g})
    # Within-site pH sensitivity targets deviations from site and campaign-position means.
    site_dummies = pd.get_dummies(cohort.site.astype(str), drop_first=True, dtype=float).to_numpy()
    n_site = np.column_stack([x[:, columns['N']], site_dummies])
    sse_site = direct_sse(n_site, y, weights)
    sse_site_ph = direct_sse(np.column_stack([n_site, x[:, k]]), y, weights)
    within = {'partial_R2_within_site': (sse_site-sse_site_ph)/sse_site,
              'common_primary_denominator_fraction': (sse_site-sse_site_ph)/primary['sse_N'],
              'nuisance_rank': int(np.linalg.matrix_rank(n_site)),
              'full_rank': int(np.linalg.matrix_rank(np.column_stack([n_site, x[:, k]])))}
    # Cases bootstrap: no row permutations, no hypothesis-test p-values.
    order = np.argsort(cohort.groupby('site').transect_km.first().loc[stats.sites].to_numpy())
    draws, intervals, resampling_design = [], [], []
    for length in [1, 5, 10, 15]:
        rng = np.random.default_rng(SEED + length)
        method = 'site_cluster' if length == 1 else f'noncircular_block_{length}'
        for replicate in range(replicates):
            multiplicity = draw_sites(rng, order, length)
            values, rank = stats.evaluate(multiplicity)
            rank_ok = identifiable(rank)
            draws.append({'method': method, 'replicate': replicate, 'unique_sites': int((multiplicity > 0).sum()),
                          'effects_identifiable': rank_ok, 'nuisance_rank': rank['N'], **values})
        selected = [row for row in draws if row['method'] == method]
        valid = [row for row in selected if row['effects_identifiable']]
        assert len(valid) >= .95*replicates
        for fraction in FRACTIONS:
            lo, median, hi = np.quantile([row[fraction] for row in valid], [.025,.5,.975])
            intervals.append({'method': method, 'fraction': fraction, 'point': primary[fraction],
                              'percentile_2_5': lo, 'bootstrap_median': median, 'percentile_97_5': hi,
                              'basic_2_5':2*primary[fraction]-hi, 'basic_97_5':2*primary[fraction]-lo,
                              'valid_draws': len(valid), 'unidentifiable_effect_draws': replicates-len(valid), 'draws_missing_nuisance_categories':sum(row['nuisance_rank']<expected['N'] for row in selected)})
        ordered_dist = cohort.groupby('site').transect_km.first().loc[stats.sites].to_numpy()[order]
        spans = ordered_dist[length-1:] - ordered_dist[:len(order)-length+1]
        resampling_design.append({'method':method,'sites_per_block':length,'median_block_span_km':float(np.median(spans)),
                                  'minimum_block_span_km':float(spans.min()),'maximum_block_span_km':float(spans.max())})
    between_draws=[]; between_intervals=[]
    for length in [1,5,10,15]:
        rng=np.random.default_rng(SEED+length)
        method='site_cluster' if length==1 else f'noncircular_block_{length}'
        values=[]
        for replicate in range(replicates):
            multiplicity=draw_sites(rng,order,length)
            value,rank=between_stats.evaluate(multiplicity)
            assert identifiable(rank)
            values.append(value)
            between_draws.append({'method':method,'replicate':replicate,**value})
        for fraction in FRACTIONS:
            lo,median,hi=np.quantile([v[fraction] for v in values],[.025,.5,.975])
            between_intervals.append({'method':method,'fraction':fraction,'point':between[fraction],
                                     'percentile_2_5':lo,'bootstrap_median':median,'percentile_97_5':hi,
                                     'basic_2_5':2*between[fraction]-hi,'basic_97_5':2*between[fraction]-lo})
    pd.DataFrame(between_draws).to_csv(output/'between_site_bootstrap_draws.tsv',sep='\t',index=False,float_format='%.12g')
    pd.DataFrame(between_intervals).to_csv(output/'between_site_bootstrap_intervals.tsv',sep='\t',index=False,float_format='%.12g')
    pd.DataFrame(draws).to_csv(output/'bootstrap_draws.tsv',sep='\t',index=False,float_format='%.12g')
    pd.DataFrame(intervals).to_csv(output/'bootstrap_intervals.tsv',sep='\t',index=False,float_format='%.12g')
    influence=[]
    for i, site in enumerate(stats.sites):
        multiplicity = np.ones(60); multiplicity[i] = 0
        value, rank = stats.evaluate(multiplicity)
        between_value,between_rank=between_stats.evaluate(multiplicity)
        assert identifiable(between_rank)
        influence.append({'excluded_site':int(site), 'full_rank':rank==expected, **value,
                          **{'between_site_'+key:val for key,val in between_value.items()}})
    pd.DataFrame(influence).to_csv(output/'leave_site_out.tsv',sep='\t',index=False,float_format='%.12g')
    # Maximum pH group removal, holding all selected features fixed; reweight surviving site rows.
    maximum = int(cohort.ph.argmax())
    keep = np.arange(len(cohort)) != maximum
    subset = cohort.loc[keep]
    subset_weights = 1/subset.groupby('site').site.transform('size').to_numpy()
    maximum_partition = partition_from_sse({key:direct_sse(x[keep][:,cols], y[keep], subset_weights) for key,cols in columns.items()})
    summary = {'seed':SEED,'bootstrap_replicates_per_method':replicates,'groups':len(cohort),'sites':60,'genera':len(taxa),
               'selected_profile_columns':profile_count,'pseudocount':.5,'minimum_group_reads':2000,
               'primary':primary,'design_ranks':ranks,'collinearity':collinearity,'within_site_sensitivity':within,
               'maximum_ph_group':cohort.iloc[maximum][['trip','site','compartment','ph','n_ph_specimens']].to_dict(),
               'maximum_ph_group_removed_partition':maximum_partition,'resampling_design':resampling_design,
               'between_site':between, 'between_site_design_ranks':between_ranks,
               'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
               'fractions':'raw descriptive; common nuisance-residual weighted SS denominator',
               'linkage':'campaign-site-compartment; physical specimen pairing unverified',
               'resampling':'percentile uncertainty under site or contiguous-block sampling assumptions; no hypothesis-test p-values'}
    (output/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    output_files = sorted(output.glob('*'))
    (output/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in output_files))
    print(json.dumps({'primary':primary,'collinearity':collinearity,'within_site':within},indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replicates',type=int,default=1999)
    args=parser.parse_args()
    if args.replicates < 20:
        parser.error('At least 20 bootstrap draws required')
    run(args.output,args.replicates)
