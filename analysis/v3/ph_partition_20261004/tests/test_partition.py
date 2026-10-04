import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from analyze import SiteAverageStatistics, SiteStatistics, partition_from_sse, direct_sse, draw_sites, load_cohort, design, verify_inputs, identifiable

class PartitionTests(unittest.TestCase):
    def test_known_partition_and_signed_shared(self):
        values=partition_from_sse({'N':100., 'P':70., 'G':60., 'PG':40.})
        self.assertAlmostEqual(values['unique_ph'],.2)
        self.assertAlmostEqual(values['unique_geography'],.3)
        self.assertAlmostEqual(values['shared'],.1)
        self.assertAlmostEqual(values['unexplained'],.4)
        self.assertLess(partition_from_sse({'N':100.,'P':90.,'G':90.,'PG':50.})['shared'],0)

    def test_site_draw_multiplicity_matches_explicit_repeated_rows(self):
        rng=np.random.default_rng(3)
        sites=np.repeat(np.arange(8),[2,3,2,4,2,3,2,4])
        x=np.column_stack([np.ones(len(sites)),rng.normal(size=(len(sites),3))])
        y=rng.normal(size=(len(sites),5))
        w=1/np.bincount(sites)[sites]
        cols={'N':np.array([0]),'P':np.array([0,1]),'G':np.array([0,2,3]),'PG':np.arange(4)}
        stats=SiteStatistics(x,y,sites,w,cols)
        draws=np.array([0,0,2,3,3,3,6,7])
        counts=np.bincount(draws,minlength=8)
        actual,_=stats.evaluate(counts)
        indices=np.concatenate([np.flatnonzero(sites==site) for site in draws])
        expected=partition_from_sse({k:direct_sse(x[indices][:,c],y[indices],w[indices]) for k,c in cols.items()})
        for k in actual:
            self.assertAlmostEqual(actual[k],expected[k],places=10)

    def test_between_site_refits_common_nuisance(self):
        rng=np.random.default_rng(41)
        sites=np.repeat(np.arange(10),3)
        x=np.column_stack([np.ones(30),np.tile([0,1,0],10),rng.normal(size=(30,3))])
        y=rng.normal(size=(30,6)); w=np.ones(30)/3
        columns={'N':np.arange(2),'P':np.arange(3),'G':np.array([0,1,3,4]),'PG':np.arange(5)}
        draws=np.array([0,0,2,3,3,5,6,7,8,9]); m=np.bincount(draws,minlength=10)
        value,rank=SiteAverageStatistics(x,y,sites,w,columns).evaluate(m)
        indices=np.concatenate([np.flatnonzero(sites==site) for site in draws])
        n=x[indices,:2]; joint=np.column_stack([y,x[:,2:]])[indices]
        residual=joint-n@np.linalg.lstsq(n,joint,rcond=None)[0]
        means=residual.reshape(10,3,-1).mean(axis=1)
        sx=np.column_stack([np.ones(10),means[:,6:]])
        cols={'N':[0],'P':[0,1],'G':[0,2,3],'PG':[0,1,2,3]}
        expected=partition_from_sse({k:direct_sse(sx[:,c],means[:,:6],np.ones(10)) for k,c in cols.items()})
        for key in value: self.assertAlmostEqual(value[key],expected[key],places=10)
        self.assertTrue(identifiable(rank))

    def test_nuisance_signal_and_global_weight_invariance(self):
        rng=np.random.default_rng(19)
        n=60
        x=np.column_stack([np.ones(n),rng.normal(size=(n,4))]); y=rng.normal(size=(n,5)); w=np.ones(n)
        cols={'N':np.array([0,1]),'P':np.array([0,1,2]),'G':np.array([0,1,3,4]),'PG':np.arange(5)}
        a=partition_from_sse({k:direct_sse(x[:,c],y,w) for k,c in cols.items()})
        shifted=y+x[:,:2]@rng.normal(size=(2,5))
        b=partition_from_sse({k:direct_sse(x[:,c],shifted,w*3) for k,c in cols.items()})
        for key in ['unique_ph','unique_geography','shared','unexplained']:
            self.assertAlmostEqual(a[key],b[key],places=11)

    def test_non_circular_blocks_and_determinism(self):
        class Stub:
            def integers(self, high, size):
                self.high=high
                return np.full(size,high-1)
        stub=Stub(); counts=draw_sites(stub,np.arange(60),10)
        self.assertEqual(stub.high,51)
        self.assertEqual(counts.sum(),60)
        np.testing.assert_equal(counts[:50],0)
        np.testing.assert_equal(counts[50:],6)
        np.testing.assert_equal(draw_sites(np.random.default_rng(12),np.arange(60),5),draw_sites(np.random.default_rng(12),np.arange(60),5))

    def test_frozen_cohort_weights_and_design(self):
        verify_inputs(); cohort,y,taxa,reads,profiles=load_cohort()
        self.assertEqual((len(cohort),cohort.site.nunique(),len(taxa),profiles),(560,60,200,1089))
        self.assertFalse(cohort[['trip','site','compartment']].duplicated().any())
        self.assertTrue((reads>=2000).all())
        np.testing.assert_allclose(y.sum(axis=1),0,atol=1e-11)
        w=1/cohort.groupby('site').site.transform('size').to_numpy()
        np.testing.assert_allclose(np.bincount(cohort.site,weights=w)[1:],1)
        x,cols=design(cohort)
        value,rank=SiteStatistics(x,y,cohort.site.to_numpy(),w,cols).evaluate()
        self.assertTrue(identifiable(rank))
        self.assertAlmostEqual(sum(value[k] for k in ['unique_ph','unique_geography','shared','unexplained']),1)
        self.assertFalse(identifiable({'N':12,'P':13,'G':14,'PG':14}))

if __name__=='__main__': unittest.main()
