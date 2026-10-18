import numpy as np
import pytest

from codex_lfe.matched_filter import extract_candidates


def site_series():
    rng=np.random.default_rng(82)
    sites={s:rng.normal(0,.03,1000) for s in ['A','B','C','D']}
    for s in sites: sites[s][500]+=.8
    sites['A'][700]=10  # must not bypass four-site support
    sites['D'][20]=np.nan
    return sites


def test_equal_site_weights_preserve_candidates_and_masks_exactly():
    sites=site_series()
    args=dict(fs=20,start_s=0,threshold=5,min_sites=4)
    a,ma=extract_candidates(sites,**args)
    b,mb=extract_candidates(sites,**args,site_weights={s:7 for s in sites})
    assert a==b and np.array_equal(ma,mb)
    assert len(a)==1 and a[0]['time_s']==25


def test_weight_scale_invariance_and_single_site_artifact_rejection():
    sites=site_series();weights={'A':3.,'B':1.,'C':.7,'D':.4}
    args=dict(fs=20,start_s=0,threshold=5,min_sites=4)
    a,ma=extract_candidates(sites,**args,site_weights=weights)
    b,mb=extract_candidates(sites,**args,site_weights={s:2*w for s,w in weights.items()})
    assert np.array_equal(ma,mb) and len(a)==len(b)==1
    assert a[0]['time_s']==b[0]['time_s']==25
    assert a[0]['score_mad']==pytest.approx(b[0]['score_mad'])


@pytest.mark.parametrize('weights',[{'A':1}, {'A':1,'B':1,'C':1,'D':0},
                                  {'A':1,'B':1,'C':1,'D':float('nan')}])
def test_reject_invalid_site_weights(weights):
    with pytest.raises(ValueError):
        extract_candidates(site_series(),20,0,5,site_weights=weights)
