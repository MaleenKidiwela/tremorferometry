import numpy as np
from codex_lfe.central_cross_station import pair_scores

def test_pair_lag_sign_and_amplitude_invariance():
    rng=np.random.default_rng(812);a=rng.normal(size=(2,2,160));b=rng.normal(size=(2,2,240))
    b[0,1,50:210]=-7*a[0,0]+3
    b[1,0,30:190]=.1*a[1,1]-2
    scores,lags,components=pair_scores(a,b)
    np.testing.assert_allclose(scores,[1,1],atol=1e-10)
    np.testing.assert_allclose(lags,[.5,-.5],atol=1e-10)
    np.testing.assert_array_equal(components,[1,2])
