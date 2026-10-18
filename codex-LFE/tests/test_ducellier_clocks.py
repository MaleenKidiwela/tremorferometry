import numpy as np
from codex_lfe.ducellier_template_clocks import lag_cc

def test_template_clock_sign_with_delayed_waveform():
    rng=np.random.default_rng(12)
    a=rng.normal(size=(3,800));b=np.pad(a,((0,0),(20,0)))
    lags=np.arange(-40,41);cc=np.array([lag_cc(a,b,int(lag)) for lag in lags])
    lag=lags[np.nanargmax(cc)]
    assert lag==20 and abs(cc.max()-1)<1e-12
    # B matches one second earlier; adding +lag/fs restores A's clock.
    assert 99+lag/20==100
