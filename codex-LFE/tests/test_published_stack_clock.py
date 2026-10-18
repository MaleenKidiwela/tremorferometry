import numpy as np
import pytest
from codex_lfe.published_stack_clock import best_s_match, phase_site_clocks
import pandas as pd


def test_transferred_phase_time_includes_both_window_anchors():
    rng=np.random.default_rng(19)
    template=rng.normal(size=(3,160))
    stack=rng.normal(scale=.01,size=(3,1400))
    stack[:,640:800]+=template
    arrival,cc,second=best_s_match(stack,template)
    assert arrival==pytest.approx(24.)
    assert cc>.99
    assert cc-second>.5


def test_two_phases_are_not_two_sites_and_disagreement_is_rejected():
    rows=pd.DataFrame({'station':['CN.PGC','CN.PGC','CN.SNB','CN.SNB'],
                       'origin_minus_published_s':[3.,3.5,3.,6.],
                       'accepted_waveform_transfer':[True]*4})
    clocks=phase_site_clocks(rows)
    assert clocks.to_dict()=={'CN.PGC':3.25}

def test_interval_identity_handles_decimal_boundary_without_nearest_guess():
    import pandas as pd
    from types import SimpleNamespace
    from codex_lfe.published_stack_clock import interval_source_identity
    s=SimpleNamespace(latitude=48.43,longitude=-123.65,depth_km=34.32)
    f=pd.DataFrame([dict(source_latitude_deg=48.425667,source_longitude_deg=-123.655,source_depth_km=34.32)],index=[47])
    assert interval_source_identity(s,f)==47
    s.depth_km=34.34
    assert interval_source_identity(s,f) is None
    s.depth_km=34.32
    assert interval_source_identity(s,pd.concat([f,f.rename(index={47:48})])) is None
