import numpy as np
import pandas as pd

from codex_lfe.miss_diagnosis import matcher_window, phase_window, support_counts


def test_phase_window_rejects_gap_and_day_wrap_and_keeps_peak_distinct():
    x=np.zeros(200);x[100]=.8
    assert phase_window(x,0,10,.1,1,np.array([]),np.array([]))[0] is False
    covered,maximum,peak,residual=phase_window(x,0,10,10,1,np.array([10.5]),np.array([.4]))
    assert covered and maximum==.8 and peak==.4 and residual==.5
    x[119]=np.nan  # outside arrival window but inside peak guard
    assert phase_window(x,0,10,10,1,np.array([]),np.array([]))[0] is False


def test_s_support_does_not_count_p_or_colocated_stations_twice():
    frame=pd.DataFrame({'station':['PGC','B011','SNB','VGZ','B926'],
        'site':['PGC','PGC','SNB','VGZ','B926'],'phase':['S','S','S','S','P'],
        'covered':[True]*5,'max_pick_score':[.5]*5,'max_probability':[.6]*5})
    counts=support_counts(frame,.3)
    assert counts['n_s_stations']==4 and counts['n_s_sites']==3


def test_matcher_joint_support_must_occur_at_same_sample():
    good=np.ones(5,bool); z=np.array([0,10,8,0,0]); support=np.array([0,3,4,0,0])
    row=matcher_window(z,support,good,2,2,4,9)
    assert row['status']=='network_threshold_passes_but_site_support_fails'
    assert row['max_supported_mad']==8
    support[1]=4
    assert matcher_window(z,support,good,2,2,4,9)['status']=='joint_threshold_and_site_support_pass'
