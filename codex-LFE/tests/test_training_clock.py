import numpy as np
import pandas as pd

from codex_lfe.training_clock import peak_offset, accepted_site_clocks


def test_histogram_recovers_offset_above_dense_background():
    rng = np.random.default_rng(12)
    ref = np.arange(100.,10100.,100.)
    picks = np.r_[ref+12.4+rng.normal(0,.015,len(ref)),rng.uniform(0,11000,10000)]
    got = peak_offset(ref,picks)
    assert got['passed_peak_check']
    assert abs(got['offset_s']-12.4)<=.25


def test_no_picks_cannot_make_clock():
    assert not peak_offset(np.arange(10),[])['passed_peak_check']


def test_colocated_clocks_do_not_make_independent_confirmation():
    rows = pd.DataFrame({'site_id':['PGC_B011','PGC_B011'],
                         'offset_s':[12.,12.25],'accepted_station':[True,True]})
    assert accepted_site_clocks(rows)[0] is None
    rows.loc[2] = ['SNB',12.5,True]
    assert accepted_site_clocks(rows)[0] is not None
    rows.loc[2,'offset_s'] = 20.
    assert accepted_site_clocks(rows)[0] is None
