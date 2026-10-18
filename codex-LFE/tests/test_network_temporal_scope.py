import pandas as pd
from codex_lfe.network_temporal_scope import exposure_summary


def test_colocated_handoff_and_instrument_boundary_do_not_double_count():
    lo = pd.Timestamp('2011-01-01', tz='UTC')
    mid = lo + pd.Timedelta(days=1)
    hi = mid + pd.Timedelta(days=1)
    frame = pd.DataFrame([
        dict(id='A', site_id='S', lo=lo, hi=mid),
        dict(id='A', site_id='S', lo=mid, hi=hi),
        dict(id='B', site_id='S', lo=mid, hi=hi),
        dict(id='C', site_id='T', lo=lo, hi=mid),
    ])
    got = exposure_summary(frame, lo, hi)
    assert got['mean_stations'] == 2
    assert got['mean_sites'] == 1.5
    assert (got['min_sites'], got['max_sites']) == (1, 2)
    assert exposure_summary(frame, mid, hi)['mean_sites'] == 1


def test_empty_region_and_partial_month_gap_have_real_zero_exposure():
    lo = pd.Timestamp('2019-01-01', tz='UTC')
    hi = lo + pd.Timedelta(days=4)
    frame = pd.DataFrame([dict(id='A', site_id='S', lo=lo, hi=lo+pd.Timedelta(days=1))])
    got = exposure_summary(frame, lo, hi)
    assert got['mean_sites'] == .25
    assert got['min_sites'] == 0
    assert exposure_summary(frame.iloc[:0], lo, hi)['mean_sites'] == 0
