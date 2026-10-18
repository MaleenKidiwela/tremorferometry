import pandas as pd
import pytest
from codex_lfe.month_units import plan_units

def epochs():
    return pd.DataFrame([dict(network='XX',station='A',channels='XX.A..HHZ',start_utc='2026-08-01T00:00:00Z',end_utc_exclusive='2026-09-02T12:00:00Z'),dict(network='XX',station='A',channels='XX.A..BHZ',start_utc='2026-09-02T12:00:00Z',end_utc_exclusive='2026-09-03T10:00:00Z')])

def test_descending_days_partition_core_and_clip_response_epochs():
    units,requests=plan_units(epochs(),'2026-09','2026-09-03T10:00:00Z')
    assert units.unit_id.tolist()==['20260903','20260902','20260901']
    spans=pd.to_datetime(requests.core_end_utc_exclusive)-pd.to_datetime(requests.core_start_utc)
    assert spans.sum().total_seconds()==58*3600
    assert len(requests[requests.unit_id=='20260902'])==2
    boundary=requests[(requests.unit_id=='20260902')&(requests.channels=='XX.A..HHZ')].iloc[0]
    assert boundary.right_padding_s==0
    assert requests.iloc[-1].request_start_utc=='2026-08-31T23:58:00+00:00'
    assert requests.iloc[0].right_padding_s==0

def test_overlapping_metadata_rejected():
    x=epochs();x.loc[1,'start_utc']='2026-09-02T11:00:00Z'
    with pytest.raises(ValueError,match='Overlapping'):plan_units(x,'2026-09','2026-09-03T10:00:00Z')

def test_partial_day_metadata_gaps_remain_explicit():
    from codex_lfe.month_units import metadata_gaps
    x=epochs().iloc[:1].copy();x.loc[0,'start_utc']='2026-09-02T06:00:00Z'
    units,requests=plan_units(x,'2026-09','2026-09-03T10:00:00Z')
    gaps=metadata_gaps(units,requests,{'XX.A','XX.B'})
    a=gaps[gaps.station=='XX.A']
    missing=(pd.to_datetime(a.core_end_utc_exclusive)-pd.to_datetime(a.core_start_utc)).sum().total_seconds()
    assert missing==52*3600
    assert len(a[a.unit_id=='20260902'])==2
    assert len(gaps[gaps.station=='XX.B'])==3

def test_month_without_epochs_keeps_request_schema_and_full_gaps():
    from codex_lfe.month_units import metadata_gaps
    units,requests=plan_units(epochs(),'2010-01','2026-09-03T10:00:00Z')
    assert requests.empty and {'network','station','unit_id'}<=set(requests.columns)
    gaps=metadata_gaps(units,requests,{'XX.A','XX.B'})
    assert len(gaps)==62
    assert (pd.to_datetime(gaps.core_end_utc_exclusive)-pd.to_datetime(gaps.core_start_utc)).sum().total_seconds()==62*86400
