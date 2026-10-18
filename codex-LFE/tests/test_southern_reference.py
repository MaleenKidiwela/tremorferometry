import pandas as pd
import pytest
from codex_lfe.southern_reference import read_release


@pytest.fixture
def release(tmp_path):
    metadata=tmp_path/'metadata.xml'
    locations='\n'.join(f'{i} 40.065 -123.620 28.1' for i in range(1,28))
    metadata.write_text('<metadata><abstract>\n'+locations+'\n</abstract></metadata>')
    catalog=tmp_path/'catalog.csv'
    pd.DataFrame([{'%year':2024,'month':1,'day':1,'hour':23,'minute':59,'second':59.99,
                   's_of_day':86399.99,'ID':'example_named_family','idnum':1}]).to_csv(catalog,index=False)
    return catalog,metadata


def test_named_family_and_subsecond_clock_preserved(release):
    events,_=read_release(*release)
    assert events.iloc[0].family_id=='shelly2025_example_named_family'
    assert events.iloc[0].time_s==pytest.approx(pd.Timestamp('2024-01-01T23:59:59.99Z').timestamp(),abs=1e-6,rel=0)
    assert events.iloc[0].time_basis=='shelly2025_template_start_first_station_estimated_S_minus_1s'
    assert events.iloc[0].latitude==40.065


def test_disagreeing_clock_rejected(release):
    catalog,_=release
    frame=pd.read_csv(catalog);frame['s_of_day']-=1;frame.to_csv(catalog,index=False)
    with pytest.raises(ValueError,match='clocks disagree'):read_release(*release)


def test_duplicate_source_row_rejected(release):
    catalog,_=release
    frame=pd.read_csv(catalog);pd.concat([frame,frame]).to_csv(catalog,index=False)
    with pytest.raises(ValueError,match='Duplicate'):read_release(*release)


def test_unmapped_family_rejected(release):
    catalog,_=release
    frame=pd.read_csv(catalog);frame['idnum']=28;frame.to_csv(catalog,index=False)
    with pytest.raises(ValueError,match='Unmapped'):read_release(*release)
