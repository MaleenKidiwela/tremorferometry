import numpy as np
import pandas as pd
import pytest

from codex_lfe.reuse_association import prepare_picks, shift_picks, finish_events, site_groups


def test_network_identity_boundaries_and_duplicate_rejection():
    f = pd.DataFrame({'net': ['CN', 'XX', 'CN'], 'sta': ['PGC']*3,
                      'phase': ['S']*3, 'prob': [.8]*3, 'thr_first': [.1]*3,
                      'time': pd.to_datetime([0, 10, 20], unit='s', utc=True)})
    picks, removed = prepare_picks(f, {'CN.PGC'}, .3, 0, 86400)
    assert removed == 1
    assert picks.source_pick_id.tolist() == [2]
    assert picks.station.tolist() == ['CN.PGC']
    with pytest.raises(ValueError, match='duplicate'):
        prepare_picks(pd.concat([f, f.iloc[[2]]]), {'CN.PGC'}, .3, 0, 86400)


@pytest.mark.parametrize('mode', ['station', 'station_hour'])
def test_controls_preserve_identity_counts_and_common_phase_shift(mode):
    p = pd.DataFrame({'station': ['CN.PGC']*4 + ['PB.B011']*2,
                      'phase': ['P', 'S']*3, 'time': [100., 110., 3700., 3710., 100., 110.],
                      'source_pick_id': np.arange(6)})
    shifted, _ = shift_picks(p, mode, 7, 0, 7200)
    shifted = shifted.sort_values('source_pick_id')
    assert set(shifted.source_pick_id) == set(p.source_pick_id)
    delta = shifted.time.to_numpy() - p.time.to_numpy()
    assert delta[0] == pytest.approx(delta[1])
    assert delta[4] == pytest.approx(delta[5])
    assert np.all((np.abs(delta) >= 300) & (np.abs(delta) <= 3000))
    assert delta[0] != delta[4]


def test_distinct_s_stations_and_edge_filter():
    e = pd.DataFrame({'idx': [0, 1, 2], 'time': [5000., 5000., 1.], 'depth': [30.]*3})
    a = pd.DataFrame({'event_idx': [0]*4 + [1]*4 + [2]*4,
                     'pick_idx': range(12), 'station': ['A', 'A', 'B', 'C'] + ['A','B','C','D']*2,
                     'phase': ['S']*12, 'residual': [.1]*12})
    events, assignments = finish_events(e, a, 4, 3000, 9000)
    assert events.idx.tolist() == [1]
    assert set(assignments.event_idx) == {1}


def test_colocated_sensors_share_control_shifts():
    stations = pd.DataFrame({'id': ['CN.PGC', 'PB.B011', 'CN.SNB'],
                             'latitude': [48.6498, 48.649544, 48.7751],
                             'longitude': [-123.4521, -123.448189, -123.1723]})
    groups = site_groups(stations, 1.)
    assert groups['CN.PGC'] == groups['PB.B011'] != groups['CN.SNB']
    p = pd.DataFrame({'station': stations.id, 'time': [100.]*3})
    for mode in ['station', 'station_hour']:
        shifted, _ = shift_picks(p, mode, 9, 0, 7200, groups)
        indexed = shifted.set_index('station')
        assert indexed.loc['CN.PGC', 'time'] == indexed.loc['PB.B011', 'time']
        assert indexed.loc['CN.PGC', 'time'] != indexed.loc['CN.SNB', 'time']
