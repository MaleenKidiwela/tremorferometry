import numpy as np
import pytest

from codex_lfe.neural import orientation_matrix, picks_from_probabilities


def test_sensor_rotation_uses_metadata_not_channel_suffix():
    matrix = orientation_matrix([{'azimuth': 0, 'dip': -90}, {'azimuth': 35, 'dip': 0}, {'azimuth': 125, 'dip': 0}])
    geo = np.random.default_rng(1).normal(size=(3, 100))
    sensor = matrix @ geo
    np.testing.assert_allclose(np.linalg.solve(matrix, sensor), geo, atol=1e-12)
    assert not np.allclose(sensor[1], geo[1])


def test_parallel_sensor_axes_rejected():
    with pytest.raises(ValueError, match='ill-conditioned'):
        orientation_matrix([{'azimuth': 0, 'dip': 0}] * 3)


def test_missing_interval_edges_do_not_generate_picks():
    x = np.zeros((2, 300), np.float32)
    x[:, :30] = np.nan
    x[1, 30] = 1
    x[1, 100] = 0.8
    picks = picks_from_probabilities(x, 1000, 20, 'XX.TEST')
    assert len(picks) == 1
    assert picks.time_s.iloc[0] == 1005
    assert picks.phase.iloc[0] == 'S'


def test_prediction_coverage_excludes_picker_edge_guard_and_gaps():
    from codex_lfe.neural import prediction_coverage
    values=np.ones((2,200));values[:,100:120]=np.nan
    result=prediction_coverage(values,1000.,10.,'CN.PGC')
    p=result[result.phase=='P']
    assert p.start_s.tolist()==[1001.,1013.]
    assert p.end_s.tolist()==[1009.,1019.]
