import numpy as np
from codex_lfe.heldout import sample_series


def test_heldout_sampling_respects_time_and_missing_boundaries():
    series = np.arange(100, dtype=float)
    series[40] = np.nan
    times = [1000, 1000.125, 1002, 999, 1005]
    result = sample_series(series, times, 1000, 20)
    np.testing.assert_allclose(result[:2], [0, 2.5])
    assert np.isnan(result[2:]).all()
