import numpy as np
import pytest
from codex_lfe.resampling import fractional_samples


@pytest.mark.parametrize('frequency', [1, 2, 5, 8])
@pytest.mark.parametrize('fraction', [0, 0.15, 0.3, 0.5, 0.9])
def test_fractional_interpolation_preserves_full_lfe_passband(frequency, fraction):
    fs = 20
    t = np.arange(5000) / fs
    x = np.sin(2 * np.pi * frequency * t)
    actual = fractional_samples(x, 200 + fraction, 4000)
    truth = np.sin(2 * np.pi * frequency * (np.arange(4000) + 200 + fraction) / fs)
    np.testing.assert_allclose(actual, truth, atol=3e-4)


def test_kernel_does_not_cross_gap_boundary():
    with pytest.raises(ValueError, match='beyond valid run'):
        fractional_samples(np.ones(100), 10.5, 50)
