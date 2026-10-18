"""Band-limited fractional sampling for UTC alignment after rate conversion."""
from __future__ import annotations

import numpy as np
from scipy.signal import fftconvolve


def fractional_samples(signal, first_index, count, half_width=32):
    """Sample x[first_index + arange(count)] with a Kaiser-windowed sinc.

    The interpolation preserves the 1-8 Hz passband at 20 Hz. Linear interpolation
    near the upper passband can attenuate an off-grid channel by roughly a factor
    of two. Caller must exclude >= half_width samples at run edges, where the
    finite kernel would otherwise reach outside valid data.
    """
    x = np.asarray(signal, float)
    integer = int(np.floor(first_index))
    fraction = float(first_index - integer)
    if fraction > 1 - 1e-8:
        integer += 1
        fraction = 0.0
    if fraction < 1e-8:
        fraction = 0.0
    if integer < half_width or integer + count > len(x) - half_width:
        raise ValueError('fractional sampling kernel extends beyond valid run')
    if fraction == 0:
        return x[integer:integer + count].copy()
    m = np.arange(-half_width, half_width + 1)
    kernel = np.sinc(m + fraction) * np.kaiser(len(m), 8.6)
    kernel /= kernel.sum()
    shifted = fftconvolve(x, kernel, mode='same')
    return shifted[integer:integer + count]
