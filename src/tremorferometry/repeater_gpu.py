"""GPU port of the all-pairs max-shifted cross-correlation used in family discovery.

The CPU version (`repeater.all_pairs_cc_max_shifted`) loops over rows in Python and
does one batched irfft per row -- ~7.9 s for a 2000x80 bin.  This collapses that into
a few large batched FFT ops on the GPU.  Algorithm is identical (FFT cross-correlation
of L2-normalized rows, max over lags in [-s, +s], zero diagonal, symmetrize), so the
output matches the CPU version to FFT roundoff (~1e-6), well below the CC>=0.8 cluster
threshold.

cupy is imported lazily so importing this module never touches the GPU.
"""
from __future__ import annotations

import numpy as np


def all_pairs_cc_max_shifted_gpu(
    X: np.ndarray, max_shift_samples: int = 20, block: int = 256,
) -> np.ndarray:
    """GPU all-pairs max-shifted CC on L2-normalized rows of X.

    Returns symmetric (n, n) float32 matrix of max CC over lags in [-s, +s],
    diagonal zeroed -- identical contract to repeater.all_pairs_cc_max_shifted.
    """
    import cupy as cp

    n, m = X.shape
    pad = 1 << int(np.ceil(np.log2(2 * m)))
    s = int(max_shift_samples)

    Xg = cp.asarray(X, dtype=cp.float32)
    F = cp.fft.rfft(Xg, n=pad, axis=1)          # (n, nf) complex64
    out = cp.zeros((n, n), dtype=cp.float32)

    for i0 in range(0, n, block):
        i1 = min(i0 + block, n)
        # prod[b, j, :] = conj(F[i0+b]) * F[j]   -> (b, n, nf)
        prod = cp.conj(F[i0:i1])[:, None, :] * F[None, :, :]
        cc_full = cp.fft.irfft(prod, n=pad, axis=2)   # (b, n, pad)
        cc_pos = cc_full[:, :, : s + 1].max(axis=2)
        cc_neg = cc_full[:, :, -s:].max(axis=2)
        out[i0:i1] = cp.maximum(cc_pos, cc_neg)
        del prod, cc_full

    cp.fill_diagonal(out, 0.0)
    out = cp.maximum(out, out.T)
    return cp.asnumpy(out)
