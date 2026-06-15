"""LFE waveform feature definitions (Stage 1, hand-crafted).

Each window is a 1-D numpy array already sliced to [t_d - PRE, t_d + POST].
Features are physically motivated to separate emergent low-frequency LFEs from
impulsive / cultural / EQ / glitch contaminants. Vertical features always;
polarization features when 3-component windows are supplied.

compute_all(z, fs, t_pre, h1=None, h2=None) -> dict
"""
import numpy as np
from scipy import signal


def _detrend(x):
    x = np.asarray(x, float)
    if x.size == 0 or not np.all(np.isfinite(x)):
        x = np.nan_to_num(x)
    return signal.detrend(x) if x.size > 1 else x


def _psd(x, fs):
    x = _detrend(x)
    nper = int(min(len(x), max(64, fs * 4)))
    f, P = signal.welch(x, fs=fs, nperseg=nper)
    return f, P + 1e-20


def _bandpass(x, fs, lo, hi):
    hi = min(hi, fs / 2.0 - 0.5)
    sos = signal.butter(4, [lo, hi], btype="band", fs=fs, output="sos")
    return signal.sosfiltfilt(sos, _detrend(x))


def spectral_features(x, fs):
    f, P = _psd(x, fs)
    tot = P.sum()
    Pn = P / tot
    centroid = float((f * Pn).sum())
    spread = float(np.sqrt(((f - centroid) ** 2 * Pn).sum()))
    flatness = float(np.exp(np.mean(np.log(P))) / np.mean(P))
    c = np.cumsum(P)
    rolloff = float(f[np.searchsorted(c, 0.85 * c[-1])])
    peak = float(f[np.argmax(P)])

    def bf(lo, hi):
        m = (f >= lo) & (f < hi)
        return float(P[m].sum() / tot)

    nyq = fs / 2.0
    feats = dict(
        sp_centroid=centroid, sp_spread=spread, sp_flatness=flatness,
        sp_rolloff85=rolloff, sp_peakfreq=peak,
        bf_1_2=bf(1, 2), bf_2_4=bf(2, 4), bf_4_8=bf(4, 8),
        bf_8_16=bf(8, 16), bf_16up=bf(16, nyq),
        hf_ratio=float(P[f >= 8].sum() / (P[(f >= 2) & (f < 8)].sum() + 1e-20)),
        sp_entropy=float(-(Pn * np.log(Pn)).sum() / np.log(len(Pn))),
    )
    m = (f >= 2) & (f < min(40, nyq))
    feats["sp_slope"] = float(np.polyfit(np.log(f[m] + 1e-9), np.log(P[m]), 1)[0]) if m.sum() > 3 else np.nan
    return feats


def envelope_features(x, fs, t_pre):
    xb = _bandpass(x, fs, 2, 8)
    env = np.abs(signal.hilbert(xb))
    pk = float(env.max()) + 1e-20
    rms = float(np.sqrt(np.mean(xb ** 2))) + 1e-20
    e = env / pk
    em = e.mean(); es = e.std() + 1e-20
    ipk = int(env.argmax())
    pre = env[: ipk + 1]
    try:
        i10 = np.where(pre >= 0.1 * pk)[0][0]
        i90 = np.where(pre >= 0.9 * pk)[0][0]
        rise = float((i90 - i10) / fs)
    except Exception:
        rise = np.nan
    n_pre = int(max(0.5, (t_pre - 5)) * fs)
    n_pre = min(n_pre, len(xb) // 3) if n_pre > 0 else len(xb) // 3
    noise = xb[:n_pre]; sig = xb[n_pre:]
    snr = (np.sqrt(np.mean(sig ** 2)) + 1e-20) / (np.sqrt(np.mean(noise ** 2)) + 1e-20)
    return dict(
        env_crest=float(pk / rms),
        env_kurt=float(((e - em) ** 4).mean() / es ** 4),
        env_skew=float(((e - em) ** 3).mean() / es ** 3),
        rise_time=rise,
        duration=float((env >= 0.5 * pk).sum() / fs),
        zcr=float(((xb[:-1] * xb[1:]) < 0).mean() * fs),
        snr=float(snr),
    )


def polarization_features(z, h1, h2, fs):
    Z = _bandpass(z, fs, 2, 8); H1 = _bandpass(h1, fs, 2, 8); H2 = _bandpass(h2, fs, 2, 8)
    ez = np.mean(Z ** 2) + 1e-20
    eh = np.mean(H1 ** 2) + np.mean(H2 ** 2)
    M = np.cov(np.vstack([Z, H1, H2]))
    w = np.sort(np.linalg.eigvalsh(M))[::-1] + 1e-20
    return dict(
        hv_ratio=float(eh / ez),
        rectilinearity=float(1 - (w[1] + w[2]) / (2 * w[0])),
        planarity=float(1 - 2 * w[2] / (w[0] + w[1])),
    )


def compute_all(z, fs, t_pre, h1=None, h2=None):
    feats = {}
    feats.update(spectral_features(z, fs))
    feats.update(envelope_features(z, fs, t_pre))
    if h1 is not None and h2 is not None and len(h1) == len(z) == len(h2) and len(z) > 8:
        try:
            feats.update(polarization_features(z, h1, h2, fs))
        except Exception:
            pass
    return feats


FEATURE_ORDER = [
    "sp_centroid", "sp_spread", "sp_flatness", "sp_rolloff85", "sp_peakfreq",
    "bf_1_2", "bf_2_4", "bf_4_8", "bf_8_16", "bf_16up", "hf_ratio", "sp_slope", "sp_entropy",
    "env_crest", "env_kurt", "env_skew", "rise_time", "duration", "zcr", "snr",
    "hv_ratio", "rectilinearity", "planarity",
]
