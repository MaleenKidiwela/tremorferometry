#!/usr/bin/env python
"""S12c step 3a: daily single-station ambient-noise auto- and cross-component correlations (0-30 s lag).

Input : raw day mseed already on disk (no download)
        2014+ : lfe_extension/traces/<NET.STA>/<year>/<doy>.mseed      (HH 100 Hz)
        2003-2013 : lfe_extension/traces_backward/<NET.STA>/<year>/<doy>.mseed (BH 40 Hz / HH)
Output: data/<NET.STA>/<year>.npz  (one file per station-year; written once at the end of the year, resumable)
        dates (ndays,) 'YYYY-MM-DD'; nh_<band> (ndays,) hours used;
        <band>_<pair> (ndays, 24, nlag) float32 HOURLY normalised correlations, NaN = hour rejected
        (hourly kept so tremor hours / noisy hours can be masked later; v0 daily-mean format is in data_v0_daily/);
        rms_<band> and rms_b2_45 (ndays, 24, 3) hourly RMS (counts) of Z,R,T: 2-4.5 Hz = tremor proxy.
        autos ZZ RR TT : lags 0..30 s ; cross ZR ZT RT : lags -30..+30 s (positive lag = second comp later).
Recipe (fixed in PLAN.md before measuring):
  merge (gaps -> 0) -> demean -> decimate to 10 Hz (anti-alias resample_poly) -> rotate N,E to R = slab dip direction
  (RF result dipdir), T = R+90 -> per band b in {0.5-2 Hz, 1-4 Hz}: zero-phase Butterworth bandpass ->
  running-absolute-mean normalisation (window 1/f1 s, weight = max over the 3 components, keeps relative amplitudes)
  -> 1-hour windows; skip a window if it holds >1 % zero-fill or its pre-normalisation RMS is >4x or <0.25x the day's
  median hourly RMS (quakes, glitches, dead) -> spectral whitening by a boxcar-smoothed amplitude spectrum
  (W = 0.4 Hz for 0.5-2, 0.8 Hz for 1-4; the boxcar averages out spectral ripple from echoes at lags > 2.5/1.25 s,
  so whitening keeps reflections) common to all 3 components (mean of the 3 smoothed spectra; keeps relative phase and
  amplitude between components) -> cosine-tapered band window -> correlations, normalised by sqrt of zero-lag autos ->
  mean over windows.
"""
from __future__ import annotations
import os, sys, glob, json, argparse, time
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import numpy as np

ROOT = "/home/jovyan/tremorferometry"
WS = f"{ROOT}/svi_science/results/s12_lvl/acf_pgc"
FS = 10.0
MAXLAG = 30.0
NL = int(MAXLAG * FS)
WIN = 3600
BANDS = {"b05_2": (0.5, 2.0, 0.4), "b1_4": (1.0, 4.0, 0.8)}
AUTOS = ["ZZ", "RR", "TT"]
CROSS = ["ZR", "ZT", "RT"]
DIPDIR = {}


def dipdir(ns):
    if ns not in DIPDIR:
        f = f"{ROOT}/svi_science/results/s12_lvl/rf_gosselin/tables/full/{ns}_result.json"
        DIPDIR[ns] = float(json.load(open(f))["dipdir"]) if os.path.exists(f) else 55.0
    return DIPDIR[ns]


def day_path(ns, y, d):
    for sub in ("traces", "traces_backward"):
        for ext in (".mseed",):
            p = f"{ROOT}/lfe_extension/{sub}/{ns}/{y}/{d:03d}{ext}"
            if os.path.exists(p) and os.path.getsize(p) > 0:
                return p
    return None


def load_day(p):
    from obspy import read, UTCDateTime
    from scipy.signal import resample_poly
    st = read(p)
    for tr in st:
        r = float(tr.stats.sampling_rate)
        if r > 0 and abs(r - round(r)) / r < 1e-4:
            tr.stats.sampling_rate = float(round(r))
        tr.data = tr.data.astype(np.float64)
        tr.data -= np.median(tr.data)
    st.merge(method=1, fill_value=None)
    t0 = min(tr.stats.starttime for tr in st) + 60
    day0 = UTCDateTime(t0.year, t0.month, t0.day)
    st.trim(day0, day0 + 86400 - 1e-6, pad=True, fill_value=None)
    comp = {}
    for tr in st:
        c = tr.stats.channel[-1]
        c = {"1": "N", "2": "E"}.get(c, c)
        if c not in "ZNE" or c in comp:
            continue
        fs0 = tr.stats.sampling_rate
        q = fs0 / FS
        if abs(q - round(q)) > 1e-6:
            return None
        x = np.ma.filled(np.ma.asarray(tr.data, dtype=np.float64), 0.0)
        z = np.ma.getmaskarray(tr.data)
        x = resample_poly(x, 1, int(round(q))) if q > 1 else x
        z = z[: len(z) // int(round(q)) * int(round(q))].reshape(-1, int(round(q))).any(1) if q > 1 else z
        n = int(86400 * FS)
        xx = np.zeros(n); zz = np.ones(n, bool)
        m = min(n, len(x), len(z)); xx[:m] = x[:m]; zz[:m] = z[:m]
        comp[c] = (xx, zz)
    if not all(c in comp for c in "ZNE"):
        return None
    return comp, str(day0.date)


def bandpass(x, f1, f2):
    from scipy.signal import butter, sosfiltfilt
    sos = butter(4, [f1 / (FS / 2), f2 / (FS / 2)], btype="band", output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def ram(X, f1):
    from scipy.ndimage import uniform_filter1d
    n = max(3, int(round(FS / f1)))
    w = uniform_filter1d(np.abs(X), size=n, axis=-1).max(0)
    w[w <= 0] = np.inf
    return X / w[None, :]


def process_day(args):
    ns, y, d = args
    p = day_path(ns, y, d)
    if p is None:
        return None
    try:
        r = load_day(p)
    except Exception as e:
        return ("ERR", f"{type(e).__name__}: {e}")
    if r is None:
        return None
    comp, date = r
    th = np.deg2rad(dipdir(ns))
    Z, N, E = comp["Z"][0], comp["N"][0], comp["E"][0]
    R = np.cos(th) * N + np.sin(th) * E
    T = -np.sin(th) * N + np.cos(th) * E
    zero = comp["Z"][1] | comp["N"][1] | comp["E"][1]
    X0 = np.vstack([Z, R, T])
    nwin = int(86400 // WIN); L = int(WIN * FS); nfft = 2 * L
    f = np.fft.rfftfreq(nfft, 1 / FS)
    out = {"date": date}
    zfrac = zero[: nwin * L].reshape(nwin, L).mean(1)
    for bn, (f1, f2, W) in BANDS.items():
        Xb = bandpass(X0, f1, f2)
        rms = np.sqrt((Xb[:, : nwin * L].reshape(3, nwin, L) ** 2).mean(-1)).max(0)
        ok = zfrac < 0.01
        if ok.sum() < 4:
            out[bn] = None; out[f"rms_{bn}"] = np.full((nwin, 3), np.nan, np.float32); continue
        med = np.median(rms[ok])
        ok &= (rms < 4 * med) & (rms > 0.25 * med)
        Xn = ram(Xb, f1)
        nbw = max(3, int(round(W / (f[1] - f[0]))))
        from scipy.ndimage import uniform_filter1d
        taper = np.zeros_like(f)
        fw = 0.1 * (f2 - f1)
        m = (f >= f1) & (f <= f2); taper[m] = 1
        lo = (f > f1 - fw) & (f < f1); taper[lo] = 0.5 * (1 - np.cos(np.pi * (f[lo] - (f1 - fw)) / fw))
        hi = (f > f2) & (f < f2 + fw); taper[hi] = 0.5 * (1 + np.cos(np.pi * (f[hi] - f2) / fw))
        acc = {k: np.full((nwin, NL + 1), np.nan, np.float32) for k in AUTOS}
        acc.update({k: np.full((nwin, 2 * NL + 1), np.nan, np.float32) for k in CROSS})
        nh = 0
        tw = np.hanning(L)  # window taper
        for i in np.where(ok)[0]:
            w = Xn[:, i * L:(i + 1) * L]
            w = (w - w.mean(1, keepdims=True)) * tw[None, :]
            S = np.fft.rfft(w, nfft, axis=-1)
            sm = uniform_filter1d(np.abs(S), size=nbw, axis=-1).mean(0)
            sm = sm + 1e-3 * sm.max()
            S = S / sm[None, :] * taper[None, :]
            idx = {"Z": 0, "R": 1, "T": 2}
            a0 = {}
            for k in AUTOS:
                c = np.fft.irfft(S[idx[k[0]]] * np.conj(S[idx[k[0]]]), nfft)
                a0[k[0]] = c[0]
                acc[k][i] = c[: NL + 1] / c[0]
            for k in CROSS:
                # positive lag = second component arrives later: c(tau) = sum a(t) b(t+tau)
                c = np.fft.irfft(np.conj(S[idx[k[0]]]) * S[idx[k[1]]], nfft)
                c = np.concatenate([c[-NL:], c[: NL + 1]]) / np.sqrt(a0[k[0]] * a0[k[1]])
                acc[k][i] = c
            nh += 1
        out[bn] = (nh, acc)
        out[f"rms_{bn}"] = np.sqrt((Xb[:, : nwin * L].reshape(3, nwin, L) ** 2).mean(-1)).T.astype(np.float32)
    # tremor proxy: hourly RMS 2-4.5 Hz per component (Z, R, T), counts (no response), NaN where gappy
    Xh = bandpass(X0, 2.0, 4.5)
    rr = np.sqrt((Xh[:, : nwin * L].reshape(3, nwin, L) ** 2).mean(-1)).T
    rr[zfrac >= 0.01] = np.nan
    out["rms_b2_45"] = rr.astype(np.float32)
    return out


def run_year(ns, y, workers):
    od = f"{WS}/data/{ns}"; os.makedirs(od, exist_ok=True)
    fo = f"{od}/{y}.npz"
    if os.path.exists(fo):
        print(f"[skip] {fo} exists", flush=True); return
    days = []
    for sub in ("traces", "traces_backward"):
        days += [int(os.path.basename(p)[:3]) for p in glob.glob(f"{ROOT}/lfe_extension/{sub}/{ns}/{y}/*.mseed")]
    days = sorted(set(days))
    t = time.time()
    from multiprocessing import Pool
    res = []
    with Pool(workers) as pool:
        for r in pool.imap(process_day, [(ns, y, d) for d in days], chunksize=2):
            if r is None: continue
            if isinstance(r, tuple):
                print(f"  [err] {ns} {y}: {r[1]}", flush=True); continue
            res.append(r)
    res.sort(key=lambda r: r["date"])
    sv = {"dates": np.array([r["date"] for r in res])}
    nw = int(86400 // WIN)
    sv["rms_b2_45"] = np.stack([r["rms_b2_45"] for r in res]) if res else np.zeros((0, nw, 3), np.float32)
    for bn in BANDS:
        sv[f"nh_{bn}"] = np.array([r[bn][0] if r[bn] else 0 for r in res], dtype=np.int16)
        sv[f"rms_{bn}"] = np.stack([r[f"rms_{bn}"] for r in res]) if res else np.zeros((0, nw, 3), np.float32)
        for k in AUTOS + CROSS:
            n = NL + 1 if k in AUTOS else 2 * NL + 1
            sv[f"{bn}_{k}"] = np.stack([r[bn][1][k] if r[bn] else np.full((nw, n), np.nan, np.float32) for r in res]) \
                if res else np.zeros((0, nw, n), np.float32)
    np.savez(fo + ".tmp.npz", **sv); os.replace(fo + ".tmp.npz", fo)
    print(f"[done] {ns} {y}: {len(res)}/{len(days)} days in {time.time() - t:.0f} s -> {fo}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", required=True)
    ap.add_argument("--years", required=True, help="e.g. 2018-2021 or 2003,2005")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--test", default=None, help="NET.STA,year,doy : one day, print timing")
    a = ap.parse_args()
    if a.test:
        ns, y, d = a.test.split(","); t = time.time()
        r = process_day((ns, int(y), int(d)))
        print(r["date"], {bn: r[bn][0] for bn in BANDS}, f"{time.time() - t:.1f} s")
        sys.exit()
    ys = a.years
    years = list(range(int(ys.split("-")[0]), int(ys.split("-")[1]) + 1)) if "-" in ys else [int(v) for v in ys.split(",")]
    for y in years:
        for ns in a.stations.split(","):
            run_year(ns, y, a.workers)
