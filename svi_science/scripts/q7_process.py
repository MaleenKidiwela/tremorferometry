#!/usr/bin/env python
"""q7_process (2026-10-07 PT) -- raw PBO gauge counts -> hourly tide- and pressure-corrected gauge strain.

Rules in svi_science/results/q7_strainmeter/PLAN.md ("Processing").  Per station:
    counts -> linearised gauge strain e = (0.2/87) * r/(1-r), r = counts/1e8   (UNAVCO GTSM formula)
    hourly means (>= 4 of 6 samples); pressure RDO hourly mean
    per station-year: high-pass (minus 72-h running median), LSQ of 8 tidal constituents + pressure admittance,
        2 passes with 5-MAD rejection; subtract fitted tides + pressure from the full series
    spikes: |Y - 5-h running median| > max(10 MAD, 20 ne) -> NaN
    steps: hourly increments (across gaps: between consecutive valid samples) > max(50 MAD, 30 ne) -> set to 0,
        series re-integrated per segment (the threshold per station-year is saved; injection uses it)
Output data/bsm/processed/<STA>_hourly.npz:  t0 (hours since 2009-01-01 UTC), Y (nh,4) float32 corrected,
    Yraw (nh,4) float32 (tide/pressure corrected, NOT de-stepped), thr (nyear,4), az (4,), qc json.

    python -P svi_science/scripts/q7_process.py B009 [B010 ...]
"""
import sys, os, json, glob, warnings
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
RAW = f"{ROOT}/data/bsm/raw"
OUT = f"{ROOT}/data/bsm/processed"
os.makedirs(OUT, exist_ok=True)
T0 = np.datetime64("2009-01-01T00:00")
T1 = np.datetime64("2026-10-07T00:00")
NH = int((T1 - T0).astype("timedelta64[h]").astype(int))
YEARS = np.arange(2009, 2027)
# tidal constituent angular speeds, degrees per hour
TIDES = dict(M2=28.9841042, S2=30.0, N2=28.4397295, K2=30.0821373, K1=15.0410686, O1=13.9430356, P1=14.9589314,
             Q1=13.3986609)
GAP, DIA = 0.2, 87.0


def hourly(tr, fs_per_h):
    """Masked trace -> hourly mean on the global grid, NaN where < 2/3 of samples."""
    x = np.ma.filled(tr.data.astype(float), np.nan) if np.ma.isMaskedArray(tr.data) else tr.data.astype(float)
    t = np.datetime64(tr.stats.starttime.datetime) + (np.arange(x.size) * tr.stats.delta * 1e6).astype("timedelta64[us]")
    h = ((t - T0) / np.timedelta64(1, "h")).astype(float)
    hi = np.floor(h).astype(int)
    ok = np.isfinite(x) & (hi >= 0) & (hi < NH)
    s = np.bincount(hi[ok], weights=x[ok], minlength=NH)
    n = np.bincount(hi[ok], minlength=NH)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = s / n
    m[n < max(1, int(np.ceil(fs_per_h * 2 / 3)))] = np.nan
    return m


def runmed(x, w):
    return pd.Series(x).rolling(w, center=True, min_periods=w // 3).median().to_numpy()


def mad(x):
    x = x[np.isfinite(x)]
    return float(np.median(np.abs(x - np.median(x)))) * 1.4826 if x.size else np.nan


def process(sta):
    import obspy
    files = sorted(glob.glob(f"{RAW}/{sta}/{sta}.*.mseed"))
    st = obspy.Stream()
    for f in files:
        try:
            st += obspy.read(f)
        except Exception as e:  # noqa
            print(sta, "read fail", f, e)
    if len(st) == 0:
        return None
    try:
        st.merge(method=0)
    except Exception:
        groups = {}
        for tr in st:
            groups.setdefault((tr.id, round(tr.stats.sampling_rate, 8)), obspy.Stream()).append(tr)
        st = obspy.Stream()
        for g in groups.values():
            g.merge(method=0); st += g
    inv_az = {}
    try:
        from obspy import read_inventory
        inv = read_inventory(f"{ROOT}/data/bsm/stationxml/{sta}.xml")
        for ch in inv[0][0]:
            inv_az[ch.code] = ch.azimuth
    except Exception as e:  # noqa
        print(sta, "no stationxml", e)
    E = np.full((NH, 4), np.nan)
    for k in range(4):
        trs = [tr for tr in st if tr.stats.channel == f"RS{k + 1}"]
        acc = np.full(NH, np.nan)
        for tr in trs:   # several location codes / epochs possible
            r = tr.copy()
            r.data = np.ma.masked_invalid(np.ma.filled(r.data.astype(float), np.nan) / 1e8) if np.ma.isMaskedArray(r.data) \
                else r.data.astype(float) / 1e8
            rr = np.ma.filled(r.data, np.nan) if np.ma.isMaskedArray(r.data) else r.data
            r.data = (GAP / DIA) * rr / (1 - rr)
            h = hourly(r, 6)
            acc = np.where(np.isfinite(acc), acc, h)
        E[:, k] = acc
    P = np.full(NH, np.nan)
    for tr in [tr for tr in st if tr.stats.channel == "RDO"]:
        h = hourly(tr, 2)
        P = np.where(np.isfinite(P), P, h)
    th = np.arange(NH, dtype=float)
    yr = (T0 + th.astype("timedelta64[h]")).astype("datetime64[Y]").astype(int) + 1970
    Ycorr = np.full_like(E, np.nan)
    qc = dict(station=sta, az={c: inv_az.get(c) for c in ["RS1", "RS2", "RS3", "RS4"]}, years={})
    thr = np.full((YEARS.size, 4), np.nan)
    for iy, y in enumerate(YEARS):
        sel = np.flatnonzero(yr == y)
        Ey, Py = E[sel], P[sel]
        nday = np.isfinite(Ey).all(1).sum() / 24
        if nday < 60:
            qc["years"][int(y)] = dict(days=float(nday), used=False)
            continue
        cols = [np.ones(sel.size)]
        for w in TIDES.values():
            a = np.deg2rad(w) * th[sel]
            cols += [np.cos(a), np.sin(a)]
        Xt = np.column_stack(cols)
        usep = np.isfinite(Py).mean() > 0.5
        Php = Py - runmed(Py, 72) if usep else None
        info = dict(days=float(nday), used=True, pressure=bool(usep), rms_resid=[], adm=[])
        for k in range(4):
            e = Ey[:, k]
            ehp = e - runmed(e, 72)
            X = np.column_stack([Xt, Php]) if usep else Xt
            ok = np.isfinite(ehp) & np.isfinite(X).all(1)
            if ok.sum() < 24 * 30:
                info["rms_resid"].append(None); info["adm"].append(None); continue
            for _ in range(3):
                c, *_ = np.linalg.lstsq(X[ok], ehp[ok], rcond=None)
                res = ehp - np.where(np.isfinite(X), X, 0) @ c
                m = mad(res[ok])
                ok = ok & (np.abs(res) < 5 * m)
            Xf = X.copy()
            if usep:
                Xf[:, -1] = np.where(np.isfinite(Py), Py - np.nanmedian(Py), np.nan)   # full pressure, not HP
            model = Xf[:, 1:] @ c[1:]                      # no constant
            Ycorr[sel, k] = e - model
            info["rms_resid"].append(float(m)); info["adm"].append(float(c[-1]) if usep else None)
        # dead-gauge QC (added 2026-10-07 ~5:40 pm PT after the diagnostic plot showed flat gauges, e.g. B926 g3):
        # a gauge-year whose high-passed residual RMS is < 0.25 x the median of the station's gauges is dead -> NaN
        rr = np.array([np.nan if v is None else v for v in info["rms_resid"]], float)
        if np.isfinite(rr).sum() >= 3:
            med = np.nanmedian(rr)
            dead = np.flatnonzero(rr < 0.25 * med)
            for k in dead:
                Ycorr[sel, k] = np.nan
            info["dead_gauges"] = [int(k) + 1 for k in dead]
        qc["years"][int(y)] = info
    # spikes
    Yr = Ycorr.copy()
    for k in range(4):
        x = Yr[:, k]
        dev = x - runmed(x, 5)
        for iy, y in enumerate(YEARS):
            sel = yr == y
            m = mad(dev[sel])
            if np.isfinite(m):
                bad = sel & (np.abs(dev) > max(10 * m, 20e-9))
                x[bad] = np.nan
        Yr[:, k] = x
    # steps
    Y = np.full_like(Yr, np.nan)
    nstep = np.zeros(4, int)
    for k in range(4):
        x = Yr[:, k]
        v = np.flatnonzero(np.isfinite(x))
        if v.size < 2:
            continue
        d = np.diff(x[v])
        yv = yr[v[1:]]
        dd = d.copy()
        for iy, y in enumerate(YEARS):
            s = yv == y
            if s.sum() < 100:
                continue
            m = mad(d[s])
            thr[iy, k] = max(50 * m, 30e-9)
            big = s & (np.abs(d) > thr[iy, k])
            dd[big] = 0.0
            nstep[k] += int(big.sum())
        out = np.empty(v.size); out[0] = x[v[0]]; out[1:] = x[v[0]] + np.cumsum(dd)
        Y[v, k] = out
    qc["n_steps_removed"] = nstep.tolist()
    np.savez_compressed(f"{OUT}/{sta}_hourly.npz", Y=Y.astype(np.float32), Yraw=Yr.astype(np.float32),
                        thr=thr, years=YEARS, az=np.array([inv_az.get(f"RS{k+1}", np.nan) for k in range(4)], float),
                        P=P.astype(np.float32))
    json.dump(qc, open(f"{OUT}/{sta}_qc.json", "w"), indent=1)
    print(sta, "done; steps removed", nstep.tolist(), "valid days per gauge",
          (np.isfinite(Y).sum(0) / 24).round(0).tolist(), flush=True)


if __name__ == "__main__":
    warnings.simplefilter("ignore")
    for s in sys.argv[1:]:
        process(s)
