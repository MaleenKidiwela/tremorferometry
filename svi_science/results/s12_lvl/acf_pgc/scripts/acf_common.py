"""S12c shared helpers: loading hourly ACFs, hour masks, stacks, windowed shift / stretch / CC measurement."""
from __future__ import annotations
import os, glob, json
import numpy as np, pandas as pd

ROOT = "/home/jovyan/tremorferometry"
WS = f"{ROOT}/svi_science/results/s12_lvl/acf_pgc"
FS = 10.0
NL = 300
LAG_AUTO = np.arange(NL + 1) / FS
LAG_X = np.arange(-NL, NL + 1) / FS
STA_LL = {"CN.PGC": (48.6478, -123.4511), "CN.VGZ": (48.4131, -123.3244), "CN.SNB": (48.7747, -123.1717)}
AUTOS = ["ZZ", "RR", "TT"]; CROSS = ["ZR", "ZT", "RT"]

# PGC windows (PLAN §4b); other stations: shifted by (their band centre - PGC band centre), see windows_for()
WIN_PGC = {
    "SS": {"pre": (13.0, 17.5), "lvl": (22.5, 26.0), "ctl_a": (9.0, 12.5), "ctl_b": (13.0, 17.5)},
    "PS": {"pre": (10.0, 13.5), "lvl": (17.5, 20.5), "ctl_a": (6.0, 9.5), "ctl_b": (10.0, 13.5)},
    "PP": {"pre": (6.0, 9.5), "lvl": (12.5, 15.5), "ctl_a": (3.0, 5.5), "ctl_b": (6.0, 9.5)},
}
PAIR_PHASE = {"ZZ": "PP", "RR": "SS", "TT": "SS", "ZR+": "PS", "ZR-": "PS", "ZT+": "PS", "ZT-": "PS"}


def windows_for(ns):
    """Shift the PGC LVL windows by the station's predicted bottom-time offset (top window kept below its top)."""
    if ns == "CN.PGC":
        return WIN_PGC
    B = pd.read_csv(f"{WS}/tables/predicted_lag_bands.csv")
    out = {}
    for ph in WIN_PGC:
        c = {}
        for side, iface in [("lvl", "bottom"), ("pre", "top")]:
            g = lambda s: B[(B.station == s) & (B.interface == iface)][[f"{ph}_lo", f"{ph}_hi"]].values.mean()
            d = g(ns) - g("CN.PGC")
            a, b = WIN_PGC[ph][side]; c[side] = (round(a + d, 1), round(b + d, 1))
        da = c["pre"][0] - WIN_PGC[ph]["pre"][0]
        for k in ("ctl_a", "ctl_b"):
            a, b = WIN_PGC[ph][k]; c[k] = (round(max(2.0, a + da), 1), round(b + da, 1))
        out[ph] = c
    return out


def load(ns, years, band, pairs=None, data_dir=None):
    """Return dates (datetime64[D]), dict pair -> (nd, 24, nlag) float32, dict rms arrays."""
    data_dir = data_dir or f"{WS}/data"
    pairs = pairs or AUTOS + CROSS
    D, H, R = [], {p: [] for p in pairs}, {"rms_b2_45": [], f"rms_{band}": []}
    for y in years:
        f = f"{data_dir}/{ns}/{y}.npz"
        if not os.path.exists(f):
            continue
        z = np.load(f)
        if len(z["dates"]) == 0:
            continue
        D.append(z["dates"].astype("datetime64[D]"))
        for p in pairs:
            H[p].append(z[f"{band}_{p}"])
        for k in R:
            R[k].append(z[k])
    D = np.concatenate(D)
    return D, {p: np.concatenate(v) for p, v in H.items()}, {k: np.concatenate(v) for k, v in R.items()}


def tremor_hours(ns, dates, r_km=60.0, pad_min=30):
    """(nd, 24) bool: PNSN-catalog tremor epicentre within r_km of the station within +-pad_min of the hour.
    Catalog starts 2010-08 (earlier days: all False; use rms_proxy_hours)."""
    lat, lon = STA_LL[ns]
    c = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    r = np.hypot((c.lat - lat) * 111.2, (c.lon - lon) * 111.2 * np.cos(np.radians(lat)))
    t = pd.to_datetime(c.time[r < r_km]).values.astype("datetime64[m]")
    M = np.zeros((len(dates), 24), bool)
    di = {d: i for i, d in enumerate(dates.astype("datetime64[D]"))}
    hs = np.unique(np.concatenate([(t - np.timedelta64(pad_min, "m")).astype("datetime64[h]"),
                                   (t + np.timedelta64(pad_min, "m")).astype("datetime64[h]")]))
    for th in hs:
        i = di.get(th.astype("datetime64[D]"))
        if i is not None:
            M[i, int((th - th.astype("datetime64[D]").astype("datetime64[h]")).astype(int))] = True
    return M


def rms_proxy_hours(dates, rms_b2_45, q=90.0, half_doy=15):
    """(nd, 24) bool: Z RMS 2-4.5 Hz above the q-th percentile of all hours within +-half_doy days of the same day of
    year (over all loaded years). Flags tremor-like (and storm/cultural) hours; works before the PNSN catalog."""
    z = rms_b2_45[:, :, 0].astype(float)
    # gain-proof (2026-10-13 ~11:14 am PT, before the full run): divide by a +-30-day running median of daily medians,
    # so sensor / gain changes between eras (BH 40 Hz -> HH 100 Hz, counts) do not decide which hours are masked
    dm = pd.Series(np.nanmedian(np.where(np.isfinite(z), z, np.nan), 1), index=pd.to_datetime(dates))
    rm = dm.rolling("61D", center=True, min_periods=5).median().values
    z = z / rm[:, None]
    doy = (dates - dates.astype("datetime64[Y]")).astype(int)
    M = np.zeros(z.shape, bool)
    for d in np.unique(doy):
        sel = np.minimum(np.abs(doy - d), 366 - np.abs(doy - d)) <= half_doy
        v = z[sel]; thr = np.nanpercentile(v, q) if np.isfinite(v).any() else np.inf
        i = doy == d
        with np.errstate(invalid="ignore"):
            M[i] = z[i] > thr
    return M


def hour_mask(ns, dates, R, mode):
    """True = keep. mode: 'none' (all hours) | 'tremorfree' (PLAN amendment A.1: catalog OR rms proxy removed)."""
    if mode == "none":
        return None
    return ~(tremor_hours(ns, dates) | rms_proxy_hours(dates, R["rms_b2_45"]))


def day_stack(Hp, hour_ok=None):
    """Mean over valid hours -> (nd, nlag); days with no valid hour -> NaN row. Returns stack and n hours."""
    X = Hp.copy()
    if hour_ok is not None:
        X[~hour_ok] = np.nan
    n = np.isfinite(X[:, :, 0]).sum(1)
    with np.errstate(invalid="ignore"):
        S = np.nansum(X, 1) / np.maximum(n, 1)[:, None]
    S[n == 0] = np.nan
    return S, n


def split_pair(S, pair):
    """For cross pairs return positive ('+') or negative ('-') lag half as a 0..30 s trace."""
    if pair.endswith("+"):
        return S[..., NL:]
    if pair.endswith("-"):
        return S[..., NL::-1]
    return S


def _up(x, k=10):
    from scipy.signal import resample
    return resample(x, len(x) * k)


def measure(ref, cur, win, maxshift=0.2, up=10, n_iter=3):
    """Shift of cur vs ref in window win (lag s); + = cur arrives later.
    Primary (2026-10-13 ~11:18 am PT, before any full-record number): LINEARISED delay (first-order Taylor,
    c(t) = a r(t - d) => d = -sum w (c/a - r) r' / sum w r'^2, 3 iterations re-shifting cur), on stacks upsampled x10.
    It cannot cycle-skip, which the CC-peak pick did at 0.5-2 Hz (ACF period ~1 s; scrambled draws jumped by
    +-100 ms in the pilot). Valid for |d| << T/4 ~ 250 ms; the target is 40-60 ms. CC-peak shift (+-maxshift) kept as
    a second estimate. Returns (d_lin, cc_after_shift, cc_zero_shift, d_ccpeak)."""
    if not (np.all(np.isfinite(ref)) and np.all(np.isfinite(cur))):
        return np.nan, np.nan, np.nan, np.nan
    fs = FS * up
    r = _up(ref, up); c = _up(cur, up)
    t = np.arange(len(r)) / fs
    m = (t >= win[0]) & (t < win[1])
    i0, i1 = np.where(m)[0][[0, -1]]
    w = np.hanning(i1 - i0 + 1)
    tw = t[i0:i1 + 1]; rw = r[i0:i1 + 1]; rd = np.gradient(r, 1 / fs)[i0:i1 + 1]
    def cc(x, y):
        x = (x - x.mean()) * w; y = (y - y.mean()) * w
        return float(np.dot(x, y) / np.sqrt(np.dot(x, x) * np.dot(y, y) + 1e-30))
    cc0 = cc(rw, c[i0:i1 + 1])
    d = 0.0
    for _ in range(n_iter):
        cw = np.interp(tw + d, t, c)          # cur advanced by d -> should match ref if cur = ref delayed by d
        amp = np.dot(w * cw, rw) / np.dot(w * rw, rw)
        if amp == 0: break
        d += -np.dot(w * (cw / amp - rw), rd) / np.dot(w * rd, rd)
        if abs(d) > 0.5: break
    cw = np.interp(tw + d, t, c)
    # CC-peak estimate
    ns_ = int(maxshift * fs); best = []
    for s_ in range(-ns_, ns_ + 1):
        if i0 + s_ < 0 or i1 + s_ >= len(c):
            best.append(-np.inf); continue
        best.append(cc(rw, c[i0 + s_:i1 + s_ + 1]))
    best = np.array(best); k = int(np.argmax(best))
    dk = 0.0
    if 0 < k < len(best) - 1:
        y0, y1, y2 = best[k - 1:k + 2]; den = y0 - 2 * y1 + y2
        dk = 0.5 * (y0 - y2) / den if den != 0 else 0.0
    return float(d), cc(rw, cw), cc0, (k - ns_ + dk) / fs


def differential(ref, cur, W):
    """Primary (amendment A.2) ddt = dt_lvl - dt_pre (additive); ddt_str = stretch-scaled form of PLAN §4b.
    Control: same between ctl_a and pre (both above the LVL top)."""
    dp, cp, cp0, dpk = measure(ref, cur, W["pre"])
    dl, cl, cl0, dlk = measure(ref, cur, W["lvl"])
    tp = np.mean(W["pre"]); tl = np.mean(W["lvl"])
    da, ca, _, _ = measure(ref, cur, W["ctl_a"])
    ta = np.mean(W["ctl_a"])
    return dict(dt_pre=dp, dt_lvl=dl, cc_pre=cp, cc_lvl=cl, cc0_pre=cp0, cc0_lvl=cl0,
                ddt=dl - dp, ddt_str=dl - dp * tl / tp, ddt_ccpk=dlk - dpk, dcc=cl0 - cp0,
                ddt_ctl=dp - da, ddt_ctl_str=dp - da * tp / ta, cc_ctl=ca)


def ets_mid():
    E = pd.read_csv(f"{ROOT}/svi_science/results/s12_lvl/rf_gosselin/tables/ets_mid_SVI_ALBH.csv")
    return pd.to_datetime(E["mid"]).values.astype("datetime64[D]")
