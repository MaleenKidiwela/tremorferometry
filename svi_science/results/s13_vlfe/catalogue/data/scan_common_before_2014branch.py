"""VLFE scan shared code (rules: ../PLAN_SCAN.md). Grid, fast synthetic templates, v3 data loading with per-run
detrend/taper/band-pass and 300-s mask dilation, earthquake masks, season/episode helpers."""
from __future__ import annotations
import os, sys, re
import numpy as np, pandas as pd
from scipy.special import erf
from scipy.signal import butter, sosfiltfilt
from scipy.ndimage import binary_dilation
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_green as G

V3 = f'{C.DATA}/disp1sps_v3'
# PLAN_SCAN Amendment S8 (2026-10-15 PT): SCAN_YEAR selects the scan year (default 2011 -> every 2011 path/value unchanged).
YEAR = int(os.environ.get('SCAN_YEAR', '2011'))
SD = f'{C.DATA}/scan' if YEAR == 2011 else f'{C.DATA}/scan_{YEAR}'
os.makedirs(SD, exist_ok=True)


def tab(name):
    """Year-specific table path: unchanged for 2011; for 2012 'scan_x'/'scan2_x' -> 'scan3_x' (Amendment S8)."""
    if YEAR == 2011:
        return f'{C.TB}/{name}'
    for pre in ('scan2_', 'scan_'):
        if name.startswith(pre):
            return f'{C.TB}/scan3_' + name[len(pre):]
    return f'{C.TB}/scan3_' + name
PILOT = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV']
NOISE_ONLY = ['CN.LZB']
LAT0, LON0 = 48.37, -123.61
FIT0, FIT1 = -60, 120
TAU = np.arange(FIT0, FIT1 + 1)            # 181 samples, origin at index 60
NT = TAU.size
MECHS = {'thrust': (315, 12, 90), 'strikeslip': (315, 90, 0), 'steepdip': (315, 80, 90)}
ELEM = [np.diag([1., -1, 0]), np.array([[0, 1., 0], [1, 0, 0], [0, 0, 0]]), np.array([[0, 0, 1.], [0, 0, 0], [1, 0, 0]]),
        np.array([[0, 0, 0], [0, 0, 1.], [0, 1, 0]]), np.diag([-1., -1, 2]) / np.sqrt(3)]
DIL = 300
TAPER = 120
MINRUN = 1800
Y0 = pd.Timestamp('2011-01-01')
NDAY = 365
ETS = (pd.Timestamp('2011-08-12'), pd.Timestamp('2011-09-04'))     # 2011-08-12 .. 2011-09-03 inclusive
HG = [('2', '2011-08-21 04:09:06', 48.3, -123.6, 54, 3.7), ('3', '2011-08-21 08:43:29', 48.4, -123.45, 52, 3.3),
      ('4', '2011-08-24 10:47:38', 48.35, -123.625, 30, 3.5), ('5', '2011-08-26 13:44:34', 48.35, -123.175, 55, 3.4)]
if YEAR == 2012:                                                    # Amendment S8.2
    Y0 = pd.Timestamp('2012-01-01')
    NDAY = 366
    ETS = (pd.Timestamp('2012-08-31'), pd.Timestamp('2012-09-29'))  # SVI major 2012-08: 2012-08-31 .. 2012-09-28 inclusive
    HG = []                                                         # H&G list no 2012 events (Table S1 not on disk)
elif YEAR != 2011:
    raise ValueError(f'SCAN_YEAR {YEAR} not configured')


def kmxy(lat, lon):
    kx = 111.19 * np.cos(np.radians(LAT0))
    return (np.asarray(lon) - LON0) * kx, (np.asarray(lat) - LAT0) * 111.19


def dist_km(la1, lo1, la2, lo2):
    x1, y1 = kmxy(la1, lo1); x2, y2 = kmxy(la2, lo2)
    return np.hypot(x1 - x2, y1 - y2)


def slab_depth(lat, lon):
    S = np.loadtxt(f'{C.ROOT}/data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz', delimiter=',')
    S = S[np.isfinite(S[:, 2])]
    lo = np.where(S[:, 0] > 180, S[:, 0] - 360, S[:, 0])
    from scipy.spatial import cKDTree
    t = cKDTree(np.c_[lo * np.cos(np.radians(LAT0)), S[:, 1]])
    d, i = t.query(np.c_[np.asarray(lon) * np.cos(np.radians(LAT0)), np.asarray(lat)])
    dep = -S[i, 2]
    return np.where(d < 0.06, dep, np.nan)


def grid():
    kx = 111.19 * np.cos(np.radians(LAT0))
    ys = np.arange((47.8 - LAT0) * 111.19, (49.3 - LAT0) * 111.19 + 1e-6, 10.0)
    xs = np.arange((-124.8 - LON0) * kx, (-122.6 - LON0) * kx + 1e-6, 10.0)
    X, Y = np.meshgrid(xs, ys)
    lat = LAT0 + Y.ravel() / 111.19; lon = LON0 + X.ravel() / kx
    dep = slab_depth(lat, lon)
    ok = np.isfinite(dep)
    N = pd.DataFrame(dict(lat=lat[ok], lon=lon[ok], dep=dep[ok]))
    N.loc[len(N)] = (48.37444982698962, -123.61085121107266, 34.66806228373703)    # S13 2011-08 pooled centroid (6b)
    N['node'] = np.arange(len(N))
    return N


def synth_fast(M, src_ned, rec_ned, t, tau=1.0, vp=G.VP, vs=G.VS, rho=G.RHO, free_surface=2.0):
    """Vectorised copy of s13_green.synth (identical formula); returns N, E, Up."""
    d = np.asarray(rec_ned, float) - np.asarray(src_ned, float)
    r = np.linalg.norm(d); g = d / r
    dt = t[1] - t[0]
    mf = lambda tt: 0.5 * (1 + erf(tt / tau))                       # noqa: E731
    mr = lambda tt: np.exp(-(tt / tau) ** 2) / (tau * np.sqrt(np.pi))  # noqa: E731
    a, b = vp, vs
    taus = np.arange(r / a, r / b + dt, dt / 4)
    I = np.trapezoid(taus[None, :] * mf(t[:, None] - taus[None, :]), taus, axis=1)
    u = np.zeros((3, t.size))
    gMg = g @ M @ g; Mg = M @ g; trM = np.trace(M)
    for n in range(3):
        AN = 15 * g[n] * gMg - 3 * g[n] * trM - 6 * Mg[n]
        AIP = 6 * g[n] * gMg - g[n] * trM - 2 * Mg[n]
        AIS = -(6 * g[n] * gMg - g[n] * trM - 3 * Mg[n])
        AFP = g[n] * gMg
        AFS = -(g[n] * gMg - Mg[n])
        u[n] = (AN * I / r ** 4 + AIP * mf(t - r / a) / (a ** 2 * r ** 2) + AIS * mf(t - r / b) / (b ** 2 * r ** 2)
                + AFP * mr(t - r / a) / (a ** 3 * r) + AFS * mr(t - r / b) / (b ** 3 * r)) / (4 * np.pi * rho)
    u = u * free_surface
    return np.vstack([u[0], u[1], -u[2]])


TS = np.arange(-400, 600, 0.25)


def template(M, src, sta_ll, tau=1.0, out_t=TAU, **kw):
    """Band-passed (0.02-0.05 Hz) displacement N,E,Z for moment tensor M (N m) at src=(lat,lon,dep_km) observed at
    sta_ll=(lat,lon), sampled at out_t (s after origin). Same recipe as s13_stack.templates."""
    n, e = G.ned_offset(src[0], src[1], sta_ll[0], sta_ll[1])
    u = G.bandpass(synth_fast(M, (0, 0, src[2] * 1e3), (n, e, 0.0), TS, tau=tau, **kw), 4.0)
    return np.vstack([np.interp(out_t, TS, c) for c in u])


def raw_disp(M, src, sta_ll, t, tau=1.0, **kw):
    """UNFILTERED displacement N,E,Z at times t (s after origin) — for injection."""
    n, e = G.ned_offset(src[0], src[1], sta_ll[0], sta_ll[1])
    return synth_fast(M, (0, 0, src[2] * 1e3), (n, e, 0.0), t, tau=tau, **kw)


def vetoes(netsta):
    V = pd.read_csv(f'{C.ROOT}/svi_science/instrument_vetoes.csv')
    V = V[(V.decision.astype(str).str.strip() == 'veto') & (V.netsta == netsta)]
    out = []
    for r in V.itertuples():
        ch = set(re.findall(r'\b[A-Z]{2}([ZNE12])\b', str(r.channels)))
        comps = set()
        for c in ch:
            comps |= {'Z': {2}, 'N': {0}, 'E': {1}, '1': {0, 1}, '2': {0, 1}}[c]
        out.append((pd.Timestamp(r.start), pd.Timestamp(r.end), sorted(comps or {0, 1, 2})))
    return out


SOS = butter(4, [0.02, 0.05], btype='band', fs=1.0, output='sos')


def runs(mask):
    m = np.concatenate([[False], mask, [False]]).astype(np.int8)
    d = np.diff(m)
    return np.where(d == 1)[0], np.where(d == -1)[0]


def bandpass_runs(x):
    """x (n,) with NaN gaps -> band-passed per contiguous finite run (detrend, 120-s cosine taper, filtfilt), NaN within
    DIL s of any NaN; runs < MINRUN dropped."""
    out = np.full(x.shape, np.nan)
    fin = np.isfinite(x)
    s0, s1 = runs(fin)
    for a, b in zip(s0, s1):
        if b - a < MINRUN:
            continue
        y = x[a:b].astype(np.float64)
        tt = np.arange(b - a)
        p = np.polyfit(tt, y, 1); y = y - np.polyval(p, tt)
        w = np.ones(b - a); k = np.arange(TAPER)
        ramp = 0.5 * (1 - np.cos(np.pi * k / TAPER))
        w[:TAPER] = ramp; w[-TAPER:] = ramp[::-1]
        out[a:b] = sosfiltfilt(SOS, y * w)
    bad = ~np.isfinite(out) | ~fin
    if bad.any():
        out[binary_dilation(bad, structure=np.ones(2 * DIL + 1, bool))] = np.nan
    return out


def load_year(netsta, root=V3, year=None, inject=None):
    """Continuous 1-sps N,E,Z for the year (3, NDAY*86400) float64, raw (NaN gaps), vetoes applied, optional
    injection callback inject(X) adding signals in place (before band-pass)."""
    year = YEAR if year is None else year
    X = np.full((3, NDAY * 86400), np.nan)
    days = pd.date_range(f'{year}-01-01', periods=NDAY, freq='D')
    VV = vetoes(netsta)
    for i, d in enumerate(days):
        f = f'{root}/{netsta}/{d.year}/{d.dayofyear:03d}.npz'
        if os.path.isfile(f):
            X[:, i * 86400:(i + 1) * 86400] = np.load(f)['d']
            for a0, a1, comps in VV:
                if a0 <= d <= a1:
                    X[comps, i * 86400:(i + 1) * 86400] = np.nan
    if inject is not None:
        inject(netsta, X)
    return X


def bandpass_year(X):
    return np.vstack([bandpass_runs(X[c]) for c in range(3)])


def mask_vector(tier='head', year=None):
    """Boolean per-second mask over the year for tier in {'head','m45','t58'} (each united with 'reg')."""
    year = YEAR if year is None else year
    E = pd.read_csv(f'{C.TB}/scan_eq_events.csv')
    y0 = pd.Timestamp(f'{year}-01-01').value / 1e9
    m = np.zeros(NDAY * 86400, bool)
    for t in (tier, 'reg'):
        a, b = E[f'{t}_t0'].values, E[f'{t}_t1'].values
        ok = np.isfinite(a) & (b > y0) & (a < y0 + len(m))
        for s, e in zip(a[ok], b[ok]):
            m[int(max(s - y0, 0)):int(min(e - y0, len(m)))] = True
    return m


def episodes_svi():
    E = pd.read_csv(f'{C.ROOT}/svi_science/results/b_ets_record/episodes_b.csv')
    return E[E.region == 'SVI']


def season(ts):
    m = pd.DatetimeIndex(ts).month
    return np.where((m >= 4) & (m <= 9), 'summer', 'winter')
