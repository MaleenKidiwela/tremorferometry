"""Grid-search centroid moment tensor inversion for hg_repro (PLAN.md row 10, 'changed': whole-space GF).

Whole-space moment-tensor displacement (same physics as ../scripts/s13_green.py: near + intermediate + far P and S
terms, free-surface factor 2, Vp 6.9, Vs 3.95 km/s, rho 3000), vectorised: for each node-station pair the five radial
time functions are computed once and the deviatoric basis (Kikuchi & Kanamori 5 elements) enters through the
coefficients. Moment function = smooth step 0.5(1+erf(t/tau)). Synthetics and data are band-passed identically
(Butterworth order 2 band, 0.02-0.05 Hz, zero phase). Window = 90 s from (origin + shift) at every station, shift
-10..+10 s at 1 s. Linear least squares per node and shift; best = max variance reduction VR = 1 - |d-Gm|^2/|d|^2.
Coarse grid 0.1 deg x 0.1 deg x 5 km (47.3-49.0 N, 124.4-122.4 W, 20-70 km), then fine 0.025 deg x 1 km (+-0.075 deg, +-3 km) around the best.
"""
import os, glob
import numpy as np
from scipy.special import erf
from scipy.signal import butter, sosfiltfilt
import hg_common as H

VP, VS, RHO, FS = 6900.0, 3950.0, 3000.0, 2.0
SOS = butter(2, H.BAND, btype='band', fs=1.0, output='sos')
TPRE, TPOST = 400, 1200                       # synthetic time axis (s) around the origin
TT = np.arange(-TPRE, TPOST, 1.0)
BASIS = np.array([[[0, 1, 0], [1, 0, 0], [0, 0, 0]],
                  [[1, 0, 0], [0, -1, 0], [0, 0, 0]],
                  [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
                  [[0, 0, 1], [0, 0, 0], [1, 0, 0]],
                  [[-1, 0, 0], [0, 0, 0], [0, 0, 1]]], float)
_COORD = None


def coords():
    """Station lat, lon (deg) from the StationXML files used by the prep."""
    global _COORD
    if _COORD is None:
        from obspy import read_inventory
        c = {}
        files = sorted(glob.glob(f'{H.S13}/data/resp/*_bb.xml')) + sorted(glob.glob(f'{H.HG}/data/resp/*.xml'))
        for f in files:
            inv = read_inventory(f)
            for n in inv:
                for s in n:
                    k = f'{n.code}.{s.code}'
                    if k in H.STAS and k not in c:
                        c[k] = (s.latitude, s.longitude)
        _COORD = c
    return _COORD


def ned(lat0, lon0, lat, lon):
    R = 6371e3
    return np.radians(lat - lat0) * R, np.radians(lon - lon0) * R * np.cos(np.radians(0.5 * (lat + lat0)))


def timefuncs(r, tau):
    """Five radial time functions (5, nt) for source-receiver distance r (m)."""
    a, b = VP, VS
    mf = lambda t: 0.5 * (1 + erf(t / tau))                          # noqa: E731
    mr = lambda t: np.exp(-(t / tau) ** 2) / (tau * np.sqrt(np.pi))   # noqa: E731
    taus = np.linspace(r / a, r / b, 64)
    I = np.trapezoid(taus[None, :] * mf(TT[:, None] - taus[None, :]), taus, axis=1)
    return np.vstack([I / r ** 4, mf(TT - r / a) / (a ** 2 * r ** 2), mf(TT - r / b) / (b ** 2 * r ** 2),
                      mr(TT - r / a) / (a ** 3 * r), mr(TT - r / b) / (b ** 3 * r)])


def green(src, rec, tau):
    """(5 basis, 3 comps N,E,Up, nt) band-passed displacement for unit basis tensors. src/rec = (n, e, down) m."""
    d = np.asarray(rec, float) - np.asarray(src, float)
    r = np.linalg.norm(d); g = d / r
    F = timefuncs(r, tau)
    out = np.zeros((5, 3, TT.size))
    for k, Mb in enumerate(BASIS):
        gMg = g @ Mb @ g; Mg = Mb @ g
        for n in range(3):
            co = np.array([15 * g[n] * gMg - 6 * Mg[n], 6 * g[n] * gMg - 2 * Mg[n], -(6 * g[n] * gMg - 3 * Mg[n]),
                           g[n] * gMg, -(g[n] * gMg - Mg[n])])
            out[k, n] = co @ F / (4 * np.pi * RHO) * FS
    out[:, 2] *= -1                                                  # down -> up
    return sosfiltfilt(SOS, out, axis=-1)


WIN0, WIN1 = -10, H.L + 10                     # stored GF / data window relative to the origin (s)


def gf_cache(nodes, tau, tag=None):
    """(nnode, nsta, 5, 3, 110) float32 GF cut to [origin-10, origin+100) for every node and every H.STAS station."""
    f = f'{H.MF}/gf_{tag}_tau{tau:g}.npy' if tag else None
    if f and os.path.isfile(f):
        return np.load(f)
    C = coords()
    out = np.zeros((len(nodes), len(H.STAS), 5, 3, WIN1 - WIN0), np.float32)
    for q, (la, lo, z) in enumerate(nodes):
        for si, s in enumerate(H.STAS):
            n_, e_ = ned(la, lo, *C[s])
            g = green((0.0, 0.0, z * 1e3), (n_, e_, 0.0), tau)
            out[q, si] = g[:, :, TPRE + WIN0:TPRE + WIN1]
    if f:
        np.save(f, out)
    return out


def invert(data, chans, nodes, G, shifts=range(-10, 11)):
    """data: (nch, 110) observed displacement on [origin-10, origin+100) (band-passed); chans 'NET.STA.C';
    G = gf_cache(nodes, tau). Returns list of dicts (one per node, best shift) sorted by VR."""
    comps = {'N': 0, 'E': 1, 'Z': 2}
    si = [H.STAS.index(c.rsplit('.', 1)[0]) for c in chans]
    ci = [comps[c[-1]] for c in chans]
    sh = np.array(list(shifts))
    D = np.stack([data[:, -WIN0 + k:-WIN0 + k + H.L].reshape(-1) for k in sh], axis=1).astype(np.float64)  # (n, nsh)
    dd = np.sum(D ** 2, axis=0)
    res = []
    for q, (la, lo, z) in enumerate(nodes):
        A = G[q][si, :, ci, -WIN0:-WIN0 + H.L]                     # (nch, 5, 90): GF on [origin+sh, +90)
        A2 = A.transpose(1, 0, 2).reshape(5, -1).T.astype(np.float64)
        Mm = np.linalg.lstsq(A2, D, rcond=None)[0]                  # (5, nsh)
        vr = 1 - np.sum((D - A2 @ Mm) ** 2, axis=0) / dd
        k = int(np.argmax(vr))
        res.append({'lat': la, 'lon': lo, 'dep': z, 'shift': int(sh[k]), 'vr': float(vr[k]), 'm': Mm[:, k]})
    return sorted(res, key=lambda d: -d['vr'])


def mt_from_m(m):
    return sum(mi * B for mi, B in zip(m, BASIS))


def describe(m):
    M = mt_from_m(m)
    w, v = np.linalg.eigh(M)
    m0 = np.sqrt(np.sum(M ** 2) / 2)
    mw = 2 / 3 * (np.log10(m0) - 9.1) if m0 > 0 else np.nan
    order = np.argsort(np.abs(w))
    eps = -w[order[0]] / abs(w[order[2]]) if w[order[2]] != 0 else np.nan
    P, T = v[:, 0], v[:, 2]                                         # NED eigenvectors (min, max)

    def axis(u):
        if u[2] < 0:
            u = -u
        return (np.degrees(np.arctan2(u[1], u[0])) % 360, np.degrees(np.arcsin(u[2])))
    return {'M0': m0, 'Mw': mw, 'clvd_pct': 200 * abs(eps), 'P_az': axis(P)[0], 'P_pl': axis(P)[1],
            'T_az': axis(T)[0], 'T_pl': axis(T)[1]}


def kagan(m1, m2):
    """Kagan angle (deg) between the double-couple frames of two moment tensors."""
    def frame(m):
        w, v = np.linalg.eigh(mt_from_m(m))
        v = v[:, [0, 1, 2]]
        if np.linalg.det(v) < 0:
            v[:, 1] *= -1
        return v
    A, B = frame(m1), frame(m2)
    best = 180.0
    for s in ([1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]):
        R = A @ np.diag(s) @ B.T
        ang = np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))
        best = min(best, ang)
    return best


def coarse_nodes():
    la = np.round(np.arange(47.3, 49.0001, 0.1), 3)
    lo = np.round(np.arange(-124.4, -122.3999, 0.1), 3)
    z = np.arange(20, 70.1, 5)
    return np.array([(a, b, c) for a in la for b in lo for c in z])


def fine_nodes(best):
    la = best['lat'] + np.arange(-0.075, 0.0751, 0.025)
    lo = best['lon'] + np.arange(-0.075, 0.0751, 0.025)
    z = best['dep'] + np.arange(-3, 3.1, 1)
    return np.array([(a, b, c) for a in la for b in lo for c in z if c >= 10])
