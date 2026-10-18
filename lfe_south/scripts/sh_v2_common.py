#!/usr/bin/env python
"""SH_V2 (2026-09-28) -- shared cleaning for the Shelly (2025) Mendocino detector v2.

Ports the Washington fixes (notes/2026-09-25_Notes.md sec. Q and sec. V) to the lfe_south 20 Hz cache
(`cache20hz/<NET.STA>/<yyyy>/<doy>.npz`, x (C, N) float, NO validity mask, gaps zero-filled).

WHY (results/SH_V2_LOG.md step 1): 664 of 10,395 row-C family-days (6.4 %) carried 98.5 % of all
detections.  Every one sits on a day where a DEAD channel (BK.BRIC HHE, BK.RVIT HHE, NC.KCR/KSM EHZ:
median |x| 1e-17..3e-13, float residue, never exactly constant) also carries a multi-million-count
glitch (peak 3.9e6-7.5e6).  The day-long cumsum S2 then has absolute round-off ~1e-3 >> the residue
energy, the AA guard (std > 1e-4 x median std) passes, and |cc| reaches 31-211.

WHAT (all per station-day, applied identically in the matched filter AND in every template stack):
  1. flat samples (exactly zero or repeating the previous sample) on the RAW cache are invalid, padded
     by FLAT_PAD_S (the 2 Hz high-pass transient a gap edge injects); a sample is valid in `m`.
  2. 2 Hz high-pass -> Shelly's 2-9 Hz band (cache is 1.5-9 Hz).
  3. BURST MUTING (lfe_array_discovery/scripts/ab_mute.py, NW_mute parameters K=15, qmax=2,
     qlong=1.3, pad 3 s, 2 s windows, +-3 windows moveout, 5 nearest other sites, sites < 3 km):
     a 2 s window of a station is muted when its energy is > K x that channel's daily median (on any
     component) while the neighbour sites are quiet (2 s Q < 2 and 60 s Q < 1.3).  The daily median
     is taken over ALL positive windows, so a dead channel's median is its residue and a glitch on it
     is muted (the WA 'truly non-zero' median would have made the glitch the channel's own level).
     Muted samples are exactly like gap samples: m False, x 0.
  4. LINE NOTCHING (ab_mute.find_lines / notch): hourly Welch, daily median prominence >= 4 over a
     +-0.25 Hz running median, narrow within +-0.05 Hz, band 2-9 Hz here; zero-phase 0.06 Hz notch
     over each valid run >= 60 s.
  5. EXACT ENERGIES + RELATIVE GUARD + DEAD FLOORS (ab_common.py window_stats_guarded / channel_refs /
     day_refs, GUARD_REL 1e-3, GUARD_STA_REL 1e-4, GUARD_MAX_FLAT 0.01, GUARD_NET_REL 1e-6,
     GUARD_MIN_LIVE_S 600): block-local window sums (no day-long cumsum); a channel-window is valid
     only if every sample is valid, <= 1 % of samples are flat, the channel has a reference energy
     (ref_c > 0: >= 10 min of truly non-zero windows) and e > 1e-3 ref_c; a station whose loudest
     ref is < 1e-6 x the day's median station level is dead.
"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t_common import *          # ROOT, WS, RES, CACHE, FS, npz_path, log

V2 = f"{RES}/sh_v2"
PREP = os.environ.get("SH_PREP_DIR", f"{V2}/prep")   # v3b (2026-09-29): overridable, default unchanged
HP_HZ = 2.0
FLAT_PAD_S = 2.0
# ---- ab_mute.py (NW_mute run, notes sec. V) --------------------------------------------------------
WIN_S = 2.0; NW_S = int(round(WIN_S * FS)); NK = int(86400 / WIN_S)
MOVE = 3; PAD_S = 3.0; LONG_S = 60.0; SITE_KM = 3.0; N_NB_SITES = 5
MUTE_K = 15.0; MUTE_QMAX = 2.0; MUTE_QLONG = 1.3
LINE_P = 4.0; LINE_NARROW_HZ = 0.05; LINE_BAND = (2.0, 9.0); NOTCH_BW_HZ = 0.06
LINE_NPER = 2048                     # 20 sps: same 0.0098 Hz resolution as ab_mute's 4096 @ 40 sps
# ---- ab_common.py gap guard (notes sec. Q) ---------------------------------------------------------
GUARD_REL = 1e-3; GUARD_STA_REL = 1e-4; GUARD_MAX_FLAT = 0.01; GUARD_NET_REL = 1e-6
GUARD_MIN_LIVE_S = 600.0

_SOS = None
def sos_hp():
    global _SOS
    if _SOS is None:
        from scipy.signal import butter
        _SOS = butter(4, HP_HZ / (FS / 2), btype="highpass", output="sos")
    return _SOS


def coords():
    f = pd.read_csv(f"{RES}/fdsn_channels_south.csv", dtype=str)
    f["ns"] = f.net + "." + f.sta
    g = f.groupby("ns").first()
    return {ns: (float(r.lat), float(r.lon)) for ns, r in g.iterrows()}


def epi_km(la1, lo1, la2, lo2):
    R = 6371.0
    p1, p2 = np.radians(la1), np.radians(la2)
    dp = p2 - p1; dl = np.radians(lo2 - lo1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def load_raw(ns, d):
    net, sta = ns.split(".")
    p = npz_path(net, sta, d)
    if not os.path.exists(p):
        return None
    try:
        z = np.load(p); x = z["x"].astype(np.float64); t0 = float(z["t0"])
        ch = [str(c) for c in z["chans"]] if "chans" in z.files else []
    except Exception:
        return None
    if x.ndim != 2 or x.shape[1] < 0.5 * 86400 * FS or not np.any(x):
        return None
    return t0, x, ch


def flat_mask(x):
    """(C, N) bool: sample exactly zero or repeating the previous sample (zero-filled gap)."""
    fl = (x == 0)
    fl[:, 1:] |= (x[:, 1:] == x[:, :-1])
    return fl


def _pad_bool(b, w):
    if w <= 0 or not b.any():
        return b
    from scipy.ndimage import maximum_filter1d
    return maximum_filter1d(b.astype(np.uint8), 2 * w + 1).astype(bool)


def base_clean(x):
    """raw (C, N) -> (X hp2 filtered with invalid samples zeroed, mc (C, N) per-channel validity).
    Flat runs >= 0.5 s are invalid and padded by FLAT_PAD_S (filter transient); isolated repeats
    (a zero crossing hitting 0.0 exactly) only mark that sample."""
    from scipy.signal import sosfiltfilt
    fl = flat_mask(x)
    C, N = x.shape
    mc = np.ones((C, N), bool)
    runmin = int(0.5 * FS); w = int(FLAT_PAD_S * FS)
    for c in range(C):
        f = fl[c]
        if not f.any():
            continue
        # runs of flat samples
        d = np.diff(np.concatenate([[0], f.astype(np.int8), [0]]))
        a = np.flatnonzero(d == 1); b = np.flatnonzero(d == -1)
        long_ = np.zeros(N, bool)
        for i, j in zip(a, b):
            if j - i >= runmin:
                long_[i:j] = True
        mc[c] = ~_pad_bool(long_, w)          # isolated exact repeats stay valid (WA: <=1 % flat)
    X = sosfiltfilt(sos_hp(), np.where(fl, 0.0, x), axis=-1)
    X = np.where(mc, X, 0.0)
    return X, mc


# ------------------------------------------------------------------ muting (ab_mute port) ----------
def energy_2s(X, mc, t0, d):
    """(C, NK) 2 s window energies on the midnight grid, NaN unless every sample valid."""
    C, N = X.shape
    out = np.full((C, NK), np.nan)
    mid = pd.Timestamp(d).timestamp()
    i0 = int(round((mid - t0) * FS))
    ks = np.arange(NK); s0 = i0 + ks * NW_S
    ok = (s0 >= 0) & (s0 + NW_S <= N)
    kk = ks[ok]
    if not len(kk):
        return out
    a, b = s0[ok][0], s0[ok][-1] + NW_S
    Xs = X[:, a:b].reshape(C, -1, NW_S)
    e = (Xs * Xs).sum(-1)
    v = mc[:, a:b].reshape(C, -1, NW_S).all(-1)
    e[~v] = np.nan; e[e <= 0] = np.nan
    out[:, kk] = e
    return out


def _maxfilt(x, w):
    from scipy.ndimage import maximum_filter1d
    y = np.where(np.isfinite(x), x, -np.inf)
    y = maximum_filter1d(y, 2 * w + 1, axis=-1, mode="nearest")
    return np.where(np.isfinite(y), y, np.nan)


def neighbour_sites(stas, C):
    n = len(stas); lab = list(range(n))
    for i in range(n):
        for j in range(i):
            if epi_km(*C[stas[i]], *C[stas[j]]) < SITE_KM:
                old, new = lab[i], lab[j]; lab = [new if x == old else x for x in lab]
    u = sorted(set(lab)); sid = np.array([u.index(x) for x in lab])
    nb = {}
    for i in range(n):
        dist = {}
        for j in range(n):
            if sid[j] == sid[i]:
                continue
            dd = epi_km(*C[stas[i]], *C[stas[j]]); dist[sid[j]] = min(dist.get(sid[j], 1e9), dd)
        nb[i] = [k for k, _ in sorted(dist.items(), key=lambda kv: kv[1])[:N_NB_SITES]]
    return sid, nb


def ratios(Elist, stas, C):
    """Elist [ (C_s, NK) ] -> R (S, NK) station ratio, Q (S, NK) 2 s array ratio, QL 60 s array ratio
    (ab_mute.ratios with a per-station channel count)."""
    from scipy.ndimage import uniform_filter1d
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        R = []
        for E in Elist:
            med = np.nanmedian(E, axis=1, keepdims=True)
            r = E / med
            Rs = np.nanmax(np.where(np.isnan(r), -np.inf, r), axis=0)
            Rs[~np.isfinite(Rs)] = np.nan
            R.append(Rs)
        R = np.stack(R)
        Mx = _maxfilt(R, MOVE)
        Mn = Mx / np.nanmedian(Mx, axis=1, keepdims=True)
        sid, nb = neighbour_sites(stas, C)
        SV = np.stack([np.nanmedian(Mn[sid == k], axis=0) if (sid == k).sum() > 1 else Mn[sid == k][0]
                       for k in range(sid.max() + 1)])
        Q = np.full_like(R, np.nan)
        for i in range(R.shape[0]):
            Q[i] = np.nanmedian(SV[nb[i]], axis=0) if nb[i] else np.nan
        w = int(round(LONG_S / WIN_S))
        Rf = np.where(np.isfinite(R), R, 0.0)
        cnt = uniform_filter1d(np.isfinite(R).astype(float), w, axis=1, mode="nearest")
        L = uniform_filter1d(Rf, w, axis=1, mode="nearest") / np.where(cnt > 0.5, cnt, np.nan)
        L = L / np.nanmedian(L, axis=1, keepdims=True)
        SL = np.stack([np.nanmedian(L[sid == k], axis=0) if (sid == k).sum() > 1 else L[sid == k][0]
                       for k in range(sid.max() + 1)])
        QL = np.full_like(R, np.nan)
        for i in range(R.shape[0]):
            QL[i] = np.nanmedian(SL[nb[i]], axis=0) if nb[i] else np.nan
    return R, Q, QL


def mute_from(R, Q, QL):
    with np.errstate(invalid="ignore"):
        mu = (R > MUTE_K) & (Q < MUTE_QMAX) & (QL < MUTE_QLONG)
    pw = int(np.ceil(PAD_S / WIN_S))
    from scipy.ndimage import maximum_filter1d
    return maximum_filter1d(mu.astype(np.uint8), 2 * pw + 1, axis=-1).astype(bool)


def mute_samples(mu, t0, d, N):
    """(NK,) 2 s mute -> (N,) sample mute."""
    ms = np.zeros(N, bool)
    if not mu.any():
        return ms
    i0 = int(round((pd.Timestamp(d).timestamp() - t0) * FS))
    k = np.flatnonzero(mu)
    a = np.clip(i0 + k * NW_S, 0, N); b = np.clip(i0 + (k + 1) * NW_S, 0, N)
    dd = np.zeros(N + 1, np.int32); np.add.at(dd, a, 1); np.add.at(dd, b, -1)
    return np.cumsum(dd[:-1]) > 0


# ------------------------------------------------------------------ lines (ab_mute port) -----------
def line_scan(X, mc):
    from scipy.signal import welch
    from scipy.ndimage import median_filter
    nh = int(3600 * FS); P = []
    ok = mc.all(0)
    for h in range(24):
        a = h * nh
        if a + nh > X.shape[1] or not ok[a:a + nh].all():
            continue
        f, p = welch(X[:, a:a + nh], fs=FS, nperseg=LINE_NPER, axis=-1)
        bg = median_filter(p, size=(1, 51), mode="nearest")
        P.append(p / np.maximum(bg, 1e-300))
    if len(P) < 12:
        return None
    return f, np.median(np.stack(P), 0)


def find_lines(f, P):
    df = f[1] - f[0]; w = max(int(round(LINE_NARROW_HZ / df)), 1); out = []
    for c in range(P.shape[0]):
        p = P[c]
        idx = np.flatnonzero((p >= LINE_P) & (f >= LINE_BAND[0]) & (f <= LINE_BAND[1]))
        for i in idx:
            if i < w or i + w >= p.size or p[i] < p[i - 1] or p[i] < p[i + 1]:
                continue
            if p[i - w] < p[i] / 2 and p[i + w] < p[i] / 2:
                out.append((c, float(f[i]), float(p[i])))
    out.sort(key=lambda r: -r[2]); keep = []
    for r in out:
        if not any(k[0] == r[0] and abs(k[1] - r[1]) < NOTCH_BW_HZ for k in keep):
            keep.append(r)
    return keep


def valid_runs(m, min_len):
    mm = np.asarray(m, bool).astype(np.int8)
    d = np.diff(np.concatenate([[0], mm, [0]]))
    a = np.flatnonzero(d == 1); b = np.flatnonzero(d == -1)
    return [(int(i), int(j)) for i, j in zip(a, b) if j - i >= min_len]


def notch(X, mc, lines):
    from scipy.signal import iirnotch, sosfiltfilt, tf2sos
    if not len(lines):
        return X
    X = X.copy()
    for c in range(X.shape[0]):
        fl = [f0 for cc, f0, _ in lines if int(cc) == c]
        if not fl:
            continue
        sos = np.vstack([tf2sos(*iirnotch(f0, f0 / NOTCH_BW_HZ, fs=FS)) for f0 in fl])
        for a, b in valid_runs(mc[c], int(60 * FS)):
            X[c, a:b] = sosfiltfilt(sos, X[c, a:b])
    return X


# ------------------------------------------------------------------ exact energies (ab_common) ----
def box_sum(v, n, starts, xp=np):
    """exact sliding sums of v (C, N) over [s, s+n) for s in starts: block-local prefix + suffix."""
    C, N = v.shape
    K = -(-N // n) + 1
    pad = xp.zeros((C, K * n))
    pad[:, :N] = v
    B = pad.reshape(C, K, n)
    pre = xp.concatenate([xp.zeros((C, K, 1)), xp.cumsum(B, axis=2)], 2)
    suf = xp.concatenate([xp.cumsum(B[:, :, ::-1], axis=2)[:, :, ::-1], xp.zeros((C, K, 1))], 2)
    k = starts // n; r = starts % n
    return suf[:, k, r] + pre[:, k + 1, r]


def channel_refs(X, mc, n, sub=20):
    """ab_common.channel_refs with per-channel validity: (C,) reference window energy, 0 = dead."""
    C, N = X.shape
    if N < n:
        return np.zeros(C)
    starts = np.arange(0, N - n + 1, sub)
    e = box_sum(X * X, n, starts)
    inv = (~mc).astype(np.int32)
    cf = np.concatenate([np.zeros((C, 1), np.int64), np.cumsum(inv, 1)], 1)
    live = ((cf[:, starts + n] - cf[:, starts]) / n) <= GUARD_MAX_FLAT
    min_live = max(int(GUARD_MIN_LIVE_S * FS / sub), 1)
    med = np.array([float(np.median(e[c, live[c]])) if live[c].sum() >= min_live else 0.0
                    for c in range(C)])
    top = float(med.max()) if med.size else 0.0
    ref = np.zeros(C)
    if top <= 0:
        return ref
    for c in range(C):
        good = live[c] & (e[c] > GUARD_STA_REL * top)
        if good.sum() >= min_live:
            ref[c] = float(np.median(e[c, good]))
    return ref


# ------------------------------------------------------------------ day prep + load ----------------
def prep_day(d, stas, C, nref):
    """Compute mute (per station, 2 s grid), lines and channel refs (window nref samples) for every
    station on date d; write PREP/<date>.npz.  Returns a summary row list."""
    d = pd.Timestamp(d)
    out = f"{PREP}/{d.date()}.npz"
    if os.path.exists(out):
        return []
    RAW = {}
    for ns in stas:
        v = load_raw(ns, d)
        if v is not None:
            RAW[ns] = v
    S = [s for s in RAW if s in C]
    rows = []; save = {}
    if not S:
        np.savez_compressed(out, stas=np.array([], dtype=str)); return rows
    CL = {}
    for ns in S:
        t0, x, ch = RAW[ns]
        X, mc = base_clean(x)
        CL[ns] = (t0, X, mc, ch)
    E = [energy_2s(CL[s][1], CL[s][2], CL[s][0], d) for s in S]
    # DEAD CHANNELS FIRST (SH_V2 step 1): a channel whose PLAIN daily median 2 s energy is below
    # GUARD_STA_REL x its station's loudest channel median is residue (BK.BRIC HHE: 1e-26 vs 1e4)
    # -- dead for the whole day, whatever glitch it carries; excluded from the burst ratio R (its
    # residue fluctuations would otherwise read as bursts) and given ref = 0.  Same for a station
    # whose loudest plain median is < GUARD_NET_REL x the day's median station level.
    with np.errstate(invalid="ignore"):
        PM = [np.nanmedian(e, axis=1) for e in E]
    DEADC = {}
    for s, pm in zip(S, PM):
        pm = np.where(np.isfinite(pm), pm, 0.0); top = pm.max()
        DEADC[s] = (pm <= GUARD_STA_REL * top) if top > 0 else np.ones(pm.size, bool)
    stop = np.array([np.nanmax(np.where(np.isfinite(pm), pm, 0.0)) for pm in PM])
    nfloor = GUARD_NET_REL * float(np.median(stop[stop > 0])) if (stop > 0).any() else np.inf
    for s, t_ in zip(S, stop):
        if t_ < nfloor: DEADC[s][:] = True
    E = [np.where(DEADC[s][:, None], np.nan, e) for s, e in zip(S, E)]
    if len(S) >= 2:
        R, Q, QL = ratios(E, S, C)
        MU = mute_from(R, Q, QL)
    else:
        MU = np.zeros((len(S), NK), bool); R = np.full((len(S), NK), np.nan)
    tops = {}
    REF = {}
    for i, ns in enumerate(S):
        t0, X, mc, ch = CL[ns]
        ms = mute_samples(MU[i], t0, d, X.shape[1])
        mc2 = mc & ~ms[None, :]
        X2 = np.where(mc2, X, 0.0)
        r = line_scan(X2, mc2)
        L = np.array(find_lines(*r), float).reshape(-1, 3) if r is not None else np.zeros((0, 3))
        X3 = notch(X2, mc2, L)
        ref = channel_refs(X3, mc2, nref)
        ref[DEADC[ns]] = 0.0
        REF[ns] = ref; tops[ns] = ref.max()
        save[f"{ns}|mute"] = MU[i]; save[f"{ns}|lines"] = L
        live = np.isfinite(R[i])
        rows.append(dict(dead_plain=int(DEADC[ns].sum()), **dict(day=str(d.date()), netsta=ns, nch=X.shape[0], chans=",".join(ch),
                         invalid_frac=float((~mc).mean()), muted_frac=float((MU[i] & live).sum() / max(live.sum(), 1)),
                         n_lines=int(len(L)), lines=";".join(f"{'%d' % c}:{f0:.2f}" for c, f0, _ in L),
                         ref=",".join(f"{v:.3g}" for v in ref))))
    tp = np.array([v for v in tops.values() if v > 0])
    floor = GUARD_NET_REL * float(np.median(tp)) if tp.size else np.inf
    for ns in S:
        if REF[ns].max() < floor:
            REF[ns] = np.zeros_like(REF[ns])
        save[f"{ns}|ref"] = REF[ns]
    for r_ in rows:
        r_["station_dead"] = bool(REF[r_["netsta"]].max() <= 0)
        r_["dead_ch"] = int((REF[r_["netsta"]] <= 0).sum())
    np.savez_compressed(out, stas=np.array(S), nref=nref, **save)
    return rows


class DayPrep:
    def __init__(self, d):
        self.d = pd.Timestamp(d)
        p = f"{PREP}/{self.d.date()}.npz"
        self.z = dict(np.load(p, allow_pickle=False)) if os.path.exists(p) else None
        self.stas = set(str(s) for s in self.z["stas"]) if self.z is not None else set()

    def load(self, ns):
        """(t0, X (C,N) cleaned+muted+notched, mc (C,N) validity, ref (C,)) or None."""
        if ns not in self.stas:
            return None
        v = load_raw(ns, self.d)
        if v is None:
            return None
        t0, x, ch = v
        X, mc = base_clean(x)
        ms = mute_samples(self.z[f"{ns}|mute"].astype(bool), t0, self.d, X.shape[1])
        mc = mc & ~ms[None, :]
        X = np.where(mc, X, 0.0)
        X = notch(X, mc, self.z[f"{ns}|lines"])
        return t0, X, mc, self.z[f"{ns}|ref"]


def window_ok(X, mc, ref, n, starts):
    """(C, len(starts)) bool: the gap guard for windows [s, s+n) (numpy, used by the stackers)."""
    C, N = X.shape
    e = box_sum(X * X, n, starts)
    inv = (~mc).astype(np.int32)
    cf = np.concatenate([np.zeros((C, 1), np.int64), np.cumsum(inv, 1)], 1)
    ninv = cf[:, starts + n] - cf[:, starts]
    return (ninv == 0) & (ref[:, None] > 0) & (e > GUARD_REL * ref[:, None])
