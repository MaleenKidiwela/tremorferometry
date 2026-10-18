#!/usr/bin/env python
"""BURST MUTING of the 40 sps cache before correlation (2026-09-26, user-approved analyst-style cleaning).

Why: single-station repeating bursts / ringing (W060, W020 ~15 s wave packets, C04A/DOSE bursts, the
B013/W030/B943 site ringing) seed noise families; 8 noise catch-alls absorbed 2,137 of the 2,695
starting NW templates.  Removing stations only moved the noise (notes/2026-09-25_Notes.md sec. U).
An analyst would cut the bursts out of the record and keep the station.  This does that.

RULE (per station-day; all numbers chosen from data, see results/ab_royer/NW_mute/mute/):
  1. energy e[s,c,k] of every channel c in 2 s windows k on a common absolute grid (midnight + 2k s);
     a window is usable only if every sample is valid in the cache mask.
  2. per-channel ratio r[s,c,k] = e / median_k(e)  (the channel's own daily median).
     station ratio R[s,k] = max_c r[s,c,k]  (a burst on any component).
  3. array ratio  Q[s,k] = median over the OTHER live stations s' of  max_{|k'-k|<=MOVE} R[s',k'] /
     median_k(max_{|k'-k|<=MOVE} R[s',k'])  -- i.e. how elevated the rest of the subarray is, allowing
     +-MOVE windows of moveout, normalised so a quiet array reads ~1.
  4. window k of station s is MUTED when R[s,k] > K  AND  Q[s,k] < Q_MAX (the array is quiet):
     a single-station transient.  Array-coherent signals (earthquakes, tremor, LFEs) raise Q and
     are kept.  Mutes are padded by PAD seconds on both sides.
  5. muted samples are handed to every stage EXACTLY like gap samples: load_day() ANDs the mute
     into the validity mask `m` (and sets the samples to 0, as a gap is; never zero without masking).
     Everything downstream (window validity, gap guard, channel refs, cut(), restack) already
     honours `m`.  The mask is per STATION (the cache mask has no channel axis): a burst on any
     component mutes the station for that span.
Spectral lines: see notes in results/ab_royer/AB_ROYER_LOG.md (line scan by --lines).

Files: <mutedir>/energy/<date>.npz  e (S, 3, NK) float32, stas
       <mutedir>/mask/<NET.STA>/<date>.npz  mute (NK,) bool on the 2 s grid (True = muted), params
       <mutedir>/params.json
Consumed by ab_common.load_day when env AB_MUTE_DIR=<mutedir> is set.

  python scripts/ab_mute.py energy --sub NW --out results/ab_royer/NW_mute/mute --workers 24
  python scripts/ab_mute.py mask   --sub NW --out results/ab_royer/NW_mute/mute --K 20 --qmax 2
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

WIN_S = 2.0
NW_S = int(round(WIN_S * A.FS))      # 80 samples
NK = int(86400 / WIN_S)              # 43200 windows per day
MOVE = 3                             # +-3 windows = +-6 s moveout tolerance (neighbour sites <~35 km)
PAD_S = 3.0
LONG_S = 60.0                        # slow array ratio window


def _raw_day(sta, date):
    """load_day WITHOUT the mute (the energy must see the bursts)."""
    old = os.environ.pop("AB_MUTE_DIR", None)
    try:
        return A.load_day(sta, date)
    finally:
        if old is not None:
            os.environ["AB_MUTE_DIR"] = old


def station_energy(sta, date):
    """(3, NK) float32 window energies on the midnight grid, NaN where not fully valid."""
    v = _raw_day(sta, date)
    out = np.full((3, NK), np.nan, np.float32)
    if v is None:
        return out
    t0, X, m = v
    mid = pd.Timestamp(date).timestamp()
    i0 = int(round((mid - t0) * A.FS))              # sample of midnight in X
    N = X.shape[1]
    ks = np.arange(NK)
    s0 = i0 + ks * NW_S
    ok = (s0 >= 0) & (s0 + NW_S <= N)
    kk = ks[ok]
    a, b = s0[ok][0], s0[ok][-1] + NW_S
    Xs = X[:, a:b].astype(np.float64)
    ms = m[a:b]
    Xs = Xs.reshape(3, -1, NW_S)
    e = (Xs * Xs).sum(-1)
    v_ok = ms.reshape(-1, NW_S).all(-1)
    e[:, ~v_ok] = np.nan
    e[e <= 0] = np.nan
    out[:, kk] = e.astype(np.float32)
    return out


def _energy_day(args):
    stas, date, out = args
    p = f"{out}/energy/{date}.npz"
    if os.path.exists(p):
        return date
    E = np.stack([station_energy(s, date) for s in stas])
    np.savez_compressed(p, e=E, stas=np.array(stas))
    return date


def line_scan(sta, date, nper=4096):
    """Hourly Welch spectra -> median over hours of PSD / running-median background (prominence)."""
    from scipy.signal import welch
    from scipy.ndimage import median_filter
    v = _raw_day(sta, date)
    if v is None:
        return None
    t0, X, m = v
    nh = int(3600 * A.FS)
    P = []
    for h in range(24):
        a = h * nh
        if a + nh > X.shape[1] or not m[a:a + nh].all():
            continue
        f, p = welch(X[:, a:a + nh].astype(np.float64), fs=A.FS, nperseg=nper, axis=-1)
        bg = median_filter(p, size=(1, 51), mode="nearest")      # +-0.25 Hz background
        P.append(p / np.maximum(bg, 1e-300))
    if len(P) < 12:
        return None
    return f, np.median(np.stack(P), 0)                         # (3, nf) daily median prominence


LINE_P = 4.0          # daily-median (over hours) prominence above the +-0.25 Hz running median
LINE_NARROW_HZ = 0.05 # the prominence must fall below half the peak within +-0.05 Hz (a LINE, not a bump)
LINE_BAND = (1.0, 8.0)
NOTCH_BW_HZ = 0.06    # -3 dB notch width


def find_lines(f, P):
    """[(channel, f0, prominence)] narrow persistent spectral lines in the 1-8 Hz band."""
    df = f[1] - f[0]
    w = max(int(round(LINE_NARROW_HZ / df)), 1)
    out = []
    for c in range(P.shape[0]):
        p = P[c]
        idx = np.flatnonzero((p >= LINE_P) & (f >= LINE_BAND[0]) & (f <= LINE_BAND[1]))
        for i in idx:
            if i < w or i + w >= p.size:
                continue
            if p[i] < p[i - 1] or p[i] < p[i + 1]:
                continue
            if p[i - w] < p[i] / 2 and p[i + w] < p[i] / 2:
                out.append((c, float(f[i]), float(p[i])))
    # merge peaks closer than the notch width on one channel (keep the strongest)
    out.sort(key=lambda r: -r[2])
    keep = []
    for r in out:
        if not any(k[0] == r[0] and abs(k[1] - r[1]) < NOTCH_BW_HZ for k in keep):
            keep.append(r)
    return keep


def _lines_day(args):
    sta, date, out = args
    p = f"{out}/lines/{A.netsta(sta)}/{date}.npz"
    if os.path.exists(p):
        return sta, date, np.load(p)["lines"]
    r = line_scan(sta, date)
    L = np.array(find_lines(*r), float).reshape(-1, 3) if r is not None else np.zeros((0, 3))
    os.makedirs(os.path.dirname(p), exist_ok=True)
    np.savez_compressed(p, lines=L)
    return sta, date, L


def notch(X, m, lines):
    """Zero-phase notch of every (channel, f0) line in `lines` (n, 3), run by run over the valid
    samples of m (gaps are not filtered through).  Returns a new float32 array."""
    from scipy.signal import iirnotch, sosfiltfilt, tf2sos
    X = X.copy()
    runs = A.valid_runs(m, int(60 * A.FS))
    for c in range(X.shape[0]):
        fl = [f0 for cc, f0, _ in lines if int(cc) == c]
        if not fl:
            continue
        sos = np.vstack([tf2sos(*iirnotch(f0, f0 / NOTCH_BW_HZ, fs=A.FS)) for f0 in fl])
        for a, b in runs:
            X[c, a:b] = sosfiltfilt(sos, X[c, a:b].astype(np.float64)).astype(np.float32)
    return X


# ------------------------------------------------------------------ mask --------------------------
def _maxfilt(x, w):
    from scipy.ndimage import maximum_filter1d
    y = np.where(np.isfinite(x), x, -np.inf)
    y = maximum_filter1d(y, 2 * w + 1, axis=-1, mode="nearest")
    return np.where(np.isfinite(y), y, np.nan)


SITE_KM = 3.0        # stations closer than this share a site (site noise is not array confirmation)
N_NB_SITES = 5       # the array test uses the nearest 5 OTHER sites


def sites(stas):
    """[site index per station]: single-linkage groups of stations < SITE_KM apart."""
    lab = list(range(len(stas)))
    for i, a in enumerate(stas):
        for j, b in enumerate(stas[:i]):
            if A.epi_km(*A.COORD[a][:2], *A.COORD[b][:2]) < SITE_KM:
                old, new = lab[i], lab[j]
                lab = [new if x == old else x for x in lab]
    u = sorted(set(lab))
    return np.array([u.index(x) for x in lab])


def neighbour_sites(stas):
    """{station index: [N_NB_SITES nearest other site ids]} by the distance to each site's closest member."""
    sid = sites(stas)
    nb = {}
    for i, a in enumerate(stas):
        dist = {}
        for j, b in enumerate(stas):
            if sid[j] == sid[i]:
                continue
            d = A.epi_km(*A.COORD[a][:2], *A.COORD[b][:2])
            dist[sid[j]] = min(dist.get(sid[j], 1e9), d)
        nb[i] = [k for k, _ in sorted(dist.items(), key=lambda kv: kv[1])[:N_NB_SITES]]
    return sid, nb


def ratios(E, stas):
    """E (S, 3, NK) -> R (S, NK) station ratio (max over components of e / the channel's daily
    median), Q (S, NK) array ratio: median over the N_NB_SITES nearest OTHER sites of the site value
    (median over the site's live stations of max_{|dk|<=MOVE} R, normalised by its own daily median)."""
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        E = E.astype(np.float64)
        med = np.nanmedian(E, axis=2, keepdims=True)
        r = E / med
        R = np.nanmax(np.where(np.isnan(r), -np.inf, r), axis=1)
        R[~np.isfinite(R)] = np.nan
        Mx = _maxfilt(R, MOVE)
        Mn = Mx / np.nanmedian(Mx, axis=1, keepdims=True)
        sid, nb = neighbour_sites(list(stas))
        SV = np.stack([np.nanmedian(Mn[sid == k], axis=0) if (sid == k).sum() > 1 else Mn[sid == k][0]
                       for k in range(sid.max() + 1)])
        Q = np.full_like(R, np.nan)
        for i in range(R.shape[0]):
            Q[i] = np.nanmedian(SV[nb[i]], axis=0)
        # slow array ratio: LONG_S running mean of each station's R, normalised by its daily median,
        # site median, median over the neighbour sites (tremor raises neighbours for minutes)
        from scipy.ndimage import uniform_filter1d
        w = int(round(LONG_S / WIN_S))
        Rf = np.where(np.isfinite(R), R, 0.0)
        cnt = uniform_filter1d(np.isfinite(R).astype(float), w, axis=1, mode="nearest")
        L = uniform_filter1d(Rf, w, axis=1, mode="nearest") / np.where(cnt > 0.5, cnt, np.nan)
        L = L / np.nanmedian(L, axis=1, keepdims=True)
        SL = np.stack([np.nanmedian(L[sid == k], axis=0) if (sid == k).sum() > 1 else L[sid == k][0]
                       for k in range(sid.max() + 1)])
        QL = np.full_like(R, np.nan)
        for i in range(R.shape[0]):
            QL[i] = np.nanmedian(SL[nb[i]], axis=0)
    return R, Q, QL


def mute_from(R, Q, K, qmax, pad_s=PAD_S, QL=None, qlong=np.inf):
    mu = (R > K) & (Q < qmax)
    if QL is not None and np.isfinite(qlong):
        mu &= (QL < qlong)
    pw = int(np.ceil(pad_s / WIN_S))
    if pw:
        from scipy.ndimage import maximum_filter1d
        mu = maximum_filter1d(mu.astype(np.uint8), 2 * pw + 1, axis=-1).astype(bool)
    return mu


def _mask_day(args):
    date, out, K, qmax, qlong, pad, use_lines = args
    z = np.load(f"{out}/energy/{date}.npz")
    E, stas = z["e"], list(z["stas"])
    R, Q, QL = ratios(E, stas)
    mu = mute_from(R, Q, K, qmax, pad, QL, qlong)
    rows = []
    for i, s in enumerate(stas):
        d = f"{out}/mask/{A.netsta(s)}"
        os.makedirs(d, exist_ok=True)
        lp = f"{out}/lines/{A.netsta(s)}/{date}.npz"
        L = np.load(lp)["lines"] if (use_lines and os.path.exists(lp)) else np.zeros((0, 3))
        np.savez_compressed(f"{d}/{date}.npz", mute=mu[i], K=K, qmax=qmax, qlong=qlong, pad=pad, win=WIN_S,
                            lines=L)
        live = np.isfinite(R[i])
        rows.append(dict(date=date, sta=s, live_win=int(live.sum()),
                         muted_win=int((mu[i] & live).sum()),
                         frac=float((mu[i] & live).sum() / max(live.sum(), 1))))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["energy", "mask", "lines"])
    ap.add_argument("--no-lines", action="store_true", help="mask files carry no notch lines")
    ap.add_argument("--sub", default="NW")
    ap.add_argument("--out", required=True)
    ap.add_argument("--dates", default="")
    ap.add_argument("--pool", default="")
    ap.add_argument("--K", type=float, default=20.0)
    ap.add_argument("--qmax", type=float, default=2.0)
    ap.add_argument("--qlong", type=float, default=np.inf)
    ap.add_argument("--pad", type=float, default=PAD_S)
    ap.add_argument("--workers", type=int, default=24)
    a = ap.parse_args()
    stas = A.pool_from_arg(a.pool, a.sub) if a.pool else A.POOL[a.sub]
    dates = a.dates.split(",") if a.dates else A.scan_days(a.sub)
    os.makedirs(f"{a.out}/energy", exist_ok=True)
    if a.mode == "energy":
        with ProcessPoolExecutor(a.workers) as ex:
            for i, d in enumerate(ex.map(_energy_day, [(stas, d, a.out) for d in dates])):
                if i % 20 == 0:
                    A.log(f"[mute energy] {i + 1}/{len(dates)} {d}")
    elif a.mode == "lines":
        rows = []
        with ProcessPoolExecutor(a.workers) as ex:
            for sta, d, L in ex.map(_lines_day, [(s, d, a.out) for s in stas for d in dates]):
                for c, f0, p in L:
                    rows.append(dict(sta=sta, date=d, ch="ZNE"[int(c)], f0=f0, prom=p))
        df = pd.DataFrame(rows)
        df.to_csv(f"{a.out}/lines_all.csv", index=False)
        if len(df):
            df["f_r"] = df.f0.round(1)
            print(df.groupby(["sta", "f_r"]).size().sort_values(ascending=False).head(40).to_string())
    elif a.mode == "mask":
        rows = []
        with ProcessPoolExecutor(a.workers) as ex:
            for r in ex.map(_mask_day, [(d, a.out, a.K, a.qmax, a.qlong, a.pad, not a.no_lines) for d in dates]):
                rows += r
        df = pd.DataFrame(rows)
        df.to_csv(f"{a.out}/muted_fraction_by_day.csv", index=False)
        g = df.groupby("sta").agg(live_win=("live_win", "sum"), muted_win=("muted_win", "sum"))
        g["frac"] = g.muted_win / g.live_win
        g.sort_values("frac", ascending=False).to_csv(f"{a.out}/muted_fraction_by_station.csv")
        json.dump(dict(K=a.K, qmax=a.qmax, qlong=a.qlong, long_s=LONG_S, pad_s=a.pad, win_s=WIN_S, move_win=MOVE, stas=stas,
                       n_dates=len(dates)), open(f"{a.out}/params.json", "w"), indent=1)
        A.log(f"[mute mask] K={a.K} qmax={a.qmax}: overall muted {g.muted_win.sum() / g.live_win.sum():.4f}")
        print(g.sort_values("frac", ascending=False).to_string())


if __name__ == "__main__":
    main()
