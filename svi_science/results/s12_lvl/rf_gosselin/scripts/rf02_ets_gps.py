#!/usr/bin/env python
"""rf02_ets_gps: ETS start / midpoint / end from GPS east (Gosselin 2020 recipe) for ALBH (SVI) and WA stations.

Recipe (paper, Materials & Methods): east component, detrend, 14-day running mean, cut a window around each
ETS, fit a 5th-order polynomial, start/end = first-derivative ~0 on either side, midpoint = second-derivative 0
(time of maximum westward rate).  Ours: UNR IGS20 tenv3; equipment / earthquake steps from UNR steps.txt are
removed (median 20 d after - 20 d before, ETS-free); linear trend + annual + semi-annual removed robustly; episodes
detected automatically as minima of the 20-day slope of the smoothed series (separated by >= 120 d) then the
polynomial is fitted in +-45 d around each minimum.  Cross-checked against repo tremor lists.
Outputs: tables/ets_gps_<STA>.csv, figures/ets_gps_<STA>.png
"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/jovyan/tremorferometry"
W = f"{ROOT}/svi_science/results/s12_lvl/rf_gosselin"
G = f"{ROOT}/data/rf_teleseis/gps"


def load(sta):
    d = pd.read_csv(f"{G}/{sta}.tenv3", sep=r"\s+")
    t = pd.to_datetime(d["__MJD"].astype(float) + 2400000.5 - 2440587.5, unit="D").dt.normalize()
    e = (d["_e0(m)"].astype(float) + d["__east(m)"].astype(float)) * 1000.0  # mm
    s = pd.Series(e.values, index=t).groupby(level=0).mean()
    s = s.reindex(pd.date_range(s.index.min(), s.index.max(), freq="D"))
    return s


def steps(sta):
    out = []
    for line in open(f"{G}/steps.txt"):
        p = line.split()
        if p and p[0] == sta:
            out.append(pd.to_datetime(p[1], format="%y%b%d"))
    return sorted(set(out))


def clean(s, sta):
    s = s.copy()
    for st in steps(sta):
        a = s[(s.index >= st - pd.Timedelta(days=20)) & (s.index < st)].median()
        b = s[(s.index > st) & (s.index <= st + pd.Timedelta(days=20))].median()
        if np.isfinite(a) and np.isfinite(b):
            s[s.index > st] -= (b - a)
    tt = (s.index - s.index[0]).days.values / 365.25
    ok = np.isfinite(s.values)
    X = np.c_[np.ones_like(tt), tt, np.sin(2 * np.pi * tt), np.cos(2 * np.pi * tt),
              np.sin(4 * np.pi * tt), np.cos(4 * np.pi * tt)]
    w = np.ones(ok.sum())
    for _ in range(5):  # IRLS (Huber-ish) so ETS reversals do not bias the trend
        c = np.linalg.lstsq(X[ok] * w[:, None], s.values[ok] * w, rcond=None)[0]
        r = s.values[ok] - X[ok] @ c
        m = 1.4826 * np.median(np.abs(r))
        w = np.minimum(1, 1.5 * m / np.maximum(np.abs(r), 1e-9))
    # keep only linear + seasonal removal; piecewise-linear secular changes are absorbed by the local windows
    return s - X @ c


def detect(r, min_sep=120, win=20):
    sm = r.rolling(14, center=True, min_periods=7).mean()
    # 20-day slope (mm/day) by centered difference of the smoothed series
    v = (sm.shift(-win // 2) - sm.shift(win // 2)) / win
    vv = v.values; idx = []
    mad = 1.4826 * np.nanmedian(np.abs(vv - np.nanmedian(vv)))
    thr = np.nanmedian(vv) - 2.5 * mad
    order = np.argsort(np.nan_to_num(vv, nan=0))
    for i in order:
        if not np.isfinite(vv[i]) or vv[i] > thr:
            break
        if all(abs(i - j) >= min_sep for j in idx):
            idx.append(i)
    return sm, v, sorted(idx), thr


def poly_times(sm, i0, half=45):
    a, b = max(0, i0 - half), min(len(sm), i0 + half + 1)
    y = sm.values[a:b]; x = np.arange(a, b, dtype=float)
    ok = np.isfinite(y)
    if ok.sum() < 30:
        return None
    xc = x - i0
    c = np.polyfit(xc[ok], y[ok], 5)
    d1 = np.polyder(c, 1); d2 = np.polyder(c, 2)
    xs = np.arange(-half, half + 0.01, 0.1)
    v1 = np.polyval(d1, xs); v2 = np.polyval(d2, xs)
    # midpoint: zero of 2nd derivative nearest to 0 where slope is negative (max westward rate)
    z2 = xs[:-1][(np.sign(v2[:-1]) != np.sign(v2[1:])) & (v1[:-1] < 0)]
    mid = z2[np.argmin(np.abs(z2))] if len(z2) else xs[np.argmin(v1)]
    # start/end: zeros of 1st derivative on either side of mid
    z1 = xs[:-1][np.sign(v1[:-1]) != np.sign(v1[1:])]
    lo = z1[z1 < mid]; hi = z1[z1 > mid]
    st = lo.max() if len(lo) else np.nan
    en = hi.min() if len(hi) else np.nan
    amp = (np.polyval(c, st) - np.polyval(c, en)) if np.isfinite(st) and np.isfinite(en) else np.nan
    return mid, st, en, amp, float(np.polyval(d1, mid))


def run(sta, t0="1993-01-01", t1="2026-12-31"):
    s = load(sta)
    r = clean(s, sta)
    sm, v, idx, thr = detect(r)
    rows = []
    for i in idx:
        res = poly_times(sm, i)
        if res is None:
            continue
        mid, st, en, amp, vmid = res
        day0 = sm.index[i]
        rows.append(dict(sta=sta, detect=day0.date(), mid=(day0 + pd.Timedelta(days=mid)).normalize().date(),
                         start=(day0 + pd.Timedelta(days=st)).date() if np.isfinite(st) else None,
                         end=(day0 + pd.Timedelta(days=en)).date() if np.isfinite(en) else None,
                         west_mm=round(float(amp), 2) if np.isfinite(amp) else np.nan,
                         rate_mid_mm_d=round(vmid, 3), v20_min=round(float(v.values[i]), 3)))
    E = pd.DataFrame(rows).sort_values("mid").reset_index(drop=True)
    E.to_csv(f"{W}/tables/ets_gps_{sta}.csv", index=False)
    fig, ax = plt.subplots(2, 1, figsize=(14, 6), sharex=True)
    ax[0].plot(r.index, r.values, ".", ms=1, color="0.7")
    ax[0].plot(sm.index, sm.values, lw=0.8, color="k")
    ax[1].plot(v.index, v.values, lw=0.7, color="k"); ax[1].axhline(thr, color="r", ls=":")
    for m in pd.to_datetime(E["mid"]):
        for a in ax:
            a.axvline(m, color="tab:blue", lw=0.6, alpha=0.6)
    ax[0].set_ylabel(f"{sta} east, detrended (mm)"); ax[1].set_ylabel("20-day slope (mm/day)")
    ax[0].set_title(f"{sta}: ETS midpoints (blue) from GPS east, Gosselin 2020 recipe; n = {len(E)}")
    fig.tight_layout(); fig.savefig(f"{W}/figures/ets_gps_{sta}.png", dpi=120); plt.close(fig)
    print(sta, len(E)); print(E.to_string())
    return E




# ---------------------------------------------------------------- anchored mode (used downstream)
PRE2003_SVI = ["1994-07-30", "1995-12-24", "1997-04-28", "1998-09-12", "1999-08-24", "2000-12-09", "2002-02-02"]


def anchored(sta, anchors, label):
    """Polynomial (paper) midpoint within +-45 d of each anchor date; anchors = list of (episode, centre)."""
    s = load(sta); r = clean(s, sta)
    sm = r.rolling(14, center=True, min_periods=7).mean()
    rows = []
    for ep, c in anchors:
        c = pd.Timestamp(c)
        if c < sm.index[0] + pd.Timedelta(days=45) or c > sm.index[-1] - pd.Timedelta(days=45):
            rows.append(dict(episode=ep, anchor=c.date(), gps_sta=sta, mid=c.date(), mid_src="anchor(no GPS)"))
            continue
        i0 = int((c - sm.index[0]).days)
        # refine centre to the steepest westward 20-d slope within +-25 d, then fit the polynomial there
        v = (sm.shift(-10) - sm.shift(10)) / 20
        seg = v.iloc[max(0, i0 - 25):i0 + 26]
        i1 = int((seg.idxmin() - sm.index[0]).days) if seg.notna().any() else i0
        res = poly_times(sm, i1)
        if res is None:
            rows.append(dict(episode=ep, anchor=c.date(), gps_sta=sta, mid=c.date(), mid_src="anchor(gap)")); continue
        mid, st, en, amp, vmid = res
        d0 = sm.index[i1]
        ok = abs((d0 + pd.Timedelta(days=mid) - c).days) <= 30 and vmid < 0
        rows.append(dict(episode=ep, anchor=c.date(), gps_sta=sta,
                         mid=((d0 + pd.Timedelta(days=mid)).normalize() if ok else c).date(),
                         mid_src="gps_poly" if ok else "anchor(gps disagrees)",
                         start=(d0 + pd.Timedelta(days=st)).date() if np.isfinite(st) else None,
                         end=(d0 + pd.Timedelta(days=en)).date() if np.isfinite(en) else None,
                         west_mm=round(float(amp), 2) if np.isfinite(amp) else np.nan,
                         rate_mid_mm_d=round(vmid, 3)))
    E = pd.DataFrame(rows)
    E.to_csv(f"{W}/tables/ets_mid_{label}.csv", index=False)
    print(label); print(E.to_string())
    return E


def svi_anchors():
    B = pd.read_csv(f"{ROOT}/svi_science/results/r23/tables/study4/big_ets_2003_2026.csv")
    a = [(f"{d[:7]}", d) for d in PRE2003_SVI]
    for r in B.itertuples():
        c = pd.Timestamp(r.start) + (pd.Timestamp(r.end) - pd.Timestamp(r.start)) / 2
        a.append((r.episode, c.normalize()))
    return a


def wa_anchors():
    B = pd.read_csv(f"{ROOT}/svi_science/results/b_ets_record/episodes_b.csv")
    B = B[(B.region == "NW") & (B.big == True)]
    # merge NW sub-episodes closer than 20 d (e.g. 2019-03/04/05) into one anchor = tremor-weighted span centre
    out = []
    for r in B.sort_values("start").itertuples():
        s, e = pd.Timestamp(r.start), pd.Timestamp(r.end)
        if out and (s - out[-1][2]).days < 20:
            out[-1] = (out[-1][0], out[-1][1], max(e, out[-1][2]))
        else:
            out.append((r.episode, s, e))
    return [(ep, (s + (e - s) / 2).normalize()) for ep, s, e in out]


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "detect":
        for sta in (sys.argv[2:] or ["ALBH", "SC04", "SEAT"]):
            run(sta)
    else:
        anchored("ALBH", svi_anchors(), "SVI_ALBH")
        anchored("SC04", wa_anchors(), "WA_SC04")
        anchored("SEAT", wa_anchors(), "WA_SEAT")
