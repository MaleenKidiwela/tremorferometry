#!/usr/bin/env python
"""S11 tomo (PLAN.md Amendment 3): pooled zone Wadati, one κ per zone with per-family intercepts.

  python s11_tomo_pooled.py --sets A,B,C,Bcl3S2,Ccl3S2
writes tables/pooled_zones.csv (+ per-zone scrambled-control distributions in pooled_scr.npz)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

T = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/tomo"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s11_tomo_wadati import along_strike  # noqa: E402

NB, NS = 1000, 1000


def pooled_fit(fid, x, y):
    """slope with per-family intercepts: demean x and y within family, OLS through origin."""
    df = pd.DataFrame(dict(f=fid, x=x, y=y))
    xm = df.groupby("f").x.transform("mean"); ym = df.groupby("f").y.transform("mean")
    xd, yd = (df.x - xm).values, (df.y - ym).values
    sxx = (xd ** 2).sum()
    if sxx <= 0:
        return np.nan, np.full(len(x), np.nan)
    b = (xd * yd).sum() / sxx
    return b, yd - b * xd


NSTA_MIN = 4


def clip_fit(fid, x, y, nfam_min=5, nsta_min=None):
    nsta_min = NSTA_MIN if nsta_min is None else nsta_min
    """V1 clipping on pooled residuals; drop families falling below nsta_min; returns slope, keep mask, rms."""
    keep = np.ones(len(x), bool)
    for _ in range(20):
        f_ok = pd.Series(fid[keep]).value_counts()
        f_ok = set(f_ok.index[f_ok >= nsta_min])
        keep &= np.isin(fid, list(f_ok))
        if len(f_ok) < nfam_min:
            return np.nan, keep, np.nan, len(f_ok)
        b, r = pooled_fit(fid[keep], x[keep], y[keep])
        full = np.full(len(x), np.nan); full[keep] = r
        # residual for all rows incl. clipped ones
        df = pd.DataFrame(dict(f=fid, x=x, y=y))
        xm = df[keep].groupby("f").x.mean(); ym = df[keep].groupby("f").y.mean()
        rr = (df.y - df.f.map(ym)) - b * (df.x - df.f.map(xm))
        rr = rr.values
        med = np.nanmedian(rr[keep]); mad = np.nanmedian(np.abs(rr[keep] - med))
        thr = max(3 * 1.4826 * mad, 0.15)
        k2 = np.isfinite(rr) & (np.abs(rr - med) <= thr) & np.isin(fid, list(f_ok))
        if (k2 == keep).all():
            break
        keep = k2
    f_ok = pd.Series(fid[keep]).value_counts()
    f_ok = set(f_ok.index[f_ok >= nsta_min]); keep &= np.isin(fid, list(f_ok))
    if len(f_ok) < nfam_min:
        return np.nan, keep, np.nan, len(f_ok)
    b, r = pooled_fit(fid[keep], x[keep], y[keep])
    return b, keep, float(np.sqrt(np.mean(r ** 2))), len(f_ok)


EST = "ols"


def kap(b):
    return 1 + b if EST == "ols" else (1 / b if b != 0 else np.nan)


def xy(tp, ts):
    """ols: (tS - tP) on tP, kappa = 1 + b;  inv (E1): tP on tS, kappa = 1/b."""
    return (tp, ts - tp) if EST == "ols" else (ts, tp)


def zone(Q, rng, label, robust=True):
    fid, tp, s, sta = Q.fam.values, Q.tP.values, Q.tS.values, Q.sta.values
    x, y = xy(tp, s)
    if robust:
        b, keep, rms, nf = clip_fit(fid, x, y)
    else:
        b, r = pooled_fit(fid, x, y); keep = np.ones(len(x), bool); rms = float(np.sqrt(np.mean(r ** 2)))
        nf = len(set(fid))
    out = dict(zone=label, n_fam_in=len(set(fid)), n_rows_in=len(x), n_fam=nf, n_rows=int(keep.sum()),
               kappa=kap(b), b=b, rms=rms, est=EST)
    if not np.isfinite(b):
        return out, None
    fk, xk, yk, sk, tk = fid[keep], x[keep], y[keep], sta[keep], s[keep]
    fams = np.unique(fk)
    out["tP_span_med"] = float(pd.Series(xk).groupby(fk).agg(lambda v: v.max() - v.min()).median())
    # bootstrap over families
    idx = {f: np.where(fk == f)[0] for f in fams}
    bb = []
    for _ in range(NB):
        pick = rng.choice(fams, len(fams))
        ii = np.concatenate([idx[f] for f in pick])
        lab = np.concatenate([[k] * len(idx[f]) for k, f in enumerate(pick)])
        bb.append(kap(pooled_fit(lab, xk[ii], yk[ii])[0]))
    out["sig_fam"] = float(np.nanstd(bb))
    # bootstrap over stations
    stas = np.unique(sk); bs = []
    for _ in range(NB):
        pick = rng.choice(stas, len(stas))
        ii = np.concatenate([np.where(sk == st)[0] for st in pick])
        lab = np.array([f"{fk[i]}" for i in ii])
        bs.append(kap(pooled_fit(lab, xk[ii], yk[ii])[0]))
    out["sig_sta"] = float(np.nanstd(bs))
    # leave one station out
    lo = {st: kap(pooled_fit(fk[sk != st], xk[sk != st], yk[sk != st])[0]) for st in stas}
    dk = {k: v - out["kappa"] for k, v in lo.items() if np.isfinite(v)}
    kmax = max(dk, key=lambda k: abs(dk[k]))
    out.update(loso_max=float(abs(dk[kmax])), loso_sta=kmax)
    # joint station terms
    st_i = {st: i for i, st in enumerate(stas)}
    df = pd.DataFrame(dict(f=fk, x=xk, y=yk, s=sk))
    # demean within family then solve y_d = b x_d + (c_s - mean_f c)
    F = pd.get_dummies(df.s).astype(float)
    Fd = F - F.groupby(df.f).transform("mean")
    xd = df.x - df.groupby("f").x.transform("mean"); yd = df.y - df.groupby("f").y.transform("mean")
    A = np.c_[xd.values, Fd.values]
    A = np.vstack([A, np.r_[0, np.ones(F.shape[1])]]); yv = np.r_[yd.values, 0]
    sol = np.linalg.lstsq(A, yv, rcond=None)[0]
    out.update(kappa_site=kap(sol[0]), d_site=float(kap(sol[0]) - kap(b)), site_sd=float(np.std(sol[1:])))
    # scrambled control: permute S picks across stations within each family (on the zone's input rows), same clipping
    sc = []
    for _ in range(NS):
        s2 = s.copy()
        for f in np.unique(fid):
            j = np.where(fid == f)[0]; s2[j] = s[rng.permutation(j)]
        xs, ys = xy(tp, s2)
        bc = clip_fit(fid, xs, ys)[0] if robust else pooled_fit(fid, xs, ys)[0]
        sc.append(bc)
    sc = np.array(sc)                      # control slopes (b scale; z on b scale, kappa columns transformed)
    out.update(scr_b_med=float(np.nanmedian(sc)), scr_b_sd=float(np.nanstd(sc)),
               scr_med=float(kap(np.nanmedian(sc))), scr_frac_measured=float(np.isfinite(sc).mean()))
    out["scr_sd"] = out["scr_b_sd"]
    out["z_vs_scr"] = (b - out["scr_b_med"]) / out["scr_b_sd"] if out["scr_b_sd"] > 0 else np.nan
    sig = max(out["sig_fam"], out["sig_sta"])
    out["sig"] = sig
    out["pass"] = bool(sig <= 0.05 and abs(out["z_vs_scr"]) > 4 and out["loso_max"] <= 0.05
                       and abs(out["d_site"]) <= 0.02)
    out["physical"] = bool(1.4 <= out["kappa"] <= 2.6)
    return out, sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", default="A,B,C,Bcl3S2,Ccl3S2")
    ap.add_argument("--rmin", type=float, default=1.2)
    ap.add_argument("--out", default="pooled_zones")
    ap.add_argument("--est", default="ols", choices=["ols", "inv"])
    ap.add_argument("--nsta-min", type=int, default=4)
    a = ap.parse_args()
    global NSTA_MIN
    NSTA_MIN = a.nsta_min
    global EST
    EST = a.est
    rng = np.random.default_rng(7)
    rows, scr = [], {}
    for st in a.sets.split(","):
        P = pd.read_csv(f"{T}/tables/picks_{st}.csv", low_memory=False)
        P["fam"] = P.fam.astype(str)
        P = P[P.tP.notna() & P.tS.notna() & (P.rP >= a.rmin) & (P.rS >= a.rmin)].copy()
        if "eye_lat" in P:
            P["lat"] = P.eye_lat.fillna(P.hyp_lat); P["lon"] = P.eye_lon.fillna(P.hyp_lon)
            P["depth_km"] = P.eye_dep.fillna(P.hyp_dep)
        P["along"] = along_strike(P.lat.values, P.lon.values)
        zones = [("all", P)]
        if "eye_class" in P:
            zones = [("eye_lfe", P[P.eye_class == "lfe"]), ("eye_noise", P[P.eye_class == "noise"])]
        else:
            lo = np.floor(P.along.min() / 30) * 30
            for z0 in np.arange(lo, P.along.max() + 1, 30):
                zones.append((f"along_{z0:+.0f}_{z0 + 30:+.0f}", P[(P.along >= z0) & (P.along < z0 + 30)]))
            for d0, d1 in ((0, 30), (30, 36), (36, 42), (42, 200)):
                zones.append((f"depth_{d0}_{d1}", P[(P.depth_km >= d0) & (P.depth_km < d1)]))
        for lab, Q in zones:
            if Q.fam.nunique() < 5:
                rows.append(dict(set=st, zone=lab, n_fam_in=Q.fam.nunique(), note="< 5 families"))
                continue
            for rob in (True, False):
                o, sc = zone(Q, rng, lab, rob)
                o.update(set=st, robust=rob, along_med=float(Q.along.median()),
                         depth_med=float(Q.depth_km.median()))
                rows.append(o)
                if sc is not None:
                    scr[f"{st}|{lab}|{rob}"] = sc
                print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in o.items()
                       if k in ("set", "zone", "robust", "n_fam", "n_rows", "kappa", "sig_fam", "sig_sta",
                                "rms", "z_vs_scr", "loso_max", "d_site", "pass")}, flush=True)
    pd.DataFrame(rows).to_csv(f"{T}/tables/{a.out}.csv", index=False)
    np.savez_compressed(f"{T}/tables/{a.out}_scr.npz", **scr)


if __name__ == "__main__":
    main()
