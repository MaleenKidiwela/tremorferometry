#!/usr/bin/env python
"""Merlin Q3/Q4 statistics for a daily-inversion product directory (res_catalog or splits/<name>):
  y_q(t)  = <pat_q, MODEL(t)> / <pat_q, pat_q>   raw patch projection (20-patch basis from res_catalog/transmission_T.npz)
  yT_q(t) = y_q(t) / T_q(t)                       divided by the directory's OWN daily transmission (if present)
  L       = mean(2022-01..2026-12) - mean(2011-01..2016-12)   two-window level difference (gauge-invariant)
Reports for deep q16 (47.2N), the leak-robust column sum q6+q16, and neighbours q15/q17, plus the deconvolved A16
as SECONDARY (own transmission_deconv.npz if present). Annual means of yT_16 printed for the shape comparison.
Usage: python tile_level_stats.py <dir> [<dir> ...]   (patterns always taken from res_catalog/transmission_T.npz)
"""
import sys, os
import numpy as np, pandas as pd

RD = "fault_tomography/inversion/res_catalog"
zt0 = np.load(f"{RD}/transmission_T.npz", allow_pickle=True)
PAT = np.asarray(zt0["patterns"], float); SELF = (PAT ** 2).sum(1)
W1 = (pd.Timestamp("2011-01-01"), pd.Timestamp("2016-12-31")); W2 = (pd.Timestamp("2022-01-01"), pd.Timestamp("2026-12-31"))


def level(s):
    a = s[(s.index >= W1[0]) & (s.index <= W1[1])].mean(); b = s[(s.index >= W2[0]) & (s.index <= W2[1])].mean()
    return b - a, a, b


def analyse(d):
    zp = np.load(f"{d}/inversion_4d_daily.npz", allow_pickle=True)
    days = pd.to_datetime(zp["days"]); MODEL = np.nan_to_num(np.asarray(zp["MODEL"], float))
    Y = (MODEL @ PAT.T) / SELF[None, :]
    T = None
    if os.path.exists(f"{d}/transmission_T.npz"):
        zt = np.load(f"{d}/transmission_T.npz", allow_pickle=True)
        assert list(zt["days"]) == list(zp["days"]) or len(zt["days"]) == len(days)
        T = np.clip(np.asarray(zt["T"], float), 0.05, None)
    YT = Y / T if T is not None else Y
    out = {}
    for name, ser in (("y16/T", YT[:, 16]), ("y6+y16 (/T)", YT[:, 6] + YT[:, 16]), ("y15/T", YT[:, 15]),
                      ("y17/T", YT[:, 17]), ("y16 raw", Y[:, 16])):
        s = pd.Series(ser, index=days)
        L, a, b = level(s)
        out[name] = dict(L=L, W1=a, W2=b)
    if os.path.exists(f"{d}/transmission_deconv.npz"):
        zd = np.load(f"{d}/transmission_deconv.npz", allow_pickle=True)
        s = pd.Series(np.asarray(zd["A"][:, 16], float), index=pd.to_datetime(zd["days"]))
        L, a, b = level(s.dropna()); out["A16 deconv (secondary)"] = dict(L=L, W1=a, W2=b)
    ann = pd.Series(YT[:, 16], index=days).resample("YE").mean()
    T16 = pd.Series(T[:, 16], index=days).resample("YE").mean() if T is not None else None
    return out, ann, T16, len(days)


if __name__ == "__main__":
    for d in sys.argv[1:]:
        out, ann, T16, nd = analyse(d)
        print(f"\n=== {d}  (frames {nd}) ===")
        for k, v in out.items():
            print(f"  {k:24s} L = {v['L']:+.3f}   (W1 {v['W1']:+.3f} -> W2 {v['W2']:+.3f})")
        print("  annual mean y16/T:", " ".join(f"{y.year}:{v:+.2f}" for y, v in ann.loc["2011":"2026"].items()))
        if T16 is not None:
            print("  annual mean T16:  ", " ".join(f"{y.year}:{v:.2f}" for y, v in T16.loc["2011":"2026"].items()))
