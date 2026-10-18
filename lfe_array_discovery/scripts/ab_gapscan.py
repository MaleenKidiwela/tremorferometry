#!/usr/bin/env python
"""STAGE 0 bug audit, part 1: how much of the Washington 20 Hz cache is zero-filled gap?

The cache is built with obspy `merge(method=1, fill_value=0)` (lfe_picker/scripts/o_cache_month.py:70)
and stores NO per-sample validity mask, so a data gap reaches every downstream correlation as a run
of exact zeros.  This script measures, per station-day of Royer & Bostock's 14 named stations over
their six episode windows:
  frac_zero       fraction of samples where all three components are exactly 0
  maxrun_s        longest such run, seconds
  n_win_zero      15 s windows (0.5 s hop) whose 3-C energy is exactly 0
  emed_zero_comp  how many of the 3 components have a MEDIAN 15-s window energy of 0
                  -> this is the trigger for the degenerate guard in k_auto.py:69-70 and
                     k_grow.py:167-169 (`good = ew > 1e-4 * emed` becomes `ew > 0`)
  hour_*          the same, per hour, since k_auto's MAD is an HOUR statistic
-> results/ab_royer/gapscan_cache20.csv  and  gapscan_hours.csv
"""
from __future__ import annotations
import glob, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
OUT = f"{STUDY}/results/ab_royer"
os.makedirs(OUT, exist_ok=True)
FS = 20.0
WIN_S, HOP_S = 15.0, 0.5

NW = ["B001", "SQM", "BS11", "W020", "W040", "GNW", "PL11"]
SW = ["N050", "S030", "W070", "N060", "PL11", "S040", "N070"]
EPS = {"NW07": ("2007-01-17", "2007-01-30"), "NW08": ("2008-05-01", "2008-05-23"),
       "NW10": ("2010-08-07", "2010-08-24"), "NW11": ("2011-08-04", "2011-08-23"),
       "SW07": ("2007-01-14", "2007-02-01"), "SW08": ("2008-05-01", "2008-05-25")}
ROOTS = [f"{STUDY}/cache20hz", "/home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz"]


def index():
    idx = {}
    for r in ROOTS:
        for p in glob.glob(f"{r}/*/*/*.npz"):
            parts = p.split(os.sep)
            sta = parts[-3].split(".")[1]
            idx.setdefault((sta, int(parts[-2]), int(parts[-1][:3])), p)
    return idx


def win_energy(x, n):
    c = np.concatenate([[0.0], np.cumsum(x.astype(np.float64) ** 2)])
    return c[n:] - c[:-n]


def maxrun(mask):
    if not mask.any():
        return 0
    d = np.diff(np.concatenate([[0], mask.view(np.int8), [0]]))
    return int((np.flatnonzero(d == -1) - np.flatnonzero(d == 1)).max())


def main():
    idx = index()
    n, hop = int(WIN_S * FS), int(HOP_S * FS)
    rows, hrows = [], []
    for ep, (a, b) in EPS.items():
        stas = NW if ep.startswith("NW") else SW
        for date in pd.date_range(a, b):
            for s in stas:
                key = (s, date.year, date.dayofyear)
                p = idx.get(key)
                if p is None:
                    rows.append(dict(ep=ep, sta=s, date=date.strftime("%Y-%m-%d"),
                                     status="absent"))
                    continue
                z = np.load(p)
                X = z["x"].astype(np.float64)
                if X.shape[0] != 3:
                    rows.append(dict(ep=ep, sta=s, date=date.strftime("%Y-%m-%d"),
                                     status=f"nchan{X.shape[0]}"))
                    continue
                dead = (np.abs(X).sum(0) == 0)
                ew = np.vstack([win_energy(X[c], n) for c in range(3)])
                emed = np.median(ew, axis=1)
                rows.append(dict(
                    ep=ep, sta=s, date=date.strftime("%Y-%m-%d"), status="ok",
                    nsamp=X.shape[1], frac_zero=float(dead.mean()),
                    maxrun_s=maxrun(dead) / FS,
                    n_win_zero=int((ew.sum(0) == 0).sum()), n_win=int(ew.shape[1]),
                    emed_zero_comp=int((emed == 0).sum()),
                    emed_min=float(emed.min()), emed_med=float(np.median(emed))))
                # ---- per-hour (k_auto's MAD block)
                t0 = float(z["t0"])
                h0 = pd.Timestamp(date).timestamp()
                for hh in range(24):
                    i0 = int(round((h0 + 3600 * hh - t0) * FS))
                    i1 = i0 + int(3600 * FS)
                    if i0 < 0 or i1 > X.shape[1]:
                        continue
                    Xh = X[:, i0:i1]
                    eh = np.vstack([win_energy(Xh[c], n) for c in range(3)])
                    em = np.median(eh, axis=1)
                    hrows.append(dict(
                        ep=ep, sta=s, date=date.strftime("%Y-%m-%d"), hour=hh,
                        frac_zero=float((np.abs(Xh).sum(0) == 0).mean()),
                        emed_zero_comp=int((em == 0).sum()),
                        frac_win_zero=float((eh.sum(0) == 0).mean()),
                        frac_win_tiny=float((eh < 1e-4 * np.maximum(em[:, None], 1e-300)).mean())))
        print(f"[gapscan] {ep} done", flush=True)
    R = pd.DataFrame(rows)
    H = pd.DataFrame(hrows)
    R.to_csv(f"{OUT}/gapscan_cache20.csv", index=False)
    H.to_csv(f"{OUT}/gapscan_hours.csv", index=False)
    ok = R[R.status == "ok"]
    print("\n=== station-days ===")
    print(f"present {len(ok)} / requested {len(R)}")
    print(f"frac_zero > 0      : {(ok.frac_zero > 0).sum()} ({(ok.frac_zero>0).mean():.1%})")
    print(f"frac_zero > 0.01   : {(ok.frac_zero > 0.01).sum()}")
    print(f"frac_zero > 0.5    : {(ok.frac_zero > 0.5).sum()}")
    print(f"any zero-energy 15s window : {(ok.n_win_zero > 0).sum()}")
    print(f"emed_zero_comp > 0 (DEGENERATE GUARD): {(ok.emed_zero_comp > 0).sum()}")
    print("\n=== station-hours ===")
    print(f"n = {len(H)}")
    print(f"frac_zero > 0                        : {(H.frac_zero > 0).sum()} ({(H.frac_zero>0).mean():.1%})")
    print(f"emed_zero_comp > 0 (DEGENERATE GUARD): {(H.emed_zero_comp > 0).sum()} ({(H.emed_zero_comp>0).mean():.2%})")
    print(f"frac_win_zero > 0.1                  : {(H.frac_win_zero > 0.1).sum()}")
    print(H.groupby("ep").agg(n=("hour", "size"), zero_hours=("frac_zero", lambda v: int((v > 0).sum())),
                              degen=("emed_zero_comp", lambda v: int((v > 0).sum()))).to_string())


if __name__ == "__main__":
    main()
