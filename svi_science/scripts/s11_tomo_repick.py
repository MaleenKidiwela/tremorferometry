#!/usr/bin/env python
"""S11 tomo (PLAN.md Amendment 5): V3 polarisation-weighted P re-pick on the same stacks; S unchanged.

  python s11_tomo_repick.py --set B        (also A, C, Bcl3S2, Ccl3S2, ...)
writes tables/picks_<set>V3.csv (tP/rP replaced; old values kept as tP_v0/rP_v0)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.signal import hilbert

ROOT = "/home/jovyan/tremorferometry"
T = f"{ROOT}/svi_science/results/s11_depth/tomo/tables"
FS = 40.0


def rotate(st):
    r = np.imag(hilbert(st, axis=-1))
    m = np.abs(r).max()
    return r / m if m > 0 else r


def pick_p_v3(st, tS, excl=0.5):
    r = rotate(st)
    ez = uniform_filter1d(st[0] ** 2, int(0.5 * FS))
    eh = uniform_filter1d(st[1] ** 2 + st[2] ** 2, int(0.5 * FS))
    w = ez / np.maximum(ez + eh, 1e-30)
    if (st[1] ** 2).sum() == 0:                       # vertical-only: no polarisation weight
        w = np.ones_like(ez)
    cf = np.abs(r[0]) * w
    kS = int(round(tS * FS))
    t0, t1 = max(kS - int(14 * FS), 0), max(kS - int(1.5 * FS), 0)
    if t1 - t0 < 4:
        return np.nan, np.nan, np.nan
    seg = cf[t0:t1]
    k = int(np.argmax(seg)); e = int(excl * FS)
    rest = np.r_[seg[:max(k - e, 0)], seg[k + e + 1:]]
    ratio = float(seg[k] / rest.max()) if rest.size and rest.max() > 0 else np.inf
    return (t0 + k) / FS, ratio, float(w[t0 + k])


def stack_loader(setname, P):
    """returns function (row) -> 3xN stack."""
    cache = {}

    def get(path):
        if path not in cache:
            cache[path] = np.load(path)
        return cache[path]
    if setname == "A":
        z = get(f"{ROOT}/lfe_array_discovery/results/ab_royer/NW/hypo_28/stacks_loc.npz")
        return lambda r: z[f"{int(r.grp)}_{r.sta}"]
    if "cl" in setname:                                 # e.g. Bcl3S2 -> clstacks_B_3_S2.npz
        base, rest = setname.split("cl")
        R, mode = rest[:-2], rest[-2:]
        z = get(f"{T}/clstacks_{base}_{R}_{mode}.npz")
        return lambda r: z[f"{r.fam[2:]}_{r.sta}"]
    return lambda r: get(f"{ROOT}/{r.src}")[f"{int(r.grp)}_{r.sta}"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True)
    a = ap.parse_args()
    P = pd.read_csv(f"{T}/picks_{a.set}.csv", low_memory=False)
    P["fam"] = P.fam.astype(str)
    get = stack_loader(a.set, P)
    P["tP_v0"], P["rP_v0"] = P.tP, P.rP
    tp, rp, wp = [], [], []
    for r in P.itertuples():
        if not np.isfinite(r.tS):
            tp.append(r.tP); rp.append(r.rP); wp.append(np.nan)     # vertical-only rows: keep v0 (no S anyway)
            continue
        x, y, w = pick_p_v3(get(r).astype(np.float64), r.tS)
        tp.append(x); rp.append(y); wp.append(w)
    P["tP"], P["rP"], P["wP"] = tp, rp, wp
    P.to_csv(f"{T}/picks_{a.set}V3.csv", index=False)
    ok0 = (P.rP_v0 >= 1.2) & (P.rS >= 1.2); ok = (P.rP >= 1.2) & (P.rS >= 1.2)
    sp0, sp = (P.tS - P.tP_v0)[ok0], (P.tS - P.tP)[ok]
    print(f"[{a.set}V3] rows rP&rS>=1.2: v0 {ok0.sum()} -> V3 {ok.sum()}; tS-tP<2.5 s share: v0 "
          f"{(sp0 < 2.5).mean():.2f} -> V3 {(sp < 2.5).mean():.2f}; changed picks "
          f"{(np.abs(P.tP - P.tP_v0) > 0.05).mean():.2f}", flush=True)


if __name__ == "__main__":
    main()
