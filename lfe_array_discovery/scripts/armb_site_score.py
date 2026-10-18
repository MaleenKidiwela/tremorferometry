#!/usr/bin/env python
"""Site-averaged impulsiveness score of Arm-B long stacks (rule fixed by the orchestrator,
2026-09-29).  Each station stack 2-8 Hz (butter 4, sosfiltfilt, fs 20), normalised by max|.|,
averaged into sites; per site Hilbert envelope summed over 3 components e; impulsive if
e.max()/median(e) > 5 and (e > e.max()/2).sum()/40 < 3.  n_imp_sites per family.

  python scripts/armb_site_score.py --stacks <npz> --out <csv>
"""
import argparse, json, re
import numpy as np, pandas as pd
from scipy.signal import butter, sosfiltfilt, hilbert

ap = argparse.ArgumentParser()
ap.add_argument("--stacks", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

def site(s):
    if s in ("B005", "B006", "B007"):
        return "xPBB005"
    if s in ("B013", "B943"):
        return "xPBB013"
    p = s[:2]
    if p in ("BH", "BS", "DR", "GC", "LC", "PA", "PL", "TB") and re.fullmatch(r"[A-Z]{2}\d+", s):
        return "x" + p
    return s

z = np.load(a.stacks, allow_pickle=True)
RAW = [(k, np.asarray(z[k], np.float64)) for k in z.files if k != "meta"]


def score(design_fs, out):
    # design_fs=20: the stated rule, a true 2-8 Hz band on the 20 sps stacks.
    # design_fs=40: the filter DESIGNED at fs 40 and applied to 20 sps data (effective 1-4 Hz);
    # this variant reproduces the reference real-Arm-B file site_score_armb_siteavg.csv exactly
    # (197/117/79/61/48, 739/739 families identical), so it is the like-for-like comparison.
    sos = butter(4, [2.0, 8.0], btype="band", output="sos", fs=design_fs)
    fam = {}
    for k, x0 in RAW:
        g, s = k.split("_", 1)
        x = sosfiltfilt(sos, x0, axis=-1)
        m = np.abs(x).max()
        if not np.isfinite(m) or m <= 0:
            continue
        fam.setdefault(int(g), {}).setdefault(site(s), []).append(x / m)
    body(fam, out, design_fs)


def body(fam, out, design_fs):
    rows = []
    for g in sorted(fam):
        n_imp, n_sites, imp = 0, 0, []
        for st, L in fam[g].items():
            A = np.mean(L, axis=0)
            e = np.abs(hilbert(A, axis=-1)).sum(0)
            n_sites += 1
            med = np.median(e)
            ok = (med > 0 and e.max() / med > 5) and ((e > e.max() / 2).sum() / 40 < 3)
            if ok:
                n_imp += 1
                imp.append(st)
        rows.append(dict(grp=g, n_sites=n_sites, n_imp_sites=n_imp,
                         imp_sites=",".join(sorted(imp))))
    T = pd.DataFrame(rows)
    T.to_csv(out, index=False)
    print(f"{a.stacks} [filter designed at fs={design_fs:g}]: {len(T)} families -> {out}")
    for k in range(2, 7):
        print(f"  n_imp_sites >= {k}: {int((T.n_imp_sites >= k).sum())}")


score(20.0, a.out)
score(40.0, a.out.replace(".csv", "_fs40design.csv"))
