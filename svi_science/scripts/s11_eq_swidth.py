#!/usr/bin/env python
"""S11 tomo/eq (PLAN.md §6 Amendment 3): S/P pulse width vs distance on the LFE stacks, onset-pick κ and ≤ 50 km κ.

  python s11_eq_swidth.py
writes tomo/eq/tables/swidth_rows.csv, swidth_slopes.csv, swidth_kappa.csv
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.signal import hilbert, resample_poly  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s11_tomo_pooled as PO  # noqa: E402
import s11_tomo_profile as PR  # noqa: E402

ROOT = "/home/jovyan/tremorferometry"
TT = f"{ROOT}/svi_science/results/s11_depth/tomo/tables"
EQ = f"{ROOT}/svi_science/results/s11_depth/tomo/eq"
FS = 40.0
SETS = {"WA_NW": "BV3", "WA_SW": "CV3", "SVI": "EallV3"}


def width(env, k, search=0.5):
    """peak within ±search s of sample k; FWHM (s) and left half-max crossing (s)."""
    n = env.size
    a, b = max(int(k - search * FS), 0), min(int(k + search * FS) + 1, n)
    if b - a < 3:
        return np.nan, np.nan, np.nan
    kp = a + int(np.argmax(env[a:b]))
    h = 0.5 * env[kp]
    lo = kp
    while lo > 0 and env[lo] > h and kp - lo < 5 * FS:
        lo -= 1
    hi = kp
    while hi < n - 1 and env[hi] > h and hi - kp < 5 * FS:
        hi += 1
    return (hi - lo) / FS, lo / FS, kp / FS


def load_rows():
    xy = PR.coords()
    out = []
    cache = {}
    for reg, ps in SETS.items():
        P = pd.read_csv(f"{TT}/picks_{ps}.csv", low_memory=False)
        P = P[P.tS.notna() & (P.rS >= 1.2)].copy()
        if reg == "SVI":
            off = P.frame_offset.values
        else:
            off = np.zeros(len(P))
            P["epi_km"] = [np.hypot((r.lon - xy[r.sta][1]) * 111.19 * np.cos(np.radians(r.lat)),
                                    (r.lat - xy[r.sta][0]) * 111.19) if r.sta in xy else np.nan for r in P.itertuples()]
        for r, o in zip(P.itertuples(), off):
            if reg == "SVI":
                p = f"{TT}/era_stacks/{r.sta}.npz"
                if p not in cache:
                    cache[p] = np.load(p)
                st = resample_poly(cache[p][str(int(r.grp))].astype(np.float64), 2, 1, axis=-1)
            else:
                p = f"{ROOT}/{r.src}"
                if p not in cache:
                    cache[p] = np.load(p)
                st = cache[p][f"{int(r.grp)}_{r.sta}"].astype(np.float64)
            eh = np.sqrt(np.abs(hilbert(st[1])) ** 2 + np.abs(hilbert(st[2])) ** 2)
            ez = np.abs(hilbert(st[0]))
            wS, onS, pkS = width(eh, (r.tS - o) * FS)
            wP, onP, pkP = (width(ez, (r.tP - o) * FS) if np.isfinite(r.tP) else (np.nan,) * 3)
            out.append(dict(reg=reg, fam=str(r.fam), sta=r.sta, epi_km=r.epi_km, depth_km=r.depth_km, snr=r.snr,
                            rP=r.rP, rS=r.rS, tP=r.tP, tS=r.tS, wS=wS, wP=wP, tS_on=onS + o, tP_on=onP + o,
                            tS_pk=pkS + o, tP_pk=pkP + o))
        print(reg, len(P), flush=True)
    return pd.DataFrame(out)


def slope_boot(d, col, rng, nb=1000):
    d = d.dropna(subset=[col, "epi_km"])
    if len(d) < 20:
        return dict(n=len(d))
    b = np.polyfit(d.epi_km, d[col], 1)[0]
    fams = d.fam.unique(); idx = {f: np.where(d.fam.values == f)[0] for f in fams}
    bs = []
    for _ in range(nb):
        ii = np.concatenate([idx[f] for f in rng.choice(fams, len(fams))])
        bs.append(np.polyfit(d.epi_km.values[ii], d[col].values[ii], 1)[0])
    rho, p = spearmanr(d.epi_km, d[col])
    near, far = d[d.epi_km <= 30][col].median(), d[d.epi_km >= 60][col].median()
    return dict(n=len(d), slope_ms_per_km=1000 * b, slope_sd=1000 * float(np.std(bs)), rho=rho, p=p,
                med_le30=near, med_ge60=far, med_all=d[col].median())


def main():
    rng = np.random.default_rng(23)
    f = f"{EQ}/tables/swidth_rows.csv"
    R = load_rows() if not os.path.exists(f) else pd.read_csv(f)
    R.to_csv(f, index=False)
    S = []
    for reg in SETS:
        d = R[R.reg == reg]
        for col, m in (("wS", d.rS >= 1.2), ("wP", d.rP >= 1.2), ("wS_minus_wP", (d.rS >= 1.2) & (d.rP >= 1.2))):
            q = d[m].copy()
            if col == "wS_minus_wP":
                q[col] = q.wS - q.wP
            S.append(dict(reg=reg, col=col, **slope_boot(q, col, rng)))
            print(S[-1], flush=True)
    pd.DataFrame(S).to_csv(f"{EQ}/tables/swidth_slopes.csv", index=False)
    PO.NB, PO.NS, PO.EST, PO.NSTA_MIN = 1000, 300, "inv", 4
    K = []
    for reg in SETS:
        d = R[(R.reg == reg) & (R.rP >= 1.2) & (R.rS >= 1.2)]
        for pick, tp, ts in (("peak (standard)", "tP", "tS"), ("halfmax onset", "tP_on", "tS_on")):
            for dmax in (40, 50, 60, 1000):
                for nmin in ((3, 4) if dmax <= 50 else (4,)):
                    q = d[d.epi_km <= dmax]
                    q = pd.DataFrame(dict(fam=q.fam, sta=q.sta, tP=q[tp], tS=q[ts])).dropna()
                    q = q[q.tS > q.tP].groupby("fam").filter(lambda g: len(g) >= nmin)
                    row = dict(reg=reg, pick=pick, dmax=dmax, nmin=nmin, n_fam_in=q.fam.nunique())
                    if q.fam.nunique() >= 5:
                        PO.NSTA_MIN = nmin
                        o, _ = PO.zone(q, rng, "w", True)
                        row.update({k: o[k] for k in ("n_fam", "n_rows", "kappa", "sig", "rms", "z_vs_scr", "loso_max",
                                                      "d_site", "pass")})
                    K.append(row)
                    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    pd.DataFrame(K).to_csv(f"{EQ}/tables/swidth_kappa.csv", index=False)


if __name__ == "__main__":
    main()
