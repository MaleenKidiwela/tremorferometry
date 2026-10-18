#!/usr/bin/env python
"""S11 tomo/eq (PLAN.md §6 Amendment 2b): synthetic emergent LFE-like records with known κ at the real WA LFE
(family, station) geometry and SNR, real station noise, the real LFE pickers -> recovered κ − true κ.

  python s11_eq_synth.py [--workers 8] [--ndraw 3]
writes tomo/eq/tables/synth_rows.csv, synth_kappa.csv
"""
import argparse
import os
import sys
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s11_eq_pickertest as PT  # noqa: E402
import s11_eq_diff as D  # noqa: E402
import s11_tomo_pooled as PO  # noqa: E402
import s11_tomo_profile as PR  # noqa: E402

A = PT.A
FS = 40.0
PRE, LEN = 4.0, 44.0
N = int(LEN * FS)
KAPPAS = (1.75, 1.90, 2.10)
REG = {"WA_NW": ("BV3", "royer_NW_WA"), "WA_SW": ("CV3", "royer_SW_WA")}


def wavelet(rng, f0=3.0, tau=0.25):
    t = np.arange(int(3.0 * FS)) / FS
    w = t ** 2 * np.exp(-t / tau) * np.cos(2 * np.pi * f0 * t + rng.uniform(0, 2 * np.pi))
    return w / np.abs(w).max()


def place(x, k0, w, amp):
    k0 = int(round(k0))
    a, b = max(k0, 0), min(k0 + len(w), x.size)
    if b > a:
        x[a:b] += amp * w[a - k0:b - k0]


def synth(tp, ts, rng, leak):
    st = np.zeros((3, N))
    wp, ws = wavelet(rng), wavelet(rng)
    az = rng.uniform(0, 2 * np.pi); baz = rng.uniform(0, 2 * np.pi)
    kp, ks = (PRE + tp) * FS, (PRE + ts) * FS
    place(st[0], kp, wp, 1.0)
    place(st[1], kp, wp, 0.3 * np.cos(baz)); place(st[2], kp, wp, 0.3 * np.sin(baz))
    place(st[1], ks, ws, 4.0 * np.cos(az)); place(st[2], ks, ws, 4.0 * np.sin(az))
    if leak > 0:
        place(st[0], ks, ws, 4.0 * leak)
    return st


def noise_bank(sta, nwin=60, seed=0):
    """real 1-8 Hz 40-Hz noise windows from up to 3 cached days of this station."""
    import glob
    rng = np.random.default_rng(seed)
    files = sorted(glob.glob(f"{A.CACHE40}/{A.netsta(sta)}/*/*.npz"))
    if not files:
        return []
    out = []
    for p in rng.choice(files, min(3, len(files)), replace=False):
        z = np.load(p)
        X = z["x"]; m = z["m"].astype(bool) if "m" in z.files else np.ones(X.shape[1], bool)
        if X.shape[0] != 3 or X.shape[1] <= N:
            continue
        got = 0
        for _ in range(400):
            i0 = int(rng.integers(0, X.shape[1] - N))
            if m[i0:i0 + N].all():
                s = X[:, i0:i0 + N].astype(np.float64)
                if np.isfinite(s).all() and (s ** 2).sum() > 0:
                    out.append(s); got += 1
            if got >= nwin:
                break
    return out[:nwin * 3]


def do_station(args):
    reg, sta, Q, ndraw, seed = args
    rng = np.random.default_rng(seed)
    bank = noise_bank(sta, seed=seed)
    rows = []
    for r in Q.itertuples():
        for kap in KAPPAS:
            ts = kap * r.tp
            for leak in (0.3, 0.0):
                for d in range(ndraw + 1):                 # d = 0: noise-free
                    st = synth(r.tp, ts, rng, leak)
                    if d > 0:
                        if not bank:
                            continue
                        st, _ = PT.noise_match(st, bank[int(rng.integers(len(bank)))], r.snr)
                    pr = PT.pick(st)
                    rows.append(dict(reg=reg, fam=r.fam, sta=sta, kappa_true=kap, leak=leak, draw=d, snr_target=r.snr,
                                     tp_true=r.tp, ts_true=ts,
                                     **{k: (v - PRE if k.startswith("t") else v) for k, v in pr.items()}))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--ndraw", type=int, default=3)
    a = ap.parse_args()
    f = f"{PT.EQ}/tables/synth_rows.csv"
    if not os.path.exists(f):
        jobs = []
        xy = PR.coords()
        for reg, (ps, prof) in REG.items():
            L = pd.read_csv(f"{PT.TT}/picks_{ps}.csv", low_memory=False)
            L = L[(L.depth_km >= 25) & (L.depth_km <= 50) & L.snr.notna() & L.sta.isin(list(xy))].copy()
            L["epi"] = [np.hypot((r.lon - xy[r.sta][1]) * 111.19 * np.cos(np.radians(r.lat)), (r.lat - xy[r.sta][0]) * 111.19)
                        for r in L.itertuples()]
            L = L[L.epi <= 100]
            pf = D.load_profile(prof)
            L["tp"] = np.nan
            for z, g in L.groupby(L.depth_km.round(0)):
                L.loc[g.index, "tp"] = D.ray_partials(pf, z, g.epi.values).sum(1)
            L["fam"] = L.fam.astype(str)
            for i, (sta, g) in enumerate(L.groupby("sta")):
                jobs.append((reg, sta, g[["fam", "tp", "snr"]], a.ndraw, 100 + i + (0 if reg == "WA_NW" else 50)))
        print("jobs", len(jobs), "rows", sum(len(j[2]) for j in jobs), flush=True)
        with Pool(a.workers) as pool:
            res = pool.map(do_station, jobs, chunksize=1)
        R = pd.DataFrame([r for rr in res for r in rr])
        R.to_csv(f, index=False)
    R = pd.read_csv(f)
    rng = np.random.default_rng(17)
    PO.NB, PO.NS, PO.EST, PO.NSTA_MIN = 300, 100, "inv", 4
    out = []
    for (reg, kap, leak, d), g in R.groupby(["reg", "kappa_true", "leak", "draw"]):
        for pk, tp, rp in (("V3", "tPv3", "rPv3"), ("frozen", "tP0", "rP0")):
            q = g[(g[rp] >= 1.2) & (g.rS >= 1.2)]
            q = pd.DataFrame(dict(fam=q.fam.astype(str), sta=q.sta, tP=q[tp], tS=q.tS)).dropna()
            q = q[q.tS > q.tP].groupby("fam").filter(lambda x: len(x) >= 4)
            row = dict(reg=reg, kappa_true=kap, leak=leak, draw=d, picker=pk, n_rows_r12=len(q), n_fam=q.fam.nunique())
            if q.fam.nunique() >= 5:
                o, _ = PO.zone(q, rng, "synth", True)
                row.update({k: o[k] for k in ("n_rows", "kappa", "sig", "rms", "loso_max", "d_site", "pass")})
                row["bias"] = row["kappa"] - kap
            dd = g.dropna(subset=[tp])
            row.update(med_dtP=float((dd[tp] - dd.tp_true).median()), frac_dtP_gt05=float(((dd[tp] - dd.tp_true).abs() > 0.5).mean()),
                       med_dtS=float((g.tS - g.ts_true).median()), frac_dtS_gt05=float(((g.tS - g.ts_true).abs() > 0.5).mean()))
            out.append(row)
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    pd.DataFrame(out).to_csv(f"{PT.EQ}/tables/synth_kappa.csv", index=False)


if __name__ == "__main__":
    main()
