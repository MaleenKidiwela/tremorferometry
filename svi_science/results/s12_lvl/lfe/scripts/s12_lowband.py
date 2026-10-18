#!/usr/bin/env python
"""s12_lowband -- PLAN §7: 1-3 Hz restack of SVI family grand stacks from the existing raw day files (no download), then
the §2 reflection scan (scrambled control, 1000 draws) with tag _lb. Per causal (family, core station) pair: up to NMAX
highest-P data repeats (P >= 0.5); window tS = detector time + template offset; s10_windows reader and processing with
the band changed to 1-3 Hz; RMS-normalised; two passes of cc realignment (+-0.3 s, S-1..S+2 s, joint 3C); S re-picked on
the horizontal envelope of the stack.
    python s12_lowband.py build [--workers 7] [--nmax 120]
    python s12_lowband.py scan [--nscr 1000]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, time, argparse, pickle
import numpy as np
import pandas as pd
ROOT = "/home/jovyan/tremorferometry"
sys.path.insert(0, f"{ROOT}/svi_science/scripts"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s10_common as S
import s10_windows as W
import s11_common as C
import s12_reflect as RF

OUT = RF.OUT; TB = RF.TB
MT = f"{ROOT}/svi_science/results/s10_dvv/multi/tables"; XT = f"{C.TB}/ext"
_SOS = {}


def sos_lb(fs):
    from scipy.signal import butter
    if fs not in _SOS:
        _SOS[fs] = butter(4, [1.0 / (fs / 2), 3.0 / (fs / 2)], btype="band", output="sos")
    return _SOS[fs]


W.sos_for = sos_lb                      # band 1-3 Hz for every window read here (module attribute patched, file unchanged)


def log(m):
    s = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"
    print(s, flush=True); open(f"{RF.FW.LOG}/lowband.log", "a").write(s + "\n")


def plan(nmax):
    ser = pd.concat([pd.read_csv(f"{C.TB}/series.csv"), pd.read_csv(f"{C.TB}/series_ext.csv")], ignore_index=True)
    ser = ser[(ser.n_data > 0) & (ser.causal25_3c > 1.5)]
    pairs = sorted(set(zip(ser.fam.astype(int), ser.netsta)))
    rows = pd.concat([pd.read_parquet(f"{MT}/rows.parquet"), pd.read_parquet(f"{XT}/rows.parquet")], ignore_index=True)
    rows = rows[(rows.series == "data") & (rows.P >= 0.5)].drop_duplicates(["fam", "t"])
    off = pd.concat([pd.read_csv(f"{MT}/offsets_sp.csv"), pd.read_csv(f"{XT}/offsets_sp.csv")]).drop_duplicates(["fam", "netsta"])
    off = off.set_index(["fam", "netsta"])
    jobs = {}
    for fam, ns in pairs:
        if (fam, ns) not in off.index:
            continue
        r = rows[rows.fam == fam].sort_values("P", ascending=False).head(nmax)
        o = off.loc[(fam, ns)]
        for t, p in zip(r.t, r.P):
            jobs.setdefault(ns, []).append((fam, float(t) + float(o.offset), float(p), float(o.sp)))
    return jobs


def station_job(arg):
    ns, lst = arg
    raw = W.Raw(ns)
    lst = sorted(lst, key=lambda x: x[1])
    out = {}
    for fam, tS, p, sp in lst:
        win = np.full((3, S.NW), np.nan)
        try:
            for c, v in enumerate(raw.window(tS)):
                if v is None:
                    continue
                y = W.process(v[0], v[1], v[2], tS)
                if y is not None:
                    win[c] = y
        except Exception:
            continue
        if not np.isfinite(win[0]).all():
            continue
        for c in range(3):
            if np.isfinite(win[c]).all():
                r = np.sqrt(np.mean(win[c] ** 2))
                win[c] = win[c] / r if r > 0 else np.nan
        out.setdefault(fam, []).append((win.astype(np.float32), p, sp))
    return ns, out


def realign(Ws, P):
    tau = S.TAU; m = (tau >= -1.0) & (tau <= 2.0); idx = np.flatnonzero(m); L = np.arange(-12, 13)
    X = np.nan_to_num(Ws.astype(np.float64))
    for it in range(2):
        G = np.einsum("n,nct->ct", P, X) / P.sum()
        sh = np.zeros(len(X), int)
        for i in range(len(X)):
            best, bl = -9, 0
            for l in L:
                v = 0.0
                for c in range(3):
                    g = G[c, idx]; x = X[i, c, idx + l]
                    v += (x @ g) / (np.linalg.norm(x) * np.linalg.norm(g) + 1e-12)
                if v > best:
                    best, bl = v, l
            sh[i] = bl
        X = np.array([np.roll(X[i], -sh[i], axis=1) for i in range(len(X))])
    return np.einsum("n,nct->ct", P, X) / P.sum()


def cmd_build(a):
    t0 = time.time()
    jobs = plan(a.nmax)
    log(f"build: {sum(len(v) for v in jobs.values())} windows at {len(jobs)} stations (nmax {a.nmax})")
    import multiprocessing as mp
    res = {}
    with mp.get_context("fork").Pool(a.workers) as p:
        for ns, out in p.imap_unordered(station_job, list(jobs.items())):
            res[ns] = out; log(f"  {ns}: {sum(len(v) for v in out.values())} windows, {len(out)} families ({time.time() - t0:.0f}s)")
    P = []
    for ns, out in res.items():
        for fam, lst in out.items():
            if len(lst) < 15:
                continue
            Wn = np.array([x[0] for x in lst]); Pw = np.array([x[1] for x in lst]); sp = lst[0][2]
            G = realign(Wn, Pw)
            from scipy.signal import hilbert
            eh = np.abs(hilbert(G[1:], axis=1)) ** 2
            eh = eh.sum(0)
            m = (S.TAU >= -1.0) & (S.TAU <= 1.0)
            tS = float(S.TAU[m][np.argmax(eh[m])])
            X = np.array([np.interp(RF.TREL + tS, S.TAU, G[c], left=np.nan, right=np.nan) for c in range(3)])
            P.append(dict(fam=int(fam), netsta=ns, X=X, n=float(len(lst)), sp=sp, epochs="lb"))
    pickle.dump(P, open(f"{TB}/lowband_pairs.pkl", "wb"))
    open(f"{OUT}/LOWBAND_BUILD_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    log(f"build done: {len(P)} pairs ({time.time() - t0:.0f}s)")


def cmd_scan(a):
    P = pickle.load(open(f"{TB}/lowband_pairs.pkl", "rb"))
    RF.TAG = "_lb"; RF.SMN, RF.CEN, RF.FL0, RF.FL1 = 13, 0.15, 0.45, 0.9
    RF.LOADER = lambda: P
    sys.argv = [sys.argv[0], "--nscr", str(a.nscr)]
    RF.main()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd"); ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--nmax", type=int, default=120); ap.add_argument("--nscr", type=int, default=1000)
    a = ap.parse_args()
    {"build": cmd_build, "scan": cmd_scan}[a.cmd](a)
