#!/usr/bin/env python
"""s12_novel_measure -- S12-LVL LFE side, coordinator tasks A-D (PLAN.md §6). One pass over the S11 round-3 cc-weighted
day stores per region; analysis in s12_novel.py.
Per (family, station, epoch) series:
  * day table (every day): date, lag to the family's front (-10..60 d), major index, days since the previous / to the
    next major front, between-ETS flag, day weight, events per day; per-day coda dv/v (1-2, 2-3, 3-4, 4-6, 2-5 s, joint 3C
    stretch vs the series grand stack, bank +-3 %) and cc; per-day source/decay scalars: ln E slope 1-6 s, S-window
    spectral ratio ln(E 5-8 Hz / E 2-4 Hz) (3C), S envelope width (s, half maximum), ln(A_P/A_S) on Z.            [B, C]
  * daily front-lag class stacks D-10..D20 (all majors superposed) with the 200-draw joint bootstrap over majors, and
    the same daily classes with fronts scrambled to other majors (200 draws).                                      [A]
  * per-major stacks of the front-lag days 0..20 (+ 50-draw day bootstrap) and the between-ETS stack (joint bootstrap). [D]
eps = -dv/v; dv/v reported in %.
    S11_REGION=SVI|WA_NW python s12_novel_measure.py [--workers 12] [--test N]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, time, argparse, pickle
import numpy as np
import pandas as pd
from scipy.signal import hilbert
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s12_epoch as E
from s12_epoch import C, OUT, SUB, TB, WINS, log, load_store, front_tables, fronts_for, classify, day_scalars, StackObs

DAYS = list(range(-10, 21))
NB, NSCR, NBD = 200, 200, 50
_G = {}


def dvv_only(SO, X):
    o = {}
    for wn, R in SO.R.items():
        e, cc = R.measure(np.nan_to_num(X), (0, 1, 2))
        o[wn] = np.where(np.abs(e) >= 0.0299, np.nan, -100 * e); o[f"cc{wn}"] = cc
    return o


def src_scalars(St, tau, tS, sp):
    tr = tau - tS
    X = np.nan_to_num(St.astype(np.float64))
    w = (tr >= -0.2) & (tr <= 0.8); n = w.sum()
    F = np.fft.rfft(X[:, :, w] * np.hanning(n)[None, None, :], n=128, axis=2)
    f = np.fft.rfftfreq(128, 1 / 40.0)
    P = (np.abs(F) ** 2).sum(1)
    lo = P[:, (f >= 2) & (f < 4)].sum(1); hi = P[:, (f >= 5) & (f <= 8)].sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        sr = np.log(hi / lo)
    E_ = np.zeros((X.shape[0], X.shape[2]))
    for c in range(3):
        E_ += np.abs(hilbert(X[:, c], axis=1))
    k = np.ones(5) / 5
    E_ = np.apply_along_axis(lambda v: np.convolve(v, k, "same"), 1, E_)
    ws = (tr >= -0.6) & (tr <= 1.2); Es = E_[:, ws]
    pk = Es.max(1, keepdims=True)
    width = (Es >= 0.5 * pk).sum(1) / 40.0
    wP = (tr >= -sp - 0.15) & (tr <= -sp + 0.6); wS = (tr >= -0.15) & (tr <= 0.6)
    with np.errstate(divide="ignore", invalid="ignore"):
        ps = np.log(np.sqrt((X[:, 0, wP] ** 2).mean(1)) / np.sqrt((X[:, 0, wS] ** 2).mean(1)))
    return dict(sr=sr, width=width, lnPS_Z_rms=ps)


def stacks(St, Ws, msk, g, ng):
    Ps = np.zeros((ng, 3, St.shape[-1])); Wg = np.zeros((ng, 3))
    for c in range(3):
        np.add.at(Ps[:, c], g[msk], Ws[msk, c][:, None] * St[msk, c]); np.add.at(Wg[:, c], g[msk], Ws[msk, c])
    return Ps, Wg


def draw_stack(K, Ps, Wg):
    num = np.einsum("bg,gct->bct", K, Ps); den = K @ Wg
    return np.where(den[:, :, None] > 0, num / np.maximum(den[:, :, None], 1e-12), 0.0), den[:, 0]


def series_job(key):
    st = _G["store"][key]; fam = int(key.split("|")[0])
    if fam not in _G["fams"]:
        return key, None
    days = pd.to_datetime(pd.Series(st["days"])).dt.normalize()
    dn = (days.astype("datetime64[ns]").astype(np.int64) // 86400_000_000_000).to_numpy()
    M = _G["majors"]; nM = len(M)
    fr = fronts_for(fam, M, _G["start"], _G["medrel"], _G["rel"])
    midn = np.array([_G["mid"][m].normalize().value // 86400_000_000_000 for m in M], dtype=np.int64)
    cf, jf, cg, jg, betw, ib = classify(dn, fr, midn)
    prev = classify.prev.copy()
    lag = dn[:, None] - fr[None, :]
    jj = np.argmin(np.where((lag >= -10) & (lag <= 60), np.abs(lag), 10 ** 6), 1)
    lj = lag[np.arange(len(dn)), jj].astype(float); lj[(lj < -10) | (lj > 60)] = np.nan
    nxt = np.where(lag < 0, -lag, 10 ** 6).min(1).astype(float)
    St = np.nan_to_num(st["St"].astype(np.float64)); Ws = st["Ws"]
    Wt = Ws.sum(0)
    ref = np.full((3, St.shape[-1]), np.nan)
    for c in range(3):
        if Wt[c] > 0:
            ref[c] = (Ws[:, c] @ St[:, c]) / Wt[c]
    SO = StackObs(ref, st["tau"], st["tS"], st["sp"], np.nan)
    # ---- day table
    dd = dvv_only(SO, St)
    ds = day_scalars(st["St"], st["tau"], st["tS"], st["sp"])
    sc = src_scalars(st["St"], st["tau"], st["tS"], st["sp"])
    day = dict(date=np.asarray(st["days"]), dn=dn, lag=lj, major=jj, prev=prev.astype(float), nxt=nxt, betw=betw,
               ws=Ws.mean(1), nd=st["Nd"] if st.get("Nd") is not None else np.full(len(dn), -1), **{f"dvv_{k}": v for k, v in dd.items()},
               slope=ds["slope"], **sc)
    # ---- daily classes (A) with the joint bootstrap
    ng = 2 * nM
    rng = np.random.default_rng(20261013)
    K = rng.multinomial(nM, np.ones(nM) / nM, size=NB).astype(float)
    Kb = rng.multinomial(nM, np.ones(nM) / nM, size=NB).astype(float)
    Kall = np.r_[np.ones((1, ng)), np.c_[K, Kb]]
    res = {}
    if betw.any():
        Ps, Wg = stacks(St, Ws, betw, nM + ib, ng)
        X, w = draw_stack(Kall, Ps, Wg); o = dvv_only(SO, X); o["w"] = w; res["betw"] = o
    for d in DAYS:
        m = lj == d
        if m.sum() == 0:
            continue
        Ps, Wg = stacks(St, Ws, m, jj, ng)
        X, w = draw_stack(Kall, Ps, Wg); o = dvv_only(SO, X); o["w"] = w; res[f"D{d}"] = o
    # ---- front-scrambled daily classes (A control): shifted fronts, same draws for every series of a family
    scr = None
    alts = _G["alts"].get(fam)
    if alts is not None and NSCR > 0:
        rs = np.random.default_rng(_G["seed"] + fam)
        Xs, lab = [], []
        for k in range(NSCR):
            shift = {m: rs.choice(alts[m]) for m in M if len(alts.get(m, [])) > 0}
            fr2 = fronts_for(fam, M, _G["start"], _G["medrel"], _G["rel"], shift)
            l2 = dn[:, None] - fr2[None, :]
            j2 = np.argmin(np.where((l2 >= -10) & (l2 <= 60), np.abs(l2), 10 ** 6), 1)
            lj2 = l2[np.arange(len(dn)), j2]
            for d in DAYS:
                m = lj2 == d
                if not m.any():
                    continue
                x = np.zeros((3, St.shape[-1]))
                for c in range(3):
                    if Ws[m, c].sum() > 0:
                        x[c] = (Ws[m, c] @ St[m, c]) / Ws[m, c].sum()
                Xs.append(x); lab.append((k, d, float(Ws[m, 0].sum())))
        if Xs:
            o = dvv_only(SO, np.array(Xs)); scr = dict(lab=np.array(lab), **o)
    # ---- per-major stacks of front-lag days 0..20 (D) with a day bootstrap
    per = {}
    for mi in range(nM):
        m = (lj >= 0) & (lj <= 20) & (jj == mi)
        if m.sum() == 0:
            continue
        idx = np.flatnonzero(m)
        Kd = np.r_[np.ones((1, len(idx))), rng.multinomial(len(idx), np.ones(len(idx)) / len(idx), size=NBD)].astype(float)
        X = np.zeros((NBD + 1, 3, St.shape[-1]))
        for c in range(3):
            num = (Kd * Ws[idx, c][None, :]) @ St[idx, c]; den = Kd @ Ws[idx, c]
            X[:, c] = np.where(den[:, None] > 0, num / np.maximum(den[:, None], 1e-12), 0.0)
        o = dvv_only(SO, X); o["w"] = float(Ws[idx, 0].sum()); o["n_days"] = len(idx)
        o["nd"] = float(np.nanmean(np.where(day["nd"][idx] > 0, day["nd"][idx], np.nan))) if np.any(day["nd"][idx] > 0) else np.nan
        per[M[mi]] = o
    return key, dict(day=day, res=res, scr=scr, per=per, tS=st["tS"], sp=st["sp"])


def init(g):
    _G.update(g)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--test", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    store = load_store()
    FT, majors, start, medrel, mid, rel = front_tables()
    alts = {}
    for f, g in FT.groupby("fam"):
        rr = dict(zip(g.major, g.rel.astype(float)))
        alts[int(f)] = {m: np.array([v for mm, v in rr.items() if mm != m]) for m in majors}
    keys = sorted(store)[: a.test or None]
    log(f"[novel] region {C.REGION}: {len(keys)} series, {len(majors)} majors")
    g = dict(store=store, majors=majors, start=start, medrel=medrel, mid=mid, rel=rel, fams=set(FT.fam.astype(int)),
             alts=alts, seed=20261013)
    import multiprocessing as mp
    out = {}
    with mp.get_context("fork").Pool(a.workers, initializer=init, initargs=(g,)) as p:
        for k, (key, r) in enumerate(p.imap_unordered(series_job, keys, chunksize=2)):
            if r is not None:
                out[key] = r
            if (k + 1) % 100 == 0 or k + 1 == len(keys):
                log(f"  [novel] {k + 1}/{len(keys)} ({len(out)} usable) {time.time() - t0:.0f}s")
    tag = "_test" if a.test else ""
    pickle.dump(dict(out=out, majors=majors, mid={m: str(v.date()) for m, v in mid.items()}),
                open(f"{TB}/novel_raw{tag}.pkl", "wb"))
    open(f"{OUT}/NOVEL_{SUB.upper()}{tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    log(f"[novel] done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
