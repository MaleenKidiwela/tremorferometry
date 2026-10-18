#!/usr/bin/env python
"""s12_epoch -- S12-LVL (LFE side), tasks 3 + 4 (PLAN.md §3-§4): superposed-epoch observables aligned on the slip-front
arrival at each family, against a between-ETS baseline, from the S11 round-3 cc-weighted day stacks (store_ccw[_ext]).

Day classes per family (front = S11 round-7 headline front date of the family in that major):
  pre (-10..-1 d), 0, 1, 2, 3, 4-6, 7-10, 11-20, tail (21-60) from the nearest major front; 'betw' = more than 60 d after
  the previous and more than 30 d before the next major front (a family without a front in a major uses that major's
  median front). A second frame (Gosselin's) uses the major's regional midpoint (median front over families): G-15..-1,
  G0..+15, G16..45.
Observables per (series, class), all as contrasts against the series' own 'betw' class:
  * coda dv/v (%) in lapse windows 1-2, 2-3, 3-4, 4-6 and 2-5 s past S: class stack (Ws-weighted day stacks, all majors
    superposed) vs the series grand stack, joint 3C S-anchored stretch, bank +-3 %;
  * coda decay: per day, ln of the noise-corrected 3C energy relative to the S peak, mean over each lapse window, and the
    slope of ln E over 1-6 s (1/s; delta(1/Qc) = -delta slope / (2 pi f), f = 5 Hz); class = Ws-weighted day mean;
  * source side: ln(A_P/A_S) on Z and on H (projection of the class stack on the reference wavelet in P / S windows),
    and delta(S-P) (ms) from the P-window lag after S alignment;
  * SrS window shift: lag (ms) of the horizontal window around the predicted LVL-base reflection (nominal model) vs the
    reference.
Pooled over series with inverse-variance weights (series sigma from a 200-draw bootstrap over majors/inter-ETS intervals);
the pooled sigma is a joint bootstrap over majors (same draw for every series). Controls: front-scrambled (each family's
front replaced by start(e) + rel(f, e') of another major e'), 200 draws; time-of-day scrambled control stacks (Sc) through
the same classes.
    S11_REGION=SVI|WA_NW python s12_epoch.py [--workers 12] [--test N] [--nscr 200]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, time, glob, argparse
import numpy as np
import pandas as pd
from scipy.signal import hilbert
ROOT = "/home/jovyan/tremorferometry"
sys.path.insert(0, f"{ROOT}/svi_science/scripts"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s11_common as C
import s11_r3_common as R3
from s11_measure import Ref3

OUT = f"{ROOT}/svi_science/results/s12_lvl/lfe"
SUB = "svi" if C.REGION == "SVI" else C.REGION.lower()
TB = f"{OUT}/tables/{SUB}"; LOGF = f"{OUT}/logs/epoch_{SUB}.log"
R7 = f"{ROOT}/svi_science/results/s11_depth/round7/{SUB}/tables"
WINS = {"1-2": (1.0, 2.0), "2-3": (2.0, 3.0), "3-4": (3.0, 4.0), "4-6": (4.0, 6.0), "2-5": (2.0, 5.0)}
EPSB = np.linspace(-0.03, 0.03, 301)
CLS_F = [("pre", -10, -1), ("0", 0, 0), ("1", 1, 1), ("2", 2, 2), ("3", 3, 3), ("4-6", 4, 6), ("7-10", 7, 10),
         ("11-20", 11, 20), ("tail", 21, 60)]
CLS_G = [("G-15..-1", -15, -1), ("G0..15", 0, 15), ("G16..45", 16, 45)]
CLASSES = [c[0] for c in CLS_F] + [c[0] for c in CLS_G] + ["betw", "betw_early", "betw_late"]
NB = 200
F0 = 5.0
_G = {}


def log(m):
    s = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"
    print(s, flush=True); open(LOGF, "a").write(s + "\n")


# ------------------------------------------------------------------------------ inputs
def load_store():
    paths = [f"{R3.TB3}/store_ccw"] + ([f"{R3.TB3}/store_ccw_ext"] if C.REGION == "SVI" and os.path.isdir(f"{R3.TB3}/store_ccw_ext") else [])
    out = {}
    for d in paths:
        for p in sorted(glob.glob(f"{d}/*.npz")):
            z = np.load(p)
            o1 = np.r_[0, np.cumsum(z["n_days"])]; o2 = np.r_[0, np.cumsum(z["n_c2"])]
            for i, k in enumerate(z["keys"]):
                out[str(k)] = dict(days=z["days"][o1[i]:o1[i + 1]].astype(str), St=z["St"][o1[i]:o1[i + 1]],
                                   Ws=z["Ws"][o1[i]:o1[i + 1]], c2_day=z["c2_day"][o2[i]:o2[i + 1]],
                                   Sc=z["Sc"][o2[i]:o2[i + 1]], Cn=z["Cn"][o2[i]:o2[i + 1]], tS=float(z["tS"][i]),
                                   sp=float(z["sp"][i]), tau=z["tau"],
                                   Nd=z["Nd"][o1[i]:o1[i + 1]] if "Nd" in z.files else None)
    return out


def front_tables():
    FT = pd.read_csv(f"{R7}/r7_front_times.csv")
    FT = FT[FT.set == "headline"].copy()
    FT["front"] = pd.to_datetime(FT.front_utc, format="mixed").dt.normalize()
    FT["start"] = pd.to_datetime(FT.start)
    majors = sorted(FT.major.unique())
    start = FT.groupby("major").start.first()
    medrel = FT.groupby("major").rel.median()
    mid = {m: start[m] + pd.Timedelta(days=float(medrel[m])) for m in majors}
    rel = {(int(r.fam), r.major): float(r.rel) for r in FT.itertuples()}
    return FT, majors, start, medrel, mid, rel


def fronts_for(fam, majors, start, medrel, rel, shift=None):
    """front day numbers (days since 1970) of family fam in every major; shift: dict major -> alternative rel"""
    out = []
    for m in majors:
        r = rel.get((fam, m), float(medrel[m]))
        if shift is not None and m in shift:
            r = shift[m]
        out.append((start[m] + pd.Timedelta(days=r)).normalize().value // 86400_000_000_000)
    return np.array(out, dtype=np.int64)


def classify(dn, fr, midn):
    """dn day numbers (n,), fr front day numbers per major (M,), midn regional midpoints (M,) ->
    class index per frame: cf (n,) index into CLS_F or -1, gF (n,) major index; cg, gG; betw mask + interval index"""
    lag = dn[:, None] - fr[None, :]                          # (n, M)
    j = np.argmin(np.where((lag >= -10) & (lag <= 60), np.abs(lag), 10 ** 6), 1)
    lj = lag[np.arange(len(dn)), j]
    cf = np.full(len(dn), -1)
    for k, (nm, a, b) in enumerate(CLS_F):
        cf[(lj >= a) & (lj <= b)] = k
    lg = dn[:, None] - midn[None, :]
    jg = np.argmin(np.where((lg >= -15) & (lg <= 45), np.abs(lg), 10 ** 6), 1)
    lgj = lg[np.arange(len(dn)), jg]
    cg = np.full(len(dn), -1)
    for k, (nm, a, b) in enumerate(CLS_G):
        cg[(lgj >= a) & (lgj <= b)] = k
    prev = np.where(lag >= 0, lag, 10 ** 6).min(1)            # days since the previous front
    nxt = np.where(lag < 0, -lag, 10 ** 6).min(1)             # days to the next front
    betw = (prev > 60) & (nxt > 30)
    ib = np.searchsorted(np.sort(fr), dn) - 1                 # inter-ETS interval = index of the previous major
    classify.prev = prev
    return cf, j, cg, jg, betw, np.clip(ib, 0, len(fr) - 1)


# ------------------------------------------------------------------------------ per-day scalar observables
def day_scalars(St, tau, tS, sp):
    """per day: ln energy (rel. S peak, noise-corrected) means per lapse window, slope 1-6 s"""
    tr = tau - tS
    X = np.nan_to_num(St.astype(np.float64))
    E = np.zeros((X.shape[0], X.shape[2]))
    for c in range(3):
        E += np.abs(hilbert(X[:, c], axis=1)) ** 2
    k = np.ones(9) / 9
    E = np.apply_along_axis(lambda v: np.convolve(v, k, "same"), 1, E)
    noise = E[:, (tr >= -13.0) & (tr <= -sp - 1.0)].mean(1, keepdims=True)
    Es = E[:, (tr >= -0.2) & (tr <= 0.5)].max(1, keepdims=True)
    Ec = np.maximum(E - noise, 1e-6 * Es)
    L = np.log(Ec / Es)
    out = {}
    for wn, (a, b) in WINS.items():
        if wn == "2-5":
            continue
        out[f"lnE_{wn}"] = L[:, (tr >= a) & (tr < b)].mean(1)
    m = (tr >= 1.0) & (tr <= 6.0)
    t = tr[m]; tt = t - t.mean()
    out["slope"] = (L[:, m] - L[:, m].mean(1, keepdims=True)) @ tt / (tt @ tt)
    out["snr_late"] = np.log(E[:, (tr >= 4) & (tr < 6)].mean(1) / noise[:, 0])
    return out


# ------------------------------------------------------------------------------ stack observables
class StackObs:
    """measures one or many class stacks (m, 3, NW) against the series reference"""
    def __init__(self, ref, tau, tS, sp, srs_lag):
        self.tau, self.tS, self.sp = tau, tS, sp
        self.R = {wn: Ref3(ref, tS, w, eps=EPSB) for wn, w in WINS.items()}
        tr = tau - tS
        self.tr = tr
        tP = -sp
        self.wP = (tr >= tP - 0.15) & (tr <= tP + 0.6)
        self.wS = (tr >= -0.15) & (tr <= 0.6)
        self.ref = ref
        self.srs = None
        if np.isfinite(srs_lag) and srs_lag > 0.6:
            self.srs = (tr >= srs_lag - 0.3) & (tr <= srs_lag + 0.3)
        self.lags = np.arange(-6, 7)                       # +-0.15 s at 40 Hz

    def proj(self, X, c, w):
        r = self.ref[c][w]
        if not np.isfinite(r).all():
            return np.full(X.shape[0], np.nan)
        return (X[:, c][:, w] @ r) / (r @ r)

    def lag(self, X, comps, w):
        """cc lag (s) of X windows vs reference, joint over comps, parabolic refinement"""
        idx = np.flatnonzero(w)
        acc = np.zeros((X.shape[0], len(self.lags))); k = 0
        for c in comps:
            r = self.ref[c][idx]
            if not np.isfinite(r).all():
                continue
            r = (r - r.mean()) / np.linalg.norm(r - r.mean())
            for j, L in enumerate(self.lags):
                x = X[:, c, idx + L]
                x = x - x.mean(1, keepdims=True); n = np.linalg.norm(x, axis=1); n[n == 0] = np.nan
                acc[:, j] += (x @ r) / n
            k += 1
        if k == 0:
            return np.full(X.shape[0], np.nan)
        acc /= k
        i = np.argmax(np.nan_to_num(acc, nan=-9), 1); ii = np.clip(i, 1, len(self.lags) - 2); rr = np.arange(len(i))
        y0, y1, y2 = acc[rr, ii - 1], acc[rr, ii], acc[rr, ii + 1]; d2 = y0 - 2 * y1 + y2
        dl = np.where((d2 < 0) & (i == ii), 0.5 * (y0 - y2) / np.where(d2 == 0, -1, d2), 0.0)
        return (self.lags[i] + np.clip(dl, -1, 1)) / 40.0

    def measure(self, X):
        X = np.nan_to_num(X)
        o = {}
        for wn, R in self.R.items():
            e, cc = R.measure(X, (0, 1, 2))
            e = np.where(np.abs(e) >= 0.0299, np.nan, e)
            o[f"dvv_{wn}"] = -100 * e; o[f"cc_{wn}"] = cc
        aPZ, aSZ = self.proj(X, 0, self.wP), self.proj(X, 0, self.wS)
        aPH = np.nansum([np.abs(self.proj(X, c, self.wP)) for c in (1, 2)], 0)
        aSH = np.nansum([np.abs(self.proj(X, c, self.wS)) for c in (1, 2)], 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            o["lnPS_Z"] = 100 * np.log(np.abs(aPZ) / np.abs(aSZ))
            o["lnPS_H"] = 100 * np.log(aPH / aSH)
        o["dSP_ms"] = -1000 * self.lag(X, (0,), self.wP)         # P later after S alignment -> S-P shorter
        o["dSrS_ms"] = 1000 * self.lag(X, (1, 2), self.srs) if self.srs is not None else np.full(X.shape[0], np.nan)
        return o


# ------------------------------------------------------------------------------ one series
def series_job(key):
    st = _G["store"][key]; fam = int(key.split("|")[0])
    if fam not in _G["fams_with_front"]:
        return key, None
    days = pd.to_datetime(pd.Series(st["days"])).dt.normalize()
    dn = (days.astype("datetime64[ns]").astype(np.int64) // 86400_000_000_000).to_numpy()
    M = _G["majors"]
    fr = fronts_for(fam, M, _G["start"], _G["medrel"], _G["rel"])
    midn = np.array([_G["mid"][m].normalize().value // 86400_000_000_000 for m in M], dtype=np.int64)
    cf, jf, cg, jg, betw, ib = classify(dn, fr, midn)
    prev_days = classify.prev
    St = st["St"].astype(np.float64); Ws = st["Ws"]
    St = np.nan_to_num(St)
    Wt = Ws.sum(0)
    if np.any(Wt <= 0):
        pass
    ref = np.full((3, St.shape[-1]), np.nan)
    for c in range(3):
        if Wt[c] > 0:
            ref[c] = (Ws[:, c] @ St[:, c]) / Wt[c]
    srs = _G["srs"].get((fam, key.split("|")[1]), np.nan)
    SO = StackObs(ref, st["tau"], st["tS"], st["sp"], srs)
    nM = len(M)
    # group label per day for the bootstrap: ETS classes -> major index, betw -> nM + inter-ETS interval
    cls_list = []
    for k in range(len(CLS_F)):
        cls_list.append((CLASSES[k], cf == k, jf))
    for k in range(len(CLS_G)):
        cls_list.append((CLASSES[len(CLS_F) + k], cg == k, jg))
    cls_list.append(("betw", betw, nM + ib))
    prev = prev_days
    cls_list.append(("betw_early", betw & (prev <= 240), nM + ib))
    cls_list.append(("betw_late", betw & (prev > 240), nM + ib))
    ng = 2 * nM
    # partial sums per (class, group)
    P = {}; Wp = {}
    for nm, msk, grp in cls_list:
        if not msk.any():
            continue
        g = grp[msk]
        Ps = np.zeros((ng, 3, St.shape[-1])); Wg = np.zeros((ng, 3))
        for c in range(3):
            np.add.at(Ps[:, c], g, Ws[msk, c][:, None] * St[msk, c])
            np.add.at(Wg[:, c], g, Ws[msk, c])
        P[nm] = Ps; Wp[nm] = Wg
    if "betw" not in P:
        return key, None
    rng = np.random.default_rng(20261013)               # same draws for every series -> joint bootstrap over majors
    K = rng.multinomial(nM, np.ones(nM) / nM, size=NB).astype(float)       # majors resampled
    Kb = rng.multinomial(nM, np.ones(nM) / nM, size=NB).astype(float)      # inter-ETS intervals resampled
    Kall = np.c_[K, Kb]                                                    # (NB, ng)
    Kall = np.r_[np.ones((1, ng)), Kall]                                   # draw 0 = the data
    res = {}
    for nm in P:
        num = np.einsum("bg,gct->bct", Kall, P[nm]); den = Kall @ Wp[nm]
        X = np.where(den[:, :, None] > 0, num / np.maximum(den[:, :, None], 1e-12), 0.0)
        o = SO.measure(X)
        o["w"] = den[:, 0]
        res[nm] = o
    # day scalars (incoherent): class value = Ws-weighted mean over days of the class (bootstrap with the same groups)
    ds = day_scalars(st["St"], st["tau"], st["tS"], st["sp"])
    wd = Ws.mean(1)
    for nm, msk, grp in cls_list:
        if nm not in res:
            continue
        g = grp[msk]
        for q, v in ds.items():
            vv = np.nan_to_num(v[msk]); ww = wd[msk] * np.isfinite(v[msk])
            sv = np.bincount(g, ww * vv, ng); sw = np.bincount(g, ww, ng)
            num = Kall @ sv; den = Kall @ sw
            res[nm][f"day_{q}"] = np.where(den > 0, num / np.maximum(den, 1e-12), np.nan)
    # front-scrambled: reclassify ETS days with shifted fronts; stacks of the ETS classes only (betw unchanged)
    scr = {}
    alts = _G["alts"].get(fam)
    if alts is not None and _G["nscr"] > 0:
        rs = np.random.default_rng(_G["seed"] + fam)          # same draws for every series of a family
        Xs = []; lab = []
        for d in range(_G["nscr"]):
            shift = {m: rs.choice(alts[m]) for m in M if len(alts.get(m, [])) > 0}
            fr2 = fronts_for(fam, M, _G["start"], _G["medrel"], _G["rel"], shift)
            cf2, jf2, _, _, _, _ = classify(dn, fr2, midn)
            for nm in ("pre", "0..20"):
                msk = (cf2 == 0) if nm == "pre" else ((cf2 >= 1) & (cf2 <= 7))
                if not msk.any():
                    Xs.append(np.zeros((3, St.shape[-1]))); lab.append((d, nm, 0.0)); continue
                w = Ws[msk]; x = np.zeros((3, St.shape[-1]))
                for c in range(3):
                    if w[:, c].sum() > 0:
                        x[c] = (w[:, c] @ St[msk, c]) / w[:, c].sum()
                Xs.append(x); lab.append((d, nm, float(w[:, 0].sum())))
        o = SO.measure(np.array(Xs))
        scr = dict(lab=lab, **{k: v for k, v in o.items() if k.startswith(("dvv_", "lnPS", "dSP", "dSrS"))})
    # aggregate data classes "0..20" (front lags 0..20) for the headline contrast
    msk = (cf >= 1) & (cf <= 7); n020 = int(msk.sum())
    if msk.any():
        nm = "0..20"; g = jf[msk]
        Ps = np.zeros((ng, 3, St.shape[-1])); Wg = np.zeros((ng, 3))
        for c in range(3):
            np.add.at(Ps[:, c], g, Ws[msk, c][:, None] * St[msk, c]); np.add.at(Wg[:, c], g, Ws[msk, c])
        num = np.einsum("bg,gct->bct", Kall, Ps); den = Kall @ Wg
        X = np.where(den[:, :, None] > 0, num / np.maximum(den[:, :, None], 1e-12), 0.0)
        res[nm] = SO.measure(X); res[nm]["w"] = den[:, 0]
        for q, v in ds.items():
            vv = np.nan_to_num(v[msk]); ww = wd[msk] * np.isfinite(v[msk])
            sv = np.bincount(g, ww * vv, ng); sw = np.bincount(g, ww, ng)
            res[nm][f"day_{q}"] = np.where(Kall @ sw > 0, (Kall @ sv) / np.maximum(Kall @ sw, 1e-12), np.nan)
    # time-of-day scrambled control stacks through the same classes (draw 0 only), against the data reference
    ctrl = {}
    if len(st["c2_day"]):
        dc = dn[st["c2_day"]]
        cfc, _, _, _, betc, _ = classify(dc, fr, midn)
        Sc = np.nan_to_num(st["Sc"].astype(np.float64)); Cn = st["Cn"]
        for nm, msk in (("betw", betc), ("0..20", (cfc >= 1) & (cfc <= 7)), ("pre", cfc == 0)):
            if msk.sum() < 2:
                continue
            x = np.zeros((3, Sc.shape[-1]))
            for c in range(3):
                if Cn[msk, c].sum() > 0:
                    x[c] = (Cn[msk, c] @ Sc[msk, c]) / Cn[msk, c].sum()
            o = SO.measure(x[None])
            ctrl[nm] = {k: float(v[0]) for k, v in o.items()}
    meta = dict(n_days=len(dn), n_betw=int(betw.sum()), n_0_20=n020, n_pre=int((cf == 0).sum()),
                tS=st["tS"], sp=st["sp"])
    return key, dict(res=res, scr=scr, ctrl=ctrl, meta=meta)


def init(g):
    _G.update(g)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--test", type=int, default=0); ap.add_argument("--nscr", type=int, default=200)
    a = ap.parse_args()
    os.makedirs(TB, exist_ok=True)
    t0 = time.time()
    store = load_store()
    FT, majors, start, medrel, mid, rel = front_tables()
    fams = set(FT.fam.astype(int))
    alts = {}
    for f, g in FT.groupby("fam"):
        rr = dict(zip(g.major, g.rel.astype(float)))
        alts[int(f)] = {m: np.array([v for mm, v in rr.items() if mm != m]) for m in majors}
    srs = {}
    if C.REGION == "SVI" and os.path.exists(f"{OUT}/tables/fwd_rays.csv"):
        Rr = pd.read_csv(f"{OUT}/tables/fwd_rays.csv")
        srs = {(int(r.fam), r.netsta): float(getattr(r, "lagS_SrS_b_H3_5_v3", np.nan)) for r in
               Rr.rename(columns=lambda c: c.replace(".", "_")).itertuples()}
    keys = sorted(store)
    if a.test:
        keys = keys[:a.test]
    log(f"region {C.REGION}: {len(keys)} series, {len(majors)} majors, {len(fams)} families with fronts; nscr {a.nscr}")
    g = dict(store=store, majors=majors, start=start, medrel=medrel, mid=mid, rel=rel, fams_with_front=fams, alts=alts,
             srs=srs, nscr=a.nscr, seed=20261013)
    import multiprocessing as mp
    out = {}
    with mp.get_context("fork").Pool(a.workers, initializer=init, initargs=(g,)) as p:
        for k, (key, r) in enumerate(p.imap_unordered(series_job, keys, chunksize=2)):
            if r is not None:
                out[key] = r
            if (k + 1) % 50 == 0 or k + 1 == len(keys):
                log(f"  {k + 1}/{len(keys)} series ({len(out)} usable) {time.time() - t0:.0f}s")
    tag = "_test" if a.test else ""
    import pickle
    pickle.dump(out, open(f"{TB}/epoch_raw{tag}.pkl", "wb"))
    open(f"{OUT}/EPOCH_{SUB.upper()}{tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    log(f"done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
