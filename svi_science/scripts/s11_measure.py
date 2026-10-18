#!/usr/bin/env python
"""s11_measure -- S11 step 4, raw path (D3) + day-stack store for the SVD path and the injection suite.
Rules: svi_science/results/s11_depth/PLAN.md. Per (family, station, response epoch) series, the S10-multi path unchanged up
to the realigned windows (s10_measure.realign, pick_S, wmean, boot_stacks imported): grand stack, S re-pick, per-event 3C
realignment, S re-pick. Then:
  * P-weighted 3C DAY stacks (data) and equal-weight scrambled-control (c2) day stacks (<= n data windows per day) are
    stored per station in tables/store/<NET.STA>.npz (the SVD path and the waveform-level injections start from them);
  * per EPISODE: P-weighted 3C stack, leave-episode-out reference, JOINT 3C S-anchored stretch (mean of per-component cc
    = cc of the concatenated z-scored windows) in 2-5 (headline), 3-6, 2-3.5, 3.5-6 s past S; Z-only and H-only at 2-5 s;
    2-parameter (shift +-50 ms x stretch) at 2-5 s; 200-draw event bootstrap sigma; pre-arrival RMS; c2 control episode
    stack through the identical path. Rail |eps| >= 0.99 % -> NaN (kept as *_raw). eps = -dv/v.

    python s11_measure.py [--workers 12] [--stations NET.STA,...] [--fams 10,87] [--tag _test]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, zlib
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s11_common as C
from s10_measure import realign, pick_S, wmean, boot_stacks, MIR0

EPS = np.linspace(-0.01, 0.01, 201)
DE = EPS[1] - EPS[0]
SHIFTS = np.arange(-0.050, 0.0501, 0.005)
EPS2 = np.linspace(-0.01, 0.01, 101)
NBOOT = 200
COMPS = {"3c": (0, 1, 2), "Z": (0,), "H": (1, 2)}


# ------------------------------------------------------------------------------ joint multi-component stretch
class Ref3:
    """S-anchored stretch bank of a 3C reference over one window: ref(tau_rel/(1+eps)) <-> eps = -dv/v"""
    def __init__(self, ref, tS, win, shifts=None, eps=EPS):
        from scipy.interpolate import CubicSpline
        tr = S.TAU - tS
        self.idx = np.flatnonzero((tr >= win[0] - 1e-9) & (tr <= win[1] + 1e-9))
        x = tr[self.idx]
        self.eps = eps
        self.ok = np.array([np.isfinite(ref[c]).all() for c in range(3)])
        self.B = []
        for c in range(3):
            if not self.ok[c]:
                self.B.append(None); continue
            sp = CubicSpline(tr, ref[c])
            if shifts is None:
                arg = x[None, :] / (1.0 + eps[:, None])
            else:
                arg = (x[None, None, :] / (1.0 + eps[None, :, None]) - shifts[:, None, None]).reshape(-1, len(x))
            B = sp(arg); B = B - B.mean(1, keepdims=True)
            self.B.append(B / np.linalg.norm(B, axis=1, keepdims=True))

    def curves(self, X, comps):
        """X (m, 3, NW) -> joint cc curves (m, n_bank) averaged over usable components"""
        X = np.atleast_3d(X) if X.ndim == 3 else X[None]
        acc = None; k = 0
        for c in comps:
            if not self.ok[c]:
                continue
            x = X[:, c, self.idx].astype(np.float64)
            if not np.isfinite(x).all():
                continue
            x = x - x.mean(1, keepdims=True); nx = np.linalg.norm(x, axis=1, keepdims=True); nx[nx == 0] = np.nan
            cc = (x / nx) @ self.B[c].T
            acc = cc if acc is None else acc + cc; k += 1
        return None if acc is None else acc / k

    def measure(self, X, comps=(0, 1, 2)):
        cc = self.curves(X, comps)
        if cc is None:
            m = (X if X.ndim == 3 else X[None]).shape[0]
            return np.full(m, np.nan), np.full(m, np.nan)
        cc = np.nan_to_num(cc, nan=-9)
        i = np.argmax(cc, 1); ii = np.clip(i, 1, len(self.eps) - 2); r = np.arange(len(i))
        y0, y1, y2 = cc[r, ii - 1], cc[r, ii], cc[r, ii + 1]; d2 = y0 - 2 * y1 + y2
        delta = np.where((d2 < 0) & (i == ii), 0.5 * (y0 - y2) / np.where(d2 == 0, -1, d2), 0.0)
        de = self.eps[1] - self.eps[0]
        e = self.eps[i] + np.clip(delta, -1, 1) * de
        c = cc[r, i]
        c = np.where(c < -1, np.nan, c)
        return e, c


def measure_2p(ref, tS, X):
    R = Ref3(ref, tS, C.WINS[C.HEAD], shifts=SHIFTS, eps=EPS2)
    cc = R.curves(X[None], (0, 1, 2))
    if cc is None:
        return np.nan, np.nan, np.nan
    cc = cc[0].reshape(len(SHIFTS), len(EPS2)); i, j = np.unravel_index(np.nanargmax(cc), cc.shape)
    return float(EPS2[j]), float(SHIFTS[i]), float(cc[i, j])


def stack3(W, w):
    """W (n, 3, NW) with NaN components, weights w -> (3, NW) weighted mean, (3,) weight sums"""
    out = np.full((3, S.NW), np.nan); ws = np.zeros(3)
    for c in range(3):
        ok = np.isfinite(W[:, c, 0])
        if ok.any():
            ws[c] = w[ok].sum(); out[c] = (w[ok] @ np.nan_to_num(W[ok, c].astype(np.float64))) / ws[c]
    return out, ws


def boot3(W, w, nb, rng):
    """(nb, 3, NW) bootstrap weighted 3C stacks (same draw for all components)"""
    n = W.shape[0]
    cnt = np.stack([np.bincount(rng.integers(0, n, n), minlength=n) for _ in range(nb)]).astype(np.float64)
    out = np.full((nb, 3, S.NW), np.nan)
    for c in range(3):
        ok = np.isfinite(W[:, c, 0])
        if ok.sum() < 2:
            continue
        ww = cnt[:, ok] * w[None, ok]; den = ww.sum(1, keepdims=True); den[den == 0] = np.nan
        out[:, c] = (ww @ np.nan_to_num(W[ok, c].astype(np.float64))) / den
    return out


def measure_cell(Sd, ref, tS, B=None, control=None, Bc=None):
    """all columns for one episode stack Sd (3,NW) against ref (3,NW); B bootstrap stacks"""
    r = {}
    for wn, win in C.WINS.items():
        R = Ref3(ref, tS, win)
        for comp, cs in COMPS.items():
            if comp != "3c" and wn != C.HEAD:
                continue
            e, c = R.measure(Sd[None], cs)
            r[f"eps{comp}_{wn}"], r[f"cc{comp}_{wn}"] = float(e[0]), float(c[0])
            if B is not None:
                eb = R.measure(B, cs)[0]; eb = eb[np.isfinite(eb) & (np.abs(eb) < C.RAIL)]
                r[f"sig{comp}_{wn}"] = float(np.std(eb, ddof=1)) if len(eb) > 10 else np.nan
            if control is not None and comp == "3c":
                e, c = R.measure(control[None], cs)
                r[f"c2_eps{comp}_{wn}"], r[f"c2_cc{comp}_{wn}"] = float(e[0]), float(c[0])
                if Bc is not None and wn == C.HEAD:
                    eb = R.measure(Bc, cs)[0]; eb = eb[np.isfinite(eb) & (np.abs(eb) < C.RAIL)]
                    r[f"c2_sig{comp}_{wn}"] = float(np.std(eb, ddof=1)) if len(eb) > 10 else np.nan
    return r


# ------------------------------------------------------------------------------ one series
def run_series(fam, ns, ep, Xd, Wd, sp, X2, W2, unit, do_boot=True):
    rng = np.random.default_rng(C.SEED + fam * 1000 + zlib.crc32(f'{ns}|{ep}'.encode()) % 997)
    okZ = np.isfinite(Wd[:, 0, 0])
    Xd, Wd = Xd[okZ].reset_index(drop=True), Wd[okZ]
    out = dict(fam=fam, netsta=ns, epoch=ep, n_data=len(Xd), n_c2=len(X2))
    if len(Xd) < 6 or Xd.date.nunique() < 3:
        return out, None
    P = Xd.P.to_numpy()
    G0 = np.stack([wmean(Wd[:, c], P)[0] for c in range(3)])
    tS0 = pick_S(G0, 0.0)
    Wd, lag_d, cc_d = realign(Wd, G0, tS0)
    ok2 = np.isfinite(W2[:, 0, 0]) if len(X2) else np.zeros(0, bool)
    X2, W2 = X2[ok2].reset_index(drop=True), W2[ok2]
    if len(X2):
        W2 = realign(W2, G0, tS0)[0]
    G1 = np.stack([wmean(Wd[:, c], P)[0] for c in range(3)])
    tS = pick_S(G1, tS0)
    SP = sp if np.isfinite(sp) else 6.0
    mend = -(SP + 1.0)
    tr = S.TAU - tS
    mir = (tr >= MIR0) & (S.TAU >= -13.5) & (tr <= mend)
    cw = (tr >= 2.0) & (tr <= 5.0)
    rms = lambda x, m: float(np.sqrt(np.mean(x[m] ** 2)))
    out.update(tS_pick=tS, sp=SP, causal25=rms(G1[0], cw) / rms(G1[0], mir),
               causal25_3c=float(np.sqrt(np.nanmean([rms(G1[c], cw) ** 2 / rms(G1[c], mir) ** 2 for c in range(3)
                                                    if np.isfinite(G1[c]).all()]))),
               realign_cc_med=float(np.median(cc_d)), lag_med_abs=float(np.median(np.abs(lag_d))),
               n_days=int(Xd.date.nunique()), horiz_ok=bool(np.isfinite(G1[1]).all() and np.isfinite(G1[2]).all()))
    # ---------------- day stacks (store)
    days = np.array(sorted(Xd.date.unique()))
    dix = [np.flatnonzero((Xd.date == d).to_numpy()) for d in days]
    St = np.zeros((len(days), 3, S.NW), np.float32); Ws = np.zeros((len(days), 3))
    for k, ix in enumerate(dix):
        s_, w_ = stack3(Wd[ix], P[ix]); St[k] = np.nan_to_num(s_); Ws[k] = w_
    c2_rows = []; Sc = []; Cn = []
    if len(X2):
        for k, d in enumerate(days):
            mc = np.flatnonzero((X2.date == d).to_numpy())
            if len(mc) > len(dix[k]):
                mc = np.sort(rng.choice(mc, len(dix[k]), replace=False))
            if len(mc) == 0:
                continue
            s_, w_ = stack3(W2[mc], np.ones(len(mc)))
            Sc.append(np.nan_to_num(s_).astype(np.float32)); Cn.append(w_); c2_rows.append(k)
    store = dict(days=days, St=St, Ws=Ws, c2_day=np.array(c2_rows, int), Sc=np.array(Sc, np.float32).reshape(-1, 3, S.NW),
                 Cn=np.array(Cn).reshape(-1, 3), G1=G1.astype(np.float32), tS=tS, sp=SP)
    # ---------------- episode stacks
    ep_row = Xd.date.map(unit).to_numpy()
    ep2 = X2.date.map(unit).to_numpy() if len(X2) else np.zeros(0, object)
    Gs = np.stack([np.nansum(Ws[:, c, None] * St[:, c], 0) for c in range(3)]); Wt = Ws.sum(0)
    pre = (S.TAU >= -13.5) & (tr <= mend)
    rows = []
    for e_ in pd.unique(ep_row):
        if not isinstance(e_, str):
            continue
        m = ep_row == e_
        ixe = np.flatnonzero(m)
        if len(ixe) < 3:
            continue
        Sd, ws = stack3(Wd[ixe], P[ixe])
        if np.any((Wt - ws)[np.isfinite(G1).all(1)] <= 0):
            continue
        ref = np.stack([(Gs[c] - ws[c] * np.nan_to_num(Sd[c])) / (Wt[c] - ws[c]) if Wt[c] - ws[c] > 0 else np.full(S.NW, np.nan)
                        for c in range(3)])
        ref[~np.isfinite(G1).all(1)] = np.nan
        r = dict(fam=fam, netsta=ns, epoch=ep, episode=e_, n=int(len(ixe)), sumP=float(P[ixe].sum()),
                 n_days=int(Xd.date[m].nunique()), date_mid=str(sorted(Xd.date[m])[len(ixe) // 2]))
        B = boot3(Wd[ixe], P[ixe], NBOOT, rng) if do_boot else None
        ctrl = Bc = None
        mc = np.flatnonzero(ep2 == e_) if len(X2) else np.zeros(0, int)
        if len(mc) > len(ixe):
            mc = np.sort(rng.choice(mc, len(ixe), replace=False))
        r["n_ctrl2"] = int(len(mc))
        if len(mc) >= 3:
            ctrl = stack3(W2[mc], np.ones(len(mc)))[0]
            Bc = boot3(W2[mc], np.ones(len(mc)), NBOOT, rng) if do_boot else None
        r.update(measure_cell(Sd, ref, tS, B, ctrl, Bc))
        e2, s2, c2 = measure_2p(ref, tS, Sd)
        r.update(eps2p_25=e2, shift2p_25=s2, cc2p_25=c2)
        r["pre_rms_ratio"] = rms(Sd[0], pre) / max(rms(Sd[0], cw), 1e-30)
        rows.append(r)
    E = pd.DataFrame(rows)
    if len(E):
        E["pre_rel"] = E.pre_rms_ratio / E.pre_rms_ratio.median()
        for col in [c for c in E.columns if c.startswith(("eps", "c2_eps"))]:
            E[f"{col}_raw"] = E[col]; rl = E[col].abs() >= C.RAIL
            E.loc[rl, col] = np.nan
            sc = col.replace("eps", "sig")
            if sc in E:
                E.loc[rl, sc] = np.nan
    return out, (E, store)


def _run(args):
    try:
        return run_series(*args)
    except Exception as e:
        import traceback
        return dict(fam=args[0], netsta=args[1], epoch=args[2], n_data=-1, error=f"{e} {traceback.format_exc()[-800:]}"), None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--stations", default="")
    ap.add_argument("--fams", default=""); ap.add_argument("--tag", default="")
    ap.add_argument("--src", default="", help="alternative input folder with rows/row_station/offsets_sp/active_days + win/ "
                                              "(Revision 2: tables/ext for the added families)")
    a = ap.parse_args()
    C.ensure(); os.makedirs(f"{C.TB}/store{a.tag}", exist_ok=True)
    log = C.Log(f"{C.LOGD}/measure{a.tag}.log")
    TBI, WIN = (a.src, f"{a.src}/win") if a.src else (C.MTB, C.MWIN)
    X = pd.read_parquet(f"{TBI}/rows.parquet")
    RS = pd.read_parquet(f"{TBI}/row_station.parquet")
    O = pd.read_csv(f"{TBI}/offsets_sp.csv").set_index(["fam", "netsta"])
    A = pd.read_csv(f"{TBI}/active_days.csv")
    unit = dict(zip(A.date, A.episode))
    stas = a.stations.split(",") if a.stations else C.STATIONS
    fams = [int(x) for x in a.fams.split(",")] if a.fams else None
    import multiprocessing as mp
    res_out, res_E = [], []
    t0 = time.time()
    for ns in stas:
        if not os.path.exists(f"{WIN}/{ns}.npz"):
            log(f"{ns}: no windows"); continue
        z = np.load(f"{WIN}/{ns}.npz"); Wall = z["W"]; pos = pd.Series(np.arange(len(z["row_id"])), index=z["row_id"])
        p2 = f"{WIN}/{ns}__c2.npz"
        if os.path.exists(p2):
            z2 = np.load(p2); W2all = z2["W"]; pos2 = pd.Series(np.arange(len(z2["row_id"])), index=z2["row_id"])
        else:
            W2all = np.zeros((0, 3, S.NW), np.float32); pos2 = pd.Series(dtype=int)
        R = RS[RS.netsta == ns].merge(X, on=["row_id", "date"])
        if fams:
            R = R[R.fam.isin(fams)]
        jobs = []
        for (fam, ep), g in R.groupby(["fam", "epoch"]):
            gd = g[(g.series == "data") & g.row_id.isin(pos.index)].reset_index(drop=True)
            gc = g[(g.series == "c1") & g.row_id.isin(pos2.index)].reset_index(drop=True)
            jobs.append((int(fam), ns, ep, gd, Wall[pos[gd.row_id].to_numpy()], float(O.sp.get((fam, ns), np.nan)),
                         gc, W2all[pos2[gc.row_id].to_numpy()] if len(gc) else np.zeros((0, 3, S.NW), np.float32), unit))
        del Wall, W2all
        jobs.sort(key=lambda j: -len(j[3]))
        log(f"{ns}: {len(jobs)} series")
        st_all = {}
        with mp.get_context("fork").Pool(a.workers, maxtasksperchild=10) as p:
            for out, res in p.imap_unordered(_run, jobs):
                res_out.append(out)
                if "error" in out:
                    log(f"ERROR {out['fam']} {out['netsta']} {out['epoch']}: {out['error']}")
                if res is not None:
                    E, st = res; res_E.append(E); st_all[f"{out['fam']}|{ns}|{out['epoch']}"] = st
                log(f"  {out['fam']} {ns} {out['epoch']}: n {out['n_data']} days {out.get('n_days', 0)} episodes "
                    f"{0 if res is None else len(res[0])} causal25 {out.get('causal25', np.nan):.2f} {time.time() - t0:.0f}s")
        del jobs
        keys = sorted(st_all)
        np.savez(f"{C.TB}/store{a.tag}/{ns}.npz", keys=np.array(keys), tau=S.TAU,
                 tS=np.array([st_all[k]["tS"] for k in keys]), sp=np.array([st_all[k]["sp"] for k in keys]),
                 G1=np.array([st_all[k]["G1"] for k in keys]).reshape(-1, 3, S.NW),
                 n_days=np.array([len(st_all[k]["days"]) for k in keys]), n_c2=np.array([len(st_all[k]["c2_day"]) for k in keys]),
                 days=np.concatenate([st_all[k]["days"].astype(str) for k in keys]) if keys else np.zeros(0, str),
                 St=np.concatenate([st_all[k]["St"] for k in keys]) if keys else np.zeros((0, 3, S.NW), np.float32),
                 Ws=np.concatenate([st_all[k]["Ws"] for k in keys]) if keys else np.zeros((0, 3)),
                 c2_day=np.concatenate([st_all[k]["c2_day"] for k in keys]) if keys else np.zeros(0, int),
                 Sc=np.concatenate([st_all[k]["Sc"] for k in keys]) if keys else np.zeros((0, 3, S.NW), np.float32),
                 Cn=np.concatenate([st_all[k]["Cn"] for k in keys]) if keys else np.zeros((0, 3)))
        log(f"{ns}: store written ({len(keys)} series)")
        del st_all
    sfx = a.tag
    pd.DataFrame(res_out).to_csv(f"{C.TB}/series{sfx}.csv", index=False)
    pd.concat(res_E, ignore_index=True).to_parquet(f"{C.TB}/episode_raw{sfx}.parquet", index=False)
    log(f"measure done {time.time() - t0:.0f}s, series {len(res_out)}")
    open(f"{C.OUT}/MEASURE{sfx}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
