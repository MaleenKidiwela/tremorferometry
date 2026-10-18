#!/usr/bin/env python
"""s11_dist_measure -- S11 distant-station lever, step 2b (distant/PLAN.md §4): the s11_measure raw path, unchanged, on the
far windows (grid -28..+10 s). The ONLY change from s11_measure.run_series: the pre-P (causality / pre-arrival RMS)
window is expressed relative to the grid start ([tS - 26.5 s, tS - S-P - 1 s] and tau >= T0W + 0.5) instead of the near
grid's hard-coded -12.5 / -13.5 s, which hold no pre-P noise at 100-150 km. S-P per pair from the 1-D model
(far_offsets.csv). Writes distant/<region>/tables/{series.csv, episode_raw.parquet, store/<NET.STA>.npz}. MEASURE_DONE.
    S11_REGION=SVI python s11_dist_measure.py --workers 10
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
import s11_dist_common as D
D.set_grid(D.FAR_GRID)
from s10_measure import realign, pick_S, wmean
import s11_measure as M

MIR0 = S.T0W + 1.5          # -26.5 s (near grid: -12.5)
PRE0 = S.T0W + 0.5          # -27.5 s (near grid: -13.5)


def unit_map():
    if D.REG == "SVI":
        import s11_paths as P
        return P.unit_map()
    A = pd.read_csv(f"{C.MTB}/active_days.csv")
    return dict(zip(A.date, A.episode))


def run_series(fam, ns, ep, Xd, Wd, sp, X2, W2, unit, do_boot=True):
    """= s11_measure.run_series with grid-relative pre-P windows (MIR0, PRE0)"""
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
    SP = sp if np.isfinite(sp) else 14.0
    mend = -(SP + 1.0)
    tr = S.TAU - tS
    mir = (tr >= MIR0) & (S.TAU >= PRE0) & (tr <= mend)
    cw = (tr >= 2.0) & (tr <= 5.0)
    rms = lambda x, m: float(np.sqrt(np.mean(x[m] ** 2))) if m.any() else np.nan
    out.update(tS_pick=tS, sp=SP, causal25=rms(G1[0], cw) / rms(G1[0], mir),
               causal25_3c=float(np.sqrt(np.nanmean([rms(G1[c], cw) ** 2 / rms(G1[c], mir) ** 2 for c in range(3)
                                                    if np.isfinite(G1[c]).all()]))),
               realign_cc_med=float(np.median(cc_d)), lag_med_abs=float(np.median(np.abs(lag_d))),
               n_days=int(Xd.date.nunique()), horiz_ok=bool(np.isfinite(G1[1]).all() and np.isfinite(G1[2]).all()),
               n_mir=int(mir.sum()))
    days = np.array(sorted(Xd.date.unique()))
    dix = [np.flatnonzero((Xd.date == d).to_numpy()) for d in days]
    St = np.zeros((len(days), 3, S.NW), np.float32); Ws = np.zeros((len(days), 3))
    for k, ix in enumerate(dix):
        s_, w_ = M.stack3(Wd[ix], P[ix]); St[k] = np.nan_to_num(s_); Ws[k] = w_
    c2_rows = []; Sc = []; Cn = []
    if len(X2):
        for k, d in enumerate(days):
            mc = np.flatnonzero((X2.date == d).to_numpy())
            if len(mc) > len(dix[k]):
                mc = np.sort(rng.choice(mc, len(dix[k]), replace=False))
            if len(mc) == 0:
                continue
            s_, w_ = M.stack3(W2[mc], np.ones(len(mc)))
            Sc.append(np.nan_to_num(s_).astype(np.float32)); Cn.append(w_); c2_rows.append(k)
    store = dict(days=days, St=St, Ws=Ws, c2_day=np.array(c2_rows, int), Sc=np.array(Sc, np.float32).reshape(-1, 3, S.NW),
                 Cn=np.array(Cn).reshape(-1, 3), G1=G1.astype(np.float32), tS=tS, sp=SP)
    ep_row = Xd.date.map(unit).to_numpy()
    ep2 = X2.date.map(unit).to_numpy() if len(X2) else np.zeros(0, object)
    Gs = np.stack([np.nansum(Ws[:, c, None] * St[:, c], 0) for c in range(3)]); Wt = Ws.sum(0)
    pre = (S.TAU >= PRE0) & (tr <= mend)
    rows = []
    for e_ in pd.unique(ep_row):
        if not isinstance(e_, str):
            continue
        m = ep_row == e_
        ixe = np.flatnonzero(m)
        if len(ixe) < 3:
            continue
        Sd, ws = M.stack3(Wd[ixe], P[ixe])
        if np.any((Wt - ws)[np.isfinite(G1).all(1)] <= 0):
            continue
        ref = np.stack([(Gs[c] - ws[c] * np.nan_to_num(Sd[c])) / (Wt[c] - ws[c]) if Wt[c] - ws[c] > 0 else np.full(S.NW, np.nan)
                        for c in range(3)])
        ref[~np.isfinite(G1).all(1)] = np.nan
        r = dict(fam=fam, netsta=ns, epoch=ep, episode=e_, n=int(len(ixe)), sumP=float(P[ixe].sum()),
                 n_days=int(Xd.date[m].nunique()), date_mid=str(sorted(Xd.date[m])[len(ixe) // 2]))
        B = M.boot3(Wd[ixe], P[ixe], M.NBOOT, rng) if do_boot else None
        ctrl = Bc = None
        mc = np.flatnonzero(ep2 == e_) if len(X2) else np.zeros(0, int)
        if len(mc) > len(ixe):
            mc = np.sort(rng.choice(mc, len(ixe), replace=False))
        r["n_ctrl2"] = int(len(mc))
        if len(mc) >= 3:
            ctrl = M.stack3(W2[mc], np.ones(len(mc)))[0]
            Bc = M.boot3(W2[mc], np.ones(len(mc)), M.NBOOT, rng) if do_boot else None
        r.update(M.measure_cell(Sd, ref, tS, B, ctrl, Bc))
        e2, s2, c2 = M.measure_2p(ref, tS, Sd)
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
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=10); ap.add_argument("--stations", default="")
    ap.add_argument("--fams", default=""); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    D.ensure(); os.makedirs(f"{D.DTB}/store{a.tag}", exist_ok=True)
    log = C.Log(f"{D.DLOG}/measure{a.tag}.log")
    RS = pd.read_parquet(f"{D.DTB}/far_rows.parquet")
    O = pd.read_csv(f"{D.DTB}/far_offsets.csv").set_index(["fam", "netsta"])
    unit = unit_map()
    stas = a.stations.split(",") if a.stations else sorted(RS.netsta.unique())
    fams = [int(x) for x in a.fams.split(",")] if a.fams else None
    import multiprocessing as mp
    res_out, res_E = [], []
    t0 = time.time()
    for ns in stas:
        p1 = f"{D.DTB}/win/{ns}.npz"
        if not os.path.exists(p1):
            log(f"{ns}: no windows"); continue
        z = np.load(p1); Wall = z["W"]; pos = pd.Series(np.arange(len(z["row_id"])), index=z["row_id"])
        assert len(z["tau"]) == S.NW, "window grid mismatch"
        p2 = f"{D.DTB}/win/{ns}__c2.npz"
        if os.path.exists(p2):
            z2 = np.load(p2); W2all = z2["W"]; pos2 = pd.Series(np.arange(len(z2["row_id"])), index=z2["row_id"])
        else:
            W2all = np.zeros((0, 3, S.NW), np.float32); pos2 = pd.Series(dtype=int)
        R = RS[RS.netsta == ns]
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
        np.savez(f"{D.DTB}/store{a.tag}/{ns}.npz", keys=np.array(keys), tau=S.TAU,
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
    pd.DataFrame(res_out).to_csv(f"{D.DTB}/series{a.tag}.csv", index=False)
    if res_E:
        pd.concat(res_E, ignore_index=True).to_parquet(f"{D.DTB}/episode_raw{a.tag}.parquet", index=False)
    log(f"measure done {time.time() - t0:.0f}s, series {len(res_out)}")
    if not a.stations and not a.fams:
        open(f"{D.DOUT}/MEASURE{a.tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
