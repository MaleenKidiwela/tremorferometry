#!/usr/bin/env python
"""s11_r4_measure -- S11 round 4, stage S1 (PLAN.md §8.0). The round-3 REDPy measurement (s11_r3_measure.run_series: cc
weights, B1 reference path, B2/B2p all-pairs graphs on >= 50-event episode nodes, 50-draw event bootstrap) with ONE change:
a scrambled-control graph node needs >= 10 control events (round 3: 50). Data values and their bootstrap are unchanged (the
data draws precede the control draws; checked against round 3).
  near: same rows/windows as round 3 (SVI main + ext; WA main) -> round4/<region>/tables/r4_episode[_ext].parquet
  far : distant windows (far grid -28..+10 s, grid-relative pre-P as s11_dist_measure), far series with causality > 1.5
        -> r4_episode_far.parquet + cc-weighted far store round4/<region>/tables/store_ccw_far/<NET.STA>.npz
eps = -dv/v.
    S11_REGION=... python s11_r4_measure.py --set near|ext|far [--workers 12] [--stations ..] [--fams ..] [--tag _test]
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

SET = "near"
if "--set" in sys.argv:
    SET = sys.argv[sys.argv.index("--set") + 1]
if SET == "far":
    D.set_grid(D.FAR_GRID)
import s11_r4_common as R4
import s11_r3_inject as J
import s10_measure as M10
from s10_measure import realign, pick_S, wmean
from s11_measure import Ref3, stack3, boot3
from s11_r3_measure import event_cc, ctrl_cc, graph_boot, NBOOT_REF, NBOOT_G

NODE_MIN, CTRL_NODE_MIN = R4.NODE_MIN, R4.CTRL_NODE_MIN
WIN = C.WINS[C.HEAD]


def prewin(tS, SP):
    tr = S.TAU - tS
    if SET == "far":
        return (tr >= S.T0W + 1.5) & (S.TAU >= S.T0W + 0.5) & (tr <= -(SP + 1.0))
    return (tr >= M10.MIR0) & (S.TAU >= -13.5) & (tr <= -(SP + 1.0))


def run_series(fam, ns, ep, Xd, Wd, sp, X2, W2, unit, nb_ref=NBOOT_REF, nb_g=NBOOT_G):
    """= s11_r3_measure.run_series (identical rng order) with CTRL_NODE_MIN for control graph nodes and grid-aware pre-P"""
    rng = np.random.default_rng(C.SEED + 3 + fam * 1000 + zlib.crc32(f'{ns}|{ep}'.encode()) % 997)
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
    SP = sp if np.isfinite(sp) else (14.0 if SET == "far" else 6.0)
    tr = S.TAU - tS
    mir = prewin(tS, SP); cw = (tr >= 2.0) & (tr <= 5.0)
    rms = lambda x, m: float(np.sqrt(np.mean(x[m] ** 2))) if m.any() else np.nan
    out.update(tS_pick=tS, causal25=rms(G1[0], cw) / rms(G1[0], mir), realign_lag_med_ms=float(1e3 * np.median(lag_d)))
    cc_e = event_cc(Wd, P, tS); wcc = np.clip(np.nan_to_num(cc_e), 0, None)
    cc_c = ctrl_cc(W2, G1, tS); wcc2 = np.clip(np.nan_to_num(cc_c), 0, None)
    out.update(cc_event_med=float(np.nanmedian(cc_e)), cc_event_q90=float(np.nanquantile(cc_e, 0.9)),
               frac_w0=float(np.mean(wcc <= 0)), cc_ctrl_med=float(np.nanmedian(cc_c)) if len(cc_c) else np.nan)
    if wcc.sum() <= 0:
        return out, None
    days = np.array(sorted(Xd.date.unique()))
    dix = [np.flatnonzero((Xd.date == d).to_numpy()) for d in days]
    St = np.zeros((len(days), 3, S.NW), np.float32); Ws = np.zeros((len(days), 3)); Nd = np.zeros(len(days), int)
    for k, ix in enumerate(dix):
        s_, w_ = stack3(Wd[ix], wcc[ix]); St[k] = np.nan_to_num(s_); Ws[k] = w_; Nd[k] = len(ix)
    c2_rows, Sc, Cn = [], [], []
    for k, d in enumerate(days):
        mc = np.flatnonzero((X2.date == d).to_numpy()) if len(X2) else np.zeros(0, int)
        if len(mc) > len(dix[k]):
            mc = np.sort(rng.choice(mc, len(dix[k]), replace=False))
        if len(mc) == 0 or wcc2[mc].sum() <= 0:
            continue
        s_, w_ = stack3(W2[mc], wcc2[mc]); Sc.append(np.nan_to_num(s_).astype(np.float32)); Cn.append(w_); c2_rows.append(k)
    store = dict(days=days, St=St, Ws=Ws, Nd=Nd, c2_day=np.array(c2_rows, int), Sc=np.array(Sc, np.float32).reshape(-1, 3, S.NW),
                 Cn=np.array(Cn).reshape(-1, 3), G1=G1.astype(np.float32), tS=tS, sp=SP)
    ep_row = Xd.date.map(unit).to_numpy()
    ep2 = X2.date.map(unit).to_numpy() if len(X2) else np.zeros(0, object)
    Gw = np.stack([np.nansum(Ws[:, c, None] * St[:, c], 0) for c in range(3)]); Wt = Ws.sum(0)
    rows = []; node_d, node_c = [], []
    lagms = 1e3 * lag_d
    for e_ in pd.unique(ep_row):
        if not isinstance(e_, str):
            continue
        ixe = np.flatnonzero(ep_row == e_)
        if len(ixe) < 3:
            continue
        r = dict(fam=fam, netsta=ns, epoch=ep, episode=e_, n=int(len(ixe)), sum_wcc=float(wcc[ixe].sum()),
                 cc_event_med=float(np.nanmedian(cc_e[ixe])), lag_med_ms=float(np.median(lagms[ixe])),
                 lag_wmed_ms=float(np.average(lagms[ixe], weights=np.maximum(wcc[ixe], 1e-9))) if wcc[ixe].sum() > 0 else np.nan)
        Sd, ws = stack3(Wd[ixe], wcc[ixe])
        if np.any(ws > 0):
            ref = np.stack([(Gw[c] - ws[c] * np.nan_to_num(Sd[c])) / (Wt[c] - ws[c]) if Wt[c] - ws[c] > 0 else np.full(S.NW, np.nan)
                            for c in range(3)])
            ref[~np.isfinite(G1).all(1)] = np.nan
            if np.isfinite(ref).all(1).any() and np.isfinite(Sd).all(1).any():
                Rf = Ref3(ref, tS, WIN)
                e1, c1 = Rf.measure(Sd[None], (0, 1, 2)); r["b1_eps"], r["b1_cc"] = float(e1[0]), float(c1[0])
                eb = Rf.measure(boot3(Wd[ixe], wcc[ixe], nb_ref, rng), (0, 1, 2))[0]; eb = eb[np.isfinite(eb) & (np.abs(eb) < C.RAIL)]
                r["b1_sig"] = float(np.std(eb, ddof=1)) if len(eb) > 10 else np.nan
                mc = np.flatnonzero(ep2 == e_) if len(X2) else np.zeros(0, int)
                if len(mc) > len(ixe):
                    mc = np.sort(rng.choice(mc, len(ixe), replace=False))
                r["n_ctrl"] = int(len(mc))
                if len(mc) >= 3 and wcc2[mc].sum() > 0:
                    Cs = stack3(W2[mc], wcc2[mc])[0]
                    e2, c2 = Rf.measure(Cs[None], (0, 1, 2)); r["b1_c2_eps"], r["b1_c2_cc"] = float(e2[0]), float(c2[0])
                    eb = Rf.measure(boot3(W2[mc], wcc2[mc], nb_ref, rng), (0, 1, 2))[0]; eb = eb[np.isfinite(eb) & (np.abs(eb) < C.RAIL)]
                    r["b1_c2_sig"] = float(np.std(eb, ddof=1)) if len(eb) > 10 else np.nan
                if len(ixe) >= NODE_MIN:
                    node_d.append((len(rows), ixe))
                    if len(mc) >= CTRL_NODE_MIN:
                        node_c.append((len(rows), mc))
        rows.append(r)
    E = pd.DataFrame(rows)
    for col in ("b1_eps", "b1_c2_eps"):
        if col in E:
            rl = E[col].abs() >= C.RAIL; E.loc[rl, col] = np.nan
            sc = col.replace("eps", "sig")
            if sc in E:
                E.loc[rl, sc] = np.nan
    for wname, wd, wc in (("g", wcc, wcc2), ("gp", P, None)):
        if len(node_d) >= 3:
            v, sb, sv, chi2, mpc = graph_boot(Wd, wd, [ix for _, ix in node_d], tS, rng, nb_g)
            for (ri, _), a, b, c, d in zip(node_d, v, sb, sv, mpc):
                E.loc[ri, f"{wname}_eps"] = a; E.loc[ri, f"{wname}_sig"] = b; E.loc[ri, f"{wname}_sigmis"] = c; E.loc[ri, f"{wname}_paircc"] = d
            E[f"{wname}_chi2"] = chi2
        if len(node_c) >= 3:
            wcn = wc if wc is not None else np.ones(len(W2))
            v, sb, sv, chi2, mpc = graph_boot(W2, wcn, [ix for _, ix in node_c], tS, rng, nb_g)
            for (ri, _), a, b in zip(node_c, v, sb):
                E.loc[ri, f"{wname}_c2_eps"] = a; E.loc[ri, f"{wname}_c2_sig"] = b
    return out, (E, store)


def _run(args):
    try:
        return run_series(*args)
    except Exception as e:
        import traceback
        return dict(fam=args[0], netsta=args[1], epoch=args[2], n_data=-1, error=f"{e} {traceback.format_exc()[-800:]}"), None


def save_store(st_all, path):
    keys = sorted(st_all)
    if not keys:
        return
    cat = lambda f: np.concatenate([st_all[k][f] for k in keys])
    np.savez(path, keys=np.array(keys), tau=S.TAU, tS=np.array([st_all[k]["tS"] for k in keys]), sp=np.array([st_all[k]["sp"] for k in keys]),
             G1=np.array([st_all[k]["G1"] for k in keys]).reshape(-1, 3, S.NW),
             n_days=np.array([len(st_all[k]["days"]) for k in keys]), n_c2=np.array([len(st_all[k]["c2_day"]) for k in keys]),
             days=np.concatenate([st_all[k]["days"].astype(str) for k in keys]), St=cat("St"), Ws=cat("Ws"), Nd=cat("Nd"),
             c2_day=cat("c2_day"), Sc=cat("Sc"), Cn=cat("Cn"))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--stations", default="")
    ap.add_argument("--fams", default=""); ap.add_argument("--tag", default=""); ap.add_argument("--set", default="near")
    a = ap.parse_args()
    sfx = {"near": "", "ext": "_ext", "far": "_far"}[a.set] + a.tag
    lg = R4.log(f"measure{sfx}"); t0 = time.time()
    if a.set == "far":
        RS = pd.read_parquet(f"{R4.DIST}/far_rows.parquet")
        O = pd.read_csv(f"{R4.DIST}/far_offsets.csv").set_index(["fam", "netsta"])
        ser = pd.read_csv(f"{R4.DIST}/series.csv")
        WIN_ = f"{R4.DIST}/win"
        if C.REGION == "SVI":
            import s11_paths as P_; unit = P_.unit_map()
        else:
            A = pd.read_csv(f"{C.MTB}/active_days.csv"); unit = dict(zip(A.date, A.episode))
        os.makedirs(f"{R4.TB4}/store_ccw_far{a.tag}", exist_ok=True)
    else:
        TBI = f"{C.TB}/ext" if a.set == "ext" else C.MTB
        WIN_ = f"{TBI}/win" if a.set == "ext" else C.MWIN
        X = pd.read_parquet(f"{TBI}/rows.parquet"); RSx = pd.read_parquet(f"{TBI}/row_station.parquet")
        RS = RSx.merge(X, on=["row_id", "date"])
        O = pd.read_csv(f"{TBI}/offsets_sp.csv").set_index(["fam", "netsta"])
        A = pd.read_csv(f"{TBI}/active_days.csv"); unit = dict(zip(A.date, A.episode))
        ser = pd.read_csv(f"{C.TB}/series{'_ext' if a.set == 'ext' else ''}.csv")
    causal = set(zip(ser.fam[ser.causal25 > 1.5].astype(int), ser.netsta[ser.causal25 > 1.5], ser.epoch[ser.causal25 > 1.5]))
    stas = a.stations.split(",") if a.stations else sorted(RS.netsta.unique())
    fams = [int(x) for x in a.fams.split(",")] if a.fams else None
    import multiprocessing as mp
    res_out, res_E = [], []
    for ns in stas:
        if not os.path.exists(f"{WIN_}/{ns}.npz"):
            continue
        R_ = RS[RS.netsta == ns]
        if fams:
            R_ = R_[R_.fam.isin(fams)]
        z = np.load(f"{WIN_}/{ns}.npz"); Wall = z["W"]; pos = pd.Series(np.arange(len(z["row_id"])), index=z["row_id"])
        assert Wall.shape[-1] == S.NW, f"grid mismatch {ns}: {Wall.shape[-1]} vs {S.NW}"
        p2 = f"{WIN_}/{ns}__c2.npz"
        if os.path.exists(p2):
            z2 = np.load(p2); W2all = z2["W"]; pos2 = pd.Series(np.arange(len(z2["row_id"])), index=z2["row_id"])
        else:
            W2all = np.zeros((0, 3, S.NW), np.float32); pos2 = pd.Series(dtype=int)
        jobs = []
        for (fam, ep), g in R_.groupby(["fam", "epoch"]):
            if (int(fam), ns, ep) not in causal:
                continue
            gd = g[(g.series == "data") & g.row_id.isin(pos.index)].reset_index(drop=True)
            gc = g[(g.series == "c1") & g.row_id.isin(pos2.index)].reset_index(drop=True)
            jobs.append((int(fam), ns, ep, gd, Wall[pos[gd.row_id].to_numpy()], float(O.sp.get((fam, ns), np.nan)),
                         gc, W2all[pos2[gc.row_id].to_numpy()] if len(gc) else np.zeros((0, 3, S.NW), np.float32), unit))
        del Wall, W2all
        if not jobs:
            continue
        jobs.sort(key=lambda j: -len(j[3]))
        st_all = {}
        with mp.get_context("fork").Pool(a.workers, maxtasksperchild=10) as p:
            for out, res in p.imap_unordered(_run, jobs):
                res_out.append(out)
                if "error" in out:
                    lg(f"ERROR {out['fam']} {out['netsta']} {out['epoch']}: {out['error']}")
                if res is not None:
                    E, st = res; res_E.append(E); st_all[f"{out['fam']}|{ns}|{out['epoch']}"] = st
        if a.set == "far":
            save_store(st_all, f"{R4.TB4}/store_ccw_far{a.tag}/{ns}.npz")
        lg(f"{ns}: {len(jobs)} causal series, {len(st_all)} measured, {time.time() - t0:.0f}s")
        del st_all, jobs
    pd.DataFrame(res_out).to_csv(f"{R4.TB4}/r4_series{sfx}.csv", index=False)
    E = pd.concat(res_E, ignore_index=True) if res_E else pd.DataFrame()
    E.to_parquet(f"{R4.TB4}/r4_episode{sfx}.parquet", index=False)
    lg(f"measure done {time.time() - t0:.0f}s, series {len(res_out)}, episodes {len(E)}; data nodes {int(E.g_eps.notna().sum()) if 'g_eps' in E else 0}, "
       f"control nodes {int(E.g_c2_eps.notna().sum()) if 'g_c2_eps' in E else 0}")
    if a.set != "far":                          # code check: data values reproduce round 3
        import s11_r3_common as R
        p3 = f"{R.TB3}/r3_episode{'_ext' if a.set == 'ext' else ''}{a.tag}.parquet"
        if os.path.exists(p3) and len(E):
            E3 = pd.read_parquet(p3)
            m = E.merge(E3, on=["fam", "netsta", "epoch", "episode"], suffixes=("", "_r3"))
            for col in ("b1_eps", "b1_sig", "g_eps", "g_sig"):
                if col in m and f"{col}_r3" in m:
                    dd = np.nanmax(np.abs(m[col] - m[f"{col}_r3"]))
                    lg(f"reproduction vs round 3: {col} max |diff| {dd:.2e} over {int(m[col].notna().sum())} values")
            if "g_c2_eps_r3" in m:
                lg(f"control nodes: round 3 {int(m.g_c2_eps_r3.notna().sum())}, round 4 {int(m.g_c2_eps.notna().sum())}")
    if not a.tag and not a.stations and not a.fams:
        open(f"{R4.R4}/MEASURE_{C.REGION}{sfx}_DONE", "w").write(R4.pt_now() + "\n")


if __name__ == "__main__":
    main()
