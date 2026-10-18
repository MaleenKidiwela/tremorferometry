#!/usr/bin/env python
"""s11_r3_measure -- S11 round 3, task 2 (PLAN §6.2): REDPy findings in the S11 raw path.
Per (family, station, epoch) series with causality > 1.5 in the existing series table: same rows, windows, 3C realignment
and S re-pick as s11_measure (imported), then
  * w_cc = max(cc_i, 0): zero-lag 3C cc of each realigned event with its leave-self-out P-weighted grand stack, S-1..S+5 s;
    scrambled-control windows: cc with the grand stack (same rule);
  * cc-weighted day-stack store (store_ccw[_ext]/<NET.STA>.npz, + events per day Nd) for the waveform-level injections;
  * B1 ccw-ref: episode stack (w_cc) vs leave-episode-out reference, joint 3C 2-5 s, 200-draw event bootstrap; control too;
  * B2 ccw-graph and B2p P-graph: nodes = episode stacks with >= 50 data events; all node pairs (joint 3C 2-5 s), REDPy solve;
    sigma = 50-draw event bootstrap (all pairs re-measured and re-solved per draw); misclosure sigma kept as a column;
    scrambled-control nodes (>= 50 control events, <= n data) through the same graph + bootstrap.
eps = -dv/v.
    S11_REGION=... python s11_r3_measure.py [--workers 12] [--ext] [--fams 10,87] [--stations ..] [--tag _test]
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
import s11_r3_common as R
import s11_r3_inject as J
from s10_measure import realign, pick_S, wmean, MIR0
from s11_measure import Ref3, stack3, boot3

NBOOT_REF, NBOOT_G, NODE_MIN = 200, 50, 50
CCW = (-1.0, 5.0)
WIN = C.WINS[C.HEAD]


def zfeat(W, tS):
    """(n,3,NW) -> per-component demeaned windows on S-1..S+5 s (n,3,m), NaN components kept NaN"""
    m = (S.TAU >= tS + CCW[0]) & (S.TAU <= tS + CCW[1])
    X = W[:, :, m].astype(np.float64)
    return X - np.nanmean(X, 2, keepdims=True)


def cc_with(X, Gw):
    """X (n,3,m) demeaned, Gw (n or 1, 3, m) demeaned refs -> mean over finite components of normalised dot"""
    num = np.nansum(X * Gw, 2)
    den = np.sqrt(np.nansum(X * X, 2) * np.nansum(Gw * Gw, 2))
    ok = np.isfinite(X[:, :, 0]) & np.isfinite(Gw[:, :, 0]) & (den > 0)
    c = np.where(ok, num / np.where(den > 0, den, 1), 0.0)
    k = ok.sum(1)
    return np.where(k > 0, c.sum(1) / np.maximum(k, 1), np.nan)


def event_cc(Wd, P, tS):
    X = zfeat(Wd, tS)
    tot = np.zeros(X.shape[1:]); wt = np.zeros(3)
    for c in range(3):
        ok = np.isfinite(X[:, c, 0])
        tot[c] = (P[ok, None] * X[ok, c]).sum(0); wt[c] = P[ok].sum()
    fin = np.isfinite(X[:, :, :1])
    loo = (tot[None] - np.where(fin, P[:, None, None] * np.nan_to_num(X), 0.0)) / np.maximum(wt[None, :, None] - np.where(fin, P[:, None, None], 0.0), 1e-12)
    loo = np.where(fin, loo, np.nan)
    return cc_with(X, loo)


def ctrl_cc(W2, G1, tS):
    if not len(W2):
        return np.zeros(0)
    X = zfeat(W2, tS); g = zfeat(G1[None], tS)
    return cc_with(X, np.broadcast_to(g, X.shape))


def graph_boot(Wd, w, node_ix, tS, rng, nb):
    """graph on node stacks + nb-draw event bootstrap -> v, sig_boot, sig_mis, chi2, median pair cc per node"""
    St = np.array([stack3(Wd[ix], w[ix])[0] for ix in node_ix])
    v, sv, chi2, mpc = J.graph_from_stacks(St, tS)
    B = [boot3(Wd[ix], w[ix], nb, rng) for ix in node_ix]          # list of (nb,3,NW)
    vb = np.full((nb, len(node_ix)), np.nan)
    for b in range(nb):
        Sb = np.array([B[k][b] for k in range(len(node_ix))])
        vb[b] = J.graph_from_stacks(Sb, tS)[0]
    sb = np.nanstd(vb, 0, ddof=1)
    return v, sb, sv, chi2, mpc


def run_series(fam, ns, ep, Xd, Wd, sp, X2, W2, unit, nb_ref=NBOOT_REF, nb_g=NBOOT_G):
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
    SP = sp if np.isfinite(sp) else 6.0
    tr = S.TAU - tS
    mir = (tr >= MIR0) & (S.TAU >= -13.5) & (tr <= -(SP + 1.0)); cw = (tr >= 2.0) & (tr <= 5.0)
    rms = lambda x, m: float(np.sqrt(np.mean(x[m] ** 2)))
    out.update(tS_pick=tS, causal25=rms(G1[0], cw) / rms(G1[0], mir))
    # ---------------- weights
    cc_e = event_cc(Wd, P, tS); wcc = np.clip(np.nan_to_num(cc_e), 0, None)
    cc_c = ctrl_cc(W2, G1, tS); wcc2 = np.clip(np.nan_to_num(cc_c), 0, None)
    out.update(cc_event_med=float(np.nanmedian(cc_e)), cc_event_q90=float(np.nanquantile(cc_e, 0.9)),
               frac_w0=float(np.mean(wcc <= 0)), cc_ctrl_med=float(np.nanmedian(cc_c)) if len(cc_c) else np.nan)
    if wcc.sum() <= 0:
        return out, None
    # ---------------- ccw day-stack store
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
    # ---------------- episodes
    ep_row = Xd.date.map(unit).to_numpy()
    ep2 = X2.date.map(unit).to_numpy() if len(X2) else np.zeros(0, object)
    Gw = np.stack([np.nansum(Ws[:, c, None] * St[:, c], 0) for c in range(3)]); Wt = Ws.sum(0)
    rows = []; node_d, node_c = [], []
    for e_ in pd.unique(ep_row):
        if not isinstance(e_, str):
            continue
        ixe = np.flatnonzero(ep_row == e_)
        if len(ixe) < 3:
            continue
        r = dict(fam=fam, netsta=ns, epoch=ep, episode=e_, n=int(len(ixe)), sum_wcc=float(wcc[ixe].sum()),
                 cc_event_med=float(np.nanmedian(cc_e[ixe])))
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
                    if len(mc) >= NODE_MIN:
                        node_c.append((len(rows), mc))
        rows.append(r)
    E = pd.DataFrame(rows)
    for col in ("b1_eps", "b1_c2_eps"):
        if col in E:
            rl = E[col].abs() >= C.RAIL; E.loc[rl, col] = np.nan
            sc = col.replace("eps", "sig")
            if sc in E:
                E.loc[rl, sc] = np.nan
    # ---------------- graphs (B2 ccw, B2p P) on >= 50-event episode nodes
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


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--stations", default="")
    ap.add_argument("--fams", default=""); ap.add_argument("--tag", default=""); ap.add_argument("--ext", action="store_true")
    a = ap.parse_args()
    sfx = ("_ext" if a.ext else "") + a.tag
    lg = R.log(f"measure{sfx}")
    TBI = f"{C.TB}/ext" if a.ext else C.MTB
    WIN_ = f"{TBI}/win" if a.ext else C.MWIN
    os.makedirs(f"{R.TB3}/store_ccw{sfx}", exist_ok=True)
    X = pd.read_parquet(f"{TBI}/rows.parquet"); RS = pd.read_parquet(f"{TBI}/row_station.parquet")
    O = pd.read_csv(f"{TBI}/offsets_sp.csv").set_index(["fam", "netsta"])
    A = pd.read_csv(f"{TBI}/active_days.csv"); unit = dict(zip(A.date, A.episode))
    ser = pd.read_csv(f"{C.TB}/series{'_ext' if a.ext else ''}.csv")
    causal = set(zip(ser.fam[ser.causal25 > 1.5].astype(int), ser.netsta[ser.causal25 > 1.5], ser.epoch[ser.causal25 > 1.5]))
    stas = a.stations.split(",") if a.stations else C.STATIONS
    fams = [int(x) for x in a.fams.split(",")] if a.fams else None
    import multiprocessing as mp
    res_out, res_E = [], []; t0 = time.time()
    for ns in stas:
        if not os.path.exists(f"{WIN_}/{ns}.npz"):
            continue
        R_ = RS[RS.netsta == ns].merge(X, on=["row_id", "date"])
        if fams:
            R_ = R_[R_.fam.isin(fams)]
        keys = [(f, e) for (f, e) in R_.groupby(["fam", "epoch"]).groups if (int(f), ns, e) in causal]
        if not keys:
            continue
        z = np.load(f"{WIN_}/{ns}.npz"); Wall = z["W"]; pos = pd.Series(np.arange(len(z["row_id"])), index=z["row_id"])
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
        jobs.sort(key=lambda j: -len(j[3]))
        lg(f"{ns}: {len(jobs)} causal series")
        st_all = {}
        with mp.get_context("fork").Pool(a.workers, maxtasksperchild=10) as p:
            for out, res in p.imap_unordered(_run, jobs):
                res_out.append(out)
                if "error" in out:
                    lg(f"ERROR {out['fam']} {out['netsta']} {out['epoch']}: {out['error']}")
                if res is not None:
                    E, st = res; res_E.append(E); st_all[f"{out['fam']}|{ns}|{out['epoch']}"] = st
        del jobs
        keys = sorted(st_all)
        if keys:
            cat = lambda f, empty: np.concatenate([st_all[k][f] for k in keys]) if keys else empty
            np.savez(f"{R.TB3}/store_ccw{sfx}/{ns}.npz", keys=np.array(keys), tau=S.TAU,
                     tS=np.array([st_all[k]["tS"] for k in keys]), sp=np.array([st_all[k]["sp"] for k in keys]),
                     G1=np.array([st_all[k]["G1"] for k in keys]).reshape(-1, 3, S.NW),
                     n_days=np.array([len(st_all[k]["days"]) for k in keys]), n_c2=np.array([len(st_all[k]["c2_day"]) for k in keys]),
                     days=np.concatenate([st_all[k]["days"].astype(str) for k in keys]), St=cat("St", None), Ws=cat("Ws", None),
                     Nd=cat("Nd", None), c2_day=cat("c2_day", None), Sc=cat("Sc", None), Cn=cat("Cn", None))
        lg(f"{ns}: done ({len(keys)} series) {time.time() - t0:.0f}s")
        del st_all
    pd.DataFrame(res_out).to_csv(f"{R.TB3}/r3_series{sfx}.csv", index=False)
    pd.concat(res_E, ignore_index=True).to_parquet(f"{R.TB3}/r3_episode{sfx}.parquet", index=False)
    lg(f"measure done {time.time() - t0:.0f}s, series {len(res_out)}")
    if not a.tag:
        open(f"{R.R3}/MEASURE{sfx}_DONE", "w").write(R.pt_now() + "\n")


if __name__ == "__main__":
    main()
