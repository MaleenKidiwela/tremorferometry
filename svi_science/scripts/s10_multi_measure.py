#!/usr/bin/env python
"""s10_multi_measure -- S10 multi step 3 (multi/PLAN.md §0, §2): per (family, station, epoch) series, the pilot path
unchanged (s10_measure: grand stack, S re-pick, cc realign, P-weighted raw day stacks, S-anchored stretch vs the
leave-day-out epoch reference in 2-4 and 3-5 s, 200-draw bootstrap sigma, mirror, P >= 0.75 stacks, scrambled control
c2 through the same path, Revision 2 rail rule) plus per-episode stacks against the leave-episode-out reference.
Not run (characterised in the pilot): injection, pairwise graph, precision vs n, horizontals. eps = -dv/v.

    python s10_multi_measure.py --workers 12
"""
from __future__ import annotations
import os, sys, argparse, time, zlib
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_multi_common as M
from s10_measure import Ref, realign, pick_S, wmean, boot_stacks, WINS, RAIL, MIR0, NBOOT, RNG_SEED


def run_series(fam, ns, ep, Xd, Wd, sp, X2, W2, unit):
    rng = np.random.default_rng(RNG_SEED + fam * 1000 + zlib.crc32(f'{ns}|{ep}'.encode()) % 997)
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
        W2, lag_2, cc_2 = realign(W2, G0, tS0)
    G1 = np.stack([wmean(Wd[:, c], P)[0] for c in range(3)])
    tS = pick_S(G1, tS0)
    SP = sp if np.isfinite(sp) else 6.0
    mend = -(SP + 1.0)
    tr = S.TAU - tS
    coda = {k: np.sqrt(np.mean(G1[0][(tr >= a) & (tr <= b)] ** 2)) for k, (a, b) in WINS.items()}
    mirr = np.sqrt(np.mean(G1[0][(tr >= MIR0) & (S.TAU >= -13.5) & (tr <= mend)] ** 2))
    out.update(tS_pick0=tS0, tS_pick=tS, sp=SP, causal24=coda["24"] / mirr, causal35=coda["35"] / mirr,
               lag_med_abs=float(np.median(np.abs(lag_d))), realign_cc_med=float(np.median(cc_d)),
               realign_cc_ctrl2_med=float(np.median(cc_2)) if len(X2) else np.nan, n_days=int(Xd.date.nunique()))
    days = np.array(sorted(Xd.date.unique()))
    Gw = (P[:, None] * Wd[:, 0]).sum(0); Ptot = P.sum()
    hi = Xd.hiP.to_numpy()
    rows = []
    for d in days:
        m = (Xd.date == d).to_numpy()
        w = P[m]; Zd = Wd[m, 0].astype(np.float64)
        Sd = (w @ Zd) / w.sum()
        ref = (Gw - w @ Zd) / (Ptot - w.sum())
        R = {k: Ref(ref, tS, win) for k, win in WINS.items()}
        RM = {k: Ref(ref, tS, win, mirror_end=mend) for k, win in WINS.items()}
        r = dict(fam=fam, netsta=ns, epoch=ep, date=d, n=int(m.sum()), sumP=float(w.sum()), n_hi=int((m & hi).sum()))
        Bs = boot_stacks(Zd, w, NBOOT, rng) if m.sum() >= 2 else None
        for k in WINS:
            e, c = R[k].measure(Sd)
            r[f"eps{k}"], r[f"cc{k}"] = float(e[0]), float(c[0])
            r[f"sig{k}"] = float(np.std(R[k].measure(Bs)[0], ddof=1)) if Bs is not None else np.nan
            em, cm = RM[k].measure(Sd)
            r[f"mir{k}"], r[f"mircc{k}"] = float(em[0]), float(cm[0])
        mh = m & hi
        if mh.sum() >= 3:
            Zh = Wd[mh, 0].astype(np.float64); Sh = Zh.mean(0); Bh = boot_stacks(Zh, np.ones(mh.sum()), NBOOT, rng)
            for k in WINS:
                e, c = R[k].measure(Sh); r[f"hi_eps{k}"] = float(e[0]); r[f"hi_sig{k}"] = float(np.std(R[k].measure(Bh)[0], ddof=1))
        if len(X2):
            mc = np.flatnonzero((X2.date == d).to_numpy())
            if len(mc) > m.sum():
                mc = np.sort(rng.choice(mc, m.sum(), replace=False))
            r["n_ctrl2"] = int(len(mc))
            if len(mc) >= 3:
                Zc = W2[mc, 0].astype(np.float64); Sc = Zc.mean(0); Bc = boot_stacks(Zc, np.ones(len(mc)), NBOOT, rng)
                for k in WINS:
                    e, c = R[k].measure(Sc); r[f"c2_eps{k}"], r[f"c2_cc{k}"] = float(e[0]), float(c[0])
                    r[f"c2_sig{k}"] = float(np.std(R[k].measure(Bc)[0], ddof=1))
        else:
            r["n_ctrl2"] = 0
        rows.append(r)
    D = pd.DataFrame(rows)
    D["episode"] = D.date.map(unit)
    for pre in ("", "hi_", "c2_"):
        for k in WINS:
            e = f"{pre}eps{k}"
            if e in D:
                D[f"{e}_raw"] = D[e]; rl = D[e].abs() >= RAIL
                D[f"{pre}rail{k}"] = rl.where(D[e].notna())
                D.loc[rl, e] = np.nan
                sg = f"{pre}sig{k}"
                if sg in D:
                    D.loc[rl, sg] = np.nan
    for k in WINS:
        D[f"mirrail{k}"] = D[f"mir{k}"].abs() >= RAIL
    # per-episode stacks (secondary): leave-episode-out reference
    erows = []
    ep_row = Xd.date.map(unit).to_numpy()
    for e_ in pd.unique(D.episode):
        m = ep_row == e_
        if m.sum() < 3 or Ptot - P[m].sum() <= 0:
            continue
        w = P[m]; Zd = Wd[m, 0].astype(np.float64)
        Sd = (w @ Zd) / w.sum(); ref = (Gw - w @ Zd) / (Ptot - w.sum())
        r = dict(fam=fam, netsta=ns, epoch=ep, episode=e_, n=int(m.sum()), sumP=float(w.sum()))
        Bs = boot_stacks(Zd, w, NBOOT, rng)
        # Revision 1 (post-hoc): scrambled-control episode stack, same days, equal weights, <= n data windows
        mc = np.flatnonzero(X2.date.map(unit).to_numpy() == e_) if len(X2) else np.zeros(0, int)
        if len(mc) > m.sum():
            mc = np.sort(rng.choice(mc, m.sum(), replace=False))
        r["n_ctrl2"] = int(len(mc))
        Sc = W2[mc, 0].astype(np.float64).mean(0) if len(mc) >= 3 else None
        Bc = boot_stacks(W2[mc, 0].astype(np.float64), np.ones(len(mc)), NBOOT, rng) if Sc is not None else None
        for k, win in WINS.items():
            R = Ref(ref, tS, win); e, c = R.measure(Sd)
            r[f"eps{k}"], r[f"cc{k}"], r[f"sig{k}"] = float(e[0]), float(c[0]), float(np.std(R.measure(Bs)[0], ddof=1))
            if Sc is not None:
                e, c = R.measure(Sc)
                r[f"c2_eps{k}"], r[f"c2_cc{k}"], r[f"c2_sig{k}"] = float(e[0]), float(c[0]), float(np.std(R.measure(Bc)[0], ddof=1))
        erows.append(r)
    E = pd.DataFrame(erows)
    for pre in ("", "c2_"):
        for k in WINS:
            if len(E) and f"{pre}eps{k}" in E:
                E[f"{pre}eps{k}_raw"] = E[f"{pre}eps{k}"]; rl = E[f"{pre}eps{k}"].abs() >= RAIL; E[f"{pre}rail{k}"] = rl
                E.loc[rl, [f"{pre}eps{k}", f"{pre}sig{k}"]] = np.nan
    return out, (D, E, G1, tS)


def _run(args):
    try:
        return run_series(*args)
    except Exception as e:
        import traceback
        return dict(fam=args[0], netsta=args[1], epoch=args[2], n_data=-1, error=f"{e} {traceback.format_exc()[-600:]}"), None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--stations", default="")
    ap.add_argument("--fams", default=""); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    M.ensure(); log = S.Log(f"{M.LOGD}/measure{a.tag}.log")
    X = pd.read_parquet(f"{M.TB}/rows.parquet")
    RS = pd.read_parquet(f"{M.TB}/row_station.parquet")
    O = pd.read_csv(f"{M.TB}/offsets_sp.csv").set_index(["fam", "netsta"])
    A = pd.read_csv(f"{M.TB}/active_days.csv")
    unit = dict(zip(A.date, A.episode))
    stas = a.stations.split(",") if a.stations else M.STATIONS
    fams = [int(x) for x in a.fams.split(",")] if a.fams else None

    def series_iter():
        for ns in stas:
            p = f"{M.WIN}/{ns}.npz"
            if not os.path.exists(p):
                log(f"{ns}: no windows"); continue
            z = np.load(p); Wall = z["W"]; pos = pd.Series(np.arange(len(z["row_id"])), index=z["row_id"])
            p2 = f"{M.WIN}/{ns}__c2.npz"
            if os.path.exists(p2):
                z2 = np.load(p2); W2all = z2["W"]; pos2 = pd.Series(np.arange(len(z2["row_id"])), index=z2["row_id"])
            else:
                W2all = np.zeros((0, 3, S.NW), np.float32); pos2 = pd.Series(dtype=int)
            R = RS[RS.netsta == ns].merge(X, on=["row_id", "date"])
            if fams:
                R = R[R.fam.isin(fams)]
            for (fam, ep), g in R.groupby(["fam", "epoch"]):
                gd = g[(g.series == "data") & g.row_id.isin(pos.index)].reset_index(drop=True)
                gc = g[(g.series == "c1") & g.row_id.isin(pos2.index)].reset_index(drop=True)
                yield (int(fam), ns, ep, gd, Wall[pos[gd.row_id].to_numpy()], float(O.sp.get((fam, ns), np.nan)),
                       gc, W2all[pos2[gc.row_id].to_numpy()] if len(gc) else np.zeros((0, 3, S.NW), np.float32), unit)
            del Wall, W2all

    import multiprocessing as mp
    res_out, res_D, res_E, res_G = [], [], [], {}
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers, maxtasksperchild=20) as p:
        for out, res in p.imap_unordered(_run, series_iter()):
            res_out.append(out)
            if res is not None:
                D, E, G1, tS = res
                res_D.append(D); res_E.append(E); res_G[f"{out['fam']}|{out['netsta']}|{out['epoch']}"] = (G1, tS)
            if "error" in out:
                log(f"ERROR {out['fam']} {out['netsta']} {out['epoch']}: {out['error']}")
            log(f"  {out['fam']} {out['netsta']} {out['epoch']}: n_data {out['n_data']} days {out.get('n_days', 0)} "
                f"causal24 {out.get('causal24', np.nan):.2f} {time.time() - t0:.0f}s")
    sfx = a.tag
    pd.DataFrame(res_out).to_csv(f"{M.TB}/series{sfx}.csv", index=False)
    pd.concat(res_D, ignore_index=True).to_parquet(f"{M.TB}/daily{sfx}.parquet", index=False)
    pd.concat(res_E, ignore_index=True).to_parquet(f"{M.TB}/episode_stacks{sfx}.parquet", index=False)
    np.savez(f"{M.TB}/grandstacks{sfx}.npz", keys=np.array(list(res_G)), G1=np.array([v[0] for v in res_G.values()]),
             tS=np.array([v[1] for v in res_G.values()]), tau=S.TAU)
    log(f"measure done {time.time() - t0:.0f}s, series {len(res_out)}, with values {len(res_D)}")
    open(f"{M.OUT}/MEASURE{sfx}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
