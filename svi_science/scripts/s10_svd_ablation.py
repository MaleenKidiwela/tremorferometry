#!/usr/bin/env python
"""s10_svd_ablation -- POST-HOC (not in SVD.md): is the SVD shrinkage set by the pre-noise column choice?
Noise columns = the pre-P mirror window (tau - tS in [-12.5, -(S-P) - 1] s, pure noise) instead of the original
[-4, -2) s (P coda here). Per series: rank, sum w, natural sigma by a 50-draw joint bootstrap, step recovery and
fixed-reference leakage at the 3 pre-registered step dates for a = +0.3 %, 2-4 s. Everything else identical."""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, time, zlib
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_measure as M
import s10_svd_measure as V

A = 0.003; WIN = (2.0, 4.0); NB = 50


def run(args):
    fam, ns, ep, X, W, X2, W2, EPI, sp = args
    rng = np.random.default_rng(V.SEED + 7 + fam * 1000 + zlib.crc32(f'{ns}|{ep}'.encode()) % 997)
    p = V.prep(fam, ns, ep, X, W, X2, W2)
    if p is None:
        return []
    Xd, Z, P, tS, days = p["Xd"], p["Z"], p["P"], p["tS"], p["days"]
    n = len(days); tr = S.TAU - tS
    idx = [np.flatnonzero((Xd.date == d).to_numpy()) for d in days]
    ws = np.array([P[ix].sum() for ix in idx])
    St = np.stack([(P[ix] @ Z[ix]) / P[ix].sum() for ix in idx])
    cl = M.clusters(days, EPI)
    bnd = [i for i in range(1, n) if cl[i] != cl[i - 1]]; chosen = []
    for q in (n / 3, n / 2, 2 * n / 3):
        cand = [b for b in bnd if b not in chosen]
        if cand:
            chosen.append(min(cand, key=lambda b: abs(b - q)))
    out = []
    SPv = sp if np.isfinite(sp) else 6.0
    for label, cols in (("orig_[-4,-2)", (tr >= -4.0) & (tr < -2.0)),
                        ("preP_mirror", (tr >= -12.5) & (S.TAU >= -13.5) & (tr <= -(SPv + 1.0)))):
        if cols.sum() < 20:          # S-P too long for a pre-P window inside the stored span (fam 17 at B926)
            continue
        Rf, info = V.svdw(St, cols, info=True)
        e_nat_raw = V.measure_all(St, ws, tS, WIN)[0]; e_nat = V.measure_all(Rf, ws, tS, WIN)[0]
        good = np.abs(e_nat_raw) < V.RAIL
        J = np.array([V.measure_all(V.svdw(np.stack([M.boot_stacks(Z[ix], P[ix], 1, rng)[0] for ix in idx]), cols), ws, tS, WIN)[0]
                      for _ in range(NB)])
        sigJ = np.nanmedian(np.where(np.abs(e_nat) < V.RAIL, J.std(0, ddof=1), np.nan))
        refs = [V.loo_ref(Rf, ws, i) for i in range(n)]
        w = good.astype(float)
        wm = lambda v, m: np.sum(w[m] * v[m]) / np.sum(w[m]) if np.sum(w[m]) > 0 else np.nan
        recs, leaks = [], []
        for b in chosen:
            post = np.arange(n) >= b; pre = ~post
            Ai = St.copy()
            for i in np.flatnonzero(post):
                Ai[i] = M.stretch_trace(St[i], tS, A)
            Aif = V.svdw(Ai, cols)
            vi = V.measure_all(Aif, ws, tS, WIN)[0]; vf = V.measure_all(Aif, ws, tS, WIN, refs=refs)[0]
            recs.append(((wm(vi, post) - wm(vi, pre)) - (wm(e_nat, post) - wm(e_nat, pre))) / A)
            leaks.append(wm(vf - e_nat, pre) / A)
        out.append(dict(fam=fam, netsta=ns, epoch=ep, noise_cols=label, n_days=n, sigma_noise=info["sigma_noise"],
                        rank_gt_nmf=info["rank_gt_nmf"], sum_w=info["sum_w"], sigJ24_med=sigJ,
                        rec_med=float(np.median(recs)), rec_min=float(np.min(recs)), leak_fixed_maxabs=float(np.max(np.abs(leaks)))))
    return out


def main():
    X = pd.read_parquet(f"{S.TB}/rows.parquet"); RS = pd.read_parquet(f"{S.TB}/row_station.parquet")
    O = pd.read_csv(f"{S.TB}/offsets_sp.csv").set_index(["fam", "netsta"])
    EPI = pd.read_csv(f"{S.ROOT}/svi_science/results/r23/tables/episodes_2003_2026.csv")
    jobs = []
    for ns in V.STAS:
        z = np.load(f"{S.WIN}/{ns}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}
        R = RS[RS.netsta == ns].merge(X, on=["row_id", "date"])
        for (fam, ep), g in R.groupby(["fam", "epoch"]):
            if fam not in S.FAMS:
                continue
            g = g[g.row_id.isin(pos)].reset_index(drop=True)
            jobs.append((int(fam), ns, ep, g, z["W"][[pos[r] for r in g.row_id]], g.iloc[:0], np.zeros((0, 3, S.NW), np.float32),
                         EPI, float(O.sp.get((fam, ns), np.nan))))
    import multiprocessing as mp
    rows = []
    with mp.get_context("fork").Pool(8) as pool:
        for r in pool.imap_unordered(run, jobs):
            rows.extend(r)
    T = pd.DataFrame(rows).sort_values(["fam", "netsta", "epoch", "noise_cols"])
    T.to_csv(f"{V.OUTD}/tables/ablation_noise_cols.csv", index=False)
    print(T.round(4).to_string())
    print(T.groupby("noise_cols")[["sigma_noise", "rank_gt_nmf", "sum_w", "sigJ24_med", "rec_med", "leak_fixed_maxabs"]].median())


if __name__ == "__main__":
    main()
