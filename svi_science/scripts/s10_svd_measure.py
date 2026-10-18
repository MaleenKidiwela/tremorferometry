#!/usr/bin/env python
"""s10_svd_measure -- SVD-Wiener denoising of the S10 daily stacks (rules: results/s10_dvv/pilot/svd/SVD.md).

Same windows / S re-pick / cc-realign / P-weighted active-day stacks as s10_measure.py (its functions are imported).
SVD-Wiener filter = scripts/dvv_roll30cal.py::_family unchanged (row peak-normalisation, sigma from pre-noise columns,
nmf = sigma (sqrt(ne) + sqrt(ns)), w = s^2 / (s^2 + nmf^2)). Per series: natural raw + SVD values (2-4, 3-5 s) against the
leave-day-out reference, per-day bootstrap with the SVD re-run in every draw, joint bootstrap, half-split check,
step / transient injection, scrambled control (c2) through the same path. eps = -dv/v.

    python s10_svd_measure.py [--workers 8] [--series 10|PB.B011|E1,...]
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
import s10_measure as M

OUTD = f"{S.OUT}/svd"
STAS = ["PB.B011", "PB.B009", "PB.B010", "PB.B926", "CN.PGC"]
WINS = M.WINS
RAIL = M.RAIL
NBOOT, NBOOT_INS, NSPLIT = 200, 100, 20
AMPS = (0.001, 0.003)
A_TR = 0.003
SEED = 20261012


# ------------------------------------------------------------------------------ SVD-Wiener (dvv_roll30cal._family)
def svdw(R, noise_cols, info=False):
    """R (ne, ns) day stacks -> filtered (same scale). Identical algebra to dvv_roll30cal.py::_family."""
    pk = np.max(np.abs(R), axis=1, keepdims=True); pk[pk == 0] = 1.0; Rn = R / pk
    ne, ns = Rn.shape; sig = Rn[:, noise_cols].std()
    try:
        U, sv, Vt = np.linalg.svd(Rn, full_matrices=False)
        nmf = sig * (np.sqrt(ne) + np.sqrt(ns)); w = sv ** 2 / (sv ** 2 + nmf ** 2); Rf = (U * (sv * w)) @ Vt
    except np.linalg.LinAlgError:
        Rf = Rn; sv = np.full(ne, np.nan); w = np.full(ne, np.nan); nmf = np.nan
    if info:
        return Rf * pk, dict(rank_gt_nmf=int(np.sum(sv > nmf)), sum_w=float(np.nansum(w)), nmf=float(nmf),
                             sigma_noise=float(sig), s1=float(sv[0]), s2=float(sv[1]) if len(sv) > 1 else np.nan,
                             w1=float(w[0]), w2=float(w[1]) if len(w) > 1 else np.nan,
                             w3=float(w[2]) if len(w) > 2 else np.nan)
    return Rf * pk


def loo_ref(Rf, ws, i):
    return (ws @ Rf - ws[i] * Rf[i]) / (ws.sum() - ws[i])


def measure_all(Rf, ws, tS, win, rows=None, refs=None):
    """eps, cc of rows of Rf against their leave-day-out references (or given refs)"""
    n = Rf.shape[0]; rows = range(n) if rows is None else rows
    e = np.full(n, np.nan); c = np.full(n, np.nan)
    G = ws @ Rf
    for i in rows:
        ref = refs[i] if refs is not None else (G - ws[i] * Rf[i]) / (ws.sum() - ws[i])
        ee, cc = M.Ref(ref, tS, win).measure(Rf[i]); e[i], c[i] = ee[0], cc[0]
    return e, c


def nan_rail(x):
    x = np.array(x, float); x[np.abs(x) >= RAIL] = np.nan; return x


def theil(x, y):
    from scipy import stats
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 5:
        return np.nan, np.nan
    return float(stats.theilslopes(y[ok], x[ok])[0]), float(np.mean(np.sign(x[ok]) == np.sign(y[ok])))


# ------------------------------------------------------------------------------ one series
def prep(fam, ns, ep, X, W, X2, W2):
    """replicates s10_measure.run_series up to the day stacks"""
    data = (X.series == "data").to_numpy(); okZ = np.isfinite(W[:, 0, 0])
    Xd, Wd = X[data & okZ].reset_index(drop=True), W[data & okZ]
    if len(Xd) < 6 or Xd.date.nunique() < 3:
        return None
    P = Xd.P.to_numpy()
    G0 = np.stack([M.wmean(Wd[:, c], P)[0] for c in range(3)])
    tS0 = M.pick_S(G0, 0.0)
    Wd = M.realign(Wd, G0, tS0)[0]
    ok2 = np.isfinite(W2[:, 0, 0]) if len(X2) else np.zeros(0, bool)
    X2, W2 = X2[ok2].reset_index(drop=True), W2[ok2]
    if len(X2):
        W2 = M.realign(W2, G0, tS0)[0]
    G1 = np.stack([M.wmean(Wd[:, c], P)[0] for c in range(3)])
    tS = M.pick_S(G1, tS0)
    days = np.array(sorted(Xd.date.unique()))
    return dict(Xd=Xd, Z=Wd[:, 0].astype(np.float64), P=P, X2=X2, Z2=W2[:, 0].astype(np.float64) if len(X2) else None,
                tS=tS, days=days)


def run_series(fam, ns, ep, X, W, X2, W2, EPI):
    t0 = time.time()
    rng = np.random.default_rng(SEED + fam * 1000 + zlib.crc32(f'{ns}|{ep}'.encode()) % 997)
    p = prep(fam, ns, ep, X, W, X2, W2)
    out = dict(fam=fam, netsta=ns, epoch=ep)
    if p is None:
        return out, None
    Xd, Z, P, tS, days = p["Xd"], p["Z"], p["P"], p["tS"], p["days"]
    n = len(days)
    out.update(n_days=n, n_rep=len(Xd), tS=tS)
    tr = S.TAU - tS
    noise = (tr >= -4.0) & (tr < -2.0)
    idx = [np.flatnonzero((Xd.date == d).to_numpy()) for d in days]
    ws = np.array([P[ix].sum() for ix in idx])
    St = np.stack([(P[ix] @ Z[ix]) / P[ix].sum() for ix in idx])
    cl = M.clusters(days, EPI)
    Rf, info = svdw(St, noise, info=True)
    out.update(info)
    D = pd.DataFrame(dict(fam=fam, netsta=ns, epoch=ep, date=days, n=[len(ix) for ix in idx], sumP=ws, cluster=cl))
    # ---- natural values
    for k, win in WINS.items():
        for tag, A in (("raw", St), ("svd", Rf)):
            e, c = measure_all(A, ws, tS, win)
            D[f"{tag}_eps{k}_raw"] = e; D[f"{tag}_eps{k}"] = nan_rail(e); D[f"{tag}_cc{k}"] = c
    # ---- per-day bootstrap (raw: fixed ref, vectorised; svd: SVD re-run in every draw)
    for k in WINS:
        D[f"raw_sig{k}"] = np.nan; D[f"svd_sig{k}"] = np.nan
    Graw = ws @ St
    for i, ix in enumerate(idx):
        if len(ix) < 2:
            continue
        Bs = M.boot_stacks(Z[ix], P[ix], NBOOT, rng)
        refr = (Graw - ws[i] * St[i]) / (ws.sum() - ws[i])
        for k, win in WINS.items():
            D.loc[i, f"raw_sig{k}"] = float(np.std(M.Ref(refr, tS, win).measure(Bs)[0], ddof=1))
        eb = {k: np.empty(NBOOT) for k in WINS}
        A = St.copy()
        for b in range(NBOOT):
            A[i] = Bs[b]
            Af = svdw(A, noise); ref = loo_ref(Af, ws, i)
            for k, win in WINS.items():
                eb[k][b] = M.Ref(ref, tS, win).measure(Af[i])[0][0]
        for k in WINS:
            D.loc[i, f"svd_sig{k}"] = float(np.std(eb[k], ddof=1))
    for tag in ("raw", "svd"):
        for k in WINS:
            D.loc[D[f"{tag}_eps{k}"].isna(), f"{tag}_sig{k}"] = np.nan   # railed -> not a measurement (Rev. 2)
    out["t_boot"] = time.time() - t0
    # ---- joint bootstrap (all days resampled at once; one SVD per draw)
    J = {f"{tag}{k}": np.full((NBOOT, n), np.nan) for tag in ("raw", "svd") for k in WINS}
    for b in range(NBOOT):
        A = np.stack([M.boot_stacks(Z[ix], P[ix], 1, rng)[0] for ix in idx])
        Af = svdw(A, noise)
        for k, win in WINS.items():
            J[f"raw{k}"][b] = measure_all(A, ws, tS, win)[0]
            J[f"svd{k}"][b] = measure_all(Af, ws, tS, win)[0]
    for key, v in J.items():
        tag, k = key[:3], key[3:]
        D[f"{tag}_sigJ{k}"] = np.where(D[f"{tag}_eps{k}"].notna(), np.std(v, axis=0, ddof=1), np.nan)
    # ---- half-split check
    hs = {f"{tag}{k}": [] for tag in ("raw", "svd") for k in WINS}
    for s in range(NSPLIT):
        keep = [i for i, ix in enumerate(idx) if len(ix) >= 2]
        if len(keep) < 4:
            break
        SA, SB, wA, wB = [], [], [], []
        for i in keep:
            pp = rng.permutation(idx[i]); h = len(pp) // 2; a_, b_ = pp[:h], pp[h:2 * h]
            SA.append((P[a_] @ Z[a_]) / P[a_].sum()); SB.append((P[b_] @ Z[b_]) / P[b_].sum())
            wA.append(P[a_].sum()); wB.append(P[b_].sum())
        SA, SB, wA, wB = np.array(SA), np.array(SB), np.array(wA), np.array(wB)
        FA, FB = svdw(SA, noise), svdw(SB, noise)
        for k, win in WINS.items():
            for tag, a_, b_ in (("raw", SA, SB), ("svd", FA, FB)):
                ea = nan_rail(measure_all(a_, wA, tS, win)[0]); eb_ = nan_rail(measure_all(b_, wB, tS, win)[0])
                d = ea - eb_; hs[f"{tag}{k}"].extend(d[np.isfinite(d)].tolist())
    for key, v in hs.items():
        tag, k = key[:3], key[3:]
        v = np.array(v)
        emp = 1.4826 * np.median(np.abs(v - np.median(v))) / np.sqrt(2) if len(v) > 10 else np.nan
        pred = np.sqrt(2) * np.nanmedian(D[f"{tag}_sig{k}"])
        out[f"{tag}_half_emp{k}"] = emp; out[f"{tag}_half_pred{k}"] = pred; out[f"{tag}_half_ratio{k}"] = emp / pred
        out[f"{tag}_half_n{k}"] = len(v)
    # ---- summary columns raw vs svd
    for tag in ("raw", "svd"):
        for k in WINS:
            out[f"{tag}_sig{k}_med"] = float(np.nanmedian(D[f"{tag}_sig{k}"]))
            out[f"{tag}_sigJ{k}_med"] = float(np.nanmedian(D[f"{tag}_sigJ{k}"]))
            out[f"{tag}_cc{k}_med"] = float(np.nanmedian(D[f"{tag}_cc{k}"]))
            out[f"{tag}_rail{k}"] = float(np.mean(D[f"{tag}_eps{k}"].isna()))
            out[f"{tag}_sd_eps{k}"] = float(np.nanstd(D[f"{tag}_eps{k}"], ddof=1))
        sl, sg = theil(D[f"{tag}_eps24"].to_numpy(), D[f"{tag}_eps35"].to_numpy())
        out[f"{tag}_lapse_slope"] = sl; out[f"{tag}_lapse_sign"] = sg
        out[f"{tag}_lapse_pass"] = None if not np.isfinite(sl) else bool(sl >= 0.6 and sg >= 0.7)
    # ---- injection (step at 3 cluster boundaries; transients)
    inj = []
    if n >= 6:
        bnd = [i for i in range(1, n) if cl[i] != cl[i - 1]]
        chosen = []
        for q in (n / 3, n / 2, 2 * n / 3):
            cand = [b for b in bnd if b not in chosen]
            if cand:
                chosen.append(min(cand, key=lambda b: abs(b - q)))
        sig_w = D["raw_sig24"].fillna(D["raw_sig24"].median()).to_numpy()
        nat = {}
        for k, win in WINS.items():
            nat[("raw", k)] = measure_all(St, ws, tS, win)[0]
            nat[("svd", k)] = measure_all(Rf, ws, tS, win)[0]
        refs_nat = {"raw": [loo_ref(St, ws, i) for i in range(n)], "svd": [loo_ref(Rf, ws, i) for i in range(n)]}
        for k in WINS:
            good = D[f"raw_eps{k}"].notna().to_numpy()
            wts = np.where(good, 1 / np.maximum(D[f"raw_sig{k}"].fillna(np.nanmedian(D[f"raw_sig{k}"])).to_numpy(), 1e-6) ** 2, 0.0)
            nat[("w", k)] = wts
        wm = lambda v, m, w: np.sum(w[m] * v[m]) / np.sum(w[m]) if np.sum(w[m]) > 0 else np.nan

        def inject(mask, a):
            A = St.copy()
            for i in np.flatnonzero(mask):
                A[i] = M.stretch_trace(St[i], tS, a)
            return A
        for b in chosen:
            post = np.arange(n) >= b; pre = ~post
            last_pre = pre & (cl == cl[b - 1]); early_pre = pre & ~last_pre
            for a in AMPS:
                Ai = inject(post, a); Aif = svdw(Ai, noise)
                for k, win in WINS.items():
                    w = nat[("w", k)]
                    for tag, A_ in (("raw", Ai), ("svd", Aif)):
                        vi = measure_all(A_, ws, tS, win)[0]; vn = nat[(tag, k)]
                        vfix = measure_all(A_, ws, tS, win, refs=refs_nat[tag])[0]
                        rec = ((wm(vi, post, w) - wm(vi, pre, w)) - (wm(vn, post, w) - wm(vn, pre, w))) / a
                        lk_fix = wm(vfix - vn, pre, w) / a
                        rec_fix = wm(vfix - vn, post, w) / a
                        lk_con = (((wm(vi, last_pre, w) - wm(vi, early_pre, w)) - (wm(vn, last_pre, w) - wm(vn, early_pre, w))) / a
                                  if early_pre.any() else np.nan)
                        inj.append(dict(fam=fam, netsta=ns, epoch=ep, kind="step", path=tag, win=k, a=a, b_idx=b,
                                        boundary=days[b], n_pre=int(pre.sum()), n_post=int(post.sum()), recovery=rec,
                                        leak_fixed=lk_fix, recovery_fixed=rec_fix, leak_contrast=lk_con))
        # transients (columns only)
        mid = n // 2
        for kind, mask in (("episode", cl == cl[mid]), ("day", np.arange(n) == mid)):
            if mask.all():
                continue
            Ai = inject(mask, A_TR); Aif = svdw(Ai, noise)
            for k, win in WINS.items():
                w = nat[("w", k)].copy(); w[mask] = np.maximum(w[mask], 1e-12)
                for tag, A_ in (("raw", Ai), ("svd", Aif)):
                    vi = measure_all(A_, ws, tS, win)[0]; vn = nat[(tag, k)]
                    rec = ((wm(vi, mask, w) - wm(vi, ~mask, w)) - (wm(vn, mask, w) - wm(vn, ~mask, w))) / A_TR
                    inj.append(dict(fam=fam, netsta=ns, epoch=ep, kind=kind, path=tag, win=k, a=A_TR, b_idx=mid,
                                    boundary=days[mid], n_pre=int((~mask).sum()), n_post=int(mask.sum()), recovery=rec))
    INJ = pd.DataFrame(inj)
    # ---- scrambled control c2 through the same path
    ctrl = []
    if p["Z2"] is not None:
        X2, Z2 = p["X2"], p["Z2"]
        cdays, cst, cix = [], [], []
        for i, d in enumerate(days):
            mc = np.flatnonzero((X2.date == d).to_numpy())
            if len(mc) > len(idx[i]):
                mc = np.sort(rng.choice(mc, len(idx[i]), replace=False))
            if len(mc) >= 3:
                cdays.append(i); cix.append(mc); cst.append(Z2[mc].mean(0))
        if cdays:
            C = np.array(cst)
            refs_raw = [loo_ref(St, ws, i) for i in cdays]; refs_svd = [loo_ref(Rf, ws, i) for i in cdays]
            Cf = svdw(C, noise) if len(cdays) >= 3 else None
            for j, i in enumerate(cdays):
                r = dict(fam=fam, netsta=ns, epoch=ep, date=days[i], n_ctrl=len(cix[j]))
                Bc = M.boot_stacks(Z2[cix[j]], np.ones(len(cix[j])), NBOOT, rng)
                # raw (pilot path)
                for k, win in WINS.items():
                    R_ = M.Ref(refs_raw[j], tS, win)
                    e, c = R_.measure(C[j]); r[f"raw_eps{k}"], r[f"raw_cc{k}"] = e[0], c[0]
                    r[f"raw_sig{k}"] = float(np.std(R_.measure(Bc)[0], ddof=1))
                # own matrix
                if Cf is not None:
                    eb = {k: np.empty(NBOOT) for k in WINS}; A = C.copy()
                    for b in range(NBOOT):
                        A[j] = Bc[b]; Af = svdw(A, noise)
                        for k, win in WINS.items():
                            eb[k][b] = M.Ref(refs_svd[j], tS, win).measure(Af[j])[0][0]
                    for k, win in WINS.items():
                        e, c = M.Ref(refs_svd[j], tS, win).measure(Cf[j])
                        r[f"own_eps{k}"], r[f"own_cc{k}"], r[f"own_sig{k}"] = e[0], c[0], float(np.std(eb[k], ddof=1))
                # inserted into the data matrix
                A = St.copy(); A[i] = C[j]; Af = svdw(A, noise); ref = loo_ref(Af, ws, i)
                for k, win in WINS.items():
                    e, c = M.Ref(ref, tS, win).measure(Af[i]); r[f"ins_eps{k}"], r[f"ins_cc{k}"] = e[0], c[0]
                eb = {k: np.empty(NBOOT_INS) for k in WINS}
                for b in range(NBOOT_INS):
                    A[i] = Bc[b]; Af = svdw(A, noise); ref = loo_ref(Af, ws, i)
                    for k, win in WINS.items():
                        eb[k][b] = M.Ref(ref, tS, win).measure(Af[i])[0][0]
                for k in WINS:
                    r[f"ins_sig{k}"] = float(np.std(eb[k], ddof=1))
                ctrl.append(r)
    CT = pd.DataFrame(ctrl)
    for tag in ("raw", "own", "ins"):
        for k in WINS:
            c_ = f"{tag}_eps{k}"
            if c_ in CT:
                CT[f"{c_}_raw"] = CT[c_]; rl = CT[c_].abs() >= RAIL
                CT.loc[rl, [c_, f"{tag}_sig{k}"]] = np.nan
    out["t_total"] = time.time() - t0
    return out, (D, INJ, CT, dict(St=St.astype(np.float32), Rf=Rf.astype(np.float32), days=days, tS=tS))


def _run(args):
    try:
        return run_series(*args)
    except Exception as e:
        import traceback
        return dict(fam=args[0], netsta=args[1], epoch=args[2], error=f"{e} {traceback.format_exc()[-800:]}"), None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--series", default="", help="fam|NET.STA|epoch,... (test subset)"); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    os.makedirs(f"{OUTD}/tables", exist_ok=True); os.makedirs(f"{OUTD}/logs", exist_ok=True)
    log = S.Log(f"{OUTD}/logs/svd_measure{a.tag}.log")
    X = pd.read_parquet(f"{S.TB}/rows.parquet"); RS = pd.read_parquet(f"{S.TB}/row_station.parquet")
    EPI = pd.read_csv(f"{S.ROOT}/svi_science/results/r23/tables/episodes_2003_2026.csv")
    want = set(a.series.split(",")) if a.series else None
    jobs = []
    for ns in STAS:
        z = np.load(f"{S.WIN}/{ns}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}
        z2 = np.load(f"{S.WIN}/{ns}__c1f.npz"); pos2 = {r: i for i, r in enumerate(z2["row_id"])}
        R = RS[RS.netsta == ns].merge(X, on=["row_id", "date"])
        for (fam, ep), g in R.groupby(["fam", "epoch"]):
            if fam not in S.FAMS or (want and f"{fam}|{ns}|{ep}" not in want):
                continue
            g = g[g.row_id.isin(pos)].reset_index(drop=True)
            W = z["W"][[pos[r] for r in g.row_id]]
            g2 = g[(g.series == "c1") & g.row_id.isin(pos2)].reset_index(drop=True)
            W2 = z2["W"][[pos2[r] for r in g2.row_id]]
            jobs.append((int(fam), ns, ep, g, W, g2, W2, EPI))
    jobs.sort(key=lambda j: -len(j[3]))
    log(f"{len(jobs)} series, {a.workers} workers")
    import multiprocessing as mp
    res_o, res_D, res_I, res_C, stk = [], [], [], [], {}
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers) as pool:
        for o, r in pool.imap_unordered(_run, jobs):
            res_o.append(o)
            if r is not None:
                D, I, C, st = r; res_D.append(D); res_I.append(I); res_C.append(C)
                stk[f"{o['fam']}|{o['netsta']}|{o['epoch']}"] = st
            log(f"  {o['fam']} {o['netsta']} {o['epoch']}: days {o.get('n_days')} rank {o.get('rank_gt_nmf')} "
                f"sumw {o.get('sum_w', np.nan):.2f} sig24 raw {o.get('raw_sig24_med', np.nan):.4f} svd "
                f"{o.get('svd_sig24_med', np.nan):.4f} {o.get('error', '')[:300]} {time.time() - t0:.0f}s")
    T = f"{OUTD}/tables"; sfx = a.tag
    pd.DataFrame(res_o).sort_values(["fam", "netsta", "epoch"]).to_csv(f"{T}/series_svd{sfx}.csv", index=False)
    pd.concat(res_D, ignore_index=True).to_parquet(f"{T}/daily_svd{sfx}.parquet", index=False)
    pd.concat(res_I, ignore_index=True).to_csv(f"{T}/injection_svd{sfx}.csv", index=False)
    pd.concat([c for c in res_C if len(c)], ignore_index=True).to_csv(f"{T}/control_svd{sfx}.csv", index=False)
    keys = list(stk)
    np.savez_compressed(f"{T}/stacks_svd{sfx}.npz", keys=np.array(keys), tau=S.TAU,
                        **{f"St_{i}": stk[k]["St"] for i, k in enumerate(keys)},
                        **{f"Rf_{i}": stk[k]["Rf"] for i, k in enumerate(keys)},
                        **{f"days_{i}": stk[k]["days"].astype(str) for i, k in enumerate(keys)},
                        tS=np.array([stk[k]["tS"] for k in keys]))
    log(f"done {time.time() - t0:.0f}s")
    open(f"{OUTD}/SVD_MEASURE{sfx}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
