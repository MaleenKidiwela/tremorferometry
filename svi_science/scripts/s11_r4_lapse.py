#!/usr/bin/env python
"""s11_r4_lapse -- S11 round 4, task 3 (PLAN.md §8.3): does the joint-lapse chi2 pick the kernel that generated the data?
Synthetic data = the smooth model fitted to the real round-3 lapse cells with each generating kernel g (prediction, offsets
free) + noise (N1 Gaussian sigma*sqrt(max(1, chi2)); N2 residual bootstrap within station x window); 300 draws per g; selection
= argmin chi2 over V6 (as round 3). Confusion matrix, uniform-prior confidence, class K = {dchi2_data <= 4}; adoption rule;
Q1/Q2 under the adopted setup with the round-3 gates (s11_r3_kerr.rereport, outputs redirected to round 4) + posterior-
weighted kernel-error recovery from round-3 kerr_rules.csv; column F1 (near + far, M/B free).
    S11_REGION=... python s11_r4_lapse.py [--ndraw 300] [--nperm 500]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import sys, argparse, time, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_r4_common as R4
import s11_r3_common as R
import s11_r3_fit as LF

V6 = R4.V6


def chi_batch(L, Dm):
    """L: s11_r3_fit.Lin; Dm (n, K) data columns -> chi2 sums (K,)"""
    Th = L.Op @ Dm
    mu = np.zeros((len(L.Wg), Dm.shape[1])); np.add.at(mu, L.codes, L.w[:, None] * Dm); mu /= L.Wg[:, None]
    Rr = Dm - mu[L.codes] - L.Gt @ Th
    return (L.w[:, None] * Rr ** 2).sum(0)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--ndraw", type=int, default=300); ap.add_argument("--nperm", type=int, default=500); ap.add_argument("--test", action="store_true")
    a = ap.parse_args()
    lg = R4.log("lapse"); t0 = time.time()
    lsel = R4.lapse_sel()
    CL = R.lapse_cells()
    okall = np.ones(len(CL), bool)
    for v in V6:
        c = R.lapse_with_setup(CL, v); okall &= np.isfinite(c.a_d.to_numpy()) & np.isfinite(c.a_s.to_numpy())
    Ls, Fs = {}, {}
    for v in V6:
        c = R.lapse_with_setup(CL, v)[okall].reset_index(drop=True)
        Ls[v] = LF.smooth(c); Fs[v] = Ls[v].solve()
    base = Ls[V6[0]]
    n = len(base.w)
    dat = np.array([float((Fs[v]["wr2"]).sum()) for v in V6])
    dchi = dat - dat.min(); vstar = V6[int(np.argmin(dat))]
    K = [v for v, x in zip(V6, dchi) if x <= 4.0]
    lg(f"region {C.REGION}: lapse cells {n}; data chi2 per variant " + ", ".join(f"{v} {x:.1f}" for v, x in zip(V6, dchi)) + f"; L-SEL {vstar}; class K {K}")
    # station x window groups for the residual bootstrap
    cells0 = Ls[V6[0]].cells
    grp = pd.factorize(cells0.netsta + "|" + cells0.win.astype(str))[0]
    rng = np.random.default_rng(C.SEED + 44)
    conf = []
    for noise in ("N1", "N2", "N0ctrl"):
        for g in V6:
            F = Fs[g]; L = Ls[g]
            s = L.Gt @ F["th"] if noise != "N0ctrl" else np.zeros(n)
            sc = np.sqrt(max(1.0, F["chi2"]))
            Dm = np.zeros((n, a.ndraw))
            for k in range(a.ndraw):
                if noise in ("N1", "N0ctrl"):
                    e = rng.standard_normal(n) * L.cells.sig.to_numpy() * sc
                else:
                    r = F["resid"]; e = np.empty(n)
                    for gi in np.unique(grp):
                        ix = np.flatnonzero(grp == gi); e[ix] = r[rng.choice(ix, len(ix), replace=True)]
                Dm[:, k] = s + e
            chis = np.array([chi_batch(Ls[v], Dm) for v in V6])          # (6, ndraw)
            sel = np.array(V6)[np.argmin(chis, 0)]
            dd = chis - chis.min(0)
            for v in V6:
                conf.append(dict(noise=noise, true=g, sel=v, p=float(np.mean(sel == v))))
            inK = np.isin(sel, K)
            conf.append(dict(noise=noise, true=g, sel="K", p=float(np.mean(inK))))
            # how often the runner-up gap is at least the data's gap (strength of evidence), for the data's v*
            gap_data = np.sort(dchi)[1]
            gaps = np.sort(dd, 0)[1]
            conf.append(dict(noise=noise, true=g, sel="gap>=data", p=float(np.mean((sel == vstar) & (gaps >= gap_data)))))
            if noise == "N0ctrl":
                break                                              # noise-only control does not depend on g
            lg(f"  {noise} true {g}: " + ", ".join(f"{v} {np.mean(sel == v):.2f}" for v in V6) + f"; in K {np.mean(inK):.2f}")
        lg(f"  {noise} done {time.time() - t0:.0f}s")
    CF = pd.DataFrame(conf); CF.to_csv(f"{R4.TB4}/lapse_confusion.csv", index=False)
    ctrl = CF[(CF.noise == "N0ctrl") & CF.sel.isin(V6)]
    lg("noise-only control selection: " + ", ".join(f"{r.sel} {r.p:.2f}" for r in ctrl.itertuples()))
    res = dict(region=C.REGION, data_dchi2=dict(zip(V6, dchi.round(3).tolist())), vstar=vstar, K=K, written=R4.pt_now())
    ad = {}
    for noise in ("N1", "N2"):
        M = CF[(CF.noise == noise) & CF.sel.isin(V6)].pivot(index="true", columns="sel", values="p").reindex(index=V6, columns=V6).fillna(0)
        diag = float(M.loc[vstar, vstar]); off = float(M.drop(index=vstar)[vstar].max())
        post = float(M.loc[vstar, vstar] / M[vstar].sum()) if M[vstar].sum() > 0 else np.nan
        postall = (M[vstar] / M[vstar].sum()).to_dict() if M[vstar].sum() > 0 else {}
        kk = CF[(CF.noise == noise) & (CF.sel == "K")].set_index("true").p
        pK_in = float(kk[K].mean()); pK_out = float(kk.drop(index=K).max()) if len(set(V6) - set(K)) else 0.0
        postK = float(kk[K].sum() / kk.sum()) if kk.sum() > 0 else np.nan
        ad[noise] = dict(diag=diag, max_off=off, posterior_vstar=post, posterior_over_true=postall, pK_in=pK_in, pK_out_max=pK_out, posterior_K=postK,
                         single_ok=bool(diag >= 0.8 and off <= 0.2), class_ok=bool(pK_in >= 0.8 and pK_out <= 0.2))
        lg(f"{noise}: P(sel=v*|v*) {diag:.2f}, max P(sel=v*|g≠v*) {off:.2f}, posterior {post:.2f}; class K: in {pK_in:.2f}, out max {pK_out:.2f}, "
           f"posterior {postK:.2f}")
    if all(ad[k]["single_ok"] for k in ad):
        adopted = vstar; setup = vstar; conf_ = min(ad[k]["posterior_vstar"] for k in ad)
    elif all(ad[k]["class_ok"] for k in ad) and len(K) > 1:
        adopted = "K:" + "+".join(K); setup = ("avg", K); conf_ = min(ad[k]["posterior_K"] for k in ad)
    elif all(ad[k]["class_ok"] for k in ad):
        adopted = vstar; setup = vstar; conf_ = min(ad[k]["posterior_K"] for k in ad)
    else:
        adopted = None; setup = None; conf_ = np.nan
    res.update(adoption=ad, adopted=adopted, confidence=conf_)
    lg(f"ADOPTED for Q1: {adopted} (confidence {conf_:.2f})")
    if setup is not None:
        import s11_r3_kerr as KR
        cells = R.head_cells(); cctrl = R.head_cells(control=True); R4.coef_table()
        old_tb3 = R.TB3; R.TB3 = R4.TB4                   # rereport writes kerr_<label>_gates.csv / _q3.csv -> round 4
        try:
            G = KR.rereport("ADOPT", setup, cells, cctrl, CL, R.majors(), a.nperm, lg, "lapse_adopted")
        finally:
            R.TB3 = old_tb3
        q = G[G.item == "Q1Q2"].iloc[0]
        res.update(Q1=float(q.Q1), Q1_sig=float(q.Q1_sig), Q1_peak=float(q.Q1_peak), Q2=float(q.Q2), Q2_sig=float(q.Q2_sig))
        res["gates"] = {r.item: (None if pd.isna(r.passed) else bool(r.passed)) for r in G.itertuples() if "passed" in G and not pd.isna(getattr(r, "passed", np.nan))}
        # posterior-weighted kernel-error recovery (round-3 rules for the matching setup)
        kr = pd.read_csv(f"{C.ROOT}/svi_science/results/s11_depth/round3/{R4.SUB}/tables/kerr_rules.csv")
        name = {"ell40": "NOM"}.get(adopted, "L-SEL" if adopted == lsel else None)
        if name is not None and name in set(kr.setup):
            k_ = kr[kr.setup == name].set_index("gen")
            post = ad["N1"]["posterior_over_true"]
            wv = pd.Series(post).reindex(V6).fillna(0)
            res["kerr_expected"] = {col: float(np.nansum(wv.to_numpy() * k_.reindex(V6)[col].to_numpy())) for col in ("R1_Q2rec", "R2_Q1leak", "R3_Bleak")}
            res["kerr_insample"] = {col: float(k_.loc[adopted, col]) for col in ("R1_Q2rec", "R2_Q1leak", "R3_Bleak", "R4_Pdrec", "R5_Psleak", "R6_B_Q3comp")} if adopted in k_.index else {}
            res["kerr_range"] = {col: [float(k_.reindex(V6)[col].min()), float(k_.reindex(V6)[col].max())] for col in ("R1_Q2rec", "R3_Bleak")}
        # F1 column (near + far B0 cells, M/B free), Q1 under the adopted setup
        try:
            import s11_dist_model as Mo
            fc = Mo.all_cells()
            fc[["a_s", "a_m", "a_d", "a_b"]] = R4.coefs(fc, setup)
            fc = fc[np.isfinite(fc.a_d) & np.isfinite(fc.a_s)].reset_index(drop=True)
            qf = Mo.Qs(Mo.fit(fc, "smooth", MB=True))
            res["F1_column"] = dict(Q1=qf["Q1"], Q1_sig=qf["Q1_sig"], Q1_peak=qf["Q1_peak"], Q2=qf["Q2"], Q2_sig=qf["Q2_sig"], n_cells=len(fc),
                                    n_far=int(fc.far.sum()))
            lg(f"F1 column under {adopted}: " + R4.fmt(res["F1_column"]))
        except Exception as e:
            import traceback
            lg(f"F1 column failed: {e} {traceback.format_exc()[-400:]}")
    R4.dump(res, f"{R4.TB4}/lapse_verdict.json")
    lg(f"lapse done {time.time() - t0:.0f}s")
    if not a.test:
        open(f"{R4.R4}/LAPSE_{C.REGION}_DONE", "w").write(R4.pt_now() + "\n")


if __name__ == "__main__":
    main()
