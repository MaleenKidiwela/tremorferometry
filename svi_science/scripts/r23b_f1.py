#!/usr/bin/env python
"""r23b_f1 (2026-10-11) -- Finding 1 on the 23-year SVI record: patch strength vs front reach (Q1), persistence /
turn-taking (S3 rule 7, Q1 H2), neighbours move together (episode_rates share test + Mantel).
Rules: svi_science/results/r23b/PLAN.md §2. Code reused by import (q1_tests, s3_patches, episode_rates), unchanged.

    python svi_science/scripts/r23b_f1.py [fast]
"""
from __future__ import annotations
import importlib.util, json, os, sys, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
SC = f"{ROOT}/svi_science/scripts"
OUT = f"{ROOT}/svi_science/results/r23b/tables/f1"
os.makedirs(OUT, exist_ok=True)
FAST = len(sys.argv) > 1 and sys.argv[1] == "fast"
sys.argv = sys.argv[:1]                       # q1_tests reads sys.argv for NPERM/NCAL
os.environ["EPR_OUT"] = OUT


def imp(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


Q = imp("q1t", f"{SC}/q1_tests.py")
S = imp("s3p", f"{SC}/s3_patches.py")
EPR = imp("epr", f"{ROOT}/lfe_extension_wa/scripts/episode_rates.py")
if FAST:
    Q.NPERM = 50; S.NPERM = 50
TB = f"{ROOT}/svi_science/results/r23b/tables"
LOGF = open(f"{ROOT}/svi_science/results/r23b/logs/r23b_f1{'_fast' if FAST else ''}.log", "w")


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"; print(s, flush=True); LOGF.write(s + "\n"); LOGF.flush()


def load():
    W = pd.read_csv(f"{TB}/windows.csv", parse_dates=["w0", "w1"])
    PE = pd.read_parquet(f"{TB}/patch_ep.parquet")
    FE = pd.read_parquet(f"{TB}/fam_ep.parquet")
    TR = pd.read_csv(f"{TB}/tremor_ep.csv")
    DT = pd.read_csv(f"{TB}/det_patch_ep.csv")
    G = pd.read_csv(f"{ROOT}/svi_science/results/s0/groups.csv")
    CL = []
    for v, col in [("ceq", "C_eq"), ("once", "n_once"), ("kc", "C_kc")]:
        CL.append(pd.DataFrame(dict(catalogue="SVI", variant=v, patchdef="group", patch=PE.group.astype(int).astype(str),
                                    episode=PE.episode, n=PE[col].astype(float), t_med=PE.t_med)))
    CL = pd.concat(CL, ignore_index=True)
    PA = pd.DataFrame(dict(catalogue="SVI", patchdef="group", patch=G.group.astype(str), as_km=G.along_strike_km))
    DT = DT.assign(catalogue="SVI", patchdef="group", patch=DT.group.astype(int).astype(str))
    GEO = W[["episode", "gnss_D_mm"]].rename(columns={"gnss_D_mm": "D"})
    D = dict(CL=CL, PA=PA, TR=TR, DT=DT, GEO=GEO, KIN=pd.DataFrame(columns=["region", "episode", "direction"]))
    return W, PE, FE, D


PER = {"2003-13": ("2003-01-01", "2013-12-31"), "2014-26": ("2014-01-01", "2026-12-31"), "2003-26": ("2003-01-01", "2026-12-31")}


def design(D, W, variant, period, pool="ALL", targets=None, hist_from=None):
    a, b = PER[period]
    H = W[(W.start >= (hist_from or a)) & (W.start <= b)].copy()
    H["catalogue"] = "SVI"
    H["era"] = "ALL" if pool == "ALL" else H.seggroup
    if targets is None:
        targets = H[H.big].episode.tolist()
    des = Q.Design(D, "SVI", variant, "group", H[["catalogue", "episode", "w0", "w1", "era", "major"]], [("SVI", e) for e in targets])
    return des


def s3_tests(C, mid):
    allc, cons, seps = S.pair_stats(C, mid)
    exc, R = S.excursion(C)
    pm = []
    for _ in range(S.NPERM):
        Cp = np.stack([S.rng.permutation(C[:, k]) for k in range(C.shape[1])], 1)
        pm.append(S.pair_stats(Cp, mid)[0].mean())
    pm = np.array(pm)
    mm, mc = [], []
    for Sm in S.multinomial_ctrl(C, 200 if FAST else 500):
        a, c, _ = S.pair_stats(Sm.astype(float), mid); mm.append(a.mean()); mc.append(c.mean() if len(c) else np.nan)
    mm = np.array(mm)
    oc, oe = [], []
    for _ in range(S.NPERM):
        o = S.rng.permutation(C.shape[1])
        _, c, _ = S.pair_stats(C[:, o], mid)
        oc.append(c.mean() if len(c) else np.nan); oe.append(S.excursion(C[:, o])[0])
    oc, oe = np.array(oc), np.array(oe)
    t = dict(n_patches=int(C.shape[0]), n_episodes=int(C.shape[1]), n_pairs_ge90d=int(len(allc)), mean_resid_corr=float(allc.mean()),
             labelshuffle_mean=float(pm.mean()), p_persist_labelshuffle=float((np.sum(pm >= allc.mean()) + 1) / (len(pm) + 1)),
             multinomial_mean=float(mm.mean()), p_persist_multinomial=float((np.sum(mm >= allc.mean()) + 1) / (len(mm) + 1)),
             n_consecutive_pairs=int(len(cons)), mean_consecutive_corr=float(cons.mean()) if len(cons) else None,
             ordershuffle_consecutive_mean=float(np.nanmean(oc)),
             p_catchup_ordershuffle=float((np.sum(oc <= cons.mean()) + 1) / (np.isfinite(oc).sum() + 1)) if len(cons) else None,
             excursion=exc, ordershuffle_excursion_q05=float(np.quantile(oe, 0.05)),
             p_excursion_low=float((np.sum(oe <= exc) + 1) / (len(oe) + 1)),
             chi2_overdispersion_per_df=float((R ** 2).sum() / ((C.shape[0] - 1) * (C.shape[1] - 1))))
    t["persistent"] = bool(t["mean_resid_corr"] > 0 and t["p_persist_labelshuffle"] < 0.01)
    t["catch_up"] = bool(len(cons) and t["mean_consecutive_corr"] < 0 and t["p_catchup_ordershuffle"] < 0.01 and exc < t["ordershuffle_excursion_q05"])
    return t


def seg_cluster(des, W, nperm=10000):
    """post-hoc (Revision 1): residual correlation of target pairs in the same vs different network segment group."""
    sg = W.set_index("episode").seggroup
    labs = np.array([sg[des.keys[j][1]] for j in des.tslots]); R = des.rL
    n = R.shape[1]; i, j = np.triu_indices(n, 1)
    C = np.array([np.corrcoef(R[:, a], R[:, b])[0, 1] for a, b in zip(i, j)])
    def stat(l):
        same = l[i] == l[j]
        return C[same].mean() - C[~same].mean() if same.any() and (~same).any() else np.nan
    obs = stat(labs); rng = np.random.default_rng(5)
    nul = np.array([stat(rng.permutation(labs)) for _ in range(nperm)])
    return dict(same_minus_diff=float(obs), p=float((np.sum(nul >= obs) + 1) / (nperm + 1)), n_pairs=int(len(C)),
                mean_same=float(C[labs[i] == labs[j]].mean()), mean_diff=float(C[labs[i] != labs[j]].mean()),
                n_out_of_reach_cells=None)


def partial_r2(des, df0, K):
    df = df0.sort_values(["p", "ti"])
    r = df.r.to_numpy(); g = df.ti.to_numpy()
    T1 = np.where(df.cens, np.nan, np.log10(df.T1))
    cols = dict(T1=T1, T2=df.T2.to_numpy(), T3=df.T3.to_numpy(), in_reach=K["inr"].reshape(-1), km_outside=K["kout"].reshape(-1),
                dir_mode=K["dm"].reshape(-1), u_along_path=K["u"].reshape(-1), detectability=K["det"].reshape(-1))
    out = []
    for k, x in cols.items():
        b_, r2, n = Q.fe_slope(r, x, g); out.append(dict(predictor=k, beta=b_, r2_within=r2, n=n))
    kin = np.column_stack([cols["in_reach"], cols["km_outside"], np.nan_to_num(cols["dir_mode"]), np.nan_to_num(cols["u_along_path"], nan=0.5)])
    bud = np.column_stack([np.nan_to_num(T1, nan=np.nanmedian(T1)), np.nan_to_num(cols["T2"]), np.nan_to_num(cols["T3"])])
    out.append(dict(predictor="kinematics_set", r2_within=Q.fe_multi_r2(r, kin, g), n=len(r)))
    out.append(dict(predictor="budget_set", r2_within=Q.fe_multi_r2(r, bud, g), n=len(r)))
    out.append(dict(predictor="all", r2_within=Q.fe_multi_r2(r, np.column_stack([kin, bud, np.nan_to_num(cols["detectability"])]), g), n=len(r)))
    return out


def run(name, des, D, full=True, per_patch=False, kin=True, s3=True, h4perm=None):
    res = dict(n_targets=int(len(des.tslots)), n_history=int(len(des.H)), n_patches=len(des.patches),
               targets=[des.keys[j][1] for j in des.tslots], N_targets=[float(des.N[j]) for j in des.tslots])
    t0 = time.time()
    b, df0 = Q.run_budget(des, per_patch=per_patch, nperm=Q.NPERM if full else min(Q.NPERM, 500))
    res["budget"] = b; res["n_censored"] = int(df0.cens.sum())
    df0["patch"] = [des.patches[p] for p in df0.p]; df0["episode"] = [des.keys[j][1] for j in df0.slot]
    df0.to_csv(f"{OUT}/rows_{name}.csv", index=False)
    if kin:
        cov = Q.covariates(des, D, "SVI", "tremor")
        kt, K = Q.kin_tests(des, cov, nperm=Q.NPERM if full else min(Q.NPERM, 500))
        res["kin"] = kt
        res["H4_L2"] = Q.h4_and_l2(des, cov, K, nperm=h4perm or (Q.NPERM if full else min(Q.NPERM, 500)))
        e, ev = Q.episode_level(des, cov, D, "SVI")
        res["episode_level"] = e; res["episode_level_values"] = ev
        res["partial_r2"] = partial_r2(des, df0, dict(K, det=cov["det"]))
        np.savez(f"{OUT}/kin_{name}.npz", **K, det=cov["det"], r=des.rL, patches=np.array(des.patches), eps=np.array(cov["eps"]), pos=des.pos)
    if s3:
        mid = (des.H.w0 + (des.H.w1 - des.H.w0) / 2).to_numpy()[des.tslots]
        res["s3"] = s3_tests(des.C[:, des.tslots].astype(float), mid)
    bb = {k: (round(v["beta"], 3) if np.isfinite(v["beta"]) else None, v["p_order"], v.get("p_perpatch")) for k, v in b.items()}
    kk = {k: (round(v["beta"], 3), round(v["r2"], 3), v["p"]) for k, v in res.get("kin", {}).items()}
    log(f"{name} [{time.time() - t0:.0f}s] n_t={res['n_targets']} budget {bb} kin {kk} "
        f"H4/L2 {({k: round(v, 3) if isinstance(v, float) else v for k, v in res.get('H4_L2', {}).items()})} "
        f"S3 persist {res.get('s3', {}).get('mean_resid_corr')} p {res.get('s3', {}).get('p_persist_labelshuffle')} "
        f"consec {res.get('s3', {}).get('mean_consecutive_corr')} p {res.get('s3', {}).get('p_catchup_ordershuffle')} "
        f"chi2/df {res.get('s3', {}).get('chi2_overdispersion_per_df')}")
    return res


def neighbours(FE, W, F):
    fams = F[["fid", "lat", "lon"]].copy()
    fids = fams.fid.astype(int).tolist()
    J = EPR.set_jaccard_svi(fids)
    out, bins_all = {}, []
    maj = W[W.major]
    blocks = [("C_eq|2003-13", "C_eq", maj[maj.start < "2014-01-01"].episode),
              ("C_eq|2014-26", "C_eq", maj[maj.start >= "2014-01-01"].episode),
              ("C_eq|2003-26", "C_eq", maj.episode),
              ("once|S15f_2010-13", "n_once", maj[maj.seggroup == "S15f"].episode),
              ("once|FWD_2014-26", "n_once", maj[maj.seggroup == "FWD"].episode),
              ("once|2003-26_crossnetwork", "n_once", maj.episode)]
    for key, col, eps in blocks:
        T = FE[FE.episode.isin(eps)].rename(columns={"fam": "fid"})
        FT, RT, dist = EPR.share_tests(T, fams, key, ncol=col, nboot=500 if FAST else 2000, seed=11)
        t = FT[FT.tested]
        mt, pairs = EPR.mantel(RT, fids, dist, J, nperm=100 if FAST else 1000)
        d = dict(n_episodes=int(len(set(eps))), n_tested=int(len(t)), n_overdispersed_q05=int((t.q_bh < 0.05).sum()),
                 median_chi2_df=float(t.chi2_df.median()) if len(t) else None,
                 coef_logdist=mt.get("coef_logdist"), p_logdist=mt.get("p_logdist"), coef_jaccard=mt.get("coef_jaccard"),
                 p_jaccard=mt.get("p_jaccard"), n_pairs=mt.get("n_pairs"))
        if pairs is not None:
            pairs["b"] = pd.cut(pairs.dist_km, [0, 3, 6, 10, 20, 40, 80, 400])
            gb = pairs.groupby("b", observed=True)["corr"].agg(["mean", "size", "std"]).reset_index()
            gb["se"] = gb["std"] / np.sqrt(gb["size"]); gb["block"] = key; gb["b"] = gb.b.astype(str)
            bins_all.append(gb)
            d["corr_bins"] = {r.b: [round(float(r["mean"]), 3), int(r["size"])] for _, r in gb.iterrows()}
            m36 = gb.set_index("b")["mean"]
            d["neighbours_move_together"] = bool(d["coef_logdist"] < 0 and d["p_logdist"] < 0.01 and m36.get("(3, 6]", -1) > 0 and m36.get("(6, 10]", -1) > 0)
            pairs.drop(columns=["b", "dist_bin"], errors="ignore").assign(block=key).to_csv(f"{OUT}/pairs_{key.replace('|', '_')}.csv", index=False)
        FT.to_csv(f"{OUT}/families_{key.replace('|', '_')}.csv", index=False)
        out[key] = d
        log(f"neighbours {key}: {json.dumps(d, default=str)}")
    pd.concat(bins_all).to_csv(f"{OUT}/neighbour_corr_bins.csv", index=False)
    return out


def main():
    W, PE, FE, D = load()
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
    R = {}
    s3big = json.load(open(f"{ROOT}/svi_science/results/s3/build_info.json"))["big"]["SVI"]
    # R0 reproduction (Q1 / S3)
    R["R0_once_S3big_2014-26"] = run("R0", design(D, W, "once", "2014-26", targets=s3big), D, full=True)
    # headline C_eq
    for per in ["2003-26", "2003-13", "2014-26"]:
        R[f"A_ceq_{per}"] = run(f"A_ceq_{per}", design(D, W, "ceq", per), D, full=True, per_patch=(per == "2003-26"))
    # columns
    for per in ["2003-26", "2003-13"]:
        R[f"C_once_segpool_{per}"] = run(f"C_once_segpool_{per}", design(D, W, "once", per, pool="seg"), D, full=False)
        R[f"C_kc_{per}"] = run(f"C_kc_{per}", design(D, W, "kc", per), D, full=False)
    big = W[W.big]
    ge300 = big[big.n_once >= 300].episode.tolist()
    R["C_ceq_ge300_2003-26"] = run("C_ceq_ge300_2003-26", design(D, W, "ceq", "2003-26", targets=ge300), D, full=False)
    for per in ["2003-13", "2003-26"]:
        des_ = design(D, W, "ceq", per); cov_ = Q.covariates(des_, D, "SVI", "tremor"); inr_, _, _, _ = Q.kin_rows(des_, cov_)
        R[f"P_segcluster_{per}"] = seg_cluster(des_, W); R[f"P_segcluster_{per}"]["n_out_of_reach_cells"] = int((inr_ == 0).sum())
        log(f"P_segcluster_{per} {R[f'P_segcluster_{per}']}")
    R["C_ceq_2014-26_S3big"] = run("C_ceq_2014-26_S3big", design(D, W, "ceq", "2014-26", targets=s3big), D, full=False)
    # env tremor for 2010-13 as well (kinematics only)
    TRe = D["TR"].copy()
    m = (TRe.episode.map(W.set_index("episode").start) >= "2010-01-06") & (TRe.episode.map(W.set_index("episode").start) < "2014-01-01")
    for c in ["n_trem", "trem_p05", "trem_p95", "trem_speed", "trem_dir"]:
        TRe.loc[m, c] = TRe.loc[m, f"{c}_env"]
    De = dict(D); De["TR"] = TRe
    for per in ["2003-13", "2003-26"]:
        des = design(De, W, "ceq", per)
        cov = Q.covariates(des, De, "SVI", "tremor"); kt, K = Q.kin_tests(des, cov, nperm=min(Q.NPERM, 500))
        e, _ = Q.episode_level(des, cov, De, "SVI")
        R[f"C_ceq_envtremor_{per}"] = dict(kin=kt, episode_level=e)
        log(f"C_ceq_envtremor_{per} kin {({k: (round(v['beta'], 3), round(v['r2'], 3), v['p']) for k, v in kt.items()})} E1 {e.get('trem_width')}")
    # calibration (no-memory multinomial) for the joint headline design
    des = design(D, W, "ceq", "2003-26")
    des.C = np.round(des.C)
    Q._recompute_resid(des)
    R["calibration_A_ceq_2003-26"] = Q.calibrate(des, nsim=20 if FAST else 100, nperm=50 if FAST else 200)
    log(f"calibration {R['calibration_A_ceq_2003-26']}")
    # neighbours
    R["neighbours"] = neighbours(FE, W, F)
    json.dump(R, open(f"{OUT}/f1_results.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    # flat table
    rows = []
    for k, v in R.items():
        if not isinstance(v, dict):
            continue
        for sec in ("budget", "kin"):
            for p, x in v.get(sec, {}).items():
                rows.append(dict(design=k, section=sec, predictor=p, beta=x.get("beta"), r2=x.get("r2"), n=x.get("n"),
                                 p=x.get("p", x.get("p_order")), p_perpatch=x.get("p_perpatch")))
        for p, x in v.get("episode_level", {}).items():
            rows.append(dict(design=k, section="E1", predictor=p, beta=x.get("rho"), n=x.get("n"), p=x.get("p")))
        if "s3" in v:
            s = v["s3"]
            rows.append(dict(design=k, section="S3", predictor="persistence_mean_resid_corr", beta=s["mean_resid_corr"], p=s["p_persist_labelshuffle"]))
            rows.append(dict(design=k, section="S3", predictor="consecutive_corr", beta=s["mean_consecutive_corr"], p=s["p_catchup_ordershuffle"]))
            rows.append(dict(design=k, section="S3", predictor="chi2_per_df", beta=s["chi2_overdispersion_per_df"]))
        if "H4_L2" in v:
            h = v["H4_L2"]
            rows.append(dict(design=k, section="H4", predictor="H4_dR2", beta=h["H4_dR2"], p=h["p_H4"], r2=h["R2_full"]))
            rows.append(dict(design=k, section="L2", predictor="L2", beta=h["L2_beta"], p=h["p_L2"]))
        for x in v.get("partial_r2", []):
            rows.append(dict(design=k, section="partialR2", predictor=x["predictor"], beta=x.get("beta"), r2=x["r2_within"], n=x["n"]))
    pd.DataFrame(rows).to_csv(f"{OUT}/f1_table.csv", index=False)
    log("done")


if __name__ == "__main__":
    main()
