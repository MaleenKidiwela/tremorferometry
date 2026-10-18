#!/usr/bin/env python
"""s11_r3_q3 -- S11 round 3, task 2 (PLAN §6.2): Q3 = DeltaD(e) with the REDPy estimators.
Configurations: B0 (raw path, P weights, n>=20; recomputed), B0n50, B1 (cc weights, reference path, n>=20), B2 (cc-weighted
all-pairs graph, >=50-event episode nodes, bootstrap sigma), B2p (P-weighted graph). Gates per configuration (§4 rules):
acceptance, scrambled control, leave-one-site-out, family split, prior x2/x0.5, deep pulse and common shallow pulse injected
into the estimator's own day stacks and run through that estimator. Adoption rule fixed in PLAN §6.2.
    S11_REGION=... python s11_r3_q3.py [--workers 12]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_r3_common as R
import s11_r3_fit as LF
import s11_r3_inject as J
import s11_paths as P

CONFIGS = ["B0", "B0n50", "B1", "B2", "B2p"]


def r3_tables():
    E = [pd.read_parquet(f"{R.TB3}/r3_episode.parquet")]
    if C.REGION == "SVI" and os.path.exists(f"{R.TB3}/r3_episode_ext.parquet"):
        E.append(pd.read_parquet(f"{R.TB3}/r3_episode_ext.parquet"))
    E = pd.concat(E, ignore_index=True)
    B = pd.read_parquet(f"{C.TB}/episode_raw{R.TAG}.parquet")[["fam", "netsta", "epoch", "episode", "cc3c_25"]]
    return E.merge(B, on=["fam", "netsta", "epoch", "episode"], how="left")


def inflate(cells, ctrl):
    """per series sqrt(max(1, mean (d/sig)^2 of its control values)), >= 3 values; else station median; else 1"""
    c = ctrl[np.isfinite(ctrl.dc) & np.isfinite(ctrl.sc) & (ctrl.sc >= 0.01)].copy()
    c["z2"] = (c.dc / c.sc) ** 2
    g = c.groupby(["fam", "netsta", "epoch"]).z2.agg(["mean", "size"])
    gs = g[g["size"] >= 3].reset_index(); st = gs.groupby("netsta")["mean"].median(); st0 = c.groupby("netsta").z2.median()
    out = []
    for f, n, e in zip(cells.fam, cells.netsta, cells.epoch):
        if (f, n, e) in g.index and g.loc[(f, n, e), "size"] >= 3:
            chi = g.loc[(f, n, e), "mean"]
        else:
            chi = st.get(n, st0.get(n, 1.0))
        out.append(float(np.sqrt(max(1.0, chi))))
    return np.array(out)


CTRL = {"B1": ("b1_c2_eps", "b1_c2_sig"), "B2": ("g_c2_eps", "g_c2_sig"), "B2p": ("gp_c2_eps", "gp_c2_sig")}


def build(E, cfg):
    """-> (cells, control cells) for a round-3 configuration (NOM kernel attached)"""
    ce, cs_ = CTRL[cfg]
    Ec = E.copy()
    Ec["dc"] = -100 * Ec[ce] if ce in Ec else np.nan
    Ec["sc"] = 100 * Ec[cs_] if cs_ in Ec else np.nan
    if cfg == "B1":
        q = Ec[(Ec.n >= R.NMIN) & (Ec.b1_cc >= R.CCMIN)].copy(); q["d"] = -100 * q.b1_eps; q["sig"] = 100 * q.b1_sig; q["cc"] = q.b1_cc
    elif cfg == "B2":
        q = Ec[(Ec.n >= 50) & (Ec.b1_cc >= R.CCMIN)].copy(); q["d"] = -100 * q.g_eps; q["sig"] = 100 * q.g_sig; q["cc"] = q.b1_cc
    elif cfg == "B2p":
        q = Ec[(Ec.n >= 50) & (Ec.cc3c_25 >= R.CCMIN)].copy(); q["d"] = -100 * q.gp_eps; q["sig"] = 100 * q.gp_sig; q["cc"] = q.cc3c_25
    else:
        raise ValueError(cfg)
    q["site"] = q.netsta.map(C.SITE)
    q = q[np.isfinite(q.d) & np.isfinite(q.sig) & (q.sig > 0) & q.site.notna()].reset_index(drop=True)
    q["sig_raw"] = q.sig; q["infl"] = inflate(q, Ec); q["sig"] = q.sig_raw * q.infl
    key = q.fam.astype(str) + "|" + q.netsta + "|" + q.epoch
    q = q[key.map(key.value_counts()) >= 3].reset_index(drop=True)
    q["series"] = q.fam.astype(str) + "|" + q.netsta + "|" + q.epoch
    q = R.with_setup(q, "ell40"); q = q[np.isfinite(q.a_d) & np.isfinite(q.a_s)].reset_index(drop=True)
    qc = q[np.isfinite(q.dc) & np.isfinite(q.sc) & (q.sc > 0)].copy()
    qc["d"] = qc.dc; qc["sig"] = qc.sc * qc.infl
    return q, qc.reset_index(drop=True)


def cells_for(cfg, E):
    if cfg == "B0":
        return R.head_cells(), R.head_cells(control=True)
    if cfg == "B0n50":
        return R.head_cells(nmin=50), R.head_cells(nmin=50, control=True)
    return build(E, cfg)


def gates(cfg, cells, ctrl, majors, lg, setup="ell40"):
    cs = R.with_setup(cells, setup) if setup != "ell40" else cells
    cs = cs[np.isfinite(cs.a_d)].reset_index(drop=True)
    F = LF.temporal(cs).solve(); T, comp = LF.q3_table(F, majors)
    if comp is None or not len(T):
        return dict(config=cfg, setup=R.setup_name(setup), n_cells=len(cs), note="no measurable major"), T
    Tm = T.set_index("episode"); mj = T[T.major]
    G = dict(config=cfg, setup=R.setup_name(setup), n_cells=len(cs), n_series=F["nser"], n_episodes=len(F["eps"]), chi2=F["chi2"],
             sig_major_median=float(mj.sig.median()), sig_major_min=float(mj.sig.min()), sig_major_max=float(mj.sig.max()),
             n_major_sig_le_0p1=int((mj.sig <= 0.10).sum()), n_majors=int(len(mj)), composite=comp["composite"],
             composite_sig=comp["composite_sig"], chi2_majors=float(np.mean((mj.dD / mj.sig) ** 2)),
             n_beyond_2sig=int((np.abs(mj.dD / mj.sig) > 2).sum()), n_beyond_3sig=int((np.abs(mj.dD / mj.sig) >= 3).sum()))
    G["acceptance"] = bool(G["n_major_sig_le_0p1"] >= 5)
    if setup != "ell40":
        return G, T
    mv = []
    for f_ in (2.0, 0.5):
        Tf, cf = LF.q3_table(LF.temporal(cs, prior=f_).solve(), majors); Tf = Tf.set_index("episode")
        e_ = [e for e in mj.episode if e in Tf.index]
        if not e_ or cf is None:
            continue
        G[f"prior_x{f_:g}"] = float(np.max(np.abs((Tf.dD[e_] - Tm.dD[e_]) / Tm.sig[e_])))
        G[f"prior_x{f_:g}_comp"] = abs(cf["composite"] - comp["composite"]) / comp["composite_sig"]
    G["prior_pass"] = bool(max(G.get("prior_x2", 9), G.get("prior_x0.5", 9), G.get("prior_x2_comp", 9), G.get("prior_x0.5_comp", 9)) < 0.1)
    cc = ctrl[np.isfinite(ctrl.a_d)].reset_index(drop=True)
    if len(cc) >= 20:
        Tc, ccomp = LF.q3_table(LF.temporal(cc).solve(), majors); Tc = Tc.set_index("episode")
        meas = [e for e in mj.episode if e in Tc.index and Tc.n_cells[e] >= 5]
        frac = float(np.mean([abs(Tc.dD[e]) < 2 * Tc.sig[e] for e in meas])) if meas else np.nan
        G.update(control_frac=frac, control_n=len(meas), control_cells=len(cc), control_comp=ccomp["composite"] if ccomp else np.nan,
                 control_comp_sig=ccomp["composite_sig"] if ccomp else np.nan,
                 control_pass=bool(frac >= 0.9 and (ccomp is None or abs(ccomp["composite"]) < 2 * ccomp["composite_sig"])) if meas else None)
        Tc.reset_index().to_csv(f"{R.TB3}/r3q3_control_{cfg}.csv", index=False)
    else:
        G.update(control_cells=len(cc), control_pass=None)
    worst, worst_site = 1.0, ""
    for s in sorted(cs.site.unique()):
        sub = cs[cs.site != s].reset_index(drop=True)
        if sub.episode.nunique() < 3:
            continue
        Ts, _ = LF.q3_table(LF.temporal(sub).solve(), majors); Ts = Ts.set_index("episode")
        ee = [e for e in mj.episode if e in Ts.index]
        if ee:
            fr = float(np.mean([abs(Ts.dD[e] - Tm.dD[e]) < 2 * Ts.sig[e] for e in ee]))
            if fr < worst:
                worst, worst_site = fr, s
    G.update(loso_worst=worst, loso_worst_site=worst_site, loso_pass=bool(worst >= 0.9))
    To, co = LF.q3_table(LF.temporal(cs[cs.fam % 2 == 1].reset_index(drop=True)).solve(), majors)
    Te, ce = LF.q3_table(LF.temporal(cs[cs.fam % 2 == 0].reset_index(drop=True)).solve(), majors)
    if co and ce:
        J_ = To.set_index("episode").join(Te.set_index("episode"), lsuffix="_o", rsuffix="_e", how="inner"); J_ = J_[J_.major_o]
        fr = float(np.mean(np.abs(J_.dD_o - J_.dD_e) < 2 * np.hypot(J_.sig_o, J_.sig_e))) if len(J_) else np.nan
        G.update(famsplit_frac=fr, famsplit_odd=co["composite"], famsplit_odd_sig=co["composite_sig"], famsplit_even=ce["composite"],
                 famsplit_even_sig=ce["composite_sig"],
                 famsplit_pass=bool(fr >= 0.9 and abs(co["composite"] - ce["composite"]) < 2 * np.hypot(co["composite_sig"], ce["composite_sig"])))
    return G, T


def pulse_runs(cs_coef, majors):
    def fn(key, st, unit):
        if key not in cs_coef:
            return {"nat": None}
        a_s, a_d = cs_coef[key]
        ep = np.array([unit.get(d, None) for d in st["days"]], dtype=object)
        runs = {"nat": None}
        for e in majors:
            m = ep == e
            if m.any():
                runs[f"Pd|{e}"] = -(0.003 * a_d) * m; runs[f"Ps|{e}"] = -(0.001 * a_s) * m
        mall = np.isin(ep, majors)
        if mall.any():
            runs["Psall"] = -(0.001 * a_s) * mall
        return runs
    return fn


def pulses(cfg, cells, vals, majors, lg):
    dn = J.cell_vector(vals, cells, "nat")
    lg(f"  [{cfg}] natural-from-store vs measured: max |diff| {np.nanmax(np.abs(dn - cells.d.to_numpy())):.1e} %, missing {int(np.isnan(dn).sum())}")
    ok = np.isfinite(dn); cs = cells[ok].reset_index(drop=True); dn = dn[ok]
    Lt = LF.temporal(cs); F0 = Lt.delta(np.zeros(len(cs)))
    rows = []
    for e in majors:
        g = LF.dvec(F0, e)
        if g is None:
            continue
        for kind, amp in (("deep", 0.3), ("shal", 0.1)):
            di = J.cell_vector(vals, cells, f"{'Pd' if kind == 'deep' else 'Ps'}|{e}")[ok]
            x = np.where(np.isfinite(di), di - dn, 0.0)
            v = float(g @ (Lt.Op @ x)); rows.append(dict(config=cfg, episode=e, kind=kind, injected=amp, response=v, ratio=v / amp))
    gc = LF.comp_vec(F0, majors)
    di = J.cell_vector(vals, cells, "Psall")[ok]; x = np.where(np.isfinite(di), di - dn, 0.0)
    if gc is not None:
        v = float(gc @ (Lt.Op @ x)); rows.append(dict(config=cfg, episode="all", kind="shal", injected=0.1, response=v, ratio=v / 0.1))
    T = pd.DataFrame(rows)
    dp = T[(T.kind == "deep")].ratio; sp = T[(T.kind == "shal") & (T.episode != "all")].ratio.abs()
    out = dict(deep_median=float(dp.median()), deep_lo=float(dp.min()), deep_hi=float(dp.max()), deep_pass=bool(0.8 <= dp.median() <= 1.1),
               shal_median=float(sp.median()), shal_lo=float(sp.min()), shal_hi=float(sp.max()), shal_pass=bool(sp.median() <= 0.35),
               shal_all_composite=float(T[T.episode == "all"].ratio.iloc[0]) if (T.episode == "all").any() else np.nan)
    lg(f"  [{cfg}] pulses: deep median {out['deep_median']:.3f} ({out['deep_lo']:.2f}-{out['deep_hi']:.2f}), shallow median {out['shal_median']:.3f}, "
       f"all-majors shallow composite {out['shal_all_composite']:.3f}")
    return out, T


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); a = ap.parse_args()
    lg = R.log("r3q3"); t0 = time.time()
    majors = R.majors()
    E = r3_tables()
    vj = f"{R.TB3}/kerr_verdict.json"
    win3 = json.load(open(vj))["winner"].get("q3") if os.path.exists(vj) else None
    setups_ = {}
    if win3:
        import s11_r3_kerr as KR
        lap = json.load(open(vj))["lapse"]
        setups_ = KR.setups(lap)
    summ, pul, allT = [], [], []
    CEL = {}
    for cfg in CONFIGS:
        cells, ctrl = cells_for(cfg, E)
        CEL[cfg] = cells
        G, T = gates(cfg, cells, ctrl, majors, lg)
        T["config"] = cfg; allT.append(T)
        if win3:
            Gw, Tw = gates(cfg, cells, ctrl, majors, lg, setup=setups_[win3])
            G.update(winner_setup=win3, winner_composite=Gw["composite"], winner_composite_sig=Gw["composite_sig"],
                     winner_sig_major_median=Gw["sig_major_median"])
            Tw["config"] = cfg + "|" + win3; allT.append(Tw)
        summ.append(G)
        lg(f"{cfg}: " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in G.items()))
    # pulses through each estimator
    unit = P.unit_map()
    ccw = [f"{R.TB3}/store_ccw"] + ([f"{R.TB3}/store_ccw_ext"] if C.REGION == "SVI" else [])
    nodes = {}
    for (f, n, e), g in E[(E.n >= 50) & E.b1_cc.notna()].groupby(["fam", "netsta", "epoch"]):
        nodes[f"{f}|{n}|{e}"] = list(g.episode)              # measurement order (pair direction matters slightly)
    est = {"B0": (J.est_raw, J.store_paths_P()), "B1": (J.est_raw, ccw), "B2": (J.est_graph, ccw), "B2p": (J.est_graph, J.store_paths_P())}
    vals_cache = {}
    for cfg in CONFIGS:
        cells = CEL[cfg]
        base = "B0" if cfg in ("B0", "B0n50") else cfg
        if base not in vals_cache:
            cb = CEL[base] if base != "B0" else pd.concat([CEL["B0"], CEL["B0n50"]]).drop_duplicates(["series", "episode"])
            coef = {k: (s_, d_) for k, s_, d_ in zip(cb.series, cb.a_s, cb.a_d)}
            fn, stores = est[base]
            vals_cache[base] = J.run(set(cb.series), pulse_runs(coef, majors), fn, stores, a.workers, lg, nodes=nodes, unit=unit)
            lg(f"  pulses {base} measured {time.time() - t0:.0f}s")
        o, T = pulses(cfg, cells, vals_cache[base], majors, lg)
        pul.append(T)
        for G in summ:
            if G["config"] == cfg:
                G.update(o)
    S = pd.DataFrame(summ)
    b0 = S.set_index("config").composite_sig["B0"]
    S["adopt"] = [bool(getattr(r, "deep_pass", False) is True and getattr(r, "control_pass", False) is True and getattr(r, "loso_pass", False) is True
                       and r.composite_sig <= 1.1 * b0) if r.config in ("B1", "B2", "B2p") else False for r in S.itertuples()]
    adopted = next((c for c in ("B2", "B1", "B2p") if S.set_index("config").adopt.get(c, False)), None)
    S.to_csv(f"{R.TB3}/r3q3_summary.csv", index=False)
    pd.concat(pul, ignore_index=True).to_csv(f"{R.TB3}/r3q3_pulses.csv", index=False)
    pd.concat(allT, ignore_index=True).to_csv(f"{R.TB3}/r3q3_dD.csv", index=False)
    R.dump(dict(adopted=adopted, b0_composite_sig=float(b0), written=R.pt_now()), f"{R.TB3}/r3q3_verdict.json")
    lg(f"ADOPTED estimator for Q3: {adopted}; done {time.time() - t0:.0f}s")
    open(f"{R.R3}/R3Q3_DONE", "w").write(R.pt_now() + "\n")


if __name__ == "__main__":
    main()
