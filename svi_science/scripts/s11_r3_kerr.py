#!/usr/bin/env python
"""s11_r3_kerr -- S11 round 3, task 1 (PLAN §6.1): a kernel-robust deep term.
Setups NOM / AVG4 / AVG6 / LOKO_k (k in V6) / L-SEL / L-W / L-STA; waveform-level injections (P store, raw path) generated
with every kernel of V6 and with AVG6; six rules R1-R6; winner by the fixed order; Q1/Q2/Q3 re-reported under the winner
with the §3-§5 gates. dv/v in %, eps = -dv/v.
    S11_REGION=... python s11_r3_kerr.py [--workers 12] [--nperm 500]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time
import numpy as np
import pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_r3_common as R
import s11_r3_fit as LF
import s11_r3_inject as J
import s11_v2 as V

STEP_DATE, PEAK = "2015-01-01", 233.0
GEN = R.V6 + ["AVG6"]


def gen_setup(g):
    return ("avg", R.V6) if g == "AVG6" else g


def setups(lapse):
    S = {"NOM": "ell40", "AVG4": ("avg", R.V4), "AVG6": ("avg", R.V6)}
    for k in R.V6:
        S[f"LOKO_{k}"] = ("avg", [v for v in R.V6 if v != k])
    S["L-SEL"] = lapse["sel"]; S["L-W"] = ("w", lapse["w"]); S["L-STA"] = ("site", lapse["site"])
    return S


def outside(sname, g):
    """generating kernels that are out of sample for a setup (PLAN §6.1 'Pass')"""
    if sname.startswith("LOKO_"):
        return g == sname[5:]
    if sname.startswith("L-"):
        return g in R.V6
    if sname == "AVG4":
        return g in ("ell10", "two3_20")
    if sname == "NOM":
        return g in R.V6 and g != "ell40"
    if sname == "AVG6":
        return False
    return False


# ------------------------------------------------------------------------------ lapse choice (data only)
def lapse_choice(CL, lg):
    rows, site_chi = [], {}
    okall = np.ones(len(CL), bool)                   # identical cell set for every variant
    for v in R.V6:
        c = R.lapse_with_setup(CL, v); okall &= np.isfinite(c.a_d.to_numpy()) & np.isfinite(c.a_s.to_numpy())
    lg(f"lapse cells usable under all six variants: {int(okall.sum())} of {len(CL)}")
    for v in R.V6:
        c = R.lapse_with_setup(CL, v)
        c = c[okall].reset_index(drop=True)
        F = LF.smooth(c).solve()
        chi = float(F["wr2"].sum())
        rows.append(dict(variant=v, chi2_sum=chi, chi2_red=F["chi2"], n_cells=len(c)))
        site_chi[v] = pd.Series(F["wr2"]).groupby(c.site.to_numpy()).sum()
    T = pd.DataFrame(rows)
    T["dchi2"] = T.chi2_sum - T.chi2_sum.min()
    T["w"] = np.exp(-T.dchi2 / 2); T["w"] /= T.w.sum()
    sel = T.variant[T.dchi2.idxmin()]
    SC = pd.DataFrame(site_chi)
    ncell = CL.groupby("site").size()
    site = {}
    for s in SC.index:
        site[s] = SC.loc[s].idxmin() if ncell.get(s, 0) >= 30 else sel
    for s in R.C.SITES:
        site.setdefault(s, sel)
    lg("lapse choice: " + ", ".join(f"{r.variant} dchi2 {r.dchi2:.1f} w {r.w:.3f}" for r in T.itertuples()) + f"; L-SEL {sel}")
    lg(f"L-STA per site: {pd.Series(site).value_counts().to_dict()}")
    SC.assign(n_lapse_cells=ncell.reindex(SC.index), chosen=[site[s] for s in SC.index]).to_csv(f"{R.TB3}/kerr_lapse_sites.csv")
    T.to_csv(f"{R.TB3}/kerr_lapse_choice.csv", index=False)
    return dict(sel=sel, w=dict(zip(T.variant, T.w)), site=site)


def station_lapse(CL, cells, lg):
    """diagnostic: per-station kernel-free in-phase annual amplitude ratio 3.5-6 / 2-3.5 s vs each variant's predictions"""
    c0 = cells.copy(); c0["dm"] = V.demeaned(c0)
    pk = C.annual_fit(c0.episode.map(V.U.mid), c0.dm.to_numpy(), c0.sig.to_numpy())["peak_doy"]
    out = []
    for ns in sorted(CL.netsta.unique()):
        a = CL[(CL.netsta == ns) & (CL.win == "2035")]; b = CL[(CL.netsta == ns) & (CL.win == "3560")]
        if len(a) < 10 or len(b) < 10:
            continue
        A1, s1 = V.inphase(a, pk); A2, s2 = V.inphase(b, pk)
        r = dict(netsta=ns, site=C.SITE[ns], n2035=len(a), n3560=len(b), A2035=A1, A2035_sig=s1, A3560=A2, A3560_sig=s2,
                 ratio=A2 / A1 if A1 != 0 else np.nan, ratio_sig=abs(A2 / A1) * np.hypot(s1 / A1, s2 / A2) if A1 != 0 and A2 != 0 else np.nan)
        for v in R.V6:
            aa = R.lapse_with_setup(pd.concat([a, b]), v)
            r[f"deep_{v}"] = aa[aa.win == "3560"].a_d.mean() / aa[aa.win == "2035"].a_d.mean()
            r[f"shal_{v}"] = aa[aa.win == "3560"].a_s.mean() / aa[aa.win == "2035"].a_s.mean()
        out.append(r)
    T = pd.DataFrame(out); T.to_csv(f"{R.TB3}/kerr_station_lapse.csv", index=False)
    lg(f"station lapse diagnostic: {len(T)} stations, pooled peak day {pk:.0f}")
    return T, pk


# ------------------------------------------------------------------------------ injections
def make_runs_fn(gcoef, majors, unit):
    def fn(key, st, unit_):
        t = pd.to_datetime(pd.Series(st["days"]))
        post = (t >= STEP_DATE).to_numpy(); ann = np.cos(2 * np.pi * (t.dt.dayofyear.to_numpy() - PEAK) / 365.25)
        ep = np.array([unit_.get(d, None) for d in st["days"]], dtype=object)
        runs = {"nat": None}
        for g in GEN:
            if key not in gcoef[g]:
                continue
            a_s, a_d = gcoef[g][key]
            if not (np.isfinite(a_s) and np.isfinite(a_d)):
                continue
            runs[f"A|{g}"] = -(0.003 * a_d) * post
            runs[f"B|{g}"] = -(0.001 * a_s) * ann
            mall = np.isin(ep, majors)
            if mall.any():
                runs[f"Psall|{g}"] = -(0.001 * a_s) * mall
            for e in majors:
                m = ep == e
                if m.any():
                    runs[f"Pd|{g}|{e}"] = -(0.003 * a_d) * m
                    runs[f"Ps|{g}|{e}"] = -(0.001 * a_s) * m
        return runs
    return fn


def evaluate(cells, vals, S, majors, lg):
    """six rules for every setup x generating kernel"""
    dn = J.cell_vector(vals, cells, "nat")
    ok0 = np.isfinite(dn)
    lg(f"natural-from-store vs measured: max |diff| {np.nanmax(np.abs(dn - cells.d.to_numpy())):.1e} %, missing {int((~ok0).sum())}")
    runs = {}

    def dv(run):
        if run not in runs:
            runs[run] = J.cell_vector(vals, cells, run)
        return runs[run]
    rows = []
    for sname, setup in S.items():
        cs = R.with_setup(cells, setup)
        good = ok0 & np.isfinite(cs.a_d.to_numpy()) & np.isfinite(cs.a_s.to_numpy())
        cs = cs[good].reset_index(drop=True); dn_ = dn[good]
        Ls, Lt = LF.smooth(cs), LF.temporal(cs)
        F0 = Lt.delta(np.zeros(len(cs)))
        gcomp = LF.comp_vec(F0, majors)
        gmaj = {e: LF.dvec(F0, e) for e in majors}
        gmaj = {e: g_ for e, g_ in gmaj.items() if g_ is not None}
        for g in GEN:
            r = dict(setup=sname, setup_def=R.setup_name(setup), gen=g, outside=outside(sname, g))

            def diff(run):
                di = dv(run)[good]
                m = np.isfinite(di)
                x = np.where(m, di - dn_, 0.0)             # cells railed under injection: response treated as 0 (counted)
                return x, int((~m).sum())
            xa, na = diff(f"A|{g}"); Fd = Ls.delta(xa)
            r["R1_Q2rec"] = float(V.q2vec(Fd) @ Fd["th"]) / 0.3; r["R2_Q1leak"] = V.amp(Fd)[0]
            xb, nb = diff(f"B|{g}"); Fb = Ls.delta(xb); r["R3_Bleak"] = V.amp(Fb)[0] / 0.1
            r["R6_B_Q3comp"] = abs(float(gcomp @ Lt.delta(xb)["th"])) / 0.1 if gcomp is not None else np.nan
            rec, lk = [], []
            for e, g_e in gmaj.items():
                x, _ = diff(f"Pd|{g}|{e}"); rec.append(float(g_e @ (Lt.Op @ x)) / 0.3)
                x, _ = diff(f"Ps|{g}|{e}"); lk.append(abs(float(g_e @ (Lt.Op @ x))) / 0.1)
            r["R4_Pdrec"] = float(np.median(rec)) if rec else np.nan
            r["R4_lo"], r["R4_hi"] = (float(np.min(rec)), float(np.max(rec))) if rec else (np.nan, np.nan)
            r["R5_Psleak"] = float(np.median(lk)) if lk else np.nan
            xs, _ = diff(f"Psall|{g}")
            r["Psall_comp_leak"] = float(gcomp @ (Lt.Op @ xs)) / 0.1 if gcomp is not None else np.nan
            r["n_railed_A"], r["n_railed_B"] = na, nb
            r["p1"] = bool(0.7 <= r["R1_Q2rec"] <= 1.3); r["p2"] = bool(r["R2_Q1leak"] < 0.03); r["p3"] = bool(r["R3_Bleak"] < 0.25)
            r["p4"] = bool(0.8 <= r["R4_Pdrec"] <= 1.1); r["p5"] = bool(r["R5_Psleak"] <= 0.35); r["p6"] = bool(r["R6_B_Q3comp"] <= 0.35)
            rows.append(r)
        lg(f"  setup {sname}: " + "; ".join(f"{x['gen']} R1 {x['R1_Q2rec']:.2f} R2 {x['R2_Q1leak']:.3f} R3 {x['R3_Bleak']:.2f} "
                                             f"R4 {x['R4_Pdrec']:.2f} R5 {x['R5_Psleak']:.2f} R6 {x['R6_B_Q3comp']:.2f}"
                                             for x in rows if x["setup"] == sname))
    return pd.DataFrame(rows)


RULES = [("p1", "R1_Q2rec", 0.10), ("p2", "R2_Q1leak", 0.01), ("p3", "R3_Bleak", 0.10), ("p4", "R4_Pdrec", 0.10),
         ("p5", "R5_Psleak", 0.10), ("p6", "R6_B_Q3comp", 0.10)]


def verdict(T, lg):
    """PLAN §6.1 + Amendment 1: a rule that already fails in-sample (NOM generated with ell40) is not a kernel error; it is then
    judged on |R_g - R_in| <= tol (0.10; R2 0.01 %). 'strict' = every rule absolute."""
    ins = T[(T.setup == "NOM") & (T.gen == "ell40")].iloc[0]
    out = {}
    for s, g in T.groupby("setup"):
        o = g[g.outside]
        if not len(o):
            out[s] = dict(q12=None, q3=None, n_out=0); continue
        crit = {}
        for pf, col, tol in RULES:
            crit[pf] = o[pf] if bool(ins[pf]) else (o[col] - ins[col]).abs() <= tol
        q12 = bool((crit["p1"] & crit["p2"] & crit["p3"]).all()); q3 = bool((crit["p4"] & crit["p5"] & crit["p6"]).all())
        out[s] = dict(q12=q12, q3=q3, n_out=len(o), q12_strict=bool((o.p1 & o.p2 & o.p3).all()), q3_strict=bool((o.p4 & o.p5 & o.p6).all()),
                      failing=sorted({f"{pf}:{gg}" for pf, _, _ in RULES for gg, ok in zip(o.gen, crit[pf]) if not ok}))
    loko = [f"LOKO_{k}" for k in R.V6]
    order = ["L-W", "L-STA", "L-SEL", "AVG4"]
    win = {}
    for part in ("q12", "q3"):
        win[part] = "AVG6" if all(out[s][part] for s in loko) else next((s for s in order if out[s][part]), None)
    insd = {col: float(ins[col]) for _, col, _ in RULES}; insd.update({f"{pf}_in_ok": bool(ins[pf]) for pf, _, _ in RULES})
    lg(f"verdict per setup: {out}; in-sample (NOM, ell40): {insd}; WINNER Q1/Q2 {win['q12']}, Q3 {win['q3']}")
    return out, win, insd


# ------------------------------------------------------------------------------ re-report under a setup
def rereport(sname, setup, cells, cctrl, CL, majors, nperm, lg, label):
    cs = R.with_setup(cells, setup); cs = cs[np.isfinite(cs.a_d) & np.isfinite(cs.a_s)].reset_index(drop=True)
    G = []
    Ls = LF.smooth(cs); F = Ls.solve(); q = LF.Qsmooth(F)
    G.append(dict(item="Q1Q2", **q))
    for f_ in (2.0, 0.5):
        qq = LF.Qsmooth(LF.smooth(cs, prior=f_).solve())
        G.append(dict(item=f"prior_x{f_:g}", dQ1_over_sig=(qq["Q1"] - q["Q1"]) / q["Q1_sig"], dQ2_over_sig=(qq["Q2"] - q["Q2"]) / q["Q2_sig"],
                      passed=bool(abs((qq["Q1"] - q["Q1"]) / q["Q1_sig"]) < 0.1 and abs((qq["Q2"] - q["Q2"]) / q["Q2_sig"]) < 0.1)))
    pk0 = q["Q1_peak"]; dpk = lambda p: abs((p - pk0 + 182.6) % 365.25 - 182.6)
    ok1 = True
    for s in sorted(cs.site.unique()):
        sub = cs[cs.site != s].reset_index(drop=True)
        if len(sub) < 20:
            continue
        qs = LF.Qsmooth(LF.smooth(sub).solve())
        if not (qs["Q1_sig"] > 0):
            continue
        w_ = abs(qs["Q1"] - q["Q1"]) < 2 * qs["Q1_sig"]; p_ = dpk(qs["Q1_peak"]) <= 60 if qs["Q1"] / qs["Q1_sig"] >= 1 else True
        ok1 &= bool(w_ and p_)
        G.append(dict(item=f"G1_without_{s}", Q1=qs["Q1"], Q1_sig=qs["Q1_sig"], Q1_peak=qs["Q1_peak"], Q2=qs["Q2"], Q2_sig=qs["Q2_sig"], passed=bool(w_ and p_)))
    G.append(dict(item="G1_all", passed=bool(ok1)))
    qo = LF.Qsmooth(LF.smooth(cs[cs.fam % 2 == 1].reset_index(drop=True)).solve())
    qe = LF.Qsmooth(LF.smooth(cs[cs.fam % 2 == 0].reset_index(drop=True)).solve())
    w2 = abs(qo["Q1"] - qe["Q1"]) < 2 * np.hypot(qo["Q1_sig"], qe["Q1_sig"])
    pko = abs((qo["Q1_peak"] - qe["Q1_peak"] + 182.6) % 365.25 - 182.6) <= 60 if min(qo["Q1"] / max(qo["Q1_sig"], 1e-12), qe["Q1"] / max(qe["Q1_sig"], 1e-12)) >= 1 else True
    G.append(dict(item="G2_odd", **qo)); G.append(dict(item="G2_even", **qe, passed=bool(w2 and pko)))
    cc = R.with_setup(cctrl, setup); cc = cc[np.isfinite(cc.a_d) & np.isfinite(cc.a_s)].reset_index(drop=True)
    qc = LF.Qsmooth(LF.smooth(cc).solve())
    G.append(dict(item="G3_scrambled_control", **qc, passed=bool(qc["Q1"] < 2 * qc["Q1_sig"])))
    rng = np.random.default_rng(C.SEED); eu = cs.episode.unique(); mids = V.U.mid.reindex(eu).to_numpy()
    perm = np.array([LF.Qsmooth(LF.smooth(cs, mid=cs.episode.map(dict(zip(eu, rng.permutation(mids))))).solve())["Q1"] for _ in range(nperm)])
    G.append(dict(item="G3col_date_permutation", value=float((1 + np.sum(perm >= q["Q1"])) / (1 + len(perm))), perm_q95=float(np.quantile(perm, 0.95))))
    # G4 station pattern (kernel-free amplitudes vs the setup's mean a_d / a_s)
    c0 = cs.copy(); c0["dm"] = V.demeaned(c0)
    pkp = C.annual_fit(c0.episode.map(V.U.mid), c0.dm.to_numpy(), c0.sig.to_numpy())["peak_doy"]
    st = []
    for ns in sorted(cs.netsta.unique()):
        qq = cs[cs.netsta == ns]
        if len(qq) < 10:
            continue
        A_, s_ = V.inphase(qq, pkp); w = 1 / qq.sig ** 2
        st.append(dict(netsta=ns, n=len(qq), A=A_, A_sig=s_, a_d=float(np.average(qq.a_d, weights=w)), a_s=float(np.average(qq.a_s, weights=w))))
    ST = pd.DataFrame(st); chi = {}
    for col in ("a_d", "a_s"):
        x, y, s_ = ST[col].to_numpy(), ST.A.to_numpy(), ST.A_sig.to_numpy()
        k = np.sum(x * y / s_ ** 2) / np.sum(x ** 2 / s_ ** 2); chi[col] = float(np.sum((y - k * x) ** 2 / s_ ** 2))
        G.append(dict(item=f"G4_chi2_prop_{col}", value=chi[col], dof=len(x) - 1, p=float(stats.chi2.sf(chi[col], len(x) - 1))))
    G.append(dict(item="G4_deep_pattern_supported", value=chi["a_s"] - chi["a_d"],
                  passed=bool(stats.chi2.sf(chi["a_d"], len(ST) - 1) > 0.05 and chi["a_s"] - chi["a_d"] > 4)))
    # G5 lapse ratio vs the setup's window coefficients
    cA, cB = CL[CL.win == "2035"], CL[CL.win == "3560"]
    A1, s1 = V.inphase(cA, pkp); A2, s2 = V.inphase(cB, pkp)
    ratio = A2 / A1; sr = abs(ratio) * np.hypot(s1 / A1, s2 / A2)
    CLs = R.lapse_with_setup(CL, setup)
    pdp = CLs[CLs.win == "3560"].a_d.mean() / CLs[CLs.win == "2035"].a_d.mean()
    psp = CLs[CLs.win == "3560"].a_s.mean() / CLs[CLs.win == "2035"].a_s.mean()
    G.append(dict(item="G5_lapse_ratio", value=float(ratio), sig=float(sr), pred_deep=float(pdp), pred_shallow=float(psp),
                  z_deep=float((ratio - pdp) / sr), z_shallow=float((ratio - psp) / sr), passed=bool(abs((ratio - psp) / sr) - abs((ratio - pdp) / sr) >= 1)))
    # Q3 under the setup
    Lt = LF.temporal(cs); Ft = Lt.solve(); T, comp = LF.q3_table(Ft, majors)
    T.to_csv(f"{R.TB3}/kerr_{label}_q3.csv", index=False)
    mj = T[T.major]
    G.append(dict(item="Q3", composite=comp["composite"], composite_sig=comp["composite_sig"], n_majors=comp["n_majors"],
                  sig_major_median=float(mj.sig.median()), n_major_sig_le_0p1=int((mj.sig <= 0.1).sum()),
                  chi2_majors=float(np.mean((mj.dD / mj.sig) ** 2)), n_beyond_2sig=int((np.abs(mj.dD / mj.sig) > 2).sum())))
    Tm = T.set_index("episode")
    for f_ in (2.0, 0.5):
        Tf, cf = LF.q3_table(LF.temporal(cs, prior=f_).solve(), majors); Tf = Tf.set_index("episode")
        e_ = [e for e in mj.episode if e in Tf.index]
        mv = float(np.max(np.abs((Tf.dD[e_] - Tm.dD[e_]) / Tm.sig[e_])))
        G.append(dict(item=f"Q3_prior_x{f_:g}", value=mv, composite=abs(cf["composite"] - comp["composite"]) / comp["composite_sig"],
                      passed=bool(mv < 0.1 and abs(cf["composite"] - comp["composite"]) / comp["composite_sig"] < 0.1)))
    Tc, ccomp = LF.q3_table(LF.temporal(cc).solve(), majors); Tc = Tc.set_index("episode")
    meas = [e for e in mj.episode if e in Tc.index and Tc.n_cells[e] >= 5]
    frac = float(np.mean([abs(Tc.dD[e]) < 2 * Tc.sig[e] for e in meas])) if meas else np.nan
    G.append(dict(item="Q3_scrambled_control", value=frac, n=len(meas), composite=ccomp["composite"] if ccomp else np.nan,
                  composite_sig=ccomp["composite_sig"] if ccomp else np.nan,
                  passed=bool(frac >= 0.9 and (ccomp is None or abs(ccomp["composite"]) < 2 * ccomp["composite_sig"])) if meas else None))
    worst = 1.0
    for s in sorted(cs.site.unique()):
        if (cs.site != s).sum() < 20:
            continue
        Ts, csx = LF.q3_table(LF.temporal(cs[cs.site != s].reset_index(drop=True)).solve(), majors); Ts = Ts.set_index("episode")
        ee = [e for e in mj.episode if e in Ts.index]
        fr = float(np.mean([abs(Ts.dD[e] - Tm.dD[e]) < 2 * Ts.sig[e] for e in ee])) if ee else np.nan
        worst = min(worst, fr) if np.isfinite(fr) else worst
        G.append(dict(item=f"Q3_loso_{s}", value=fr, composite=csx["composite"] if csx else np.nan, composite_sig=csx["composite_sig"] if csx else np.nan,
                      passed=bool(fr >= 0.9) if ee else None))
    G.append(dict(item="Q3_loso_all", value=worst, passed=bool(worst >= 0.9)))
    To, co = LF.q3_table(LF.temporal(cs[cs.fam % 2 == 1].reset_index(drop=True)).solve(), majors)
    Te, ce = LF.q3_table(LF.temporal(cs[cs.fam % 2 == 0].reset_index(drop=True)).solve(), majors)
    if co and ce:
        J_ = To.set_index("episode").join(Te.set_index("episode"), lsuffix="_o", rsuffix="_e", how="inner")
        J_ = J_[J_.major_o]
        fr = float(np.mean(np.abs(J_.dD_o - J_.dD_e) < 2 * np.hypot(J_.sig_o, J_.sig_e))) if len(J_) else np.nan
        G.append(dict(item="Q3_family_split", value=fr, composite=co["composite"], composite_sig=co["composite_sig"],
                      composite_even=ce["composite"], composite_even_sig=ce["composite_sig"],
                      passed=bool(fr >= 0.9 and abs(co["composite"] - ce["composite"]) < 2 * np.hypot(co["composite_sig"], ce["composite_sig"]))))
    GT = pd.DataFrame(G); GT.insert(0, "setup", sname)
    GT.to_csv(f"{R.TB3}/kerr_{label}_gates.csv", index=False)
    for r in G:
        lg(f"  [{label}] " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    return GT


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--nperm", type=int, default=500)
    ap.add_argument("--limit", type=int, default=0, help="test: only this many series")
    a = ap.parse_args()
    lg = R.log("kerr"); t0 = time.time()
    majors = R.majors()
    cells = R.head_cells()
    ref = pd.read_csv(f"{C.TB}/q3_cells_c_115fam_n20.csv" if C.REGION == "SVI" else f"{C.TB}/v2_cells.csv")
    lg(f"region {C.REGION}: {len(cells)} headline cells (stored file {len(ref)}; max |d diff| "
       f"{np.max(np.abs(np.sort(cells.d.to_numpy()) - np.sort(ref.d.to_numpy()))) if len(ref) == len(cells) else 'n/a'})")
    if a.limit:                                     # test run: outputs to tables/test/
        have = set(R.coefs("ell10", C.HEAD).index)    # loads the kernel table before the redirect
        R.TB3 = f"{R.TB3}/test"; os.makedirs(R.TB3, exist_ok=True)
        keep = [k for k in dict.fromkeys(cells.series) if k in have][:a.limit]; cells = cells[cells.series.isin(keep)].reset_index(drop=True)
    LF.check_against_original(cells, lg)
    cctrl = R.head_cells(control=True)
    CL = R.lapse_cells()
    lg(f"control cells {len(cctrl)}, lapse cells 2035 {int((CL.win == '2035').sum())} 3560 {int((CL.win == '3560').sum())}")
    lap = lapse_choice(CL, lg)
    station_lapse(CL, cells, lg)
    S = setups(lap)
    # values of every setup on the data (Q1, Q2, Q3 composite)
    vrows = []
    for sname, setup in S.items():
        cs = R.with_setup(cells, setup); cs = cs[np.isfinite(cs.a_d) & np.isfinite(cs.a_s)].reset_index(drop=True)
        q = LF.Qsmooth(LF.smooth(cs).solve())
        T, comp = LF.q3_table(LF.temporal(cs).solve(), majors); mj = T[T.major]
        vrows.append(dict(setup=sname, setup_def=R.setup_name(setup), **q, Q3_composite=comp["composite"], Q3_composite_sig=comp["composite_sig"],
                          Q3_sig_major_median=float(mj.sig.median()), a_s_med=float(cs.a_s.median()), a_d_med=float(cs.a_d.median())))
    pd.DataFrame(vrows).to_csv(f"{R.TB3}/kerr_values.csv", index=False)
    for r in vrows:
        lg("value " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    # generating coefficients per series key (2-5 s)
    sk = cells.drop_duplicates("series")[["fam", "netsta", "epoch", "series", "site"]].reset_index(drop=True)
    gcoef = {}
    for g in GEN:
        s_, d_ = R.setup_coefs(gen_setup(g), C.HEAD, sk)
        gcoef[g] = {k: (x, y) for k, x, y in zip(sk.series, s_, d_)}
    unit = __import__("s11_paths").unit_map()
    vals = J.run(set(sk.series), make_runs_fn(gcoef, majors, unit), J.est_raw, J.store_paths_P(), a.workers, lg, unit=unit)
    lg(f"injections done {time.time() - t0:.0f}s ({len(vals)} series)")
    T = evaluate(cells, vals, S, majors, lg)
    T.to_csv(f"{R.TB3}/kerr_rules.csv", index=False)
    out, win, ins = verdict(T, lg)
    R.dump(dict(per_setup=out, winner=win, in_sample=ins, setups={k: R.setup_name(v) for k, v in S.items()}, lapse=lap,
                written=R.pt_now()), f"{R.TB3}/kerr_verdict.json")
    # re-report under the winner(s); AVG6 and NOM always as columns
    done = set()
    for label, sname in (("NOM", "NOM"), ("AVG6", "AVG6"), ("win_q12", win["q12"]), ("win_q3", win["q3"])):
        if sname is None or sname in done:
            continue
        rereport(sname, S[sname], cells, cctrl, CL, majors, a.nperm if label.startswith("win") or sname == "AVG6" else 50, lg, sname)
        done.add(sname)
    if C.REGION == "SVI":                          # SVI v2 column (76 families, n >= 50) for Q1/Q2
        c50 = R.I.prepare("raw", C.HEAD, "", "ell40", nmin=50, ccmin=R.CCMIN)
        rr = []
        for sname in ("NOM", "AVG6", win["q12"]):
            if sname is None:
                continue
            cs = R.with_setup(c50, S[sname]); cs = cs[np.isfinite(cs.a_d)].reset_index(drop=True)
            rr.append(dict(setup=sname, cells="v2_76fam_n50", **LF.Qsmooth(LF.smooth(cs).solve())))
        pd.DataFrame(rr).to_csv(f"{R.TB3}/kerr_svi_v2_column.csv", index=False)
        for r in rr:
            lg("SVI v2 column " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    lg(f"kerr done {time.time() - t0:.0f}s")
    if not a.limit:
        open(f"{R.R3}/KERR_DONE", "w").write(R.pt_now() + "\n")


if __name__ == "__main__":
    main()
