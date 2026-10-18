#!/usr/bin/env python
"""s11_r6_model -- S11 round 6, task 1 (PLAN.md §10.1): why WA L-SEL deep-pulse recovery T1 is 0.75 under E4, and the fix.
  --stage forms : per region, forms E3 / E4 (round-5 reproductions), E4y (D ⊥ c1,s1), E4l (⊥ lin,c1,s1), E4h (⊥ without the
                  end hinges), E4y+Scs2 under NOM (ell40) and L-SEL (two3_20): natural Q3, every gate, T1-T5, per-major
                  cause analysis (constraint leverage, model-consistent pulse, waveform pulse), Q3 kernel band.
  --stage adopt : both regions -> the adopted form E4' (first of E4y, E4l, E4h, E4y+Scs2 passing in all 4 combinations)
                  -> round6/adopted.json
    S11_REGION=... python s11_r6_model.py --stage forms ;  python s11_r6_model.py --stage adopt
Reuses round-5 cells H, the round-4 waveform-level injections and the round-5 model class (no new measurement)."""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_v2 as V
import s11_r3_fit as LF
import s11_r4_common as R4
import s11_r5_common as R5

R6 = f"{R5.BASE}/round6"
OUT6 = f"{R6}/{R5.SUB}"
TB6, LOG6, FIG6 = f"{OUT6}/tables", f"{OUT6}/logs", f"{R6}/figures"
for _d in (TB6, LOG6, FIG6):
    os.makedirs(_d, exist_ok=True)
SETUPS = ("ell40", "two3_20")
PROJ_ALL = list(R5.PROJ)                                   # lin h2007 h2011 h2015 h2019 h2023 c1 s1 (round-5 Dp)
PROJS = {"Dp": PROJ_ALL, "Dy": ["c1", "s1"], "Dl": ["lin", "c1", "s1"],
         "Dh": ["lin", "h2011", "h2015", "h2019", "c1", "s1"], "D0": []}
FORMS6 = {   # name -> (round-5 kwargs, projection set)
    "E3": (dict(unit="one", Dproj=False), "D0"),
    "E4": (dict(unit="one", Dproj=True), "Dp"),
    "E4y": (dict(unit="one", Dproj=True), "Dy"),
    "E4l": (dict(unit="one", Dproj=True), "Dl"),
    "E4h": (dict(unit="one", Dproj=True), "Dh"),
    "E4y+Scs2": (dict(unit="one", Dproj=True, Sc=("spline", 2)), "Dy"),
}
ORDER = ["E4y", "E4l", "E4h", "E4y+Scs2"]                  # PLAN §10.1 adoption order
# §10 Amendment 1 (POST-RESULT for task 1): the shared slow directions removed on both sides (unit basis full / small)
FORMS6["E4lu"] = (dict(unit="one", Dproj=True), "Dl")
FORMS6["E4yu"] = (dict(unit="one", Dproj=True), "Dy")
UTERMS = {"E4lu": (["c1", "s1", "c2", "s2", "lin"], ["c1", "s1", "lin"]), "E4yu": (["c1", "s1", "c2", "s2"], ["c1", "s1"])}
MAIN = ["E3", "E4", "E4y", "E4l", "E4h", "E4y+Scs2"]
AMEND = ["E4lu", "E4yu"]
ORDER_A = ["E4lu", "E4yu"]
LIMIT = [0]


def log(name):
    return C.Log(f"{LOG6}/{name}.log")


def factory6(kw, proj, uterm=None):
    """round-5 Model5 with the D constraint set `proj` (R5.PROJ is read inside Model5 at construction) and optionally the
    unit basis (full, small) replaced (R5.STERMS / R5.SMALL, read by units_of and Model5 at construction)"""
    pl = PROJS[proj]

    def make(cells, **o):
        old = (R5.PROJ, R5.STERMS, R5.SMALL)
        R5.PROJ = pl
        if uterm is not None:
            R5.STERMS, R5.SMALL = list(uterm[0]), list(uterm[1])
        try:
            return R5.Model5(cells, **{**kw, **o})
        finally:
            R5.PROJ, R5.STERMS, R5.SMALL = old
    return make


def ctrl_b1(cs):
    import s11_r4_combined as CB
    return CB._ctrl_b1(cs, "ell40")


def leverage_pred(eps_, proj, e):
    """projection-only ΔD recovery of a unit pulse at e: the D constraint basis [1, proj] at the episode mids"""
    U = R5.U
    em = U.mid.reindex(eps_); te, de = V.tdec(em.to_numpy())
    cols = [np.ones(len(eps_))]
    if proj:
        bb = V.basis(te, de, PROJS[proj] if isinstance(proj, str) else proj)
        cols += [bb[t] for t in (PROJS[proj] if isinstance(proj, str) else proj) if np.any(bb[t] != 0)]
    B = np.column_stack(cols)
    P = B @ np.linalg.pinv(B)
    i = eps_.index(e)
    d = np.zeros(len(eps_)); d[i] = 1.0
    v = d - P @ d
    nb = [j for j, x in enumerate(eps_) if x != e and abs((U.mid[x] - U.mid[e]).days) <= LF.WIN3]
    pred = v[i] - (np.mean(v[nb]) if nb else np.nan)
    return float(P[i, i]), float(pred)


def per_major(H, inj, setup, majors, fac, proj, label):
    cs, good = R5.setup_cells(H, setup)
    good = good & np.isfinite(inj["nat"])
    cs = cs[good].reset_index(drop=True)
    L = fac(cs); F0 = L.delta(np.zeros(len(cs)))
    dn = inj["nat"][good]
    ep = cs.episode.to_numpy(); ad = cs.a_d.to_numpy(float)
    rows = []
    for e in majors:
        if e not in L.eps:
            continue
        r = dict(form=label, setup=setup, episode=e, n_cells=int((ep == e).sum()))
        h, pred = leverage_pred(list(L.eps), proj if proj != "D0" else None, e)
        r.update(h_ee=h, R_proj_pred=pred)
        x = np.where(ep == e, 0.3 * ad, 0.0)
        rr, _ = R5.resp(L, F0, x, [e]); r["R_struct"] = rr.get(e, np.nan) / 0.3
        k = f"Pd|{setup}|{e}"
        if k in inj:
            v = inj[k][good]; X = np.where(np.isfinite(v), v - dn, 0.0)
            rr, _ = R5.resp(L, F0, X, [e]); r["R_wave"] = rr.get(e, np.nan) / 0.3
        else:
            r["R_wave"] = np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def run_form(H, inj, majors, name, setup, lg):
    kw, proj = FORMS6[name]
    fac = factory6(kw, proj, UTERMS.get(name))
    out = f"{TB6}/forms_{name}_{setup}"
    G, T, sw, F, L, cs = R5.gates(H, setup, majors, fac, name, lg, exact=True, ctrl_b1_fn=ctrl_b1, out=out)
    row = dict(region=C.REGION, form=name, proj=proj, **G)
    if "2019-11" in set(T.episode):
        x = T.set_index("episode").loc["2019-11"]; row.update(dD_2019_11=float(x.dD), sig_2019_11=float(x.sig))
    if "2023-03" in set(T.episode):
        x = T.set_index("episode").loc["2023-03"]; row.update(dD_2023_03=float(x.dD), sig_2023_03=float(x.sig))
    TR, O = R5.evaluate(H, inj, setup, majors, fac, sw, name, lg)
    v = R5.rules(TR, O, setup)
    row.update({f"rule_{k}": (json.dumps(x) if isinstance(x, list) else x) for k, x in v.items()})
    row["G_ctrl"] = row.get("ctrl_g_pass"); row["G_loso"] = row.get("G_loso_pass"); row["G_fam"] = row.get("G_fam_pass")
    row["G_prior"] = row.get("G_prior_pass")
    for k in ("T1", "T2", "T3", "T4", "T5"):
        row[k] = row.get(f"rule_{k}")
    items = ("T1", "T2", "T3", "T4", "T5", "G_prior", "G_ctrl", "G_loso", "G_fam", "acceptance")
    row["fail_items"] = ",".join(k for k in items if not bool(row.get(k)))
    row["pass_all"] = row["fail_items"] == ""
    row["n_pass"] = int(sum(bool(row.get(k)) for k in items))
    PM = per_major(H, inj, setup, majors, fac, proj, name)
    return row, TR.assign(region=C.REGION, form=name), O.assign(region=C.REGION, form=name), T.assign(form=name, setup=setup), PM


def stage_forms(lg, names=None, pre="r6"):
    names = MAIN if names is None else names
    majors = R4.majors()
    H, inj = R5.load_H(lg, LIMIT[0])
    rows, TRs, Os, Ts, PMs = [], [], [], [], []
    for name in names:
        for setup in SETUPS:
            t0 = time.time()
            row, TR, O, T, PM = run_form(H, inj, majors, name, setup, lg)
            rows.append(row); TRs.append(TR); Os.append(O); Ts.append(T); PMs.append(PM)
            lg(f"FORM {name}|{setup}: σ/major {row['sig_major_median']:.4f} ({row['n_major_sig_le_0p1']}/{row['n_majors']} ≤ 0.1), "
               f"composite {row['composite']:+.4f} ± {row['composite_sig']:.4f}, T1 {row['rule_T1_val']:.3f}, leak {row['rule_T2_val']:.3f}, "
               f"T4 {row['rule_T4_all_comp']:.3f}/{row['rule_T4_all_major']:.3f}, prior {row['prior_max']:.3f}σ, "
               f"fail [{row['fail_items']}] ({time.time() - t0:.0f}s)")
            lg(f"  per-major T1 (wave) {PM.R_wave.min():.2f}-{PM.R_wave.max():.2f} median {PM.R_wave.median():.3f}; "
               f"struct median {PM.R_struct.median():.3f}; proj-pred median {PM.R_proj_pred.median():.3f}; "
               f"ρ(h_ee, R_wave) {PM[['h_ee', 'R_wave']].corr('spearman').iloc[0, 1]:.2f}")
    D = pd.DataFrame(rows)
    PM = pd.concat(PMs, ignore_index=True)
    # reproduction of round-5 E3 / E4
    if not LIMIT[0] and pre == "r6":
        r5 = pd.read_csv(f"{R5.TB5}/r5_forms.csv")
        chk = []
        for name in ("E3", "E4"):
            for setup in SETUPS:
                a = D[(D.form == name) & (D.setup == setup)].iloc[0]; b = r5[(r5.form == name) & (r5.setup == setup)].iloc[0]
                chk.append(dict(form=name, setup=setup, d_sig_major=float(a.sig_major_median - b.sig_major_median),
                                d_composite=float(a.composite - b.composite), d_composite_sig=float(a.composite_sig - b.composite_sig),
                                d_T1=float(a.rule_T1_val - b.rule_T1_val), d_prior=float(a.prior_max - b.prior_max)))
        CK = pd.DataFrame(chk); CK.to_csv(f"{TB6}/forms_reproduction_check.csv", index=False)
        mx = float(np.nanmax(np.abs(CK.drop(columns=["form", "setup"]).to_numpy())))
        lg(f"reproduction of round-5 E3 / E4 (σ/major, composite, σ, T1, prior): max |Δ| {mx:.1e} -> {'OK' if mx <= 1e-6 else 'MISMATCH'}")
        D["repro_max_abs"] = mx
    # cause summary: E4 - E3 shortfall split
    cz = []
    for setup in SETUPS:
        a = PM[(PM.setup == setup)].pivot(index="episode", columns="form", values=["R_wave", "R_struct", "R_proj_pred", "h_ee"])
        for f in names:
            if ("R_wave", f) not in a.columns:
                continue
            cz.append(dict(region=C.REGION, setup=setup, form=f, R_wave_med=float(a["R_wave"][f].median()),
                           R_struct_med=float(a["R_struct"][f].median()), R_proj_pred_med=float(a["R_proj_pred"][f].median()),
                           h_ee_med=float(a["h_ee"][f].median()), h_ee_max=float(a["h_ee"][f].max()),
                           rho_h_Rwave=float(pd.concat([a["h_ee"][f], a["R_wave"][f]], axis=1).corr("spearman").iloc[0, 1]),
                           **({} if ("R_wave", "E3") not in a.columns else dict(
                               shortfall_vs_E3_wave=float((a["R_wave"]["E3"] - a["R_wave"][f]).median()),
                               shortfall_vs_E3_struct=float((a["R_struct"]["E3"] - a["R_struct"][f]).median()),
                               shortfall_vs_E3_proj=float((a["R_proj_pred"]["E3"] - a["R_proj_pred"][f]).median())))))
    CZ = pd.DataFrame(cz); CZ.to_csv(f"{TB6}/{pre}_cause.csv", index=False)
    for r in CZ.itertuples():
        lg(f"CAUSE {r.form}|{r.setup}: R_wave {r.R_wave_med:.3f}, R_struct {r.R_struct_med:.3f}, proj-only {r.R_proj_pred_med:.3f}, "
           f"h_ee med/max {r.h_ee_med:.3f}/{r.h_ee_max:.3f}; ρ(h, R_wave) {r.rho_h_Rwave:.2f}")
    # Q3 band over the fit kernels (task-5 definition) for every candidate form
    band = []
    for name in names:
        kw, proj = FORMS6[name]; fac = factory6(kw, proj, UTERMS.get(name))
        for v in R5.GEN:
            cs, good = R5.setup_cells(H, v); cs = cs[good].reset_index(drop=True)
            F = fac(cs).solve(); T, comp = R5.q3(F, majors); mj = T[T.major]
            band.append(dict(region=C.REGION, form=name, fit_setup=v, plausible=v in R5.PLAUS[C.REGION], composite=comp["composite"],
                             composite_sig=comp["composite_sig"], sig_major_median=float(mj.sig.median()),
                             n_major_sig_le_0p1=int((mj.sig <= 0.1).sum()), n_majors_ge3sig=int((np.abs(mj.dD / mj.sig) >= 3).sum())))
    B = pd.DataFrame(band); B.to_csv(f"{TB6}/{pre}_band_q3.csv", index=False)
    TRa = pd.concat(TRs, ignore_index=True)
    for name in names:
        b = B[(B.form == name) & B.plausible]
        bound = float(np.max(np.abs(b.composite) + 2 * b.composite_sig))
        t1 = TRa[(TRa.form == name) & (TRa.setup == "two3_20") & TRa.gen.isin(R5.PLAUS[C.REGION])].T1_Pd_rec
        D.loc[D.form == name, "band_bound2sig"] = bound
        D.loc[D.form == name, "band_T1_min"] = float(t1.min()); D.loc[D.form == name, "band_T1_max"] = float(t1.max())
        D.loc[D.form == name, "band_within2sig"] = bool((np.abs(b.composite) < 2 * b.composite_sig).all())
        lg(f"BAND {name}: composite over P {b.composite.min():+.4f} … {b.composite.max():+.4f}; 2σ bound {bound:.3f} %; "
           f"T1 (L-SEL fit) over P {t1.min():.2f}–{t1.max():.2f} -> size band {bound / t1.max():.3f}–{bound / t1.min():.3f} %")
    D.to_csv(f"{TB6}/{pre}_forms.csv", index=False)
    TRa.to_csv(f"{TB6}/{pre}_forms_rules.csv", index=False)
    pd.concat(Os, ignore_index=True).to_csv(f"{TB6}/{pre}_forms_onesite.csv", index=False)
    pd.concat(Ts, ignore_index=True).to_csv(f"{TB6}/{pre}_forms_dD.csv", index=False)
    PM.to_csv(f"{TB6}/{pre}_t1_permajor.csv", index=False)


def stage_amend(lg):
    stage_forms(lg, AMEND, "r6a")


def stage_adopt_amend(lg):
    D = pd.concat([pd.read_csv(f"{R6}/{s}/tables/r6a_forms.csv") for s in ("svi", "wa_nw")], ignore_index=True)
    ok = D.groupby("form").pass_all.all(); npass = D.groupby("form").n_pass.sum()
    adopted = next((f for f in ORDER_A if bool(ok.get(f, False))), None)
    res = dict(adopted=adopted, pass_all=ok.to_dict(), n_pass=npass.to_dict(), written=R5.pt_now(),
               fails_by_form={f: {f"{r.region}|{r.setup}": r.fail_items for r in D[D.form == f].itertuples()} for f in AMEND},
               note="§10 Amendment 1 (post-result for task 1); task-2 column only")
    R5.dump(res, f"{R6}/adopted_amend.json")
    lg(f"AMENDMENT ADOPTED {adopted}; pass_all {ok.to_dict()}; n_pass {npass.to_dict()}")


def stage_adopt(lg):
    D = pd.concat([pd.read_csv(f"{R6}/{s}/tables/r6_forms.csv") for s in ("svi", "wa_nw")], ignore_index=True)
    ok = D.groupby("form").pass_all.all(); npass = D.groupby("form").n_pass.sum()
    adopted = next((f for f in ORDER if bool(ok.get(f, False))), None)
    res = dict(adopted=adopted, pass_all=ok.to_dict(), n_pass=npass.to_dict(), written=R5.pt_now())
    if adopted is None:
        best = sorted(ORDER, key=lambda f: (-npass[f], ORDER.index(f)))[0]
        res["best_not_passing"] = best
        res["fail_items"] = {f"{r.region}|{r.setup}": r.fail_items for r in D[D.form == best].itertuples()}
    res["fails_by_form"] = {f: {f"{r.region}|{r.setup}": r.fail_items for r in D[D.form == f].itertuples()} for f in MAIN}
    res["task2_source"] = adopted if adopted else "E4"
    R5.dump(res, f"{R6}/adopted.json")
    lg(f"ADOPTED {adopted}; pass_all {ok.to_dict()}; n_pass {npass.to_dict()}; task-2 ΔD source {res['task2_source']}")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--stage", required=True); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    global TB6
    if a.limit:
        LIMIT[0] = a.limit; TB6 = f"{TB6}/test"; os.makedirs(TB6, exist_ok=True)
    lg = log(f"model_{a.stage}") if a.stage in ("forms", "amend") else C.Log(f"{R6}/model_{a.stage}.log")
    t0 = time.time()
    {"forms": stage_forms, "adopt": stage_adopt, "amend": stage_amend, "adopt_amend": stage_adopt_amend}[a.stage](lg)
    lg(f"stage {a.stage} done {time.time() - t0:.0f}s")
    if a.limit:
        return
    mk = {"forms": f"MODEL_{C.REGION}_DONE", "adopt": "ADOPT_DONE", "amend": f"AMEND_{C.REGION}_DONE",
          "adopt_amend": "ADOPT_AMEND_DONE"}[a.stage]
    open(f"{R6}/{mk}", "w").write(R5.pt_now() + "\n")


if __name__ == "__main__":
    main()
