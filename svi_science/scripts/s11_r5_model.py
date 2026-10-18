#!/usr/bin/env python
"""s11_r5_model -- S11 round 5, tasks 1, 2 and the Q3 part of task 5 (PLAN.md §9.1, §9.2, §9.5).
  --stage forms : per region, forms r4_HnoSc / r4_H (reproductions) and E1-E4 (one shallow term per unit, exact basis) under
                  NOM (ell40) and L-SEL (two3_20): natural Q3, gates (G-prior r5, control, loso, family, site), T1-T5 from the
                  round-4 injections.
  --stage select: both regions -> the selected form (fixed pass-count rule) -> round5/selected.json
  --stage sc    : per region, the selected form + smooth regional Sc (B-spline, knots 2/3/4 yr) with the full suite; Pcal
                  column; Q3 systematic band over the fit kernels (task 5).
  --stage adopt : both regions -> Sc-spline adoption (goal leak <= 0.35 and sigma/major <= round-4 HnoSc) -> round5/sc_adoption.json
    S11_REGION=... python s11_r5_model.py --stage forms|sc ;  python s11_r5_model.py --stage select|adopt
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
import s11_r4_common as R4
import s11_r5_common as R5

SETUPS = ("ell40", "two3_20")
LIMIT = [0]
RULE_ITEMS = ("T1", "T2", "T3", "T5", "G_prior", "G_ctrl", "G_loso", "G_fam")


def ctrl_b1(cs):
    import s11_r4_combined as CB
    return CB._ctrl_b1(cs, "ell40")


def run_form(H, inj, majors, name, kw, setup, lg, tag, exact=True, ridge_kw=None, suite=True):
    fac = R5.factory(kw)
    out = f"{R5.TB5}/{tag}_{name}_{setup}"
    G, T, sw, F, L, cs = R5.gates(H, setup, majors, fac, name, lg, exact=exact, ridge_kw=ridge_kw, ctrl_b1_fn=ctrl_b1, out=out)
    row = dict(region=C.REGION, form=name, **G)
    if "2019-11" in set(T.episode):
        x = T.set_index("episode").loc["2019-11"]; row.update(dD_2019_11=float(x.dD), sig_2019_11=float(x.sig))
    TR = O = None
    if suite:
        TR, O = R5.evaluate(H, inj, setup, majors, fac, sw, name, lg)
        v = R5.rules(TR, O, setup)
        row.update({f"rule_{k}": (json.dumps(x) if isinstance(x, list) else x) for k, x in v.items()})
        TR = TR.assign(region=C.REGION, form=name); O = O.assign(region=C.REGION, form=name)
    row["G_ctrl"] = row.get("ctrl_g_pass"); row["G_loso"] = row.get("G_loso_pass"); row["G_fam"] = row.get("G_fam_pass")
    row["G_prior"] = row.get("G_prior_pass")
    for k in ("T1", "T2", "T3", "T5"):
        row[k] = row.get(f"rule_{k}")
    row["n_pass"] = int(sum(bool(row.get(k)) for k in RULE_ITEMS))
    return row, TR, O, T.assign(form=name, setup=setup), F, L, cs


def stage_forms(lg):
    majors = R4.majors()
    H, inj = R5.load_H(lg, LIMIT[0])
    rows, TRs, Os, Ts = [], [], [], []
    for name in ("r4_HnoSc", "r4_H", "E1", "E2", "E3", "E4", "E2b", "E4b", "E2m"):
        kw = R5.FORMS[name]
        for setup in SETUPS:
            t0 = time.time()
            exact = kw.get("ridge") is None
            row, TR, O, T, F, L, cs = run_form(H, inj, majors, name, kw, setup, lg, "forms", exact=exact,
                                              ridge_kw=None if exact else dict(ridge=kw["ridge"]))
            rows.append(row); TRs.append(TR); Os.append(O); Ts.append(T)
            lg(f"FORM {name}|{setup}: σ/major {row['sig_major_median']:.4f} ({row['n_major_sig_le_0p1']}/{row['n_majors']} ≤ 0.1), composite "
               f"{row['composite']:+.4f} ± {row['composite_sig']:.4f}, prior max {row['prior_max']:.3f}σ, cond {row['cond']:.1e}, "
               f"pass {row['n_pass']}/8 [{', '.join(k for k in RULE_ITEMS if not row.get(k))} fail] ({time.time() - t0:.0f}s)")
    D = pd.DataFrame(rows)
    # reproduction check vs round 4 (HnoSc and H natural values)
    g4 = pd.read_csv(f"{R4.TB4}/r4_gates.csv")
    chk = []
    for name, lab in (("r4_HnoSc", "HnoSc"), ("r4_H", "H")):
        for setup in SETUPS:
            a = D[(D.form == name) & (D.setup == setup)].iloc[0]; b = g4[(g4.model == lab) & (g4.setup == setup)]
            if len(b):
                b = b.iloc[0]
                chk.append(dict(form=name, setup=setup, d_sig_major=float(a.sig_major_median - b.sig_major_median),
                                d_composite=float(a.composite - b.composite), d_composite_sig=float(a.composite_sig - b.composite_sig)))
    CK = pd.DataFrame(chk); CK.to_csv(f"{R5.TB5}/forms_reproduction_check.csv", index=False)
    mx = float(np.nanmax(np.abs(CK[["d_sig_major", "d_composite", "d_composite_sig"]].to_numpy()))) if len(CK) else np.nan
    lg(f"reproduction of round-4 H / HnoSc (σ/major, composite, σ): max |Δ| {mx:.1e} -> {'OK' if mx <= 1e-6 else 'MISMATCH'}")
    D["repro_max_abs"] = mx
    D.to_csv(f"{R5.TB5}/r5_forms.csv", index=False)
    pd.concat(TRs, ignore_index=True).to_csv(f"{R5.TB5}/r5_forms_rules.csv", index=False)
    pd.concat(Os, ignore_index=True).to_csv(f"{R5.TB5}/r5_forms_onesite.csv", index=False)
    pd.concat(Ts, ignore_index=True).to_csv(f"{R5.TB5}/r5_forms_dD.csv", index=False)


def stage_select(lg):
    parts = []
    for sub in ("svi", "wa_nw"):
        p = f"{R5.R5}/{sub}/tables/r5_forms.csv"
        parts.append(pd.read_csv(p))
    D = pd.concat(parts, ignore_index=True)
    E = D[D.form.isin(R5.ELIGIBLE)]
    sc = E.groupby("form").n_pass.sum()
    sig = E[E.setup == "two3_20"].groupby("form").sig_major_median.mean()
    order = sorted(sc.index, key=lambda f: (-sc[f], round(sig[f], 6), R5.TIE.index(f)))
    sel = order[0]
    res = dict(selected=sel, score=sc.to_dict(), mean_lsel_sig_major=sig.to_dict(), order=order, written=R5.pt_now())
    t = E[E.form == sel]
    res["task1_pass"] = bool(all(bool(r.G_prior) and all(bool(getattr(r, k)) for k in ("T1", "T2", "T3", "T5", "G_ctrl", "G_loso", "G_fam"))
                                 for r in t[t.setup == "two3_20"].itertuples()))
    res["G_prior_by_combo"] = {f"{r.region}|{r.setup}": dict(prior_max=float(r.prior_max), passed=bool(r.G_prior)) for r in t.itertuples()}
    R5.dump(res, f"{R5.R5}/selected.json")
    lg(f"SELECTED form {sel}: scores {sc.to_dict()}, mean L-SEL σ/major {sig.round(4).to_dict()}; task-1 pass {res['task1_pass']}")


def selected():
    if LIMIT[0]:
        return "E2"                              # code test only
    return json.load(open(f"{R5.R5}/selected.json"))["selected"]


def stage_sc(lg):
    majors = R4.majors()
    H, inj = R5.load_H(lg, LIMIT[0])
    sel = selected(); kw = R5.FORMS[sel]
    lg(f"selected form {sel}: {kw}")
    rows, TRs, Os, Ts = [], [], [], []
    for k in (2, 3, 4):
        for setup in SETUPS:
            t0 = time.time()
            row, TR, O, T, F, L, cs = run_form(H, inj, majors, f"{sel}+Scs{k}", {**kw, "Sc": ("spline", k)}, setup, lg, "sc")
            row["knot_yr"] = k; row["hnosc_r4_sig"] = R5.HNOSC_R4[(C.REGION, setup)]
            row["goal_leak"] = bool(row["rule_T2_val"] <= 0.35); row["goal_sig"] = bool(row["sig_major_median"] <= row["hnosc_r4_sig"] + 1e-12)
            row["goal"] = bool(row["goal_leak"] and row["goal_sig"])
            rows.append(row); TRs.append(TR); Os.append(O); Ts.append(T)
            lg(f"SCS k={k}|{setup}: σ/major {row['sig_major_median']:.4f} (r4 HnoSc {row['hnosc_r4_sig']:.3f}), leak {row['rule_T2_val']:.3f}, "
               f"composite {row['composite']:+.4f} ± {row['composite_sig']:.4f}, spline cols {row.get('spline_kept')}, goal {row['goal']} ({time.time() - t0:.0f}s)")
    # Pcal column (selected exact form; ridge per group from the data and scrambled-control exact fits)
    for setup in SETUPS:
        cs, good = R5.setup_cells(H, setup); cs = cs[good].reset_index(drop=True)
        fac = R5.factory(kw)
        Fd = fac(cs).solve()
        cc, _ = R4.identifiable(R4.control_cells(cs))
        rd = R5.coef_rms(Fd); rc = R5.coef_rms(fac(cc).solve()) if len(cc) >= 20 and cc.episode.nunique() >= 3 else {}
        ridge = {g: max(rd.get(g, 0.0), rc.get(g, 0.0)) for g in set(rd) | set(rc)}
        ridge = {g: v for g, v in ridge.items() if v > 0}
        lg(f"Pcal ridge [{setup}] (data RMS {rd}, control RMS {rc}) -> {ridge}")
        row, TR, O, T, F, L, cs2 = run_form(H, inj, majors, f"{sel}+Pcal", {**kw, "ridge": ridge}, setup, lg, "sc", exact=False,
                                            ridge_kw=dict(ridge=ridge))
        row["pcal_ridge"] = json.dumps(ridge); rows.append(row); TRs.append(TR); Os.append(O); Ts.append(T)
        lg(f"PCAL|{setup}: σ/major {row['sig_major_median']:.4f}, composite {row['composite']:+.4f} ± {row['composite_sig']:.4f}, "
           f"prior ×2/×0.5 max {row['prior_max']:.3f}σ")
    D = pd.DataFrame(rows); D.to_csv(f"{R5.TB5}/r5_sc.csv", index=False)
    pd.concat(TRs, ignore_index=True).to_csv(f"{R5.TB5}/r5_sc_rules.csv", index=False)
    pd.concat(Os, ignore_index=True).to_csv(f"{R5.TB5}/r5_sc_onesite.csv", index=False)
    pd.concat(Ts, ignore_index=True).to_csv(f"{R5.TB5}/r5_sc_dD.csv", index=False)
    # ---- task 5: Q3 systematic band (natural values under every fit setup) for the selected form and each Sc-spline k
    band = []
    for name, fkw in [(sel, kw)] + [(f"{sel}+Scs{k}", {**kw, "Sc": ("spline", k)}) for k in (2, 3, 4)]:
        fac = R5.factory(fkw)
        for v in R5.GEN:
            cs, good = R5.setup_cells(H, v); cs = cs[good].reset_index(drop=True)
            F = fac(cs).solve(); T, comp = R5.q3(F, majors); mj = T[T.major]
            r = dict(region=C.REGION, form=name, fit_setup=v, plausible=v in R5.PLAUS[C.REGION], composite=comp["composite"],
                     composite_sig=comp["composite_sig"], sig_major_median=float(mj.sig.median()), n_major_sig_le_0p1=int((mj.sig <= 0.1).sum()),
                     n_majors_ge2sig=int((np.abs(mj.dD / mj.sig) >= 2).sum()), n_majors_ge3sig=int((np.abs(mj.dD / mj.sig) >= 3).sum()),
                     chi2_majors=float(np.mean((mj.dD / mj.sig) ** 2)))
            if "2019-11" in set(T.episode):
                x = T.set_index("episode").loc["2019-11"]; r.update(dD_2019_11=float(x.dD), sig_2019_11=float(x.sig))
            band.append(r)
            T.assign(form=name, fit_setup=v).to_csv(f"{R5.TB5}/band_dD_{name}_{v}.csv", index=False)
    B = pd.DataFrame(band)
    # deep-pulse recovery under the headline setup (L-SEL) for the generators in P (from the suite above / forms stage)
    fr = pd.read_csv(f"{R5.TB5}/r5_forms_rules.csv")
    sr = pd.read_csv(f"{R5.TB5}/r5_sc_rules.csv")
    allr = pd.concat([fr, sr], ignore_index=True)
    rec = allr[(allr.setup == "two3_20")].groupby(["form", "gen"]).T1_Pd_rec.first().unstack()
    B.to_csv(f"{R5.TB5}/r5_band_q3.csv", index=False)
    rec.to_csv(f"{R5.TB5}/r5_band_T1_lsel.csv")
    for name in B.form.unique():
        b = B[(B.form == name) & B.plausible]
        within = bool((np.abs(b.composite) < 2 * b.composite_sig).all())
        bound = float(np.max(np.abs(b.composite) + 2 * b.composite_sig))
        t1p = rec.loc[name, R5.PLAUS[C.REGION]] if name in rec.index else pd.Series(dtype=float)
        lg(f"BAND {name}: composite over P {b.composite.min():+.4f} … {b.composite.max():+.4f} (σ {b.composite_sig.min():.4f}–"
           f"{b.composite_sig.max():.4f}); all within 2σ of 0: {within}; 2σ bound {bound:.3f} %; T1 (L-SEL fit) over P "
           f"{t1p.min():.2f}–{t1p.max():.2f} -> size band {bound / max(t1p.max(), 1e-9):.3f}–{bound / max(t1p.min(), 1e-9):.3f} %")


def stage_adopt(lg):
    parts = [pd.read_csv(f"{R5.R5}/{sub}/tables/r5_sc.csv") for sub in ("svi", "wa_nw")]
    D = pd.concat(parts, ignore_index=True); D = D[D.knot_yr.notna()]
    ok = D.groupby("knot_yr").goal.all()
    sig = D[D.setup == "two3_20"].groupby("knot_yr").sig_major_median.mean()
    cand = [k for k in ok.index if ok[k]]
    adopted = int(min(cand, key=lambda k: sig[k])) if cand else None
    fail = {int(k): [f"{r.region}|{r.setup}: leak {r.rule_T2_val:.3f}, σ {r.sig_major_median:.3f} vs {r.hnosc_r4_sig:.3f}"
                     for r in D[(D.knot_yr == k) & ~D.goal].itertuples()] for k in ok.index}
    res = dict(adopted_knot_yr=adopted, goal_all=ok.to_dict(), mean_lsel_sig=sig.to_dict(), failing=fail, written=R5.pt_now())
    R5.dump(res, f"{R5.R5}/sc_adoption.json")
    lg(f"Sc-spline adoption: {adopted}; goal per k {ok.to_dict()}; failing {fail}")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--stage", required=True); ap.add_argument("--limit", type=int, default=0); a = ap.parse_args()
    if a.limit:                                  # code test: outputs to tables/test, no marker
        LIMIT[0] = a.limit; R5.TB5 = f"{R5.TB5}/test"; os.makedirs(R5.TB5, exist_ok=True)
    lg = R5.log(f"model_{a.stage}") if a.stage in ("forms", "sc") else C.Log(f"{R5.R5}/model_{a.stage}.log")
    t0 = time.time()
    {"forms": stage_forms, "select": stage_select, "sc": stage_sc, "adopt": stage_adopt}[a.stage](lg)
    lg(f"stage {a.stage} done {time.time() - t0:.0f}s")
    if a.limit:
        return
    mk = {"forms": f"MODEL_FORMS_{C.REGION}_DONE", "sc": f"MODEL_SC_{C.REGION}_DONE", "select": "SELECT_DONE", "adopt": "ADOPT_DONE"}[a.stage]
    open(f"{R5.R5}/{mk}", "w").write(R5.pt_now() + "\n")


if __name__ == "__main__":
    main()
