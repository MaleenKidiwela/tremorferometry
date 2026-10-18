#!/usr/bin/env python
"""s11_r5_q1 -- S11 round 5, task 3 and the Q1 part of task 5 (PLAN.md §9.3, §9.5, Amendment 1).
The yearly deep-weighted cycle Q1 refitted with ONE receiver-type term per unit (station with >= 20 cells, else site-rest;
coefficient 1; basis by cell count; exact basis) under two3_20; cell sets (a) near head cells (round-4 Q1 cells) and
(b) near + far (the F1 set, which contains the PGC-site stations as far stations). Groups P (PGC site), B (other PB
boreholes), O (others): Q1 per group, contribution shares, kernel-free lapse ratio per group vs deep / shallow / receiver
predictions. Gates G1-G3, date permutation, G-prior r5. Verdict DEEP / MIXED / SHALLOW / UNDECIDED by the fixed rule.
Task 5: Q1-r4 and Q1-R under every fit kernel; verdict recomputed under each plausible kernel.
    S11_REGION=... python s11_r5_q1.py [--nperm 500] [--test]
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
import s11_v2 as V
import s11_r3_common as R
import s11_r3_fit as LF
import s11_r4_common as R4
import s11_r5_common as R5

PGC = {"PB.B009", "PB.B010", "PB.B011", "CN.PGC"}
GROUPS = ("P", "B", "O")
GNAME = {"P": "PGC site (B009, B010, B011, PGC)", "B": "other PB boreholes", "O": "other stations"}
OUT = [R5.TB5]


def grp(ns):
    return "P" if ns in PGC else ("B" if ns.startswith("PB.") else "O")


def dpk(a, b):
    return abs((a - b + 182.625) % 365.25 - 182.625)


# ------------------------------------------------------------------------------ cell sets
def set_cells(which, setup, control=False):
    if which == "a":
        c = R.head_cells(control=control)
        c = R4.with_setup(c, setup) if setup != "ell40" else c
    else:
        import s11_dist_model as Mo
        c = Mo.all_cells(control=control)
        c["series"] = c.fam.astype(str) + "|" + c.netsta + "|" + c.epoch
        c[["a_s", "a_m", "a_d", "a_b"]] = R4.coefs(c, setup)
    c = c[np.isfinite(c.a_d) & np.isfinite(c.a_s) & np.isfinite(c.a_m) & np.isfinite(c.a_b)].reset_index(drop=True)
    c["grp"] = c.netsta.map(grp)
    return c


def lapse_set(which, setup):
    out = []
    for w in ("2035", "3560"):
        if which == "a":
            c = R.head_cells(win=w)
        else:
            import s11_dist_model as Mo
            c = Mo.all_cells(win=w)
            c["series"] = c.fam.astype(str) + "|" + c.netsta + "|" + c.epoch
        c = c.copy(); c[["a_s", "a_m", "a_d", "a_b"]] = R4.coefs(c, setup, win=w)
        c["win"] = w; c["epoch_orig"] = c.epoch; c["epoch"] = c.epoch + f"_w{w}"
        out.append(c)
    c = pd.concat(out, ignore_index=True)
    c = c[np.isfinite(c.a_d) & np.isfinite(c.a_s)].reset_index(drop=True)
    c["grp"] = c.netsta.map(grp)
    return c


# ------------------------------------------------------------------------------ models
def fit_r4(c):
    return LF.Qsmooth(LF.smooth(c).solve())


def q1_model(c, form="R", ridge=None, mid=None, fold_b=False):
    kw = dict(kind="smooth", MB=False, ridge=ridge, mid=mid)
    if form == "R":
        kw["unit"] = "one"
    elif form == "Ua":
        kw["unit"] = "as"
    elif form == "RMB":
        kw.update(unit="one", MB=True, fold_b=fold_b)
    return R5.Model5(c, **kw)


def Qof(L, d=None):
    F = L.solve(d)
    q = V.Q(F)
    return q, F


def q1_perm(c, form, nperm, seed, fold_b=False):
    """date permutation: episode mid-dates shuffled among the episodes (as round 3/4)"""
    rng = np.random.default_rng(seed)
    eu = c.episode.unique(); mids = V.U.mid.reindex(eu).to_numpy()
    out = []
    for _ in range(nperm):
        perm = dict(zip(eu, rng.permutation(mids)))
        out.append(Qof(q1_model(c, form, mid=c.episode.map(perm).to_numpy(), fold_b=fold_b))[0]["Q1"])
    return np.array(out)


def call_lapse(r, sr, pdp, psp, prc=1.0):
    zd, zs, zr = abs((r - pdp) / sr), abs((r - psp) / sr), abs((r - prc) / sr)
    if zd + 1 <= min(zs, zr):
        return "deep-like", zd, zs, zr
    if min(zs, zr) + 1 <= zd:
        return "shallow-like", zd, zs, zr
    return "undecided", zd, zs, zr


def lapse_test(CL, pk, mask=None):
    c = CL if mask is None else CL[mask]
    a, b = c[c.win == "2035"], c[c.win == "3560"]
    if len(a) < 10 or len(b) < 10:
        return None
    A1, s1 = V.inphase(a, pk); A2, s2 = V.inphase(b, pk)
    if A1 == 0 or not np.isfinite(A1):
        return None
    r = A2 / A1; sr = abs(r) * np.hypot(s1 / A1, s2 / A2) if A2 != 0 else abs(s2 / A1)
    pdp = b.a_d.mean() / a.a_d.mean(); psp = b.a_s.mean() / a.a_s.mean()
    call, zd, zs, zr = call_lapse(r, sr, pdp, psp)
    return dict(A2035=A1, A2035_sig=s1, A3560=A2, A3560_sig=s2, ratio=r, ratio_sig=sr, pred_deep=pdp, pred_shallow=psp, pred_receiver=1.0,
                z_deep=zd, z_shallow=zs, z_receiver=zr, call=call, n2035=len(a), n3560=len(b))


# ------------------------------------------------------------------------------ one analysis (cell set x setup)
def analyse(which, setup, nperm, lg, fold_b, full=True):
    c = set_cells(which, setup); cc = set_cells(which, setup, control=True)
    res = dict(region=C.REGION, set=which, setup=setup, n_cells=len(c))
    q4 = fit_r4(c); res.update({f"r4_{k}": v for k, v in q4.items() if k in ("Q1", "Q1_sig", "Q1_peak", "Q2", "Q2_sig")})
    L = q1_model(c, "R"); q, F = Qof(L)
    res.update(Q1=q["Q1"], Q1_sig=q["Q1_sig"], Q1_peak=q["Q1_peak"], Q1_z=q["Q1"] / q["Q1_sig"], Q2=q["Q2"], Q2_sig=q["Q2_sig"],
               cond=L.cond, singular=L.singular, n_units=len([u for u, t in L.unit_terms.items() if t]))
    for form in ("Ua", "RMB"):
        qq, _ = Qof(q1_model(c, form, fold_b=fold_b))
        res.update({f"{form}_Q1": qq["Q1"], f"{form}_Q1_sig": qq["Q1_sig"], f"{form}_Q1_peak": qq["Q1_peak"]})
    lg(f"[{which}|{setup}] Q1-r4 {q4['Q1']:.4f} ± {q4['Q1_sig']:.4f} (day {q4['Q1_peak']:.0f}); Q1-R {q['Q1']:.4f} ± {q['Q1_sig']:.4f} "
       f"(day {q['Q1_peak']:.0f}, z {q['Q1'] / q['Q1_sig']:.2f}); Q1-Ua {res['Ua_Q1']:.4f} ± {res['Ua_Q1_sig']:.4f}; Q1-R+MB {res['RMB_Q1']:.4f} ± {res['RMB_Q1_sig']:.4f}")
    if which == "b":
        try:
            import s11_dist_model as Mo
            qf = Mo.Qs(Mo.fit(c, "smooth", MB=True)); res.update(F1_MB_Q1=qf["Q1"], F1_MB_Q1_sig=qf["Q1_sig"])
        except Exception as e:
            lg(f"F1 reference failed: {e}")
    # G-prior (r5): ridge 1 / 2 % on the unit terms
    for sp in (1.0, 2.0):
        qq, _ = Qof(q1_model(c, "R", ridge=sp)); res[f"prior_ridge{sp:g}"] = abs(qq["Q1"] - q["Q1"]) / q["Q1_sig"]
    res["G_prior"] = bool(max(res["prior_ridge1"], res["prior_ridge2"]) < 0.1)
    # G1 leave-one-site-out
    ok1 = True; g1 = []
    for s in sorted(c.site.unique()):
        sub = c[c.site != s].reset_index(drop=True)
        if len(sub) < 20:
            continue
        qs, _ = Qof(q1_model(sub, "R"))
        if not (qs["Q1_sig"] > 0):
            continue
        w_ = abs(qs["Q1"] - q["Q1"]) < 2 * qs["Q1_sig"]; p_ = dpk(qs["Q1_peak"], q["Q1_peak"]) <= 60 if qs["Q1"] / qs["Q1_sig"] >= 1 else True
        ok1 &= bool(w_ and p_); g1.append(dict(left_out=s, Q1=qs["Q1"], Q1_sig=qs["Q1_sig"], Q1_peak=qs["Q1_peak"], passed=bool(w_ and p_)))
    res["G1"] = bool(ok1)
    # G2 family split
    qo, _ = Qof(q1_model(c[c.fam % 2 == 1].reset_index(drop=True), "R")); qe, _ = Qof(q1_model(c[c.fam % 2 == 0].reset_index(drop=True), "R"))
    w2 = abs(qo["Q1"] - qe["Q1"]) < 2 * np.hypot(qo["Q1_sig"], qe["Q1_sig"])
    pko = dpk(qo["Q1_peak"], qe["Q1_peak"]) <= 60 if min(qo["Q1"] / max(qo["Q1_sig"], 1e-12), qe["Q1"] / max(qe["Q1_sig"], 1e-12)) >= 1 else True
    res.update(G2=bool(w2 and pko), G2_odd=qo["Q1"], G2_odd_sig=qo["Q1_sig"], G2_even=qe["Q1"], G2_even_sig=qe["Q1_sig"])
    # G3 scrambled control
    qc, _ = Qof(q1_model(cc, "R")); res.update(G3=bool(qc["Q1"] < 2 * qc["Q1_sig"]), G3_Q1=qc["Q1"], G3_Q1_sig=qc["Q1_sig"])
    if nperm:
        perm = q1_perm(c, "R", nperm, C.SEED + 55)
        res["perm_p"] = float((1 + np.sum(perm >= q["Q1"])) / (1 + len(perm))); res["perm_q95"] = float(np.quantile(perm, 0.95))
    # groups: Q1-R per group, contribution shares of the in-phase Q1
    ph = 2 * np.pi * (q["Q1_peak"] - 1) / 365.25
    J = np.zeros(len(F["th"])); J[F["idx"][("D", "", "c1")]] = np.cos(ph); J[F["idx"][("D", "", "s1")]] = np.sin(ph)
    contrib = (J @ L.Op) * c.d.to_numpy(); tot = contrib.sum()
    gt = {}
    for g in GROUPS:
        m = (c.grp == g).to_numpy()
        r = dict(group=g, n_cells=int(m.sum()), n_stations=int(c[m].netsta.nunique()), share=float(contrib[m].sum() / tot) if tot != 0 else np.nan)
        sub = c[m].reset_index(drop=True)
        if len(sub) >= 30 and sub.episode.nunique() >= 5:
            try:
                qg, _ = Qof(q1_model(sub, "R")); r.update(Q1=qg["Q1"], Q1_sig=qg["Q1_sig"], Q1_peak=qg["Q1_peak"])
                q4g = fit_r4(sub); r.update(r4_Q1=q4g["Q1"], r4_Q1_sig=q4g["Q1_sig"], r4_Q1_peak=q4g["Q1_peak"])
            except Exception as e:
                r["err"] = str(e)[:100]
        gt[g] = r
    # lapse per group (kernel-free in-phase amplitudes at the pooled Q1-r4 peak)
    CL = lapse_set(which, setup)
    pk = q4["Q1_peak"]
    pooled = lapse_test(CL, pk)
    res["lapse_pooled"] = pooled
    for g in GROUPS:
        lt = lapse_test(CL, pk, (CL.grp == g).to_numpy())
        if lt:
            gt[g].update({f"lapse_{k}": v for k, v in lt.items()})
    res["groups"] = gt
    # Q1-RL: the receiver-term model on the joint-lapse cells
    try:
        qL, _ = Qof(q1_model(CL, "R")); res.update(RL_Q1=qL["Q1"], RL_Q1_sig=qL["Q1_sig"], RL_Q1_peak=qL["Q1_peak"])
    except Exception as e:
        lg(f"Q1-RL failed: {e}")
    # per-station kernel-free amplitudes (figure / column)
    st = []
    for ns in sorted(CL.netsta.unique()):
        a, b = CL[(CL.netsta == ns) & (CL.win == "2035")], CL[(CL.netsta == ns) & (CL.win == "3560")]
        if len(a) < 10 or len(b) < 10:
            continue
        A1, s1 = V.inphase(a, pk); A2, s2 = V.inphase(b, pk)
        st.append(dict(netsta=ns, grp=grp(ns), n2035=len(a), n3560=len(b), A2035=A1, A2035_sig=s1, A3560=A2, A3560_sig=s2,
                       a_d_ratio=b.a_d.mean() / a.a_d.mean(), a_s_ratio=b.a_s.mean() / a.a_s.mean()))
    res["verdict"], res["verdict_why"] = verdict(res)
    lg(f"[{which}|{setup}] gates G1 {res['G1']} G2 {res['G2']} G3 {res['G3']} ({qc['Q1']:.3f} ± {qc['Q1_sig']:.3f}) prior {res['prior_ridge1']:.3f}/"
       f"{res['prior_ridge2']:.3f}σ perm p {res.get('perm_p', np.nan):.3f}; pooled lapse {pooled['call'] if pooled else None} "
       f"({pooled['ratio']:.2f} ± {pooled['ratio_sig']:.2f}; deep {pooled['pred_deep']:.2f}, shallow {pooled['pred_shallow']:.2f})" if pooled else "")
    for g, r in gt.items():
        lg(f"   group {g} ({GNAME[g]}): cells {r['n_cells']}, share {r['share']:+.2f}, Q1-R {r.get('Q1', np.nan):.4f} ± {r.get('Q1_sig', np.nan):.4f} "
           f"(day {r.get('Q1_peak', np.nan):.0f}); lapse {r.get('lapse_ratio', np.nan):.2f} ± {r.get('lapse_ratio_sig', np.nan):.2f} "
           f"(deep {r.get('lapse_pred_deep', np.nan):.2f}, shallow {r.get('lapse_pred_shallow', np.nan):.2f}) -> {r.get('lapse_call')}")
    lg(f"[{which}|{setup}] VERDICT {res['verdict']}: {res['verdict_why']}")
    return res, pd.DataFrame(g1), pd.DataFrame(st)


def verdict(res):
    z = res["Q1_z"]; gt = res["groups"]; pooled = res.get("lapse_pooled") or {}
    pc = pooled.get("call")
    big = [g for g, r in gt.items() if np.isfinite(r.get("share", np.nan)) and abs(r["share"]) >= 0.25]
    calls = {g: gt[g].get("lapse_call") for g in gt}
    present = [g for g, r in gt.items() if np.isfinite(r.get("Q1", np.nan)) and r.get("Q1_sig", 0) > 0 and r["Q1"] / r["Q1_sig"] >= 2
               and dpk(r["Q1_peak"], res["Q1_peak"]) <= 60]
    why = (f"Q1-R/σ {z:.2f}; G1 {res['G1']} G2 {res['G2']} G3 {res['G3']}; pooled lapse {pc}; groups ≥ 25 % share {big} with lapse calls "
           f"{ {g: calls[g] for g in big} }; cycle ≥ 2σ (peak within 60 d) in groups {present}")
    deep = (z >= 3 and res["G1"] and res["G2"] and res["G3"] and pc != "shallow-like"
            and not any(calls[g] == "shallow-like" for g in big) and len(present) >= 2)
    if deep:
        return "DEEP", why
    dl = [g for g in gt if calls[g] == "deep-like" and (g in big or g in present)]
    sl = [g for g in gt if calls[g] == "shallow-like" and (g in big or g in present)]
    if (dl and sl and set(dl) != set(sl)) or (z >= 3 and pc == "shallow-like"):
        return "MIXED", why
    if (z < 2 or pc == "shallow-like") and not any(calls[g] == "deep-like" for g in big):
        return "SHALLOW", why
    return "UNDECIDED", why


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--nperm", type=int, default=500); ap.add_argument("--test", action="store_true")
    a = ap.parse_args()
    if a.test:
        R5.TB5 = f"{R5.TB5}/test"; os.makedirs(R5.TB5, exist_ok=True); a.nperm = min(a.nperm, 5)
    lg = R5.log("q1" + ("_test" if a.test else "")); t0 = time.time()
    fold_b = False
    p = f"{R5.R5}/selected.json"
    if os.path.exists(p):
        fold_b = bool(R5.FORMS[json.load(open(p))["selected"]].get("fold_b", False))
    lg(f"region {C.REGION}; Q1-R+MB uses fold_b={fold_b}")
    allres, band = [], []
    sets = ("a", "b")
    for which in sets:
        for v in R5.GEN:
            plaus = v in R5.PLAUS[C.REGION]
            if a.test and v not in ("two3_20",):
                continue
            t1 = time.time()
            if plaus:
                res, g1, st = analyse(which, v, a.nperm if v == "two3_20" else 0, lg, fold_b)
                allres.append(res)
                g1.to_csv(f"{R5.TB5}/q1_G1_{which}_{v}.csv", index=False); st.to_csv(f"{R5.TB5}/q1_stations_{which}_{v}.csv", index=False)
                pd.DataFrame(list(res["groups"].values())).to_csv(
                    f"{R5.TB5}/q1_groups_{which}_{v}.csv", index=False)
            else:
                c = set_cells(which, v)
                q4 = fit_r4(c); q, _ = Qof(q1_model(c, "R"))
                res = dict(region=C.REGION, set=which, setup=v, Q1=q["Q1"], Q1_sig=q["Q1_sig"], Q1_peak=q["Q1_peak"], Q1_z=q["Q1"] / q["Q1_sig"],
                           r4_Q1=q4["Q1"], r4_Q1_sig=q4["Q1_sig"], r4_Q1_peak=q4["Q1_peak"], verdict=None)
            band.append(dict(region=C.REGION, set=which, setup=v, plausible=plaus, Q1_R=res["Q1"], Q1_R_sig=res["Q1_sig"], Q1_R_peak=res["Q1_peak"],
                             Q1_R_z=res["Q1_z"], Q1_r4=res["r4_Q1"], Q1_r4_sig=res["r4_Q1_sig"], Q1_r4_peak=res["r4_Q1_peak"],
                             Q1_r4_z=res["r4_Q1"] / res["r4_Q1_sig"], verdict=res.get("verdict")))
            lg(f"  band {which}|{v}: Q1-r4 {res['r4_Q1']:.4f} ± {res['r4_Q1_sig']:.4f}, Q1-R {res['Q1']:.4f} ± {res['Q1_sig']:.4f}, verdict {res.get('verdict')} "
               f"({time.time() - t1:.0f}s)")
    B = pd.DataFrame(band); B.to_csv(f"{R5.TB5}/q1_band.csv", index=False)
    head = {r["set"]: r for r in allres if r["setup"] == "two3_20"}
    out = dict(region=C.REGION, headline_set="b" if C.REGION == "WA_NW" else "a", results=allres, written=R5.pt_now())
    for which in sets:
        b = B[(B.set == which) & B.plausible]
        out[f"band_{which}"] = dict(Q1_R_min=float(b.Q1_R.min()), Q1_R_max=float(b.Q1_R.max()), z_min=float(b.Q1_R_z.min()), z_max=float(b.Q1_R_z.max()),
                                    Q1_r4_min=float(b.Q1_r4.min()), Q1_r4_max=float(b.Q1_r4.max()), z_r4_min=float(b.Q1_r4_z.min()),
                                    z_r4_max=float(b.Q1_r4_z.max()), verdicts=b.verdict.tolist(),
                                    holds=bool(len(set(b.verdict)) == 1 and (b.Q1_R_z >= 3).nunique() == 1))
        lg(f"BAND set {which}: Q1-R {b.Q1_R.min():.4f}–{b.Q1_R.max():.4f} (z {b.Q1_R_z.min():.2f}–{b.Q1_R_z.max():.2f}); Q1-r4 "
           f"{b.Q1_r4.min():.4f}–{b.Q1_r4.max():.4f}; verdicts {b.verdict.tolist()} -> holds {out[f'band_{which}']['holds']}")
    R5.dump(out, f"{R5.TB5}/q1_verdict.json")
    lg(f"q1 done {time.time() - t0:.0f}s")
    if not a.test:
        open(f"{R5.R5}/Q1_{C.REGION}_DONE", "w").write(R5.pt_now() + "\n")


if __name__ == "__main__":
    main()
