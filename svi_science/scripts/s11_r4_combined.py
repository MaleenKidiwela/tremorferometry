#!/usr/bin/env python
"""s11_r4_combined -- S11 round 4, task 1 (PLAN.md §8.1): the pre-registered combined model
  near B2 + far B2 cells; D(e) + regional Sc(e) per episode with exact sum-to-zero (no prior); S per site smooth; receiver
  term R per station (coefficient 1; annual + linear); regional M, B free.
Natural values for H and the columns; waveform-level injection suite through the B2 graph on the near (round 3) and far
(round 4) cc-weighted stores, generated with every kernel of V6 + AVG6 (Pd, Ps, Psall, B, A per kernel; Rann receiver
annual); rules T1-T6 for NOM (headline) and L-SEL (column); gates G-ctrl, G-loso, G-fam, G-prior, G-site.
    S11_REGION=... python s11_r4_combined.py [--workers 12] [--limit N]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s11_common as C
import s11_dist_common as D
import s11_r4_common as R4
import s11_r3_common as R
import s11_r3_inject as J
import s11_r3_fit as LF
import s11_paths as P

GEN = R4.V6 + ["AVG6"]
STEP_DATE, PEAK = "2015-01-01", 233.0
_G = {}
OUT = [R4.TB4]


# ------------------------------------------------------------------------------ cells
def load_H(lg):
    nr = R4.b2_cells(False); fr = R4.b2_cells(True)
    H = pd.concat([nr, fr], ignore_index=True) if len(fr) else nr
    H, nd = R4.identifiable(H)
    lg(f"cells: near B2 {len(nr)}, far B2 {len(fr)}; H {len(H)} after the >= 2-site episode rule ({nd} dropped); "
       f"episodes {H.episode.nunique()}, sites {H.site.nunique()}, stations {H.netsta.nunique()} (far {H[H.far].netsta.nunique()})")
    return H


# ------------------------------------------------------------------------------ injection engine (B2 graph, near + far stores)
def runs_fn(key, st, unit):
    gc = _G["gcoef"]; majors = _G["majors"]
    t = pd.to_datetime(pd.Series(st["days"]))
    post = (t >= STEP_DATE).to_numpy().astype(float)
    ann = np.cos(2 * np.pi * (t.dt.dayofyear.to_numpy() - PEAK) / 365.25)
    ep = np.array([unit.get(d, None) for d in st["days"]], dtype=object)
    runs = {"nat": None, "Rann": -0.001 * ann}
    mall = np.isin(ep, majors).astype(float)
    for g in GEN:
        if key not in gc[g]:
            continue
        a_s, a_d = gc[g][key]
        if not (np.isfinite(a_s) and np.isfinite(a_d)):
            continue
        runs[f"A|{g}"] = -(0.003 * a_d) * post
        runs[f"B|{g}"] = -(0.001 * a_s) * ann
        if mall.any():
            runs[f"Psall|{g}"] = -(0.001 * a_s) * mall
        for e in majors:
            m = (ep == e).astype(float)
            if m.any():
                runs[f"Pd|{g}|{e}"] = -(0.003 * a_d) * m
                runs[f"Ps|{g}|{e}"] = -(0.001 * a_s) * m
    return runs


def _worker(key):
    try:
        st = _G["store"][key]; nodes = _G["nodes"].get(key, [])
        if len(nodes) < 3:
            return key, None
        out = {}
        for name, amp in runs_fn(key, st, _G["unit"]).items():
            Sx = J.episode_stacks_from_days(st, _G["unit"], amp, nodes)
            out[name] = J.graph_from_stacks(Sx, st["tS"])[0]
        return key, (list(nodes), out)
    except Exception as e:
        import traceback
        print("ERROR", key, e, traceback.format_exc()[-500:], flush=True)
        return key, None


def run_inj(H, majors, workers, lg):
    import multiprocessing as mp
    unit = P.unit_map() if C.REGION == "SVI" else dict(zip(pd.read_csv(f"{C.MTB}/active_days.csv").date,
                                                            pd.read_csv(f"{C.MTB}/active_days.csv").episode))
    sk = H.drop_duplicates("series")[["fam", "netsta", "epoch", "series", "far"]].reset_index(drop=True)
    gcoef = {}
    for g in GEN:
        A = R4.coefs(sk, g)
        gcoef[g] = {k: (a[0], a[2]) for k, a in zip(sk.series, A)}
    vals = {}
    for far in (False, True):
        q = sk[sk.far == far]
        if not len(q):
            continue
        E = R4.r4_episode(far)
        nodes = {}
        for (f, n, e), g in E[(E.n >= R4.NODE_MIN) & E.b1_cc.notna()].groupby(["fam", "netsta", "epoch"]):
            nodes[f"{f}|{n}|{e}"] = list(g.episode)
        D.set_grid(D.FAR_GRID if far else D.NEAR_GRID)
        paths = [f"{R4.TB4}/store_ccw_far"] if far else ([f"{R.TB3}/store_ccw"] + ([f"{R.TB3}/store_ccw_ext"] if C.REGION == "SVI" else []))
        for ns in sorted(q.netsta.unique()):
            store = {}
            for sp in paths:
                f_ = f"{sp}/{ns}.npz"
                if os.path.exists(f_):
                    z = np.load(f_)
                    if len(z["tau"]) != S.NW:
                        raise RuntimeError(f"grid mismatch {f_}: {len(z['tau'])} vs {S.NW}")
                    store.update(J.load_store_file(f_))
            keys = [k for k in store if k in set(q.series)]
            if not keys:
                continue
            _G.clear(); _G.update(store=store, unit=unit, nodes=nodes, gcoef=gcoef, majors=majors)
            t0 = time.time()
            with mp.get_context("fork").Pool(workers) as pool:
                for k, o in pool.imap_unordered(_worker, keys, chunksize=1):
                    if o is not None:
                        vals[k] = o
            lg(f"  inject {'far' if far else 'near'} {ns}: {len(keys)} series {time.time() - t0:.0f}s")
    D.set_grid(D.NEAR_GRID)
    return vals


def vectors(vals, H):
    """-> dict run -> d (dv/v %) per cell of H (run missing for a series -> natural)"""
    key = H.series.to_numpy(); ep = H.episode.to_numpy()
    pos = np.full(len(H), -1); ok = np.zeros(len(H), bool)
    for i, (k, e) in enumerate(zip(key, ep)):
        if k in vals and e in vals[k][0]:
            pos[i] = vals[k][0].index(e); ok[i] = True
    runs = sorted({r for v in vals.values() for r in v[1]})
    out = {}
    for run in runs:
        d = np.full(len(H), np.nan)
        for i in np.flatnonzero(ok):
            o = vals[key[i]][1]
            v = o.get(run, o["nat"])[pos[i]]
            if np.isfinite(v) and abs(v) < C.RAIL:
                d[i] = -100 * float(v)
        out[run] = d
    return out


# ------------------------------------------------------------------------------ evaluation
def resp(L, F0, x, majors):
    th = L.Op @ x
    gm = {e: LF.dvec(F0, e) for e in majors}
    gm = {e: g for e, g in gm.items() if g is not None}
    gc = LF.comp_vec(F0, majors)
    return {e: float(g @ th) for e, g in gm.items()}, (float(gc @ th) if gc is not None else np.nan)


def evaluate(H, inj, setup, majors, lg, label="H", sites_w=None, model_kw=None):
    """rules T1-T6 for one setup (generating kernels GEN); x = injected - natural, railed -> 0"""
    model_kw = model_kw or {}
    cs = R4.with_setup(H, setup) if setup != "ell40" else H
    good = np.isfinite(cs.a_d.to_numpy()) & np.isfinite(cs.a_s.to_numpy()) & np.isfinite(inj["nat"])
    cs = cs[good].reset_index(drop=True)
    L = R4.Comb(cs, **model_kw); F0 = L.delta(np.zeros(len(cs)))
    dn = inj["nat"][good]
    X = lambda run, mask=None: (lambda v: np.where(np.isfinite(v), v - dn, 0.0) * (1.0 if mask is None else mask))(inj.get(run, inj["nat"])[good])
    site = cs.site.to_numpy()
    rows = []
    for g in GEN:
        rec = []; lk = []
        for e in majors:
            if f"Pd|{g}|{e}" not in inj:
                continue
            r_, _ = resp(L, F0, X(f"Pd|{g}|{e}"), [e])
            if e in r_:
                rec.append(r_[e] / 0.3)
            r_, _ = resp(L, F0, X(f"Ps|{g}|{e}"), [e])
            if e in r_:
                lk.append(abs(r_[e]) / 0.1)
        _, cps = resp(L, F0, X(f"Psall|{g}"), majors)
        rb, cb = resp(L, F0, X(f"B|{g}"), majors)
        ra, ca = resp(L, F0, X(f"A|{g}"), majors)
        r = dict(label=label, setup=R.setup_name(setup) if not isinstance(setup, str) else setup, gen=g,
                 T1_Pd_rec=float(np.median(rec)) if rec else np.nan, T1_lo=float(np.min(rec)) if rec else np.nan,
                 T1_hi=float(np.max(rec)) if rec else np.nan, T2_Ps_leak=float(np.median(lk)) if lk else np.nan,
                 T2_hi=float(np.max(lk)) if lk else np.nan, Psall_comp_leak=cps / 0.1, T5_B_comp=abs(cb) / 0.1,
                 B_major_med=float(np.median(np.abs(list(rb.values())))) / 0.1 if rb else np.nan, A_comp=ca / 0.3,
                 A_major_med=float(np.median(list(ra.values()))) / 0.3 if ra else np.nan, n_major=len(rec))
        rows.append(r)
    T = pd.DataFrame(rows)
    # one-site (T3) and receiver annual (T4), generated with the in-sample kernel of the setup
    gin = setup if isinstance(setup, str) and setup in GEN else "ell40"
    ons = []
    sites = sites_w.index[sites_w >= 0.05] if sites_w is not None else []
    for s in sites:
        m = (site == s).astype(float)
        lk = []
        for e in majors:
            if f"Ps|{gin}|{e}" not in inj:
                continue
            r_, _ = resp(L, F0, X(f"Ps|{gin}|{e}", m), [e])
            if e in r_:
                lk.append(abs(r_[e]) / 0.1)
        rr, cr = resp(L, F0, X("Rann", m), majors)
        ons.append(dict(label=label, setup=r["setup"], site=s, weight_share=float(sites_w[s]), T3_onesite_leak=float(np.median(lk)) if lk else np.nan,
                        T3_hi=float(np.max(lk)) if lk else np.nan, T4_rann_comp=abs(cr) / 0.1,
                        T4_rann_major_med=float(np.median(np.abs(list(rr.values())))) / 0.1 if rr else np.nan))
    rr, cr = resp(L, F0, X("Rann"), majors)
    ons.append(dict(label=label, setup=r["setup"], site="ALL", weight_share=1.0, T3_onesite_leak=np.nan, T4_rann_comp=abs(cr) / 0.1,
                    T4_rann_major_med=float(np.median(np.abs(list(rr.values())))) / 0.1 if rr else np.nan))
    O = pd.DataFrame(ons)
    for x in rows:
        lg(f"  [{label}|{x['setup']}] gen {x['gen']}: T1 {x['T1_Pd_rec']:.3f} ({x['T1_lo']:.2f}-{x['T1_hi']:.2f}) T2 {x['T2_Ps_leak']:.3f} "
           f"Psall {x['Psall_comp_leak']:+.3f} T5 {x['T5_B_comp']:.3f} A comp {x['A_comp']:+.3f}")
    for x in ons:
        lg(f"  [{label}|{x['setup']}] site {x['site']} (w {x['weight_share']:.2f}): T3 {x['T3_onesite_leak']:.3f} T4 comp {x['T4_rann_comp']:.3f} "
           f"major {x['T4_rann_major_med']:.3f}")
    return T, O


def rule_verdict(T, O, gin):
    ins = T[T.gen == gin].iloc[0]
    v = dict(T1=bool(0.8 <= ins.T1_Pd_rec <= 1.1), T2=bool(ins.T2_Ps_leak <= 0.35), T2_fixed=bool(ins.T2_Ps_leak <= 0.10),
             T5=bool(ins.T5_B_comp <= 0.35))
    o = O[O.site != "ALL"]
    v["T3"] = bool((o.T3_onesite_leak <= 0.35).all()) if len(o) else None
    v["T4"] = bool(((O.T4_rann_comp <= 0.35) & (O.T4_rann_major_med <= 0.35)).all())
    fail = []
    for r in T.itertuples():
        for col, lo, hi, key in (("T1_Pd_rec", 0.8, 1.1, "T1"), ("T2_Ps_leak", -1, 0.35, "T2"), ("T5_B_comp", -1, 0.35, "T5")):
            val = getattr(r, col)
            ok = (lo <= val <= hi) if v[key] else abs(val - getattr(ins, col)) <= 0.10
            if not ok:
                fail.append(f"{key}:{r.gen}")
    v["T6"] = not fail; v["T6_failing"] = fail
    v["T1_range_over_gen"] = [float(T.T1_Pd_rec.min()), float(T.T1_Pd_rec.max())]
    return v


def gates(H, setup, majors, lg, label, model_kw=None):
    model_kw = model_kw or {}
    cs = R4.with_setup(H, setup) if setup != "ell40" else H
    cs = cs[np.isfinite(cs.a_d) & np.isfinite(cs.a_s)].reset_index(drop=True)
    F = R4.Comb(cs, **model_kw).solve(); T, comp = R4.q3(F, majors); Tm = T.set_index("episode"); mj = T[T.major]
    G = R4.summary(F, T, comp, label)
    G["acceptance"] = bool(G["n_major_sig_le_0p1"] >= 5)
    w = 1 / cs.sig ** 2; sw = (w.groupby(cs.site).sum() / w.sum()).sort_values(ascending=False)
    G["top_site"], G["top_site_share"] = sw.index[0], float(sw.iloc[0]); G["G_site_pass"] = bool(sw.iloc[0] <= 0.5)
    for f_ in (2.0, 0.5):
        Tf, cf = R4.q3(R4.Comb(cs, prior=f_, **model_kw).solve(), majors); Tf = Tf.set_index("episode")
        e_ = [e for e in mj.episode if e in Tf.index]
        G[f"prior_x{f_:g}"] = float(np.max(np.abs((Tf.dD[e_] - Tm.dD[e_]) / Tm.sig[e_])))
        G[f"prior_x{f_:g}_comp"] = abs(cf["composite"] - comp["composite"]) / comp["composite_sig"]
    G["G_prior_pass"] = bool(max(G["prior_x2"], G["prior_x0.5"], G["prior_x2_comp"], G["prior_x0.5_comp"]) < 0.1)
    for kind in ("g", "b1"):
        cc = R4.control_cells(cs if kind == "g" else _ctrl_b1(cs, setup))
        cc, _ = R4.identifiable(cc)
        if len(cc) < 20 or cc.episode.nunique() < 3:
            G[f"ctrl_{kind}_pass"] = None; G[f"ctrl_{kind}_cells"] = len(cc); continue
        Fc = R4.Comb(cc, **model_kw).solve(); Tc, ccomp = R4.q3(Fc, majors); Tc = Tc.set_index("episode")
        meas = [e for e in mj.episode if e in Tc.index and Tc.n_cells[e] >= 5]
        frac = float(np.mean([abs(Tc.dD[e]) < 2 * Tc.sig[e] for e in meas])) if meas else np.nan
        G.update({f"ctrl_{kind}_cells": len(cc), f"ctrl_{kind}_n": len(meas), f"ctrl_{kind}_frac": frac,
                  f"ctrl_{kind}_comp": ccomp["composite"] if ccomp else np.nan, f"ctrl_{kind}_comp_sig": ccomp["composite_sig"] if ccomp else np.nan,
                  f"ctrl_{kind}_pass": bool(frac >= 0.9 and (ccomp is None or abs(ccomp["composite"]) < 2 * ccomp["composite_sig"])) if meas else None})
        Tc.reset_index().to_csv(f"{OUT[0]}/r4_control_{kind}_{label}.csv", index=False)
    worst, ws_ = 1.0, ""; lo = []
    for s in sorted(cs.site.unique()):
        sub, _ = R4.identifiable(cs[cs.site != s])
        if sub.episode.nunique() < 3:
            continue
        Ts, csx = R4.q3(R4.Comb(sub, **model_kw).solve(), majors); Ts = Ts.set_index("episode")
        ee = [e for e in mj.episode if e in Ts.index]
        if not ee:
            continue
        fr = float(np.mean([abs(Ts.dD[e] - Tm.dD[e]) < 2 * Ts.sig[e] for e in ee]))
        lo.append(Ts.reset_index().assign(left_out=s))
        if fr < worst:
            worst, ws_ = fr, s
    if lo:
        pd.concat(lo, ignore_index=True).to_csv(f"{OUT[0]}/r4_loso_{label}.csv", index=False)
    G.update(loso_worst=worst, loso_worst_site=ws_, G_loso_pass=bool(worst >= 0.9))
    oo, ee_ = (R4.identifiable(cs[cs.fam % 2 == k])[0] for k in (1, 0))
    To, co = R4.q3(R4.Comb(oo, **model_kw).solve(), majors); Te, ce = R4.q3(R4.Comb(ee_, **model_kw).solve(), majors)
    if co and ce:
        J_ = To.set_index("episode").join(Te.set_index("episode"), lsuffix="_o", rsuffix="_e", how="inner"); J_ = J_[J_.major_o]
        fr = float(np.mean(np.abs(J_.dD_o - J_.dD_e) < 2 * np.hypot(J_.sig_o, J_.sig_e))) if len(J_) else np.nan
        G.update(famsplit_frac=fr, famsplit_odd=co["composite"], famsplit_odd_sig=co["composite_sig"], famsplit_even=ce["composite"],
                 famsplit_even_sig=ce["composite_sig"],
                 G_fam_pass=bool(fr >= 0.9 and abs(co["composite"] - ce["composite"]) < 2 * np.hypot(co["composite_sig"], ce["composite_sig"])))
    lg(f"gates [{label}|{R.setup_name(setup) if not isinstance(setup, str) else setup}]: " + R4.fmt(G))
    T["label"] = label
    return G, T, sw


_B1 = {}


def _ctrl_b1(cs, setup):
    """same cells with the B1 reference-path control values"""
    if "b1" not in _B1:
        b = pd.concat([R4.b2_cells(False, ctrl="b1"), R4.b2_cells(True, ctrl="b1")], ignore_index=True)
        _B1["b1"] = b.set_index(["series", "episode"])[["dc", "sc"]]
    m = cs.drop(columns=["dc", "sc"]).join(_B1["b1"], on=["series", "episode"])
    return m


# ------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    lg = R4.log("combined"); t0 = time.time()
    majors = R4.majors(); lsel = R4.lapse_sel()
    H = load_H(lg)
    if a.limit:
        keep = list(dict.fromkeys(H.series))[:a.limit]; H = H[H.series.isin(keep)].reset_index(drop=True)
        OUT[0] = f"{R4.TB4}/test"; os.makedirs(OUT[0], exist_ok=True)
    H.to_csv(f"{OUT[0]}/r4_cells_H.csv", index=False)
    # ---------------- natural values: H and columns
    rows, Ts = [], []
    nr = H[~H.far].reset_index(drop=True)
    b0f = R4.b0_far_cells()
    cols = {"H": (H, {}), "near_only": (R4.identifiable(nr)[0], {}), "noR": (H, dict(Rr=False)), "noSc": (H, dict(Sc=False)),
            "prior_not_constraint": (H, dict(constraint=False, prior_ep=1.0)), "noMB": (H, dict(MB=False))}
    if len(b0f):
        nb = pd.concat([nr, b0f[[c for c in b0f.columns if c in nr.columns]]], ignore_index=True)
        cols["near_B2_far_B0"] = (R4.identifiable(nb)[0], {})
    for setup in ("ell40", lsel):
        for name, (cc, kw) in cols.items():
            if setup != "ell40" and name not in ("H", "near_only"):
                continue
            cs = R4.with_setup(cc, setup) if setup != "ell40" else cc
            cs = cs[np.isfinite(cs.a_d) & np.isfinite(cs.a_s)].reset_index(drop=True)
            L = R4.Comb(cs, **kw); F = L.solve(); T, comp = R4.q3(F, majors)
            r = R4.summary(F, T, comp, name); r.update(setup=setup, cond=L.cond, singular=L.singular, n_R_stations=len(L.stas_R))
            if len(T) and "2019-11" in set(T.episode):
                x = T.set_index("episode").loc["2019-11"]; r.update(dD_2019_11=float(x.dD), sig_2019_11=float(x.sig))
            rows.append(r); Ts.append(T.assign(label=name, setup=setup))
            lg(f"natural {name}/{setup}: " + R4.fmt(r))
    pd.DataFrame(rows).to_csv(f"{OUT[0]}/r4_models.csv", index=False)
    pd.concat(Ts, ignore_index=True).to_csv(f"{OUT[0]}/r4_dD.csv", index=False)
    # ---------------- injections
    vals = run_inj(H, majors, a.workers, lg)
    inj = vectors(vals, H)
    dn = inj["nat"]; okn = np.isfinite(dn)
    lg(f"injections measured {time.time() - t0:.0f}s: {len(inj)} runs; natural-from-store vs measured max |diff| "
       f"{np.nanmax(np.abs(dn[okn] - H.d.to_numpy()[okn])):.1e} %, missing {int((~okn).sum())} of {len(H)}")
    np.savez_compressed(f"{OUT[0]}/r4_injected.npz", runs=np.array(sorted(inj)), d=np.array([inj[k] for k in sorted(inj)]),
                        series=H.series.to_numpy(), episode=H.episode.to_numpy())
    # ---------------- gates + rules: H under NOM (headline) and L-SEL (column); in-sample rules for the columns
    allR, allO, allG, allT = [], [], [], []
    verdict = {}
    for lab, setup, gin, kw in (("H", "ell40", "ell40", {}), ("H", lsel, lsel, {}), ("HnoSc", "ell40", "ell40", dict(Sc=False)),
                                ("HnoSc", lsel, lsel, dict(Sc=False))):            # HnoSc: PLAN §8 Amendment 1
        G, T, sw = gates(H, setup, majors, lg, f"{lab}_{setup}", kw)
        allG.append(dict(model=lab, setup=setup, **G)); allT.append(T)
        TR, O = evaluate(H, inj, setup, majors, lg, lab, sw, kw)
        allR.append(TR); allO.append(O)
        v = rule_verdict(TR, O, gin)
        gp = dict(G_ctrl=G.get("ctrl_g_pass"), G_loso=G.get("G_loso_pass"), G_fam=G.get("G_fam_pass"), G_prior=G.get("G_prior_pass"),
                  G_site=G.get("G_site_pass"), acceptance=G.get("acceptance"))
        mj = T[T.major].set_index("episode")
        named = [e for e in mj.index if abs(mj.dD[e] / mj.sig[e]) >= 3]
        allpass = all(v[k] for k in ("T1", "T2", "T3", "T4", "T5") if v[k] is not None) and all(bool(gp[k]) for k in ("G_ctrl", "G_loso", "G_fam", "G_prior"))
        verdict[f"{lab}_{setup}"] = dict(rules=v, gates=gp, majors_ge_3sig=named, wording_all_pass=bool(allpass),
                                         deep_weighted_episode_changes=named if allpass else [], composite=G["composite"],
                                         composite_sig=G["composite_sig"], sig_major_median=G["sig_major_median"])
        lg(f"VERDICT {lab}_{setup}: {verdict[f'{lab}_{setup}']}")
    for name in ("near_only", "noR", "prior_not_constraint", "noMB"):
        cc, kw = cols[name]
        if name == "near_only":
            sub = H[~H.far]; inj_s = {k: v[~H.far.to_numpy()] for k, v in inj.items()}
            sub = sub.reset_index(drop=True); sub, _ = R4.identifiable(sub)
            keep = set(zip(sub.series, sub.episode)); m = np.array([(s, e) in keep for s, e in zip(H.series[~H.far], H.episode[~H.far])])
            inj_s = {k: v[m] for k, v in inj_s.items()}
        else:
            sub, inj_s = H, inj
        TR, O = evaluate(sub, inj_s, "ell40", majors, lg, name, None, kw)
        allR.append(TR[TR.gen == "ell40"]); allO.append(O)
    pd.concat(allR, ignore_index=True).to_csv(f"{OUT[0]}/r4_rules.csv", index=False)
    pd.concat(allO, ignore_index=True).to_csv(f"{OUT[0]}/r4_onesite_receiver.csv", index=False)
    pd.DataFrame(allG).to_csv(f"{OUT[0]}/r4_gates.csv", index=False)
    pd.concat(allT, ignore_index=True).to_csv(f"{OUT[0]}/r4_dD_H.csv", index=False)
    R4.dump(dict(verdict=verdict, lsel=lsel, written=R4.pt_now(), n_cells=len(H)), f"{OUT[0]}/r4_verdict.json")
    lg(f"combined done {time.time() - t0:.0f}s")
    if not a.limit:
        open(f"{R4.R4}/COMBINED_{C.REGION}_DONE", "w").write(R4.pt_now() + "\n")


if __name__ == "__main__":
    main()
