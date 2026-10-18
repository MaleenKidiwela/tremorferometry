#!/usr/bin/env python
"""s11_r7_joint -- S11 round 7, task 3 (PLAN.md §11.3): joint arrival-time (ballistic S-P, ray-traced) + coda (E4 rows)
linear inversion for shallow (0-5 km), mid-crust (5-25 km) and deep (25-45 km) dv/v per episode and per decade.
Fits: coda-only, ballistic-only, joint; columns Vs-only, + far-pair dS rows, + source-locus nuisance.
    S11_REGION=... python s11_r7_joint.py [--test N]
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
import s11_r7_common as Q
import s11_r5_common as R5
import s11_r4_common as R4

BALL = f"{Q.BASE}/ballistic/round3/{Q.SUB}/tables"
GRP = ["s", "m", "d"]
STEP = pd.Timestamp("2015-01-01")


# ------------------------------------------------------------------------------ ballistic rows
def sp_scale(rays):
    """§11 Amendment 1: f = observed S-P / traced (tS - tP) per (family, station), median over epochs (1 where missing)"""
    geo = pd.read_csv(f"{BALL}/series_geom.csv")[["fam", "netsta", "SP"]]
    m = geo.merge(rays[["fam", "netsta", "tS", "tP"]], on=["fam", "netsta"], how="inner")
    m["f"] = m.SP / (m.tS - m.tP)
    return m.groupby(["fam", "netsta"]).f.median()


def scaled_rays(rays, on=True):
    r = rays.copy()
    if not on:
        return r
    f = sp_scale(rays).reindex(list(zip(r.fam, r.netsta))).to_numpy()
    f = np.where(np.isfinite(f) & (f > 0), f, 1.0)
    for c in [c for c in r.columns if c.startswith("tS") or c.startswith("tP")]:
        r[c] = r[c] * f
    return r


def ball_rows(kind="EP", eps_=None, rays=None, vs_only=False, lg=None):
    if kind == "EP":
        T = pd.read_parquet(f"{BALL}/lags_station.parquet"); T = T[(T.bin == "EP") & (T.kind == "data")]
        T = T[T.lab.isin(eps_)]
    else:
        T = pd.read_parquet(f"{BALL}/lags_station_dh.parquet"); T = T[(T.bin == "DH") & (T.kind == "data") & T.lab.isin(["E", "L"])]
        T = T[T.netsta != "PB.B013"]
    T = T[~T.railS.astype(bool) & ~T.railP3.astype(bool) & np.isfinite(T.lagS) & np.isfinite(T.lagP3) & np.isfinite(T.sigSmP3) & (T.sigSmP3 > 0)].copy()
    T["y"] = T.lagS - T.lagP3; T["sig"] = T.sigSmP3
    T["series"] = "b:" + T.fam.astype(str) + "|" + T.netsta + "|" + T.epoch
    T = T.merge(rays, on=["fam", "netsta"], how="inner")
    for g in ["s", "m", "d", "b"]:
        T[f"K_{g}"] = -10.0 * (T[f"tS_{g}"] if vs_only else (T[f"tS_{g}"] - T[f"tP_{g}"]))
    if kind == "EP":
        T["episode"] = T.lab
        n = T.groupby("series").episode.nunique(); T = T[T.series.map(n) >= 3]
    else:
        n = T.groupby("series").lab.nunique(); T = T[T.series.map(n) == 2]
    if lg:
        lg(f"ballistic {kind} rows: {len(T)} ({T.series.nunique()} series, {T.fam.nunique()} families, {T.netsta.nunique()} stations)")
    return T.reset_index(drop=True)


def pair_rows(eps_, rays, lg=None):
    P = pd.read_parquet(f"{BALL}/lags_pair.parquet"); P = P[(P.bin == "EP") & (P.kind == "data") & P.lab.isin(eps_)]
    site = C.SITE
    P = P[[site.get(a, a) != site.get(b, b) for a, b in zip(P.A, P.B)]]
    P = P[~P.railS.astype(bool) & np.isfinite(P.dS) & np.isfinite(P.sig_dS) & (P.sig_dS > 0)].copy()
    ra = rays.set_index(["fam", "netsta"])
    ka = ra.reindex(list(zip(P.fam, P.A))); kb = ra.reindex(list(zip(P.fam, P.B)))
    for g in ["s", "m", "d", "b"]:
        P[f"K_{g}"] = -10.0 * (kb[f"tS_{g}"].to_numpy() - ka[f"tS_{g}"].to_numpy())
    P = P[np.isfinite(P.K_d)]
    P["y"] = P.dS; P["sig"] = P.sig_dS; P["episode"] = P.lab
    P["series"] = "p:" + P.fam.astype(str) + "|" + P.pair + "|" + P.epA + "|" + P.epB
    n = P.groupby("series").episode.nunique(); P = P[P.series.map(n) >= 3]
    if lg:
        lg(f"far-pair dS rows: {len(P)} ({P.series.nunique()} series, pairs {sorted(P.pair.unique())[:8]}...)")
    return P.reset_index(drop=True)


# ------------------------------------------------------------------------------ design
def episode_design(cs, B, eps_, use_coda=True, use_ball=True, locus=False):
    """rows: coda cells (E4 + S(e), M(e), D(e)) and/or ballistic rows; shared S/M/D episode unknowns"""
    blocks_G, series, w, y, rowtype = [], [], [], [], []
    names = [(g, "", e) for g in GRP for e in eps_]
    ie = {e: i for i, e in enumerate(eps_)}
    extra_names = []; extra_cols = []
    nrows = 0
    parts = []
    if use_coda:
        c = cs.reset_index(drop=True)
        cols = [(g, "", c[f"a_{g}"].to_numpy(float)) for g in GRP]
        Gc, nmc, _ = Q.e4_blocks(c, cols)
        # e4_blocks orders episode columns by the episodes of c; map to eps_
        epc = Q.episodes_of(c)
        assert epc == [e for e in eps_ if e in set(epc)]
        P = len(names)
        Gfull = np.zeros((len(c), P))
        k = 0
        for g in GRP:
            for e in epc:
                Gfull[:, names.index((g, "", e))] = Gc[:, k]; k += 1
        rest = Gc[:, k:]; rest_names = nmc[k:]
        parts.append(("coda", Gfull, rest, rest_names, c.series.to_numpy(), 1 / c.sig.to_numpy() ** 2, c.d.to_numpy(float)))
    if use_ball and B is not None and len(B):
        b = B.reset_index(drop=True)
        P = len(names); Gb = np.zeros((len(b), P))
        rr = np.arange(len(b))
        for g in GRP:
            Gb[rr, [names.index((g, "", e)) for e in b.episode]] = b[f"K_{g}"].to_numpy()
        mid = b.episode.map(Q.U.mid); t, doy = V.tdec(mid); bR = V.basis(t, doy, R5.RT)
        restb = np.column_stack([b.K_b.to_numpy() * bR[tm] for tm in R5.RT]); rnb = [("B", "", tm) for tm in R5.RT]
        if locus:
            sp = b.series.str.startswith("b:").to_numpy()
            key = b.fam.astype(str) + "#" + b.episode
            cnt = b[sp].groupby(key[sp]).netsta.nunique()
            lk = sorted(cnt.index[cnt >= 2])
            Lm = np.zeros((len(b), len(lk))); li = {k_: i for i, k_ in enumerate(lk)}
            for i, (k_, s_) in enumerate(zip(key, sp)):
                if s_ and k_ in li:
                    Lm[i, li[k_]] = 1.0
            restb = np.column_stack([restb, Lm]); rnb = rnb + [("L", k_, "") for k_ in lk]
        parts.append(("ball", Gb, restb, rnb, b.series.to_numpy(), 1 / b.sig.to_numpy() ** 2, b.y.to_numpy(float)))
    # merge the non-episode columns by name
    allrest = []
    for p in parts:
        for nm in p[3]:
            if nm not in allrest:
                allrest.append(nm)
    ri = {nm: i for i, nm in enumerate(allrest)}
    Gs, S_, W_, Y_, T_ = [], [], [], [], []
    for kind, Ge, Rr, rn, ser, ww, yy in parts:
        Rm = np.zeros((len(Ge), len(allrest)))
        for j, nm in enumerate(rn):
            Rm[:, ri[nm]] += Rr[:, j]
        Gs.append(np.column_stack([Ge, Rm])); S_.append(ser); W_.append(ww); Y_.append(yy); T_ += [kind] * len(Ge)
    G = np.vstack(Gs); allnames = names + allrest
    idx = {nm: i for i, nm in enumerate(allnames)}; Pn = len(allnames)
    # keep only episodes with data in at least one row
    A = []
    used = np.abs(G).sum(0) > 0
    eps_used = [e for e in eps_ if any(used[idx[(g, "", e)]] for g in GRP)]
    for g in GRP:
        ee = [e for e in eps_used if used[idx[(g, "", e)]]]
        if len(ee) >= 2:
            A += Q.proj_rows(idx, Pn, g, "", ee)
        for e in eps_:
            if not used[idx[(g, "", e)]]:
                a = np.zeros(Pn); a[idx[(g, "", e)]] = 1.0; A.append(a)
    L = Q.LinJ(G, allnames, np.concatenate(S_), np.concatenate(W_), A)
    L.eps = eps_used; L.rowtype = np.array(T_)
    return L, np.concatenate(Y_)


def decade_design(cs, B, use_coda=True, use_ball=True):
    parts = []
    names = [("X", g, "step") for g in ["s", "m", "d", "b"]]
    if use_coda:
        c = cs.reset_index(drop=True)
        st = (c.mid >= STEP).to_numpy().astype(float)
        Ge = np.column_stack([c[f"a_{g}"].to_numpy(float) * st for g in ["s", "m", "d", "b"]])
        Gu, nmu, _ = Q.e4_blocks(c, [], units=True, MB=False, uterms=["c1", "s1", "c2", "s2"], usmall=["c1", "s1"])
        parts.append((Ge, Gu, nmu, c.series.to_numpy(), 1 / c.sig.to_numpy() ** 2, c.d.to_numpy(float)))
    if use_ball and B is not None and len(B):
        b = B.reset_index(drop=True); st = (b.lab == "L").to_numpy().astype(float)
        Ge = np.column_stack([b[f"K_{g}"].to_numpy() * st for g in ["s", "m", "d", "b"]])
        parts.append((Ge, np.zeros((len(b), 0)), [], b.series.to_numpy(), 1 / b.sig.to_numpy() ** 2, b.y.to_numpy(float)))
    allrest = []
    for p in parts:
        for nm in p[2]:
            if nm not in allrest:
                allrest.append(nm)
    ri = {nm: i for i, nm in enumerate(allrest)}
    Gs, S_, W_, Y_ = [], [], [], []
    for Ge, Rr, rn, ser, ww, yy in parts:
        Rm = np.zeros((len(Ge), len(allrest)))
        for j, nm in enumerate(rn):
            Rm[:, ri[nm]] += Rr[:, j]
        Gs.append(np.column_stack([Ge, Rm])); S_.append(ser); W_.append(ww); Y_.append(yy)
    G = np.vstack(Gs); allnames = names + allrest
    keep = np.abs(G).sum(0) > 0
    G = G[:, keep]; allnames = [n for n, k in zip(allnames, keep) if k]
    L = Q.LinJ(G, allnames, np.concatenate(S_), np.concatenate(W_), None)
    return L, np.concatenate(Y_)


# ------------------------------------------------------------------------------ summaries
def ep_summary(L, F, majors, label):
    rows = []
    for g in GRP:
        gs = {}
        for e in majors:
            v = Q.dvec_g(F, L.eps, e, g)
            if v is not None:
                gs[e] = v
        for e, v in gs.items():
            val, sg = Q.val_sig(F, v); rows.append(dict(fit=label, group=g, episode=e, value=val, sig=sg))
    T = pd.DataFrame(rows)
    S = []
    for g in GRP:
        t = T[T.group == g] if len(T) else T
        if not len(t):
            continue
        gc = np.mean([Q.dvec_g(F, L.eps, e, g) for e in t.episode], 0)
        cv, cs_ = Q.val_sig(F, gc)
        S.append(dict(fit=label, group=g, n_majors=len(t), sig_major_med=float(t.sig.median()), composite=cv, composite_sig=cs_,
                      n_3sig=int((np.abs(t.value) >= 3 * t.sig).sum()), chi2_majors=float(np.mean((t.value / t.sig) ** 2))))
    return T, pd.DataFrame(S)


def inject_ep(L, rows_info, majors, F0):
    """model-consistent pulses through the operator: D +0.3 % at each major; S +0.1 % at each major"""
    out = []
    for kind, gsrc, amp in (("Dpulse", "d", 0.3), ("Spulse", "s", 0.1)):
        rec = {g: [] for g in GRP}
        for e in majors:
            if Q.dvec_g(F0, L.eps, e, "d") is None:
                continue
            x = rows_info["coef"][gsrc] * amp * (rows_info["episode"] == e)
            th = L.Op @ x
            Fx = dict(F0, th=th)
            for g in GRP:
                v = Q.dvec_g(F0, L.eps, e, g)
                rec[g].append(float(v @ th) / amp)
        out.append(dict(pulse=kind, **{f"to_{g}_med": float(np.median(rec[g])) if rec[g] else np.nan for g in GRP},
                        **{f"to_{g}_absmax": float(np.max(np.abs(rec[g]))) if rec[g] else np.nan for g in GRP}))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--test", type=int, default=0)
    a = ap.parse_args(); test = a.test > 0
    TB = Q.tbdir(test); lg = Q.log("joint_test" if test else "joint"); t0 = time.time()
    H, inj = R5.load_H(lg, limit=a.test)
    majors = Q.majors()
    cs0, good = R5.setup_cells(H, Q.LSEL); good = good & np.isfinite(inj["nat"])
    cs = cs0[good].reset_index(drop=True)
    eps_ = Q.episodes_of(cs)
    # rays
    LB = pd.read_parquet(f"{BALL}/lags_station.parquet")[["fam", "netsta"]].drop_duplicates()
    LP = pd.read_parquet(f"{BALL}/lags_pair.parquet")[["fam", "A", "B"]]
    pairs = pd.concat([LB, LP.rename(columns={"A": "netsta"})[["fam", "netsta"]], LP.rename(columns={"B": "netsta"})[["fam", "netsta"]]]).drop_duplicates()
    rays = Q.ray_table(pairs, lg)
    rays.to_csv(f"{TB}/r7_joint_rays.csv", index=False)
    geo = pd.read_csv(f"{BALL}/series_geom.csv")
    m = geo.merge(rays, on=["fam", "netsta"], how="inner", suffixes=("_obs", ""))
    if len(m):
        rr = (m.tS - m.tP) / m.SP
        lg(f"ray check: traced tS−tP / observed SP median {np.median(rr):.3f} (q10 {np.quantile(rr, 0.1):.3f}, q90 {np.quantile(rr, 0.9):.3f}) over {len(m)} series")
    raw_rays = rays
    rays = scaled_rays(raw_rays, True)                     # §11 Amendment 1 (post-code-test): kernel total = observed S-P
    fsc = sp_scale(raw_rays); lg(f"S-P kernel scale f (observed / traced): median {fsc.median():.3f}, q10 {fsc.quantile(0.1):.3f}, q90 {fsc.quantile(0.9):.3f}")
    BA = ball_rows("EP", eps_, rays, lg=lg); BV = ball_rows("EP", eps_, rays, vs_only=True)
    BR = ball_rows("EP", eps_, raw_rays)
    PR = pair_rows(eps_, rays, lg)
    if test:                                    # code test: ballistic rows of the subset's families only
        fs = set(H.fam.unique()); BA = BA[BA.fam.isin(fs)]; BV = BV[BV.fam.isin(fs)]; PR = PR[PR.fam.isin(fs)]; BR = BR[BR.fam.isin(fs)]
    # ballistic chi2 factor (ballistic-only episode fit, constant Vp/Vs)
    Lb, yb = episode_design(cs, BA, eps_, use_coda=False)
    Fb = Lb.solve(yb, scale_chi2=False); fac = float(np.sqrt(max(1.0, Fb["chi2"])))
    Lp_, yp_ = episode_design(cs, PR, eps_, use_coda=False) if len(PR) else (None, None)
    facp = float(np.sqrt(max(1.0, Lp_.solve(yp_, scale_chi2=False)["chi2"]))) if Lp_ is not None else 1.0
    lg(f"ballistic error factors: S−P {fac:.2f}, dS {facp:.2f}")
    BA = BA.assign(sig=BA.sig * fac); BV = BV.assign(sig=BV.sig * fac); PR = PR.assign(sig=PR.sig * facp); BR = BR.assign(sig=BR.sig * fac)
    fits = [("coda", dict(use_coda=True, use_ball=False), BA), ("ballistic", dict(use_coda=False, use_ball=True), BA),
            ("joint", dict(use_coda=True, use_ball=True), BA), ("ballistic_Vs", dict(use_coda=False, use_ball=True), BV),
            ("joint_Vs", dict(use_coda=True, use_ball=True), BV),
            ("ballistic_dS", dict(use_coda=False, use_ball=True), pd.concat([BA, PR], ignore_index=True)),
            ("joint_dS", dict(use_coda=True, use_ball=True), pd.concat([BA, PR], ignore_index=True)),
            ("ballistic_locus", dict(use_coda=False, use_ball=True, locus=True), BA),
            ("joint_locus", dict(use_coda=True, use_ball=True, locus=True), BA),
            ("ballistic_raw", dict(use_coda=False, use_ball=True), BR), ("joint_raw", dict(use_coda=True, use_ball=True), BR)]
    TE, SU, IN = [], [], []
    # E4 reference (D only)
    LE = Q.e4_model(cs, [("D", "", cs.a_d.to_numpy())]); FE = LE.solve(cs.d.to_numpy())
    rows = []
    for e in majors:
        v = Q.dvec_g(FE, LE.eps, e, "D")
        if v is not None:
            rows.append(dict(fit="E4_Donly", group="d", episode=e, value=Q.val_sig(FE, v)[0], sig=Q.val_sig(FE, v)[1]))
    t = pd.DataFrame(rows); TE.append(t)
    gc = np.mean([Q.dvec_g(FE, LE.eps, e, "D") for e in t.episode], 0)
    SU.append(pd.DataFrame([dict(fit="E4_Donly", group="d", n_majors=len(t), sig_major_med=float(t.sig.median()),
                                 composite=Q.val_sig(FE, gc)[0], composite_sig=Q.val_sig(FE, gc)[1])]))
    for label, kw, Bx in fits:
        try:
            L, y = episode_design(cs, Bx, eps_, **kw)
        except Exception as ex:
            lg(f"[{label}] failed: {ex}"); continue
        F = L.solve(y)
        T, S = ep_summary(L, F, majors, label)
        if len(S):
            S["chi2_fit"] = F["chi2"]; S["n_rows"] = len(y); S["n_coda"] = int((L.rowtype == "coda").sum()); S["n_ball"] = int((L.rowtype == "ball").sum())
        TE.append(T); SU.append(S)
        # injections
        n_c = int((L.rowtype == "coda").sum())
        coef = {g: np.zeros(len(y)) for g in GRP}; epis = np.empty(len(y), dtype=object)
        if kw.get("use_coda"):
            for g in GRP:
                coef[g][:n_c] = cs[f"a_{g}"].to_numpy(float)
            epis[:n_c] = cs.episode.to_numpy()
        if kw.get("use_ball"):
            for g in GRP:
                coef[g][n_c:] = Bx[f"K_{g}"].to_numpy()
            epis[n_c:] = Bx.episode.to_numpy()
        for r in inject_ep(L, dict(coef=coef, episode=epis), majors, F):
            r["fit"] = label; IN.append(r)
        lg(f"[{label}] rows {len(y)} chi2 {F['chi2']:.2f}: " + "; ".join(f"{r.group}: σ/major {r.sig_major_med:.4f}, composite {r.composite:+.4f} ± {r.composite_sig:.4f}"
                                                                  for r in S.itertuples()))
    TE = pd.concat(TE, ignore_index=True); SU = pd.concat(SU, ignore_index=True)
    TE.to_csv(f"{TB}/r7_joint_episode.csv", index=False); SU.to_csv(f"{TB}/r7_joint_summary.csv", index=False)
    pd.DataFrame(IN).to_csv(f"{TB}/r7_joint_inject.csv", index=False)
    # decade
    BD = ball_rows("DH", None, rays, lg=lg)
    if test:
        BD = BD[BD.fam.isin(set(H.fam.unique()))]
    Lbd, ybd = decade_design(cs, BD, use_coda=False)
    facd = float(np.sqrt(max(1.0, Lbd.solve(ybd, scale_chi2=False)["chi2"]))); BD = BD.assign(sig=BD.sig * facd)
    lg(f"ballistic decade error factor {facd:.2f}")
    DR = []
    for label, kw in (("coda", dict(use_coda=True, use_ball=False)), ("ballistic", dict(use_coda=False, use_ball=True)),
                      ("joint", dict(use_coda=True, use_ball=True))):
        L, y = decade_design(cs, BD, **kw); F = L.solve(y)
        for g in ["s", "m", "d", "b"]:
            nm = ("X", g, "step")
            if nm not in F["idx"]:
                continue
            gv = np.zeros(len(F["th"])); gv[F["idx"][nm]] = 1.0; v, s = Q.val_sig(F, gv)
            # injection: +0.1 % deep step
            n_c = len(cs) if kw["use_coda"] else 0
            x = np.zeros(len(y))
            if kw["use_coda"]:
                x[:n_c] = cs.a_d.to_numpy() * 0.1 * (cs.mid >= STEP).to_numpy()
            if kw["use_ball"]:
                x[n_c:] = BD.K_d.to_numpy() * 0.1 * (BD.lab == "L").to_numpy()
            th = L.Op @ x
            DR.append(dict(fit=label, group=g, value=v, sig=s, bound2=abs(v) + 2 * s, chi2_fit=F["chi2"], n_rows=len(y),
                           deep_step_to_this=float(gv @ th) / 0.1))
        lg(f"[decade {label}] " + "; ".join(f"{r['group']} {r['value']:+.4f} ± {r['sig']:.4f}" for r in DR if r["fit"] == label))
    DR = pd.DataFrame(DR); DR.to_csv(f"{TB}/r7_joint_decade.csv", index=False)
    # reading
    def get(fit, g, col):
        t = SU[(SU.fit == fit) & (SU.group == g)]
        return float(t.iloc[0][col]) if len(t) else np.nan
    V_ = dict(region=C.REGION, ballistic_factor=fac, decade_factor=facd, written=Q.pt_now())
    for g in GRP:
        sc, sb, sj = get("coda", g, "sig_major_med"), get("ballistic", g, "sig_major_med"), get("joint", g, "sig_major_med")
        cc, cb, cj = (get(f, g, "composite") for f in ("coda", "ballistic", "joint"))
        ccs, cbs, cjs = (get(f, g, "composite_sig") for f in ("coda", "ballistic", "joint"))
        V_[f"{g}_sig_major"] = dict(coda=sc, ballistic=sb, joint=sj, joint_over_coda=sj / sc, joint_over_ball=sj / sb,
                                    ball_info_share=1 - (sj / sc) ** 2 if np.isfinite(sc) else np.nan)
        V_[f"{g}_composite"] = dict(coda=[cc, ccs], ballistic=[cb, cbs], joint=[cj, cjs],
                                    detected=bool(abs(cj) >= 3 * cjs and np.sign(cc) == np.sign(cb)), bound2=abs(cj) + 2 * cjs,
                                    ballistic_adds=bool(cjs <= 0.8 * ccs), coda_adds=bool(cjs <= 0.8 * cbs))
        d = DR[DR.group == g].set_index("fit")
        if {"coda", "ballistic", "joint"} <= set(d.index):
            V_[f"{g}_decade"] = dict(coda=[d.loc["coda", "value"], d.loc["coda", "sig"]], ballistic=[d.loc["ballistic", "value"], d.loc["ballistic", "sig"]],
                                     joint=[d.loc["joint", "value"], d.loc["joint", "sig"]],
                                     detected=bool(abs(d.loc["joint", "value"]) >= 3 * d.loc["joint", "sig"] and np.sign(d.loc["coda", "value"]) == np.sign(d.loc["ballistic", "value"])),
                                     bound2=d.loc["joint", "bound2"], ballistic_adds=bool(d.loc["joint", "sig"] <= 0.8 * d.loc["coda", "sig"]),
                                     coda_adds=bool(d.loc["joint", "sig"] <= 0.8 * d.loc["ballistic", "sig"]))
    Q.dump(V_, f"{TB}/r7_joint_verdict.json")
    lg("VERDICT " + json.dumps(V_, default=str))
    lg(f"done {time.time() - t0:.0f}s")
    if not test:
        open(f"{Q.R7}/JOINT_{C.REGION}_DONE", "w").write(Q.pt_now() + "\n")


if __name__ == "__main__":
    main()
