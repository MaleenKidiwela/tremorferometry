#!/usr/bin/env python
"""s11_ball3_phase -- S11 ballistic round 3, §2 ETS-phase composite, both regions (round3/PLAN.md + Amendment 2).
PH labels BTW (bursts) / ARR (<= 1 d after the family onset) / AFT (> 1 d). Per series contrasts C1 = ARR - BTW,
C2 = AFT - ARR, C3 = AFT - BTW; items: S - P3 per station, dS per (reference, station) (uncorrected headline, corrected
column), absolute S and P3 per station (columns). Pooling 1/sigma^2 with chi2 inflation; region-wide in dv/v units with the
random-effects error across stations. Controls: within-series label permutation (10,000), episode-date permutation (PHp,
200 draws), scrambled control (c2). PH5 / PHE family-common curves as columns.

    python s11_ball3_phase.py [--tag _test]
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_ball3_lib as L

NPERM_LAB = 10000
SEED = 20261013
LABS = ["BTW", "ARR", "AFT"]
CONTR = {"C1": ("ARR", "BTW"), "C2": ("AFT", "ARR"), "C3": ("AFT", "BTW")}
PH_REFS = {"SVI": {"BH": ["PB.B009", "PB.B010", "PB.B011"]}, "WA_NW": {"G1BH": ["PB.B005", "PB.B006", "PB.B007"],
                                                                       "G2": ["PB.B013", "PB.B943"]}}
PERMSAVE = {}
COLOC_X = {"SVI": {"BH": ["CN.PGC"]}, "WA_NW": {"G1BH": ["UW.LRIV"], "G2": []}}


# ------------------------------------------------------------------------------------------ series matrices
def to_matrix(df, keys, val, sig, labs=LABS):
    """df rows (keys..., lab, val, sig) -> series frame S (keys), V (n,3), E (n,3) variances"""
    d = df[df.lab.isin(labs) & np.isfinite(df[val]) & np.isfinite(df[sig])]
    piv = d.pivot_table(index=keys, columns="lab", values=[val, sig], aggfunc="first")
    V = np.column_stack([piv[(val, l)].to_numpy() if (val, l) in piv else np.full(len(piv), np.nan) for l in labs])
    E = np.column_stack([piv[(sig, l)].to_numpy() if (sig, l) in piv else np.full(len(piv), np.nan) for l in labs]) ** 2
    S = piv.index.to_frame(index=False)
    keep = np.isfinite(V).sum(1) >= 2
    return S[keep].reset_index(drop=True), V[keep], E[keep]


def contrasts(V, E):
    out = {}
    for c, (a, b) in CONTR.items():
        ia, ib = LABS.index(a), LABS.index(b)
        out[c] = (V[..., ia] - V[..., ib], E[..., ia] + E[..., ib])
    return out


def pooled(C, Vr, scale=None):
    """vectorised 1/var mean along last axis; C, Vr (..., n) -> mean (...)"""
    sc = 1.0 if scale is None else scale
    x = C * sc; v = Vr * sc ** 2
    w = np.where(np.isfinite(x) & np.isfinite(v) & (v > 0), 1 / np.where(v > 0, v, 1), 0.0)
    return (np.nan_to_num(x) * w).sum(-1) / np.maximum(w.sum(-1), 1e-30)


def label_perm(V, E, rng, n=NPERM_LAB, scale=None):
    """within-series permutation of (value, variance) among present labels -> dict contrast -> array (n,) pooled"""
    pres = np.isfinite(V)
    Ppos = np.argsort(~pres, axis=1, kind="stable")
    res = {c: np.empty(n) for c in CONTR}
    r_ = np.arange(len(V))[:, None]
    for k in range(n):
        keys = rng.random(V.shape); keys[~pres] = 2.0
        order = np.argsort(keys, axis=1)
        Vp = np.full_like(V, np.nan); Ep = np.full_like(E, np.nan)
        Vp[r_, Ppos] = V[r_, order]; Ep[r_, Ppos] = E[r_, order]
        for c, (C, Vr) in contrasts(Vp, Ep).items():
            res[c][k] = pooled(C, Vr, scale)
    return res


def summarize(C, Vr, groups=None, scale=None):
    sc = np.ones(len(C)) if scale is None else scale
    ok = np.isfinite(C) & np.isfinite(Vr) & np.isfinite(sc)
    return L.pool_re((C * sc)[ok], (np.sqrt(Vr) * np.abs(sc))[ok], None if groups is None else np.asarray(groups)[ok])


# ------------------------------------------------------------------------------------------ region items
def region_items(reg, tag):
    _, TB, _, _ = L.dirs(reg)
    st = L.load(reg, "lags_station", tag); pr = L.load(reg, "lags_pair", tag)
    SP = L.sp_table(reg, tag)
    items = {}
    for kind in ("data", "ctrl"):
        s = st[(st.kind == kind) & (st.bin == "PH")].copy()
        s["SmP"] = np.where(s.railS | s.railP3, np.nan, s.lagS - s.lagP3)
        s["sSmP"] = s.sigSmP3
        s["aS"] = np.where(s.railS, np.nan, s.lagS); s["aP"] = np.where(s.railP3, np.nan, s.lagP3)
        for obs, val, sig in (("S-P3", "SmP", "sSmP"), ("absS", "aS", "sigS"), ("absP3", "aP", "sigP3")):
            Sx, V, E = to_matrix(s, ["fam", "netsta", "epoch"], val, sig)
            Sx["X"] = Sx.netsta; Sx["SP"] = [SP.get((f, n), np.nan) for f, n in zip(Sx.fam, Sx.netsta)]
            items[(obs, "-", kind)] = (Sx, V, E)
        # pairs dS: reference frames (uncorrected headline; Amendment 2)
        p = pr[(pr.kind == kind) & (pr.bin == "PH")].copy()
        p.loc[p.railS, "dS"] = np.nan; p["sig"] = np.maximum(p.sig_dS, 0.5); p["r"] = p.dS
        p = p[np.isfinite(p.r)]
        if not len(p):
            continue
        Lg = L.oriented(p, cols=("r",))
        for rname, rs in PH_REFS[reg].items():
            F = L.ref_frame(Lg, rs, "r")
            Sx, V, E = to_matrix(F.assign(sigv=F.sig), ["fam", "X", "epX"], "v", "sigv")
            Sx["SP"] = [SP.get((f, n), np.nan) for f, n in zip(Sx.fam, Sx.X)]
            items[("dS", rname, kind)] = (Sx, V, E)
    # corrected dS column (data): kcorr from pair_corrected where present
    pcp = f"{TB}/pair_corrected{tag}.parquet"
    if os.path.exists(pcp):
        Dc = pd.read_parquet(pcp); Dc = Dc[Dc.bin == "PH"]
        p = pr[(pr.kind == "data") & (pr.bin == "PH")].copy(); p.loc[p.railS, "dS"] = np.nan
        p = p.merge(Dc[["fam", "pair", "epA", "epB", "lab", "kcorr"]], on=["fam", "pair", "epA", "epB", "lab"], how="left")
        p["r"] = p.dS + p.kcorr.fillna(0.0); p["sig"] = np.maximum(p.sig_dS, 0.5); p = p[np.isfinite(p.r)]
        Lg = L.oriented(p, cols=("r",))
        for rname, rs in PH_REFS[reg].items():
            F = L.ref_frame(Lg, rs, "r")
            Sx, V, E = to_matrix(F.assign(sigv=F.sig), ["fam", "X", "epX"], "v", "sigv")
            Sx["SP"] = [SP.get((f, n), np.nan) for f, n in zip(Sx.fam, Sx.X)]
            items[("dScorr", rname, "data")] = (Sx, V, E)
    return items


def perm_matrices(reg, tag, items):
    """date-permutation draws: for each item, (NPERM, n_series, 3) value matrices aligned to the data series (BTW real)"""
    _, TB, _, _ = L.dirs(reg)
    ps = pd.read_parquet(f"{TB}/perm_station{tag}.parquet"); pp = pd.read_parquet(f"{TB}/perm_pair{tag}.parquet")
    ndraw = int(max(ps.draw.max(), pp.draw.max() if len(pp) else 0)) + 1
    out = {}
    ps["SmP"] = np.where(ps.railS | ps.railP3, np.nan, ps.lagS - ps.lagP3)
    ps["aS"] = np.where(ps.railS, np.nan, ps.lagS)
    for obs, val in (("S-P3", "SmP"), ("absS", "aS")):
        Sx, V, E = items[(obs, "-", "data")]
        key = pd.MultiIndex.from_frame(Sx[["fam", "netsta", "epoch"]])
        M = np.repeat(V[None], ndraw, 0).copy(); M[:, :, 1:] = np.nan
        for code in (1, 2):
            q = ps[ps.code == code]
            ix = key.get_indexer(pd.MultiIndex.from_frame(q[["fam", "netsta", "epoch"]]))
            ok = ix >= 0
            M[q.draw.to_numpy()[ok], ix[ok], code] = q[val].to_numpy()[ok]
        out[(obs, "-")] = M
    if len(pp):
        pp = pp.copy(); pp["A"] = pp.pair.str.split("-").str[0]; pp["B"] = pp.pair.str.split("-").str[1]
        pp.loc[pp.railS, "dS"] = np.nan
        for rname, rs in PH_REFS[reg].items():
            Sx, V, E = items[("dS", rname, "data")]
            a = pp[pp.A.isin(rs) & ~pp.B.isin(rs)][["fam", "B", "epB", "draw", "code", "dS"]].rename(columns={"B": "X", "epB": "epX"})
            b = pp[pp.B.isin(rs) & ~pp.A.isin(rs)][["fam", "A", "epA", "draw", "code", "dS"]].rename(columns={"A": "X", "epA": "epX"})
            b["dS"] = -b.dS
            q = pd.concat([a, b]).groupby(["fam", "X", "epX", "draw", "code"]).dS.mean().reset_index()
            key = pd.MultiIndex.from_frame(Sx[["fam", "X", "epX"]])
            M = np.repeat(V[None], ndraw, 0).copy(); M[:, :, 1:] = np.nan
            ix = key.get_indexer(pd.MultiIndex.from_frame(q[["fam", "X", "epX"]])); ok = ix >= 0
            M[q.draw.to_numpy()[ok], ix[ok], q.code.to_numpy()[ok]] = q.dS.to_numpy()[ok]
            out[("dS", rname)] = M
    return out, ndraw


def scale_for(obs, Sx, unit):
    """per-series multiplier converting a contrast in ms to percent dv/v"""
    sp = Sx.SP.to_numpy(dtype=float)
    if obs == "S-P3":
        T = {"constVpVs": sp, "VsOnly": L.TS(sp), "VpOnly": -L.TP(sp)}[unit]
    elif obs in ("dS", "dScorr"):
        T = L.TS(sp)
    else:
        T = L.TS(sp)
    return -100.0 / (T * 1e3)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default=""); ap.add_argument("--regions", default="SVI,WA_NW")
    a = ap.parse_args(); tag = a.tag
    rng = np.random.default_rng(SEED)
    rows, curves = [], []
    for reg in a.regions.split(","):
        items = region_items(reg, tag)
        PM, ndraw = perm_matrices(reg, tag, items)
        for (obs, rname, kind), (Sx, V, E) in items.items():
            CT = contrasts(V, E)
            coloc = Sx.X.isin(COLOC_X[reg].get(rname, [])).to_numpy() if obs in ("dS", "dScorr") else np.zeros(len(Sx), bool)
            units = ["ms"] + (["constVpVs", "VsOnly", "VpOnly"] if obs == "S-P3" else ["dVs_Xleg"] if obs in ("dS", "dScorr") else [])
            # per station (X) and region-wide (all non-co-located X)
            groups_ = [("X", x) for x in sorted(Sx.X.unique())] + [("REGION", "all")]
            for gk, gv in groups_:
                m = (Sx.X == gv).to_numpy() if gk == "X" else ~coloc
                if m.sum() < 3:
                    continue
                for unit in units:
                    if gk == "X" and unit not in ("ms", "constVpVs", "dVs_Xleg"):
                        continue
                    sc = None if unit == "ms" else scale_for(obs, Sx, unit)[m]
                    r = dict(region=reg, obs=obs, ref=rname, kind=kind, X=gv, unit=unit, n_series=int(m.sum()),
                             n_fam=int(Sx.fam[m].nunique()))
                    for c, (C, Vr) in CT.items():
                        sm = summarize(C[m], Vr[m], Sx.X[m].to_numpy() if gk == "REGION" else None, sc)
                        r.update({f"{c}": sm["mean"], f"{c}_sig": sm.get("sig"), f"{c}_chi2dof": sm.get("chi2dof"),
                                  f"{c}_sig_re": sm.get("sig_re"), f"{c}_n": sm["n"]})
                    # label permutation (pooled mean statistic, same scaling)
                    nl = NPERM_LAB if (gk == "REGION" or obs == "S-P3") else 2000
                    lp = label_perm(V[m], E[m], rng, nl, sc)
                    for c in CONTR:
                        obsv = r[c]
                        r[f"{c}_labp"] = float((1 + (np.abs(lp[c]) >= abs(obsv)).sum()) / (nl + 1)) if np.isfinite(obsv) else np.nan
                    # date permutation (data only, headline items)
                    if kind == "data" and (obs, rname) in PM:
                        Mm = PM[(obs, rname)][:, m]                # (ndraw, ns, 3)
                        Ed = np.repeat(E[m][None], ndraw, 0)
                        for c, (aa, bb) in CONTR.items():
                            ia, ib = LABS.index(aa), LABS.index(bb)
                            Cd = Mm[..., ia] - Mm[..., ib]; Vd = Ed[..., ia] + Ed[..., ib]
                            stat = pooled(Cd, Vd, sc)
                            nok = np.isfinite(Cd).sum(1)
                            stat = stat[nok >= 3]
                            r[f"{c}_datep"] = float((1 + (np.abs(stat) >= abs(r[c])).sum()) / (len(stat) + 1)) if np.isfinite(r[c]) else np.nan
                            r[f"{c}_date_sd"] = float(np.std(stat)) if len(stat) else np.nan
                            r[f"{c}_date_mean"] = float(np.mean(stat)) if len(stat) else np.nan
                            if gk == "REGION":
                                PERMSAVE[f"{reg}|{obs}|{rname}|{unit}|{c}|date"] = stat
                    if gk == "REGION":
                        for c in CONTR:
                            PERMSAVE[f"{reg}|{obs}|{rname}|{unit}|{c}|{kind}|label"] = lp[c]
                    rows.append(r)
            # family-common phase curve (demeaned per series), region-wide, ms and constVpVs
            for unit in (["ms", "constVpVs"] if obs == "S-P3" else ["ms"]):
                m = ~coloc
                sc = np.ones(m.sum()) if unit == "ms" else scale_for(obs, Sx, unit)[m]
                Vm = V[m] * sc[:, None]; Em = E[m] * sc[:, None] ** 2
                Vd = Vm - np.nanmean(Vm, 1, keepdims=True)
                for j, lab in enumerate(LABS):
                    x = Vd[:, j]; w = 1 / Em[:, j]; ok = np.isfinite(x) & np.isfinite(w)
                    if ok.sum():
                        curves.append(dict(region=reg, obs=obs, ref=rname, kind=kind, scheme="PH", unit=unit, lab=lab,
                                           value=float((x[ok] * w[ok]).sum() / w[ok].sum()), n=int(ok.sum())))
        # PH5 / PHE curves (data and control): S-P3 and dS (uncorrected), region-wide, ms
        st = L.load(reg, "lags_station", tag); SPt = L.sp_table(reg, tag)
        for scheme, labs in (("PH5", ["BPRE", "BMID", "BPOST", "EARLY", "ARR0", "AFT1", "AFT3"]), ("PHE", ["BTW", "E1", "E2", "E3"])):
            for kind in ("data", "ctrl"):
                s = st[(st.kind == kind) & (st.bin == scheme)].copy()
                s["SmP"] = np.where(s.railS | s.railP3, np.nan, s.lagS - s.lagP3)
                Sx, V, E = to_matrix(s, ["fam", "netsta", "epoch"], "SmP", "sigSmP3", labs)
                if not len(Sx):
                    continue
                sp = np.array([SPt.get((f, n), np.nan) for f, n in zip(Sx.fam, Sx.netsta)])
                for unit in ("ms", "constVpVs"):
                    sc = np.ones(len(sp)) if unit == "ms" else -100.0 / (sp * 1e3)
                    Vm = V * sc[:, None]; Em = E * sc[:, None] ** 2
                    Vd = Vm - np.nanmean(Vm, 1, keepdims=True)
                    for j, lab in enumerate(labs):
                        x = Vd[:, j]; w = 1 / Em[:, j]; ok = np.isfinite(x) & np.isfinite(w)
                        if ok.sum():
                            curves.append(dict(region=reg, obs="S-P3", ref="-", kind=kind, scheme=scheme, unit=unit, lab=lab,
                                               value=float((x[ok] * w[ok]).sum() / w[ok].sum()), n=int(ok.sum()),
                                               sig=float(1 / np.sqrt(w[ok].sum()))))
    np.savez_compressed(f"{L.ROOT3}/tables/phase_perm_dists{tag}.npz", **PERMSAVE)
    R = pd.DataFrame(rows); Cv = pd.DataFrame(curves)
    R.to_csv(f"{L.ROOT3}/tables/phase_contrasts{tag}.csv", index=False); Cv.to_csv(f"{L.ROOT3}/tables/phase_curves{tag}.csv", index=False)
    calls(R, tag)


def calls(R, tag):
    """detection rule (PLAN §2) on region-wide items; candidate = all but the cross-region clause"""
    out = []
    H = R[(R.X == "all") & (R.kind == "data")]
    for _, r in H.iterrows():
        ctrl = R[(R.region == r.region) & (R.obs == r.obs) & (R.ref == r.ref) & (R.kind == "ctrl") & (R.X == "all") & (R.unit == r.unit)]
        for c in CONTR:
            cp = float(ctrl[f"{c}_labp"].iloc[0]) if len(ctrl) else np.nan
            ok = (np.isfinite(r[c]) and abs(r[c]) >= 3 * r[f"{c}_sig"] and r[f"{c}_labp"] < 0.01 and
                  r.get(f"{c}_datep", 1.0) < 0.01 and not (np.isfinite(cp) and cp < 0.01))
            other = H[(H.region != r.region) & (H.obs == r.obs) & (H.unit == r.unit)]
            cross = bool(len(other) and any((np.sign(o[c]) == np.sign(r[c])) and (abs(o[c]) >= 1.96 * o[f"{c}_sig"]) for _, o in other.iterrows()))
            out.append(dict(region=r.region, obs=r.obs, ref=r.ref, unit=r.unit, contrast=c, value=r[c], sig=r[f"{c}_sig"],
                            sig_re=r[f"{c}_sig_re"], chi2dof=r[f"{c}_chi2dof"], labp=r[f"{c}_labp"], datep=r.get(f"{c}_datep"),
                            ctrl_labp=cp, n_series=r.n_series, bound2s=abs(r[c]) + 2 * r[f"{c}_sig"],
                            call="DETECTION" if (ok and cross) else ("candidate" if ok else "bound")))
    K = pd.DataFrame(out); K.to_csv(f"{L.ROOT3}/tables/phase_calls{tag}.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_rows", 300):
        print(K.round(4).to_string())


if __name__ == "__main__":
    main()
