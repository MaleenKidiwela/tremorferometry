#!/usr/bin/env python
"""s11_ball3_dh -- S11 ballistic round 3, Amendment 5: decade-half (DH) stacks and 4-year (Q4) time course, both regions.
Per series contrast C = [lag(L) - lag(E)] / 2 (Amendment 5a), var = sig_E^2 + sig_L^2. Group clock terms on within-group
pair contrasts (median gauge, leave-family-out over >= 3 other families); reference means; 1/var pooling with chi2
inflation; 10,000 random sign flips; scrambled control. Rule 1 re-judged on DH; region-wide S - P; drift separation;
joint (a, b); site-joint S - P; Q4 curves; figure.

    python s11_ball3_dh.py [--tag _dh]
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_ball3_lib as L

NFLIP = 10000
SEED = 20261014
rng = np.random.default_rng(SEED)


def flip_p(C, v, scale=None, n=NFLIP):
    sc = np.ones(len(C)) if scale is None else scale
    x = C * sc; w = 1 / (v * sc ** 2); ok = np.isfinite(x) & np.isfinite(w)
    x, w = x[ok], w[ok]
    if len(x) < 2:
        return np.nan
    m0 = (w * x).sum() / w.sum()
    S = rng.choice([-1.0, 1.0], size=(n, len(x)))
    mp = (S * (w * x)).sum(1) / w.sum()
    return float((1 + (np.abs(mp) >= abs(m0)).sum()) / (n + 1))


def pooled(C, v, groups=None, scale=None):
    sc = np.ones(len(C)) if scale is None else np.asarray(scale, float)
    ok = np.isfinite(C) & np.isfinite(v) & np.isfinite(sc) & (v > 0)
    if ok.sum() == 0:
        return dict(mean=np.nan, sig=np.nan, n=0, flip_p=np.nan)
    pr = L.pool_re((C * sc)[ok], (np.sqrt(v) * np.abs(sc))[ok], None if groups is None else np.asarray(groups)[ok])
    pr["flip_p"] = flip_p(C[ok], v[ok], sc[ok])
    return pr


def contrasts(df, keys, val, sig, rail):
    d = df[(df.bin == "DH") & ~df[rail]].copy()
    d = d[np.isfinite(d[val]) & np.isfinite(d[sig])]
    piv = d.pivot_table(index=keys, columns="lab", values=[val, sig], aggfunc="first")
    if ("%s" % val, "E") not in piv or (val, "L") not in piv:
        return pd.DataFrame(columns=keys + ["C", "v"])
    out = piv.index.to_frame(index=False)
    out["C"] = (piv[(val, "L")].to_numpy() - piv[(val, "E")].to_numpy()) / 2
    out["v"] = np.maximum(piv[(sig, "L")].to_numpy(), 0.5) ** 2 + np.maximum(piv[(sig, "E")].to_numpy(), 0.5) ** 2
    return out[np.isfinite(out.C) & np.isfinite(out.v)].reset_index(drop=True)


def clock_dh(Cp, region):
    """group clock terms on within-group pair contrasts -> dict (fam, node) -> leave-family-out k"""
    ng = L.node_group(region)
    Cp = Cp.assign(nodeA=Cp.A + "|" + Cp.epA, nodeB=Cp.B + "|" + Cp.epB, gA=Cp.A.map(ng), gB=Cp.B.map(ng))
    Q = Cp[Cp.gA.notna() & (Cp.gA == Cp.gB)]
    rows = []
    for (f, g_), g in Q.groupby(["fam", "gA"]):
        nodes = sorted(set(g.nodeA) | set(g.nodeB))
        if len(nodes) < 2:
            continue
        ix = {n: i for i, n in enumerate(nodes)}
        Am = np.zeros((len(g), len(nodes))); w = 1 / np.sqrt(g.v.to_numpy())
        Am[np.arange(len(g)), [ix[n] for n in g.nodeB]] = 1; Am[np.arange(len(g)), [ix[n] for n in g.nodeA]] = -1
        k = np.linalg.lstsq(Am * w[:, None], g.C.to_numpy() * w, rcond=None)[0]; k = k - np.median(k)
        prec = np.abs(Am).T @ (w ** 2)
        for n in nodes:
            rows.append(dict(fam=f, node=n, k=k[ix[n]], var=1 / max(prec[ix[n]], 1e-9)))
    T = pd.DataFrame(rows)
    if not len(T):
        return {}, T
    T["w"] = 1 / T["var"]
    g = T.groupby("node")
    T["W"] = g.w.transform("sum"); T["WK"] = (T.w * T.k).groupby(T.node).transform("sum"); T["nf"] = g.k.transform("size")
    T["K_lfo"] = (T.WK - T.w * T.k) / (T.W - T.w); T.loc[T.nf - 1 < 3, "K_lfo"] = np.nan
    return {(f, n): k for f, n, k in zip(T.fam, T.node, T.K_lfo)}, T


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default="_dh"); a = ap.parse_args(); tag = a.tag
    pairs_out, st_out, reg_out, joint_out, q4_out, sj_out, clk_out = [], [], [], [], [], [], []
    verdict = {}
    for reg in ("SVI", "WA_NW"):
        TB = L.dirs(reg)[1]
        if not os.path.exists(f"{TB}/lags_pair{tag}.parquet"):
            continue
        st = L.load(reg, "lags_station", tag); pr = L.load(reg, "lags_pair", tag); SP = L.sp_table(reg, tag)
        ng = L.node_group(reg)
        res = {}
        for kind in ("data", "ctrl"):
            p = pr[pr.kind == kind]
            for obs, val, sig, rail in (("dS", "dS", "sig_dS", "railS"), ("dP3", "dP3", "sig_dP3", "railP3")):
                Cp = contrasts(p, ["fam", "pair", "A", "B", "epA", "epB"], val, sig, rail)
                kmap, T = clock_dh(Cp, reg) if (kind == "data" and obs == "dS") else (res.get(("kmap",)), None)
                if kind == "data" and obs == "dS":
                    res[("kmap",)] = kmap
                    if T is not None and len(T):
                        clk_out.append(T.assign(region=reg))
                kmap = kmap or {}
                kA = np.array([kmap.get((f, a_ + "|" + e), np.nan) for f, a_, e in zip(Cp.fam, Cp.A, Cp.epA)])
                kB = np.array([kmap.get((f, b_ + "|" + e), np.nan) for f, b_, e in zip(Cp.fam, Cp.B, Cp.epB)])
                inA = Cp.A.map(ng).notna().to_numpy(); inB = Cp.B.map(ng).notna().to_numpy()
                corr = np.where(inA, np.nan_to_num(kA), 0) - np.where(inB, np.nan_to_num(kB), 0)
                Cp["Cc"] = Cp.C + (corr if kind == "data" else 0.0)
                res[(obs, kind)] = Cp
                # oriented and reference means
                a1 = Cp.assign(X=Cp.B, epX=Cp.epB, M=Cp.A, cv=Cp.Cc); a2 = Cp.assign(X=Cp.A, epX=Cp.epA, M=Cp.B, cv=-Cp.Cc)
                Lg = pd.concat([a1, a2], ignore_index=True)
                for rname, rs in L.REFSETS[reg].items():
                    F = Lg[Lg.M.isin(rs) & ~Lg.X.isin(rs)].groupby(["fam", "X", "epX"]).agg(
                        C=("cv", "mean"), v=("v", lambda s: s.mean() / len(s)), nM=("M", "size")).reset_index()
                    for X, g in F.groupby("X"):
                        tsx = np.array([L.TS(SP.get((f, X), np.nan)) for f in g.fam])
                        tsr = np.array([np.nanmean([L.TS(SP.get((f, m), np.nan)) for m in rs]) for f in g.fam])
                        pr_ = pooled(g.C.to_numpy(), g.v.to_numpy())
                        leg = pooled(g.C.to_numpy(), g.v.to_numpy(), scale=-100 / (tsx * 1e3))
                        dT = np.nanmean(tsx - tsr)
                        uni = pooled(g.C.to_numpy(), g.v.to_numpy(), scale=-100 / ((tsx - tsr) * 1e3)) if abs(dT) >= 2 else {}
                        pairs_out.append(dict(region=reg, obs=obs, kind=kind, ref=rname, X=X, n_series=pr_["n"], C_ms=pr_["mean"],
                                              sig_ms=pr_["sig"], chi2dof=pr_.get("chi2dof"), flip_p=pr_["flip_p"], TS_X=np.nanmean(tsx),
                                              TS_ref=np.nanmean(tsr), leg_pct=leg["mean"], leg_sig_pct=leg["sig"],
                                              leg_bound2s_pct=abs(leg["mean"]) + 2 * leg["sig"] if leg["n"] else np.nan,
                                              uni_pct=uni.get("mean", np.nan), uni_bound2s_pct=(abs(uni["mean"]) + 2 * uni["sig"]) if uni.get("n") else np.nan))
            # co-located (uncorrected and corrected)
            Cp = res[("dS", kind)]
            W = Cp[Cp.A.map(ng).notna() & (Cp.A.map(ng) == Cp.B.map(ng))]
            for pair, g in W.groupby("pair"):
                for col in ("C", "Cc"):
                    pr_ = pooled(g[col].to_numpy(), g.v.to_numpy())
                    pairs_out.append(dict(region=reg, obs="coloc" + ("" if col == "C" else "_corr"), kind=kind, ref=g.A.map(ng).iloc[0],
                                          X=pair, n_series=pr_["n"], C_ms=pr_["mean"], sig_ms=pr_["sig"], chi2dof=pr_.get("chi2dof"),
                                          flip_p=pr_["flip_p"]))
            # station S-P, absS, absP3
            s = st[st.kind == kind].copy()
            s["SmP"] = s.lagS - s.lagP3; s["rSmP"] = s.railS | s.railP3
            s["SmP3c"] = s.lagS3c - s.lagP3; s["rSmP3c"] = s.railS3c | s.railP3
            items = {"S-P3": ("SmP", "sigSmP3", "rSmP"), "absS": ("lagS", "sigS", "railS"), "absP3": ("lagP3", "sigP3", "railP3")}
            if kind == "data":
                items["S3c-P3"] = ("SmP3c", "sigS3cmP3", "rSmP3c")
            for obs, (val, sig, rail) in items.items():
                Cs = contrasts(s, ["fam", "netsta", "epoch"], val, sig, rail)
                if not len(Cs):
                    continue
                Cs["sp"] = [SP.get((f, n), np.nan) for f, n in zip(Cs.fam, Cs.netsta)]
                units = {"ms": None}
                if obs.startswith(("S-", "S3c")):
                    units.update({"constVpVs": -100 / (Cs.sp.to_numpy() * 1e3), "VsOnly": -100 / (L.TS(Cs.sp.to_numpy()) * 1e3),
                                  "VpOnly": 100 / (L.TP(Cs.sp.to_numpy()) * 1e3)})
                for ns, g in Cs.groupby("netsta"):
                    pr_ = pooled(g.C.to_numpy(), g.v.to_numpy())
                    pc_ = pooled(g.C.to_numpy(), g.v.to_numpy(), scale=-100 / (g.sp.to_numpy() * 1e3)) if obs.startswith(("S-", "S3c")) else {}
                    st_out.append(dict(region=reg, obs=obs, kind=kind, X=ns, n_series=pr_["n"], C_ms=pr_["mean"], sig_ms=pr_["sig"],
                                       chi2dof=pr_.get("chi2dof"), flip_p=pr_["flip_p"], constVpVs_pct=pc_.get("mean"), constVpVs_sig=pc_.get("sig")))
                for unit, sc in units.items():
                    for excl in ("none", "B013") if (reg == "WA_NW" and obs.startswith("S-")) else ("none",):
                        m = np.ones(len(Cs), bool) if excl == "none" else (Cs.netsta != "PB.B013").to_numpy()
                        pr_ = pooled(Cs.C.to_numpy()[m], Cs.v.to_numpy()[m], Cs.netsta.to_numpy()[m], None if sc is None else sc[m])
                        reg_out.append(dict(region=reg, obs=obs, kind=kind, unit=unit, excluded=excl, n_series=pr_["n"],
                                            n_st=pr_.get("n_groups"), mean=pr_["mean"], sig=pr_["sig"], sig_formal=pr_.get("sig_formal"),
                                            sig_re=pr_.get("sig_re"), chi2dof=pr_.get("chi2dof"), flip_p=pr_["flip_p"],
                                            bound2s=abs(pr_["mean"]) + 2 * pr_["sig"] if pr_["n"] else np.nan))
                if obs == "S-P3":
                    res[("SmP", kind)] = Cs
        # joint (a, b) on DH contrasts, headline reference
        for jname, (rname, rs) in (("headline", (L.HEAD_REFS[reg][0], L.REFSETS[reg][L.HEAD_REFS[reg][0]])),
                                   ("column", (L.HEAD_REFS[reg][1], L.REFSETS[reg][L.HEAD_REFS[reg][1]]))):
            for excl in ("none", "B013"):
                if excl == "B013" and reg != "WA_NW":
                    continue
                Cs = res[("SmP", "data")]; Cs = Cs[Cs.netsta != "PB.B013"] if excl == "B013" else Cs
                C_l, V_l, D_l = [Cs.C.to_numpy() / 1e3], [Cs.v.to_numpy() / 1e6], [np.column_stack([-L.TS(Cs.sp.to_numpy()), L.TP(Cs.sp.to_numpy())])]
                for obs in ("dS", "dP3"):
                    Cp = res[(obs, "data")]
                    if excl == "B013":
                        Cp = Cp[(Cp.A != "PB.B013") & (Cp.B != "PB.B013")]
                    a1 = Cp.assign(X=Cp.B, M=Cp.A, cv=Cp.Cc); a2 = Cp.assign(X=Cp.A, M=Cp.B, cv=-Cp.Cc)
                    Lg = pd.concat([a1, a2]); F = Lg[Lg.M.isin(rs) & ~Lg.X.isin(rs)].groupby(["fam", "X"]).agg(
                        C=("cv", "mean"), v=("v", lambda s: s.mean() / len(s))).reset_index()
                    tsx = np.array([L.TS(SP.get((f, x), np.nan)) for f, x in zip(F.fam, F.X)])
                    tsr = np.array([np.nanmean([L.TS(SP.get((f, m), np.nan)) for m in rs]) for f in F.fam])
                    dT = tsx - tsr
                    D_l.append(np.column_stack([-dT, 0 * dT]) if obs == "dS" else np.column_stack([0 * dT, -dT / L.VPVS]))
                    C_l.append(F.C.to_numpy() / 1e3); V_l.append(F.v.to_numpy() / 1e6)
                C = np.concatenate(C_l); V = np.concatenate(V_l); Dm = np.concatenate(D_l)
                ok = np.isfinite(C) & np.isfinite(V) & np.isfinite(Dm).all(1); C, V, Dm = C[ok], V[ok], Dm[ok]; w = 1 / V
                N = (Dm * w[:, None]).T @ Dm; x = np.linalg.solve(N, (Dm * w[:, None]).T @ C); cov = np.linalg.inv(N)
                chi = float((w * (C - Dm @ x) ** 2).sum()) / max(len(C) - 2, 1); infl = max(1, np.sqrt(chi))
                S_ = rng.choice([-1.0, 1.0], size=(2000, len(C)))
                xp = np.linalg.solve(N, ((Dm * w[:, None]).T @ (S_ * C).T)).T
                sa, sb = np.sqrt(np.diag(cov)) * infl; sab = np.sqrt(cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]) * infl
                joint_out.append(dict(region=reg, fit=jname, ref=rname, excluded=excl, n=len(C), chi2dof=chi,
                                      a_dVs_pct=x[0] * 100, a_sig=sa * 100, a_flip_p=float((1 + (np.abs(xp[:, 0]) >= abs(x[0])).sum()) / 2001),
                                      b_dVp_pct=x[1] * 100, b_sig=sb * 100, b_flip_p=float((1 + (np.abs(xp[:, 1]) >= abs(x[1])).sum()) / 2001),
                                      bma_pct=(x[1] - x[0]) * 100, bma_sig=sab * 100,
                                      a_bound2s=(abs(x[0]) + 2 * sa) * 100, b_bound2s=(abs(x[1]) + 2 * sb) * 100,
                                      bma_bound2s=(abs(x[1] - x[0]) + 2 * sab) * 100, corr_ab=cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])))
        # site-joint S-P (DH)
        p = f"{TB}/site_joint{tag}.parquet"
        if os.path.exists(p):
            J = pd.read_parquet(p); J = J[(J.kind == "data") & (J.bin == "DH")].copy()
            J["SmP"] = J.lagS - J.lagP3; J.loc[J.railS | J.railP3, "SmP"] = np.nan
            for grp, g in J.groupby("group"):
                piv = g.pivot_table(index="fam", columns="lab", values=["SmP", "sigSmP3"])
                if ("SmP", "L") in piv and ("SmP", "E") in piv:
                    C = ((piv[("SmP", "L")] - piv[("SmP", "E")]) / 2).to_numpy(); v = (piv[("sigSmP3", "L")] ** 2 + piv[("sigSmP3", "E")] ** 2).to_numpy()
                    pr_ = pooled(C, v); sj_out.append(dict(region=reg, group=grp, estimator="joint", n_series=pr_["n"], C_ms=pr_["mean"], sig_ms=pr_["sig"], flip_p=pr_["flip_p"]))
                for c in g.columns:
                    if c.startswith("m_") and c.endswith("_lagSmP3"):
                        nm = c[2:-8]; pv = g.pivot_table(index="fam", columns="lab", values=[c, c.replace("lagSmP3", "sigSmP3")])
                        if (c, "L") in pv and (c, "E") in pv:
                            C = ((pv[(c, "L")] - pv[(c, "E")]) / 2).to_numpy()
                            v = (pv[(c.replace("lagSmP3", "sigSmP3"), "L")] ** 2 + pv[(c.replace("lagSmP3", "sigSmP3"), "E")] ** 2).to_numpy()
                            pr_ = pooled(C, v); sj_out.append(dict(region=reg, group=grp, estimator=nm, n_series=pr_["n"], C_ms=pr_["mean"], sig_ms=pr_["sig"], flip_p=pr_["flip_p"]))
        # Q4 curves (data): S-P per station (selected), coloc raw, refs
        s = st[(st.kind == "data") & (st.bin == "Q4")].copy(); s["SmP"] = np.where(s.railS | s.railP3, np.nan, s.lagS - s.lagP3)
        s["sp"] = [SP.get((f, n), np.nan) for f, n in zip(s.fam, s.netsta)]
        s["x"] = -100 * s.SmP / (s.sp * 1e3); s["sx"] = 100 * np.maximum(s.sigSmP3, 0.5) / (s.sp * 1e3)
        for nm, sel in [("S-P3 region (constVpVs %)", np.ones(len(s), bool))] + [(f"S-P3 {x} (ms)", (s.netsta == x).to_numpy())
                                                                                for x in ("PB.B013", "PB.B943", "PB.B001", "PB.B005", "PB.B011", "PB.B010", "CN.PGC")]:
            col, scol = ("x", "sx") if "region" in nm else ("SmP", "sigSmP3")
            d = s[sel & np.isfinite(s[col])].copy()
            d["K"] = d.groupby(["fam", "netsta", "epoch"])[col].transform("size"); d = d[d.K >= 3]
            if not len(d):
                continue
            d["r"] = d[col] - d.groupby(["fam", "netsta", "epoch"])[col].transform("mean"); d["w"] = 1 / np.maximum(d[scol], 1e-6) ** 2
            for lab, g in d.groupby("lab"):
                q4_out.append(dict(region=reg, item=nm, lab=lab, value=np.average(g.r, weights=g.w), n=len(g), cc=g.ccP3.median()))
        Pp = pr[(pr.kind == "data") & (pr.bin == "Q4") & ~pr.railS].copy()
        Lg = pd.concat([Pp.assign(X=Pp.B, epX=Pp.epB, M=Pp.A, cv=Pp.dS), Pp.assign(X=Pp.A, epX=Pp.epA, M=Pp.B, cv=-Pp.dS)])
        Lg["ccX"] = np.minimum(Pp.ccS_A, Pp.ccS_B).tolist() * 2
        for rname in L.HEAD_REFS[reg][:2]:
            rs = L.REFSETS[reg][rname]
            F = Lg[Lg.M.isin(rs) & ~Lg.X.isin(rs)].groupby(["fam", "X", "epX", "lab"]).agg(v=("cv", "mean"), s=("sig_dS", "median"), cc=("ccX", "median")).reset_index()
            for X in F.X.unique():
                d = F[F.X == X].copy(); d["K"] = d.groupby(["fam", "epX"]).v.transform("size"); d = d[d.K >= 3]
                if d.fam.nunique() < 5:
                    continue
                d["r"] = d.v - d.groupby(["fam", "epX"]).v.transform("mean"); d["w"] = 1 / np.maximum(d.s, 0.5) ** 2
                for lab, g in d.groupby("lab"):
                    q4_out.append(dict(region=reg, item=f"dS {rname}-{X} (ms, uncorrected)", lab=lab, value=np.average(g.r, weights=g.w), n=len(g), cc=g.cc.median()))
    P = pd.DataFrame(pairs_out); Sx = pd.DataFrame(st_out); Rg = pd.DataFrame(reg_out); Jt = pd.DataFrame(joint_out)
    Q = pd.DataFrame(q4_out); SJ = pd.DataFrame(sj_out)
    T3 = f"{L.ROOT3}/tables"
    for nm, df in (("dh_pairs", P), ("dh_stations", Sx), ("dh_region", Rg), ("dh_joint_ab", Jt), ("dh_q4_curves", Q), ("dh_sitejoint", SJ)):
        df.to_csv(f"{T3}/{nm}.csv", index=False)
    if clk_out:
        pd.concat(clk_out).to_csv(f"{T3}/dh_clock_terms.csv", index=False)
    # ---- rule 1 on DH (WA), detection clauses
    v = {}
    for reg in ("SVI", "WA_NW"):
        head = L.HEAD_REFS[reg][:2]
        D = P[(P.region == reg) & (P.obs == "dS") & (P.kind == "data")]
        Cc = P[(P.region == reg) & (P.obs == "dS") & (P.kind == "ctrl")].set_index(["ref", "X"])
        co = P[(P.region == reg) & (P.obs == "coloc") & (P.kind == "data")]
        legs, dets = [], []
        for X, g in D[D.ref.isin(head)].groupby("X"):
            g = g.set_index("ref")
            if not all(r in g.index for r in head):
                continue
            b = [g.loc[r, "leg_bound2s_pct"] for r in head]; n = [g.loc[r, "n_series"] for r in head]
            if all(np.isfinite(b)) and max(b) <= 0.03 and min(n) >= 20:
                legs.append(dict(X=X, bounds=b, n=n))
            for r in head:
                x = g.loc[r]; cp = Cc.flip_p.get((r, X), np.nan)
                single = abs(x.C_ms) >= 3 * x.sig_ms and x.flip_p < 0.01 and not (np.isfinite(cp) and cp < 0.01)
                oth = [g.loc[o] for o in head if o != r]
                if single and all(np.sign(o.C_ms) == np.sign(x.C_ms) and o.flip_p < 0.05 for o in oth):
                    dets.append(dict(ref=r, X=X, C_ms=x.C_ms, sig=x.sig_ms, p=x.flip_p))
        v[reg] = dict(coloc={r.X: dict(C=r.C_ms, sig=r.sig_ms, p=r.flip_p, within3s=bool(abs(r.C_ms) < 3 * r.sig_ms)) for _, r in co.iterrows()},
                      legs_le_003=legs, detections=dets)
        v[reg]["rule_a"] = bool(all(x["within3s"] for x in v[reg]["coloc"].values()))
        v[reg]["PASS"] = bool(v[reg]["rule_a"] and (len(legs) > 0 or len(dets) > 0))
    json.dump(v, open(f"{T3}/rule1_dh.json", "w"), indent=1, default=float)
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        print(P[P.kind == "data"].round(4).to_string()); print(Sx.round(4).to_string()); print(Rg.round(4).to_string())
        print(Jt.round(4).T.to_string()); print(SJ.round(3).to_string()); print(json.dumps(v, indent=1, default=float))
    figure(P, Q)


def figure(P, Q):
    cm = plt.get_cmap("viridis"); nm = Normalize(0.3, 1.0)
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(2, 1, 1)
    D = P[(P.region == "WA_NW") & (P.obs == "dS") & (P.kind == "data") & P.ref.isin(["G1BH", "G2"]) & (P.n_series >= 5)]
    order = D.groupby("X").TS_X.mean().sort_values().index.tolist(); xi = {x: i for i, x in enumerate(order)}
    for ref, mk, col in (("G1BH", "o", "C0"), ("G2", "s", "C3")):
        q = D[D.ref == ref]; ax.plot([xi[x] for x in q.X], q.leg_bound2s_pct, mk, color=col, label=f"reference {ref}")
    ax.axhline(0.03, color="r", ls="--", lw=0.8, label="rule (0.03 %)"); ax.axhline(0.014, color="k", ls=":", lw=0.8, label="SVI B926 leg (0.014 %)")
    ax.set_yscale("log"); ax.set_xticks(range(len(order))); ax.set_xticklabels([x.split(".")[1] for x in order], rotation=90, fontsize=7)
    ax.set_ylabel("2σ bound on δVs/Vs on the leg (%)"); ax.legend(fontsize=7)
    ax.set_title("WA-NW, decade-half stacks (≥2015 vs ≤2014): bound per station (sorted by S travel time)", fontsize=9, loc="left")
    items = [("WA_NW", "dS G1BH-PB.B013 (ms, uncorrected)"), ("WA_NW", "dS G1BH-PB.B943 (ms, uncorrected)"), ("WA_NW", "dS G2-PB.B001 (ms, uncorrected)"),
             ("WA_NW", "S-P3 PB.B013 (ms)"), ("WA_NW", "S-P3 PB.B943 (ms)"), ("WA_NW", "S-P3 PB.B001 (ms)"),
             ("SVI", "S-P3 PB.B010 (ms)"), ("SVI", "S-P3 PB.B011 (ms)"), ("SVI", "dS BH-PB.B926 (ms, uncorrected)")]
    for i, (reg, it) in enumerate(items):
        ax = fig.add_subplot(4, 5, 11 + i)
        q = Q[(Q.region == reg) & (Q["item"] == it)].sort_values("lab")
        if len(q):
            ax.plot(range(len(q)), q.value, color="0.4", lw=1); ax.scatter(range(len(q)), q.value, c=q.cc, cmap=cm, norm=nm, s=25, zorder=3)
            ax.set_xticks(range(len(q))); ax.set_xticklabels([l[1:] for l in q.lab], rotation=60, fontsize=6)
        ax.axhline(0, color="0.6", lw=0.5); ax.set_title(f"{reg} {it}".replace("PB.", "").replace("CN.", ""), fontsize=7)
    fig.tight_layout(); fig.savefig(f"{L.ROOT3}/figures/fig_dh_decade.png", dpi=105); plt.close(fig)


if __name__ == "__main__":
    main()
