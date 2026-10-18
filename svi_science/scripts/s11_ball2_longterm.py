#!/usr/bin/env python
"""s11_ball2_longterm -- S11 ballistic round 2, §3 pooled long-term constraints (round2/PLAN.md + Amendment 2).
Year (YR, headline) and two-year (Y2) bins. Items: S-P per station; far-station dS and dP3 with borehole-mean (headline),
B011, B009+B010, PGC and site-mean references (clock-corrected); co-located dS. Per item: family-common component per
year (one-way decomposition + within-series year-label permutation), decade contrast (>= 2015 minus <= 2014), split
contrast, linear trend, all with the same permutation control; scrambled control; conversion to dVs/Vs, dVp/Vp.

    python s11_ball2_longterm.py [--tag _test]
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from s11_ball_analysis import decomp
import s11_ball2_clock as CK

OUT = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/ballistic/round2"
TB, FIG = f"{OUT}/tables", f"{OUT}/figures"
STATIONS = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011", "PB.B926", "CN.VGZ", "CN.SNB"]
FAR = ["PB.B926", "CN.SNB", "CN.VGZ"]
REFS = {"BH": ["PB.B009", "PB.B010", "PB.B011"], "B011": ["PB.B011"], "B0910": ["PB.B009", "PB.B010"],
        "PGC": ["CN.PGC"], "SITE": ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]}
NPERM, SEED, SPLIT = 10000, 20261012, 2015
VPVS = 1.75


def year_of(lab):
    return int(str(lab)[:4])


def series_frame(df, val, sig, keys, cc):
    d = df[np.isfinite(df[val]) & np.isfinite(df[sig])].copy()
    d["K"] = d.groupby(keys)[val].transform("size"); d = d[d.K >= 3].copy()
    d["r"] = d[val] - d.groupby(keys)[val].transform("mean")
    d["sid"] = d.groupby(keys).ngroup(); d["sig"] = np.maximum(d[sig], 0.5)
    d["t"] = d.lab.map(year_of).astype(float); d["episode"] = d.lab; d["cc"] = cc(d) if callable(cc) else d[cc]
    return d.reset_index(drop=True)


def contrast_stats(d, rng, nperm=NPERM):
    """decade contrast, split contrast and trend with within-series permutation p"""
    out = {}
    if len(d) == 0:
        return out
    sid = d.sid.to_numpy(); y = d.r.to_numpy(); s2 = d.sig.to_numpy() ** 2; t = d.t.to_numpy()
    ns = sid.max() + 1
    late = t >= SPLIT
    med = d.groupby("sid").t.transform("median").to_numpy(); late2 = t > med
    defs = {}
    for nm, L in (("decade", late), ("split", late2)):
        nl = np.bincount(sid, L, ns); ne = np.bincount(sid, ~L, ns)
        ok = (nl >= 2) & (ne >= 2)
        if nm == "split":
            ok = (nl >= 1) & (ne >= 1)
        v = np.bincount(sid, s2 * L, ns) / np.maximum(nl, 1) ** 2 + np.bincount(sid, s2 * ~L, ns) / np.maximum(ne, 1) ** 2
        w = np.where(ok, 1 / np.where(v > 0, v, np.inf), 0.0)
        defs[nm] = (L, nl, ne, ok, w)
    w_i = 1 / s2
    tb = np.bincount(sid, w_i * t, ns) / np.bincount(sid, w_i, ns); tc = t - tb[sid]
    Stt = (w_i * tc * tc).sum()

    def stats(yv):
        res = {}
        for nm, (L, nl, ne, ok, w) in defs.items():
            C_s = np.bincount(sid, yv * L, ns) / np.maximum(nl, 1) - np.bincount(sid, yv * ~L, ns) / np.maximum(ne, 1)
            res[nm] = (w * C_s).sum() / w.sum() if w.sum() > 0 else np.nan
            res[nm + "_Cs"] = C_s
        res["trend"] = (w_i * tc * yv).sum() / Stt if Stt > 0 else np.nan
        return res
    obs = stats(y)
    order = np.argsort(sid, kind="stable"); sid_s = sid[order]
    perm = {k: np.empty(nperm) for k in ("decade", "split", "trend")}
    for i in range(nperm):
        key = sid_s + rng.random(len(sid_s)); pr = order[np.argsort(key, kind="stable")]
        yp = np.empty_like(y); yp[order] = y[pr]
        st = stats(yp)
        for k in perm:
            perm[k][i] = st[k]
    for nm in ("decade", "split"):
        L, nl, ne, ok, w = defs[nm]
        C = obs[nm]; out[f"{nm}_n_series"] = int(ok.sum())
        if ok.sum() == 0:
            continue
        sig = 1 / np.sqrt(w.sum()); chi2 = float((w * (obs[nm + "_Cs"] - C) ** 2)[ok].sum()); dof = max(ok.sum() - 1, 1)
        infl = max(1.0, np.sqrt(chi2 / dof))
        out.update({f"{nm}_ms": C, f"{nm}_sig_ms": sig * infl, f"{nm}_chi2dof": chi2 / dof,
                    f"{nm}_perm_p": float((1 + (np.abs(perm[nm]) >= abs(C)).sum()) / (nperm + 1)),
                    f"{nm}_perm_sd_ms": float(np.nanstd(perm[nm]))})
    resid = y - obs["trend"] * tc
    chi2t = float((w_i * resid ** 2).sum()); doft = max(len(y) - ns - 1, 1)
    out.update(trend_ms_per_yr=obs["trend"], trend_sig=1 / np.sqrt(Stt) * max(1, np.sqrt(chi2t / doft)),
               trend_perm_p=float((1 + (np.abs(perm["trend"]) >= abs(obs["trend"])).sum()) / (nperm + 1)))
    return out


def far_ref_frame(D, refset, valcol, far_list=FAR):
    F = D[D.B.isin(far_list) & D.A.isin(refset)]
    g = F.groupby(["fam", "B", "epB", "bin", "lab"])
    X = g.agg(v=(valcol, "mean"), sig=("sig", lambda s: np.sqrt((s ** 2).mean() / len(s))), date=("date", "min"),
              ccX=("ccX", "median"), nA=("A", "size")).reset_index()
    return X.rename(columns={"B": "far", "epB": "epoch"})


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default=""); a = ap.parse_args(); tag = a.tag
    rng = np.random.default_rng(SEED)
    st = pd.read_parquet(f"{TB}/lags_station{tag}.parquet")
    pr = pd.read_parquet(f"{TB}/lags_pair{tag}.parquet")
    geo = pd.read_csv(f"{TB}/series_geom{tag}.csv")
    SPm = geo.groupby(["fam", "netsta"]).SmP_pick.median()
    ver = json.load(open(f"{TB}/P_verdict{tag}.json"))
    if ver["P3"]["PASS"]:
        PE = "P3"
    else:
        H = pd.read_csv(f"{TB}/P_splithalf{tag}.csv").set_index("netsta")
        PE = max(["P3", "P3j", "Pz", "Pzr"], key=lambda k: H[f"half_cc{k}"].median())
    # ---- pair frames (data: clock-corrected via the S-based correction; P gets the same correction)
    Dc = pd.read_parquet(f"{TB}/pair_corrected{tag}.parquet")
    keys = ["fam", "pair", "epA", "epB", "bin", "lab"]
    DP = CK.prep(pr, "dP3", "sig_dP3", "railP3").merge(Dc[keys + ["kcorr", "has_corr"]], on=keys, how="left")
    DP["kcorr"] = DP.kcorr.fillna(0.0); DP["rc0"] = DP.r + DP.kcorr
    DP["rc"] = DP.rc0 - DP.groupby(CK.PKEY).rc0.transform("mean")
    for D_ in (Dc, DP):
        D_["ccX"] = np.minimum(D_.ccS_A, D_.ccS_B) if "ccS_A" in D_ else np.nan
    DP["ccX"] = np.minimum(DP.ccP3_A, DP.ccP3_B)
    # control pairs (uncorrected)
    pc = pr[pr.kind == "ctrl"].copy()
    rows, curves = [], []

    def run_item(name, d, conv, kind, bk, extra=None):
        o, ept = decomp(d.assign(date=d.date), "sig", name, rng) if len(d) else ({}, None)
        r = dict(item=name, bin=bk, kind=kind, n_values=len(d), n_series=d.sid.nunique() if len(d) else 0,
                 sd_year=np.sqrt(max(o.get("var_ep", np.nan), 0)) if np.isfinite(o.get("var_ep", np.nan)) else np.nan,
                 year_perm_p=o.get("perm_p", np.nan), scatter_ms=o.get("scatter_ms", np.nan))
        r.update(contrast_stats(d, rng, NPERM if kind == "data" else 2000) if len(d) else {})
        if extra:
            r.update(extra)
        for k_, T in conv.items():                               # conversions of the decade contrast
            if T and np.isfinite(T) and "decade_ms" in r:
                r[f"dec_{k_}_pct"] = -r["decade_ms"] / 1e3 / T * 100 if not k_.startswith("+") else r["decade_ms"] / 1e3 / T * 100
                r[f"dec_{k_}_sig_pct"] = r["decade_sig_ms"] / 1e3 / abs(T) * 100
                r[f"dec_{k_}_bound2s_pct"] = (abs(r["decade_ms"]) + 2 * r["decade_sig_ms"]) / 1e3 / abs(T) * 100
                sgc = max(r["decade_sig_ms"], r.get("decade_perm_sd_ms", 0.0))       # conservative: larger of formal / permutation SD
                r[f"dec_{k_}_bound2s_cons_pct"] = (abs(r["decade_ms"]) + 2 * sgc) / 1e3 / abs(T) * 100
        rows.append(r)
        if kind == "data" and len(d):
            cy = d.groupby("lab").apply(lambda x: pd.Series(dict(common=np.average(x.r, weights=1 / x.sig ** 2),
                                                                 n=len(x)))).reset_index()
            cy["item"], cy["bin"] = name, bk; curves.append(cy)
        return r

    for bk in ("YR", "Y2"):
        # ---- S-P per station
        for kind in ("data", "ctrl"):
            s = st[(st.kind == kind) & (st.bin == bk)].copy()
            s["SmP"] = np.where(s.railS | s[f"rail{PE}"], np.nan, s.lagS - s[f"lag{PE}"])
            s["sigSmP"] = s[f"sigSm{PE}"] if kind == "data" else np.hypot(s.sigS, s[f"sig{PE}"])
            for ns in STATIONS:
                d = series_frame(s[s.netsta == ns], "SmP", "sigSmP", ["fam", "netsta", "epoch", "bin"], f"cc{PE}")
                SP = float(np.nanmean([SPm.get((f, ns), np.nan) for f in d.fam.unique()])) if len(d) else np.nan
                run_item(f"S-P {ns}", d, {"constVpVs": SP, "VsOnly": SP * VPVS / (VPVS - 1), "+VpOnly": SP / (VPVS - 1)},
                         kind, bk, dict(P_est=PE, SP_s=SP))
        # ---- far dS / dP3 with references
        for obs, D_, cclab in (("dS", Dc, "ccS"), ("dP3", DP, "ccP3")):
            Db = D_[D_.bin == bk]
            for ref, rs in REFS.items():
                for valcol, corr in (("rc", "corrected"), ("r", "uncorrected")):
                    if corr == "uncorrected" and ref not in ("BH", "B011"):
                        continue
                    X = far_ref_frame(Db, rs, valcol)
                    for far in FAR:
                        d = series_frame(X[X.far == far].assign(lab=lambda x: x.lab), "v", "sig",
                                         ["fam", "far", "epoch", "bin"], "ccX")
                        TSX = float(np.nanmean([SPm.get((f, far), np.nan) for f in d.fam.unique()]) * VPVS / (VPVS - 1)) if len(d) else np.nan
                        TSA = float(np.nanmean([SPm.get((f, "PB.B011"), np.nan) for f in d.fam.unique()]) * VPVS / (VPVS - 1)) if len(d) else np.nan
                        if obs == "dS":
                            conv = {"dVs_Xleg": TSX, "dVs_uniform": TSX - TSA}
                        else:
                            conv = {"dVp_Xleg": TSX / VPVS, "dVp_uniform": (TSX - TSA) / VPVS}          # T_P = T_S / (Vp/Vs)
                        run_item(f"{obs} {ref}-{far} {corr}", d, conv, "data", bk, dict(TS_X=TSX, TS_ref=TSA))
                    # scrambled control (uncorrected, same reference averaging)
                    if corr == "corrected" and obs == "dS":
                        c = pc[(pc.bin == bk)].copy(); c.loc[c.railS, "dS"] = np.nan
                        c["sig"] = np.maximum(c.sig_dS, 0.5); c["ccX"] = np.minimum(c.ccS_A, c.ccS_B)
                        cpd = []
                        for _, g in c.groupby(["fam", "pair", "epA", "epB"]):
                            g = g[np.isfinite(g.dS)]
                            if len(g) >= 3:
                                cpd.append(g.assign(r=g.dS - g.dS.mean()))
                        if cpd:
                            C_ = pd.concat(cpd)
                            Xc = far_ref_frame(C_, rs, "r")
                            for far in FAR:
                                d = series_frame(Xc[Xc.far == far], "v", "sig", ["fam", "far", "epoch", "bin"], "ccX")
                                run_item(f"{obs} {ref}-{far} corrected", d, {}, "ctrl", bk)
        # ---- co-located (corrected) at year scale
        Db = Dc[Dc.bin == bk]
        for ref, pl in CK.FLOOR_SETS.items():
            d = Db[Db.pair.isin(pl)].copy()
            d = series_frame(d.assign(v=d.rc), "v", "sig", CK.PKEY, "ccX")
            run_item(f"coloc dS {ref}-set corrected", d, {}, "data", bk)
    R = pd.DataFrame(rows); R.to_csv(f"{TB}/longterm{tag}.csv", index=False)
    Cv = pd.concat(curves, ignore_index=True) if curves else pd.DataFrame(); Cv.to_csv(f"{TB}/longterm_curves{tag}.csv", index=False)
    # ---- detection / bound calls (headline YR, decade)
    calls = []
    for _, r in R[(R.bin == "YR") & (R.kind == "data")].iterrows():
        it = r["item"]
        if "uncorrected" in it or it.startswith("coloc"):
            continue
        if it.startswith(("dS", "dP3")) and not it.split()[1].startswith("BH-"):
            continue
        cn = R[(R["item"] == it) & (R.bin == "YR") & (R.kind == "ctrl")]
        ctrl_p = float(cn.decade_perm_p.iloc[0]) if len(cn) and "decade_perm_p" in cn else np.nan
        det = (abs(r.get("decade_ms", 0)) >= 3 * r.get("decade_sig_ms", np.inf)) and (r.get("decade_perm_p", 1) < 0.01) \
            and not (np.isfinite(ctrl_p) and ctrl_p < 0.01)
        if det and it.startswith(("dS", "dP3")):
            obs, rf = it.split()[0], it.split()[1]; far = rf.split("-", 1)[1]
            ok2 = True
            for alt in ("B011", "B0910"):
                q = R[(R["item"] == f"{obs} {alt}-{far} corrected") & (R.bin == "YR") & (R.kind == "data")]
                ok2 &= bool(len(q) and np.sign(q.decade_ms.iloc[0]) == np.sign(r["decade_ms"]) and q.decade_perm_p.iloc[0] < 0.05)
            det = det and ok2
        b = {k: r[k] for k in r.index if str(k).startswith("dec_") and str(k).endswith("bound2s_pct")}
        calls.append(dict(item=it, decade_ms=r.get("decade_ms"), decade_sig_ms=r.get("decade_sig_ms"),
                          perm_p=r.get("decade_perm_p"), ctrl_perm_p=ctrl_p, detection=bool(det),
                          trend_ms_per_yr=r.get("trend_ms_per_yr"), trend_sig=r.get("trend_sig"),
                          trend_perm_p=r.get("trend_perm_p"), year_perm_p=r.get("year_perm_p"), sd_year=r.get("sd_year"), **b))
    K = pd.DataFrame(calls); K.to_csv(f"{TB}/longterm_calls{tag}.csv", index=False)
    with pd.option_context("display.width", 300, "display.max_rows", 400, "display.max_columns", 40):
        print("P estimator for S-P:", PE)
        print(K.round(4).to_string())
        cols = ["item", "bin", "kind", "n_values", "n_series", "sd_year", "year_perm_p", "decade_ms", "decade_sig_ms",
                "decade_perm_p", "decade_n_series", "split_ms", "split_sig_ms", "split_perm_p", "trend_ms_per_yr", "trend_sig",
                "trend_perm_p"]
        print(R[[c for c in cols if c in R]].round(3).to_string())
    figure(Dc, DP, st, PE, R, tag)


def figure(Dc, DP, st, PE, R, tag):
    from matplotlib.colors import Normalize
    norm = Normalize(0.3, 1.0); cmap = plt.get_cmap("viridis")
    panels = [("S-P", ns) for ns in STATIONS] + [("dS", f) for f in FAR] + [("dP3", f) for f in FAR]
    fig, axs = plt.subplots(len(panels), 1, figsize=(12, 2.0 * len(panels)), sharex=True)
    s = st[(st.kind == "data") & (st.bin == "YR")].copy()
    s["SmP"] = np.where(s.railS | s[f"rail{PE}"], np.nan, s.lagS - s[f"lag{PE}"])
    for ax, (obs, ns) in zip(axs, panels):
        if obs == "S-P":
            d = series_frame(s[s.netsta == ns], "SmP", f"sigSm{PE}", ["fam", "netsta", "epoch", "bin"], f"cc{PE}")
            name = f"S-P {ns}"
        else:
            D_ = Dc if obs == "dS" else DP
            X = far_ref_frame(D_[D_.bin == "YR"], REFS["BH"], "rc")
            d = series_frame(X[X.far == ns], "v", "sig", ["fam", "far", "epoch", "bin"], "ccX")
            name = f"{obs} BH-{ns} corrected"
        if not len(d):
            continue
        d["yr"] = d.t + 0.5
        for _, g in d.groupby("sid"):
            g = g.sort_values("yr"); ax.plot(g.yr, g.r, color="0.82", lw=0.5, zorder=1)
        sc = ax.scatter(d.yr, d.r, c=d.cc, cmap=cmap, norm=norm, s=8, lw=0, zorder=2)
        cy = d.groupby("yr").apply(lambda x: np.average(x.r, weights=1 / x.sig ** 2)).sort_index()
        ax.plot(cy.index, cy.values, color="k", lw=1.6, zorder=3)
        lim = np.nanpercentile(np.abs(d.r), 97) * 1.2
        ax.set_ylim(-lim, lim); ax.axhline(0, color="0.5", lw=0.5); ax.axvline(2015, color="r", lw=0.6, ls=":")
        row = R[(R["item"] == name) & (R.bin == "YR") & (R.kind == "data")]
        txt = ""
        if len(row) and "decade_ms" in row and np.isfinite(row.decade_ms.iloc[0]):
            r0 = row.iloc[0]; txt = f"  decade contrast {r0.decade_ms:+.1f} ± {r0.decade_sig_ms:.1f} ms (perm p {r0.decade_perm_p:.3f})"
        ax.set_title(name.replace("CN.", "").replace("PB.", "") + txt, fontsize=8, loc="left")
        ax.set_ylabel("ms", fontsize=8)
    fig.colorbar(sc, ax=axs, shrink=0.25, label="cc")
    fig.suptitle("Year stacks: grey = one family series (demeaned), black = family-common (1/σ² mean); red dotted = 2015 split",
                 fontsize=10)
    fig.savefig(f"{FIG}/fig_longterm_years{tag}.png", dpi=105, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
