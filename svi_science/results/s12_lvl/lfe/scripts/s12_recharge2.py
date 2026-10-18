#!/usr/bin/env python
"""s12_recharge2 -- PLAN §7: confirmation tests of the between-ETS recharge slope from b_days.parquet (s12_b_measure).
    S11_REGION=WA_NW|SVI python s12_recharge2.py
Outputs tables/{sub}/r2_*.csv, r2_verdict.json; figure figures/fig_recharge_{sub}.png.
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import pandas as pd
from scipy.stats import binomtest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s12_epoch import C, OUT, SUB, TB, log
from s12_novel import within, wls

rng = np.random.default_rng(77)
WL = ["1-2", "2-3", "3-4", "4-6", "2-5"]
PAT = {"LVL": [0.191, 0.223, 0.216, 0.201], "shallow 0-5 km": [0.166, 0.263, 0.304, 0.321], "25-45 km": [0.649, 0.543, 0.505, 0.445]}
NB = int(os.environ.get("S12_NB", 200))


def prep(F, q, extra=None):
    wn = q.split("_")[1] if q.startswith("dvv") else None
    d = F[np.isfinite(F[q])].copy()
    if wn:
        d = d[d[f"dvv_cc{wn}"] > 0.3]
        d["w"] = d[f"dvv_cc{wn}"] ** 2
    else:
        d["w"] = 1.0
    if extra is not None:
        d = d[extra(d)]
    return d


def annual_cols(d, by="site", semi=False):
    cols = []
    for s in sorted(d[by].unique()):
        m = (d[by] == s).astype(float)
        d[f"c_{s}"] = d.cs * m; d[f"s_{s}"] = d.sn * m; cols += [f"c_{s}", f"s_{s}"]
        if semi:
            d[f"c2_{s}"] = d.cs2 * m; d[f"s2_{s}"] = d.sn2 * m; cols += [f"c2_{s}", f"s2_{s}"]
    return cols


def fit(d, q, reg=("ty",), by="site", semi=False, extra_cols=(), boot=True, ret_all=False):
    d = d.copy()
    cols = list(reg) + list(extra_cols) + annual_cols(d, by, semi)
    w = d.w.to_numpy()
    G = within(d, [q] + cols, pd.Series(w, index=d.index))
    b = wls(G[q].to_numpy(), G[cols].to_numpy(), w)
    out = dict(b=float(b[0]), n_days=len(d), n_series=d.key.nunique(), n_fam=d.fam.nunique(), n_iv=d.prev_major.nunique())
    if boot:
        uc = d.prev_major.unique(); cl = d.prev_major.to_numpy(); bs = []
        for k in range(NB):
            m = dict(zip(uc, rng.multinomial(len(uc), np.ones(len(uc)) / len(uc))))
            wm = w * np.array([m[c] for c in cl])
            Gm = within(d, [q] + cols, pd.Series(wm, index=d.index))
            bs.append(wls(Gm[q].to_numpy(), Gm[cols].to_numpy(), wm)[0])
        out["sig"] = float(np.nanstd(bs)); out["z"] = out["b"] / out["sig"] if out["sig"] > 0 else np.nan
    if ret_all:
        return out, G, cols, b
    return out


def scramble(d, q, groups, n=500, by="site", semi=False):
    d = d.copy(); b0 = fit(d, q, by=by, semi=semi, boot=False)["b"]; sc = []
    for k in range(n):
        e = d.copy()
        e["ty"] = e.groupby(groups).ty.transform(lambda x: rng.permutation(x.to_numpy()))
        sc.append(fit(e, q, by=by, semi=semi, boot=False)["b"])
    sc = np.array(sc)
    return dict(b=b0, scr_sd=float(sc.std()), scr_p=float((np.sum(np.abs(sc) >= abs(b0)) + 1) / (n + 1)))


def main():
    t0 = time.time()
    F = pd.read_parquet(f"{TB}/b_days.parquet")
    F["date"] = pd.to_datetime(F.date)
    F["site"] = [C.SITE.get(n, n) for n in F.netsta]
    doy = F.date.dt.dayofyear.to_numpy(); ph = 2 * np.pi * doy / 365.25
    F["cs"], F["sn"], F["cs2"], F["sn2"] = np.cos(ph), np.sin(ph), np.cos(2 * ph), np.sin(2 * ph)
    F["ty"] = F.prev / 365.25; F["year"] = F.date.dt.year
    F["quarter"] = F.date.dt.quarter
    fam_tab = pd.read_csv(f"{C.MTB}/families.csv") if os.path.exists(f"{C.MTB}/families.csv") else None
    B = F[F.betw & (F.prev < 2000)].copy()
    R = {}; rows = []
    q = "dvv_2-5"
    d = prep(B, q)
    # ---- headline (all families, every store series)
    h = fit(d, q); h.update(scramble(d, q, ["key"], n=300)); R["headline_all"] = h; log(f"[r2 {SUB}] headline all {h}")
    # ---- C1 families
    R["C1_with_front"] = fit(d[d.has_front], q) if d.has_front.any() else {}
    R["C1_without_front"] = fit(d[~d.has_front], q) if (~d.has_front).sum() > 30 else {}
    # per-family slopes on season-removed residuals
    o, G, cols, b = fit(d, q, boot=False, ret_all=True)
    res = G[q].to_numpy() - G[cols[1:]].to_numpy() @ b[1:]
    d["res"] = res; d["tyw"] = G["ty"].to_numpy()
    pf = []
    for f, g in d.groupby("fam"):
        if len(g) < 15 or g.ty.std() < 0.05:
            continue
        sl = float(np.sum(g.w * g.tyw * g.res) / np.sum(g.w * g.tyw ** 2))
        pf.append(dict(fam=f, slope=sl, n=len(g), has_front=bool(g.has_front.iloc[0]), ty_span=float(g.ty.max() - g.ty.min())))
    PF = pd.DataFrame(pf); PF.to_csv(f"{TB}/r2_per_family.csv", index=False)
    nneg = int((PF.slope < 0).sum()) if len(PF) else 0
    R["C1_per_family"] = dict(n=len(PF), n_neg=nneg, frac_neg=nneg / max(len(PF), 1),
                              binom_p=float(binomtest(nneg, len(PF), 0.5, alternative="greater").pvalue) if len(PF) else np.nan)
    # ---- C2 source / family mix
    both = d.groupby("fam").prev.agg(lambda x: (x < 182).any() and (x >= 182).any())
    R["C2a_balanced"] = fit(d[d.fam.isin(both[both].index)], q)
    dd = d[np.isfinite(d.dSP_ms) & (np.abs(d.dSP_ms) < 140)].copy()
    R["C2b_with_dSP"] = fit(dd, q, extra_cols=("dSP_ms",))
    R["C2b_same_days_no_dSP"] = fit(dd, q)
    ds = prep(B, "dSP_ms", extra=lambda x: np.abs(x.dSP_ms) < 140)
    R["C2c_dSP_slope_ms_per_yr"] = fit(ds, "dSP_ms")
    if fam_tab is not None and {"lat", "lon"} <= set(fam_tab.columns):
        pos = fam_tab.set_index("fam")[["lat", "lon"]]
        e = d.join(pos, on="fam")
        for c in ("lat", "lon"):
            x = e[c].to_numpy(); t = e.ty.to_numpy(); w = e.w.to_numpy(); ok = np.isfinite(x)
            sl = np.cov(t[ok], x[ok], aweights=w[ok])[0, 1] / np.cov(t[ok], t[ok], aweights=w[ok])[0, 1]
            bs = []
            uc = e.prev_major.unique(); cl = e.prev_major.to_numpy()
            for k in range(NB):
                m = dict(zip(uc, rng.multinomial(len(uc), np.ones(len(uc)) / len(uc))))
                wm = w * np.array([m[cc] for cc in cl]); o2 = ok & (wm > 0)
                bs.append(np.cov(t[o2], x[o2], aweights=wm[o2])[0, 1] / np.cov(t[o2], t[o2], aweights=wm[o2])[0, 1])
            R[f"C2c_family_{c}_drift_deg_per_yr"] = dict(b=float(sl), sig=float(np.nanstd(bs)), km_per_yr=float(sl * 111.19 * (np.cos(np.radians(47.5)) if c == "lon" else 1)))
    # family-mix by early/late: share of families
    tab = pd.crosstab(d.fam, d.prev >= 182)
    R["C2c_family_mix_n_fam_early_only"] = int(((tab.get(True, 0) == 0)).sum()); R["C2c_family_mix_n_fam_late_only"] = int(((tab.get(False, 0) == 0)).sum())
    # ---- C3 stations
    bh = d.netsta.str.startswith("PB.")
    R["C3_borehole"] = fit(d[bh], q) if bh.sum() > 30 else {}
    R["C3_surface"] = fit(d[~bh], q) if (~bh).sum() > 30 else {}
    ps = []
    for s, g in d.groupby("netsta"):
        if len(g) >= 20 and g.ty.std() > 0.05:
            r = fit(g, q, by="netsta", boot=False); r["netsta"] = s; ps.append(r)
    pd.DataFrame(ps).to_csv(f"{TB}/r2_per_station.csv", index=False)
    # ---- C4 season
    R["C4_station_annual_semi"] = fit(d, q, by="netsta", semi=True)
    R["C4_season_matched_scramble"] = scramble(d, q, ["key", "quarter"], n=500)
    R["C4_no_annual"] = {}
    # ---- lapse
    lap = {}
    for wn in WL:
        dw = prep(B, f"dvv_{wn}")
        lap[wn] = fit(dw, f"dvv_{wn}")
    R["lapse"] = lap
    y = np.array([lap[w]["b"] for w in WL[:4]]); s = np.array([lap[w]["sig"] for w in WL[:4]])
    pat = {}
    for nm, p in PAT.items():
        p = np.array(p); k = np.sum(y * p / s ** 2) / np.sum(p ** 2 / s ** 2)
        pat[nm] = dict(k=float(k), chi2=float(np.sum((y - k * p) ** 2 / s ** 2)), implied_change_pct_per_yr=float(k))
    R["lapse_patterns"] = pat
    # ---- shape
    shp = []
    dS = prep(F[(F.prev >= 21) & (F.prev < 2000) & (F.nxt > 30)], q)
    for lab, dz in (("betw (>=60 d)", d), ("with 21-60 d days", dS)):
        for tau in (20, 30, 45, 60, 90, 120, 180, 270, 365, 540):
            dz = dz.copy(); dz["ex"] = np.exp(-dz.prev / tau)
            o, G, cols, b = fit(dz, q, reg=("ex",), boot=False, ret_all=True)
            r_ = G[q].to_numpy() - G[cols].to_numpy() @ b
            sse = float(np.sum(dz.w.to_numpy() * r_ ** 2)); n = len(dz)
            shp.append(dict(set=lab, form="exp", tau=tau, amp=float(b[0]), wsse=sse, n=n))
        o, G, cols, b = fit(dz, q, boot=False, ret_all=True)
        r_ = G[q].to_numpy() - G[cols].to_numpy() @ b
        shp.append(dict(set=lab, form="linear", tau=np.nan, amp=float(b[0]), wsse=float(np.sum(dz.w.to_numpy() * r_ ** 2)), n=len(dz)))
    SH = pd.DataFrame(shp)
    SH["aic"] = SH.n * np.log(SH.wsse / SH.n) + 2 * np.where(SH.form == "exp", 2, 1)
    SH.to_csv(f"{TB}/r2_shape.csv", index=False)
    for lab, g in SH.groupby("set"):
        be = g[g.form == "exp"].sort_values("aic").iloc[0]; li = g[g.form == "linear"].iloc[0]
        R[f"shape_{lab}"] = dict(best_tau=float(be.tau), amp=float(be.amp), dAIC_exp_minus_linear=float(be.aic - li.aic), linear_slope=float(li.amp))
    # ---- era split
    if SUB == "svi":
        for era, m in (("2003-2013", d.date.dt.year <= 2013), ("2014-2026", d.date.dt.year >= 2014)):
            R[f"era_{era}"] = fit(d[m], q) if m.sum() > 30 else {}
    # ---- tie-in: early level per interval vs size of the preceding major
    sp = f"{C.ROOT}/svi_science/results/s11_depth/round6/{SUB}/tables/r6_sizes.csv"
    if os.path.exists(sp):
        S = pd.read_csv(sp).set_index("episode")
        e = d[(d.prev >= 60) & (d.prev <= 180)]
        lev = e.groupby("prev_major").apply(lambda g: pd.Series(dict(v=np.sum(g.w * g.res) / g.w.sum(), w=g.w.sum(), n=len(g))),
                                            include_groups=False)
        lev = lev[lev.n >= 5].join(S[[c for c in ("Du", "T", "Fs") if c in S]])
        tie = {}
        for c in ("Du", "T", "Fs"):
            g = lev[np.isfinite(lev[c])] if c in lev else pd.DataFrame()
            if len(g) >= 5:
                X = np.column_stack([np.ones(len(g)), g[c]]); w = g.w.to_numpy()
                bb = wls(g.v.to_numpy(), X, w)[1]
                bs = [wls(g.v.to_numpy()[i], X[i], w[i])[1] for i in (rng.integers(0, len(g), len(g)) for _ in range(500))]
                tie[c] = dict(slope=float(bb), sig=float(np.nanstd(bs)), n_intervals=len(g))
        lev.to_csv(f"{TB}/r2_interval_levels.csv")
        R["tie_in_size"] = tie
    # ---- binned curves for the figure
    cur = []
    for lab, m in (("all", np.ones(len(d), bool)), ("with front (11 of §6)", d.has_front.to_numpy()), ("without front", ~d.has_front.to_numpy()),
                   ("borehole", bh.to_numpy()), ("surface", ~bh.to_numpy())):
        mo = d.prev.to_numpy() / 30.44
        for a, b_ in ((2, 4), (4, 6), (6, 9), (9, 12), (12, 18), (18, 30)):
            mm = m & (mo >= a) & (mo < b_)
            if mm.sum() >= 8:
                cur.append(dict(set=lab, mid=(a + b_) / 2, value=float(np.sum(d.w[mm] * d.res[mm]) / d.w[mm].sum()), weight=float(d.w[mm].sum()), n=int(mm.sum())))
    pd.DataFrame(cur).to_csv(f"{TB}/r2_curve.csv", index=False)
    # ---- verdict (PLAN §7)
    c1 = (R["C1_without_front"].get("z", 0) <= -2) or (R["C1_per_family"]["frac_neg"] >= 0.65 and R["C1_per_family"]["binom_p"] < 0.05)
    hb = R["headline_all"]["b"]
    c2 = (np.sign(R["C2a_balanced"]["b"]) == np.sign(hb) and abs(R["C2a_balanced"]["b"]) >= 0.5 * abs(hb)
          and abs(R["C2b_with_dSP"]["b"]) >= 0.5 * abs(R["C2b_same_days_no_dSP"]["b"]) and np.sign(R["C2b_with_dSP"]["b"]) == np.sign(hb))
    c3 = bool(R["C3_borehole"] and R["C3_surface"] and np.sign(R["C3_borehole"]["b"]) == np.sign(R["C3_surface"]["b"]) == np.sign(hb))
    c4 = R["C4_station_annual_semi"]["z"] <= -2 and R["C4_season_matched_scramble"]["scr_p"] < 0.01
    R["verdict"] = dict(C1=bool(c1), C2=bool(c2), C3=bool(c3), C4=bool(c4), confirmed=bool(c1 and c2 and c3 and c4 and hb < 0),
                        written=C.pt_now() if hasattr(C, "pt_now") else time.strftime("%Y-%m-%d %H:%M"))
    json.dump(R, open(f"{TB}/r2_verdict.json", "w"), indent=1, default=float)
    log(f"[r2 {SUB}] verdict {R['verdict']} ({time.time() - t0:.0f}s)\n" + json.dumps(R, indent=1, default=float)[:6000])
    open(f"{OUT}/R2_{SUB.upper()}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
