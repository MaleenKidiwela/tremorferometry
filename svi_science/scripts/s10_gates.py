#!/usr/bin/env python
"""s10_gates -- S10 pilot step 4 (PLAN §5-7, §9): gates table, precision vs number of repeats, change tests,
POLARIS coverage / agreement and the KELB->KLNB check. Reads the measure tables; writes tables/*.csv and result.json.

    python s10_gates.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_measure as M

log = S.Log(f"{S.LOGD}/gates.log")
RES = {}


def theil_sen(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 5:
        return np.nan, np.nan, np.nan
    r = stats.theilslopes(y[ok], x[ok])
    return float(r[0]), float(r[2]), float(r[3])


def wfit(t, y, s):
    """weighted LS slope of y on t with sigma scaled by sqrt(chi2_red) when > 1"""
    ok = np.isfinite(t) & np.isfinite(y) & np.isfinite(s) & (s > 0)
    if ok.sum() < 4:
        return np.nan, np.nan, np.nan
    t, y, w = t[ok], y[ok], 1 / s[ok] ** 2
    A = np.c_[np.ones_like(t), t]
    N = A.T @ (A * w[:, None]); b = A.T @ (w * y)
    if np.linalg.cond(N) > 1e12:
        return np.nan, np.nan, np.nan
    p = np.linalg.solve(N, b); r = y - A @ p
    chi2 = float(np.sum(w * r ** 2) / max(ok.sum() - 2, 1))
    cov = np.linalg.inv(N) * max(chi2, 1.0)
    return float(p[1]), float(np.sqrt(cov[1, 1])), chi2


def gates(D, ser, INJ):
    rows = []
    for (fam, ns, ep), g in D.groupby(["fam", "netsta", "epoch"]):
        s_ = ser[(ser.fam == fam) & (ser.netsta == ns) & (ser.epoch == ep)].iloc[0]
        r = dict(fam=fam, netsta=ns, epoch=ep, n_days=len(g), n_win_median=float(g.n.median()),
                 sig24_median=float(g.sig24.median()), sig35_median=float(g.sig35.median()),
                 g_sig24_median=float(g.g_sig24.median()), cc24_median=float(g.cc24.median()),
                 date0=g.date.min(), date1=g.date.max(), rail24_frac=float(g.rail24.mean()),
                 n_days_measured=int(g.eps24.notna().sum()))
        r["causal24"], r["causal35"] = s_.causal24, s_.causal35
        r["a_pass"] = bool(s_.causal24 > 1.5) if np.isfinite(s_.causal24) else None
        for pre, nm in (("c_", "b1"), ("c2_", "b2")):
            z = (g[f"{pre}eps24"] / g[f"{pre}sig24"]).to_numpy() if f"{pre}eps24" in g else np.array([])
            z = z[np.isfinite(z)]
            r[f"{nm}_n"] = len(z)
            r[f"{nm}_frac_absz_lt2"] = float(np.mean(np.abs(z) < 2)) if len(z) else np.nan
            r[f"{nm}_mean_z"] = float(np.mean(z)) if len(z) else np.nan
            r[f"{nm}_cc24_median"] = float(g[f"{pre}cc24"].median()) if f"{pre}cc24" in g else np.nan
            r[f"{nm}_sd_eps24"] = float(g[f"{pre}eps24"].std()) if f"{pre}eps24" in g else np.nan
        mm = np.isfinite(g.mir24) & np.isfinite(g.eps24)
        r["mirror_frac_lt_coda"] = float(np.mean(np.abs(g.mir24[mm]) < np.abs(g.eps24[mm]))) if mm.any() else np.nan
        r["mirror_sd24"] = float(g.mir24.std()); r["data_sd24"] = float(g.eps24.std())
        r["b_pass"] = (None if (r["b2_n"] < 3 or not np.isfinite(r["mirror_frac_lt_coda"])) else
                       bool(r["b2_frac_absz_lt2"] >= 0.9 and abs(r["b2_mean_z"]) <= 0.5 and r["mirror_frac_lt_coda"] >= 0.7))
        r["b_pass_ctrl_only"] = None if r["b2_n"] < 3 else bool(r["b2_frac_absz_lt2"] >= 0.9 and abs(r["b2_mean_z"]) <= 0.5)
        # (c)
        if "hi_eps24" in g:
            dz = (g.hi_eps24 - g.eps24) / np.sqrt(np.abs(g.hi_sig24 ** 2 - g.sig24 ** 2))
            dz = dz[np.isfinite(dz)]
            r["c_n"] = len(dz); r["c_frac_absz_lt2"] = float(np.mean(np.abs(dz) < 2)) if len(dz) else np.nan
            r["c_pass"] = None if len(dz) < 3 else bool(r["c_frac_absz_lt2"] >= 0.9)
        # (d)
        sl, lo, hi = theil_sen(g.eps24.to_numpy(), g.eps35.to_numpy())
        sg = np.sign(g.eps24) == np.sign(g.eps35)
        r.update(d_slope=sl, d_slope_lo=lo, d_slope_hi=hi, d_sign_agree=float(sg.mean()))
        big = np.abs(g.eps24 / g.sig24) > 2
        r["d_sign_agree_z2"] = float(sg[big].mean()) if big.any() else np.nan
        r["d_n_z2"] = int(big.sum())
        r["d_pass"] = None if not np.isfinite(sl) else bool(sl >= 0.6 and r["d_sign_agree"] >= 0.7)
        # (e)
        I = INJ[(INJ.fam == fam) & (INJ.netsta == ns) & (INJ.epoch == ep) & (INJ.win == 24)]
        for est in ("daily", "graph"):
            J = I[I.est == est]
            if len(J):
                r[f"e_{est}_rec_min"] = float(J.recovery.min()); r[f"e_{est}_rec_max"] = float(J.recovery.max())
                r[f"e_{est}_leak_maxabs"] = float(J.leakage.abs().max())
                ok = J.recovery.between(0.8, 1.2).all() and (J.leakage.abs() < 0.2).all()
                r[f"e_{est}_pass"] = bool(ok) if J.leakage.notna().all() else None
        r["e_pass"] = (None if r.get("e_daily_pass") is None or r.get("e_graph_pass") is None
                       else bool(r["e_daily_pass"] and r["e_graph_pass"]))
        r["prepost_cc24"] = s_.get("prepost_cc24", np.nan)
        # (f)
        inside = [f"{d} {lab}" for d, lab in S.CHANGES if g.date.min() < d < g.date.max()]
        r["f_overlay"] = "drawn"; r["f_changes_inside"] = "; ".join(inside)
        ps = [r.get(k) for k in ("a_pass", "b_pass", "c_pass", "d_pass", "e_pass")]
        r["n_gates_pass"] = int(sum(p is True for p in ps)); r["n_gates_fail"] = int(sum(p is False for p in ps))
        r["n_gates_na"] = int(sum(p is None for p in ps))
        r["all_pass"] = bool(all(p is True for p in ps))
        rows.append(r)
    return pd.DataFrame(rows)


def precision(P):
    """median sigma(n) over busy days, per window and station kind; fit sigma = c n^b on the series-days that reach
    n >= 100 (so every n level is the same set of days), and per-day fits summarised by their median."""
    out = {}
    P = P.copy(); P["key"] = P.fam.astype(str) + P.netsta + P.epoch + P.day
    full = P[P.n == 100].key.unique()
    P = P[P.key.isin(full)]
    P = P.copy(); P["kind"] = np.where(P.netsta.str.startswith("PB."), "borehole", "surface")
    for (win, kind), g in list(P.groupby(["win", "kind"])) + [((w, "all"), g) for w, g in P.groupby("win")]:
        med = g.groupby("n").sigma.median()
        n = med.index.to_numpy(float); s = med.to_numpy()
        b, a = np.polyfit(np.log(n), np.log(s), 1)
        need = {lvl: float(np.exp((np.log(lvl) - a) / b)) for lvl in (0.001, 0.0005)}
        rf = g.groupby("n").rail_frac.mean()
        okn = rf.index[rf < 0.05].to_numpy(float)
        if len(okn) >= 2:
            b2, a2 = np.polyfit(np.log(okn), np.log(med.loc[okn].to_numpy()), 1)
            need2 = {lvl: float(np.exp((np.log(lvl) - a2) / b2)) for lvl in (0.001, 0.0005)}
        else:
            b2, need2 = np.nan, {0.001: np.nan, 0.0005: np.nan}
        out[f"{win}_{kind}"] = dict(sigma_by_n={int(k): float(v) for k, v in med.items()},
                                    rail_frac_by_n={int(k): float(v) for k, v in rf.items()},
                                    unrailed_fit_slope=float(b2), unrailed_n_for_0p1pct=need2[0.001],
                                    unrailed_n_for_0p05pct=need2[0.0005],
                                    n_series=int(g.groupby(["fam", "netsta", "epoch"]).ngroups), slope=float(b),
                                    n_for_0p1pct=need[0.001], n_for_0p05pct=need[0.0005])
    return out


def front_dates():
    X = pd.read_parquet(f"{S.TB}/rows.parquet")
    E = pd.read_csv(f"{S.ROOT}/svi_science/results/r23/tables/episodes_2003_2026.csv")
    H = X[(X.series == "data") & X.hiP]
    rows = []
    for e in E.itertuples():
        a = (pd.Timestamp(e.start) - pd.Timedelta(days=2)).strftime("%Y-%m-%d")
        b = (pd.Timestamp(e.end) + pd.Timedelta(days=2)).strftime("%Y-%m-%d")
        for f in S.FAMS:
            h = H[(H.fam == f) & (H.date >= a) & (H.date <= b)]
            if len(h):
                rows.append(dict(episode=e.episode, fam=f, front=h.date.min(), big=bool(e.big), major=bool(e.major)))
    return pd.DataFrame(rows)


def within_episode(D, FR):
    D = D.merge(FR, left_on=["cluster", "fam"], right_on=["episode", "fam"], how="inner")
    D["dsf"] = (pd.to_datetime(D.date) - pd.to_datetime(D.front)).dt.days
    out, pooled = [], []
    for key, g in D.groupby(["fam", "netsta", "epoch", "cluster"]):
        if len(g) < 2:
            continue
        for pre in ("", "g_", "c2_"):
            e, s = f"{pre}eps24", f"{pre}sig24"
            if e not in g:
                continue
            ok = np.isfinite(g[e]) & np.isfinite(g[s])
            if ok.sum() < 2:
                continue
            w = 1 / g[s][ok] ** 2
            dm = g[e][ok] - np.sum(w * g[e][ok]) / w.sum()
            pooled.append(pd.DataFrame(dict(fam=key[0], netsta=key[1], epoch=key[2], cluster=key[3], est=pre or "d_",
                                            dsf=g.dsf[ok], de=dm, sig=g[s][ok])))
    P = pd.concat(pooled, ignore_index=True)
    for (fam, est), g in list(P.groupby(["fam", "est"])) + [(("all", est), g) for est, g in P.groupby("est")]:
        sl, se, chi2 = wfit(g.dsf.to_numpy(float), g.de.to_numpy(), g.sig.to_numpy())
        out.append(dict(fam=fam, est=est, n=len(g), slope_pct_per_day=100 * sl, se=100 * se, z=sl / se, chi2=chi2,
                        early_minus_late_pct=100 * (np.average(g.de[g.dsf <= 3], weights=1 / g.sig[g.dsf <= 3] ** 2)
                                                    - np.average(g.de[g.dsf > 3], weights=1 / g.sig[g.dsf > 3] ** 2))
                        if ((g.dsf <= 3).any() and (g.dsf > 3).any()) else np.nan))
    return pd.DataFrame(out), P


def longterm(D, Ep):
    rows = []
    for (fam, ns, ep), g in D.groupby(["fam", "netsta", "epoch"]):
        t = (pd.to_datetime(g.date) - pd.Timestamp("2000-01-01")).dt.days.to_numpy() / 365.25
        if t.max() - t.min() < 2:
            continue
        r = dict(fam=fam, netsta=ns, epoch=ep, span_yr=float(t.max() - t.min()), n=len(g))
        for pre in ("", "g_", "c2_"):
            e, s = f"{pre}eps24", f"{pre}sig24"
            if e in g:
                sl, se, chi2 = wfit(t, g[e].to_numpy(), g[s].to_numpy())
                r[f"{pre or 'd_'}slope_pct_per_yr"] = 100 * sl; r[f"{pre or 'd_'}se"] = 100 * se
                r[f"{pre or 'd_'}z"] = sl / se if se else np.nan
        sl, lo, hi = theil_sen(t, g.mir24.to_numpy())
        r["mirror_theilsen_pct_per_yr"] = 100 * sl if np.isfinite(sl) else np.nan
        e_ = Ep[(Ep.fam == fam) & (Ep.netsta == ns) & (Ep.epoch == ep)]
        if len(e_) >= 4:
            te = (pd.to_datetime(e_.date_mid) - pd.Timestamp("2000-01-01")).dt.days.to_numpy() / 365.25
            sl, se, chi2 = wfit(te, e_.eps24.to_numpy(), e_.sig24.to_numpy())
            r.update(ep_slope_pct_per_yr=100 * sl, ep_se=100 * se, ep_chi2=chi2, ep_n=len(e_))
        rows.append(r)
    return pd.DataFrame(rows)


def polaris(D):
    rows, agr = [], []
    ref = ["CN.PGC", "CN.VGZ", "CN.SNB"]
    for fam in S.FAMS:
        Df = D[(D.fam == fam) & (D.date < "2007-01-01")]
        refdays = set(Df[Df.netsta.isin(ref)].date)
        alldays = set(Df.date)
        for ns in S.POLARIS:
            g = Df[Df.netsta == ns]
            rows.append(dict(fam=fam, netsta=ns, days_with_stack=g.date.nunique(),
                             days_without_PGC_VGZ_SNB=len(set(g.date) - refdays),
                             family_days_2003_06=len({d for d in alldays}), epochs=",".join(sorted(g.epoch.unique()))))
            for (ep, ), gp in g.groupby(["epoch"]):
                for rn in ref:
                    for (ep2, ), gr in Df[Df.netsta == rn].groupby(["epoch"]):
                        m = gp.merge(gr, on="date", suffixes=("_a", "_b"))
                        m = m[np.isfinite(m.eps24_a) & np.isfinite(m.eps24_b)]
                        if len(m) < 4:
                            continue
                        s2 = m.sig24_a ** 2 + m.sig24_b ** 2; w = 1 / s2
                        dd = m.eps24_a - m.eps24_b; off = np.sum(w * dd) / w.sum()
                        chi2 = float(np.sum(w * (dd - off) ** 2)); dof = len(m) - 1
                        wc = np.cov(m.eps24_a, m.eps24_b, aweights=1 / (m.sig24_a * m.sig24_b))
                        agr.append(dict(fam=fam, polaris=f"{ns}:{ep}", ref=f"{rn}:{ep2}", n_common=len(m), chi2=chi2,
                                        dof=dof, p_chi2=float(stats.chi2.sf(chi2, dof)), agree=bool(stats.chi2.sf(chi2, dof) > 0.01),
                                        wcorr=float(wc[0, 1] / np.sqrt(wc[0, 0] * wc[1, 1])),
                                        sd_diff_pct=100 * float(np.std(dd - off)), med_sig_pair_pct=100 * float(np.median(np.sqrt(s2)))))
    return pd.DataFrame(rows), pd.DataFrame(agr)


def core_agreement(D):
    """same-day agreement between core station series (ε 2-4 s): chi2 after removing the mean offset, weighted corr"""
    out = []
    for fam in S.FAMS:
        Df = D[(D.fam == fam) & D.netsta.isin(S.CORE)]
        keys = sorted(Df.groupby(["netsta", "epoch"]).groups)
        for i, (na, ea) in enumerate(keys):
            for nb, eb in keys[i + 1:]:
                if na == nb:
                    continue
                m = Df[(Df.netsta == na) & (Df.epoch == ea)].merge(Df[(Df.netsta == nb) & (Df.epoch == eb)], on="date", suffixes=("_a", "_b"))
                m = m[np.isfinite(m.eps24_a) & np.isfinite(m.eps24_b)]
                if len(m) < 8:
                    continue
                s2 = m.sig24_a ** 2 + m.sig24_b ** 2; w = 1 / s2; dd = m.eps24_a - m.eps24_b
                off = np.sum(w * dd) / w.sum(); chi2 = float(np.sum(w * (dd - off) ** 2)); dof = len(m) - 1
                r, p = stats.spearmanr(m.eps24_a, m.eps24_b)
                out.append(dict(fam=fam, a=f"{na}:{ea}", b=f"{nb}:{eb}", n_common=len(m), chi2_red=chi2 / dof,
                                p_chi2=float(stats.chi2.sf(chi2, dof)), spearman=float(r), p_spearman=float(p),
                                sd_a_pct=100 * float(m.eps24_a.std()), sd_b_pct=100 * float(m.eps24_b.std()),
                                med_sig_a_pct=100 * float(m.sig24_a.median()), med_sig_b_pct=100 * float(m.sig24_b.median())))
    return pd.DataFrame(out)


def kelb_check():
    """PLAN §9: KELB(b) vs KLNB gain-normalised check per family (Z)."""
    X = pd.read_parquet(f"{S.TB}/rows.parquet"); RS = pd.read_parquet(f"{S.TB}/row_station.parquet")
    ser = pd.read_csv(f"{S.TB}/series.csv")
    rng = np.random.default_rng(7)
    out = []
    G = {}
    for ns, ep in (("PO.KELB", "Eb"), ("PO.KLNB", "E1")):
        z = np.load(f"{S.WIN}/{ns}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}
        R = RS[(RS.netsta == ns) & (RS.epoch == ep)].merge(X, on=["row_id", "date"])
        R = R[(R.series == "data") & R.row_id.isin(pos)]
        for fam, g in R.groupby("fam"):
            W = z["W"][[pos[r] for r in g.row_id]]; ok = np.isfinite(W[:, 0, 0]); g = g[ok]; W = W[ok]
            if len(g) < 20:
                continue
            P = g.P.to_numpy()
            G0 = np.stack([M.wmean(W[:, c], P)[0] for c in range(3)]); tS0 = M.pick_S(G0)
            W, _, _ = M.realign(W, G0, tS0)
            G1 = (P @ W[:, 0]) / P.sum(); tS = M.pick_S(np.stack([M.wmean(W[:, c], P)[0] for c in range(3)]), tS0)
            hs = []
            for _ in range(50):
                pp = rng.permutation(len(g)); h = len(g) // 2
                A_, B_ = (P[pp[:h]] @ W[pp[:h], 0]) / P[pp[:h]].sum(), (P[pp[h:]] @ W[pp[h:], 0]) / P[pp[h:]].sum()
                hs.append(M.Ref(A_, tS, (2, 4)).measure(B_)[1][0])
            G[(ns, fam)] = (G1, tS, W[:, 0], P, float(np.median(hs)))
    for fam in S.FAMS:
        if ("PO.KELB", fam) not in G or ("PO.KLNB", fam) not in G:
            continue
        Ga, ta, Wa, Pa, hsa = G[("PO.KELB", fam)]; Gb, tb, Wb, Pb, hsb = G[("PO.KLNB", fam)]
        # put KELB on the KLNB S time (sub-sample shift)
        sh = ta - tb
        f = np.fft.rfftfreq(S.NW, 1 / S.FS)
        shift = lambda x: np.fft.irfft(np.fft.rfft(x) * np.exp(2j * np.pi * f * sh), n=S.NW)
        R = M.Ref(Gb, tb, (2, 4))
        e, cc = R.measure(shift(Ga))
        boots = []
        for _ in range(200):
            ia = rng.integers(0, len(Pa), len(Pa)); ib = rng.integers(0, len(Pb), len(Pb))
            ga = (Pa[ia] @ Wa[ia]) / Pa[ia].sum(); gb = (Pb[ib] @ Wb[ib]) / Pb[ib].sum()
            boots.append(M.Ref(gb, tb, (2, 4)).measure(shift(ga))[0][0])
        sig = float(np.std(boots, ddof=1))
        sa = ser[(ser.fam == fam) & (ser.netsta == "PO.KELB") & (ser.epoch == "Eb")]
        sb = ser[(ser.fam == fam) & (ser.netsta == "PO.KLNB") & (ser.epoch == "E1")]
        ra = float((sa.rms_coda24 / sa.rms_S).iloc[0]) if len(sa) else np.nan
        rb = float((sb.rms_coda24 / sb.rms_S).iloc[0]) if len(sb) else np.nan
        crit_cc = min(hsa, hsb) - 0.05
        merge = bool(cc[0] >= crit_cc and abs(e[0]) < 2 * sig and abs(ra / rb - 1) < 0.2)
        out.append(dict(fam=fam, cc_kelb_klnb_24=float(cc[0]), halfsplit_cc_kelb=hsa, halfsplit_cc_klnb=hsb,
                        cc_criterion=crit_cc, stretch_pct=100 * float(e[0]), stretch_sig_pct=100 * sig,
                        coda_S_ratio_kelb=ra, coda_S_ratio_klnb=rb, ratio_of_ratios=ra / rb, merge=merge,
                        n_kelb=len(Pa), n_klnb=len(Pb)))
    return pd.DataFrame(out)


def main():
    D = pd.read_parquet(f"{S.TB}/daily.parquet")
    Ep = pd.read_parquet(f"{S.TB}/episodes.parquet")
    INJ = pd.read_csv(f"{S.TB}/injection.csv"); INJ["win"] = INJ.win.astype(int)
    ser = pd.read_csv(f"{S.TB}/series.csv")
    P = pd.read_csv(f"{S.TB}/precision.csv")
    GT = gates(D, ser, INJ)
    GT.to_csv(f"{S.TB}/gates.csv", index=False)
    log("gates per series:\n" + GT[["fam", "netsta", "epoch", "n_days", "sig24_median", "causal24", "a_pass", "b2_frac_absz_lt2",
                                   "b2_mean_z", "mirror_frac_lt_coda", "b_pass", "c_pass", "d_slope", "d_sign_agree", "d_pass",
                                   "e_daily_rec_min", "e_graph_rec_min", "e_pass", "all_pass"]].round(3).to_string())
    RES["precision"] = precision(P)
    log("precision: " + json.dumps(RES["precision"], indent=1))
    FR = front_dates(); FR.to_csv(f"{S.TB}/front_dates.csv", index=False)
    WE, WP = within_episode(D, FR); WE.to_csv(f"{S.TB}/within_episode.csv", index=False); WP.to_parquet(f"{S.TB}/within_episode_points.parquet")
    log("within-episode:\n" + WE.round(4).to_string())
    LT = longterm(D, Ep); LT.to_csv(f"{S.TB}/longterm.csv", index=False)
    log("long-term:\n" + LT.round(4).to_string())
    PC, PA = polaris(D); PC.to_csv(f"{S.TB}/polaris_coverage.csv", index=False); PA.to_csv(f"{S.TB}/polaris_agreement.csv", index=False)
    log("POLARIS coverage:\n" + PC.to_string()); log("POLARIS agreement:\n" + PA.round(3).to_string())
    CA = core_agreement(D); CA.to_csv(f"{S.TB}/core_agreement.csv", index=False)
    log("core same-day agreement:\n" + CA.round(3).to_string())
    KB = kelb_check(); KB.to_csv(f"{S.TB}/kelb_klnb_check.csv", index=False)
    log("KELB/KLNB:\n" + KB.round(4).to_string())
    core = GT[GT.netsta.isin(S.CORE)]
    RES["gates_summary"] = {k: dict(core_pass=int((core[k] == True).sum()), core_fail=int((core[k] == False).sum()),
                                    all_pass=int((GT[k] == True).sum()), all_fail=int((GT[k] == False).sum()))
                            for k in ("a_pass", "b_pass", "b_pass_ctrl_only", "c_pass", "d_pass", "e_daily_pass", "e_graph_pass", "e_pass", "all_pass")}
    RES["daily_sigma24_median_core"] = float(D[D.netsta.isin(S.CORE)].sig24.median())
    RES["daily_n_median_core"] = float(D[D.netsta.isin(S.CORE)].n.median())
    S.Log  # noqa
    json.dump(RES, open(f"{S.TB}/result.json", "w"), indent=1, default=str)
    log(json.dumps(RES["gates_summary"], indent=1))
    open(f"{S.OUT}/GATES_DONE", "w").write("done\n")


if __name__ == "__main__":
    main()
