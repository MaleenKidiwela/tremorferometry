#!/usr/bin/env python
"""s12_epoch_pool -- pool the per-series superposed-epoch values of s12_epoch.py (PLAN.md §3-§4) and apply the rules.
contrast = class - betw (same series), per draw (draw 0 = data, 1..200 = joint bootstrap over majors / inter-ETS
intervals). Pooled = inverse-variance mean over series (sigma floor = 25th percentile of series sigmas); pooled sigma =
std over the joint draws. Also per site, odd/even families, the front-scrambled control (pre - 0..20) and the time-of-day
scrambled control (0..20 - betw). For coda dv/v, the LVL change implied by the data: delta = sum w a x / sum w a^2 with
a = the series' forward LVL coefficient (ell40, nominal slab; SVI) or the SVI median a (WA, labelled).
    S11_REGION=SVI|WA_NW python s12_epoch_pool.py [--test]
"""
from __future__ import annotations
import os, sys, json, pickle, time, argparse
import numpy as np
import pandas as pd
ROOT = "/home/jovyan/tremorferometry"
sys.path.insert(0, f"{ROOT}/svi_science/scripts"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s11_common as C
from s12_epoch import OUT, SUB, TB, CLASSES, log

OBS = ["dvv_1-2", "dvv_2-3", "dvv_3-4", "dvv_4-6", "dvv_2-5", "day_lnE_1-2", "day_lnE_2-3", "day_lnE_3-4", "day_lnE_4-6",
       "day_slope", "lnPS_Z", "lnPS_H", "dSP_ms", "dSrS_ms"]
DELTA = 0.02


def a_lvl_table():
    K = pd.read_parquet(f"{OUT}/tables/fwd_coef.parquet")
    K = K[(K.variant == "ell40") & (K.comp == "3c")]
    return {(r.series, r.win): r.a_lvl_nom for r in K.itertuples()}, K.groupby("win").a_lvl_nom.median().to_dict()


def pool(X, S, keep):
    """X (n_series, n_draw) contrasts, S (n_series,) sigmas -> pooled per draw"""
    s = np.where(keep, S, np.nan)
    fl = np.nanpercentile(s, 25) if np.isfinite(s).any() else np.nan
    w = np.where(keep, 1.0 / np.maximum(S, fl) ** 2, 0.0)
    W = w[:, None] * np.isfinite(X)
    return np.nansum(W * np.nan_to_num(X), 0) / np.maximum(W.sum(0), 1e-30), w


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--test", action="store_true"); a = ap.parse_args()
    tag = "_test" if a.test else ""
    D = pickle.load(open(f"{TB}/epoch_raw{tag}.pkl", "rb"))
    aT, aMed = a_lvl_table() if os.path.exists(f"{OUT}/tables/fwd_coef.parquet") else ({}, {})
    FOOT = {}
    if C.REGION == "SVI" and os.path.exists(f"{OUT}/tables/fwd_coef_foot.parquet"):
        Kf = pd.read_parquet(f"{OUT}/tables/fwd_coef_foot.parquet"); Kf = Kf[(Kf.variant == "ell40") & (Kf.comp == "3c")]
        for pl in ("foot15", "foot25"):
            FOOT[pl] = ({(r.series, r.win): getattr(r, f"a_lvl_{pl}") for r in Kf.itertuples()}, None)
    if C.REGION != "SVI":
        aT = {}
    keys = sorted(D)
    site = np.array([C.SITE.get(k.split("|")[1], k.split("|")[1]) for k in keys])
    fam = np.array([int(k.split("|")[0]) for k in keys])
    rows, srows = [], []
    for cl in [c for c in CLASSES if c != "betw"] + ["0..20"]:
        for q in OBS:
            X = np.full((len(keys), 201), np.nan)
            for i, k in enumerate(keys):
                r = D[k]["res"]
                if cl in r and "betw" in r and q in r[cl]:
                    X[i] = np.asarray(r[cl][q], float) - np.asarray(r["betw"][q], float)
            fin = np.isfinite(X[:, 1:]).mean(1)
            S = np.nanstd(X[:, 1:], 1)
            keep = np.isfinite(X[:, 0]) & (fin >= 0.5) & (S > 0)
            if keep.sum() < 3:
                continue
            P, w = pool(X, S, keep)
            row = dict(cls=cl, obs=q, value=float(P[0]), sigma=float(np.nanstd(P[1:])), n_series=int(keep.sum()),
                       n_fam=int(len(set(fam[keep]))), weight=float(w.sum()))
            for sname in sorted(set(site)):
                m = keep & (site == sname)
                if m.sum() >= 2:
                    Ps, _ = pool(X[m], S[m], np.ones(m.sum(), bool))
                    row[f"site_{sname}"] = float(Ps[0]); row[f"site_{sname}_sig"] = float(np.nanstd(Ps[1:]))
                    srows.append(dict(cls=cl, obs=q, site=sname, value=float(Ps[0]), sigma=float(np.nanstd(Ps[1:])),
                                      n_series=int(m.sum())))
            for par, m in (("odd", keep & (fam % 2 == 1)), ("even", keep & (fam % 2 == 0))):
                if m.sum() >= 2:
                    Pp, _ = pool(X[m], S[m], np.ones(m.sum(), bool))
                    row[f"{par}"] = float(Pp[0]); row[f"{par}_sig"] = float(np.nanstd(Pp[1:]))
            if q.startswith("dvv_") and aT:
                wn = q.split("_")[1]
                av = np.array([aT.get((k, wn), aMed.get(wn, np.nan)) for k in keys])
                ok = keep & np.isfinite(av)
                wa = w * ok
                num = np.nansum((wa * av)[:, None] * np.nan_to_num(X) * np.isfinite(X), 0)
                den = np.nansum(((wa * av ** 2)[:, None]) * np.isfinite(X), 0)
                dl = num / np.maximum(den, 1e-30)              # implied LVL dVs/Vs in % (X is in %)
                row["delta_lvl_pct"] = float(dl[0]); row["delta_lvl_sig"] = float(np.nanstd(dl[1:]))
                row["target_pct"] = 100 * DELTA * float(aMed.get(wn, np.nan))
                for pl, (aTp, aMp) in FOOT.items():                  # Gosselin footprint: LVL change only near PGC/SNB/VGZ
                    av = np.array([aTp.get((k, wn), np.nan) for k in keys]); ok = keep & np.isfinite(av)
                    wa = w * ok
                    num = np.nansum((wa * av)[:, None] * np.nan_to_num(X) * np.isfinite(X), 0)
                    den = np.nansum(((wa * av ** 2)[:, None]) * np.isfinite(X), 0)
                    dl = num / np.maximum(den, 1e-30)
                    row[f"delta_{pl}_pct"] = float(dl[0]); row[f"delta_{pl}_sig"] = float(np.nanstd(dl[1:]))
            rows.append(row)
    T = pd.DataFrame(rows); T.to_csv(f"{TB}/epoch_pooled{tag}.csv", index=False)
    pd.DataFrame(srows).to_csv(f"{TB}/epoch_sites{tag}.csv", index=False)
    # ---- front-scrambled control: pre - 0..20
    scr_rows = []
    for q in [o for o in OBS if not o.startswith("day_")]:
        Xd = np.full(len(keys), np.nan); Xs = None; Sg = np.full(len(keys), np.nan)
        for i, k in enumerate(keys):
            r = D[k]["res"]; sc = D[k]["scr"]
            if "pre" in r and "0..20" in r and q in r["pre"]:
                x = np.asarray(r["pre"][q]) - np.asarray(r["0..20"][q])
                Xd[i] = x[0]; Sg[i] = np.nanstd(x[1:])
                if sc and q in sc:
                    lab = sc["lab"]; v = np.asarray(sc[q])
                    nd = max(l[0] for l in lab) + 1
                    pre = np.array([v[j] for j, l in enumerate(lab) if l[1] == "pre"])
                    po = np.array([v[j] for j, l in enumerate(lab) if l[1] == "0..20"])
                    wpre = np.array([l[2] for l in lab if l[1] == "pre"])
                    xs = np.where(wpre > 0, pre - po, np.nan)
                    if Xs is None:
                        Xs = np.full((len(keys), nd), np.nan)
                    Xs[i, :len(xs)] = xs
        keep = np.isfinite(Xd) & np.isfinite(Sg) & (Sg > 0)
        if keep.sum() < 3:
            continue
        fl = np.percentile(Sg[keep], 25); w = np.where(keep, 1 / np.maximum(Sg, fl) ** 2, 0)
        val = float(np.sum(w[keep] * Xd[keep]) / w[keep].sum())
        row = dict(obs=q, pre_minus_020=val, n_series=int(keep.sum()))
        if Xs is not None:
            W = w[:, None] * np.isfinite(Xs)
            ps = np.nansum(W * np.nan_to_num(Xs), 0) / np.maximum(W.sum(0), 1e-30)
            row.update(scr_sd=float(np.std(ps)), scr_mean=float(np.mean(ps)),
                       scr_p=float((np.sum(np.abs(ps - np.mean(ps)) >= abs(val - np.mean(ps))) + 1) / (len(ps) + 1)))
        scr_rows.append(row)
    pd.DataFrame(scr_rows).to_csv(f"{TB}/epoch_frontscr{tag}.csv", index=False)
    # ---- time-of-day scrambled control: 0..20 - betw
    crow = []
    for q in [o for o in OBS if not o.startswith("day_")]:
        xs = []; ws = []
        for i, k in enumerate(keys):
            c = D[k]["ctrl"]
            if "0..20" in c and "betw" in c and np.isfinite(c["0..20"].get(q, np.nan)) and np.isfinite(c["betw"].get(q, np.nan)):
                xs.append(c["0..20"][q] - c["betw"][q])
        xs = np.array(xs)
        if len(xs) >= 3:
            crow.append(dict(obs=q, ctrl_mean=float(np.mean(xs)), ctrl_sig=float(np.std(xs, ddof=1) / np.sqrt(len(xs))),
                             n=len(xs)))
    pd.DataFrame(crow).to_csv(f"{TB}/epoch_ctrl{tag}.csv", index=False)
    # ---- verdicts (PLAN §3-§4)
    V = {}
    hd = T[T.cls == "0..20"].set_index("obs")
    C_ = pd.DataFrame(crow).set_index("obs") if crow else pd.DataFrame()
    sites = sorted(set(site))
    for q in hd.index:
        r = hd.loc[q]
        z = r.value / r.sigma if r.sigma > 0 else np.nan
        sgn_sites = [np.sign(r.get(f"site_{s}", np.nan)) for s in sites if np.isfinite(r.get(f"site_{s}", np.nan))]
        same = int(np.sum(np.array(sgn_sites) == np.sign(r.value)))
        ctrl_ok = (q in C_.index) and abs(C_.loc[q, "ctrl_mean"]) <= 2 * C_.loc[q, "ctrl_sig"] if len(C_) else None
        V[q] = dict(value=float(r.value), sigma=float(r.sigma), z=float(z), sites_same_sign=same, n_sites=len(sgn_sites),
                    ctrl_ok=None if ctrl_ok is None else bool(ctrl_ok), bound2=float(abs(r.value) + 2 * r.sigma))
        if q.startswith("dvv_") and "delta_lvl_pct" in r:
            V[q].update(delta_lvl_pct=float(r.delta_lvl_pct), delta_lvl_sig=float(r.delta_lvl_sig),
                        target_pct=float(r.target_pct), bound2_over_target=float((abs(r.value) + 2 * r.sigma) / r.target_pct))
        if q.startswith("dvv_"):
            V[q]["detect"] = bool(abs(z) >= 3 and r.value > 0 and same >= 3 and bool(ctrl_ok))
    V["_written"] = C.pt_now() if hasattr(C, "pt_now") else time.strftime("%Y-%m-%d %H:%M")
    json.dump(V, open(f"{TB}/epoch_verdict{tag}.json", "w"), indent=1, default=float)
    show = T[T.cls.isin(["0..20", "pre", "0", "4-6", "11-20", "tail", "G0..15", "betw_early"])]
    log("pooled (class - betw):\n" + show[["cls", "obs", "value", "sigma", "n_series"] +
                                         [c for c in ("delta_lvl_pct", "delta_lvl_sig", "target_pct") if c in T]].round(4).to_string())
    log("front-scrambled pre-0..20:\n" + pd.DataFrame(scr_rows).round(4).to_string())
    log("time-of-day control 0..20-betw:\n" + pd.DataFrame(crow).round(4).to_string())
    open(f"{OUT}/POOL_{SUB.upper()}{tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
