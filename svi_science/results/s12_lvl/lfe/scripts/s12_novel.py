#!/usr/bin/env python
"""s12_novel -- analysis of coordinator tasks A-D (PLAN.md §6) from s12_novel_measure.py output.
    S11_REGION=SVI|WA_NW python s12_novel.py [--test]
Outputs tables/{sub}/novel_{A,A_curve,B,C,D,D_majors}.csv, novel_verdict.json; figures/fig_novel_{sub}.png.
"""
from __future__ import annotations
import os, sys, json, pickle, time, argparse
import numpy as np
import pandas as pd
from scipy.stats import norm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s12_epoch import C, OUT, SUB, TB, log

WN = ["2-5", "4-6"]
DELTA = 0.02
SIG_AS = 7.0
rng = np.random.default_rng(11)


def fshape(d, v=10.0):
    return norm.cdf(v * np.asarray(d, float) / SIG_AS)


def a_tables():
    p = f"{OUT}/tables/fwd_coef.parquet"
    K = pd.read_parquet(p); K = K[(K.variant == "ell40") & (K.comp == "3c")]
    return {(r.series, r.win): r.a_lvl_nom for r in K.itertuples()}, K.groupby("win").a_lvl_nom.median().to_dict()


def site_of(key):
    ns = key.split("|")[1]
    return C.SITE.get(ns, ns)


# ------------------------------------------------------------------------------ A
def task_A(D, aT, aM):
    keys = sorted(k for k in D if D[k]["res"].get("betw") is not None)
    days = list(range(-10, 21))
    rows, curve = [], []
    out = {}
    for wn in WN:
        for v in (10.0, 5.0):
            f = fshape(days, v)
            X = np.full((len(keys), len(days), 201), np.nan); a_ = np.zeros(len(keys))
            for i, k in enumerate(keys):
                r = D[k]["res"]; b = np.asarray(r["betw"][wn], float)
                for j, d in enumerate(days):
                    if f"D{d}" in r:
                        X[i, j] = np.asarray(r[f"D{d}"][wn], float) - b
                a_[i] = aT.get((k, wn), aM.get(wn, np.nan))
            S = np.nanstd(X[:, :, 1:], 2)
            fl = np.nanpercentile(S[S > 0], 25)
            W = np.where(np.isfinite(X[:, :, 0]) & (S > 0), 1 / np.maximum(S, fl) ** 2, 0.0)

            def fitA(Y, W):
                ok = np.isfinite(Y) & (W > 0)
                lo_ = (ok & (f[None, :] < 0.3)).any(1); hi_ = (ok & (f[None, :] > 0.7)).sum(1) >= 1
                ok &= (lo_ & hi_)[:, None]
                Wk = np.where(ok, W, 0); Yk = np.where(ok, Y, 0)
                sw = Wk.sum(1, keepdims=True); sw[sw == 0] = np.nan
                fb = (Wk * f[None, :]).sum(1, keepdims=True) / sw; yb = (Wk * Yk).sum(1, keepdims=True) / sw
                df = np.where(ok, f[None, :] - fb, 0); dy = np.where(ok, Yk - yb, 0)
                A = np.nansum(Wk * df * dy) / np.nansum(Wk * df ** 2)
                aa = np.where(np.isfinite(a_), a_, 0)[:, None]
                dl = np.nansum(Wk * aa * df * dy) / np.nansum(Wk * (aa * df) ** 2)
                return A, dl
            A0, d0 = fitA(X[:, :, 0], W)
            Ab = np.array([fitA(X[:, :, b], W) for b in range(1, 201)])
            # front-scrambled
            As = []
            for kdraw in range(200):
                Ys = np.full((len(keys), len(days)), np.nan)
                for i, k in enumerate(keys):
                    sc = D[k]["scr"]
                    if sc is None:
                        continue
                    lab = sc["lab"]; m = lab[:, 0] == kdraw
                    if not m.any():
                        continue
                    b0 = float(D[k]["res"]["betw"][wn][0])
                    for d, val in zip(lab[m, 1].astype(int), np.asarray(sc[wn])[m]):
                        Ys[i, days.index(d)] = val - b0
                As.append(fitA(Ys, W)[0])
            As = np.array(As)
            tgt = 100 * DELTA * aM.get(wn, np.nan)
            r = dict(win=wn, v_kmd=v, A=A0, A_sig=float(np.nanstd(Ab[:, 0])), delta_lvl_pct=100 * d0 / 100,
                     delta_sig=float(np.nanstd(Ab[:, 1])), A_target=tgt, scr_sd=float(np.nanstd(As)),
                     scr_p=float((np.sum(np.abs(As) >= abs(A0)) + 1) / (np.isfinite(As).sum() + 1)),
                     n_series=int((W.sum(1) > 0).sum()))
            r["bound2_over_target"] = (abs(A0) + 2 * r["A_sig"]) / tgt
            rows.append(r)
            if v == 10.0:
                # pooled daily curve (inverse-variance over series, no offsets) + 20 scrambled curves
                Pc = np.nansum(W[:, :, None] * np.nan_to_num(X), 0) / np.maximum((W[:, :, None] * np.isfinite(X)).sum(0), 1e-30)
                for j, d in enumerate(days):
                    curve.append(dict(win=wn, day=d, value=Pc[j, 0], weight=W[:, j].sum(),
                                      boot_sd=float(np.nanstd(Pc[j, 1:])), target_curve=tgt * f[j]))
    return pd.DataFrame(rows), pd.DataFrame(curve)


# ------------------------------------------------------------------------------ day table helpers
def day_frame(D):
    fr = []
    for k, v in D.items():
        d = pd.DataFrame(v["day"]); d["key"] = k; d["site"] = site_of(k); d["fam"] = int(k.split("|")[0])
        fr.append(d)
    F = pd.concat(fr, ignore_index=True)
    F["date"] = pd.to_datetime(F.date)
    doy = F.date.dt.dayofyear.to_numpy()
    F["cs"] = np.cos(2 * np.pi * doy / 365.25); F["sn"] = np.sin(2 * np.pi * doy / 365.25)
    return F


def within(F, cols, w):
    """demean columns within series (weighted)"""
    G = F[cols].copy()
    for c in cols:
        num = (F[c] * w).groupby(F.key).transform("sum"); den = w.groupby(F.key).transform("sum")
        G[c] = F[c] - num / den.where(den > 0)
    return G.fillna(0.0)


def wls(y, X, w):
    Xw = X * w[:, None]
    A = X.T @ Xw
    try:
        return np.linalg.lstsq(A, Xw.T @ y, rcond=None)[0]
    except np.linalg.LinAlgError:
        return np.full(X.shape[1], np.nan)


# ------------------------------------------------------------------------------ B
def task_B(F):
    rows = []
    Fb = F[F.betw & (F.prev < 2000)].copy()
    Fb["ty"] = Fb.prev / 365.25; Fb["ex"] = np.exp(-Fb.prev / 45.0)
    Fb["clu"] = (Fb.major.astype(str) + "_" + (Fb.prev // 1).astype(int).astype(str)).str.len()  # placeholder
    # cluster = inter-ETS interval: date of the previous front ~ date - prev
    Fb["clu"] = ((Fb.date - pd.to_timedelta(Fb.prev, "D")).dt.year * 100 + (Fb.date - pd.to_timedelta(Fb.prev, "D")).dt.month)
    for q in ["dvv_2-5", "dvv_4-6", "slope"]:
        for form in ("linear", "exp45"):
            for ann in (True, False):
                d = Fb[np.isfinite(Fb[q])].copy()
                if q.startswith("dvv"):
                    cc = d[f"dvv_cc{q.split('_')[1]}"] if f"dvv_cc{q.split('_')[1]}" in d else None
                    d = d[(cc > 0.3)] if cc is not None else d
                    w = (d[f"dvv_cc{q.split('_')[1]}"] ** 2).to_numpy()
                else:
                    w = np.ones(len(d))
                if len(d) < 30:
                    continue
                reg = "ty" if form == "linear" else "ex"
                cols = [reg]
                if ann:
                    for s in sorted(d.site.unique()):
                        d[f"c_{s}"] = d.cs * (d.site == s); d[f"s_{s}"] = d.sn * (d.site == s); cols += [f"c_{s}", f"s_{s}"]
                ww = pd.Series(w, index=d.index)
                G = within(d, [q] + cols, ww)
                y = G[q].to_numpy(); X = G[cols].to_numpy()
                b0 = wls(y, X, w)[0]
                cl = d.clu.to_numpy(); uc = np.unique(cl)
                bs = []
                for k in range(200):
                    m = rng.multinomial(len(uc), np.ones(len(uc)) / len(uc))
                    mult = dict(zip(uc, m)); wm = w * np.array([mult[c] for c in cl])
                    if wm.sum() == 0:
                        continue
                    Gm = within(d, [q] + cols, pd.Series(wm, index=d.index))
                    try:
                        bs.append(wls(Gm[q].to_numpy(), Gm[cols].to_numpy(), wm)[0])
                    except Exception:
                        pass
                r = dict(obs=q, form=form, annual=ann, b=b0, b_sig=float(np.nanstd(bs)), n_days=len(d), n_series=d.key.nunique(),
                         n_intervals=len(uc))
                if ann and form == "linear":
                    for s in sorted(d.site.unique()):
                        m = (d.site == s).to_numpy()
                        if m.sum() > 30:
                            Gs = within(d[m], [q] + cols, ww[m]); cs_ = [reg] + [c for c in cols if c.endswith(f"_{s}")]
                            r[f"site_{s}"] = float(wls(Gs[q].to_numpy(), Gs[cs_].to_numpy(), w[m])[0])
                    for par, m in (("odd", (d.fam % 2 == 1).to_numpy()), ("even", (d.fam % 2 == 0).to_numpy())):
                        if m.sum() > 30:
                            Gs = within(d[m], [q] + cols, ww[m])
                            r[par] = float(wls(Gs[q].to_numpy(), Gs[cols].to_numpy(), w[m])[0])
                rows.append(r)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ C
def task_C(F):
    F = F.copy()
    cls = np.full(len(F), "", dtype=object)
    lag = F.lag.to_numpy()
    cls[(lag >= -10) & (lag <= -1)] = "pre"; cls[lag == 0] = "0"; cls[lag == 1] = "1"
    cls[(lag >= 2) & (lag <= 3)] = "2-3"; cls[(lag >= 4) & (lag <= 20)] = "4-20"; cls[F.betw.to_numpy()] = "betw"
    F["cls"] = cls; F = F[F.cls != ""].copy()
    F["lnws"] = np.log(np.maximum(F.ws, 1e-3)); F["lnnd"] = np.log(np.where(F.nd > 0, F.nd, np.nan))
    F["clu"] = np.where(F.cls == "betw", 1000 + F.major, F.major)
    rows = []
    levels = ["pre", "0", "1", "2-3", "4-20"]
    for q in ["sr", "width", "lnPS_Z_rms", "slope"]:
        for adj in (True, False):
            d = F[np.isfinite(F[q])].copy()
            if adj and np.isfinite(d.lnnd).mean() > 0.8:
                d = d[np.isfinite(d.lnnd)]
            for L in levels:
                d[f"k_{L}"] = (d.cls == L).astype(float)
            cols = [f"k_{L}" for L in levels] + (["lnws"] + (["lnnd"] if np.isfinite(d.lnnd).all() else []) if adj else [])
            w = np.ones(len(d)); ww = pd.Series(w, index=d.index)
            G = within(d, [q] + cols, ww)
            b0 = wls(G[q].to_numpy(), G[cols].to_numpy(), w)
            cl = d.clu.to_numpy(); uc = np.unique(cl); bs = []
            for k in range(200):
                m = rng.multinomial(len(uc), np.ones(len(uc)) / len(uc)); mult = dict(zip(uc, m))
                wm = np.array([mult[c] for c in cl], float)
                Gm = within(d, [q] + cols, pd.Series(wm, index=d.index))
                bs.append(wls(Gm[q].to_numpy(), Gm[cols].to_numpy(), wm))
            bs = [b_ for b_ in bs if np.all(np.isfinite(b_))] if k == 199 else bs
            bs = np.array(bs)
            for j, L in enumerate(levels):
                r = dict(obs=q, adjusted=adj, cls=L, minus="betw", value=float(b0[j]), sigma=float(np.nanstd(bs[:, j])),
                         n_days=int((d.cls == L).sum()))
                if adj:
                    for s in sorted(d.site.unique()):
                        m = (d.site == s).to_numpy()
                        Gs = within(d[m], [q] + cols, ww[m])
                        try:
                            r[f"site_{s}"] = float(wls(Gs[q].to_numpy(), Gs[cols].to_numpy(), w[m])[j])
                        except Exception:
                            pass
                rows.append(r)
            # behind - at front: (2-3, 4-20) vs (0)
            r = dict(obs=q, adjusted=adj, cls="2-20 minus 0", minus="0",
                     value=float(0.5 * (b0[3] + b0[4]) - b0[1]), sigma=float(np.nanstd(0.5 * (bs[:, 3] + bs[:, 4]) - bs[:, 1])))
            rows.append(r)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ D
def task_D(D, majors, aM):
    sizes = pd.read_csv(f"{C.ROOT}/svi_science/results/s11_depth/round6/{SUB}/tables/r6_sizes.csv").set_index("episode")
    rows = []
    for wn in WN:
        per = {}
        for k, v in D.items():
            b = v["res"].get("betw")
            if b is None:
                continue
            bv = float(b[wn][0]); bs = float(np.nanstd(np.asarray(b[wn][1:], float)))
            fam = int(k.split("|")[0])
            for m, o in v["per"].items():
                x = float(o[wn][0]) - bv; s = float(np.sqrt(np.nanstd(np.asarray(o[wn][1:], float)) ** 2 + bs ** 2))
                if np.isfinite(x) and s > 0:
                    per.setdefault(m, []).append((fam, x, s, o["n_days"], o.get("nd", np.nan)))
        for m, lst in per.items():
            a = pd.DataFrame(lst, columns=["fam", "x", "s", "n_days", "nd"])
            fl = a.s.quantile(.25); a["w"] = 1 / np.maximum(a.s, fl) ** 2
            fm = a.groupby("fam").apply(lambda g: pd.Series(dict(x=np.sum(g.w * g.x) / g.w.sum(), w=g.w.sum(), nd=np.nanmean(g.nd))),
                                        include_groups=False)
            val = float(np.sum(fm.w * fm.x) / fm.w.sum())
            sig = float(np.sqrt(1 / fm.w.sum()) * max(1.0, np.sqrt(np.sum(fm.w * (fm.x - val) ** 2) / max(len(fm) - 1, 1))))
            rows.append(dict(win=wn, major=m, value=val, sigma=sig, n_fam=len(fm), n_series=len(a), lognd=float(np.log(np.nanmean(fm.nd)))))
    Mj = pd.DataFrame(rows)
    Mj = Mj.join(sizes[[c for c in ("Du", "Fs", "T", "L", "M0") if c in sizes]], on="major")
    out = []
    for wn in WN:
        g = Mj[(Mj.win == wn)].copy()
        tgt = 100 * DELTA * aM.get(wn, np.nan)
        for sz in [c for c in ("Du", "Fs", "T", "L", "M0") if c in g]:
            h = g[np.isfinite(g[sz]) & np.isfinite(g.value)]
            if len(h) < 5:
                continue
            for form in ("size_only", "size_lognd"):
                X = np.column_stack([np.ones(len(h)), h[sz]] + ([h.lognd] if form == "size_lognd" else []))
                w = 1 / h.sigma.to_numpy() ** 2
                b = wls(h.value.to_numpy(), X, w)
                bb = []
                for k in range(1000):
                    i = rng.integers(0, len(h), len(h))
                    try:
                        bb.append(wls(h.value.to_numpy()[i], X[i], w[i])[1])
                    except Exception:
                        pass
                out.append(dict(win=wn, size=sz, form=form, slope=float(b[1]), slope_sig=float(np.nanstd(bb)), n_majors=len(h),
                                mean_size=float(h[sz].mean()), target_slope=tgt / float(h[sz].mean()),
                                mean_contrast=float(np.sum(w * h.value) / w.sum())))
    return Mj, pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--test", action="store_true"); a = ap.parse_args()
    tag = "_test" if a.test else ""
    R = pickle.load(open(f"{TB}/novel_raw{tag}.pkl", "rb")); D = R["out"]
    aT, aM = a_tables()
    if C.REGION != "SVI":
        aT = {}
    t0 = time.time()
    A, Ac = task_A(D, aT, aM); A.to_csv(f"{TB}/novel_A{tag}.csv", index=False); Ac.to_csv(f"{TB}/novel_A_curve{tag}.csv", index=False)
    log(f"[A] {time.time() - t0:.0f}s\n" + A.round(4).to_string())
    F = day_frame(D)
    F.to_parquet(f"{TB}/novel_days{tag}.parquet", index=False)
    B = task_B(F); B.to_csv(f"{TB}/novel_B{tag}.csv", index=False)
    log(f"[B] {time.time() - t0:.0f}s\n" + B.round(4).to_string())
    Cc = task_C(F); Cc.to_csv(f"{TB}/novel_C{tag}.csv", index=False)
    log(f"[C] {time.time() - t0:.0f}s\n" + Cc[Cc.adjusted].round(4).to_string())
    Mj, Dd = task_D(D, R["majors"], aM); Mj.to_csv(f"{TB}/novel_D_majors{tag}.csv", index=False); Dd.to_csv(f"{TB}/novel_D{tag}.csv", index=False)
    log(f"[D] {time.time() - t0:.0f}s\n" + Dd.round(4).to_string())
    open(f"{OUT}/NOVELAN_{SUB.upper()}{tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
