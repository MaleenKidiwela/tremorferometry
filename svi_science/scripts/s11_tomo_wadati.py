#!/usr/bin/env python
"""S11 tomo pilot: per-family Wadati κ + gates (tomo/PLAN.md §1-§3).

  python s11_tomo_wadati.py --set B [--rmin 1.2] [--tag ""]
Reads tomo/tables/picks_<set>.csv (columns fam, sta, tP, rP, tS, rS, + metadata), writes
  tables/kappa_<set><tag>.csv  (one row per family, ALL families incl. unmeasured)
  tables/gates_<set><tag>.json (gate statistics)
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy import stats

T = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/tomo"
NB, NS = 1000, 1000
SIG_OK, RMS_OK = 0.03, 0.25


def ols(x, y):
    """vectorised OLS along last axis: returns slope, intercept, rms."""
    xm, ym = x.mean(-1, keepdims=True), y.mean(-1, keepdims=True)
    sxx = ((x - xm) ** 2).sum(-1)
    b = ((x - xm) * (y - ym)).sum(-1) / np.where(sxx > 0, sxx, np.nan)
    a = ym[..., 0] - b * xm[..., 0]
    r = y - (a[..., None] + b[..., None] * x)
    return b, a, np.sqrt((r ** 2).mean(-1))


def tls_slope(x, y):
    xm, ym = x - x.mean(), y - y.mean()
    _, _, vt = np.linalg.svd(np.c_[xm, ym], full_matrices=False)
    v = vt[-1]
    return -v[0] / v[1]


def robust_mask(x, y, nmin):
    """V1: iterative clipping |res| > max(3*1.4826*MAD, 0.15 s); returns keep mask (None if < nmin kept)."""
    m = np.ones(len(x), bool)
    for _ in range(20):
        if m.sum() < max(nmin, 3):
            return None
        b, a, _ = ols(x[m], y[m])
        r = y - (a + b * x)
        thr = max(3 * 1.4826 * np.median(np.abs(r[m] - np.median(r[m]))), 0.15)
        m2 = np.abs(r - np.median(r[m])) <= thr
        if (m2 == m).all():
            break
        m = m2
    return m if m.sum() >= nmin else None


EST = "ols"


def XY(tp, ts):
    """ols: (tS - tP) on tP, kappa = 1 + b;  inv (E1): tP on tS, kappa = 1/b."""
    return (tp, ts - tp) if EST == "ols" else (ts, tp)


def KAP(b):
    b = np.asarray(b, float)
    return 1 + b if EST == "ols" else np.where(b != 0, 1 / np.where(b != 0, b, 1), np.nan)


def family(q, rng, robust=False, nmin=6):
    qfull = q
    tp, s = q.tP.values, q.tS.values
    x, y = XY(tp, s)
    if robust:
        k = robust_mask(x, y, nmin)
        if k is None:
            return None
        q, tp, s, x, y = q[k], tp[k], s[k], x[k], y[k]
        nclip = int((~k).sum())
    n = len(x)
    b, a, rms = ols(x, y)
    out = dict(n_sta=n, kappa=float(KAP(b)), b=float(b), icpt=a, rms=rms, tP_range=float(tp.max() - tp.min()),
               n_clipped=nclip if robust else 0, kept=";".join(q.sta), est=EST)
    # bootstrap over stations
    idx = rng.integers(0, n, (NB, n))
    ok = np.array([len(np.unique(i)) >= 3 for i in idx])
    bb, _, _ = ols(x[idx[ok]], y[idx[ok]])
    kb = KAP(bb[np.isfinite(bb)]); kb = kb[np.isfinite(kb)]
    out.update(sig_k=float(np.std(kb)), hw68=float(np.diff(np.percentile(kb, [16, 84]))[0] / 2))
    # leave one station out
    lo = np.array([float(KAP(ols(x[np.arange(n) != i], y[np.arange(n) != i])[0])) for i in range(n)])
    out.update(loso_max=float(np.nanmax(np.abs(lo - out["kappa"]))),
               loso_sta=q.sta.values[int(np.nanargmax(np.abs(lo - out["kappa"])))])
    # scrambled control: permute S picks across stations (break P-S pairing); robust form re-clips from the full set
    tp0, s0 = (qfull.tP.values, qfull.tS.values) if robust else (tp, s)
    n0 = len(tp0)
    perm0 = np.argsort(rng.random((NS, n0)), axis=1)
    sb, srms = np.full(NS, np.nan), np.full(NS, np.inf)
    kept_first = None
    for i in range(NS):
        xx, yy = XY(tp0, s0[perm0[i]])
        kk = robust_mask(xx, yy, nmin) if robust else np.ones(n0, bool)
        if kk is not None:
            sb[i], _, srms[i] = ols(xx[kk], yy[kk])
            if kept_first is None:
                kept_first = (xx[kk], yy[kk], i)
    out.update(scr_kappa=KAP(sb), scr_rms=srms, p_rms=float((srms <= rms).mean()),
               scr_kappa_med=float(np.nanmedian(KAP(sb))), scr_rms_med=float(np.median(srms)),
               scr_b_med=float(np.nanmedian(sb)), scr_b_sd=float(np.nanstd(sb)),
               scr_frac_measurable=float(np.isfinite(srms).mean()))
    out["z_b_vs_scr"] = (b - out["scr_b_med"]) / out["scr_b_sd"] if out["scr_b_sd"] > 0 else np.nan
    if kept_first is None:
        out.update(scr1_kappa=np.nan, scr1_rms=np.inf, scr1_sig=np.inf)
    else:
        xk, yk, i = kept_first
        ib = rng.integers(0, len(xk), (NB, len(xk)))
        okb = np.array([len(np.unique(j)) >= 3 for j in ib])
        bb1, _, _ = ols(xk[ib[okb]], yk[ib[okb]])
        k1 = KAP(bb1); k1 = k1[np.isfinite(k1)]
        out.update(scr1_kappa=float(KAP(sb[i])), scr1_rms=float(srms[i]), scr1_sig=float(np.std(k1)) if len(k1) else np.inf)
    # columns
    out["kappa_tS_on_tP"] = float(ols(tp, s)[0])
    out["kappa_tls"] = float(tls_slope(tp, s))
    out["kappa_inv_col"] = float(1 / ols(s, tp)[0])
    return out


def spearman_perm(a, b, rng, n=10000):
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if len(a) < 5:
        return np.nan, np.nan, len(a)
    r = stats.spearmanr(a, b)[0]
    rb = np.array([stats.spearmanr(a, rng.permutation(b))[0] for _ in range(n)])
    return float(r), float((np.abs(rb) >= abs(r)).mean()), int(len(a))


def along_strike(lat, lon, az=330.0, lat0=47.0, lon0=-123.0):
    x = (lon - lon0) * 111.19 * np.cos(np.radians(lat0))
    y = (lat - lat0) * 111.19
    return x * np.sin(np.radians(az)) + y * np.cos(np.radians(az))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True)
    ap.add_argument("--rmin", type=float, default=1.2)
    ap.add_argument("--nmin", type=int, default=6)
    ap.add_argument("--tag", default="")
    ap.add_argument("--robust", action="store_true", help="V1 iterative clipping")
    ap.add_argument("--maxdist", type=float, default=0, help="V2: epicentral distance limit (km)")
    ap.add_argument("--est", default="ols", choices=["ols", "inv"])
    a = ap.parse_args()
    global EST
    EST = a.est
    rng = np.random.default_rng(12)
    P = pd.read_csv(f"{T}/tables/picks_{a.set}.csv", low_memory=False)
    P["fam"] = P.fam.astype(str)
    if a.maxdist > 0:
        import sys, os
        sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_array_discovery/scripts")
        cwd = os.getcwd(); os.chdir("/home/jovyan/tremorferometry/lfe_array_discovery")
        import a4f_common as AC
        os.chdir(cwd)
        if "eye_lat" in P:
            flat, flon = P.eye_lat.fillna(P.hyp_lat), P.eye_lon.fillna(P.hyp_lon)
        else:
            flat, flon = P.lat, P.lon
        slat = P.sta.map(lambda s_: AC.COORD[s_][0]); slon = P.sta.map(lambda s_: AC.COORD[s_][1])
        P["epi_km"] = np.hypot((slon - flon) * 111.19 * np.cos(np.radians(flat)), (slat - flat) * 111.19)
        P = P[P.epi_km <= a.maxdist]
    okrow = P.tP.notna() & P.tS.notna() & (P.rP >= a.rmin) & (P.rS >= a.rmin)
    meta_cols = [c for c in ("eye_class", "hyp_lat", "hyp_lon", "hyp_dep", "eye_lat", "eye_lon", "eye_dep",
                             "lat", "lon", "depth_km", "n_detections", "era_start", "era_end",
                             "location_quality", "grp", "src") if c in P]
    M = P.groupby("fam")[meta_cols].first()
    rows, scr = [], {}
    for f, q in P[okrow].groupby("fam"):
        if len(q) < a.nmin:
            continue
        o = family(q, rng, a.robust, a.nmin)
        if o is None:
            continue
        scr[f] = (o.pop("scr_kappa"), o.pop("scr_rms"))
        o["fam"] = f
        o["stations"] = o.pop("kept")
        rows.append(o)
    K = pd.DataFrame(rows).set_index("fam")
    K = M.join(K, how="left")
    K["n_ok_rows"] = P[okrow].groupby("fam").size().reindex(K.index).fillna(0).astype(int)
    K["measured"] = K.kappa.notna()
    K["good"] = K.measured & (K.sig_k <= SIG_OK) & (K.rms <= RMS_OK)
    # location for mapping
    if "eye_lat" in K:                       # set A: eye-checked relocation where present, else run-1 Hyp2000
        K["map_lat"] = K.eye_lat.fillna(K.hyp_lat); K["map_lon"] = K.eye_lon.fillna(K.hyp_lon)
        K["map_dep"] = K.eye_dep.fillna(K.hyp_dep)
    else:
        K["map_lat"], K["map_lon"], K["map_dep"] = K.lat, K.lon, K.depth_km
    K["along_km"] = along_strike(K.map_lat.values, K.map_lon.values)
    if "era_start" in K and K.era_start.notna().any():
        days = (pd.to_datetime(K.era_end) - pd.to_datetime(K.era_start)).dt.days.clip(lower=1)
        K["rate_per_day"] = K.n_detections / days
    G = {"set": a.set, "rmin": a.rmin, "nmin": a.nmin, "n_families": int(len(K)),
         "n_measured": int(K.measured.sum()), "n_good": int(K.good.sum()),
         "median_rms_measured": float(K.rms.median()), "median_sig_measured": float(K.sig_k.median())}
    if "eye_class" in K:
        for c in ("lfe", "noise"):
            k = K[K.eye_class == c]
            G[f"{c}_measured"], G[f"{c}_good"] = int(k.measured.sum()), int(k.good.sum())
    G["G1_pass"] = G["n_good"] >= 50
    good = K[K.good]
    gf = list(good.index)
    G["kappa_good_median"] = float(good.kappa.median()) if len(good) else np.nan
    G["kappa_good_p16_p84"] = [float(v) for v in np.percentile(good.kappa, [16, 84])] if len(good) else []
    # control good count (one scramble)
    meas = K[K.measured]
    G["scr_good_count"] = int(((meas.scr1_sig <= SIG_OK) & (meas.scr1_rms <= RMS_OK)).sum())
    if len(gf) >= 3:
        sp_real = 1.4826 * stats.median_abs_deviation(good.kappa)
        SK = np.array([scr[f][0] for f in gf])        # (nfam, NS)
        sp_ctrl = 1.4826 * stats.median_abs_deviation(SK, axis=0, nan_policy="omit")
        G["G2a_spread_real"] = float(sp_real)
        G["G2a_spread_ctrl_p05_p50_p95"] = [float(v) for v in np.percentile(sp_ctrl, [5, 50, 95])]
        G["G2a_literal_pass"] = bool(sp_real > np.percentile(sp_ctrl, 95))
        G["G2a_real_below_ctrl_p05"] = bool(sp_real < np.percentile(sp_ctrl, 5))
        G["G2b_frac_good_p_lt_05"] = float((good.p_rms < 0.05).mean())
        G["G2b_pass"] = bool(G["G2b_frac_good_p_lt_05"] >= 0.8 and G["scr_good_count"] <= 0.05 * len(gf))
        w = 1 / good.sig_k.clip(lower=1e-3) ** 2
        km = float((w * good.kappa).sum() / w.sum())
        chi2 = float((((good.kappa - km) / good.sig_k.clip(lower=1e-3)) ** 2).sum())
        nu = len(gf) - 1
        G.update(G2c_wmean=km, G2c_chi2_nu=chi2 / nu, G2c_p=float(stats.chi2.sf(chi2, nu)),
                 G2c_pass=bool(stats.chi2.sf(chi2, nu) < 0.05))
        G["G3_frac_loso_le_005"] = float((good.loso_max <= 0.05).mean())
        G["G3_pass"] = bool(G["G3_frac_loso_le_005"] >= 0.8)
    # split-half (families with >= 10 stations), real vs scrambled
    big = [f for f in K.index[K.n_ok_rows >= 10]]
    if len(big) >= 5:
        Pk = P[okrow & P.fam.isin(big)]
        cr, cc = [], []
        for _ in range(200):
            h1, h2, s1, s2 = [], [], [], []
            for f, q in Pk.groupby("fam"):
                x, s = q.tP.values, q.tS.values
                i = rng.permutation(len(x)); A_, B_ = i[: len(x) // 2], i[len(x) // 2:]
                h1.append(ols(x[A_], s[A_] - x[A_])[0]); h2.append(ols(x[B_], s[B_] - x[B_])[0])
                sp = s[rng.permutation(len(s))]
                s1.append(ols(x[A_], sp[A_] - x[A_])[0]); s2.append(ols(x[B_], sp[B_] - x[B_])[0])
            cr.append(stats.spearmanr(h1, h2)[0]); cc.append(stats.spearmanr(s1, s2)[0])
        G.update(split_n=len(big), split_rho_real_med=float(np.nanmedian(cr)),
                 split_rho_ctrl_med=float(np.nanmedian(cc)),
                 split_rho_ctrl_p95=float(np.nanpercentile(cc, 95)))
    # G4 site terms: joint fit over measured families
    if len(meas) >= 5:
        mf = list(meas.index)
        Q = P[okrow & P.fam.isin(mf)].copy()
        Q = Q[[st in K.loc[f_, "stations"].split(";") for f_, st in zip(Q.fam, Q.sta)]]
        stas = sorted(Q.sta.unique()); si = {s: i for i, s in enumerate(stas)}
        fi = {f: i for i, f in enumerate(mf)}
        nf, ns = len(mf), len(stas)
        Amat = np.zeros((len(Q) + 1, 2 * nf + ns)); yv = np.zeros(len(Q) + 1)
        for r, row in enumerate(Q.itertuples()):
            j = fi[row.fam]
            xx_, yy_ = XY(row.tP, row.tS)
            Amat[r, 2 * j], Amat[r, 2 * j + 1], Amat[r, 2 * nf + si[row.sta]] = 1, xx_, 1
            yv[r] = yy_
        Amat[-1, 2 * nf:] = 1e3
        sol = np.linalg.lstsq(Amat, yv, rcond=None)[0]
        K.loc[mf, "kappa_joint"] = KAP(sol[1:2 * nf:2])
        site = pd.Series(sol[2 * nf:], index=stas, name="site_term_s")
        site.to_csv(f"{T}/tables/site_terms_{a.set}{a.tag}.csv")
        if len(gf) >= 5:
            d = (K.loc[gf, "kappa_joint"] - K.loc[gf, "kappa"]).abs()
            rho = stats.spearmanr(K.loc[gf, "kappa_joint"], K.loc[gf, "kappa"])[0]
            G.update(G4i_med_abs_dk=float(d.median()), G4i_spearman=float(rho),
                     G4i_pass=bool(d.median() <= 0.02 and rho >= 0.8),
                     site_term_sd=float(site.std()))
            # presence regression
            pres = pd.DataFrame({s: [s in good.loc[f, "stations"].split(";") for f in gf] for s in stas},
                                index=gf).astype(float)
            pres = pres.loc[:, (pres.sum() >= 5) & (pres.sum() <= len(gf) - 1)]
            if pres.shape[1] >= 1 and len(gf) > pres.shape[1] + 2:
                X = np.c_[np.ones(len(gf)), pres.values]
                yk = good.kappa.values

                def r2(yy):
                    c = np.linalg.lstsq(X, yy, rcond=None)[0]
                    return 1 - ((yy - X @ c) ** 2).sum() / ((yy - yy.mean()) ** 2).sum()
                R2 = r2(yk); R2p = np.array([r2(rng.permutation(yk)) for _ in range(1000)])
                G.update(G4ii_n_stations=int(pres.shape[1]), G4ii_R2=float(R2),
                         G4ii_R2_perm_p95=float(np.percentile(R2p, 95)),
                         G4ii_p=float((R2p >= R2).mean()), G4ii_pass=bool(R2 < np.percentile(R2p, 95)))
            else:
                G["G4ii_note"] = f"too few families ({len(gf)}) for {pres.shape[1]} station indicators"
    # mapping correlations
    for sub, KK in (("good", K[K.good]), ("measured", K[K.measured])):
        for col in ("map_dep", "along_km", "rate_per_day"):
            if col in KK:
                r, p, n = spearman_perm(KK.kappa.values, KK[col].values.astype(float), rng, 2000)
                G[f"map_{sub}_{col}"] = dict(rho=r, p=p, n=n)
    K.to_csv(f"{T}/tables/kappa_{a.set}{a.tag}.csv")
    json.dump(G, open(f"{T}/tables/gates_{a.set}{a.tag}.json", "w"), indent=1, default=float)
    print(json.dumps(G, indent=1, default=float))


if __name__ == "__main__":
    main()
