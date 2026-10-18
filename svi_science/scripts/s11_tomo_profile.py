#!/usr/bin/env python
"""S11 tomo (PLAN.md Amendment 7): start-free inlier-profile Wadati per zone.

  python s11_tomo_profile.py --sets AV3,BV3,CV3,S1V3,A,B,C
writes tables/profile_zones.csv, tables/profile_curves.npz, tables/profile_inliers.csv
"""
import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
T = f"{ROOT}/svi_science/results/s11_depth/tomo"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s11_tomo_wadati import along_strike  # noqa: E402

KG = np.round(np.arange(1.30, 2.6001, 0.005), 3)
C = 1 / KG
THR = 0.25
NB, NSC = 1000, 200


def coords():
    sys.path.insert(0, f"{ROOT}/lfe_array_discovery/scripts")
    cwd = os.getcwd(); os.chdir(f"{ROOT}/lfe_array_discovery")
    import a4f_common as AC
    os.chdir(cwd)
    xy = {s: (v[0], v[1]) for s, v in AC.COORD.items()}
    E = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_stations.csv")
    for r in E.itertuples():
        xy.setdefault(r.sta, (r.lat, r.lon))
    return xy


def fam_counts(tp, ts):
    """per-family inlier count over the kappa grid; tp, ts arrays of one family."""
    r = tp[:, None] - C[None, :] * ts[:, None]
    a = np.median(r, axis=0)
    return (np.abs(r - a[None, :]) <= THR).sum(0), r - a[None, :]


def pooled_slope(fid, x, y):
    df = pd.DataFrame(dict(f=fid, x=x, y=y))
    xd = (df.x - df.groupby("f").x.transform("mean")).values
    yd = (df.y - df.groupby("f").y.transform("mean")).values
    return xd, yd, (xd * yd).sum() / (xd ** 2).sum()


def deming(xd, yd, delta):
    sxx, syy, sxy = (xd ** 2).sum(), (yd ** 2).sum(), (xd * yd).sum()
    return (syy - delta * sxx + np.sqrt((syy - delta * sxx) ** 2 + 4 * delta * sxy ** 2)) / (2 * sxy)


def colocated_sigma(P, xy):
    stas = sorted(P.sta.unique())
    out = {"P": [], "S": []}
    ok = P[(P.rP >= 1.2) & (P.rS >= 1.2)].set_index(["fam", "sta"])
    for a, b in itertools.combinations(stas, 2):
        if a not in xy or b not in xy:
            continue
        (la, lo), (lb, lob) = xy[a], xy[b]
        if np.hypot((lo - lob) * 111.19 * np.cos(np.radians(la)), (la - lb) * 111.19) >= 1.0:
            continue
        try:
            j = ok.xs(a, level="sta").join(ok.xs(b, level="sta"), lsuffix="_a", rsuffix="_b", how="inner")
        except KeyError:
            continue
        for ph in "PS":
            d = (j[f"t{ph}_a"] - j[f"t{ph}_b"]).values
            if len(d) >= 5:
                out[ph].append(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2))
    sp = float(np.median(out["P"])) if out["P"] else np.nan
    ss = float(np.median(out["S"])) if out["S"] else np.nan
    if not (np.isfinite(sp) and np.isfinite(ss)):      # no co-located pairs: fall back to the WA-NW (B013 site) values
        sp, ss = 0.039, 0.092
        print(f"  no co-located pairs with >= 5 common families: sigma_P/S fallback {sp}/{ss} s", flush=True)
    return max(sp, 0.025), max(ss, 0.025), len(out["S"])


MINROWS = 4


def zone(Q, rng, sigP, sigS):
    Q = Q.groupby("fam").filter(lambda g: len(g) >= MINROWS)
    fams = Q.fam.unique()
    if len(fams) < 5:
        return dict(n_fam=len(fams), note="< 5 families"), None
    grp = [Q[Q.fam == f] for f in fams]
    cnt = np.array([fam_counts(g.tP.values, g.tS.values)[0] for g in grp])     # (nf, K)
    score = cnt.sum(0)
    k = int(np.argmax(score)); kap = float(KG[k])
    # half-width of the peak at 90 % of the max
    hi = KG[score >= 0.9 * score[k]]
    out = dict(n_fam=len(fams), n_rows=len(Q), kappa=kap, score=int(score[k]), frac_inlier=score[k] / len(Q),
               peak90_lo=float(hi.min()), peak90_hi=float(hi.max()))
    bk = [KG[np.argmax(cnt[rng.integers(0, len(fams), len(fams))].sum(0))] for _ in range(NB)]
    out["sig_boot"] = float(np.std(bk))
    # scrambled control
    sc_peak, sc_k = [], []
    for _ in range(NSC):
        cs = np.array([fam_counts(g.tP.values, g.tS.values[rng.permutation(len(g))])[0] for g in grp]).sum(0)
        sc_peak.append(cs.max()); sc_k.append(KG[np.argmax(cs)])
    out.update(scr_peak_p99=float(np.percentile(sc_peak, 99)), scr_peak_med=float(np.median(sc_peak)),
               scr_kappa_sd=float(np.std(sc_k)), ctrl_pass=bool(score[k] > np.percentile(sc_peak, 99)))
    # LOSO
    lo = {}
    for st in Q.sta.unique():
        cs = np.array([fam_counts(g.tP.values[g.sta.values != st], g.tS.values[g.sta.values != st])[0]
                       if (g.sta.values != st).sum() >= 3 else np.zeros(len(KG)) for g in grp]).sum(0)
        lo[st] = KG[np.argmax(cs)] - kap
    smax = max(lo, key=lambda s: abs(lo[s]))
    out.update(loso_max=float(abs(lo[smax])), loso_sta=smax)
    # inlier rows at the peak
    rows = []
    for g in grp:
        _, r = fam_counts(g.tP.values, g.tS.values)
        rows.append(g[np.abs(r[:, k]) <= THR])
    I = pd.concat(rows)
    I = I.groupby("fam").filter(lambda g: len(g) >= 3)
    xd, yd, b = pooled_slope(I.fam.values, I.tS.values, I.tP.values)
    out["kappa_E1"] = float(1 / b)
    var_ts = float(np.mean(xd ** 2))
    RS = 1 - sigS ** 2 / var_ts
    out.update(R_S=RS, kappa_E3=float(RS / b), kappa_E2=float(1 / deming(xd, yd, sigP ** 2 / sigS ** 2)),
               sigP=sigP, sigS=sigS, n_inlier_rows=len(I), n_inlier_fam=I.fam.nunique(),
               n_inlier_sta=I.sta.nunique())
    # E1 bootstrap over families on inliers
    fl = I.fam.unique(); idx = {f: np.where(I.fam.values == f)[0] for f in fl}
    be = []
    for _ in range(NB):
        pick = rng.choice(fl, len(fl))
        ii = np.concatenate([idx[f] for f in pick]); lab = np.concatenate([[j] * len(idx[f]) for j, f in enumerate(pick)])
        be.append(1 / pooled_slope(lab, I.tS.values[ii], I.tP.values[ii])[2])
    out["sig_E1"] = float(np.std(be))
    # station terms on inliers
    F = pd.get_dummies(I.sta).astype(float)
    Fd = (F - F.groupby(I.fam.values).transform("mean")).values
    A = np.vstack([np.c_[xd, Fd], np.r_[0, np.ones(F.shape[1])]])
    sol = np.linalg.lstsq(A, np.r_[yd, 0], rcond=None)[0]
    out.update(kappa_site=float(1 / sol[0]), d_site=float(1 / sol[0] - 1 / b), site_sd=float(np.std(sol[1:])))
    sig = max(out["sig_boot"], out["sig_E1"])
    out["sig"] = sig
    agree = all(abs(out[c] - kap) <= 2 * max(sig, 0.01) for c in ("kappa_E1", "kappa_E2", "kappa_E3"))
    out["estimators_agree"] = bool(agree)
    out["pass"] = bool(out["ctrl_pass"] and sig <= 0.05 and out["loso_max"] <= 0.05 and abs(out["d_site"]) <= 0.02
                       and agree)
    return out, (score, np.array(sc_peak), I)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", default="AV3,BV3,CV3,S1V3,A,B,C")
    ap.add_argument("--rmin", type=float, default=1.2)
    ap.add_argument("--out", default="profile")
    ap.add_argument("--minrows", type=int, default=4)
    a = ap.parse_args()
    global MINROWS
    MINROWS = a.minrows
    rng = np.random.default_rng(31)
    xy = coords()
    rows, curves, inl = [], {}, []
    for st in a.sets.split(","):
        P = pd.read_csv(f"{T}/tables/picks_{st}.csv", low_memory=False)
        P["fam"] = P.fam.astype(str)
        sigP, sigS, npair = colocated_sigma(P, xy)
        P = P[P.tP.notna() & P.tS.notna() & (P.rP >= a.rmin) & (P.rS >= a.rmin)].copy()
        if "eye_lat" in P:
            P["lat"] = P.eye_lat.fillna(P.hyp_lat); P["lon"] = P.eye_lon.fillna(P.hyp_lon)
            P["depth_km"] = P.eye_dep.fillna(P.hyp_dep)
        P["along"] = along_strike(P.lat.values, P.lon.values)
        if "eye_class" in P:
            zones = [("eye_lfe", P[P.eye_class == "lfe"]), ("eye_noise", P[P.eye_class == "noise"])]
        else:
            zones = [("all", P)]
            if P.lat.notna().any():
                lo = np.floor(P.along.min() / 30) * 30
                for z0 in np.arange(lo, P.along.max() + 1, 30):
                    zones.append((f"along_{z0:+.0f}_{z0 + 30:+.0f}", P[(P.along >= z0) & (P.along < z0 + 30)]))
                for d0, d1 in ((0, 30), (30, 36), (36, 42), (42, 200)):
                    zones.append((f"depth_{d0}_{d1}", P[(P.depth_km >= d0) & (P.depth_km < d1)]))
        for lab, Q in zones:
            o, extra = zone(Q, rng, sigP, sigS)
            o.update(set=st, zone=lab, n_coloc_pairs=npair, along_med=float(Q.along.median()) if len(Q) else np.nan,
                     depth_med=float(Q.depth_km.median()) if len(Q) else np.nan)
            rows.append(o)
            if extra is not None:
                curves[f"{st}|{lab}"] = extra[0]; curves[f"{st}|{lab}|scr"] = extra[1]
                inl.append(extra[2].assign(set=st, zone=lab))
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in o.items()
                   if k in ("set", "zone", "n_fam", "n_rows", "kappa", "sig", "frac_inlier", "score", "scr_peak_p99",
                            "kappa_E1", "kappa_E2", "kappa_E3", "loso_max", "d_site", "pass", "note")}, flush=True)
    pd.DataFrame(rows).to_csv(f"{T}/tables/{a.out}_zones.csv", index=False)
    np.savez_compressed(f"{T}/tables/{a.out}_curves.npz", kgrid=KG, **curves)
    if inl:
        pd.concat(inl).to_csv(f"{T}/tables/{a.out}_inliers.csv", index=False)


if __name__ == "__main__":
    main()
