#!/usr/bin/env python
"""S11 tomo/eq (PLAN.md §6.5): differential Wadati — κ of the layer between 25 km and the LFE sources.

  python s11_eq_diff.py [--kappa eq_kappa]
Inputs: tomo/eq/tables/<kappa>_zones.csv + <kappa>_rows.csv (earthquake κ), tomo/tables (LFE κ, LFE picks for the
geometry), CVM v1.7 regional profiles. Writes tomo/eq/tables/diff_layer.csv, diff_crosscheck.csv, diff_model.json.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
TT = f"{ROOT}/svi_science/results/s11_depth/tomo"
EQ = f"{TT}/eq"
PROF = f"{ROOT}/data/velocity_models/usgs_cascadia_v1.7/profiles"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s11_tomo_wadati import along_strike  # noqa: E402

Z_CRUST = 25.0
REG = {"SVI": ("south_VI", "EallV3"), "WA_NW": ("royer_NW_WA", "BV3"), "WA_SW": ("royer_SW_WA", "CV3")}
# LFE κ headline values (tomo/RESULTS.md; clipped tP-on-tS, all families, passes Amendment 3) + columns
LFE = {"SVI": [("era clipped tP-on-tS (headline)", 1.773, 0.010), ("era inlier profile", 1.825, 0.042)],
       "WA_NW": [("V3 clipped tP-on-tS (headline)", 1.950, 0.029), ("frozen-P inlier profile", 1.815, 0.047),
                 ("eye-real 33 clipped", 1.909, 0.039)],
       "WA_SW": [("V3 clipped tP-on-tS (headline)", 1.904, 0.015), ("V3 inlier profile", 1.910, 0.024),
                 ("frozen-P clipped", 1.880, 0.024)]}
KAPPA_RF = (2.35, 2.8)
KAPPA_SLAB = (1.73, 1.80)


def load_profile(name):
    d = pd.read_csv(f"{PROF}/{name}.csv")
    d["zbot_km"] = d.zbot_km.replace(np.inf, 400.0).astype(float)
    return d


def ray_partials(prof, zs, X, zcuts=(Z_CRUST,)):
    """Direct up-going P ray, flat layers: per-row travel time in [0, zcut0), [zcut0, zs) (and below for deeper cuts).
    Returns array (len(X), len(zcuts)+1): P time in each depth band; band edges 0, zcuts..., zs."""
    edges = sorted(set([0.0, zs] + list(np.unique(np.r_[prof.ztop_km.values, prof.zbot_km.values]))))
    edges = [e for e in edges if e <= zs] + ([zs] if zs not in edges else [])
    edges = sorted(set(edges + [c for c in zcuts if c < zs]))
    h = np.diff(edges)
    zm = 0.5 * (np.array(edges[:-1]) + np.array(edges[1:]))
    v = np.array([prof.vp.values[(prof.ztop_km.values <= z) & (prof.zbot_km.values > z)][0] for z in zm])
    th = np.radians(np.linspace(0, 89.9, 6000))
    p = np.sin(th) / v.max()
    s = np.clip(p[:, None] * v[None, :], 0, 0.999999)
    cos = np.sqrt(1 - s ** 2)
    xl = h[None, :] * s / cos
    tl = h[None, :] / (v[None, :] * cos)
    xt = xl.sum(1)
    band = np.digitize(zm, list(zcuts))           # 0: above first cut, ...
    out = np.zeros((len(X), len(zcuts) + 1))
    for b in range(len(zcuts) + 1):
        tb = tl[:, band == b].sum(1)
        out[:, b] = np.interp(X, xt, tb)
    return out


def e1_kappa(fid, tp, ts):
    df = pd.DataFrame(dict(f=fid, x=ts, y=tp))
    xd = df.x - df.groupby("f").x.transform("mean"); yd = df.y - df.groupby("f").y.transform("mean")
    return float((xd ** 2).sum() / (xd * yd).sum())


_DM = {}


def _demeaned(parts, fid):
    """within-family demeaned band partials (cached by object identity of parts)."""
    key = (id(parts), len(parts))
    if key not in _DM:
        df = pd.DataFrame(parts)
        Pd = (df - df.groupby(np.asarray(fid)).transform("mean")).values
        _DM[key] = (parts, Pd)          # keep parts alive so id() stays unique
    return _DM[key][1]


def model_kappa(parts, fid, kap_bands):
    """E1 (tP on tS, family intercepts) κ of model times tP = Σ parts, tS = Σ κ_b parts_b."""
    Pd = _demeaned(parts, fid)
    tpd = Pd.sum(1)
    tsd = Pd @ np.asarray(kap_bands, float)
    return float((tsd @ tsd) / (tpd @ tsd))


def solve_layer(parts, fid, kc, kobs, kbelow=None):
    """κ_layer such that the model E1 κ on this geometry equals kobs (bands: crust, layer[, below])."""
    def f(kl):
        kb = [kc, kl] + ([kbelow] if kbelow is not None else [])
        return model_kappa(parts, fid, kb) - kobs
    lo, hi = 0.8, 5.0
    flo, fhi = f(lo), f(hi)
    if np.sign(flo) == np.sign(fhi):
        return np.nan
    for _ in range(40):
        m = 0.5 * (lo + hi); fm = f(m)
        if np.sign(fm) == np.sign(flo):
            lo, flo = m, fm
        else:
            hi = m
    return 0.5 * (lo + hi)


def lfe_geometry(reg):
    _, pset = REG[reg]
    P = pd.read_csv(f"{TT}/tables/picks_{pset}.csv", low_memory=False)
    P = P[P.tP.notna() & P.tS.notna() & (P.rP >= 1.2) & (P.rS >= 1.2)].copy()
    P = P[(P.depth_km >= 25) & (P.depth_km <= 50)]
    if "epi_km" not in P or P.epi_km.isna().all():
        import s11_tomo_profile as PR
        xy = PR.coords()
        la = P.sta.map(lambda s: xy.get(s, (np.nan, np.nan))[0]); lo = P.sta.map(lambda s: xy.get(s, (np.nan, np.nan))[1])
        P["epi_km"] = np.hypot((P.lon - lo) * 111.19 * np.cos(np.radians(P.lat)), (P.lat - la) * 111.19)
    P = P[P.epi_km.notna()]
    P["along"] = along_strike(P.lat.values, P.lon.values)
    P["zone"] = (np.floor(P.along / 30) * 30).astype(int)
    return P[["fam", "sta", "depth_km", "epi_km", "zone", "along"]].assign(fam=lambda d: d.fam.astype(str))


def parts_for(prof, geo, zcol="depth_km", xcol="epi_km", zcuts=(Z_CRUST,)):
    out = np.zeros((len(geo), len(zcuts) + 1))
    for zs, idx in geo.groupby(np.round(geo[zcol].values, 0)).groups.items():
        ii = geo.index.get_indexer(idx)
        out[ii] = ray_partials(prof, max(zs, max(zcuts) + 0.5), geo[xcol].values[ii], zcuts)
    return out


def get_eq(Z, kind, reg, zone="all", est="clip_inv", variant="d100"):
    # zone labels in eq_kappa_zones are "all" or signed ints ("+30", "-30", "+0"); fixed 2026-10-13 (was str(zone))
    lab = zone if isinstance(zone, str) else f"{int(zone):+d}"
    r = Z[(Z.kind == kind) & (Z.reg == reg) & (Z.zone.astype(str) == lab) & (Z.est == est) & (Z.variant == variant)]
    return r.iloc[0] if len(r) else None


def lfe_zone(reg, zone):
    _, pset = REG[reg]
    out = []
    for f, est in ((f"pooled_{pset}_inv.csv", "clip"), ):
        d = pd.read_csv(f"{TT}/tables/{f}")
        if "robust" in d:
            d = d[d.robust == True]  # noqa: E712
        lab = f"along_{zone:+d}_{zone + 30:+d}"
        r = d[d.zone == lab]
        if len(r):
            r = r.iloc[0]
            out.append((est, r.kappa, r.sig, bool(r["pass"])))
    pf = "profile_E_zones.csv" if reg == "SVI" else "profile_V3wa_zones.csv"
    d = pd.read_csv(f"{TT}/tables/{pf}")
    d = d[(d.set == ("EallV3" if reg == "SVI" else pset)) & (d.zone == f"along_{zone:+d}_{zone + 30:+d}")]
    if len(d) and np.isfinite(d.iloc[0].kappa):
        out.append(("profile", d.iloc[0].kappa, d.iloc[0].sig, bool(d.iloc[0]["pass"])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", default="eq_kappa")
    ap.add_argument("--ndraw", type=int, default=2000)
    ap.add_argument("--prefix", default="", help="output file prefix (tests)")
    a = ap.parse_args()
    rng = np.random.default_rng(5)
    Z = pd.read_csv(f"{EQ}/tables/{a.kappa}_zones.csv")
    R = pd.read_csv(f"{EQ}/tables/{a.kappa}_rows.csv", low_memory=False)
    sel = pd.read_csv(f"{EQ}/selection.csv")
    R = R.merge(sel[["cc_id", "slab2_km"]].drop_duplicates("cc_id"), on="cc_id", how="left")
    rows, cross, meta = [], [], {}
    for reg, (pname, pset) in REG.items():
        prof = load_profile(pname)
        geo = lfe_geometry(reg).reset_index(drop=True)
        parts = parts_for(prof, geo)
        zs_med = float(geo.depth_km.median())
        # vertical-ray fraction below 25 km (simple layered column)
        pv = ray_partials(prof, zs_med, np.array([0.0]))[0]
        f_vert = pv[1] / pv.sum()
        # geometry-weighted sensitivity: d κ_model / d κ_layer at the CVM values
        kc0 = float(np.average(prof.vpvs, weights=np.clip(np.minimum(prof.zbot_km, Z_CRUST) - prof.ztop_km, 0, None)))
        g = (model_kappa(parts, geo.fam.values, [kc0, 1.85]) - model_kappa(parts, geo.fam.values, [kc0, 1.75])) / 0.1
        meta[reg] = dict(profile=pname, n_lfe_rows=len(geo), n_lfe_fam=geo.fam.nunique(), zs_med=zs_med,
                         epi_med=float(geo.epi_km.median()), f_vertical=float(f_vert), w_layer_wadati=float(g),
                         kappa_crust_cvm=kc0)
        for scope in ["all"] + sorted(geo.zone.unique().tolist()):
            if scope == "all":
                gi = np.ones(len(geo), bool); lfe_vals = LFE[reg]
            else:
                gi = (geo.zone == scope).values
                lfe_vals = [(f"zone {e}", k, s, p) for e, k, s, p in lfe_zone(reg, scope)]
            if gi.sum() < 10 or geo.fam[gi].nunique() < 3:
                continue
            Pg, Fg = parts[gi], geo.fam.values[gi]
            for lv in lfe_vals:
                lab, kl_obs, sl = lv[0], lv[1], lv[2]
                lfe_pass = lv[3] if len(lv) > 3 else lab.endswith("(headline)") or "clipped" in lab
                for est in ("clip_inv", "profile_0.25", "profile_0.5", "clip_ols"):
                    e = get_eq(Z, "crust", reg, scope, est)
                    if e is None or not np.isfinite(e.get("kappa", np.nan)):
                        continue
                    kc, sc = float(e.kappa), float(e.sig) if np.isfinite(e.sig) else 0.05
                    kl = solve_layer(Pg, Fg, kc, kl_obs)
                    dr = [solve_layer(Pg, Fg, kc + sc * rng.standard_normal(),
                                      kl_obs + sl * rng.standard_normal()) for _ in range(a.ndraw if scope == "all" else 400)]
                    dr = np.array(dr)
                    simple = (kl_obs - (1 - f_vert) * kc) / f_vert
                    kbg = float(np.average(prof.vpvs, weights=np.clip(np.minimum(prof.zbot_km, zs_med) -
                                                                      np.maximum(prof.ztop_km, Z_CRUST), 0, None)))
                    hfrac = [(kl - kbg) / (krf - kbg) for krf in KAPPA_RF]
                    rows.append(dict(reg=reg, scope=str(scope), lfe_est=lab, kappa_lfe=kl_obs, sig_lfe=sl,
                                     lfe_pass=lfe_pass, eq_est=est, kappa_crust=kc, sig_crust=sc,
                                     crust_pass=bool(e["pass"]) if "pass" in e and pd.notna(e["pass"]) else False,
                                     n_eq=e.get("n_fam"), kappa_layer=kl, kl_p16=np.nanpercentile(dr, 16),
                                     kl_p84=np.nanpercentile(dr, 84), kl_sd=np.nanstd(dr),
                                     kappa_layer_simple=simple, f_vertical=f_vert, w_layer_wadati=g,
                                     layer_top=Z_CRUST, layer_bot=float(geo.depth_km[gi].median()),
                                     kappa_bg_cvm=kbg, h_frac_rf235=hfrac[0], h_frac_rf280=hfrac[1],
                                     h_km_rf235=hfrac[0] * (float(geo.depth_km[gi].median()) - Z_CRUST),
                                     h_km_rf280=hfrac[1] * (float(geo.depth_km[gi].median()) - Z_CRUST)))
        # cross-checks: intraslab and layer events, predicted from (κ_crust, κ_layer, κ_slab) on their own geometry
        eqc = get_eq(Z, "crust", reg)
        head = [r for r in rows if r["reg"] == reg and r["scope"] == "all" and r["eq_est"] == "clip_inv"
                and r["lfe_est"].endswith("(headline)")]
        if eqc is None or not head:
            continue
        kc, kl = float(eqc.kappa), head[0]["kappa_layer"]
        for kind in ("slab", "layer"):
            q = R[(R.kind == kind) & (R.reg == reg) & (R.dist_km <= 100)].copy()
            q = q.groupby("fam").filter(lambda g: len(g) >= 4).reset_index(drop=True)
            eo = get_eq(Z, kind, reg)
            if len(q) < 20 or eo is None:
                continue
            q["zev"] = q.o_depth_km.fillna(q.cat_depth_km)
            if kind == "slab":
                q["zint"] = q.slab2_km.clip(lower=Z_CRUST + 1)
                preds = {}
                for ks in KAPPA_SLAB:
                    tp_all, ts_all = [], []
                    for i, r in q.iterrows():
                        pr = ray_partials(prof, r.zev, np.array([r.dist_km]), (Z_CRUST, min(r.zint, r.zev - 0.5)))[0]
                        tp_all.append(pr.sum()); ts_all.append(pr[0] * kc + pr[1] * kl + pr[2] * ks)
                    preds[ks] = e1_kappa(q.fam.values, np.array(tp_all), np.array(ts_all))
                cross.append(dict(reg=reg, kind=kind, n_ev=q.fam.nunique(), kappa_obs=eo.kappa, sig_obs=eo.sig,
                                  obs_pass=eo["pass"], pred_kslab_173=preds[KAPPA_SLAB[0]],
                                  pred_kslab_180=preds[KAPPA_SLAB[1]], kappa_crust=kc, kappa_layer=kl,
                                  zev_med=float(q.zev.median()), zint_med=float(q.zint.median())))
            else:
                pr = parts_for(prof, q, "zev", "dist_km")
                pred = model_kappa(pr, q.fam.values, [kc, kl])
                cross.append(dict(reg=reg, kind=kind, n_ev=q.fam.nunique(), kappa_obs=eo.kappa, sig_obs=eo.sig,
                                  obs_pass=eo["pass"], pred=pred, kappa_lfe=head[0]["kappa_lfe"], kappa_crust=kc,
                                  kappa_layer=kl, zev_med=float(q.zev.median())))
    D = pd.DataFrame(rows)
    D.to_csv(f"{EQ}/tables/{a.prefix}diff_layer.csv", index=False)
    pd.DataFrame(cross).to_csv(f"{EQ}/tables/{a.prefix}diff_crosscheck.csv", index=False)
    json.dump(meta, open(f"{EQ}/tables/{a.prefix}diff_model.json", "w"), indent=1)
    with pd.option_context("display.width", 250, "display.max_rows", 400):
        print(json.dumps(meta, indent=1))
        if len(D):
            print(D[["reg", "scope", "lfe_est", "kappa_lfe", "eq_est", "kappa_crust", "crust_pass", "kappa_layer",
                     "kl_p16", "kl_p84", "kappa_layer_simple", "h_km_rf235", "h_km_rf280"]].round(3).to_string())
        print(pd.DataFrame(cross).round(3).to_string())


if __name__ == "__main__":
    main()
