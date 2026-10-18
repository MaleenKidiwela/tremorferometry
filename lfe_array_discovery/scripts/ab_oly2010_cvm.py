#!/usr/bin/env python
"""Relocate the 2010 N-Olympic arrays-as-arrays candidate families (real + whole-chain site-shift
null) with Hyp2000, models USGS_NW (USGS Cascadia CVM v1.7, royer_NW_WA profile, LVZ merged,
Vp/Vs 1.823) and GSC (Vp/Vs 1.73), pick ratio R 1.5, free depth -- and score both beside the
adopted L9_xcdiff locations (results/u_locate/wa/wa_{real,null}_ranked.csv).  2026-09-28.

Inputs are read-only: the 44 s / 12 s-pre-node 1-8 Hz 20-sps long stacks that u_locate built from
the frozen det_it0 detections (results/u_locate/wa_{real,null}_stacks.npz); site coordinates
results/o_oly2010/sites.csv (beam = array centroid).  No detection or growth is re-run.

Picks = ab_hypo.make_picks (Hilbert-rotated, zero-phase extrema; FS set to 20).  Hyp2000 inputs =
ab_hypo.write_hyp_inputs / run_hyp (station codes mapped to <=5 chars).

Scoring (every family keeps its row; no cut):
  * u_locate's Chestler scorer VERBATIM (uu_wa.summarise: nearest-coordinate distance, families
    within 5/10/20 km, coordinates hit within 5/10/20 km), real beside the site-shift null;
  * one-to-one (Hungarian) matching of the 69 Chestler-derived coordinates, 3-D and epicentral,
    fraction within 3/5/10 km, plus the rotated+shifted-location null of ab_truth_score;
  * Royer & Bostock (2014) NW digitized Hyp2000 templates via ab_truth_score.score (subarray NW).
Rows: real[all], real[280 rich: dedupe_n_partners>=20], null[all], null[top 280 by partners].

  python scripts/ab_oly2010_cvm.py --out results/ab_royer/OLY2010_cvm
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_hypo_cvm                                                          # noqa: E402,F401
import ab_hypo as H                                                         # noqa: E402
import ab_common as A                                                       # noqa: E402
import uu_common as U                                                       # noqa: E402
import uu_wa as W                                                           # noqa: E402
import ab_truth_score as TS                                                 # noqa: E402

RES = U.RES
UL = f"{RES}/u_locate"
MODELS = [("USGS_NW", 1.823), ("GSC", 1.73)]
RMIN = 1.5
NNULL = 1000


def sites():
    S = pd.read_csv(f"{RES}/o_oly2010/sites.csv")
    C = pd.read_csv(f"{RES}/channels_2006_2012.csv").dropna(subset=["lat", "lon"])
    el = C.groupby("sta").elev.median()
    code = {}
    for r in S.itertuples():
        c = r.site.upper() if len(r.site) <= 5 else ("XPB" + r.site[-2:])
        code[r.site] = c
        mem = str(r.members).split(",")
        e = float(np.nanmedian([el.get(m, np.nan) for m in mem])) if r.kind == "ARRAY" \
            else float(el.get(r.site, 0.0))
        A.COORD[c] = (float(r.lat), float(r.lon), 0.0 if not np.isfinite(e) else e)
        A.NET[c] = r.net
    assert len(set(code.values())) == len(code)
    return code


def load_stacks(p, code):
    z = np.load(p)
    ST = {}
    for k in z.files:
        if k == "meta":
            continue
        g, s = k.split("_", 1)
        ST.setdefault(int(g), {})[code[s]] = z[k].astype(np.float64)
    return ST


def one_to_one(L, C, rng):
    """Hungarian match of the published coordinates to located families; 3-D and epi."""
    ok = L[L.located.astype(bool)]
    out = {}
    if not len(ok):
        return out
    lat0, lon0 = C.latitude.mean(), C.longitude.mean()
    tx, ty = TS.to_xy(C.longitude.values, C.latitude.values, lon0, lat0)
    ox, oy = TS.to_xy(ok.lon.values, ok.lat.values, lon0, lat0)
    tz, oz = C.depth.values, ok.dep.values

    def m(ox, oy):
        De = np.hypot(tx[:, None] - ox[None], ty[:, None] - oy[None])
        D3 = np.sqrt(De ** 2 + (tz[:, None] - oz[None]) ** 2)
        r = {}
        for tag, D in (("3d", D3), ("epi", De)):
            i, j = linear_sum_assignment(D)
            d = np.full(len(tx), np.inf)
            d[i] = D[i, j]
            r[tag] = d
        return r

    real = m(ox, oy)
    R = np.sqrt(np.mean((tx - tx.mean()) ** 2 + (ty - ty.mean()) ** 2))
    mx, my = ox.mean(), oy.mean()
    nul = {t: np.zeros((NNULL, 3)) for t in ("3d", "epi")}
    for k in range(NNULL):
        th, rr, ph = rng.uniform(0, 2 * np.pi), R * np.sqrt(rng.uniform()), rng.uniform(0, 2 * np.pi)
        x0, y0 = ox - mx, oy - my
        nx = np.cos(th) * x0 - np.sin(th) * y0 + mx + rr * np.cos(ph)
        ny = np.sin(th) * x0 + np.cos(th) * y0 + my + rr * np.sin(ph)
        dn = m(nx, ny)
        for t in nul:
            nul[t][k] = [np.sum(dn[t] <= r) for r in (3, 5, 10)]
    for t in ("3d", "epi"):
        for i, r in enumerate((3, 5, 10)):
            v = int(np.sum(real[t] <= r))
            out[f"o2o_{t}_{r}km"] = v
            out[f"o2o_{t}_{r}km_rotnull_mean"] = float(nul[t][:, i].mean())
            out[f"o2o_{t}_{r}km_rotnull_p"] = float((1 + np.sum(nul[t][:, i] >= v)) / (NNULL + 1))
        out[f"o2o_{t}_median_km"] = float(np.median(real[t]))
    return out


def depth_stats(L):
    ok = L[L.located.astype(bool)]
    d = dict(n_located=int(len(ok)))
    if len(ok):
        d.update(dep_med=float(ok.dep.median()),
                 frac_25_45=float(ok.dep.between(25, 45).mean()),
                 dep_p10=float(ok.dep.quantile(.1)), dep_p90=float(ok.dep.quantile(.9)))
        for c in ("erh", "erz", "rms", "gap", "num"):
            if c in ok:
                d[f"{c}_med"] = float(ok[c].median())
    return d


def royer_nw(L):
    ok = L[L.located.astype(bool)].rename(columns={"dep": "depth"})
    T = pd.read_csv(TS.TRUTH)
    T = T[T.locator == "Hyp2000"]
    if not len(ok):
        return {}
    r = TS.score(T, ok[["lon", "lat", "depth"]], ["NW"], NNULL, 0).iloc[0].to_dict()
    return {f"royerNW_{k}": v for k, v in r.items() if k != "subarray"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{RES}/ab_royer/OLY2010_cvm")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    A.FS = 20.0                                    # u_locate long stacks are 20 sps
    code = sites()
    json.dump(code, open(f"{a.out}/station_codes.json", "w"), indent=1)
    stas = sorted(set(code.values()))
    C = pd.read_csv(U.CHESTLER)
    allP, locs = [], {}
    for pop in ("real", "null"):
        ST = load_stacks(f"{UL}/wa_{pop}_stacks.npz", code)
        P = H.make_picks(ST)
        P.insert(0, "population", pop)
        allP.append(P)
        U.log(f"[{pop}] {len(ST)} families, {len(P)} station stacks, "
              f"P {P.tP.notna().sum()} S {P.tS.notna().sum()} "
              f"(R>={RMIN}: P {(P.rP >= RMIN).sum()} S {(P.rS >= RMIN).sum()})")
        base = pd.read_csv(f"{UL}/wa/wa_{pop}_ranked.csv")
        keep = ["grp", "dedupe_n_partners", "sf_top", "n_sta", "located", "lat", "lon", "dep",
                "edge", "boot_km", "boot_dz"]
        B = base[[c for c in keep if c in base]].rename(
            columns={c: f"L9_{c}" for c in ("n_sta", "located", "lat", "lon", "dep", "edge",
                                            "boot_km", "boot_dz")})
        locs[(pop, "L9")] = B.rename(columns={"L9_located": "located", "L9_lat": "lat",
                                              "L9_lon": "lon", "L9_dep": "dep",
                                              "L9_edge": "edge"})
        for model, vpvs in MODELS:
            wd = f"{a.out}/run_{pop}_{model}"
            ids = H.write_hyp_inputs(wd, P.drop(columns="population"), stas, vpvs, RMIN,
                                     model=model)
            Hy, lg = H.run_hyp(wd)
            U.log(f"[{pop} {model} {vpvs}] Hyp2000 {len(Hy)}/{len(ids)} located"
                  + ("" if len(Hy) else "\n" + lg))
            M = B[["grp", "dedupe_n_partners", "sf_top"]].merge(Hy, on="grp", how="left")
            M["located"] = M.lat.notna()
            M["n_input"] = M.grp.isin(ids)
            M["edge"] = False
            M["model"], M["vpvs"], M["rmin"], M["population"] = model, vpvs, RMIN, pop
            locs[(pop, model)] = M
            if pop == "real":
                M.to_csv(f"{a.out}/loc_{model}.csv", index=False)
            else:
                M.to_csv(f"{a.out}/loc_{model}_null.csv", index=False)
    pd.concat(allP, ignore_index=True).to_csv(f"{a.out}/picks.csv", index=False)
    # compare.csv -- per real family, all three locations side by side (plus the null file)
    for pop in ("real", "null"):
        Cm = locs[(pop, "L9")].copy()
        for model, _ in MODELS:
            M = locs[(pop, model)][["grp", "lat", "lon", "dep", "erh", "erz", "rms", "num",
                                    "gap", "located"]]
            Cm = Cm.merge(M.rename(columns={c: f"{model}_{c}" for c in M.columns if c != "grp"}),
                          on="grp", how="left")
        for model, _ in MODELS:
            ok = Cm.located.astype(bool) & Cm[f"{model}_located"].astype(bool)
            Cm[f"{model}_vs_L9_epi_km"] = np.where(ok, U.epi_km(
                Cm.lat, Cm.lon, Cm[f"{model}_lat"], Cm[f"{model}_lon"]), np.nan)
            Cm[f"{model}_vs_L9_dz_km"] = np.where(ok, Cm[f"{model}_dep"] - Cm.dep, np.nan)
        Cm.to_csv(f"{a.out}/compare{'' if pop == 'real' else '_null'}.csv", index=False)
    # score table
    rng = np.random.default_rng(0)
    rows = []
    for loc in ("L9", "GSC", "USGS_NW"):
        R, N = locs[("real", loc)], locs[("null", loc)]
        n280 = int((R.sf_top == True).sum())                               # noqa: E712
        subs = {"real[all]": R, "real[280 rich]": R[R.sf_top == True],     # noqa: E712
                "null[all]": N,
                "null[top 280 by partners]": N.sort_values("dedupe_n_partners",
                                                           ascending=False).head(n280)}
        for tag, T in subs.items():
            T = T.copy()
            T["located"] = T.located.fillna(False).astype(bool)
            T["edge"] = T.edge.fillna(False).astype(bool)
            d, _ = W.summarise(tag, T, C)
            d = {("u_" + k if k not in ("pop", "n_families", "n_located") else k): v
                 for k, v in d.items()}
            row = dict(locator=loc, **d)
            row.update(depth_stats(T))
            row.update(one_to_one(T, C, rng))
            row.update(royer_nw(T))
            rows.append(row)
            U.log(f"[score] {loc} {tag}: n_loc {row['n_located']} "
                  f"hit10 {row.get('u_coords_hit_10km')} o2o3d_10 {row.get('o2o_3d_10km')}")
    S = pd.DataFrame(rows)
    S.to_csv(f"{a.out}/score_summary.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 80):
        U.log("\n" + S.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
