#!/usr/bin/env python
"""SECOND location column for the Royer & Bostock (2014) reproduction: Hyp2000 with the depth
FIXED to the plate interface under the epicentre (the free-depth ab_hypo.py headline stays the
main result; this is a column beside it).

Inputs (unchanged from the headline): <hypo>/picks.csv, the headline model GSC onshore 1-D,
Vp/Vs 1.73, pick ratio R 1.5, same station pool, same Hyp2000 control file (ab_hypo.write_hyp_inputs
is reused as-is; only the per-event terminator line is rewritten).

Method -- fixed-point iteration (NOT a grid search):
  start epicentre = the free-depth headline epicentre (else the centroid of the picked stations);
  k = 0, 1, ...: z_k = slab depth under epicentre_k; locate with Hyp2000, depth fixed at z_k
  (terminator col 35 '-', trial lat/lon/depth in cols 15-34; verified on the USGS test event:
  depth held exactly, remark '-', same answer from two different trial epicentres);
  epicentre_{k+1} = the Hyp2000 solution T(x_k) for the first 10 iterations, then the relaxed
  x_k + 0.5 (T(x_k) - x_k) (breaks 2-cycles); converged when |T(x_k) - x_k| < 0.1 km (max 60 its);
  the reported location is T(x_k), i.e. a fixed point of 'locate at the slab depth below me'.
  Families whose epicentre leaves the slab grid are reported as off_grid.
Hyp2000 MIN 4 weighted phases is kept.  Free depth has 4 unknowns, so a family with exactly 4
phases is only exactly determined there; with depth fixed (3 unknowns) it is over-determined.
Those families are flagged free_zero_dof and scored separately.  Families with <= 3 phases are
counted but cannot be located (Hyp2000 MIN 4).

Slab models (data/slab_models/README.md): slab2 = Slab2 cas_slab2_dep_02.24.18.grd (headline of
this column); mccrory = McCrory et al. 2012 / USGS DS633 Modeled_slab_XYZ_Grid (5').

  python scripts/ab_hypo_slab.py --sub NW
-> results/ab_royer/<SUB>/hypo_slab/loc_slab_{slab2,mccrory}.csv, score_*.csv, summary.txt
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import RegularGridInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402
import ab_hypo as HY                                                        # noqa: E402
import ab_truth_score as TS                                                 # noqa: E402

SLAB_DIR = os.path.join(HERE, "../../data/slab_models")
RES = os.path.join(HERE, "../results/ab_royer")
MODEL, VPVS, RMIN = HY.HEAD
TOL_KM, MAXIT, N_UNDAMPED, ALPHA = 0.1, 60, 10, 0.5


def slab_interp(name):
    """Return f(lat, lon) -> slab-top depth (km, positive down), NaN off the grid."""
    if name == "slab2":
        d = xr.open_dataset(f"{SLAB_DIR}/slab2/cas_slab2_dep_02.24.18.grd")
        x, y, z = d.x.values, d.y.values, -d.z.values.astype(float)
        fi = RegularGridInterpolator((y, x), z, bounds_error=False, fill_value=np.nan)
        return lambda la, lo: fi(np.c_[np.atleast_1d(la), np.mod(np.atleast_1d(lo), 360)])
    if name == "mccrory":
        t = pd.read_csv(f"{SLAB_DIR}/mccrory2012/Modeled_slab_XYZ_Grid_GCS_WGS84.txt")
        t["lo"], t["la"] = t.Longitude.round(4), t.Latitude.round(4)
        g = t.pivot_table(index="la", columns="lo", values="Depth_km")
        fi = RegularGridInterpolator((g.index.values, g.columns.values), -g.values,
                                     bounds_error=False, fill_value=np.nan)
        return lambda la, lo: fi(np.c_[np.atleast_1d(la), np.atleast_1d(lo)])
    raise ValueError(name)


def terminator(g, la, lo, z):
    """Hypoinverse terminator (H71 terminator format 1): trial lat/lon/depth, depth FIXED ('-')."""
    t = [" "] * 72
    la_h = int(round(abs(la) * 6000))           # hundredths of a minute
    lo_h = int(round(abs(lo) * 6000))
    t[14:16] = f"{la_h // 6000:2d}"
    t[17:21] = f"{la_h % 6000:4d}"
    t[21:24] = f"{lo_h // 6000:3d}"
    t[25:29] = f"{lo_h % 6000:4d}"
    t[29:34] = f"{z:5.2f}"
    t[34] = "-"
    t[62:72] = f"{int(g):10d}"
    return "".join(t)


def locate_fixed(wd, P, stas, trial):
    """One Hyp2000 pass; trial = {grp: (lat, lon, z)}. Returns the out.sum table."""
    Pa = P[P.grp.isin(list(trial))]
    HY.write_hyp_inputs(wd, Pa, stas, VPVS, RMIN, model=MODEL)
    lines = open(f"{wd}/pha.arc").read().rstrip("\n").split("\n")
    out = []
    for ln in lines:
        s = ln.strip()
        if ln.startswith(" " * 62) and s.isdigit() and int(s) in trial:
            la, lo, z = trial[int(s)]
            ln = terminator(int(s), la, lo, z)
        out.append(ln)
    open(f"{wd}/pha.arc", "w").write("\n".join(out) + "\n")
    for f in ("out.sum", "out.prt", "out.arc"):
        if os.path.exists(f"{wd}/{f}"):
            os.remove(f"{wd}/{f}")
    H, _ = HY.run_hyp(wd)
    return H


def run_model(name, P, stas, H0, start, out):
    slab = slab_interp(name)
    wd = f"{out}/run_{name}"
    cur = dict(start)                                  # grp -> (lat, lon)
    done, hist = {}, []
    active = set(cur)
    for it in range(MAXIT):
        trial = {}
        for g in sorted(active):
            z = float(slab(*cur[g])[0])
            if np.isfinite(z):
                trial[g] = (cur[g][0], cur[g][1], z)
            else:
                done[g] = dict(grp=g, status="off_grid", n_iter=it)
        if not trial:
            break
        H = locate_fixed(wd, P, stas, trial).set_index("grp") if trial else pd.DataFrame()
        for g, (la, lo, z) in trial.items():
            if g not in H.index:
                done[g] = dict(grp=g, status="hyp_failed", n_iter=it + 1)
                continue
            h = H.loc[g]
            mv = float(A.epi_km(h.lat, h.lon, la, lo))
            rec = dict(grp=g, lat=h.lat, lon=h.lon, dep_fixed=z, hyp_dep=h.dep, erh=h.erh,
                       erz=h.erz, rms=h.rms, n_picks=int(h.num), gap=int(h.gap), ot_s=h.ot_s,
                       last_move_km=mv, n_iter=it + 1)
            # plain fixed-point step for the first N_UNDAMPED iterations, then relaxed
            # (x <- x + ALPHA (T(x) - x)) to break 2-cycles; convergence is always judged on the
            # UNDAMPED Hyp2000 move |T(x) - x| < TOL_KM and the Hyp2000 solution T(x) is reported
            al = 1.0 if it < N_UNDAMPED else ALPHA
            cur[g] = (la + al * (h.lat - la), lo + al * (h.lon - lo))
            if mv < TOL_KM:
                z2 = float(slab(h.lat, h.lon)[0])
                rec.update(status="converged", slab_dep_at_final=z2)
                done[g] = rec
            else:
                done[g] = dict(rec, status="not_converged")
        active = {g for g, r in done.items() if r["status"] == "not_converged"}
        hist.append((it + 1, len(trial), len(active)))
        A.log(f"[slab {name}] it {it + 1}: located {len(trial)}, still moving {len(active)}")
        if not active:
            break
    L = pd.DataFrame(list(done.values()))
    if "slab_dep_at_final" not in L:
        L["slab_dep_at_final"] = np.nan
    ok = L.lat.notna() if "lat" in L else pd.Series(False, index=L.index)
    L["slab_dep_at_final"] = np.where(L.status == "converged", L.slab_dep_at_final,
                                      [float(slab(a, b)[0]) if o else np.nan
                                       for a, b, o in zip(L.get("lat", np.nan), L.get("lon", np.nan), ok)])
    Hf = H0.set_index("grp")
    L["free_lat"] = [Hf.lat.get(g, np.nan) for g in L.grp]
    L["free_lon"] = [Hf.lon.get(g, np.nan) for g in L.grp]
    L["free_dep"] = [Hf.dep.get(g, np.nan) for g in L.grp]
    L["free_erh"] = [Hf.erh.get(g, np.nan) for g in L.grp]
    L["shift_vs_free_km"] = [float(A.epi_km(a, b, c, d)) if np.isfinite([a, b, c, d]).all()
                             else np.nan for a, b, c, d in zip(L.get("lat", np.nan),
                                                               L.get("lon", np.nan),
                                                               L.free_lat, L.free_lon)]
    L["slab_model"] = name
    return L.sort_values("grp"), hist


def score(df, sub, nnull=1000):
    truth = pd.read_csv(TS.TRUTH)
    truth = truth[truth.locator == "Hyp2000"]
    return TS.score(truth, df[["lat", "lon", "depth"]].dropna(subset=["lat"]), [sub], nnull, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", required=True, choices=["NW", "SW", "NVI"])
    ap.add_argument("--hypo", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--models", default="slab2,mccrory")
    ap.add_argument("--pool", default="")
    ap.add_argument("--nnull", type=int, default=1000)
    a = ap.parse_args()
    hypo = a.hypo or f"{RES}/{a.sub}/hypo"
    out = a.out or f"{RES}/{a.sub}/hypo_slab"
    os.makedirs(out, exist_ok=True)
    stas = A.pool_from_arg(a.pool, a.sub)
    P = pd.read_csv(f"{hypo}/picks.csv")
    H0 = pd.read_csv(f"{hypo}/loc_headline.csv")
    okP = P.tP.notna() & (P.rP >= RMIN)
    okS = P.tS.notna() & (P.rS >= RMIN)
    nph = (okP.astype(int) + okS.astype(int)).groupby(P.grp).sum()
    Hi = H0.set_index("grp")
    start = {}
    for g, n in nph.items():
        if n < 4:
            continue
        if g in Hi.index:
            start[int(g)] = (float(Hi.lat[g]), float(Hi.lon[g]))
        else:                                   # centroid of the stations with a usable pick
            ss = P.sta[(P.grp == g) & (okP | okS)]
            start[int(g)] = (float(np.mean([A.COORD[s][0] for s in ss])),
                             float(np.mean([A.COORD[s][1] for s in ss])))
    A.log(f"[slab {a.sub}] families {len(nph)}; >=4 phases {len(start)}; "
          f"==4 phases {int((nph == 4).sum())}; <=3 phases {int((nph < 4).sum())}; "
          f"free-depth located {len(H0)}")
    free = H0.rename(columns={"dep": "depth"})
    S_free = score(free, a.sub, a.nnull)
    rows = [dict(set="free-depth headline (GSC 1.73 R1.5)", n_located=len(H0),
                 median_erh=float(H0.erh.median()), median_rms=float(H0.rms.median()),
                 median_shift_km=0.0, **S_free.iloc[0].to_dict())]
    extra = []
    for name in a.models.split(","):
        L, hist = run_model(name, P, stas, H0, start, out)
        L["n_phases"] = [int(nph.get(g, 0)) for g in L.grp]
        L["free_zero_dof"] = L.n_phases == 4
        L.to_csv(f"{out}/loc_slab_{name}.csv", index=False)
        ok = L[L.status.isin(["converged", "not_converged"])].copy()
        okc = ok.rename(columns={"dep_fixed": "depth"})
        Sc = score(okc, a.sub, a.nnull)
        Sc.to_csv(f"{out}/score_{name}.csv", index=False)
        rows.append(dict(set=f"slab-fixed {name} (all >=4-phase)", n_located=len(ok),
                         median_erh=float(ok.erh.median()), median_rms=float(ok.rms.median()),
                         median_shift_km=float(ok.shift_vs_free_km.median()),
                         **Sc.iloc[0].to_dict()))
        od = okc[~okc.free_zero_dof]
        if len(od):
            So = score(od, a.sub, a.nnull)
            rows.append(dict(set=f"slab-fixed {name} (>=5 phases only)", n_located=len(od),
                             median_erh=float(od.erh.median()), median_rms=float(od.rms.median()),
                             median_shift_km=float(od.shift_vs_free_km.median()),
                             **So.iloc[0].to_dict()))
        z4 = ok[ok.free_zero_dof]
        extra.append(dict(model=name, n_status=L.status.value_counts().to_dict(),
                          iters=hist, n4_located=len(z4),
                          n4_median_erh=float(z4.erh.median()) if len(z4) else np.nan,
                          n4_median_erh_free=float(z4.free_erh.median()) if len(z4) else np.nan,
                          n4_median_shift=float(z4.shift_vs_free_km.median()) if len(z4) else np.nan,
                          median_fixed_dep=float(ok.dep_fixed.median()),
                          median_free_dep=float(ok.free_dep.median())))
    free_od = free[free.grp.map(nph) >= 5]
    if len(free_od):
        So = score(free_od, a.sub, a.nnull)
        rows.insert(1, dict(set="free-depth headline (>=5 phases only)", n_located=len(free_od),
                            median_erh=float(free_od.erh.median()),
                            median_rms=float(free_od.rms.median()), median_shift_km=0.0,
                            **So.iloc[0].to_dict()))
    T = pd.DataFrame(rows)
    T.to_csv(f"{out}/summary.csv", index=False)
    cols = ["set", "n_located", "median_erh", "median_rms", "median_shift_km", "n_ours_in_region",
            "median_match_km"] + [c for r in TS.RADII for c in
                                  (f"frac_{r:g}km", f"null_mean_{r:g}km", f"p_{r:g}km")]
    pd.set_option("display.width", 300)
    txt = [f"ab_hypo_slab.py {a.sub}: Hyp2000 depth FIXED to the slab top under the epicentre "
           f"(fixed-point iteration, tol {TOL_KM} km, max {MAXIT}); model {MODEL} Vp/Vs {VPVS} "
           f"R {RMIN}; truth = Royer 2014 Hyp2000 templates; null = our set rotated+shifted "
           f"({a.nnull}x), null_mean = the rotated-positions expectation.",
           f"families {len(nph)}; >=4 phases {len(start)} (== 4: {int((nph == 4).sum())}, "
           f"free-depth zero-dof); <=3 phases {int((nph < 4).sum())} (not locatable, MIN 4)",
           "", T[cols].round(3).to_string(index=False), ""]
    for e in extra:
        txt.append(str(e))
    open(f"{out}/summary.txt", "w").write("\n".join(txt) + "\n")
    print("\n".join(txt))


if __name__ == "__main__":
    main()
