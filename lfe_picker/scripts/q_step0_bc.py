#!/usr/bin/env python
"""P3 STEP 0(b) + 0(c) — the corrected event-level null, and the loud-day collapse table.

(b) Merlin's flaw #2.  The pre-registered moveout-shuffle null compared COUNTS of detections.  The
    honest statistic for a back-projection detector is the shuffled run's OWN chance-corrected
    Bostock recall, divided by the true run's: a wrong travel-time table that still recovers most of
    the events is not carrying geometry.  Reported alongside the uniform-over-grid expectation for
    the "fraction of detections > 5 km from every known family" number, which is what that 0.61
    actually equals.

(c) Merlin's flaw #3.  Per-cell per-day two-sided MAD normalisation collapses on the loudest ETS
    days.  Table: Bostock events vs P2 detections vs frozen-template detections, per day, on the
    loud days (2011 doy 239-242; 2012 doy 253, 255, 267).

Reads stored products only; nothing is re-run, nothing is tuned, no referee budget is consumed.

    python scripts/q_step0_bc.py
"""
from __future__ import annotations
import os, sys, json, glob
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

P3 = f"{P.RES}/p3"
OUT = f"{P.RES}/p2"
SHIFT = 307.0
TOL = 2.0


def one2one(det_t, ev_t, tol=TOL):
    """identical to p2_score.one2one"""
    out = np.full(len(det_t), -1, np.int64)
    if len(det_t) == 0 or len(ev_t) == 0:
        return out
    order = np.argsort(det_t)
    used = np.zeros(len(ev_t), bool)
    ev_sorted = np.argsort(ev_t)
    ev_t_s = ev_t[ev_sorted]
    for i in order:
        j = np.searchsorted(ev_t_s, det_t[i])
        best, bd = -1, tol
        for k in range(max(0, j - 3), min(len(ev_t_s), j + 3)):
            if used[k]:
                continue
            d = abs(ev_t_s[k] - det_t[i])
            if d <= bd:
                best, bd = k, d
        if best >= 0:
            used[best] = True
            out[i] = ev_sorted[best]
    return out


def recall_cc(det_t, ev_t, nev):
    m = one2one(det_t, ev_t)
    mp = one2one(det_t + SHIFT, ev_t)
    mm = one2one(det_t - SHIFT, ev_t)
    rec = (m >= 0).sum() / nev
    ch = 0.5 * ((mp >= 0).sum() + (mm >= 0).sum()) / nev
    return float(rec), float(ch), float((rec - ch) / max(1e-9, 1 - ch))


def main():
    os.makedirs(P3, exist_ok=True)
    B = P.bostock()
    FAM = P.families()

    # ---------------------------------------------------------------- (b) corrected null
    rows = []
    for year, tag, null in [(2011, "p2nnpu_2011", "p2nnpu_2011_shuf"),
                            (2012, "p2nnpu_2012", "p2nnpu_2012_shuf")]:
        sc = json.load(open(f"{OUT}/score_{tag}.json"))
        ED = sc["ets_days"]; QD = sc["quiet_days"]
        thr = {float(k): v for k, v in sc["thresholds"].items()}
        D = pd.read_parquet(f"{OUT}/det_{tag}.parquet")
        N = pd.read_parquet(f"{OUT}/det_{null}.parquet")
        Bx = B[B.year == year]
        ev = Bx[Bx.doy.isin(ED)].t_origin.to_numpy()
        nev = len(ev)
        for r, z in sorted(thr.items()):
            d = D[(D.z >= z) & D.doy.isin(ED)]
            n = N[(N.z >= z) & N.doy.isin(ED)]
            # the null at its OWN quiet-day FA calibration (fair threshold)
            nq = np.sort(N[N.doy.isin(QD)].z.values)[::-1]
            k = int(round(r * max(1, len(QD))))
            zn = float(nq[min(k, len(nq)) - 1]) if len(nq) else np.inf
            n2 = N[(N.z >= zn) & N.doy.isin(ED)]
            tr, tc, tcc = recall_cc(d.t_origin.values, ev, nev)
            sr, scn, scc = recall_cc(n.t_origin.values, ev, nev)
            s2r, s2c, s2cc = recall_cc(n2.t_origin.values, ev, nev)
            rows.append(dict(year=year, fa_per_day=r, z_true=z, z_null_own_fa=zn,
                             n_events=nev, n_det_true=len(d), n_det_null=len(n),
                             n_det_null_own_fa=len(n2),
                             recall_true=tr, chance_true=tc, recall_cc_true=tcc,
                             recall_null=sr, chance_null=scn, recall_cc_null=scc,
                             recall_cc_null_own_fa=s2cc,
                             count_ratio=len(n) / max(1, len(d)),
                             RECALL_RATIO_cc=scc / max(1e-9, tcc),
                             RECALL_RATIO_cc_own_fa=s2cc / max(1e-9, tcc)))
            print(f"[{year} {r}/day] true cc {tcc:.4f}  shuffled cc {scc:.4f} "
                  f"-> ratio {scc/max(1e-9,tcc):.3f}   (count ratio {len(n)/max(1,len(d)):.3f})",
                  flush=True)
    NUL = pd.DataFrame(rows)
    NUL.to_csv(f"{P3}/p2_corrected_event_null.csv", index=False)

    # ---------------------------------------------------------------- uniform-grid expectation
    G = pd.read_csv(f"{OUT}/grid_cells.csv")
    fa = FAM[["lat", "lon"]].to_numpy()
    dmin = np.min(np.stack([P.epi_km(G.lat.values, G.lon.values, fa[i, 0], fa[i, 1])
                            for i in range(len(fa))], 1), 1)
    uni = dict(n_cells=int(len(G)),
               frac_gt5km=float((dmin > 5).mean()), frac_gt10km=float((dmin > 10).mean()),
               median_km=float(np.median(dmin)))
    # observed, at 6/day on 2011 and 2012
    for year, tag in [(2011, "p2nnpu_2011"), (2012, "p2nnpu_2012")]:
        sc = json.load(open(f"{OUT}/score_{tag}.json"))
        uni[f"observed_{year}_frac_gt5km"] = sc["rate_6.0"]["frac_gt5km_from_every_family"]
        uni[f"observed_{year}_frac_gt10km"] = sc["rate_6.0"]["frac_gt10km_from_every_family"]
        dd = pd.read_parquet(f"{OUT}/det_{tag}.parquet")
        dd = dd[(dd.z >= sc["thresholds"]["6.0"]) & dd.doy.isin(sc["ets_days"])]
        uni[f"observed_{year}_n_argmax_cells"] = int(dd.cell.nunique())
    json.dump(uni, open(f"{P3}/p2_uniform_grid_expectation.json", "w"), indent=1)
    print("\n[uniform grid]", json.dumps(uni, indent=1))

    # ---------------------------------------------------------------- (c) loud-day collapse
    LOUD = {2011: [239, 240, 241, 242], 2012: [253, 255, 267]}
    tmpl = {}
    t11 = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f"{P.RES}/trackO/detections_2011-*.parquet"))],
                    ignore_index=True)
    tmpl[2011] = t11[t11.cat_core].groupby("doy").size()
    t12 = pd.read_parquet(f"{P.EXT}/results/e5r_detections_a1_mergeB_L12.parquet")
    t12 = t12[(t12.L == 12.0) & (t12.net_mad >= 8.75)]
    t12["doy"] = pd.to_datetime(t12.t, unit="s").dt.dayofyear
    t12["year"] = pd.to_datetime(t12.t, unit="s").dt.year
    tmpl[2012] = t12[t12.year == 2012].groupby("doy").size()
    crows = []
    for year in (2011, 2012):
        sc = json.load(open(f"{OUT}/score_p2nnpu_{year}.json"))
        z6 = sc["thresholds"]["6.0"]
        D = pd.read_parquet(f"{OUT}/det_p2nnpu_{year}.parquet")
        d6 = D[D.z >= z6].groupby("doy").size()
        Bx = B[B.year == year]
        nb = Bx.groupby("doy").size()
        ED = sc["ets_days"]
        for doy in ED:
            crows.append(dict(year=year, doy=doy, loud=doy in LOUD[year],
                              bostock=int(nb.get(doy, 0)),
                              p2_det_6perday=int(d6.get(doy, 0)),
                              template_det_core=int(tmpl[year].get(doy, 0))))
    C = pd.DataFrame(crows)
    C["p2_per_1000_bostock"] = 1000 * C.p2_det_6perday / C.bostock.clip(lower=1)
    C["template_per_1000_bostock"] = 1000 * C.template_det_core / C.bostock.clip(lower=1)
    C.to_csv(f"{P3}/loud_day_collapse.csv", index=False)
    print("\n[loud-day collapse] loud days")
    print(C[C.loud].to_string(index=False))
    print("\n[loud-day collapse] summary loud vs rest")
    print(C.groupby(["year", "loud"])[["bostock", "p2_det_6perday", "template_det_core",
                                       "p2_per_1000_bostock", "template_per_1000_bostock"]]
          .mean().round(2).to_string())
    print(f"\n[done] -> {P3}/")


if __name__ == "__main__":
    main()
