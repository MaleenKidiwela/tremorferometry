#!/usr/bin/env python
"""P2 scoring — everything Merlin asked for, from stored back-projection detections.

  python scripts/p2_score.py --tag p2core_2012 --null p2core_2012_shuf --year 2012

Reported (all at a threshold set by QUIET-DAY false events/day ONLY):
  * own recall vs Bostock, +-2 s one-to-one, family-agnostic, with a TIME-SHIFTED chance control
    subtracted (a family-agnostic detector matching a 1,400-event ETS day by accident is the first
    thing that can fake this number);
  * network-wide false events/day on quiet days -> the operating points;
  * moveout-shuffle null ratio on ETS days;
  * recall by label provenance / n_sites / magnitude / family;
  * location error vs the matched event's Bostock family, and the CENTROID-COLLAPSE diagnostic:
    the distribution of back-projected cells against the 130 family positions and the fraction of
    detections whose best cell is > 5 km from every known family;
  * calibration: isotonic p_real fitted on the dev year, checked on the gate year against (a) the
    Bostock match rate and (b) an independent WITHHELD-STATION energy excess (B926/NLLB/LRIV, which
    never enter the back-projection), against matched random-time controls.
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

OUT = f"{P.RES}/p2"
FS = P.FS


# --------------------------------------------------------------------------- matching
def one2one(det_t, ev_t, tol=2.0):
    """greedy one-to-one, strongest-first is not available here so we go by time (Track A step 0b).

    Returns (matched index into ev for each detection, or -1)."""
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


def quiet_days(year, labels):
    q = labels[(labels.kind == "neg_certain_day") & (labels.year == year)]
    return sorted(q.doy.unique())


def ets_days(year, B):
    c = B[B.year == year].groupby("doy").size()
    return sorted(c[c >= 100].index)


# --------------------------------------------------------------------------- withheld-station check
def withheld_excess(det, stas, store, tt_cells, rng, e_win=6.0, e_pre=20.0, ncontrol=200):
    """log10 E(S-1..S+5) / E(same window 20 s earlier) in that station-day's MAD units, at the
    predicted S of the detection's cell, for stations that never entered the back-projection."""
    nE, nPre = int(e_win * FS), int(e_pre * FS)
    rows = []
    for (year, doy), g in det.groupby(["year", "doy"], sort=False):
        for s in stas:
            if not store.has(s, year, doy):
                continue
            i = store.idx(s, year, doy)
            n = int(store.nvalid[i])
            x = np.asarray(store.get_i(i)[:, :n], np.float32)
            e = (x ** 2).sum(0)
            c = np.concatenate([[0.0], np.cumsum(e.astype(np.float64))])
            W = c[nE:] - c[:-nE]
            v = np.full(len(W), np.nan)
            ok = slice(nPre, len(W))
            with np.errstate(divide="ignore", invalid="ignore"):
                v[ok] = np.log10(W[nPre:] / W[:len(W) - nPre])
            fin = np.isfinite(v)
            if fin.sum() < 1000:
                continue
            md = np.median(v[fin])
            sd = np.median(np.abs(v[fin] - md)) * 1.4826
            if sd <= 0:
                continue
            z = (v - md) / sd
            ts = float(store.t0[i])
            j = np.rint((g.t_origin.values + tt_cells[s][g.cell.values] - 1.0 - ts) * FS).astype(int)
            good = (j >= 0) & (j < len(z))
            zz = np.full(len(g), np.nan)
            zz[good] = z[j[good]]
            jc = rng.integers(nPre, len(z), ncontrol)
            rows.append(pd.DataFrame(dict(idx=g.index.values, sta=s, z=zz,
                                          ctl_p95=np.nanquantile(z[jc], 0.95),
                                          ctl_p99=np.nanquantile(z[jc], 0.99))))
    if not rows:
        return pd.DataFrame(columns=["idx", "sta", "z", "ctl_p95", "ctl_p99"])
    return pd.concat(rows, ignore_index=True)


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True, help="det_<tag>.parquet")
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--null", default="", help="tag of the moveout-shuffle run")
    ap.add_argument("--rates", default="6,21", help="network-wide quiet-day false events/day")
    ap.add_argument("--tol", type=float, default=2.0)
    ap.add_argument("--shift", type=float, default=307.0, help="chance-control time shift, s")
    ap.add_argument("--calib-from", default="", help="det tag of the DEV year for the isotonic fit")
    ap.add_argument("--calib-year", type=int, default=2012)
    ap.add_argument("--withheld", default="B926,NLLB,LRIV")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    rates = [float(x) for x in a.rates.split(",")]
    L = pd.read_parquet(f"{P.RES}/labels_v0.parquet")
    B = P.bostock()
    FAM = P.families()
    hid, hsta = P.holdouts()
    D = pd.read_parquet(f"{OUT}/det_{a.tag}.parquet")
    meta = json.load(open(f"{OUT}/det_{a.tag}_meta.json"))
    QD, ED = quiet_days(a.year, L), ets_days(a.year, B)
    print(f"[scope] {a.year}: {len(D)} peaks, {len(QD)} quiet days, {len(ED)} ETS days, "
          f"stations={meta['stations']}")

    # -------- operating points from QUIET days only
    q = D[D.doy.isin(QD)]
    thr = {}
    for r in rates:
        k = int(round(r * max(1, len(QD))))
        v = np.sort(q.z.values)[::-1]
        thr[r] = float(v[min(k, len(v)) - 1]) if len(v) else np.inf
    print(f"[fa] thresholds from {len(QD)} quiet days: " +
          "  ".join(f"{r}/day -> z={thr[r]:.2f}" for r in rates))

    res = dict(tag=a.tag, year=a.year, meta=meta, quiet_days=QD, ets_days=ED,
               thresholds={str(r): thr[r] for r in rates})
    Bx = B[B.year == a.year]
    ev_t = Bx.t_origin.to_numpy()
    tabs = []
    for r in rates:
        d = D[(D.z >= thr[r]) & D.doy.isin(ED)].reset_index(drop=True)
        m = one2one(d.t_origin.values, ev_t, tol=a.tol)
        ms = one2one(d.t_origin.values + a.shift, ev_t, tol=a.tol)
        ms2 = one2one(d.t_origin.values - a.shift, ev_t, tol=a.tol)
        nev = int((Bx.doy.isin(ED)).sum())
        rec = (m >= 0).sum() / nev
        chance = 0.5 * ((ms >= 0).sum() + (ms2 >= 0).sum()) / nev
        d["match"] = m
        d["bkey"] = np.where(m >= 0, Bx.bkey.to_numpy()[np.clip(m, 0, None)], -1)
        d["mag"] = np.where(m >= 0, Bx.mag.to_numpy()[np.clip(m, 0, None)], np.nan)
        d["hidden"] = np.where(m >= 0, Bx.hidden.to_numpy()[np.clip(m, 0, None)], False)
        tabs.append((r, d))
        # location error for matched detections
        ok = d.match >= 0
        if ok.any():
            fl = FAM.loc[d.bkey[ok]]
            derr = P.epi_km(d.lat[ok].values, d.lon[ok].values, fl.lat.values, fl.lon.values)
        else:
            derr = np.array([np.nan])
        # centroid collapse: distance from every detection's cell to the NEAREST known family
        dmin = np.full(len(d), np.nan)
        if len(d):
            fa = FAM[["lat", "lon"]].to_numpy()
            dd = np.stack([P.epi_km(d.lat.values, d.lon.values, fa[i, 0], fa[i, 1])
                           for i in range(len(fa))], 1)
            dmin = dd.min(1)
        res[f"rate_{r}"] = dict(
            n_det_ets=int(len(d)), n_events=nev, recall=float(rec), chance=float(chance),
            recall_chance_corrected=float((rec - chance) / max(1e-9, 1 - chance)),
            precision_vs_bostock=float((m >= 0).mean()) if len(d) else np.nan,
            loc_err_km_median=float(np.nanmedian(derr)), loc_err_km_p90=float(np.nanquantile(derr, .9)),
            dist_to_nearest_family_km_median=float(np.nanmedian(dmin)),
            frac_gt5km_from_every_family=float(np.nanmean(dmin > 5.0)),
            frac_gt10km_from_every_family=float(np.nanmean(dmin > 10.0)))
        print(f"[rate {r}/day] z>={thr[r]:.2f}  ETS detections {len(d)}  recall {rec:.4f} "
              f"(chance {chance:.4f} -> corrected {(rec-chance)/max(1e-9,1-chance):.4f})  "
              f"loc err {np.nanmedian(derr):.1f} km  >5 km from every family {np.nanmean(dmin>5):.3f}")

    # -------- strata at the first rate
    r0, d0 = tabs[0]
    Bs = Bx[Bx.doy.isin(ED)].reset_index(drop=True)
    mall = one2one(D[(D.z >= thr[r0]) & D.doy.isin(ED)].t_origin.values, Bs.t_origin.values, a.tol)
    Bs["found"] = False
    Bs.loc[mall[mall >= 0], "found"] = True
    tt_lab = P.fam_tt_labels()
    prov = {b: ("picked" if any(tt_lab.get((b, s), (0, "x"))[1] == "picked" for s in P.S15)
                else "predicted") for b in FAM.index}
    Bs["fam_prov"] = Bs.bkey.map(prov)
    Bs["magbin"] = pd.cut(Bs.mag, [-9, 1.0, 1.4, 1.8, 2.2, 9])
    strata = {
        "by_hidden": Bs.groupby("hidden").found.agg(["size", "mean"]).to_dict(),
        "by_mag": Bs.groupby("magbin", observed=True).found.agg(["size", "mean"]).to_dict(),
        "by_family_provenance": Bs.groupby("fam_prov").found.agg(["size", "mean"]).to_dict(),
    }
    Bs.groupby("bkey").found.agg(["size", "mean"]).reset_index().to_csv(
        f"{OUT}/score_{a.tag}_by_family.csv", index=False)
    nsit = D[(D.z >= thr[r0]) & D.doy.isin(ED)].groupby("n_sites").size().to_dict()
    res["strata"] = json.loads(pd.Series(strata).to_json())
    res["det_by_n_sites"] = {str(k): int(v) for k, v in nsit.items()}
    print("[strata] recall by hidden family:",
          Bs.groupby("hidden").found.mean().round(4).to_dict())
    print("[strata] recall by family label provenance:",
          Bs.groupby("fam_prov").found.mean().round(4).to_dict())

    # -------- moveout-shuffle null
    if a.null:
        N = pd.read_parquet(f"{OUT}/det_{a.null}.parquet")
        for r in rates:
            nq = N[N.doy.isin(QD)]
            v = np.sort(nq.z.values)[::-1]
            k = int(round(r * max(1, len(QD))))
            tn = float(v[min(k, len(v)) - 1]) if len(v) else np.inf
            n_true = int(((D.z >= thr[r]) & D.doy.isin(ED)).sum())
            n_null_same = int(((N.z >= thr[r]) & N.doy.isin(ED)).sum())
            n_null_own = int(((N.z >= tn) & N.doy.isin(ED)).sum())
            res[f"shuffle_{r}"] = dict(n_true=n_true, n_null_at_true_thr=n_null_same,
                                       ratio_same_threshold=n_null_same / max(1, n_true),
                                       null_thr_own_fa=tn, n_null_at_own_fa=n_null_own,
                                       ratio_matched_fa=n_null_own / max(1, n_true))
            print(f"[null {r}/day] true {n_true}  shuffled {n_null_same} "
                  f"(same z) -> ratio {n_null_same/max(1,n_true):.3f}; "
                  f"at its own quiet-day FA {n_null_own} -> {n_null_own/max(1,n_true):.3f}")

    # -------- calibration
    if a.calib_from:
        from sklearn.isotonic import IsotonicRegression
        Dc = pd.read_parquet(f"{OUT}/det_{a.calib_from}.parquet")
        Bc = B[B.year == a.calib_year]
        EDc = ets_days(a.calib_year, B)
        dc = Dc[Dc.doy.isin(EDc)]
        mc = one2one(dc.t_origin.values, Bc.t_origin.to_numpy(), a.tol)
        iso = IsotonicRegression(out_of_bounds="clip", y_min=0, y_max=1)
        iso.fit(dc.z.values, (mc >= 0).astype(float))
        dd = D[D.doy.isin(ED)].reset_index(drop=True)
        dd["p_real"] = iso.predict(dd.z.values)
        mm = one2one(dd.t_origin.values, ev_t, a.tol)
        dd["bostock_match"] = mm >= 0
        dd["bkey"] = np.where(mm >= 0, Bx.bkey.to_numpy()[np.clip(mm, 0, None)], -1)
        # withheld-station realness (stations never used in the back-projection)
        wst = [s for s in a.withheld.split(",") if s and s not in meta["stations"]]
        store = P.DayStore()
        z = np.load(f"{OUT}/grid_tt.npz")
        gs = [str(s) for s in z["stas"]]
        ttc = {s: z["tt"][:, gs.index(s)] for s in wst}
        top = dd[dd.z >= thr[min(rates)]].copy()
        W = withheld_excess(top, wst, store, ttc, np.random.default_rng(3))
        if len(W):
            agg = W.groupby("idx").apply(
                lambda g: pd.Series(dict(w_hit=float(np.nanmean(g.z > g.ctl_p95)),
                                         w_z=float(np.nanmean(g.z)))), include_groups=False)
            top = top.join(agg)
        top["bin"] = pd.cut(top.p_real, [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0])
        cal = top.groupby("bin", observed=True).agg(
            n=("p_real", "size"), p_mean=("p_real", "mean"),
            bostock=("bostock_match", "mean"),
            withheld=("w_hit", "mean") if "w_hit" in top else ("p_real", "size")).reset_index()
        cal["abs_dev_bostock"] = (cal.p_mean - cal.bostock).abs()
        cal.to_csv(f"{OUT}/calib_{a.tag}.csv", index=False)
        top.drop(columns=["bin"]).to_parquet(f"{OUT}/det_scored_{a.tag}.parquet",
                                             index=False)
        res["calibration"] = dict(
            max_abs_dev_vs_bostock=float(cal.abs_dev_bostock.max()),
            within_0p1=bool(cal.abs_dev_bostock.max() <= 0.1),
            withheld_stations=wst,
            withheld_hit_rate_overall=float(np.nanmean(top.w_hit)) if "w_hit" in top else None,
            withheld_control_rate=0.05)
        print("\n[calibration] (p_real fitted on " + a.calib_from + ")")
        print(cal.to_string(index=False))
        print(f"  max |p_real - Bostock match rate| = {cal.abs_dev_bostock.max():.3f} "
              f"-> within +-0.1: {cal.abs_dev_bostock.max()<=0.1}")

    out = a.out or f"{OUT}/score_{a.tag}.json"
    json.dump(res, open(out, "w"), indent=1, default=str)
    print(f"\n[done] -> {out}")


if __name__ == "__main__":
    main()
