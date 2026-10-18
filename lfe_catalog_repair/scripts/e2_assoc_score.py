#!/usr/bin/env python
"""E2 step 3a-score (added 2026-09-10, resumed run) — score PyOcto associated events.

For each associated-event file:
  * events/day, station-count histogram, depth distribution;
  * ETS vs quiet event rate (PNSN tremor inside the association box, 30-d rolling, >p90 / <=p25 --
    the D1 convention);
  * recall / precision against Lin et al. (2024) ORIGIN TIMES and (where the window allows) the
    Bostock et al. (2015) catalogue times, with
      - a fixed-tolerance count next to a shifted-time null (+-600 s, +-1800 s, one random shift),
      - a chance-free PEAK-EXCESS estimate: the histogram of (associated time - catalogue time) is
        a narrow peak on a flat chance background; the excess area over +-1 s around the measured
        peak divided by n_gt (recall) or n_events (precision).
  * Bostock time convention: the Slow-Earthquake-DB record is NOT an origin time on our time base;
    its offset is measured at PGC in e2_pgc2010_offsets.csv (+12.4 s w.r.t. t + C + ttS) and applied
    here as an equivalent origin time t + C + 12.4 s.  The offset is distance-dependent at the +-5 s
    level (see NOTES), so the Bostock numbers are lower bounds.
  * spatial comparison of matched events with the Lin epicentre (Lin depths are unreliable).
Writes results/e2_assoc_score_<tag>.csv (one row per gt x tolerance) and appends the summary row to
results/e2_assoc_scores.csv.
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_resid_hist as RH  # noqa: E402

TOLS = (2.0, 5.0)
NULL_SHIFTS = (-600.0, 600.0, -1800.0, 1800.0)
BOX = dict(lat=(47.7, 49.6), lon=(-125.4, -122.3))
R_KM_DEG = 111.195


def epoch(s):
    return pd.to_datetime(s).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9


def match_frac(a, b, tol):
    if len(a) == 0 or len(b) == 0:
        return 0.0, 0.0
    a = np.sort(a); b = np.sort(b)

    def frac(x, y):
        i = np.searchsorted(y, x)
        d = np.full(len(x), np.inf)
        for off in (-1, 0):
            j = np.clip(i + off, 0, len(y) - 1)
            d = np.minimum(d, np.abs(y[j] - x))
        return float((d <= tol).mean())
    return frac(a, b), frac(b, a)


def peak_excess(gt, ev, maxlag=40.0, binw=0.25, half=1.0, half2=2.0):
    d = RH.pairs_within(np.sort(gt), np.sort(ev), maxlag)
    e = np.arange(-maxlag, maxlag + binw / 2, binw)
    h, _ = np.histogram(d, bins=e)
    c = 0.5 * (e[1:] + e[:-1])
    k = int(np.argmax(h))
    bgm = np.abs(c - c[k]) > 5.0
    bg = float(np.median(h[bgm]))
    sbg = float(h[bgm].std(ddof=1))
    out = {}
    for lab, w in (("", half), ("2s_", half2)):
        m = np.abs(c - c[k]) <= w
        exc = float(h[m].sum() - bg * m.sum())
        sig = float(np.sqrt(max(h[m].sum(), 1) + m.sum() ** 2 * sbg ** 2 / max(bgm.sum(), 1)))
        out[f"{lab}excess"] = exc
        out[f"{lab}excess_snr"] = exc / max(sig, 1e-9)
    out.update(peak_lag_s=float(c[k]), peak_over_bg=h[k] / max(bg, 1e-9), bg_per_bin=bg)
    return out


def bostock_offset():
    """equivalent-origin-time offset for the Slow-Earthquake-DB Bostock times (see docstring).

    E3 2026-09-11: the `C` term came from `e2_bostock_lin_dt.csv`, whose estimator is
    chance-dominated (Merlin F7 / E2 bug 3) and which e2_groundtruth.py no longer writes; it was
    +0.086 s, i.e. numerically irrelevant.  If the file is absent C is taken as 0 and the offset is
    the measured PGC peak alone.  PREFER the per-family table (`--bostock-family-offsets`).
    """
    f = f"{RES}/e2_bostock_lin_dt.csv"
    C = float(np.median(pd.read_csv(f).dt)) if os.path.exists(f) else 0.0
    o = pd.read_csv(f"{RES}/e2_pgc2010_offsets.csv").set_index("gt")
    return C + float(o.loc["bostock_stack", "peak_lag_s"]), C


def bostock_times(bo, famoff=None):
    """Equivalent ORIGIN times for the Bostock table.

    famoff = path to `e3_bostock_family_offsets.csv` (E3 step 5).  With it, each family gets its
    own measured offset (accepted families only; the rest fall back to the median accepted
    offset).  Without it, the single E2 constant is used -- the old behaviour.
    Returns (times, label, n_family_specific).
    """
    t0 = epoch(bo.time)
    if not famoff:
        Cb, _ = bostock_offset()
        return t0 + Cb, f"single constant {Cb:+.2f} s", 0
    fo = pd.read_csv(famoff)
    acc = fo[fo.accept.astype(bool)]
    med = float(acc.offset_s.median())
    m = dict(zip(acc.bkey.astype(int), acc.offset_s.astype(float)))
    off = bo.bkey.astype(int).map(m)
    n_spec = int(off.notna().sum())
    off = off.fillna(med).to_numpy()
    return t0 + off, (f"per-family (n={len(m)} families, median {med:+.2f} s, "
                      f"sd {acc.offset_s.std():.2f}; {n_spec}/{len(bo)} events family-specific)"), n_spec


def tremor_daily(box):
    tr = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time).astype("datetime64[ns]")
    tr = tr[(tr.lat.between(*box["lat"])) & (tr.lon.between(*box["lon"]))]
    return tr.groupby(tr.time.dt.date).size().rename("tremor")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True)
    ap.add_argument("--min-stations", type=int, default=0)
    ap.add_argument("--bostock-family-offsets", default="",
                    help="E3 step 5: path to e3_bostock_family_offsets.csv -> use the measured "
                         "PER-FAMILY Bostock time offsets instead of the single E2 constant")
    ap.add_argument("--tag-suffix", default="",
                    help="appended to the output tag so an E3 re-score cannot overwrite E2 files")
    args = ap.parse_args()
    tag = os.path.basename(args.events)[len("e2_events_"):-4]
    m = re.match(r"(vi\d+)_(\w+?)_t([\d.]+)", tag)
    setname, model, thr = m.group(1), m.group(2), float(m.group(3))

    ev = pd.read_csv(args.events)
    if args.min_stations:
        ev = ev[ev["n_stations"] >= args.min_stations]
        tag += f"_n{args.min_stations}"
    tag += args.tag_suffix
    et = ev.time.to_numpy()                          # PyOcto origin time, epoch s
    t0, t1 = et.min(), et.max()
    days = (t1 - t0) / 86400.0
    print(f"[{tag}] {len(ev)} events over {days:.1f} d = {len(ev)/days:.1f} events/day; "
          f"stations {ev.n_stations.value_counts().sort_index().to_dict()}")
    print(f"[{tag}] depth km: median {ev.depth.median():.1f}, "
          f"frac<5 km {float((ev.depth<5).mean()):.3f}, frac 25-45 km "
          f"{float(ev.depth.between(25,45).mean()):.3f}; rms resid median {ev.rms_resid.median():.2f} s")

    # ---- ground truth
    lin = pd.read_parquet(f"{RES}/e2_gt_lin.parquet")
    # NOT filtered on the `window` column: e2_groundtruth.py de-duplicates across windows, so the
    # vi2010 rows were dropped as duplicates of the (wider) pgc2010 rows and window=="vi2010" is
    # EMPTY.  Selecting on the event time itself is equivalent and window-name independent.
    lin = lin.copy()
    lin["t"] = epoch(lin.OT)
    lin = lin[(lin.t >= t0 - 60) & (lin.t <= t1 + 60)]
    gts = {"lin_all": lin,
           "lin_hq": lin[(lin.residual < 0.5) & (lin.N >= 4)],
           "lin_N4": lin[lin.N >= 4]}
    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    bo["t"], blab, _ = bostock_times(bo, args.bostock_family_offsets)
    bo = bo[(bo.t >= t0 - 60) & (bo.t <= t1 + 60)]
    if len(bo):
        gts["bostock"] = bo
        print(f"[gt] Bostock equivalent-origin offset: {blab}")
    print("[gt] " + ", ".join(f"{k}={len(v)}" for k, v in gts.items()))

    rng = np.random.default_rng(0)
    rows = []
    for gk, g in gts.items():
        if len(g) < 50:      # e.g. the vi2017 window: the Lin catalogue ends 2017-02-21, so only
            continue         # 4-7 Lin events overlap -> no usable ground truth there
        gt = g.t.to_numpy()
        pe = peak_excess(gt, et)
        for tol in TOLS:
            rec, prec = match_frac(gt, et, tol)
            nr, np_ = [], []
            for sh in NULL_SHIFTS:
                r, p = match_frac(gt + sh, et, tol)
                nr.append(r); np_.append(p)
            r, p = match_frac(gt + rng.uniform(3000, 40000), et, tol)
            nr.append(r); np_.append(p)
            rows.append(dict(set=setname, model=model, threshold=thr, tag=tag, gt=gk,
                             tol_s=tol, n_gt=len(gt), n_events=len(ev),
                             events_per_day=len(ev) / days,
                             recall=rec, recall_null=float(np.mean(nr)),
                             precision=prec, precision_null=float(np.mean(np_)),
                             recall_excess=pe["excess"] / len(gt),
                             precision_excess=pe["excess"] / len(ev),
                             recall_excess_2s=pe["2s_excess"] / len(gt),
                             precision_excess_2s=pe["2s_excess"] / len(ev),
                             **pe))
    sc = pd.DataFrame(rows)
    sc.to_csv(f"{RES}/e2_assoc_score_{tag}.csv", index=False)
    if len(sc):
        print(sc[["gt", "tol_s", "n_gt", "recall", "recall_null", "precision", "precision_null",
              "peak_lag_s", "peak_over_bg", "excess_snr", "recall_excess", "precision_excess",
                  "recall_excess_2s", "precision_excess_2s"]].round(4).to_string(index=False))
    else:
        print("[score] no ground truth with >= 50 events in this window -- rates only")

    # ---- ETS vs quiet
    trd = tremor_daily(BOX)
    ev["date"] = pd.to_datetime(ev.time_utc).dt.date
    c = ev.groupby("date").size().rename("events")
    df = pd.concat([c, trd], axis=1)
    df = df[(df.index >= pd.to_datetime(t0, unit="s").date())
            & (df.index <= pd.to_datetime(t1, unit="s").date())].fillna(0.0).sort_index()
    roll = df.tremor.rolling(30, min_periods=5, center=True).sum()
    hi, lo = np.nanpercentile(roll, 90), np.nanpercentile(roll, 25)
    from scipy.stats import spearmanr
    rho, pv = spearmanr(df.events, df.tremor)
    ets, qui = df.events[roll >= hi].mean(), df.events[roll <= lo].mean()
    print(f"[rate] ETS {ets:.1f} ev/d, quiet {qui:.1f} ev/d, ratio {ets/max(qui,1e-9):.2f}; "
          f"spearman(events, daily tremor in box) {rho:+.3f} (p={pv:.1e}, n={len(df)} d)")
    df.to_csv(f"{RES}/e2_assoc_daily_{tag}.csv")

    # ---- location comparison for time-matched Lin events (epicentre only; Lin depths unreliable)
    if "lin_hq" in gts and len(gts["lin_hq"]):
        g = gts["lin_hq"]
        gt = g.t.to_numpy()
        i = np.searchsorted(np.sort(et), gt)
        se = np.sort(et)
        order = np.argsort(et)
        best = np.full(len(gt), np.nan, dtype=int)
        bd = np.full(len(gt), np.inf)
        for off in (-1, 0):
            j = np.clip(i + off, 0, len(se) - 1)
            d = np.abs(se[j] - gt)
            m = d < bd
            bd, best = np.where(m, d, bd), np.where(m, order[j], best)
        ok = bd <= 2.0
        if ok.sum() >= 10:
            e_ok = ev.iloc[best[ok]]
            dx = (e_ok.longitude.values - g.lon.values[ok]) * R_KM_DEG * np.cos(np.radians(48.6))
            dy = (e_ok.latitude.values - g.lat.values[ok]) * R_KM_DEG
            dh = np.hypot(dx, dy)
            dz = e_ok.depth.values - np.abs(g.depth.values[ok])
            print(f"[loc vs Lin] {int(ok.sum())} matched (<=2 s): horizontal offset median "
                  f"{np.median(dh):.1f} km, IQR {np.percentile(dh,25):.1f}-{np.percentile(dh,75):.1f}; "
                  f"depth diff median {np.median(dz):+.1f} km, IQR {np.percentile(dz,25):+.1f}.."
                  f"{np.percentile(dz,75):+.1f}")
            pd.DataFrame(dict(dt=bd[ok], dh_km=dh, dz_km=dz,
                              lat_ev=e_ok.latitude.values, lon_ev=e_ok.longitude.values,
                              dep_ev=e_ok.depth.values, lat_lin=g.lat.values[ok],
                              lon_lin=g.lon.values[ok], dep_lin=g.depth.values[ok],
                              n_sta=e_ok.n_stations.values)
                         ).to_csv(f"{RES}/e2_assoc_locmatch_{tag}.csv", index=False)
    # E3 2026-09-11: an E3 re-score must not append to the E2 master file (it did once, 24 rows,
    # and they had to be split back out).  Any run with an E3 option writes its own master.
    allf = (f"{RES}/e3_assoc_scores_famoff.csv"
            if (args.bostock_family_offsets or args.tag_suffix)
            else f"{RES}/e2_assoc_scores.csv")
    sc.to_csv(allf, mode="a", header=not os.path.exists(allf), index=False)
    print(f"[write] {RES}/e2_assoc_score_{tag}.csv (+ appended to {allf})")


if __name__ == "__main__":
    main()
