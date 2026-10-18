#!/usr/bin/env python
"""E4 step 4 — score the tight-tolerance association re-null.

Observed vs shifted-pick null event rates at pick_match_tolerance 1.2 s, edt_pick_std 0.5 s,
zlim 20-50 km, n_picks run DIRECTLY at 4 and 5, thresholds 0.3 and 0.5; 5 constant-shift seeds and
2 hourly-shift seeds per configuration.  Plus Bostock recall by magnitude using the per-family
time-base offsets, and a cross-match of the clean catalog against any regional earthquake/blast
catalog found on disk.

Writes results/e4_assoc_null_table.csv, e4_assoc_bostock_recall.csv, e4_assoc_eqcheck.txt
"""
from __future__ import annotations
import glob, os, sys
import numpy as np, pandas as pd

RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
ROOT = "/home/jovyan/tremorferometry"


def rate(path):
    if not os.path.exists(path):
        return np.nan, 0
    e = pd.read_csv(path)
    if not len(e):
        return 0.0, 0
    t = e.time.to_numpy()
    span_d = (t.max() - t.min()) / 86400.0
    # use the same 34-day window span for every configuration
    return len(e) / 34.0, len(e)


def main():
    rows = []
    for thr in ("0.3", "0.5"):
        for n in (4, 5):
            obs, nobs = rate(f"{RES}/e4_assoc_vi2010_cascadia_t{thr}_n{n}.csv")
            for mode, seeds in (("station", [1, 2, 3, 4, 5]), ("station_hour", [1, 2])):
                v = []
                for s in seeds:
                    r, _ = rate(f"{RES}/e4_assocnull_vi2010_cascadia_t{thr}_n{n}_{mode}_s{s}.csv")
                    if np.isfinite(r):
                        v.append(r)
                if not v:
                    continue
                rows.append(dict(thr=float(thr), n_stations=n, shift=mode, obs_per_day=obs,
                                 n_obs=nobs, null_per_day=float(np.mean(v)),
                                 null_sd=float(np.std(v)), n_seeds=len(v),
                                 null_over_obs=float(np.mean(v) / obs) if obs else np.nan))
    T = pd.DataFrame(rows)
    T.to_csv(f"{RES}/e4_assoc_null_table.csv", index=False)
    print("=== tight-tolerance association null (tol 1.2 s, std 0.5 s, z 20-50 km) ===")
    print(T.round(4).to_string(index=False))
    # E3 comparison (wide tolerances), for the record
    old = pd.read_csv(f"{RES}/e3_nullevents_rates.csv") if \
        os.path.exists(f"{RES}/e3_nullevents_rates.csv") else None
    if old is not None:
        print("\n[E3 wide-tolerance reference]")
        print(old.round(3).to_string(index=False))

    # ---------------- Bostock recall by magnitude, per-family offsets
    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    fo = pd.read_csv(f"{RES}/e3_bostock_family_offsets.csv")
    acc = fo[fo.accept.astype(bool)]
    omap = dict(zip(acc.bkey.astype(int), acc.offset_s.astype(float)))
    bo = bo[bo.bkey.astype(int).isin(omap)].copy()
    bo["t"] = (bo.time.astype("datetime64[ns]").astype("int64") / 1e9
               + bo.bkey.astype(int).map(omap).to_numpy())
    MB = [(-9, 1.5, "M<1.5"), (1.5, 1.8, "M1.5-1.8"), (1.8, 2.1, "M1.8-2.1"), (2.1, 9, "M>=2.1")]
    rr = []
    for thr in ("0.3", "0.5"):
        for n in (4, 5):
            f = f"{RES}/e4_assoc_vi2010_cascadia_t{thr}_n{n}.csv"
            if not os.path.exists(f):
                continue
            e = pd.read_csv(f)
            if not len(e):
                continue
            et = np.sort(e.time.to_numpy())
            # chance floor: the same matching against 20 shifted copies of the event list
            rng = np.random.default_rng(3)
            for tol in (2.0, 5.0):
                for lo, hi, lab in MB + [(-9, 9, "all")]:
                    h = bo[(bo.mag > lo) & (bo.mag <= hi)]
                    if not len(h):
                        continue
                    idx = np.searchsorted(et, h.t.to_numpy())
                    d = np.full(len(h), np.inf)
                    for k in (-1, 0):
                        j = np.clip(idx + k, 0, len(et) - 1)
                        d = np.minimum(d, np.abs(et[j] - h.t.to_numpy()))
                    rec = float((d <= tol).mean())
                    ch = []
                    for _ in range(10):
                        sh = np.sort(et + rng.uniform(300, 3000) * rng.choice([-1, 1]))
                        i2 = np.searchsorted(sh, h.t.to_numpy())
                        d2 = np.full(len(h), np.inf)
                        for k in (-1, 0):
                            j = np.clip(i2 + k, 0, len(sh) - 1)
                            d2 = np.minimum(d2, np.abs(sh[j] - h.t.to_numpy()))
                        ch.append(float((d2 <= tol).mean()))
                    rr.append(dict(thr=float(thr), n_stations=n, tol_s=tol, magbin=lab,
                                   n_bostock=len(h), recall=rec, chance=float(np.mean(ch)),
                                   recall_excess=rec - float(np.mean(ch)),
                                   n_events=len(e), ev_per_day=len(e) / 34.0))
    RR = pd.DataFrame(rr)
    RR.to_csv(f"{RES}/e4_assoc_bostock_recall.csv", index=False)
    print("\n=== Bostock recall by magnitude (per-family offsets, +-2 s) ===")
    print(RR[RR.tol_s == 2.0].round(4).to_string(index=False))

    # ---------------- cross-match the clean catalog against the regional earthquake catalogs
    xrows = []
    for cat in ("eq_catalog_pgc", "eq_catalog_b011"):
        f = f"{ROOT}/data/{cat}.csv"
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        te = pd.to_datetime(d.ot, format="mixed", utc=True).dt.tz_localize(None) \
               .astype("datetime64[ns]").astype("int64").to_numpy() / 1e9
        for thr in ("0.3", "0.5"):
            for n in (4, 5):
                fp = f"{RES}/e4_assoc_vi2010_cascadia_t{thr}_n{n}.csv"
                if not os.path.exists(fp):
                    continue
                e = pd.read_csv(fp)
                if not len(e):
                    continue
                et = np.sort(e.time.to_numpy())
                tw = np.sort(te[(te >= et.min() - 120) & (te <= et.max() + 120)])
                if not len(tw):
                    xrows.append(dict(catalog=cat, thr=float(thr), n_stations=n, n_eq=0,
                                      n_events=len(e), n_matched=0, frac=0.0))
                    continue
                i = np.searchsorted(tw, et)
                dmin = np.full(len(et), np.inf)
                for k in (-1, 0):
                    j = np.clip(i + k, 0, len(tw) - 1)
                    dmin = np.minimum(dmin, np.abs(tw[j] - et))
                # an LFE "event" within 60 s of a catalogued regular earthquake in the region is a
                # candidate contamination (the earthquake's own P/S picks can associate)
                xrows.append(dict(catalog=cat, thr=float(thr), n_stations=n, n_eq=int(len(tw)),
                                  n_events=len(e), n_matched=int((dmin <= 60).sum()),
                                  frac=float((dmin <= 60).mean())))
    X = pd.DataFrame(xrows)
    X.to_csv(f"{RES}/e4_assoc_eq_crossmatch.csv", index=False)
    print("\n=== cross-match vs regional earthquake catalogs (|dt| <= 60 s) ===")
    print(X.round(4).to_string(index=False))

    # ---------------- regional earthquake / blast catalogs on disk?
    pats = ["**/*comcat*", "**/*anss*", "**/*ANSS*", "**/*usgs*", "**/*earthquake*",
            "**/*eq_catalog*", "**/*pnsn*event*", "**/*nrcan*"]
    found = []
    for base in (f"{ROOT}/catalogs", f"{ROOT}/data", f"{RES}"):
        for p in pats:
            found += [f for f in glob.glob(f"{base}/{p}", recursive=True)
                      if os.path.isfile(f) and os.path.getsize(f) > 1000]
    found = sorted(set(found))
    with open(f"{RES}/e4_assoc_eqcheck.txt", "w") as fh:
        if not found:
            msg = ("E4 step 4: NO regional earthquake/blast catalog (ANSS/ComCat/NRCan event list) "
                   "was found on disk under catalogs/, data/ or results/ -- searched for "
                   "comcat/anss/usgs/earthquake/eq_catalog/nrcan patterns. The clean catalog was "
                   "therefore NOT screened against regular earthquakes or quarry blasts. Fetching "
                   "one is a small metadata download (FDSN event service), not a waveform pull.")
        else:
            msg = "E4 step 4: candidate earthquake/blast catalogs on disk:\n" + \
                  "\n".join(f"  {f} ({os.path.getsize(f)/1e6:.1f} MB)" for f in found)
        fh.write(msg + "\n")
    print("\n" + msg)


if __name__ == "__main__":
    main()
