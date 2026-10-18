#!/usr/bin/env python
"""E3 step 6a — the off-ground-truth attribution, recomputed (Merlin F6).

E2 read "WISH and YACT are flat through ETS while the borehole B943 responds" as evidence about
GENERALISATION and about surface-vs-borehole instruments.  It is neither: it is distance to the
tremor.  Three measurements settle it, all from data on disk.

  1. How much tremor is actually near each off-GT station: the number of days in 2017-2018 with any
     PNSN tremor within 40 / 60 / 80 / 100 km.
  2. Whether SURFACE stations respond when they ARE near the source: the daily S-pick rate vs
     tremor for the eight vi2017 stations, ETS/quiet defined on the FULL PNSN record within 40 km
     of each station (the F10 fix), so surface CN stations inside the 2017 Vancouver Island source
     zone can be compared with the co-located borehole B011.
  3. The direct surface/borehole control: CN.PGC and PB.B011 are 0.6 km apart; their
     single-station recall against Lin-hq in the vi2010 window, by the chance-free +-1 s
     peak-excess estimator.

Writes results/e3_offgt_attribution.csv, results/e3_surface_vs_borehole.csv
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
TREMOR = "catalogs/pnsn_tremor_cascadia_full.csv"
R_KM_DEG = 111.195
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_resid_hist as RH  # noqa: E402
import e3_null_score as NS  # noqa: E402


def epoch(s):
    return pd.to_datetime(s).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9


def main():
    st = pd.read_csv(f"{RES}/e2_stations.csv").set_index("sta")
    tr = pd.read_csv(TREMOR, usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time).astype("datetime64[ns]")

    # ---- 1. tremor days near each station, 2017-2018
    t2 = tr[(tr.time >= "2017-01-01") & (tr.time < "2019-01-01")]
    rows = []
    for s in ["WISH", "YACT", "B943", "PGC", "B011", "SYMB", "SNB", "VGZ", "CLRS", "PTRF",
              "LRIV"]:
        if s not in st.index:
            continue
        la, lo = float(st.loc[s, "lat"]), float(st.loc[s, "lon"])
        d = np.hypot((t2.lon.values - lo) * R_KM_DEG * np.cos(np.radians(la)),
                     (t2.lat.values - la) * R_KM_DEG)
        r = dict(sta=s, lat=la, lon=lo)
        for R in (40, 60, 80, 100):
            r[f"tremor_days_{R}km"] = int(t2.time[d <= R].dt.normalize().nunique())
            r[f"tremor_count_{R}km"] = int((d <= R).sum())
        rows.append(r)
    A = pd.DataFrame(rows)
    A.to_csv(f"{RES}/e3_offgt_attribution.csv", index=False)
    print("[1] PNSN tremor days in 2017-2018 within R km of each station")
    print(A[["sta"] + [f"tremor_days_{R}km" for R in (40, 60, 80, 100)]].to_string(index=False))

    # ---- 2. per-station ETS/quiet inside the 2017 VI source zone, FULL-record definition
    from scipy.stats import spearmanr
    rows = []
    for setname, pf in (("vi2017", f"{RES}/e2_picks_vi2017_cascadia.parquet"),
                        ("offgt", None)):
        files = ([pf] if pf else
                 [f"{RES}/e2_picks_offgt_{s}_cascadia.parquet" for s in ("WISH", "B943", "YACT")])
        for f in files:
            if not os.path.exists(f):
                continue
            pk = pd.read_parquet(f)
            pk = pk[pk.phase == "S"]
            for s, g in pk.groupby("sta"):
                if s not in st.index:
                    continue
                lab, cnt, roll, hi, lo = NS.ets_quiet_labels(
                    float(st.loc[s, "lat"]), float(st.loc[s, "lon"]), radius_km=40.0)
                for thr in (0.1, 0.3, 0.5):
                    sel = g[(g.thr_first <= thr) & (g.prob >= thr)]
                    c = sel.groupby(pd.to_datetime(sel.day)).size()
                    days = pd.date_range(pd.to_datetime(g.day).min(),
                                         pd.to_datetime(g.day).max(), freq="D")
                    c = c.reindex(days, fill_value=0)
                    L = lab.reindex(days).fillna("na")
                    trd = cnt.reindex(days).fillna(0)
                    ne, nq = int((L == "ets").sum()), int((L == "quiet").sum())
                    e = float(c[L == "ets"].mean()) if ne else np.nan
                    q = float(c[L == "quiet"].mean()) if nq else np.nan
                    rho, pv = spearmanr(c.values, trd.values)
                    rows.append(dict(set=setname, sta=s, threshold=thr, n_days=len(days),
                                     n_ets_days=ne, n_quiet_days=nq,
                                     picks_per_day=float(c.mean()),
                                     ets_per_day=e, quiet_per_day=q,
                                     ets_quiet=e / q if q else np.nan,
                                     spearman=rho, p=pv,
                                     tremor_days_40km=int((trd > 0).sum())))
    E = pd.DataFrame(rows)
    E.to_csv(f"{RES}/e3_offgt_ets_quiet.csv", index=False)
    print("\n[2] daily S-pick rate vs tremor within 40 km of EACH station "
          "(ETS/quiet from the FULL PNSN record, cascadia thr 0.3)")
    print(E[E.threshold == 0.3][["set", "sta", "n_days", "tremor_days_40km", "n_ets_days",
                                 "n_quiet_days", "picks_per_day", "ets_quiet", "spearman"]]
          .round(3).to_string(index=False))

    # ---- 3. PGC vs B011, single-station recall against Lin-hq (vi2010)
    lin = pd.read_parquet(f"{RES}/e2_gt_lin.parquet")
    linq = lin[(lin.residual < 0.5) & (lin.N >= 4)]
    rows = []
    for model in ("cascadia", "jcms"):
        pk = pd.read_parquet(f"{RES}/e2_picks_vi2010_{model}.parquet")
        pk = pk[pk.phase == "S"]
        for s in ("PGC", "B011", "SNB", "VGZ", "LZB", "PFB", "YOUB", "NLLB", "B926", "LRIV"):
            g = pk[pk.sta == s]
            if len(g) == 0:
                continue
            t0, t1 = epoch(g.time).min(), epoch(g.time).max()
            gt = linq[(epoch(linq.OT) >= t0 - 60) & (epoch(linq.OT) <= t1 + 60)]
            pred = epoch(gt.OT) + gt[f"ttS_{s}"].to_numpy()
            for thr in (0.1, 0.3):
                sel = np.sort(epoch(g.time[(g.thr_first <= thr) & (g.prob >= thr)]))
                if len(sel) == 0:
                    continue
                d = RH.pairs_within(np.sort(pred), sel, 30.0)
                e = np.arange(-30, 30.001, 0.25)
                h, _ = np.histogram(d, bins=e)
                c = 0.5 * (e[1:] + e[:-1])
                bg = float(np.median(h[np.abs(c) > 10]))
                m = np.abs(c) <= 1.0
                exc = float(h[m].sum() - bg * m.sum())
                rows.append(dict(model=model, sta=s, threshold=thr,
                                 is_borehole=s.startswith("B") and s != "B926" or s in
                                 ("B011", "B926"),
                                 n_gt=len(gt), n_picks=len(sel),
                                 picks_per_day=len(sel) / ((t1 - t0) / 86400),
                                 excess=exc, recall_excess=exc / len(gt),
                                 precision_excess=exc / len(sel),
                                 sigma=exc / np.sqrt(max(h[m].sum(), 1))))
    B = pd.DataFrame(rows)
    B.to_csv(f"{RES}/e3_surface_vs_borehole.csv", index=False)
    print("\n[3] single-station recall vs Lin-hq (+-1 s peak excess), vi2010")
    print(B[B.model == "cascadia"][["sta", "threshold", "is_borehole", "picks_per_day",
                                    "recall_excess", "precision_excess", "sigma"]]
          .round(4).to_string(index=False))
    for thr in (0.1, 0.3):
        p = B[(B.model == "cascadia") & (B.threshold == thr) & (B.sta == "PGC")].recall_excess
        b = B[(B.model == "cascadia") & (B.threshold == thr) & (B.sta == "B011")].recall_excess
        if len(p) and len(b):
            print(f"  PGC (surface) {float(p.iloc[0]):.3f} vs B011 (borehole, 0.6 km away) "
                  f"{float(b.iloc[0]):.3f} at thr {thr} -> ratio "
                  f"{float(b.iloc[0])/max(float(p.iloc[0]),1e-9):.2f}x")
    print(f"\n[write] {RES}/e3_offgt_attribution.csv, e3_offgt_ets_quiet.csv, "
          f"e3_surface_vs_borehole.csv")


if __name__ == "__main__":
    main()
