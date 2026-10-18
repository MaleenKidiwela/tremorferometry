#!/usr/bin/env python
"""E2 step 4 — off-ground-truth behaviour of the picker at WISH, B943 and YACT, 2017-2018.

No LFE catalogue exists at these stations, so the only available references are
  (i)  PNSN tremor within 40 km (ETS vs quiet daily pick rate, Spearman correlation),
  (ii) WISH x B943 same-source coincidence: the two stations are 104 km apart and share 160
       discovery source cells, so a real LFE in a shared cell must appear at both with a
       differential S time set by the geometry.  We cross-correlate the two S-pick series in 0.5 s
       bins over +-40 s of lag and ask whether the peak sits inside the range of differential
       travel times predicted by the Crosson 1-D model for the shared cells, and how big it is
       relative to the +-600 s shifted background.
  (iii) YACT, the station D1 measured as 0 % responsive, as the negative control.

Writes `results/e2_offgt_summary.csv` and `results/e2_offgt_daily.csv`.
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d3_vel1d as v1d  # noqa: E402

RES = "lfe_catalog_repair/results"
STA = {"WISH": (47.11596, -123.77143), "B943": (47.813202, -122.911301),
       "YACT": (45.9325, -122.4193)}
DISC = {"WISH": "data/wish_disc.summary.csv",
        "B943": "data/b943_disc_p70_2010_2026_m3.summary.csv",
        "YACT": "data/yact_disc.summary.csv"}
THRS = (0.1, 0.3, 0.5)


def load_picks():
    fr = []
    for p in sorted(glob.glob(f"{RES}/e2_picks_offgt_*.parquet")):
        fr.append(pd.read_parquet(p))
    df = pd.concat(fr, ignore_index=True)
    df["t"] = pd.to_datetime(df.time).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9
    df["sec_of_day"] = df.t % 86400.0
    return df


# DAY-BOUNDARY ARTIFACT (found 2026-09-10, resumed run).  The `cascadia` weights fire on the first
# window of a day file: 0.18-0.67 picks per station-day land in the first 5 s of the day, 5-10x the
# rate expected for a 60 s slice, and 86 % of the WISH x B943 "zero-lag coincidences" were exactly
# those.  They are a per-day edge effect shared by every station cut on the same UTC boundary, so
# they manufacture a one-bin spike at lag 0 (5.9x background at thr 0.3) where the geometry of the
# shared source cells predicts -17.4 .. +0.6 s.  The coincidence test therefore drops the first and
# last 60 s of every day file.  Daily rates are unaffected (the artifact is 0.3-0.7 % of picks).
EDGE_S = 60.0


def drop_day_edges(df):
    return df[(df.sec_of_day >= EDGE_S) & (df.sec_of_day <= 86400.0 - EDGE_S)]


def tremor_daily(lat, lon, radius_km=40.0):
    tr = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time).astype("datetime64[ns]")
    d = np.hypot((tr.lat - lat) * 111.195, (tr.lon - lon) * 111.195 * np.cos(np.radians(lat)))
    tr = tr[(d < radius_km) & (tr.time >= "2017-01-01") & (tr.time < "2019-01-01")]
    return tr.groupby(tr.time.dt.date).size().rename("tremor")


def shared_cells():
    """0.05-deg source cells present in both the WISH and the B943 discovery family lists."""
    def cells(p):
        s = pd.read_csv(p)
        return set(zip(np.round(s.lat / 0.05) * 0.05, np.round(s.lon / 0.05) * 0.05))
    a, b = cells(DISC["WISH"]), cells(DISC["B943"])
    sh = sorted(a & b)
    return pd.DataFrame(sh, columns=["lat", "lon"])


def diff_tt_range(cells):
    (la, lo), (lb, lob) = STA["WISH"], STA["B943"]
    z = v1d.slab2_depth(cells.lat.values, cells.lon.values)
    def epi(lat, lon, s):
        return np.hypot((cells.lat - lat) * 111.195,
                        (cells.lon - lon) * 111.195 * np.cos(np.radians(lat)))
    import e2_groundtruth as gtm
    dA = np.hypot((cells.lat - la) * 111.195, (cells.lon - lo) * 111.195 * np.cos(np.radians(la)))
    dB = np.hypot((cells.lat - lb) * 111.195, (cells.lon - lob) * 111.195 * np.cos(np.radians(lb)))
    tA = gtm.tt_vec(dA.values, z, "S")
    tB = gtm.tt_vec(dB.values, z, "S")
    return tB - tA, z


def xcorr_picks(tA, tB, maxlag=40.0, dt=0.5):
    """binned cross-correlation of two pick-time series; returns lags and counts."""
    if len(tA) == 0 or len(tB) == 0:
        return np.array([]), np.array([])
    t0 = min(tA.min(), tB.min())
    nb = int((max(tA.max(), tB.max()) - t0) / dt) + 2
    a = np.bincount(((tA - t0) / dt).astype(np.int64), minlength=nb).astype(np.float32)
    b = np.bincount(((tB - t0) / dt).astype(np.int64), minlength=nb).astype(np.float32)
    nl = int(maxlag / dt)
    lags = np.arange(-nl, nl + 1) * dt
    out = np.empty(len(lags))
    for i, k in enumerate(range(-nl, nl + 1)):
        if k >= 0:
            out[i] = float(np.dot(a[:nb - k], b[k:]))
        else:
            out[i] = float(np.dot(a[-k:], b[:nb + k]))
    return lags, out


def main():
    pk = load_picks()
    print(f"[picks] {len(pk)} rows, models {sorted(pk.model.unique())}, "
          f"stations {sorted(pk.sta.unique())}")
    from scipy.stats import spearmanr
    rows, dailies = [], []
    for sta in ("WISH", "B943", "YACT"):
        trd = tremor_daily(*STA[sta])
        for model in sorted(pk.model.unique()):
            for thr in THRS:
                s = pk[(pk.sta == sta) & (pk.model == model) & (pk.phase == "S")
                       & (pk.thr_first <= thr) & (pk.prob >= thr)]
                if len(s) == 0:
                    continue
                c = s.groupby(s.day).size().rename("picks")
                df = pd.concat([c, trd], axis=1)
                df["picks"] = df.picks.fillna(0.0)
                df["tremor"] = df.tremor.fillna(0.0)
                df = df[df.picks > 0]  # only days actually processed
                df = df.sort_index()
                roll = df.tremor.rolling(30, min_periods=10, center=True).sum()
                hi, lo = np.nanpercentile(roll, 90), np.nanpercentile(roll, 25)
                ets = df.picks[roll >= hi].mean()
                qui = df.picks[roll <= lo].mean()  # BUGFIX 2026-09-10: same empty-quiet-set bug
                rho, pv = spearmanr(df.picks, df.tremor)
                edge = float(((s.sec_of_day < 60) | (s.sec_of_day > 86340)).mean())
                rows.append(dict(station=sta, model=model, threshold=thr, n_days=len(df),
                                 frac_in_day_edge_60s=edge,
                                 edge_over_expected=edge / (120.0 / 86400.0),
                                 picks_per_day=df.picks.mean(),
                                 ets_picks_per_day=ets, quiet_picks_per_day=qui,
                                 ets_quiet_ratio=ets / max(qui, 1e-9),
                                 spearman_rho=rho, spearman_p=pv,
                                 tremor_days=int((df.tremor > 0).sum())))
                d = df.reset_index().rename(columns={"index": "date", "day": "date"})
                d["station"], d["model"], d["threshold"] = sta, model, thr
                dailies.append(d)
    summ = pd.DataFrame(rows)
    summ.to_csv(f"{RES}/e2_offgt_summary.csv", index=False)
    pd.concat(dailies, ignore_index=True).to_csv(f"{RES}/e2_offgt_daily.csv", index=False)
    print(summ.round(3).to_string(index=False))

    # ---- WISH x B943 same-source coincidence
    cells = shared_cells()
    dtt, z = diff_tt_range(cells)
    print(f"\n[shared] {len(cells)} shared 0.05-deg source cells, Slab2 depth "
          f"{np.nanmin(z):.1f}-{np.nanmax(z):.1f} km; predicted S(B943)-S(WISH) "
          f"{np.nanpercentile(dtt,5):+.1f} .. {np.nanpercentile(dtt,95):+.1f} s "
          f"(median {np.nanmedian(dtt):+.1f} s)")
    co = []
    for model in sorted(pk.model.unique()):
        for thr in THRS:
            pke = drop_day_edges(pk)
            a = np.sort(pke[(pke.sta == "WISH") & (pke.model == model) & (pke.phase == "S")
                            & (pke.thr_first <= thr) & (pke.prob >= thr)].t.to_numpy())
            b = np.sort(pke[(pke.sta == "B943") & (pke.model == model) & (pke.phase == "S")
                            & (pke.thr_first <= thr) & (pke.prob >= thr)].t.to_numpy())
            lags, cc = xcorr_picks(a, b)
            if len(lags) == 0:
                continue
            bg = np.median(cc)
            k = int(np.argmax(cc))
            inband = (lags >= np.nanpercentile(dtt, 5)) & (lags <= np.nanpercentile(dtt, 95))
            kb = int(np.argmax(np.where(inband, cc, -np.inf)))
            co.append(dict(model=model, threshold=thr, n_WISH=len(a), n_B943=len(b),
                           day_edges_dropped_s=EDGE_S,
                           peak_lag_s=lags[k], peak_over_background=cc[k] / max(bg, 1e-9),
                           best_inband_lag_s=lags[kb],
                           inband_peak_over_background=cc[kb] / max(bg, 1e-9),
                           background_pairs=bg))
            np.savez(f"{RES}/e2_offgt_xcorr_{model}_t{thr}.npz", lags=lags, cc=cc,
                     dtt=dtt, slab_z=z)
    co = pd.DataFrame(co)
    co.to_csv(f"{RES}/e2_offgt_coincidence.csv", index=False)
    print("\n[coincidence WISH x B943]")
    print(co.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
