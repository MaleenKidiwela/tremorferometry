#!/usr/bin/env python
"""Step 0 (c) -- WASHINGTON station corrections from CATALOGUED INTRASLAB EARTHQUAKES.

The `a5_eq_corrections.py` recipe (Merlin Q4; lfe_extension), with its hard-coded northern `STAS`
list replaced by this study's stations and its trace roots replaced by `traces_eq/` (the +-3 min
event windows downloaded as tranche T6):

  * catalogue: `results/intraslab_eq_2006_2012.csv` (USGS M >= 2, 25-50 km depth, in the WA box,
    2006-06..2012-06);
  * picks at the NATIVE sample rate, 2-15 Hz Butterworth-4 zero-phase (the 20 Hz / 1-8 Hz LFE cache
    is too narrow for an M2 P onset), P on Z within +-3 s and S on sqrt(H1^2+H2^2) within +-4 s of
    the Crosson-1-D prediction, envelope onset + AIC, accepted at peak/noise-rms >= --min-snr;
  * residual = observed - predicted; the correction is the MEDIAN residual over events.

USABILITY RULE (frozen before use, Merlin Q4 / the Step-0 spec): a station's correction is used
only if it has n >= 5 accepted S picks AND the residual MAD is <= 0.5 s.  Otherwise the correction
is 0.0 and the station is FLAGGED `usable=False` -- never silently imputed from a neighbour.

    python scripts/x_corrections.py [--min-snr 4.0] [--workers 12]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import x_common as X                                                        # noqa: E402

MIN_PICKS = 5
MAX_MAD_S = 0.5


def _one(args):
    sta, netsta, slat, slon, E, min_snr = args
    import g_common as G
    import s_common as S
    from a5_env import env_onset
    from obspy import read
    from scipy.signal import butter, sosfiltfilt
    rows = []
    for e in E.itertuples():
        p = f"{X.TRACES_EQ}/{netsta}/{e.evid}.mseed"
        if not (os.path.exists(p) and os.path.getsize(p) > 0):
            continue
        try:
            st = read(p)
            st.merge(method=1, fill_value=0)
            if len(st) < 3:
                continue
            fs = float(st[0].stats.sampling_rate)
            n = min(tr.stats.npts for tr in st)
            if n < int(40 * fs):
                continue
            # channel order: vertical first, the two horizontals after (no rotation is needed --
            # the S pick is made on sqrt(H1^2 + H2^2), which is rotation-invariant)
            tr_z = [tr for tr in st if tr.stats.channel.endswith("Z")]
            tr_h = [tr for tr in st if not tr.stats.channel.endswith("Z")]
            if len(tr_z) != 1 or len(tr_h) < 2:
                continue
            Xw = np.vstack([tr_z[0].data[:n].astype(float),
                            tr_h[0].data[:n].astype(float), tr_h[1].data[:n].astype(float)])
            t0 = float(tr_z[0].stats.starttime.timestamp)
            f2 = min(15.0, 0.45 * fs)
            sos = butter(4, [2.0 / (fs / 2), f2 / (fs / 2)], btype="band", output="sos")
            Xw = sosfiltfilt(sos, Xw, axis=-1)
            d = float(S.epi_km(e.lat, e.lon, slat, slon))
            ttS, ttP = G.ttSP(d, float(e.depth_km))
            tor = float(e.t.timestamp())
            iP0 = int(round((tor + ttP - t0) * fs))
            iS0 = int(round((tor + ttS - t0) * fs))
            nnoise = max(int(10 * fs), 50)
            if iP0 - int(3 * fs) <= nnoise or iS0 + int(4 * fs) >= n:
                continue
            iP, snrP = env_onset(Xw[0], iP0 - int(3 * fs), iP0 + int(3 * fs), nnoise, fs)
            H = np.sqrt(Xw[1] ** 2 + Xw[2] ** 2)
            iS, snrS = env_onset(H, iS0 - int(4 * fs), iS0 + int(4 * fs), nnoise, fs)
            rP = (t0 + iP / fs) - (tor + ttP) if iP is not None else np.nan
            rS = (t0 + iS / fs) - (tor + ttS) if iS is not None else np.nan
            rows.append(dict(evid=e.evid, sta=sta, dist_km=d, dep=float(e.depth_km),
                             mag=float(e.mag), fs=fs, ttP=ttP, ttS=ttS, snrP=snrP, snrS=snrS,
                             resP=rP if snrP >= min_snr else np.nan,
                             resS=rS if snrS >= min_snr else np.nan,
                             resP_raw=rP, resS_raw=rS))
        except Exception:                                                   # noqa: BLE001
            continue
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-snr", type=float, default=4.0)
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()
    from multiprocessing import Pool
    C = X.cfg()
    ST = X.wa_station_table()
    stas = sorted(set([s for q in C["site_order"] for s in C["site_members"][q]]
                      + list(C["withheld_stations"])))
    stas = [s for s in stas if s in ST.index]
    E = pd.read_csv(f"{X.RES}/intraslab_eq_2006_2012.csv")
    E["t"] = pd.to_datetime(E.origin_utc, errors="coerce", format="mixed", utc=True)
    E = E.dropna(subset=["t"]).reset_index(drop=True)
    X.log(f"[corr] {len(stas)} stations x {len(E)} intraslab earthquakes")
    jobs = [(s, f"{ST.net[s]}.{s}", float(ST.lat[s]), float(ST.lon[s]), E, a.min_snr)
            for s in stas]
    rows = []
    with Pool(min(a.workers, 24)) as pool:
        for k, r in enumerate(pool.imap_unordered(_one, jobs, chunksize=1), 1):
            rows += r
            if k % 10 == 0:
                X.log(f"[corr] {k}/{len(jobs)} stations, {len(rows)} station-events")
    P = pd.DataFrame(rows)
    P.to_csv(f"{X.OUT}/eq_picks.csv", index=False)
    out = []
    for s in stas:
        g = P[P.sta == s] if len(P) else P
        nS = int(g.resS.notna().sum()) if len(g) else 0
        nP = int(g.resP.notna().sum()) if len(g) else 0
        cS = float(g.resS.median()) if nS else np.nan
        cP = float(g.resP.median()) if nP else np.nan
        mS = float((g.resS - cS).abs().median()) if nS else np.nan
        mP = float((g.resP - cP).abs().median()) if nP else np.nan
        usable = bool(nS >= MIN_PICKS and np.isfinite(mS) and mS <= MAX_MAD_S)
        out.append(dict(sta=s, n_P=nP, n_S=nS, corr_P=cP if usable else 0.0,
                        corr_S=cS if usable else 0.0, raw_corr_P=cP, raw_corr_S=cS,
                        mad_P=mP, mad_S=mS, usable=usable,
                        reason=("" if usable else
                                ("n_S<5" if nS < MIN_PICKS else "MAD_S>0.5 s"))))
    Cc = pd.DataFrame(out)
    Cc.to_csv(X.CORR_CSV, index=False)
    # REPORTED, NOT USED: the same statistic pooled over each site's member stations.  The frozen
    # rule is per STATION (that is what `g_common.station_corrections` consumes), but the location
    # takes a SITE median pick, so the site-level MAD is the number that says whether the array
    # members agree; it is written here as evidence for the write-up and never fed to the engine.
    srow = []
    for q in C["site_order"]:
        mm = C["site_members"][q]
        g = P[P.sta.isin(mm)] if len(P) else P
        nS = int(g.resS.notna().sum()) if len(g) else 0
        cS = float(g.resS.median()) if nS else np.nan
        mS = float((g.resS - cS).abs().median()) if nS else np.nan
        srow.append(dict(site=q, n_members=len(mm), n_S=nS, corr_S=cS, mad_S=mS,
                         would_be_usable=bool(nS >= MIN_PICKS and np.isfinite(mS)
                                              and mS <= MAX_MAD_S)))
    pd.DataFrame(srow).to_csv(f"{X.OUT}/eq_corrections_by_site_REPORTED_NOT_USED.csv", index=False)
    X.log(f"[corr] {int(Cc.usable.sum())}/{len(Cc)} stations usable "
          f"(n_S >= {MIN_PICKS} and MAD <= {MAX_MAD_S} s); "
          f"median |corr_S| {float(Cc[Cc.usable].corr_S.abs().median()) if Cc.usable.any() else float('nan'):.2f} s")
    X.log(f"[corr] wrote {X.CORR_CSV}")


if __name__ == "__main__":
    main()
