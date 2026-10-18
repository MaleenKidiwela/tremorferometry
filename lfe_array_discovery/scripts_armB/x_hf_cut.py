#!/usr/bin/env python
"""Step 0 (b) -- the WASHINGTON HF-class VETO cut for Gate 0, derived BEFORE it is used.

Merlin F10: the frozen Gate-0 HF veto (`g_gates.GateConfig.g0_hf_station = ("PGC","B011")`,
cut 3.5) was calibrated at PGC on the NORTH.  It must move to a WA 100-Hz borehole, and that
station's own quiet-day HF-ratio distribution must be reported before the cut is applied.

The statistic is `g_common.hf_ratio_members` VERBATIM: per 12-s window on the NATIVE trace,
E(10-19 Hz) / E(2-8 Hz) summed over the three components.  Nothing here is a coherence measurement,
so random windows on quiet days are exactly the right sample for the site term.

THE RULE (written here, frozen, and not retuned afterwards)

    cut_WA = 3.5 x  median(HF ratio, quiet-day random windows, WA borehole)
                  / median(HF ratio, quiet-day random windows, PGC)

  The statistic is a per-site SPECTRAL RATIO, so the whole difference between two stations on
  noise is the site+instrument response term; the north cut is transferred through it.  It is NOT
  re-derived from a WA LFE population, because none exists yet -- that is the entire point of the
  run.  Equivalently: the frozen cut sits at 3.5 / median(PGC quiet) = ~1.5 x PGC's own quiet-noise
  level, and the WA cut is placed at the SAME multiple of the WA borehole's quiet-noise level.

  BINDING clause (the queue stops if it fails): the cut must be finite, positive, and ABOVE the WA
  borehole's quiet-day median -- a cut below the site's own noise level would veto everything.

  REPORTED, NOT USED: the same statistic on the catalogued intraslab earthquakes at the same
  station (`results/intraslab_eq_2006_2012.csv`, windows at the Crosson-predicted S).  Measured
  here at 0.04 -- i.e. M>=2 intraslab events at 25-50 km depth and ~60 km distance arrive HF-
  DEPLETED at this borehole and are NOT the class this veto separates.  The number is recorded so
  that the veto is never described as an "earthquake rejector" on this bed.

    python scripts/x_hf_cut.py [--n-per-day 40]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import x_common as X                                                        # noqa: E402

NORTH_CUT = 3.5                      # the frozen north veto (Track A, calibrated at PGC)
NORTH_REF_STA = "PGC"
NORTH_TRACES = ["/home/jovyan/tremorferometry/lfe_catalog_repair/traces",
                "/home/jovyan/tremorferometry/lfe_picker/traces",
                "/home/jovyan/tremorferometry/lfe_extension/traces"]


def _ratios(nstore, netsta, days, n_per_day, seed):
    """HF ratios of `n_per_day` random 12-s windows on each of `days` (list of Timestamps)."""
    import g_common as G
    rng = np.random.default_rng(seed)
    by = {}
    for d in days:
        t0 = float(pd.Timestamp(d).timestamp())
        ts = t0 + rng.uniform(300.0, 86000.0, n_per_day)
        by[(int(pd.Timestamp(d).year), int(pd.Timestamp(d).dayofyear))] = [float(t) + 1.0
                                                                          for t in ts]
    return G.hf_ratio_members(nstore, netsta, by)


def north_quiet_days(nstore, sta=NORTH_REF_STA, want=10):
    """Days of the north PGC trace tree with no PNSN tremor within 100 km and >= 3 d from one."""
    import g_common as G
    import s_common as S
    T = pd.read_csv(X.PNSN, usecols=["time", "lat", "lon"])
    T["time"] = pd.to_datetime(T.time, errors="coerce").astype("datetime64[ns]")
    T = T.dropna(subset=["time"])
    st = G.stations()
    la, lo = float(st.lat[sta]), float(st.lon[sta])
    d = S.epi_km(T.lat.to_numpy(), T.lon.to_numpy(), la, lo)
    td = np.array(sorted({pd.Timestamp(x).value for x in T.time[d <= 100.0].dt.floor("D")}),
                  dtype="int64")
    out = []
    for root in NORTH_TRACES:
        for yr in sorted(os.listdir(f"{root}/CN.{sta}") if os.path.isdir(f"{root}/CN.{sta}")
                         else []):
            if not yr.isdigit():
                continue
            for f in sorted(os.listdir(f"{root}/CN.{sta}/{yr}")):
                if not f.endswith(".mseed"):
                    continue
                t = pd.Timestamp(year=int(yr), month=1, day=1) + pd.Timedelta(days=int(f[:3]) - 1)
                gap = np.abs(td - t.value).min() / 86400e9 if td.size else 1e9
                if gap >= 3.0:
                    out.append(t)
        if len(out) >= want:
            break
    return sorted(set(out))[:want]


def eq_hf_ratios(E, netsta, slat, slon, win_s=12.0):
    """Same HF statistic on the +-3 min INTRASLAB-EARTHQUAKE windows (`traces_eq/<netsta>/<evid>`),
    read at the native rate and cut 1 s before the Crosson-predicted S."""
    import g_common as G
    import s_common as S
    from obspy import read
    out = []
    for e in E.itertuples():
        p = f"{X.TRACES_EQ}/{netsta}/{e.evid}.mseed"
        if not (os.path.exists(p) and os.path.getsize(p) > 0):
            continue
        try:
            st = read(p)
            st.merge(method=1, fill_value=0)
            fs = float(st[0].stats.sampling_rate)
            n = min(tr.stats.npts for tr in st)
            if len(st) < 3 or n < int(win_s * fs) + 10:
                continue
            Xw = np.vstack([tr.data[:n].astype(float) for tr in st[:3]])
            t0 = float(st[0].stats.starttime.timestamp)
            d = float(S.epi_km(e.lat, e.lon, slat, slon))
            ttS, _ = G.ttSP(d, float(e.depth_km))
            i0 = int(round((float(e.t.timestamp()) + ttS - 1.0 - t0) * fs))
            nn = int(round(win_s * fs))
            if i0 < 0 or i0 + nn > n:
                continue
            Xh = G._bp(Xw, fs, 10.0, 19.0)
            Xl = G._bp(Xw, fs, 2.0, 8.0)
            eh = float((Xh[:, i0:i0 + nn] ** 2).sum())
            el = float((Xl[:, i0:i0 + nn] ** 2).sum())
            if el > 0:
                out.append(eh / el)
        except Exception:                                                   # noqa: BLE001
            continue
    return np.asarray(out, float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-day", type=int, default=40)
    ap.add_argument("--n-eq", type=int, default=60)
    a = ap.parse_args()
    import g_common as G
    C = X.cfg()
    sta = C["hf_station"]
    ST = X.wa_station_table()
    netsta = f"{ST.net[sta]}.{sta}"

    # ---- WA quiet-day noise at the veto station
    nsW = G.NativeStore(trace_dirs=[X.TRACES])
    qd = [pd.Timestamp(d) for d in C["quiet_days"]]
    rw = _ratios(nsW, netsta, qd, a.n_per_day, 20260920)
    X.log(f"[hf] {netsta}: {rw.size} quiet-day windows over {len(qd)} days, "
          f"median {np.median(rw):.4f}")

    # ---- north reference (PGC) on the SAME statistic
    nsN = G.NativeStore(trace_dirs=NORTH_TRACES)
    nd = north_quiet_days(nsN)
    rn = _ratios(nsN, f"CN.{NORTH_REF_STA}", nd, a.n_per_day, 20260921)
    if rn.size < 100:
        raise SystemExit(f"north reference failed: only {rn.size} PGC quiet windows")
    X.log(f"[hf] CN.{NORTH_REF_STA}: {rn.size} quiet-day windows over {len(nd)} days, "
          f"median {np.median(rn):.4f}")

    cut = float(NORTH_CUT * np.median(rw) / np.median(rn))

    # ---- sanity clause (ii): the catalogued intraslab earthquakes at the same station
    E = pd.read_csv(f"{X.RES}/intraslab_eq_2006_2012.csv")
    E["t"] = pd.to_datetime(E.origin_utc, errors="coerce", format="mixed", utc=True)
    E = E.dropna(subset=["t"]).head(a.n_eq)
    import s_common as S
    req = eq_hf_ratios(E, netsta, float(ST.lat[sta]), float(ST.lon[sta]))
    eq_med = float(np.median(req)) if req.size >= 5 else np.nan
    X.log(f"[hf] intraslab-earthquake windows at {netsta}: n={req.size}, "
          f"median {eq_med if np.isfinite(eq_med) else float('nan'):.4f}")

    rep = dict(
        station=netsta, statistic="g_common.hf_ratio_members (12 s, E(10-19)/E(2-8), native rate)",
        rule=("cut_WA = 3.5 x median(quiet-day HF ratio, WA borehole) / median(quiet-day HF ratio, "
              "PGC); the north cut transferred through the per-site spectral-ratio term, because "
              "no WA LFE population exists to recalibrate on"),
        north_cut=NORTH_CUT, north_station=f"CN.{NORTH_REF_STA}",
        north_n=int(rn.size), north_days=[str(d.date()) for d in nd],
        north_median=float(np.median(rn)), north_p25=float(np.quantile(rn, .25)),
        north_p75=float(np.quantile(rn, .75)), north_p95=float(np.quantile(rn, .95)),
        wa_n=int(rw.size), wa_days=[str(d.date()) for d in qd],
        wa_median=float(np.median(rw)), wa_p25=float(np.quantile(rw, .25)),
        wa_p75=float(np.quantile(rw, .75)), wa_p95=float(np.quantile(rw, .95)),
        wa_p99=float(np.quantile(rw, .99)),
        eq_n=int(req.size), eq_median=eq_med,
        eq_p25=float(np.quantile(req, .25)) if req.size >= 5 else np.nan,
        cut=cut,
        cut_over_quiet_median=float(cut / np.median(rw)),
        north_cut_over_quiet_median=float(NORTH_CUT / np.median(rn)),
        binding_above_quiet_median=bool(np.isfinite(cut) and cut > 0
                                        and cut > float(np.median(rw))),
        reported_eq_median_below_cut=bool(np.isfinite(eq_med) and eq_med < cut),
        reported_eq_note=("intraslab M>=2 events at 25-50 km depth are HF-DEPLETED at this "
                          "borehole; reported for context, NOT a clause"),
    )
    json.dump(rep, open(f"{X.OUT}/hf_cut.json", "w"), indent=1)
    pd.DataFrame(dict(ratio=rw)).to_csv(f"{X.OUT}/hf_quiet_ratios_wa.csv", index=False)
    C["hf_veto"] = cut
    C["hf_veto_note"] = rep["rule"]
    C["hf_report"] = {k: v for k, v in rep.items() if not isinstance(v, list)}
    json.dump(C, open(X.BED_JSON, "w"), indent=1)
    X.log(f"[hf] CUT = {cut:.3f} (north 3.5 x {np.median(rw):.4f}/{np.median(rn):.4f}); "
          f"= {rep['cut_over_quiet_median']:.2f} x the WA quiet median "
          f"(north cut = {rep['north_cut_over_quiet_median']:.2f} x PGC's); "
          f"BINDING clause above-quiet-median {rep['binding_above_quiet_median']}; "
          f"reported intraslab-EQ median {eq_med:.3f}")
    if not rep["binding_above_quiet_median"]:
        raise SystemExit("HF cut fails the binding clause (not above the site quiet median)")


if __name__ == "__main__":
    main()
