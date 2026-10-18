#!/usr/bin/env python
"""WA preparation (1) -- STATION-HOUR HEALTH CENSUS on the Washington caches.  Prefix `b_`.

Merlin's Step-1 audit defect F1: `C8.SHVB`'s horizontals carried a 4.76 Hz tone on 16 of 24 ETS
days and sat in every sum undetected, and "the WA bed has the same gap".  His Step-2 rule, applied
here verbatim to every station-hour of the WA caches (quiet and tremor alike):

    veto if ANY component is flat, OR the horizontal/vertical rms ratio is outside [1/30, 30],
    OR the 1-8 Hz Welch PSD peak/median exceeds 1e4 (tone).

The cache is the frozen 20 Hz 1-8 Hz chain, so the whole Welch band IS the passband and
peak/median is a clean tone statistic.

Outputs
  results/wa_station_hour_veto_<year>.csv   one row per station-hour with the reason codes
  results/wa_station_veto_summary_<year>.csv  per station: n_hours, % vetoed, reasons
  results/wa_array_veto_summary_<year>.csv    per array (station-code prefix): % vetoed

Usage
  python scripts/b_wa_health.py --year 2010 --workers 24
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
CACHE = f"{STUDY}/cache20hz"
RES = f"{STUDY}/results"
FS = 20.0
HOUR_N = int(3600 * FS)
HV_LO, HV_HI = 1.0 / 30.0, 30.0
TONE_RATIO = 1e4
FLAT_CONST_FRAC = 0.5          # > 50 % of the hour equal to the median = flat / gap-filled


def array_of(netsta):
    """Array label: XG/XU sub-array code (BH, CL, DR, GC, LC, PA, SL, TB, BS, PL, W0xx), else net."""
    net, sta = netsta.split(".")
    if net in ("XG", "XU"):
        if sta[:2].isalpha() and sta[0] != "W":
            return f"{net}.{sta[:2]}"
        return f"{net}.{sta[0]}"
    return net


def one_file(p):
    try:
        d = np.load(p)
        x = d["x"]
        order = [str(o) for o in d["order"]]
    except Exception as e:                                                  # noqa: BLE001
        return [dict(netsta=p.split("/")[-3], date="", hour=-1, error=f"{type(e).__name__}")]
    from scipy.signal import welch
    netsta = p.split("/")[-3]
    year = int(p.split("/")[-2])
    doy = int(os.path.basename(p).split(".")[0])
    date = str((pd.Timestamp(year=year, month=1, day=1) + pd.Timedelta(days=doy - 1)).date())
    iz = 0 if not order else int(np.argmax([o.endswith("Z") for o in order]))
    rows = []
    n = x.shape[1]
    for h in range(24):
        a, b = h * HOUR_N, min((h + 1) * HOUR_N, n)
        if b - a < HOUR_N // 2:
            rows.append(dict(netsta=netsta, date=date, hour=h, n_samples=int(max(b - a, 0)),
                             n_flat=3, hv_ratio=np.nan, psd_peak_med=np.nan,
                             veto=True, reason="short"))
            continue
        seg = np.asarray(x[:, a:b], dtype=np.float64)
        bad = ~np.isfinite(seg)
        if bad.any():
            seg = np.where(bad, 0.0, seg)
        n_flat = 0
        for c in range(seg.shape[0]):
            v = seg[c]
            if not np.any(v) or float(v.std()) == 0.0:
                n_flat += 1
                continue
            med = np.median(v)
            if float(np.mean(v == med)) > FLAT_CONST_FRAC:
                n_flat += 1
        rms = np.sqrt((seg ** 2).mean(axis=1))
        ih = [c for c in range(seg.shape[0]) if c != iz]
        rz = float(rms[iz])
        rh = float(np.sqrt(np.mean(rms[ih] ** 2))) if ih else np.nan
        hv = (rh / rz) if rz > 0 else np.inf
        pk = np.nan
        try:
            f, P = welch(seg, fs=FS, nperseg=1024, axis=-1)
            m = (f >= 1.0) & (f <= 8.0)
            with np.errstate(divide="ignore", invalid="ignore"):
                r = P[:, m].max(axis=1) / np.maximum(np.median(P[:, m], axis=1), 1e-300)
            pk = float(np.nanmax(r))
        except Exception:                                                   # noqa: BLE001
            pass
        reasons = []
        if n_flat:
            reasons.append("flat")
        if not np.isfinite(hv) or hv < HV_LO or hv > HV_HI:
            reasons.append("hv")
        if np.isfinite(pk) and pk > TONE_RATIO:
            reasons.append("tone")
        rows.append(dict(netsta=netsta, date=date, hour=h, n_samples=int(b - a), n_flat=n_flat,
                         hv_ratio=float(hv) if np.isfinite(hv) else np.inf,
                         psd_peak_med=pk, veto=bool(reasons), reason="+".join(reasons)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    files = sorted(glob.glob(f"{CACHE}/*/{a.year}/*.npz"))
    if a.limit:
        files = files[:a.limit]
    print(f"[health] {a.year}: {len(files)} station-days, {a.workers} workers", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for k, r in enumerate(ex.map(one_file, files, chunksize=4)):
            rows += r
            if (k + 1) % 500 == 0:
                print(f"   {k+1}/{len(files)} files", flush=True)
    D = pd.DataFrame(rows)
    D["array"] = [array_of(s) for s in D.netsta]
    out = f"{RES}/wa_station_hour_veto_{a.year}.csv"
    D.to_csv(out, index=False)
    print(f"[health] {len(D)} station-hours -> {out}")

    def summ(g):
        return pd.Series(dict(
            n_hours=len(g), n_veto=int(g.veto.sum()), pct_veto=100.0 * float(g.veto.mean()),
            pct_flat=100.0 * float(g.reason.str.contains("flat").mean()),
            pct_hv=100.0 * float(g.reason.str.contains("hv").mean()),
            pct_tone=100.0 * float(g.reason.str.contains("tone").mean()),
            pct_short=100.0 * float(g.reason.str.contains("short").mean())))

    S = D.groupby("netsta").apply(summ, include_groups=False).reset_index()
    S["array"] = [array_of(s) for s in S.netsta]
    S = S.sort_values("pct_veto", ascending=False)
    S.to_csv(f"{RES}/wa_station_veto_summary_{a.year}.csv", index=False)
    A = D.groupby("array").apply(summ, include_groups=False).reset_index().sort_values(
        "pct_veto", ascending=False)
    A["n_stations"] = [int((S.array == x).sum()) for x in A.array]
    A["n_stations_gt30pct"] = [int(((S.array == x) & (S.pct_veto > 30)).sum()) for x in A.array]
    A.to_csv(f"{RES}/wa_array_veto_summary_{a.year}.csv", index=False)
    print(A.round(2).to_string(index=False))
    hot = A[A.pct_veto > 20]
    if len(hot):
        print("\n[FLAG] arrays with > 20 % vetoed station-hours (Merlin: the site-probe admission "
              "was blind to this and the coverage numbers are overstated):")
        print(hot[["array", "pct_veto", "n_stations", "n_stations_gt30pct"]].round(2)
              .to_string(index=False))
    print("\nstations with > 30 % vetoed hours (to be excluded from the coverage recompute):")
    print(S[S.pct_veto > 30][["netsta", "array", "n_hours", "pct_veto", "pct_flat", "pct_hv",
                              "pct_tone"]].round(1).to_string(index=False))


if __name__ == "__main__":
    main()
