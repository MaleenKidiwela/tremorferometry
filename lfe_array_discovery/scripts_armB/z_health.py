#!/usr/bin/env python
"""STEP 1c' STAGE 2 -- STATION-HOUR HEALTH VETO (Merlin F1).

C8.SHVB's horizontals carry a 4.76 Hz tone on most of the 2005 ETS window; it sat in every sum of
Step 1c because neither the engine nor the bed assembly has a station-health screen.  This stage
builds one, for EVERY station-hour of the bed (tremor days AND quiet days), from the frozen 20 Hz
1-8 Hz cache:

  flat   any component whose hour is constant, or more than FLAT_ZERO_FRAC exactly-zero samples;
  hv     horizontal / vertical rms ratio outside [1/30, 30];
  tone   1-8 Hz Welch PSD peak / median above 1e4 on any component.

Writes results/step1c_prime/station_hour_veto.csv (sta, date, hour, stats, reason, veto) and a
census json.  ACCEPT (pre-registered): SHVB flagged on doy 241-251 and 255-261, the other dense
sites unflagged, and no dense site losing more than 30 % of its hours.  A dense site above 30 %
means the bed shrinks -> the queue stops and the spread rule is recomputed and reported.
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
import w_common as W                                                        # noqa: E402

FS = W.FS
HV_LO, HV_HI = 1.0 / 30.0, 30.0
TONE_RATIO = 1e4
FLAT_ZERO_FRAC = 0.5
BAND = (1.0, 8.0)
DENSE_MAX_VETO_FRAC = 0.30
# the sites Merlin names: SHVB is EXPECTED to be flagged (it is the defect), the others are the
# bed's working dense sites and must not be
EXPECT_FLAGGED = ["SHVB"]
EXPECT_CLEAN = ["GLBC", "TWKB", "MGCB", "TSJB", "TWBB"]
SHVB_BAD_DOYS = list(range(241, 252)) + list(range(255, 262))


def station_hours(sta, days):
    from scipy.signal import welch
    rows = []
    for d in days:
        ts = pd.Timestamp(d)
        v = W.load_day(sta, int(ts.year), int(ts.dayofyear))
        if v is None:
            continue
        t0, X = v
        nH = int(round(3600 * FS))
        for h in range(24):
            tgt = float(ts.timestamp()) + h * 3600.0
            i0 = int(round((tgt - t0) * FS))
            if i0 < 0 or i0 + nH > X.shape[1]:
                continue
            x = X[:, i0:i0 + nH].astype(np.float64)
            if not np.isfinite(x).all():
                x = np.nan_to_num(x)
            sd = x.std(1)
            zf = (x == 0).mean(1)
            flat = bool((sd <= 0).any() or (zf > FLAT_ZERO_FRAC).any())
            rz = float(np.sqrt((x[0] ** 2).mean()))
            rh = float(np.sqrt((x[1] ** 2 + x[2] ** 2).mean()))
            hv = rh / rz if rz > 0 else np.inf
            pk = 0.0
            for c in range(3):
                if sd[c] <= 0:
                    continue
                f, p = welch(x[c], fs=FS, nperseg=2048)
                m = (f >= BAND[0]) & (f <= BAND[1])
                med = float(np.median(p[m]))
                if med > 0:
                    pk = max(pk, float(p[m].max() / med))
            reason = []
            if flat:
                reason.append("flat")
            if not (HV_LO <= hv <= HV_HI):
                reason.append("hv")
            if pk > TONE_RATIO:
                reason.append("tone")
            rows.append(dict(sta=sta, date=ts.strftime("%Y-%m-%d"), doy=int(ts.dayofyear), hour=h,
                             rms_z=rz, rms_h=rh, hv_ratio=float(hv), psd_peak_over_median=pk,
                             frac_zero=float(zf.max()), flat=flat,
                             reason="+".join(reason), veto=bool(reason)))
    return rows


def _work(args):
    sta, days = args
    try:
        return station_hours(sta, days)
    except Exception as e:                                                  # noqa: BLE001
        W.log(f"[health] {sta} FAILED: {e}")
        return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{W.RES}/step1c_prime/station_hour_veto.csv")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    mem, of, order, det = W.sites()
    stas = sorted(set(det) | set(W.WITHHELD))
    days = sorted({d.strftime("%Y-%m-%d") for d in W.bed_days()}
                  | set(W.STEP1C_QUIET_DAYS))
    W.log(f"[health] {len(stas)} stations x {len(days)} days")
    import multiprocessing as mp
    with mp.Pool(min(a.workers, 30)) as pool:
        parts = pool.map(_work, [(s, days) for s in stas])
    R = pd.DataFrame([r for p in parts for r in p])
    R = R.sort_values(["sta", "date", "hour"]).reset_index(drop=True)
    R.to_csv(a.out, index=False)

    cen = R.groupby("sta").agg(n_hours=("veto", "size"), n_veto=("veto", "sum")).reset_index()
    cen["frac_veto"] = cen.n_veto / cen.n_hours.clip(lower=1)
    cen.to_csv(a.out.replace(".csv", "_census.csv"), index=False)
    fr = {r.sta: float(r.frac_veto) for r in cen.itertuples()}
    # SHVB on the doys Merlin names
    sh = R[(R.sta == "SHVB") & (R.doy.isin(SHVB_BAD_DOYS))]
    shvb_hit = float(sh.veto.mean()) if len(sh) else np.nan
    dense = [s for s in EXPECT_CLEAN if s in fr]
    over = {s: fr[s] for s in dense if fr[s] > DENSE_MAX_VETO_FRAC}
    reasons = R[R.veto].groupby(["sta", "reason"]).size().rename("n").reset_index()
    out = dict(
        n_station_hours=int(len(R)), n_vetoed=int(R.veto.sum()),
        frac_vetoed=float(R.veto.mean()) if len(R) else 0.0,
        frac_veto_by_station=fr,
        reasons=reasons.to_dict("records"),
        shvb_veto_frac_on_bad_doys=shvb_hit,
        shvb_flagged_as_expected=bool(np.isfinite(shvb_hit) and shvb_hit >= 0.5),
        dense_sites_checked=dense, dense_sites_over_30pc=over,
        thresholds=dict(hv_range=[HV_LO, HV_HI], tone_psd_peak_over_median=TONE_RATIO,
                        flat_zero_frac=FLAT_ZERO_FRAC, band=list(BAND),
                        dense_max_veto_frac=DENSE_MAX_VETO_FRAC),
        expect_flagged=EXPECT_FLAGGED, expect_clean=EXPECT_CLEAN,
        verdict=("ABANDON_BED_SHRINKS" if over else "ACCEPT"),
        reason=("dense site(s) lose more than 30 % of their hours: "
                f"{over} -- recompute the moveout-spread rule and report" if over else ""))
    json.dump(out, open(a.out.replace(".csv", ".json"), "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "reasons"}, indent=1, default=str))


if __name__ == "__main__":
    main()
