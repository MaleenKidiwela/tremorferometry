#!/usr/bin/env python
"""STAGE 6 -- locate the Royer templates with the REBUILT locator, and compare with the paper.

The locator is the one calibrated on Bostock's published coordinates in
`results/WA_LOCATION_REBUILD.md`: the adopted variant **L9_xcdiff** (4.0 km epicentral median, 6.4 km
3-D median, rail 0.09 against Bostock's 130 published family coordinates, versus 27.6 km for the old
`k_locate`) plus **L12_comb**, which locates more families at a larger p90.  Both are run, every
family keeps its row, and each carries its station-bootstrap radius (median 8.0 km in the north,
correlation 0.81 with the true error) -- `rank-dont-discard`.

Royer's own values, which are the REFERENCE LINES and never gates (sec.5):
  Hyp2000 then hypoDD (Table S1: IDAT 3, 5 x NITER 8 = 40 iterations, DAMP 50,50,40,40,40;
  WTCCP/WTCCS 0.01,0.01,1,1,1; WTCTP/WTCTS 1,1,0.01,0.01,0.01; WRCC -9,-9,-9,5,5;
  WDCC -9,-9,2,2,0.5; WRCT -9,-9,6,6,6; WDCT -9,8,8,8,8)
  depths 31-46 km (NW), 32-49 km (SW); mean nominal Hyp2000 uncertainty 2.0/1.3 km (NW),
  3.7/3.5 km (SW); hypoDD retains 77 % (NW) and 72 % (SW) of the events.

WHY THE INPUT IS DECIMATED TO 20 Hz.  `uu_common.FS` is 20.0 and the whole ladder -- envelope
smoothing, the cross-correlation differential-time window, the pick-quality columns -- was tuned and
measured at that rate.  Feeding it 40 sps data would be an uncalibrated change to a calibrated
instrument, so the Royer 40 sps cache is decimated to 20 Hz for the location stage only (the
detection stages stay at the paper's 40 sps).  A 40 Hz run is available as a reported ablation
(--fs 40) and is never the headline.

  python scripts/ab_locate.py --sub NW --det results/ab_royer/NW/grow/det_it3.parquet \
      --tpl results/ab_royer/NW/grow/tpl_it3.npz --out results/ab_royer/NW/locate
-> <out>/stacks.npz, <out>/picks.parquet, <out>/loc_L9_xcdiff.csv, <out>/loc_L12_comb.csv
-> <out>/depth_summary.csv   located families, depth percentiles, fraction inside Royer's band
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import resample_poly

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402
import uu_common as U                                                       # noqa: E402

WA_GRID = dict(lat=(46.30, 49.00), lon=(-125.60, -121.60), dep=(6.0, 76.0))
VARIANTS = "L9_xcdiff,L12_comb"


class RoyerCache:
    """ab_common.DayCache resampled to the locator's own rate, mask-aware.

    Only complete-and-valid segments are returned, exactly as ab_common.DayCache.cut does; the mask
    is carried by cutting at 40 sps first and decimating the segment, so an invalid sample can never
    be interpolated into a valid one."""

    def __init__(self, stas, fs_out=U.FS):
        self.dc = A.DayCache(stas)
        self.fs_out = float(fs_out)
        self.q = int(round(A.FS / self.fs_out))
        if abs(A.FS / self.fs_out - self.q) > 1e-9:
            raise SystemExit(f"{A.FS} -> {fs_out} is not an integer decimation")

    def day(self, date):
        return self.dc.day(date)

    def cut(self, sta, t_abs, len_s, date=None):
        seg = self.dc.cut(sta, t_abs, len_s, date)
        if seg is None:
            return None
        if self.q == 1:
            return seg
        out = resample_poly(seg, 1, self.q, axis=-1)
        n = int(round(len_s * self.fs_out))
        return out[:, :n] if out.shape[1] >= n else None


class RoyerMod:
    """The duck-typed module `uu_common.Bed` expects."""

    def __init__(self, sub, fs_out):
        self.sub = sub
        self.fs_out = fs_out

    def activate(self):
        return None

    def station_table(self):
        return A.station_table()

    def stations(self):
        return list(A.SUBARRAY[self.sub]["stas"])

    def bed_days(self):
        return A.days(self.sub)

    def DayCache(self, stas):
        return RoyerCache(stas, self.fs_out)

    def station_corrections(self):
        """SMOKE-TEST FIX: uu_common.station_terms_published calls this on non-a5 beds and the
        duck-typed module lacked it (AttributeError).  Same convention as l_common: the intraslab
        S terms of results/x_station_corrections.csv; missing/unusable -> 0.0 and FLAGGED."""
        p = f"{A.STUDY}/results/x_station_corrections.csv"
        C = pd.read_csv(p).set_index("sta") if os.path.exists(p) else pd.DataFrame()
        corr, flagged = {}, []
        for s in A.COORD:
            ok = s in C.index and bool(C.loc[s, "usable"])
            v = float(C.loc[s, "corr_S"]) if ok else 0.0
            corr[s] = v if np.isfinite(v) else 0.0
            if not ok:
                flagged.append(s)
        return corr, flagged


def install_bed(sub, fs_out):
    """Register a `royer` bed on uu_common WITHOUT editing the calibrated locator code."""
    mod = RoyerMod(sub, fs_out)
    orig = U.bed

    def bed(name="north"):
        if name == "royer":
            return U.Bed(mod, "royer", WA_GRID, "intraslab")
        return orig(name)

    U.bed = bed
    return mod


def depth_summary(sub, rows):
    lo, hi = A.SUBARRAY[sub]["depth_pub"]
    out = []
    for name, R in rows:
        L = R[R.located == True]                                            # noqa: E712
        d = L.dep.to_numpy(float) if len(L) else np.array([])
        out.append(dict(
            variant=name, subarray=sub, n_families=int(len(R)), n_located=int(len(L)),
            rail=float(L.edge.mean()) if len(L) else np.nan,
            dep_p05=float(np.percentile(d, 5)) if d.size else np.nan,
            dep_med=float(np.median(d)) if d.size else np.nan,
            dep_p95=float(np.percentile(d, 95)) if d.size else np.nan,
            frac_in_published_band=float(((d >= lo) & (d <= hi)).mean()) if d.size else np.nan,
            published_band=f"{lo}-{hi}",
            boot_km_med=float(L.boot_km.median()) if "boot_km" in L and len(L) else np.nan,
            boot_km_p90=float(L.boot_km.quantile(0.9)) if "boot_km" in L and len(L) else np.nan,
            published_sigma_h=A.SUBARRAY[sub]["sig_h_pub"],
            published_sigma_z=A.SUBARRAY[sub]["sig_z_pub"],
            published_n_templates=A.SUBARRAY[sub]["n_tpl_pub"]))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", required=True, choices=["NW", "SW", "NVI"])
    ap.add_argument("--det", required=True, help="det_it<K>.parquet from ab_grow")
    ap.add_argument("--tpl", required=True, help="tpl_it<K>.npz from ab_grow (station set)")
    ap.add_argument("--out", default="")
    ap.add_argument("--fs", type=float, default=U.FS, help="20 = the calibrated rate (headline)")
    ap.add_argument("--boot", type=int, default=60)
    ap.add_argument("--variants", default=VARIANTS)
    a = ap.parse_args()
    out = a.out or f"{A.RES}/{a.sub}/locate"
    os.makedirs(out, exist_ok=True)
    install_bed(a.sub, a.fs)
    import uu_locate, uu_picks, uu_stacks                                   # noqa: E402

    stacks = f"{out}/stacks.npz"
    picks = f"{out}/picks.parquet"
    if not os.path.exists(stacks):
        sys.argv = ["uu_stacks", "--bed", "royer", "--det", a.det, "--tpl", a.tpl,
                    "--out", stacks]
        uu_stacks.main()
    if not os.path.exists(picks):
        sys.argv = ["uu_picks", "--stacks", stacks, "--out", picks, "--jobs", "12"]
        uu_picks.main()
    sys.argv = ["uu_locate", "--bed", "royer", "--picks", picks, "--stacks", stacks,
                "--outdir", out, "--only", a.variants, "--boot", str(a.boot)]
    uu_locate.main()
    rows = []
    for name in a.variants.split(","):
        p = f"{out}/loc_{name}.csv"
        if os.path.exists(p):
            rows.append((name, pd.read_csv(p)))
    if rows:
        S = depth_summary(a.sub, rows)
        S.to_csv(f"{out}/depth_summary.csv", index=False)
        A.log("\n" + S.to_string(index=False))


if __name__ == "__main__":
    main()
