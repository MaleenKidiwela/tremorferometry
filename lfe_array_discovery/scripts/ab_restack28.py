"""Restack the converged NW families at 2-8 Hz (2026-09-26, user-approved).

Same arithmetic as ab_grow.restack / ab_grow_cont.restack_parallel (every pool station, 44 s window
at LONG_PRE, unit-max-3C per detection, mean over detections), except each detection segment is
band-passed 2-8 Hz (4th-order Butterworth, zero-phase) BEFORE the unit-max normalisation.
Reason: the 1-8 Hz stacks peak at 1.1-1.3 Hz (noise at the band edge) and carry narrow tones;
see notes/2026-09-25_Notes.md.  The 1-8 Hz stacks stay as the literal comparison.

  python scripts/ab_restack28.py --sub NW --det results/ab_royer/NW/grow/det_it14.parquet \
      --out results/ab_royer/NW/hypo_28/stacks_loc.npz
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt

sys.path.insert(0, os.path.dirname(__file__))
import ab_common as A  # noqa: E402
import ab_grow as G  # noqa: E402

SOS = butter(4, [2.0 / (A.FS / 2), 8.0 / (A.FS / 2)], "band", output="sos")
_orig = G.unit_max3


def unit_max3_28(seg):
    return _orig(sosfiltfilt(SOS, seg, axis=-1))


G.unit_max3 = unit_max3_28          # patched before the worker pool forks
import ab_grow_cont as GC  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", default="NW")
    ap.add_argument("--det", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args()
    det = pd.read_parquet(a.det)
    stas = A.POOL[a.sub]
    A.log(f"[restack28] {det.grp.nunique()} families, {len(det)} detections, {len(stas)} stations")
    ST = GC.restack_parallel(stas, det, workers=a.workers)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    np.savez_compressed(a.out, **{f"{g}_{s}": ST[g][s].astype(np.float32)
                                  for g in ST for s in ST[g]})
    A.log(f"[restack28] wrote {a.out} ({len(ST)} families)")


if __name__ == "__main__":
    main()
