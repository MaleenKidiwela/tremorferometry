#!/usr/bin/env python
"""Adaptive candidate P-threshold (fleet-wide policy): the B011-trained picker's P scores drift LOW at
other stations, so a fixed --thr 0.7 STARVES sparse/low-calibration stations (unrecoverable without re-run).
Since causality certifies reliability at the END, cast a wider net: pick the P threshold giving ~TARGET
candidates, floored 0.2, capped 0.7 -> thr = max(0.2, min(0.7, P_that_gives_TARGET)). Re-thresholds from the
already-scored baseline parquet (no re-scoring). Floor lowered 0.5->0.2 (2026-07-13): 0.5 starved low-
calibration stations to 0 families (B045 got 1046 cand, B030 got 34) even though they have real GOLD LFEs;
the P>=0.5 floor was calibrated for B011-like sites. The 300 hard-cap bounds densify cost regardless.
Usage: python adaptive_cand_threshold.py <STA> [TARGET=15000]  -> overwrites data/<sta>_cand_filtered.parquet
"""
import sys, numpy as np, pandas as pd

sta = sys.argv[1].lower()
TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 15000
base = pd.read_parquet(f"data/{sta}_cand_baseline.parquet")
p = np.sort(base.p_lfe.values)[::-1]
thr_raw = 0.0 if len(p) <= TARGET else float(p[TARGET - 1])   # P of the TARGET-th highest-scoring candidate
thr = max(0.05, min(0.7, thr_raw))  # floor 0.05 = GARBAGE bound only; selection is rank-based by design.
# History: 0.5 (starved low-calibration boreholes) -> 0.2 (2026-07-13) -> 0.05 (2026-07-15, Merlin): the 0.2
# floor was calibrated on BROADBAND score magnitudes and was already inert there (SHB 0.19, CLRS 0.20), but
# SHORT-PERIOD EHZ scores compress below it (SMW P_for_30k=0.10) -> 0.2 starved 208 EHZ fleet stations to
# ~0 candidates. Absolute P does NOT transfer across instrument classes (isotonic out-of-domain); RANKS do.
# Causality (coda/mirror>1.5) + the frozen >=20/>=15% gate is the referee; a junk station FLAGs there, cheaply.
filt = base[base.p_lfe >= thr]
filt.to_parquet(f"data/{sta}_cand_filtered.parquet", index=False)
print(f"{sta}: baseline {len(base)} cand, P_for_{TARGET}={thr_raw:.3f} -> thr={thr:.3f} "
      f"-> {len(filt)} candidates -> data/{sta}_cand_filtered.parquet")
