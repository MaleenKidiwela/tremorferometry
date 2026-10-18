# NOTES — rf_gosselin (times PT)

## 2026-10-13
- ~9:45 am: read paper + Materials & Methods. rfpy 0.1.2 (Audet's package, the paper's group) installed with
  `pip --user --no-deps`; its rotate/SNR/Wiener deconvolution and harmonics code are copied for our RFs.
- Event list: ComCat M>=5.8 1993-01-01..2026-10-13 (~7,100 events); per station 30-90° arrivals with BH/HH epochs:
  PGC 4059, SNB 2566, VGZ 2860, LZB 1784, PFB 2278; extension NLLB 2856, YOUB 1595, MGB 1632, GOBB 971, GNW 3803,
  SP2 3382, LON 3955, OFR 1499, WISH 2800, NLWA 2527.
- Broadband epochs (NRCAN = IRIS metadata): PGC BH 1993-07..2017-08 then HH; SNB 2000-05..2021-03; VGZ 2003-02..;
  LZB 2001-01..2014-06; PFB 2001-05..2018-12. Before 2000/2003 SNB and VGZ only have short-period EHZ at NRCAN, so
  SNB and VGZ cannot reach back to 1993 (the paper's 1993-2019 is PGC's span).
- Bug 1 (download): the EarthScope copy of CN returns no data for most pre-2010 windows (204). Switched CN to the
  NRCAN/CNDC FDSN (the paper's source). NRCAN answers HTTP 404 for "no data" -> treated as empty. The first
  (empty) EarthScope chunk markers were moved to `_chunks_iris_empty_20261013/` (not deleted).
- Bug 2 (geometry): rf00 stored the event->station azimuth as 'baz' (gps2dist_azimuth index [2]); symptom: SV RFs
  had −0.6×P at zero lag at VGZ. Fixed in rf00 (index [1]) and recomputed in rf03 (`baz_true`, old kept as
  `baz_rf00`). After the fix SV(0)/P(0) = −0.01 and VGZ shows the LVL pair (negative ~3.8 s, positive ~4.8 s).
- Near-surface velocities for the P-SV-SH rotation: rfpy defaults vp 6.0, vs 3.5 (residual P on SV −7 % at VGZ;
  vp 5.0 / vs 2.9 gives 0 %). Kept the defaults.
- Response: only sensitivity + orientation applied (RF = spectral ratio of components of one sensor, so the
  response cancels); PLAN said "response removed" — this is the change, logged here.
- ETS midpoints (rf02, paper recipe on UNR IGS20 ALBH east): anchors = GPS detections 1994-2002 (picked to keep
  the ~14-16 month cadence: 1994-07, 1995-12, 1997-04, 1998-09, 1999-08, 2000-12, 2002-02; rejected 1995-01,
  1996-06, 2012-02 as off-cadence) + repo SVI big-ETS list 2003-2026 (r23 study4). 1994-01..2019-06: **21 episodes,
  the same count as the paper**. 30 episodes 1994-2025. WA: NW big list (b_ets_record) midpoints refined on SEAT.
- Smoke test (VGZ, partial download 2003-2009, 300 draws): ΔCC_SV (15 d) = −0.19, at the 8th percentile of the
  scrambled control; not a result, only a code check.
- 10:06 am: download running (3 workers), processing queue `scripts/rf_queue.sh` launched (waits for markers).
- ~10:40 am: core results in (RESULTS.md). PGC reproduces ΔCC_SV −0.12 (pct 0.002) and replicates out of sample 2019–26
  (−0.094, pct 0.04); VGZ same sign, not significant; SNB weak; dVs size not reproduced (scale noisy). target.json written
  from core stages; recovery/extension/slip-front still in the queue.

## 2026-10-14
- 1:13 pm PT: queue finished (RF_EXT, EXT_SVI, EXT_WA, REC_EXT_SVI, REC_EXT_WA, SLIP_SVI, SLIP_WA); REC_REPRO/REC_FULL
  (done 10-13) also written up. RESULTS §6-11 added; §5 superseded by §9.
- US.NLWA (re-downloaded after the '00'->'0.0' loc fix): 2,312/2,527 events, 1,343 RFs > 5 dB, median CC_SV 0.77: sensible.
- LVL picks failed (window edge / no negative pulse) at SP2, LON, NLWA, OFR, WISH -> their dVs columns are not LVL values.
- Quick check (not a queue stage): `scripts/rf06_align_check.py` -> `tables/align_check_15d.csv` (re-bins existing per-RF CC
  aligned on midpoint / median / first local tremor). PGC: -0.133 / -0.128 / -0.013.
- target.json updated from the recovery + extension + slip-front stages (previous version kept as target_v1_20261013.json).
- Open (not run): a season-preserving scrambled control (shuffle within calendar month); ETS midpoints cluster by season.
- ~2:45 pm PT. PLAN Amendment 2 written (same-month scrambled control, seasonal regression, rules S1–S3) + pointer to the
  acf_pgc §10 C2 pre-window test, before computing.
- ~2:47 pm PT. `scripts/rf07_samemonth.py` written; smoke run (500 draws, PGC + LON, `tables/samemonth/smoke_*`). The first full run
  crashed on oos2019 (no recovery table for that tag); the `touch SAMEMONTH_DONE` in that command line ran anyway, so the
  marker existed for ~2 min before the good run (fixed the missing-file case; rerun, marker re-touched). Then the seasonal panel's
  month axis was shifted by half a month (plot only); rerun, tables byte-identical.
- Results (RESULTS §12–13): PGC passes the same-month control (pct 0.003 repro / 0.0004 full); LON survives and is not seasonal
  (pct 0.0016, season-predicted Δ 0.000, no seasonal cycle at LON); GNW's opposite sign drops to two-sided 0.088; seasonal cycles
  exist (PGC annual 0.03, F p 0.004) but predict ≤ 0.012 of any 15-d Δ. Pre-ETS H on unused stations: fail (p 0.90).
