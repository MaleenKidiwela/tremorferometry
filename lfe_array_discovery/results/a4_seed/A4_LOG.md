# A4 -- location-seeded positive control (and A6) -- running log

Agent started 2026-09-25 ~15:40Z.  Independent of the literal chain (queue pid 1928163, not touched);
only its 40 sps 1-8 Hz day caches (`cache40hz/`) are read.

## Code (all new; nothing of the chain edited)
- `scripts/a4f_common.py`, `a4f_grow.py`, `a4f_hypo.py`: FROZEN copies of `ab_common.py`
  (md5 16b14970af2792903571375a6e3fc02d), `ab_grow.py` (89b62582ae35bccab21cb0096a429ff4),
  `ab_hypo.py` (48eacf36e9785e391b7ac475c0f587f0) taken 2026-09-25 15:45Z, only the import names
  changed, so edits the running chain makes cannot change this queue mid-run.  Used: station
  coords/pools/episode days, cache paths, `window_stats`, the GSC model, `make_picks`,
  `write_hyp_inputs`, `run_hyp` (Hyp2000 v1.40, GSC, Vp/Vs 1.73, R 1.5 = the chain headline).
- `scripts/a4_make_locs.py` -> `locs_all.csv` (223 real + 223 nullR), `locs_pilot.csv` (10 per subarray + images).
- `scripts/a4_seed.py` the per-(subarray, condition) pipeline (docstring = full method).
- `scripts/a4_report.py` separation tables/curves/figure/REPORT.md.
- `scripts/a4_queue.sh` detached queue (flock, heartbeat, EXIT trap -> FAILED marker).
- `scripts/a6_prior_vs_truth.py` A6.

## Method choices (documented; each is a ranking/size, not a gate on candidates)
- Truth: Fig 4 Hyp2000 rows (NW 122, SW 54, NVI 47).  2 NW rows with no depth -> subarray median
  depth, flagged `depth_filled`; `depth_uncertain` = depth_spread_km > 1 or filled (13 rows).  Kept as columns.
- Days: the sec.3.1 autocorrelation episode days only (NW 75, SW 44, NVI 13).
- Stations: pool stations (Fig.1D pools) with cache on >= 1 episode day; nearest 10 carry the
  template (paper: 7-10), nearest 16 get diagnostic stacks for Hyp2000.  Elevation ignored in the
  predicted times.
- Travel times: GSC onshore 1-D (0-1 5.01, 1-6 6.01, 6-30 6.71, 30-46 7.11, >46 7.76), Vp/Vs 1.73,
  first arrival of direct ray / head waves (own implementation, `tt_1d`).
- Beam: local log-STA/LTA envelope (0.5 s smoothing, 30 s running mean of log env), S on
  horizontals + P on vertical at predicted moveout, mean over valid stations (>= 4), per-day robust
  z, 20 peaks/day 15 s apart, top 60 over all days.  NOTE: "tremor hours" were not separately
  selected -- all episode hours are ETS days; the beam ranking decides.
- Seed: pairwise network similarity (per-channel normalised cc, common lag +-1 s, mean over shared
  channels, >= 3 shared); seed = best mean-of-top-5; stack of it + 9 most similar (unit max 3C).
- Template window per station: [tP - 1 s, tP - 1 s + L], L = max(tS - tP) + 6 s rounded up to 2 s
  (so a few MF length groups).  Moveout preserved via per-station offsets.
- MF: per-channel normalised cc / valid channel count, INFINITE denominator on invalid windows (no
  clamp), 8 x MAD per day over samples with >= 4 valid stations, declustered by L, 3 passes
  (seed -> MF0 -> restack -> MF1 -> restack -> MF2 -> final 32 s stacks).  Restack = all detections.
- nullT: each station reads another episode day (seeded distinct whole-day shift within the
  subarray's day list) AND is circularly rolled in-day by (600 + 997 i) s -> every station has its
  own time map.  nullR: the whole template cloud rotated about its centroid + shifted within the
  RMS radius (the ab_truth_score null), depths kept; its points lie a median 2.0 (NVI), 4.9 (NW),
  5.5 (SW) km from SOME real template, so nullR is a location-specificity control, not family-free.

## Log
- 15:48Z smoke (NVI, 10 pilot locations, 2 days: 2006-09-07, 2007-06-15; `results/a4_seed/smoke/`),
  real vs nullT: Hyp2000 located 6/10 real (4 within 2.6 km of the seed: 0.4, 0.8, 1.6, 2.6 km) vs
  0/10 nullT; median stack SNR 2.74 vs 1.89 (AUC 0.93); frac of stations with S within 1 s of
  predicted 0.34 vs 0.13; detections/day 106 vs 32 -- but the null ALSO grows tens of detections per
  day at 8 x MAD, so detection counts alone are not evidence.  Code path works end to end.
- 15:53Z pilot v1 launched (raw-beam ranking).  16:07Z inspection of NW seeds: **7/10 real and
  4/10 nullR NW seeds were the SAME burst (2011-08-19 05:47:12 +-2 s)**, i.e. the raw envelope beam
  has too little lateral resolution: the brightest regional bursts top every location's list and
  the most-mutually-similar seed is shared.  Not in the USGS catalogue (no local M>~1 event), so
  probably a strong LFE/tremor burst.  Pilot v1 STOPPED by me (kill of my own process group
  1934010), partial output kept in `pilot_v1_stopped/`.
- 16:10Z tried a ring contrast (z minus the max z of 6 points 12 km away within +-2 s; smoke2 on 3 NW
  days): seeds became location-specific but the contrast (~1.3-1.7) is at noise level because the
  +-2 s pooling absorbs the uniform part of the moveout change; seed similarity fell to null level.
  Abandoned.
- 16:13Z adopted: BURST ATTRIBUTION -- beams for every location of the subarray's reference set
  (locs_all.csv, independent of the subset run), pooled +-2 s; a burst belongs to the location with
  the highest pooled z (argmax; margin kept as `contrast`).  Both seed modes run in the SAME pass
  (`seed_mode` column): `top` = raw beam ranking (v1), `attr` = attributed bursts only.  smoke3
  (NVI, 2 days, real/nullT/nullR): `top` reproduces the v1 smoke exactly (6/10 located, 4 within
  3 km; nullT 0/10); `attr` 3/10 located (2 within 3 km), nullT 0/10, nullR 1/10.
- 16:18Z pilot v2 launched (all episode days; 10 locations/subarray; real/nullT/nullR; PAR 4 x 3 threads).
