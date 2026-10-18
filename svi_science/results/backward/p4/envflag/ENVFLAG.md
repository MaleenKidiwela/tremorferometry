# P4 quiet flag, attempt 2 — envelope cross-correlation tremor days (Wech & Creager 2008 style)

Replaces the CN4 hourly-energy flag of `../P4.md` ("Quiet flag (option 3)"), which failed its own calibration
(specificity 0.28; AUC 0.51 vs PNSN quiet; the statistic tracked season). One attempt, approved by the user.
Rules below written **Thu 2026-10-08, 6:40 am PT, before any envelope was computed**.
Scripts `svi_science/scripts/bw_p4env_*.py` / `.sh`; outputs in this folder. CPU only, ≤ 8 workers.
"Scrambled control" = the same detector on the same day with each station's envelope shifted to another time of that day.

Existing code reused: `lfe_array_discovery/scripts/ab_wecc.py` (WECC for the WA Royer subarray; kappa 0.68 vs PNSN
tremor days there). Same structure; differences listed below.

## Part A — rules

### Envelopes (per station-day)
- Input: the 20 Hz caches (`cache20hz_backward` first, then `cache20hz_r1`, `cache20hz`, `lfe_catalog_repair/cache20hz`,
  `lfe_picker/cache20hz` — the order used by P5 and the energy flag). Samples where all 3 components are exactly 0 = gap.
- 2–8 Hz 4-pole Butterworth, zero phase, per component; envelope = sqrt(Σ over the 3 components of |Hilbert|²);
  mean in 1-s bins (a bin touching a gap, ± 10 s, is NaN); 2-pole 0.1 Hz zero-phase low-pass (WC2008); 1 sample/s.

### Windows and grid search
- 300-s windows, 150-s step (575 per UTC day). A station is valid in a window if all 300 samples are finite.
  Per window and station: demean, unit-norm. Full cross-correlation for every valid pair (FFT).
- Source grid: lat 47.50–49.55, lon −125.10 to −122.30, step 0.05° (2,394 nodes), fixed depth 35 km (median depth of the
  130 SVI families, 35.3 km). S travel times from the USGS CVM v1.7 `south_VI` 1-D S profile
  (`data/velocity_models/usgs_cascadia_v1.7/profiles/south_VI.csv`), direct/refracted ray by `lfe_picker/scripts/p_vel1d.py`.
  The lag bound is therefore the S moveout over the grid (pairs are read only at predicted lags, rounded to 1 s).
- Coherence(node) = mean over valid pairs of cc at the node's predicted lag. Location = argmax node; `coh` = its value.
- Per station at the located node: mean cc with the other valid stations; `s3` = the 3rd largest of these
  (the largest C at which ≥ 3 stations are coherent). Needs ≥ 3 valid stations, else the window is not usable.
- In-area: located epicentre ≤ 30 km from the nearest of the 130 SVI families (`e5r_families_2012.csv`).
- **Tremor window (threshold C)** = usable AND coh ≥ C AND s3 ≥ C AND in-area.
- Diagnostic columns only (not used by the flag): impulsiveness = max / median of the network-mean envelope in the window;
  the scrambled control (each station's day envelope circularly shifted by a distinct offset, station i by
  i × 4 h 13 min + 17 min, so no two stations stay aligned) run through the identical detector.

### Station sets
- **Main (CN4):** PGC, VGZ, NLLB, SNB only, every era.
- **Variant (all):** every roster station with a cache that day (CN4, LZB, PO ×6, C8 ×5, PB ×10, LRIV), one station per
  co-located site (PGC site = PGC > B011 > B010 > B009; LRIV site = LRIV > B005 > B006 > B007; `c_common.SITES`). Calibrated
  separately on the same 2010–11 days with the stations that existed then (no POLARIS/C8 in 2010–11, so its pre-2010
  transfer is untested; reported separately, never the main flag).

### Day statistic and flag (same structure as the energy flag, so the two compare like for like)
- N(d) = number of tremor windows that UTC day. Coverage(d) = fraction of the 575 windows that are usable;
  coverage < 0.5 → N unknown (never quiet; unknown neighbours ignored, as for Z).
- Active day: N(d) ≥ k. Env-quiet: N known AND no active day within ± 3 d.
- **Flag (main) `quiet_env3` = env-quiet AND not within ± 3 d of a Bostock window** (`truth_episodes.csv`), i.e. option 3
  with the energy statistic replaced. The env-only flag (`quiet_env`) is reported beside it.

### Calibration (once, frozen) — 2010-01-10..2011-12-31 vs the PNSN quiet flag (`bw_p4_common.pnsn_quiet`)
- Grid C ∈ {0.30, 0.32, …, 0.90}, k ∈ {1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 48}. Choose (C, k) maximising
  sensitivity + specificity of `quiet_env3`; ties → larger C, then larger k. Same class definitions and day set as P4.md
  (positive = PNSN not quiet; every PNSN-labelled day in the period counts, a day with N unknown counts as "not quiet").
  Also reported on days with N known only.
- Reported beside it, to compare with the energy flag (its values in brackets): AUC of N(d) at the chosen C for PNSN not
  quiet [0.51]; for PNSN tremor that day [0.55]; for big days (≥ 50 PNSN windows) vs PNSN quiet days [0.63]; AUC of the
  ± 3-d maximum of N (the statistic the flag actually uses) for PNSN not quiet.
- If the bars (spec ≥ 0.9 AND sens ≥ 0.8) are infeasible on the grid, the best Youden point and the best-sens point at
  spec ≥ 0.9 are reported too (diagnostic; the frozen choice is the rule above).

### Validation (once, q frozen) — 2012-01-01..2013-12-31, pass = all three (identical to P4.md)
1. specificity ≥ 0.9; 2. sensitivity ≥ 0.8;
3. control-subtracted quiet ΣP (`bw_p4_qflag.sigma_p`: frozen forward `WAVE_same` on every data and scrambled-control
   candidate of `p5f_2012`, `p5f_2013`, no P cut): env-flag rate / PNSN rate within 0.90–1.10 pooled AND within 0.90–1.10
   in ≥ 80 % of the calendar months with ≥ 5 quiet days under both flags.
- Runs once, when `P2_2012_DONE`, `P2_2013_DONE`, `P5F_2012_DONE`, `P5F_2013_DONE` all exist. No re-tuning after it.
  AUCs as above also reported on 2012–13. Variant validated the same way, reported separately.
- Fail → the env flag, the energy flag and the Bostock-only flags are all written side by side (as P4.md prescribes).

### Sanity checks (reported, not pass rules)
- N(d) on Bostock ETS window days vs days > 30 d from every window, per year (median N, fraction active), with the
  scrambled-control N beside it.
- Located tremor-window epicentres for the 2010-08-16..09-08 ETS (calibration period) against PNSN tremor epicentres of the
  same days (map + median distance to that day's PNSN median epicentre); 2012-09 added when its caches exist.

### Pre-2010
- Every year 2003–2009 is computed when `P2_<Y>_DONE` exists (nohup queue `bw_p4env_q.sh`, no polling by agents);
  the flag file is rewritten after each year with the frozen (C, k).
- Era caveat, known before computing: CN4 geometry (all 4 stations east/north of the family area) constrains E-W location
  poorly; the in-area test is therefore lenient in the west. CC is amplitude-free, so instrument epochs (VGZ, SNB) do not
  change the threshold directly, only through noise level.

## Part B — results (Thu 2026-10-08, 6:58 am PT)

**Check progress (no waiting):** `ls p4/envflag/P4ENV_*`; `tail -3 p4/envflag/logs/queue.log`; queue PID in
`p4/envflag/Q_ENV_PID`. Flag file (every day 2003–2013, rewritten after each year): `tables/quiet_flag_env.csv`
(columns `N_CN4`, `active_CN4`, `quiet_env_CN4`, `quiet_env3_CN4`, the same for `_all`, scrambled-control `Nc1_*`,
PNSN `pnsn_quiet`, the energy flag `quiet_energy_opt3`, the Bostock-only flags). Runtime ~8 min per year on 8 cores.

| step | status | marker |
|---|---|---|
| 2010–11 window tables + calibration | **done**, frozen 6:56 am PT | `tables/envflag_calibration_{CN4,all}.json` |
| 2012, 2013, then 2003 → 2009 | queued behind each `P2_<Y>_DONE` (`bw_p4env_q.sh`) | `P4ENV_<Y>_DONE` / `_FAILED` |
| 2012–13 validation (once) | runs at the first flag run with P2_2012/13 + P5F_2012/13 | `tables/envflag_validation_{CN4,all}.json`, `P4ENV_VAL_DONE` |

### Calibration 2010-01-10..2011-12-31 (721 PNSN-labelled days; 330 PNSN quiet, 391 not quiet)

| flag | C | k | sensitivity | specificity | Youden J | AUC N vs not-quiet | AUC ±3-d max N vs not-quiet | AUC N vs tremor that day | AUC big days vs quiet |
|---|---|---|---|---|---|---|---|---|---|
| energy (option 3, P4.md) | – | – | 0.81 | 0.28 | 0.09 | 0.51 | – | 0.55 | 0.63 |
| **env-CC, CN4 (main)** | 0.42 | 12 | **0.84** | **0.79** | 0.63 | 0.71 | **0.87** | **0.91** | **0.999** |
| env-CC, all stations (variant) | 0.32 | 16 | 0.74 | 0.95 | 0.69 | 0.70 | 0.87 | 0.89 | 0.996 |

- The two bars together (spec ≥ 0.9 AND sens ≥ 0.8) are **not reachable** on the grid for either set. Best sensitivity at
  spec ≥ 0.9: CN4 0.69 (C 0.56, k 12, spec 0.92); all 0.76 (C 0.34, k 12, spec 0.91). So the 2012–13 validation will fail
  bar 1 or bar 2 for the main flag as frozen (spec 0.79); it still runs once, as registered.
- What the flag misses: of the 64 PNSN-not-quiet days it calls quiet (CN4), the median PNSN count that day is 0 (max 23):
  they are days that are only "not quiet" because a weak tremor day lies within ± 3 d. Big days are separated almost
  perfectly (AUC 0.999). Of the 68 PNSN-quiet days it calls not quiet, only 12 are active
  themselves (N 12–51, median 16.5); the other 56 (median N 3) are within ± 3 d of such a day. Those 12 days are weak
  coherent envelopes in the box with no PNSN tremor — tremor below PNSN's catalogue or other coherent sources; unresolved.
- Confusion, CN4, 2010–11: PNSN quiet & flag quiet 262, PNSN quiet & flag active 68, PNSN active & flag quiet 64,
  both active 327 (energy flag: 93 / 237 / 74 / 317).
- No season: the monthly fraction of active days follows the PNSN tremor months (e.g. Apr 2010 0.00, Oct 2010 0.00,
  Aug 2010 0.71), unlike the energy statistic.

### Sanity checks
- Scrambled control: CN4 control N per day is ≤ 3 on every 2010–11 day (99 % ≤ 2), all-stations control 0 — the control
  never reaches k, so a flag-active day is never produced by chance alignment of the four envelopes.
- Bostock ETS window days vs days > 30 d from every window (`tables/sanity_ets_vs_far.csv`), CN4: 2010 median N 235 vs 3
  (active 96 % vs 13 %); 2011 168 vs 3 (77 % vs 16 %). Variant: 295 vs 2, 192 vs 2.
- 2010-08-16..09-08 ETS locations (`figures/ets_locations_2010-08.png`, `tables/ets_locations_2010-08.json`): the
  located windows reproduce the PNSN tremor migration from the SE (Olympic Peninsula / Juan de Fuca strait, day ~228)
  to the NW up Vancouver Island (day ~250). Median distance to that day's PNSN median epicentre 15 km (CN4; 80 % within
  27 km) and 11 km (all stations; 80 % within 20 km). CN4 scatters west of −124.4° (no stations there), as expected.
  2012-09 is added automatically when 2012 is computed.
