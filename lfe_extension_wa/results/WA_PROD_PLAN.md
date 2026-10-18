# WA production plan: Washington LFE families on permanent stations, data on disk first

Written 2026-10-02 (read-only planning: nothing run on GPU, downloaded, deleted or committed). For Merlin's review, then the
build agent. Paths are relative to `/home/jovyan/tremorferometry/`; `results/` = `lfe_extension_wa/results/`.
"Scrambled control" = the same detector run with each station site's data shifted to a different time of the same day,
so nothing lines up across sites except by chance. Inventory files behind section 2: `results/wa_prod_plan/`
(`inventory_by_station.csv`, `download_candidates_2012_2026.csv`, `inventory_scan.py`).

## 0. In short

- Goal: for every Washington family, (1) repeat times with a confidence column P, (2) firing rate through time, (3) raw
  −15…+30 s windows at all stations for coda dv/v. Locations are not needed. Every ≥ 8×MAD candidate is kept with its P.
- Method so far (`notes/decision_tree/royer_path.html` §ext; T1 all templates, confidence on array truth, member
  probability, T2B new-station templates, daily stacks) was measured on array-era days (2008–2011) only.
- On disk: 21 permanent stations have near-continuous raw 2014 → 2026 (2.4 TB for the roster), and ~14–22 stations
  have near-continuous raw 2007 → 2013. Most UW short-period stations and US.NLWA after 2014 are NOT on disk.
- Plan: audits first (completeness, station health, swap checks), then a production matched filter 2007 → 2026 on a
  fixed "rate network", confidence refitted to the SVI WAVE standard with a one-look sealed episode (NW 2009-05 +
  2007-01, never used so far), a small PNET gate, rates with step checks, the cut-off table, window cut, delete list.
- Cost (on-disk part): GPU ≈ 27 h (+ ≈ 9 h if PNET passes), CPU ≈ 25–35 h overlapped, wall ≈ 3–4 days. Downloads
  (separate, for approval): gap fill 0.1 TB; NW/Olympic UW 0.55 TB; SW carriers 0.61 TB; new 3-C UW ≤ 10 km 0.9 TB.
- Biggest open item: Olympic Arm A single repeats are not separable from the scrambled control (§1); expect no P ≥ 0.75
  repeats there. It can still be followed by stacking, and its new-station templates recover well (0.85).

## 1. Current state per sub-region (all numbers from the files named; held-out = V2008/V2010 test half + A2011)

| | NW (Royer) | SW (Royer) | Olympic Arm A | Olympic Arm B (Chestler-seeded) |
|---|---|---|---|---|
| catalogue v1.2 rows / representatives (v1.1 reps used in all WA runs) | 141 / 129 (132) | 81 / 79 (79) | 84 / 75 (75) | 47 / 46 (43) |
| located | 140 | 81 | 54 | 43 |
| array truth months (detections) | 38 months, 2007-01 → 2011-12 | 2007-01/02/05/07/08/12, 2008-03/05/07 | 2010-07 → 09 only | 2010-07 → 09 only |
| T1 observable families (all templates) | 132 | 73 | 28 | 41 |
| T1 recall own, all scanned truth / precision own (control chance) | 0.119 / 0.254 (0.004) | 0.010 / 0.837 (0.000) | 0.226 / 0.052 (0.001) | 0.002 / 0.224 (0.002) |
| confidence AUC own, held out (`CONFIDENCE_RESULTS.md`) | 0.85–0.87 | 0.89 (296 cands) | 0.90 | 0.82 (86 pos.; slope 1.40) |
| held-out candidates P_own ≥ 0.75: data / mean of 4 controls | 4,838 / 137 | 187 / 0 | 6 / 0 | 4 / 0 |
| same at P_nbr3 (≤ 3 km) ≥ 0.75 | 5,752 / 200 | 179 / 0.5 | 1,043 / 603 | 21 / 0 |
| T2B new-station template cc (ceiling) / pairs with ≥ 1 window | 0.54 (0.68) / 97 % | 0.33 (0.67) / 65 % | 0.85 (0.87) / 91 % | 0.42 (0.59) / 45 % |
| daily stack AUC, median per family (2008+2010) | 0.69 (best 20 in 2011: 0.84) | 0.57 | 0.56 | 0.57 |
| clear-template families on the 21 on-disk 2015+ stations (≥ 1 / ≥ 2 sites) | 108 / 45 | 14 / 0 | 21 / 1 | 7 / 1 |
| median share of a family's template weight on those 21 | 0.44 | 0.59 | 0.55 | 0.41 |

- **NW** is the strong region: single repeats separate from the scrambled control (data/control ≈ 35× at P_own ≥ 0.75,
  episode days), stacks peak at the 2012, 2013 and 2014 ETS in the 2012-06 → 2014-12 series (20 families, 12 stations).
- **SW**: few but clean repeats (precision 0.84); its step-1 probabilities failed in 2011 (control = data); its carriers
  (CPW, SMW, RRHS, TA D03D/D04D/D04E) are mostly not on disk after 2014. Only B017, B018, GNW, HDW are.
- **Olympic Arm A**: 21 unlocated families sit at the B013/B943 approximate position (a 23-member position-rule group);
  they make 6,719 of the ~9,200 candidates per day and fire on borehole noise (OLY:938: 238 vs 238 per day). Even at
  the ≤ 3-km level the control is ~58 % of the data at P ≥ 0.75. Flagged, kept (rank, do not discard).
- **Arm B**: weak and thin (precision 0.22, slope 1.40). **Others**: v1.2 has no separate Puget or other WA region; the
  Puget-area (XU 2008) families are in Royer SW/NW.
- **Not yet done** (vs the southern standard, `lfe_extension/results/SVI_REFINE_PLAN.md`, `SVI_REF_RESULTS.md`):
  WAVE model (WA confidence uses busy-hour features dens_own/dens_grp, region code and step-1 P); station health and
  swap checks; > 3 % data-error guard; station-set strength on truth days; a sealed check; rates with step checks;
  cut-off table; windows; nothing after 2014 beyond ETS trial windows.
- v1.2 moved three duplicate-group representatives from NW to Arm B (RNW:mute:2435 → OLB:576, RNW:run1:2135 →
  OLB:352, RNW:qc:2268 → OLB:298). Same families, new labels (decision D8).

## 2. Data on disk (scan of 9 raw and 7 cache roots, symlinks resolved, 2026-10-02)

Roster = the 94 stations of `results/step4b/setup/stations.csv` (42 trial + 52 step-4b). Unique station-days:

| era | raw station-days (TB) | 20 Hz cache station-days (TB) | 40 Hz cache |
|---|---|---|---|
| 2006–2011 | 32,330 (0.55) | 9,431 (0.21) | 9,782 (0.44) |
| 2012–2013 | 20,722 (0.35) | 13,676 (0.31) | – |
| 2014–2026 | 112,694 (2.39) | 71,124 (1.46) | – |
| non-roster WA stations, 2014–2026 | 21,025 (0.56; 100 stations) | 11,337 (0.24) | – |

- **Continuous 2014 → 2026 (≥ 3,900 days each, 21 stations):** CN.PGC, CN.VGZ, PB.B001, B003–B007, B009–B011, B926,
  UW.LRIV (in `lfe_extension/traces`, shared with SVI, to 2026-09-15) and PB.B013, B014, B017, B018, B943, UW.DOSE,
  GNW, HDW (real files in `lfe_array_discovery/traces`, to 2026-05-31). Gaps: B004 2020 146 d, B014 2023/2024 137/117 d.
- **Partial 2014 → 2026:** CN.CLRS (2019–26 full), CN.GOBB (2019–26, gappy), CN.PFB (to 2018-12-20), CN.SNB (to
  2021-03-27), UW.STW and UW.UWFH (2019–26). New 3-C UW stations on disk in `lfe_picker/traces`: OSQM (2021–26, 1.2 km
  from a family), RNWY (2021–, 3.9 km), TURTL (2023–, 4.8 km), ODEER (2025–, 4.3 km), BUCKS, LOPEZ, OHOH, ODUC; also CN.SYMB, PTRF.
- **Not on disk after 2014:** every other UW station (BLN, SQM, SMW, CPW, GMW, OSD, SVOH, OTR, BHW, HTW, RRHS, PCFR, VVHS,
  WISH after 2018, SP2 …), US.NLWA, PB.B019, B941. They have only the 37 T2014 days + 2–6 sampled days in 2023–2025.
- **2007 → 2013:** 14–22 permanent stations nearly every day (the array-discovery all-days pass), 39–41 in the ETS
  windows (V2008, V2010, A2011, T2012, T2013, T2014). Raw ≤ 2011 for WA stations, array stations included: ≈ 1.0 TB.
- **Caches:** WA's own 20 Hz cache (0.69 TB) covers the trial windows only. `lfe_extension/cache20hz` has 15 shared
  stations 2014–2026 (1.31 TB), but its format differs: no validity mask, 1,728,001 samples, PB t0 +8.4 ms; channel
  correlation with the WA cache 1.0000 on the 2 station-days compared (D2). The 8 WA-only continuous stations have only
  43–381 cached days after 2014.
- **Other readers:** SVI holds 2.76 TB (18 shared stations' raw + its cache) for WA (`lfe_extension/results/refine/r9/`);
  `lfe_array_discovery` owns its 0.71 TB of 2014–26 raw, `lfe_picker` its 0.21 TB.
- **Truth episodes with permanent raw on disk** (array detections of the WA families; used so far: 2008-05, 2010-08, 2011-08):

| month | NW dets | SW dets | trial stations on disk | status |
|---|---|---|---|---|
| 2007-01 | 18,750 | 18,156 | 17 | unused → proposed sealed (with 2009-05) |
| 2007-10 | 16,230 | 0 | 15 | unused → fit |
| 2009-02 / 04 | 7,314 / 5,024 | 0 | 18 | unused → fit |
| 2009-05 | 28,452 | 0 | 21 | unused → proposed sealed |
| 2010-03 | 14,568 | 0 | 19 | unused → fit |
| 2011-04 / 12 | 3,925 / 3,140 | 0 | 20 | unused → fit |

  Olympic A/B have truth only in 2010-07..09, so no independent episode exists for them. The scanned-day list of each
  month must be read from the grow runs before use (truth exists only on scanned days).

## 3. Steps on data already on disk (one shell queue, launch-and-exit, markers `results/prod/WA_PROD_{PID,HEARTBEAT,DONE,FAILED}`)

**W0 Freeze + W1 completeness** (SVI R0/R1 patterns; CPU < 30 min). Hash-snapshot existing `results/`; new outputs only
under `results/prod/`; candidate key (fid, t, series, run). Per station-day join raw → cache → MF live station; a
station-month losing > 3 % at any stage is rebuilt first (SVI R1 found 136 silently lost days). Every data-reading
stage fails if > 3 % of its station-days error.

**W2 Station health and swap checks, first** (SVI R2 + `lfe_south/results/SH_KMRFIX_RESULTS.md`; CPU 2–3 h).
- Raw, per channel per day: DC drift, constant runs, flat fraction, 2–8 and 0.1–0.5 Hz levels, line comb 2–9 Hz, rate and
  codes (reuse SVI's daily table for shared stations). Per station-month with WA families: stack cc at the detector's
  placement vs the station read at other times (scrambled control), lag ±2 s, polarity, channel permutation, rotation.
- Swap checks at every node change: PGC 2015-01-12 and 2017-08-10, VGZ 2017-03-03, SNB 2017-06-08, DOSE and GNW
  BH40→HH100 2019-05-07, LRIV 2019-05-30, HDW 2020-12-04 and HH 3-C 2022-09-01; plus B004 2020, B014 2023–24, B926 2025.
- Output `r2_health_daily.parquet`, `r2_health_monthly.csv`, `bad_for_dvv.csv` (KMR style). Must finish before any raw
  appears on a delete list (the KMR check could not be completed once the raw was gone).

**W3 20 Hz cache, WA chain** (`scripts/wa_cache.py`; CPU 3–5 h at 8–10 station-days/s). Production stations
2007-01 → 2026-09 (~20 stations × 1,800 days + 21 × 5,370 + partials, minus what exists: ~130k station-days), ≈ 2.5–3.5 TB.

**W4 Templates and nodes** (CPU < 1 h). Frozen T1 templates: 11,358 build-half stacks on the 42 trial stations +
`tpl_step4b_a` for the others. No template iteration.
- Band-only changes (same sensor, BH40 → HH100: DOSE, GNW, LRIV): carry the template over if the W2 swap check passes.
- New nodes or stations with no array-era template (HDW HH 2022, CLRS, UWFH 2019+, OSQM, RNWY, TURTL, ODEER, BUCKS,
  LOPEZ, OHOH, ODUC, SYMB, PTRF): T2B rule, i.e. stack at candidate times with P ≥ 0.3 (WAVE model, W7) from the rest
  of the network in the node's first two ETS, then cut a template. Keep every one with its cc against a
  scrambled-control stack as a column. The days used to build a node's template are in-sample for that node.
- `tpl_step4b_a` contains A2011 events. Rebuild those stations' templates from the V2010 build half only, or flag A2011
  as in-sample for them (decision D1).

**W5 Truth runs with the production station sets** (SVI D3 analogue; GPU ≈ 6.6 h). Matched filter on every truth day
(table above + V2008, V2010, A2011) restricted to each production network: S_0 2007–11 core, S_A 2012–14 core, S_B 2015 → 2019-05
(21 + PFB/SNB), S_C 2019-05 → 2022-09 (new nodes absent), S_D 2022-09 →. Four scrambled controls on truth days. The
sealed months are run but not read.

**W6 Production matched filter** (`wa_mf.py --all-templates --drop-dead 0.05`; GPU ≈ 18 h for data + 2 controls at
~9 s/day, 2007-01-01 → 2026-09-15; Olympic set to 2026-05-31).
- Rate network: a FIXED core per epoch (the 21 continuous stations, minus W2 exclusions). It changes only at the W2
  dates, so stations that appear only in ETS windows or on sampled days cannot create activity steps
  (`lfe_south/results/duc_fin/D2011.md`: catalogue rules depend on live channel count).
- Rich pass (≈ 2.5 h): days with extra stations (ETS windows 2008–14, partial and new stations). It is used for identity,
  member probability and windows, never for rates.
- Every ≥ 8×MAD candidate is stored. No tier thresholds: the per-epoch false-alarm tiers pushed WA to 18–21×MAD and
  left recall ~0.003. `cat_core`/`cat_ext` are kept as columns.

**W7 Confidence to the WAVE standard** (copy `lfe_extension/scripts/svi_ref_fit.py`; CPU ≈ 7 h features + < 1 h fit).
- Features (waveform only): net_mad, net_cc, n_sta, n_sites, live channels / n_sta_day, wsum, site_cc_max/2nd,
  site_z_max/2nd/min, n_z2/3/5, frac_z3, dom_top, hfr_rel, site_noise_mad_med, per-station z on the rate network.
  Never: dens_own/dens_grp, family or region id, step-1 P (it contains density), step-2 attribution, tremor, season, date.
- Labels: own family, ≤ 3-km group (nbr3), any family, each within ±2 s. One model per label. HGB + isotonic, GroupKFold
  by day. Ablations WAVE / FULL / NO_PNET / MF_ONLY.
- Day roles, fixed before fitting: in-sample (never fitted or scored) = V2008/V2010 build halves (and A2011 if D1 says
  so); fit out-of-fold = V2008/V2010 test halves, A2011, 2007-10, 2009-02/04, 2010-03, 2011-04/12; **sealed one look =
  2009-05 + 2007-01** (`SEALED_OPENED.json` written first, refuses to run twice).
- Bars, reported and not tuned: AUC vs net_mad; slope 0.8–1.2; matched/ΣP 0.8–1.2; control P > 0.5 ≤ 1 %. Expect the last
  to fail on tremor days: there the scrambled control also lands on real LFEs.
- Calibration in the large is checked per epoch (A2011 matched/ΣP was 0.76: base-rate shift). P = probability of matching
  the array catalogue, a lower bound on "real"; between ETS windows it is not a measure of real.

**W8 PNET gate (small)** (GPU ≈ 0.5 h, truth days only). Use the frozen WA fine-tune (`results/pnet/models/`), no new
training. Adopt only if WITH_PNET − WAVE ΔAUC(own) ≥ 0.01 with its interval above 0 on the sealed months, AND ≥ 0 at the
production stations PNET never trained on (B017, B018, GNW, HDW); judged in the same single sealed look. Earlier WA gain
was station-specific (0.870 vs 0.778 at held-out stations). If adopted: production CF ≈ 9 GPU-h, CF files temporary.

**W9 Member probability** (`scripts/step2_attrib.py`, CPU < 1 h). P(member) for every candidate in a multi-member group
(65 groups, 179 families), as a column. It picks a member just as confidently on the scrambled control (0.77 vs 0.76),
so it says which member only if the repeat is real.

**W10 Cut-off table, before fixing any region's cut-off** (`lfe_south/scripts/cutoff_check.py` pattern). Per region and
label: P 0.5 / 0.6 / 0.7 / 0.75 / 0.8 / 0.9 × days with / without PNSN tremor within 30 and 60 km of the family
(`catalogs/pnsn_tremor_cascadia_full.csv`, from 2010-01-06), data vs scrambled control, false share, families covered.
Output `results/prod/cutoff_check/false_share_by_cutoff.csv`. The default stays P ≥ 0.75 for single repeats and rates.

**W11 Rates with step checks** (SVI R7; CPU < 1 h).
- Per family per day: ΣP and count at the region cut-off for own / nbr3 / any, with fold spread; scrambled-control ΣP
  and count subtracted and shown; k = truth/ΣP per family and station set from W5.
- Step checks (±21 d and ±45 d against the same calendar windows in other years) at every W2 date and health flag.
  A step next to an SSE cannot be separated from activity (SH_KMRFIX §4).
- Daily-stack activity A (`stack_activity.py`) on the rate network as a separate column for weak families. Inter-ETS
  rates use the SVI validation (withheld-station S energy, control ΣP over time, station-set invariance); if it fails,
  report them as a control-subtracted bracket.

**W12 Windows for dv/v** (`lfe_south/scripts/cut_windows.py`, new profile `wa_ext`; CPU ≈ 9 h).
−15…+30 s around the expected S (template offset; else CVM v1.7 moveout from the location or approximate position; else
0, recorded) at every WA station with a raw day file that day; raw counts, native rate, all channels, HDF5 per month +
parquet index joined to P. Which candidates: decision D5; a pilot month gives the real size first. Verify 20 windows
sample-identical vs the miniSEED; stack check per region vs template and scrambled control. Truth episodes included
(ties the dv/v series back to 2007).

**W13 Delete-candidate list, joint with SVI** (`svi_ref_delete.py` pattern, needs-check incl. symlink readers).
Nothing is deleted. The list goes to the user only after W2 + W12 are verified:
raw 2014–26 of the 18 SVI-shared stations (≈ 1.4 TB, now held by SVI); SVI cache (1.31 TB) and WA production caches
(2.5–3.5 TB), both regenerable; 2012–13 raw (0.35 TB); `lfe_array_discovery` / `lfe_picker` raw (0.71 / 0.21 TB) listed
only, as those projects' decision. Keep: truth raw 2006–2011 (≈ 1.0 TB), the sealed months, all results, cut store, logs.

## 4. Needs downloads (separate; costed for the user's approval; nothing here runs without it)

| tier | what | station-days | raw | why |
|---|---|---|---|---|
| A | gaps of the 21 core stations 2012–26 (incl. the 8 Olympic-set stations 2026-06 → 09) | 5,231 | 0.11 TB | complete the rate network |
| B | NW/Olympic UW + US 2012–26 while live: BLN (→ 2016-07), SQM (→ 2019-09), NLWA (→ 2018-06), OSD, SVOH, GMW, OTR, BHW, HTW, STW 2012–18, SP2, WISH, STOR | 53,268 | 0.55 TB | BLN 45 and NLWA 11 clear NW templates |
| C | SW carriers: CPW, RRHS, PCFR, SMW (→ 2021-12), VVHS, PB.B019, B941, RER | 39,679 | 0.61 TB | SW has 0 families with 2 clear sites on disk |
| D | new 3-C UW stations from 2017 within 10 km of a family (22: EQUIL, NATEM, CHIMA, KNGST, DOORS, TKEY, PAN4H, DOTY, …); 20 km: 43 stations | 32,585 (67,234) | 0.90 (1.73) TB | T2B templates, more sites after 2019 |

- Rate 1–2.3 station-days/s (measured in lfe_array_discovery): A–C (98k station-days, 1.26 TB) ≈ 12–27 h; D 4–9 h more.
- Streaming, month by month, oldest first: fetch month M → > 3 % guard → cache → W2 raw health → MF (B/C stations join
  the rate network as new epochs, step checks at their start/end) → features → P → windows → verify 20/20 → month M's
  exact delete list (downloaded raw + cache only). Raw buffer cap ~300 GB; the queue pauses when full until the user
  approves the listed months. With 66 TB free, A–C also fit with no deletion and one list at the end.
- Re-fit: P for the downloaded epochs needs those stations in W5 truth runs; they are on disk for 2007–2014 truth days.

## 5. Compute, time, storage (on-disk part)

| stage | resource | estimate | basis |
|---|---|---|---|
| W1 | CPU | < 30 min | joins |
| W2 | CPU 24 workers | 2–3 h | ~46k WA-only station-days; SVI 69.8k in 2–4 h |
| W3 | CPU 12–16 workers | 3–5 h; +2.5–3.5 TB | WA cache logs: 8–10 station-days/s |
| W4 | CPU | < 1 h | T2B 227 s for 1,003 pairs |
| W5 | GPU | ≈ 6.6 h | ~330 truth days × 4 sets × 18 s (4 controls) |
| W6 | GPU | ≈ 18 h core + 2.5 h rich | mt1_all: 97 days in 14.6 min with 2 controls |
| W7 | CPU | ≈ 7 h + < 1 h fit; features 20–40 GB | step-1 features ~1.2 s per day-series |
| W8 | GPU | 0.5 h gate (+ ≈ 9 h if adopted) | SH_CONF 73 min per 21k station-days |
| W9–W11 | CPU | < 3 h | |
| W12 | CPU 6–12 workers | ≈ 9 h; store per D5 | Ducellier 708 GB in 86 min |

- Totals: GPU ≈ 27 h (≈ 36 h with PNET); CPU stages overlap; wall ≈ 3–4 days. At most 3 matched-filter passes on the
  H100 at once; GPU stages wait for any other GPU queue lock.
- Candidate volume: ~9,200 data candidates/day today (Olympic A 6,719, NW 1,510, Arm B 235, SW 18). That is ~65 M over
  2007–26, ~200 M rows with controls. Expect a 30–60 GB parquet.
- Storage peak: caches 2.5–3.5 TB + cut store 0.5–2.5 TB (D5) + features ≈ 3–6 TB of the 66 TB free.

## 6. Open decisions for the user (each with a recommendation)

- **D1 Sealed check.** Recommend sealing NW 2009-05 (28,452 detections, 21 stations) and all of 2007-01 (NW 18,750 + SW
  18,156). Fit out-of-fold on the other episodes, A2011 included. Olympic A/B get no sealed check (2010 is their only
  truth), so they are reported out-of-fold only. Also rebuild the step-4b templates without A2011 events so that A2011
  stays leak-free. Alternative: put SW 2007-01 into the fit, since SW has only 296 test-half candidates.
- **D2 Cache.** Recommend rebuilding all production caches on the WA chain. The SVI cache differs (no mask, PB t0
  +8.4 ms), and an adapter would need its own equivalence test. Cost ≈ 1.3 TB of duplicate storage, ~2 h CPU.
- **D3/D4.** Recommend keeping caches to the end (regenerable, then on the delete list), and a fixed rate network per
  epoch for rates with extra stations in a rich pass for identity and windows only.
- **D5 Which candidates get windows.** The −15…+30 s union of every candidate covers ~57 % of each day, so the store
  would be about half the raw. Measured union (FULL model, max over labels): P ≥ 0.1 → 25 %, ≥ 0.2 → 18 %,
  ≥ 0.3 → 13 %, ≥ 0.5 → 8 %. Recommend P_max ≥ 0.2 under WAVE, plus every candidate on truth days, re-measured on
  the pilot. Every candidate keeps its P row either way. Risk: a later refit can raise a P above 0.2, and that
  window would then be missing.
- **D6 PNET.** Recommend the W8 gate; no production PNET unless it passes.
- **D7 Labels per region.** Recommend: produce the W10 table first. Tentatively P_own ≥ 0.75 for NW and SW;
  P_nbr3 + member probability for Arm B. Olympic A (incl. the 21 unlocated families at B013/B943, kept and flagged) is
  followed by stacks and template recovery only; no rate is claimed unless its quiet-day false share is < 10 %.
- **D8 Representatives.** Recommend detecting with the 329 v1.1-representative templates, which exist and were
  validated, and reporting under v1.2 uid + dup_group. No re-templating for the three moved representatives.
- **D9 Array-era production 2007–2011** (≈ 5 GPU-h of W6). Recommend yes: rates and dv/v windows from 2007, and k
  factors from many episodes.
- **D10 LRIV after 2019-05-30.** Recommend: follow SVI (out of the rate network), keep it in features and windows
  unless the W2 WA-family check puts it above the control.
- **D11 Downloads.** Recommend approving tier A after W1–W2. Decide on B and C after the on-disk W10 table shows how
  many NW/SW families reach P ≥ 0.75 without them. Hold D until T2B on the on-disk new stations (OSQM, RNWY, TURTL,
  ODEER) shows the gain.

Not in scope here: relocation, template iteration (rejected on held-out data in SVI), subspace detectors, slip fronts.
