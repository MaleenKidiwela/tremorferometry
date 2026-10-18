# S13 NOTES — dated decision log (times PT)

## 2026-10-13 (Tue)

### ~10:20 am: inventory (`s13_inventory.py`; 2,424 root × station × year rows)
- Every trace root holds RAW-COUNT continuous day files (int32 STEIM2, native 40 Hz BH or 100 Hz HH), not decimated
  or bandpassed. The 20–50 s band is intact; only the instrument response limits it.
- The 20 Hz caches (`cache20hz*`) are derived and are not used here.
- Excluded: PB boreholes (2 Hz geophones) and EH short-period channels.
- Broadband days on disk (merged roots, `tables/inventory_bb_near.csv`):

  | station | days | years |
  |---|---|---|
  | VGZ | 8,293 | 2003–26 |
  | PGC | 8,143 | 2003–26 |
  | NLLB | 7,950 | 2003–26 |
  | GNW | 6,866 | 2006–26 |
  | SNB | 6,263 | 2003–21 |
  | LRIV | 5,870 | 2009–26 |
  | DOSE | 5,701 | 2010–26 |
  | WISH | 5,344 | 2007–26 |
  | SP2 | 5,136 | 2008–26 |
  | PFB | 2,920 | 2009–18 |
  | LZB | 2,033 | 2003–14 |

  Also the 2003–06 POLARIS C8 stations, the XU CAFE array 2006–10, TA 2006–16, and about 20 UW Cascadia-compact HH
  stations from 2021 on.
- Download estimate:
  - CN has no LH archived at EarthScope.
  - The IRIS timeseries service (server-side decimation) is retired ("This service has been retired"). Any download
    is therefore raw BH/HH: about 7 MB/day at 40 Hz and 17 MB/day at 100 Hz, about 1 MB/day after our 1-sps prep.
  - Filling the SVI episode windows (±2 d) for the stations with gaps (LZB 2015+, SNB 2022+, PFB pre-2009/2019+,
    YOUB, GOBB 2017+) is about 2,500 station-days, roughly 35–45 GB raw (about 3 GB after prep). That is over the
    ~20 GB OK → NOT fetched.
  - Not needed for the first pass: the long-record near stations are already on disk.
  - Fetched only the 2012-08 pilot gap: LZB and YOUB, 74 day files, 1.2 GB. PFB 2012 has no data at EarthScope
    (channel epoch open, 404 for every day).

### ~10:30 am: expected size (`s13_green.py`)
- A 1e12 N m thrust at 35 km gives 0.4–1.2 nm peak at 20–50 s within 40 km.
- Measured band noise per window: PGC Z 25 nm, N 58, E 109 nm; VGZ Z 13 nm; LRIV Z 16 nm. That is far above the NLNM:
  the coastal microseism tail plus tilt on the horizontals.
- If the LFEs carried only their own Mw 1.7 moment, the pooled episode stack would reach about 1–3σ.
- The pilot found ~12× more (see RESULTS).
- Green's function checks:
  - an explosion gives an outward static;
  - a thrust gives uplift above the source with SW (trenchward) horizontal motion, as for an ETS.

### ~10:40 am: pilot
- A test run with 20 draws (`tables/test_*`) was looked at before the 500-draw run. No rule was changed after it.
- 2012-08 PGC lost 99 % of its samples: gaps every few minutes, and the prep did not implement the PLAN's "interpolate
  ≤ 5 s" step. Recorded as Amendment 2. The pilot verdict stands as run.
- Verdict: works (P1, P3 both episodes; P2 2011-08; P4 4/5 families). See RESULTS.md.
- No Merlin review (Fable quota rule). Recommended before the scale-up results are trusted.

### ~10:55 am: scale-up launched (Amendments 1–2)
- 39 SVI majors (13 of the 52 have < 30 LFEs at P ≥ 0.75, mostly 2014+ episodes, whose forward catalogue is thin at
  this cut).
- 56 WA-NW majors (13 below 30).
- Between-ETS sets.
- `s13_stack2.py` reproduces the pilot pooled 2011-08 moment exactly (5.33e12) on the same stations: code check
  `tables/test_s2_2011-08_*` (20 draws).
- ~11:15 am PT: gap fill (Amendment 3) launched: 1,884 station-days ≈ 28 GB raw, then scalef_* restacks after the scale-up; log logs/fill_queue.log.

## 2026-10-14 (Wed)

### ~5:45–6:15 pm PT: scale-up write-up (before gap fill)
- Summary covers all 40 SVI and 57 WA stacks; no NaN rows. No `scalef_*` table existed yet, so these are pre-fill.
- **Broken 1: flat channels.** UW.LRIV N/Z have control SD exactly 0 in 2019-04 (and partly 2019-05); UW.OPC is flat
  in WA 2009-04.
  - The zero SD gives infinite weight, so the moment collapses to ~0. 27 of 1,623 units are affected.
  - The as-run A2 "weighted mean" (5.9 N m SVI, 707 N m WA) is exactly the 2019-04 value. The WA 2019-04
    "detection" (z 3.9) is the same artefact.
  - Fix needed before the fill summary: drop channels with σ ≤ a floor (e.g. < 1 % of the station's median), or
    drop units with control SD < 1e10 N m.
  - The write-up recomputes A2/A4 without these units (`scripts/s13_writeup_checks.py` →
    `tables/writeup_bgf_*_clean.csv`).
- **Broken 2: UW.SP2 E is near-dead 2012-04 → 2019-03.** σ is 0.04–0.35 nm vs 7–20 nm on SP2 N, and it holds 33–70 %
  of the WA fit weight in 19 episodes and both WA between sets.
  - Re-fit without it from the stored stacks (`tables/writeup_bgf_refit.csv`; the all-channel re-fit reproduces the
    stored moments exactly) raises those moments 1.2–2.6×.
  - No new z: control draws are not stored.
- **Minor:** the scale-up `n_nb10` counts each LFE's own copies in its half/zone/family units (it reads 2–4, not the
  pilot's 0.16–0.19), so it is not meaningful.
- Result in RESULTS.md "Scale-up (before gap fill)":
  - SVI: 20/38 detected, 38/38 positive. Constant 5.3e12 N m per LFE 2008–26 (χ²/dof 1.03); 2003–05 lower.
  - SVI between-ETS LFEs carry 0.79× the ETS value (z 3.7).
  - SVI total vs GNSS ρ 0.81, carried by the LFE count.
  - WA: 13/55 detected, with the signal in zone 1 (SE) only. Zones 2–3 slightly negative (suspect template/location
    there). WA ±1-day control width 1.3–1.5 (SVI 0.9).
- New figure: `figures/s13_scale_both_m0_per_lfe.png` (`scripts/s13_scale_fig_both.py`).
- No Merlin review of the scale-up yet.

### ~5:50–6:00 pm PT: channel guard (PLAN.md Amendment 4)
- Rule written into PLAN.md at ~5:52 pm PT, before any guarded number: per unit and component, drop a station-component
  with control σ ≤ 1e-3 × that component's median σ over the unit's stations (or σ not finite); then cap any single
  station-component at 0.3 of the total weight (water-filling; < 4 live components → equal weights). Same weights for
  data, scrambled-control draws, ±1-day, P-weighted and free-MT fits.
- Code: `s13_common.guard_weights` / `guard_violation` (new), used in `s13_stack2.py` (stacomp gains `w_share`,
  `guard_dropped`; units gain `guard_dropped`, `guard_capped`; `n_sta` = live stations). Installed with os.replace at
  5:53 pm PT. The fill queue was still in prep then (no scalef stack had started), so every `scalef_*` stack uses
  the guarded code. Pre-guard copy: `scripts/.s13_stack2.py.pre_guard_2026-10-14`.
- Test (20 draws, SVI 2019-04, output moved to `../_quarantine_2026-10-14_envfail/s13_pre_guard/`): LRIV N/Z dropped;
  pooled moment 3.5e12 instead of ~0.
- The rule bites much more widely than the two known cases: in the first-pass tables **847 of 1,623 units** break it,
  in 56 of 97 row files (`tables/guard_violations_scale.csv`). Flat channels: LRIV N/Z (SVI/WA 2019-04, WA 2019-05,
  SVI 2023-04 Z), OPC N/E (WA 2009-04). Over the cap: UW.SP2 E (WA 2012–2019, up to 98 %), but also the best quiet
  channels, CN.PGC Z (up to 88 %, SVI 2003-02), CN.NLLB Z, XU.BS11 Z, CN.VGZ Z/E. So the cap also down-weights good
  quiet stations; expect somewhat lower z in SVI, as the price of not resting a unit on one channel.
- Re-runs: 29 of the 56 row files are re-stacked by the fill queue (`scalef_*`, guarded). The other 27 are re-run as
  `refit_*` (`scripts/s13_refit_queue.sh`, list `tables/refit_list.txt`, 4 parallel, OMP 1, same station lists and
  data dir as the first pass; their days touch no fetched day, so inputs are identical apart from the guard).
  Started 5:54 pm PT; svi_2004-07 is the long one (~1 h).
- `s13_scale_summary.py`: table choice scalef_ > refit_ > scale_; a unit from an unguarded table that breaks the rule is
  set to NaN and listed in `tables/scale_guard_excluded.csv` (count in `scale_verdict.json` → `guard_excluded_units`).

### ~6:15–6:30 pm PT: guard v1 (cap) withdrawn → v2 (Merlin); WA withdrawn
- Coordinator relayed Merlin's review: (1) replace the 0.3 cap with a defect-only rule (v2, PLAN Amendment 4 rewritten);
  (2) **WA is broken at the root (Merlin F1): `s13_units.wa_shift` moves every WA origin ~16 s early** (the WA
  catalogue `time` is already origin-like; `s_off` is a within-template offset). WA numbers are withdrawn; WA will be
  re-run from scratch after the origin fix. No WA refits.
- Stopped my v1 refit queue at ~6:18 pm PT (11 v1 refits had finished, 4 killed). All v1 outputs (`refit_*`, and the
  fill queue's SVI `scalef_*`, which ran 6:00–6:14 pm PT with v1 code) moved to
  `../_quarantine_2026-10-14_envfail/s13_pre_guard/v1_cap_rule/` (212 files). The fill queue keeps running its WA
  `scalef_*` (v1 code before 6:22 pm PT, v2 after); they will be superseded by the WA origin fix.
- v2 rule (`s13_common.guard2`): drop if n_used < 5, σ ≤ 0, σ < 0.05 × median σ of the station's other two
  components, or σ < 0.02 × median σ over all channels in the unit; no cap. `s13_stack2.py` now also writes
  `max_share` (largest (tpl/σ)² share), leave-one-station-out M0 range, drop reasons, and stores the weights in the pkl.
  Installed 6:22 pm PT.
- Tests (20 draws; outputs in `.../s13_pre_guard/tests/`):
  - SVI 2019-04 pooled 3.35e12 (was ~0): LRIV N/Z dropped (σ = 0). WA 2009-04 pooled 1.56e12, z 5.9 (OPC N/E dropped,
    n_used < 5 and σ = 0). WA 2019-04 pooled 3.85e12 (LRIV N/Z and SP2 E dropped). Pass.
  - Reproduction: on the same data the v2 stack of the first unit is bit-identical to the first pass and the
    all-channel refit gives the stored M0 to 2e-16 (later units differ only through the draw count in the test).
- **Side effect to flag:** the station-ratio clause (σ < 0.05 × the station's other two) also drops healthy quiet
  verticals whose horizontals are tilt-noisy at 20–50 s: CN.SNB Z (σ ≈ 0.3–6 nm vs N/E 10–440 nm) in 266 SVI units,
  C8.KHVB Z 177, CN.VGZ Z 101, C8.SHDB Z 33 (first-pass tables, `tables/guard2_drops_scale_svi.csv`; 641 drops in 503
  of 804 SVI units). These are not dead channels; dropping them costs sensitivity but should not bias M0
  (SVI 2015-12 pooled 5.42e12 vs 5.62e12 with SNB Z). Applied as written; not changed after seeing it.
- SVI re-runs (`scripts/s13_refit_queue.sh`, v2, 4 parallel, started 6:23 pm PT): 10 non-fill row files with a v2 drop
  → `refit_*` (2003-02, 2004-07, 2012-08, 2013-09, 2014-11, 2015-12, 2017-03, 2020-05, 2020-11, 2023-04); 29
  fill-affected row files → `scalef_*` on the filled data. svi_2005-09 has no v2 drop, so its first-pass table stands.
- `s13_scale_summary.py`: SVI only (WA marked withdrawn in `scale_verdict.json`); v2 tables used as they are; other
  tables only for units v2 would not change.

### ~6:25–6:40 pm PT: Merlin steps 2–6 + Q1c launched (PLAN.md Amendment 5)
- Amendment 5 written ~6:30 pm PT before any number. New code only: `scripts/s13_units_v2.py` (→ `data/rows_v2/`),
  `scripts/s13_stack2_v2.py`, `scripts/s13_v2_queue.sh`, `scripts/s13_v2_summary.py`. `data/rows/`, `s13_units.py`,
  `s13_stack2.py` untouched (fill and refit queues read them).
- **Check of the WA origin before running.** WA template offsets are not travel times: s_off vs hypocentral distance has
  slope 0.005 s/km over 274 families (all within −5…+20 s), so `wa_shift` (median s_off − r/Vs ≈ −16.6 s) was wrong, as
  Merlin said. But the 417 "clear" template arrivals put S at `time` + 2.7 s at 20–40 km, where S needs ~8–10 s:
  the physical origin is ≈ 5–8 s *before* `time`. So `time + 0.5` (the prescribed origin, used as the headline) is
  probably ~6 s late; the best lag should land near −5…−7 s. Recorded in the amendment as a prior; extra lag-grid
  columns `c6_` (−14…+2 s) and `wide_` (±40 s) carry their own scrambled-control z.
- **v2 stacker check:** with `--kmax 1` on `svi_2025-09` it reproduces the current (guard v2) `s13_stack2.py` on every
  unit: M0 to 3e-16 relative, z to 2e-15, same lags and drops. `--kmax 5` path tested on WA 2022-08 (runs; 78 s).
  Test outputs in `../_quarantine_2026-10-14_envfail/s13_v2_test/`.
- **WA rows v2:** 37,061 selected NW rows → 28,537 after the ±2 s cross-family dedupe (77 %; post_v2 variant-C flags
  1,765 of the input). 32 of 132 families had depth < 20 km or missing → Slab2. 54 of 69 majors built (first pass 56).
  WA_BE 220 → 220, WA_BE_broad 934 → 918 after dedupe. Prep days ±6 d (for the ±5-day controls): 2,249 days,
  11,839 new station-days to prep (no download; raw on disk).
- **SVI step 4 rows:** 2004-07 bins 5,715 / 7,130 / 2,847 (P .75–.90 / .90–.97 / .97–1), top-3,000, random-3,000;
  2005-09 4,572 / 4,363 / 1,471; 2012-08 1,504 / 774 / 586; 2013-09 771 / 351 / 237.
- **SVI step 5 rows:** BE_broad re-built from `s7/between_broad_events.parquet` with P_same (origins identical to the
  first pass, max diff 0 s). Top-5 days: 2020-09-13 (77), 2015-10-08 (64), 2015-10-10 (42), 2014-03-23 (33),
  2020-09-14 (32) = 248 of 688 (36 %). Units: all 688; minus 2020-09-13/14 579; minus 2015-10-08/10 573;
  minus top-5 440; tau_next ≥ 30 d 488; 2014–19 357; 2020–26 331.
- **Queue** `scripts/s13_v2_queue.sh` started 6:29 pm PT (PID in `logs/v2_queue.pid`, log `logs/v2_queue.log`):
  A) `pbin_svi_*` ×4 + `svibe2` (≤ 2 parallel, kmax 1); B) prep WA days (12 workers, ~1 h); C) `wa2_*` 54 episodes +
  `wa2_between` (kmax 5; ≤ 4 parallel while fill/refit queues run, ≤ 6 after). Done = `V2_QUEUE_DONE`.
- Early read of the PGC-out refit on the SVI tables as they stood at 6:33 pm PT (refit/scalef v2 still running, so
  not final): 2008–2026 weighted mean 5.39e12 all stations vs 5.67e12 PGC-out (+5 %). The white-noise σ formula
  under-states the control SD by ×6.2 (median), so single-channel and PGC-out σ use that measured ratio.
- Fixed n_nb (pooled rows, self excluded) for SVI from the row files: n_nb10 median 0.06; 2004-07 0.38, 2005-09 0.29,
  2003-02 0.19 (matches Merlin). Kernel-weighted sum is small and negative (the 20–50 s kernel's side lobe at
  ±10–20 s), so catalogued neighbours do not inflate the moment.
- **Resume:** when `V2_QUEUE_DONE` exists and the refit queue has finished (`REFIT_QUEUE_DONE` or no
  `s13_refit_queue.sh` process), run `python scripts/s13_v2_summary.py` (writes `tables/v2_verdict.json`,
  `tables/v2_*.csv`, `figures/s13_v2_wa_episodes.png`, `figures/s13_v2_svi_pbin.png`), then score W1–W5 / step 4 /
  step 5 and write RESULTS.md. If the queue died: re-run the same nohup command (idempotent, DONE markers).

### ~6:30–7:45 pm PT: SVI v2 re-runs done, summary, RESULTS rewritten
- Refit queue (v2) ran 6:22–7:23 pm PT, no failures: 10 `refit_svi_*` + 29 `scalef_svi_*` (gap-filled data, fill station
  list with GOBB). Marker `REFIT_QUEUE_DONE`. The fill queue ended 6:45 pm PT with "13 scalef stacks missing" (the SVI
  ones I had moved to quarantine and was re-running). Once all 54 `STACK_scalef_*_DONE` existed and the summary ran
  cleanly, I touched `FILL_QUEUE_DONE` by hand (same condition as the queue; logged in `logs/fill_queue.log`).
- Bug of mine, fixed: the first summary patch cut out `perm_spearman`, so the fill queue's own 6:45 pm summary run
  failed (`logs/scale_summary_fill.out`). Restored; final run `logs/scale_summary_g2.out`.
- Reproduction pass: in the 16 refit units with no v2 drop, M0 and z equal the first pass exactly.
- Summary picks 28 `scalef_` / 10 `refit_` / 1 `scale_` (2005-09, no v2 drop) episodes. 0 units excluded.
- SVI corrected numbers:
  - 21/39 detected, 0 at z ≤ −3, 38/39 positive; 2019-04 is now 3.4e12, z 1.2.
  - A2 all 39: 4.4e12 (χ²/dof 2.7, p 9e-8). **2008–26: 5.4e12, χ²/dof 1.09, p 0.34.**
  - BE_broad: 4.4e12, z 3.2 (0.81× ETS).
  - Total vs GNSS ρ 0.81 (n 15); families 29/111 at z ≥ 3.
  - ±1-day SD 1.00 / 0.94; odd/even ρ 0.70.
- Check column `scripts/s13_guard_variant.py` → `tables/guard2_zvariant.csv`: re-fit with the healthy verticals that v2
  drops added back. Median ratio 1.00 (40 pooled/between units, 33 within 10 %). It differs in tiny episodes:
  2012-04 pooled −2.7e12 → +1.0e12 (116 LFEs, z −0.6).
- `max_share` median 0.55; 2003-02 rests on PGC Z (0.89; leaving PGC out gives −1.3e13).
- **For the WA re-run (not done now):** v2 drops SP2 E in only 240 of the 315 WA first-pass units where SP2 E is
  near-dead (E/N < 0.1). In the other 75 (mostly 2012-09, 2018-05, 2019-03, 2016-08), E/N is 0.03–0.065 but E/Z is
  0.10–0.27. The median of the other two is then above 20 × σ_E. So the Merlin pass check "SP2 E dropped 2012–2019"
  is only partly met. Needs a decision before the WA re-run.
- RESULTS.md scale-up section rewritten as "SVI (after gap fill, channel guard v2)"; WA section withdrawn (old text in
  `../_quarantine_2026-10-14_envfail/s13_pre_guard/RESULTS_pre_guard_2026-10-14.md`).
- Stale on disk: `figures/s13_scale_both_m0_per_lfe.png`, `figures/s13_scale_wa_*.png`, `tables/scale_wa_*` (summary
  no longer writes WA), `tables/writeup_bgf_*` (pre-guard).

### ~7:15–7:40 pm PT: UW.SP2 E instrument veto (coordinator decision; PLAN Amendment 5 addendum)
- Data-driven span 2011-09-05 → 2019-05-22 inclusive (daily 20–50 s E/N RMS: ≥ 0.7 before, 0.38 on 09-05,
  0.004–0.04 after; normal from 2019-05-23). Row added to `svi_science/instrument_vetoes.csv`.
- `s13_stack2_v2.py` patched (blank vetoed component-days before the band-pass; logged), os.replace 7:29:58 pm PT;
  summary refits honour it too. wa2_between and wa2_2012-02 (started on old code) killed 7:30 pm PT (queue log
  shows them as FAILED); `scripts/s13_v2_rerun.sh` (nohup, PID in ps) waits for `s13_v2_queue.sh` to end, re-runs them,
  then writes `V2_QUEUE_DONE`. The 8 stacks done/started before the patch (2007-01 … 2011-07) touch no vetoed day.
- Side effect: s11/s14/s12 veto readers apply rows station-wide, so UW.SP2 is vetoed there too in that window.
- Note for the write-up (no action): guard2's 0.05 station-ratio clause drops healthy verticals when horizontals
  are tilt-noisy (SNB/KHVB/VGZ Z); M0 unchanged (`tables/guard2_zvariant.csv`) — mention in RESULTS.
- SVI step-4 (pbin_*) and step-5 (svibe2) stacks finished 6:37–6:53 pm PT (not yet summarised). WA prep done 7:16 pm
  PT (8,425 new station-days).

### ~8:10–8:15 pm PT: Amendment 5 summary + RESULTS
- `V2_QUEUE_DONE` 8:08 pm PT (wa2_2012-02, wa2_between re-run with the SP2 veto; 27 SP2 veto log lines). Summary bug
  of mine fixed (itertuples renames `zk-5` columns → indexed access); `logs/v2_summary.out`.
- **WA W1 FAIL** (headline ±8 s about time + 0.5): pooled lags bimodal, 30/54 at −4…−8 s with M0 > 0, 13/54 pinned at
  +8 s all with M0 < 0 (opposite side lobe); mode +8, 33 % at ±8. Failure branch: grand-stack ±40 s scan peaks −7 s
  (zone 1, 4.0e12) / −5 s (pooled, 2.6e12) → residual origin offset −5…−7 s, as the pre-run near-station prior said.
  WA stays withdrawn as headline. W2–W5 pass (nb2 0.00, ±k-day SD 1.14, zones 2–3 positive 20/24 18/24, ρ 0.45).
- WA column c6 (−14…+2 s, own controls): 21/54 z ≥ 3, 0 ≤ −3, 51/54 positive, 2.7e12 per LFE (χ²/dof 1.5), all
  zones positive (4.1 / 2.0 / 3.3e12), odd/even 0.57, total vs tremor ρ 0.84, BE_broad 2.4e12 z 4.2. Decision for
  coordinator/Merlin: adopt origin = time − 5.5 s (stacks exist; summary pass only). New figure
  `figures/s13_v2_wa_c6.png` (`scripts/s13_v2_wa_c6_fig.py`), grand scans `tables/v2_wa_*_grand_lagscan.csv`.
- WA after dedupe still has 0.81 neighbours per LFE within ±10 s (kernel-weighted +0.20): cross-family near-duplicates
  2–10 s apart likely remain.
- SVI step 4: pre-set outcome (c) (high/low P ratio 2.4 in 2004-07, 2.3 in 2012-08); (a) fails (top-P 2004-07 3.0–3.2e12);
  (b) misses only on 3.03e12 > 3e12. At fixed P 2003–05 still lower (high bin 3.0/4.4 vs 7.1/8.8e12); 2004-07 is not
  low-P-heavy. Single-channel Z: the 2003–05 drop rests mainly on PGC Z; TWBB 5.0 ± 1.9 and VGZ 4.0 ± 2.6 in 2004-07.
- SVI step 5 PASS: all variants positive, minus-top-5 4.4e12 z 2.7, halves 3.4e12 / 4.0e12; ratio 0.80 ± 0.25;
  tau_next ≥ 30 d weaker (2.6e12, z 1.6).
- PGC-out: 2008–26 +5 % (5.65 vs 5.39e12), all +9 %; single PGC-carried episodes flip (2003-02, 2008-05).
- RESULTS.md: new section "Amendment 5 results" (SVI provisional statement per Merlin step 7 incl. guard2 vertical-drop
  note; steps 4–6; WA corrected run with rule table); earlier SVI plain summary re-worded (pre-set "changes" first,
  2008–26 constancy post hoc). Pre-edit copy `../_quarantine_2026-10-14_envfail/RESULTS_before_amendment5_2026-10-14_2030PT.md`.

### ~8:45–9:05 pm PT: Amendment 6 (Merlin WA-origin review) — steps 1–3 done, stacks launched
- PLAN.md Amendment 6 written ~8:50 pm PT before any offset/row/stack number; Addendum 6a ~9:05 pm PT after the first
  catalogue-only step-1 pass, before any row or stack (reasons in the addendum).
- **Step 1 (catalogue only, `scripts/s13_offsets_v3.py` → `tables/wa_family_offsets.csv`, `wa3_offsets_verdict.json`).**
  - Every co-detected family pair is an exact template-alignment constant: selected rows, 325 ordered pairs with
    n ≥ 8, 302 of them IQR ≤ 0.1 s; network residual rms 0.011 s. fam38–fam75 −5.65 s (1,824 matches, IQR 0.05 s),
    fam26–fam117 +9.25 s, fam58–fam77 −4.25 s (Merlin's numbers reproduced).
  - Epoch test: fam38–fam75 = 5.65 s in 2007–11, 2012–17 and 2018–24 (drift 0.00 s); the 10 biggest pairs drift ≤ 0.05 s.
    The per-family-epoch rule does not fire.
  - Selected rows alone: 26 components, 11 anchorable, 59 families (33 % of rows) on the default offset. Addendum 6a:
    pairs counted on NW rows with P ≥ 0.5 (90,640) and IQR ≤ 0.5 s → 110 families, 15 components (44/19/18/5/4/2…),
    7 anchored with constants −6.6 … −9.9 s, residual rms 0.06 s. Flags: network 89 families (89 % of rows),
    clear_only 7, default 36 (offset −8.47 s).
  - Offsets (network families): range −14.5 … −0.3 s relative to `time + 0.5`, SD 3.5 s (row-weighted 3.2), IQR
    −10.4 / −8.5 / −4.6 s. So the WA catalogue `time` carries a family constant of several seconds (Merlin's point).
  - Anchor check: 24 anchor families; leave-one-out SD (22 anchors in components with ≥ 2 anchors) **3.16 s** — rule
    < 3 s missed, > 4 s stop not reached → proceed, WA headline capped at "present, size provisional". One anchor
    dominates: fam49 (+11.8 s; B quality, 139 rows); without it the LOO SD is ≈ 1.9 s (post hoc, not a rule).
  - Vs: 3.95 km/s straight ray (same as the SVI origins). CVM v1.7 NW column: every offset −1.99 s (lags +1.99 s).
- **Step 2 (`tables/wa3_location_flags.csv`).** 81 of 132 families have ≥ 30 % of their rows co-detected with a family
  > 15 km away; 32 are shallower than 20 km or have no depth (91 suspect in all). Duplicate clusters: 30 clusters hold
  98 families (93 % of rows); 27 clusters span > 15 km (up to 136 km). 66 families (65 % of rows) take another member's
  location for the Green's functions (median move 23 km, 44 moved > 15 km). fam75 and fam121 → fam38's location;
  fam77 → fam58's. Plainly: the 132 WA-NW families are about 64 distinct sources, and many catalogue locations disagree
  by tens of km for the same events.
- **Step 3 (`scripts/s13_units_v3.py` → `data/rows_v3/`, `tables/wa3_rows_verdict.json`).** 37,061 → 20,604 rows after
  the ±2 s dedupe on corrected times: **44 % removed** (rule 25–35 %: miss → headline cap), n_nb2 0, n_nb10 0.144
  (rule < 0.15: pass), kernel-weighted neighbours −0.01 (v2 +0.20). 51 majors built (20 with ≥ 300 LFEs; v2 54).
  WA_BE 220 → 195, WA_BE_broad 934 → 629. No new prep days needed.
- **Independent check, preview (v2 stacks, which never saw the offsets):** the catalogue-only offset + 2 s predicts the
  v2 family grand-stack lag for fam38 (−8.3 vs −8), fam106 (−6.7 vs −7), fam126 (−3.6 vs −5), fam37 (−0.9 vs 0),
  fam52 (−2.6 vs −2); fam75 misses by 11 s (its v2 template used its own mislocated position). Full table comes from
  the summary.
- **Step 5 rows (`scripts/s13_units_v3_svi.py`):** P ≥ 0.9 net_mad terciles, 100 % of pooled rows matched to the
  catalogue net_mad; n per tercile 3,326 / 1,945 / 454 / 196 (2004-07 / 2005-09 / 2012-08 / 2013-09). P still rises
  with net_mad inside P ≥ 0.9 (ρ 0.55–0.81; tercile P medians 0.93 → 0.98), noted as a composition caveat.
- Stacker: `s13_stack2_v2.py` got `--grids` (default unchanged); smoke test on wa3_2017-10 ran (outputs moved to
  `../_quarantine_2026-10-14_envfail/s13_v3_test/`, not used).
- **Queue `scripts/s13_v3_queue.sh` started 8:57 pm PT** (PID `logs/v3_queue.pid`, log `logs/v3_queue.log`): 51
  `wa3_*` + `wa3_between` (kmax 5, grids ±8 / wide / fix2) and 4 `nmt_svi_*` (kmax 1), ≤ 6 parallel. Expected ~45–60 min.
  **Done marker: `V3_QUEUE_DONE`**; the queue then runs `python scripts/s13_v3_summary.py` (→ `logs/v3_summary.out`,
  `tables/v3_verdict.json`, figures `s13_v3_*`).
- **Resume:** if `V3_QUEUE_DONE` exists, check `logs/v3_summary.out` (re-run `python scripts/s13_v3_summary.py` if it
  failed), then write RESULTS.md "Amendment 6 results". If the queue died without the marker: re-run
  `cd svi_science/results/s13_vlfe && nohup bash scripts/s13_v3_queue.sh > logs/v3_queue.out 2>&1 &` (idempotent).

### ~9:20–9:35 pm PT: Amendment 6 stacks done, summary, RESULTS
- `V3_QUEUE_DONE` 9:20 pm PT (56 stacks, none failed); the summary ran from the queue at 9:21 pm PT (`logs/v3_summary.out`,
  `tables/v3_verdict.json`). RESULTS.md "Amendment 6 results" written; pre-edit copy
  `../_quarantine_2026-10-14_envfail/RESULTS_before_amendment6_2026-10-14_2130PT.md`.
- **WA v3:** W1 pass (mode +3 s, median +2, 4 % at ±8). W2–W5 pass (W3 SD 1.01; zones 18/18, 17/20, 19/20;
  odd/even 0.61). 4.1e12 N m per LFE (χ²/dof 0.96) = 0.77 ± 0.04 × SVI. 25/51 detected, 50/51 positive.
  - Fixed +2 s column: 3.9e12, 31/51 detected.
  - Grand stack peaks +3 s pooled (3.9e12) and +2 s zone 1 (5.5e12), up from 2.6 / 4.0e12 in v2.
  - Status "present, size provisional" (caps: anchor leave-one-out SD 3.16 s, dedupe 44 %).
  - Between ETS: WA_BE 3.9e12 (z 3.0), BE_broad 2.4e12 (z 3.4).
- **Family check:** the catalogue-only offsets predict 5/6 strong families' v2 lags within 2 s (median residual
  −0.04 s). The miss is fam75 (v2 used its mislocated position). In v3 they sit at +2 s.
- **SVI net_mad terciles (P ≥ 0.9):** high/low 1.9 / 2.0 / 2.3 / 2.7 → "scales with detection amplitude" (4/4;
  combined 2.9σ). P still correlates with net_mad (caveat). 2004-07 stays low at the same tercile.

## 2026-10-14 (Wed) ~10:20–11:05 pm PT — between-ETS VLFE scan pilot (SVI 2011) started (`PLAN_SCAN.md`, Merlin `MERLIN_VLFE_SCAN.md`)
- Read the design (steps 1–8, F1–F9), RESULTS/PLAN/NOTES, the two Merlin reviews and the scripts. H&G 2018 PDF on disk:
  Table 1 transcribed (13 templates) to `tables/hg2018_table1.csv` with source line. The paper gives no numeric
  between-ETS rate (counts only in Figs 6–7); 2011 SVI templates 2–5 are Mw 3.7/3.3/3.5/3.4 at 54/52/30/55 km.
- Step 1 (`scripts/scan_truth.py`): footprint polygon (111 of 115 SVI families have ≥ 30 P ≥ 0.75 LFEs; hull + 10 km
  = 7,898 km²), `tables/station_epochs.csv`, earthquake events + masks. ComCat pulls (event metadata only, 17 MB in
  `data/eq/`): global M ≥ 4.5 2003–2026 per year, M ≥ 2.5 within 300 km, M ≥ 4 within 1,500 km (needed the sandbox off
  for DNS; nothing else downloaded).
  - Merlin's literal M ≥ 4.5 teleseism mask covers 69 % of 2011 (53–69 % every year). Amendment S1: headline mask keeps
    M ≥ 6 with his windows and M < 6 only if predicted 20-s amplitude ≥ 20 nm, in the surface-wave group window →
    34 % of 2011 (21–34 % per year). Literal masks kept as columns.
- Step 2: `scripts/scan_prep_v3.py` (fixed 600-s taper, remove_response taper=False, gap intervals saved). Acceptance
  PGC 2011-08-20: band RMS v3/v2 = 1.0000 mid-day on all components; edge ratio 2.4–3.0 in the first 300 s. PASS.
  2011 prep: 2,920 station-days in ~12 min (`SCAN_PREP_V3_DONE`; YOUB and LZB only 58 days in 2011).
- Step 3 (`scripts/scan_calib.py`): 290-node grid (289 + S13 centroid), templates in 43 s.
  - First gain run (n-weighted, guard2 drops respected) failed: LRIV Z g 19.6 (LRIV BHZ is broken 2015-12…2018-06,
    σ 89–411 nm, inside the same channel epoch) and SNB Z left with 1 episode (guard2 `<0.05xsta` drops healthy
    verticals). Amendment S2: inverse-variance episode weights, ignore `<0.05xsta`-only drops for calibration.
  - Second run: 12/21 kept (rule ≥ 12: pass). Z gains PGC 0.80, VGZ 1.30, LRIV 1.46 (rule [0.7, 1.4]: LRIV just
    outside → acceptance missed by LRIV; Merlin's failure branch "all ≈ 0.5 or ≈ 2" not met → continue). NLLB Z is
    negative (g −0.62, 4.9 σ) and PGC E is 2.4× the synthetic (30 σ): both excluded by the rule; this is the F5 amplitude
    pattern. SNB Jan–Mar (BH epoch) fully excluded. Kept: PGC N,Z; VGZ E,Z; PFB Z; YOUB Z; LRIV N,E,Z; SNB(HH) N,E,Z.
  - Noise table: SNB Z 130–155 nm Jan–Oct 2011 (excluded, as F7 said), PGC Mar/May/Jun excluded (Z 37–46 nm vs its
    all-days median 11.6 nm), LZB excluded (≈ 1,000 nm).
- Step 4 (`scripts/scan_bp.py`, `scan_detect.py` cupy, `scan_columns.py`). A 4-day code test (Aug 20/21/24/26) found:
  (1) z blow-ups at an edge node when availability changes → Amendment S3 (z from num/√den); (2) instrument glitches
  (PGC E ~10 mm at every UTC midnight; YOUB Z 2e7 nm bursts) → Amendment S4 glitch mask (> 200 σ_day, ±600 s; ≤ 2.5 %
  of each component's time); (3) every event has an opposite-sign side lobe 13–17 s away (antisymmetric band-passed
  step), and polarity-scrambled templates keep part of real signals → Amendment S5 (joint-sign clusters; |ρ| ≤ 0.3
  draws; filter validation on injections with the 85 % rule; max_share fallback = excess over expected share).
  - Code-test only (not a result): all four H&G 2011 times had a positive peak within ±4 s (z 18.0, 4.4, 18.8, 6.6;
    Mw 3.42, 3.18, 3.44, 3.21 vs H&G 3.7, 3.3, 3.5, 3.4). At 04:09 the 2–8 Hz envelope jumps 15–33× the day median
    9–24 s after the VLFE origin at all five stations (no catalogued quake) — worth a look in step 8 / test 1.
- Year queue `scripts/scan_queue.sh` started 10:57 pm PT (PID `logs/scan_queue.pid`, log `logs/scan_queue.log`): real
  + 50 polarity draws, 5 station-day-shift controls, injection (320 events), then columns. Marker `SCAN_QUEUE_DONE`.
  Resume: `cd svi_science/results/s13_vlfe && nohup bash scripts/scan_queue.sh > logs/scan_queue.out 2>&1 &`
  (idempotent), then `python scripts/scan_loc6e.py` and `python scripts/scan_validate.py`.

## 2026-10-14 (Wed) ~11:35–11:50 pm PT — scan pilot year results, verdict NO-GO
- `SCAN_QUEUE_DONE` 11:35 pm PT. First `scan_validate.py` pass: 0/140 "20-s" and 1/40 perturbed injections recovered →
  bug: τ = D in the erf moment function (spectrum exp(−(πfτ)²) = 0.03 at 0.03 Hz for τ = 20). Amendment S6
  (τ = D/2.33); injections re-run (`scripts/scan_inj2.sh`, marker `SCAN_INJ2_DONE`, ~6 min); old run kept (`runs/inj`,
  `inject_2011_tauD.parquet`, `tables/*_pass1.*`).
- Verdict (`tables/scan_verdict.json`): NO-GO. Controls never reach (a) ≤ 1/yr, (b) ≤ 2/yr: at z ≥ 7 (a) 714/yr, (b)
  462/yr between ETS vs real 515/yr. 6a 1.6× (rule 10×); ETS sign 33 vs 21; 6b 3.84e12 (0.73×; 5.18e12 with g = 1;
  within noise, not a bug); 6c 2/4 (H&G 2 and 4 at z 18, Mw 3.42/3.44; 3 at 44 km; 5 at z 6.6); 6d Mw50/80 summer
  3.42/3.59, winter 3.45/3.73, 82 % at 3.6 summer, 24 % at 3.3; location median 5.9 km at 3.6; perturbed GF 5 %;
  6e pass (17/13 km).
- Filter validation (S5.3): max_share keeps 22 % of injections → share_excess fallback (99 %); all others ≥ 93 %.
- Post-hoc diagnosis (`scripts/scan_tail.py`, `scan_loso.py`, `scan_posthoc_combo.py`): false alarms at high z are
  single-station (max share 0.99) at nodes under LRIV/PFB/VGZ; leave-dominant-station-out kills control (a) (714 → 4/yr
  at 7) but real = (b) remains (69 vs 66/yr); 60 % of those are inside the literal M ≥ 4.5 mask, the rest mostly at grid
  edges; with literal mask + footprint nodes, real 7.6/yr vs (b) 12.7/yr — no excess.
- Step 8 not run (NO-GO). Figures redrawn by `scripts/scan_figs.py`. RESULTS.md "Scan pilot (2011)" written
  (pre-edit copy `../_quarantine_2026-10-14_envfail/RESULTS_before_scan_pilot_2026-10-15_0010PT.md`).
- Next decision (coordinator/Merlin): re-design per NO-GO branch (ii) (leave-one-station-out requirement, literal
  mask or data-driven surface-wave detection, footprint-limited grid, 25-s corner, more stations) or the
  stack-at-candidate-times fallback.

## 2026-10-15 (Thu) ~0:10–1:30 am PT — scan pilot 2 (SVI 2011): Merlin NO-GO review steps 2–7 (`MERLIN_SCAN_NOGO_REVIEW.md`)
- Read the review, PLAN_SCAN S1–S6, RESULTS "Scan pilot (2011)", the scan scripts and run tables. Wrote PLAN_SCAN
  **Amendment S7** (0:15–0:50 am PT) with every rule before any second-pilot number.
- **References (S7.1).** Prep v3 (`scan_prep_ref.py`, wrapper of `scan_prep_v3.py`; responses from data/resp, else
  `lfe_array_discovery/results/stationxml`): GNW + DOSE all 2011, ETS-days set DOY 211–239; 991 station-days ok
  (OBSR 8 and D03D 6 days without data). Band-pass + S4 glitch mask (`scan_bp_ref.py`, = scan_bp.py).
  - Reference QC (`scan_refqc.py`): teleseism (M ≥ 6.8, 2011) 20–50 s Z peak ratio vs PGC. GNW 1.00 (27 events,
    0.71–1.25), DOSE 1.07. ELW 0.46 and WISH 0.48 → excluded by the [0.5, 2] rule; LEBA 0.54, RATT 0.53, SP2 0.55,
    NLWA 0.54 pass narrowly (a factor ~2 low — probably response metadata; noted, not acted on). No month excluded by the
    noise rule; GNW monthly Z noise 8–19 nm (2011 median 11.5), 96 % finite, glitch mask 1–3 % per month.
  - Code note: the QC amplitudes are computed on band-passed data WITHOUT the S4 glitch mask (`bp_window`), because the
    200-σ mask removes the surface waves of every M ≥ 6.8 teleseism (first try: 0 usable events).
- **Masks (S7.6, `scan_masks2.py`).** 2011 fractions: head 33.4 %, head2 (Ms = mb − 0.3) 46.8 %, bw5 8.6 %, reg3 4.6 %
  (≤ 15 % → in v2), reg 0.7 %; **mask v2 = 52.4 %** of 2011 (v1 33.9 %). ComCat M ≥ 3 within 2,500 km 2011 fetched
  (1,216 events, metadata only, `data/eq/comcat_reg2500_M30_2011.csv`).
- **Known body-wave list (S7.4, `scan_bwlist.py`).** The rule gives **8** of the 47 pilot-1 survivors (Merlin by eye: 11):
  all teleseisms M 5.1–5.5 at 8,800–15,500 km, 220–1,210 s after P (most 100–250 s after S). `tables/scan2_bw_known_list.csv`,
  all 47 with matches in `tables/scan2_survivors47.csv`. 80 % rule → ≥ 7 of 8 must be flagged.
- **Decay (S7.2, `scan_decay.py`; stacks `scan2dec_svi_2011-08` (v3) and `scan2dec_svi_2012-08` (v2), pooled unit).** GNW
  is not detected in either episode; upper bound R_up = 0.35 (2011-08) and 0.17 (2012-08) ≤ 0.5 → PASS; predicted (whole
  space) 0.19 / 0.13. Caveat: PGC Z itself is only ~1.7 σ in a single-station stack (A 1.26 nm vs σ 0.75 nm in 2011-08),
  so the ratios are noise-limited upper bounds. Figure `figures/scan2_decay_stacks.png`.
- **Step cleaner (S7.5, `scan_detect2.py`).** Code fix during the code test: ±600-s dilation by cumulative sums (scipy
  binary_dilation with a 1,201-sample structure took ~20 s per station-year). Cleaned runs real_c (+50 draws), ctrlA_s*_c,
  inj2_c, inj3_c done (~1 s/day + 100 s cleaner pass). It removes 35 % of LRIV's year, 13 % VGZ, 12 % PFB, 10 % PGC.
  **Acceptance FAILED:** Mw ≥ 3.6 injection recovery 0.23 vs 0.82 in pilot 1 (keep 29 %, rule ≥ 85 %) and control (a) at
  z ≥ 7 = 1,061/yr (rule ≤ 20; cluster definition now positive-only, so not strictly comparable with pilot 1's 714). The
  rule flags the station that carries a real nearby event (others < 2.5 at that station's best node). → cleaner NOT used
  (S7.5), fallback to the uncleaned pilot-1 runs + `inj3_u` (S7.12, written ~1:05 am PT before any fallback number).
  S7.12 also states that fp_interior falls under the 85 % rule (cleaned runs: 0.66).
- First `scan_validate2.py` pass crashed (node dtype) right after the cleaner and filter tables; no veto number existed.
- **Far references (S7.13).** User approval relayed ~1:15 am PT. EarthScope LHZ 2011 (`scan_fetch_ref2011.py`,
  `logs/ref_download.log`): IU.COR, US.NEW, US.HAWA + substitutes UW.OMAK (for CN.PNT) and TA.I02D (for CN.BBB) — CN.PNT
  and CN.BBB have only 40-Hz BH at NRCan/EarthScope (~2 GB/station-year of BHZ); no LH station exists north of the
  network. 285 MB downloaded (more than the ~150 MB estimate; LHZ only), 0 missing months; day files on disk 954 MB
  (obspy rewrite). Prep (LH path, Z only), band-pass, QC: all five pass (median ratio vs PGC 0.77–0.93, 29 events).
  Veto columns recomputed with the far set; `ref_excess_nf`/`ref_ratio_nf` keep the without-far version.
- Running: `scripts/scan2_fallback_queue.sh` (columns on real, ctrlA_s*, inj2, inj3_u, then validate2), marker
  `SCAN2_FALLBACK_DONE`, log `logs/scan2_fallback.log`, validation output `logs/scan2_validate.out`.

### ~1:20–1:55 am PT: fallback runs, validation, verdict NO-GO
- `scan2_fallback_queue.sh` (columns with 17 references on real, ctrlA_s*, inj2, inj3_u; validate2) done 1:38 am PT.
  Bug found and fixed in validate2: the body-wave rows were merged with all 51 variants of the pilot-1 loso table
  (18 rows instead of 8); re-run 1:45 am PT.
- Filter re-validation on the fallback injections (Mw ≥ 3.4, n = 220): all criteria ≥ 0.91 except max_share 0.26 →
  share_excess (0.98). fp_interior keeps 0.92 of family injections Mw ≥ 3.6 → used.
- Veto validation: V1/V2/V3 flag 3 of 8 body-wave clusters with the far set (0 of 8 without); H&G 2/4 never flagged;
  injections flagged 0–2 %. → no variant validated. Post-hoc: ref_ratio separates (body waves 0.67–1.11, injections mostly
  < 0.5, H&G 2/4 0.34/0.46); the excess is limited because M 5 teleseismic S is only 1–7 σ at each far station.
- Threshold grid: no pair meets (a) ≤ 0.3 and (b) ≤ 0.7. Lowest (a)+(b) with ≥ 50 % Mw-3.8 recovery: z ≥ 20, z_loso ≥ 10 →
  (a) 1.4, (b) 3.9/yr, recovery 58 % (n 48), location 6.4 km; real 0 detections (ETS 0, between 0); Poisson p = 1.
  At z ≥ 10, z_loso ≥ 10: real 5 clusters (H&G 2, H&G 4, 2011-08-28 01:08:29 UTC ETS Mw 3.6 HF-loud ref_ratio 0.36,
  09-04 15:15, 10-22 16:11 between ETS with ref_excess 4.4 / ratio 1.0), (b) 7.8/yr.
- Classification: 08-18 21:19 UTC "unresolved" by rule (ref_excess 3.9, no validated veto); far east stations at 7–12 σ,
  plane wave 0.29 s/km from ~283° → most likely a wave from offshore WNW, not a northern VLFE. 08-22 triple "local
  candidate" by rule but ref_ratio 1.2–1.5, far stations 3–6.5 σ, edge nodes outside the footprint.
- RESULTS.md "Scan pilot 2 (2011)" written (pre-edit copy `../_quarantine_2026-10-14_envfail/RESULTS_before_scan_pilot2_2026-10-15_0150PT.md`).

## 2026-10-15 (Thu) ~1:50–2:40 am PT — two parallel tracks after the pilot-2 NO-GO (coordinator instruction)
- Read RESULTS "Scan pilot (2011)" / "Scan pilot 2 (2011)", PLAN_SCAN S1–S7.13, MERLIN_SCAN_NOGO_REVIEW (fallback
  paragraph), MERLIN_VLFE_SCAN, the scan and S13 scripts. H&G 2018 PDF: "Supporting Information S1, Table S1" is only
  referenced (not on disk); the paper's 2011-07 → 2012-12 catalogue appears only as counts in Figs 6–7, so no 2012 H&G event
  exists to test against.
- **Track A = PLAN_SCAN Amendment S8** (written ~1:55–2:15 am PT, before anything ran on 2012). Ratio veto R1 fixed from the
  2011 tables: flag if n_ref ≥ 2 and ref_ratio > 0.551 (geometric midpoint of H&G 4 0.455 and the smallest 2011 body-wave
  ratio 0.668). 2011 in-sample: 8/8 body waves, 0/2 H&G, 3.8 % of family injections Mw ≥ 3.6.
  - Code: `scan_common.py` reads `SCAN_YEAR` (default 2011 = unchanged; 2012 → `data/scan_2012/`, NDAY 366, ETS
    2012-08-31 → 09-28, tables `scan3_*` via `S.tab()`); same switch in `scan_masks2.py`, `scan_refqc.py`, `scan_refs.py`
    (bp_exists dir, 2012 StationXML), `scan_inject*.py`, `scan_detect.py` (noise table), `scan_prep_ref.py` (raw_ref2012).
    New: `scan_fetch_ref2012.py`, `scan3_noise.py`, `scan3_queue.sh`, `scan_validate3.py`.
  - Downloads: EarthScope LHZ 2012 for IU.COR, US.NEW, US.HAWA, UW.OMAK, TA.I02D = 286 MB, 0 missing months
    (`logs/ref_download2012.log`, marker `REF2012_DOWNLOAD_DONE`); ComCat M ≥ 3 within 2,500 km, 2012 (metadata, 1,353 events).
  - Queue `scripts/scan3_queue.sh` started 1:54 am PT (≤ 5 workers; log `logs/scan3_queue.log`; per-step markers
    `data/scan_2012/STEP_*_DONE`; final marker `SCAN3_QUEUE_DONE`). Then `SCAN_YEAR=2012 python scripts/scan_validate3.py`.
- **Track B = PLAN.md Amendment 7** (written ~2:15–2:35 am PT, before any row/prep/stack). Day sample (catalogue only):
  7,527 between-ETS days, 1,118 eligible 7-day blocks, 224 drawn (20 %) → 1,488 days. Candidates P ≥ 0.1 between ETS
  110,007 → 106,857 after the ±2 s cross-family dedupe; drawn per class: C075 859, C050 1,706, C030 3,273, C010 11,046
  (C075 is clustered on few days, so the 20 % day sample holds only 8 % of it). PNSN epochs in the footprint between ETS:
  12,656 (10,926 with no P ≥ 0.75 row within ±5 min, 6,523 with no P ≥ 0.1 row) on 403 days. Prep 2,523 days × 5 stations
  (8,774 station-days missing from v2).
  - Code test of the envelope origins on 2012-02-04 (5 stations, 102 epochs): sensible peaks (2–3× the day median; one
    epoch at 21× — probably a local earthquake; window QC handles 20–50 s outliers).
  - Queue `scripts/s13_be3_queue.sh` started 1:57 am PT (≤ 3 workers; log `logs/be3_queue.log`; markers `PREP_be3_DONE`,
    `STACK_be3_{hi,c030,c010,pnsn}_DONE`, `BE3_QUEUE_DONE`). Then `python scripts/s13_be3_summary.py`.
- ~2:06 am PT: track B rows built (`logs/be3_rows.out`): C075 859, C050 1,706, C030 3,273, C010 11,046; PNSN_noLFE 7,333 and
  PNSN_noCand 4,758 windows after the ±60 s dedupe (10,925 of 10,926 epochs had ≥ 2 envelopes); BE_broad 688. Envelope
  origin minus PNSN epoch time: median +102 s (5–95 %: −127 … +283 s), so the peak search window matters (stated, not changed).
  Prep of 8,774 station-days running at 3 workers (~3 h). Track A 2012 test of `scan_masks2.py` (code check, rerun
  by the queue identically): mask v2 = 46 % of 2012.
- Stopped here (stages > 20 min). Resume: check markers below; queues are idempotent and can be restarted with the same
  nohup commands. Once `SCAN3_QUEUE_DONE`: `SCAN_YEAR=2012 python scripts/scan_validate3.py > logs/scan3_validate.out 2>&1`.
  Once `BE3_QUEUE_DONE`: `python scripts/s13_be3_summary.py > logs/be3_summary.out 2>&1`. When track A ends, track B can
  be relaunched with more workers (`BE3_WORKERS=6`, kill the running queue first; prep skips existing files).

## 2026-10-15 (Thu) ~3:35–4:05 am PT — scan pilot 3 (2012) results; track B restarted with 8 workers
- `SCAN3_QUEUE_DONE` 3:35 am PT (prep 1:55–2:30, band-pass/refs/QC/injections to 2:41, detector runs ~1.3 min each,
  columns 3:06–3:35). Reference QC 2012: all 17 pass (WISH 0.87 now). Noise rule: PGC Oct–Dec excluded. Mask v2 46 %.
- `SCAN_YEAR=2012 python scripts/scan_validate3.py` (58 s; log `logs/scan3_validate.out`). One fix after the first run:
  the separation figure had an empty third panel (z_loso set dedup) → panels z_loso 5/7/10; numbers unchanged.
- Verdict NO-GO (`tables/scan3_verdict.json`):
  - R1 flags 7/7 known body-wave clusters but 29.5 % of family Mw ≥ 3.6 injections → unvalidated. The cause is the
    reference mix: five 100–150 km references are live all year in 2012 (6 refs at injection times in 2011, ~11 in 2012),
    which raises the injection ratio median from 0.17 to 0.36.
  - Column R1_fixed6 (GNW + far set only) passes in 2012 (7/7, 94.7 % kept); a post-hoc lesson, not a test.
  - Best pair with R1: z ≥ 5, z_loso ≥ 10: (a) 0.76, (b) 0.0/yr between ETS (the pair rule fails on (a) > 0.3);
    Mw-3.8 recovery 51 % (z_loso ≥ 10 alone keeps 80 %), location 6.4 km.
  - 1 ETS detection (2012-09-19 19:18:46 PDT, Mw 3.45, ratio 0.28, HF-loud, 16 LFEs within ±10 min vs 1.7
    expected), 0 between ETS; Poisson p 0.048.
  - Without a veto: real between 7.3/yr = (b) 7.1/yr at z ≥ 20, z_loso ≥ 10 (as in 2011). The 3 between-ETS clusters
    that pass without a veto (06-26, 12-05, 12-31) have reference ratios 0.70–1.05 and no LFEs.
- RESULTS "Scan pilot 3 (2012, ratio veto)" written (pre-edit copy
  `../_quarantine_2026-10-14_envfail/RESULTS_before_scan_pilot3_2026-10-15_0400PT.md`).
- Track B: prep was 51 % done at 3 workers (6,379 of 12,615 station-day jobs). Killed 3:40 am PT and relaunched with
  `BE3_WORKERS=8` (track A finished; prep skips done files, atomic writes, no partial files found). First-run log kept as
  `logs/prep_be3_run1_3workers.out`. Markers unchanged (`PREP_be3_DONE`, `STACK_be3_*_DONE`, `BE3_QUEUE_DONE`), then
  `python scripts/s13_be3_summary.py`.

## 2026-10-15 (Thu) ~4:40–5:00 am PT — track B results (between-ETS stacks at candidate times)
- Queue after the 8-worker restart: prep end 4:11 am PT (3,628 new ok, 7,405 existing), stacks be3_hi 4:22, c030 4:25,
  pnsn 4:30, c010 4:38 (`BE3_QUEUE_DONE`). `python scripts/s13_be3_summary.py` → `tables/be3_summary.csv`,
  `be3_verdict.json`, `figures/s13_be3_classes.png`, `figures/s13_be3_stacks_*.png`.
- Numbers (M0 per candidate, z; ρ = /5.39e12):
  - C075 4.0e12, z 1.8, ρ 0.74 ± 0.41; C050 3.1e12, z 1.7, ρ 0.57 ± 0.33;
  - C030 −0.7e12, ρ −0.14 ± 0.18; C010 −0.4e12, ρ −0.07 ± 0.09;
  - PNSN_noLFE 1.0e12, z 1.8, ρ 0.19 ± 0.10; PNSN_noCand 0.9e12, z 1.4;
  - BE_broad (5 stations) 3.8e12, z 2.6, ρ 0.71 ± 0.28;
  - no unit reaches z ≥ 3; ±1-day |z| ≤ 2.5 everywhere.
- Mixture reading: C030 and C010 sit 2.7–2.9σ below their mean P → low-P candidates carry less than P predicts. Tremor-day
  candidates carry more than quiet-day ones in every class (column).
- Control SD per candidate is ~2× the ETS-stack value (5 stations, all seasons).
- RESULTS "Between-ETS stacks at candidate times" written (pre-edit copy
  `../_quarantine_2026-10-14_envfail/RESULTS_before_be3_2026-10-15_0450PT.md`), with the hg_repro relation. WA not run.
- One correction while writing: the scaling remark for the full C075 set is z ≈ 6 (√12 × 1.8), not 3.5.
