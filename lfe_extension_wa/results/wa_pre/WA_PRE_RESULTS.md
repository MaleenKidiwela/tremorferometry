# Washington pre-production steps (while the downloads run)

Written 2026-10-02 ~22:50Z. Everything is under `lfe_extension_wa/results/wa_pre/`. Scripts: `scripts/pre_*.py`, queue
`scripts/wa_pre_queue.sh` (frozen copy `wa_pre_queue.launched.sh`, script hashes `scripts_at_launch.sha256`), log
`logs/wa_pre_queue.log`, markers `WA_PRE_{PID,HEARTBEAT,DONE,FAILED}` and `stage_<name>.done`. No downloads, no deletions,
no git. Nothing that existed was overwritten. The download roots `traces_dl/` and `traces_dl2/` were never read; neither
was `lfe_array_discovery/traces_sealed`. The sealed months 2007-01 and 2009-05 and the excluded month 2011-04 were not
used in any run. "Scrambled control" = the same detector with each station site's data shifted to a different time of
the same day, so nothing lines up between sites except by chance.

## In short

1. **Freeze, completeness and day roles: done.** 139 input files hashed (1.39 GB). Of the 147,445 expected core
   station-days (2007-01-01 to 2026-09-15), raw files exist for 131,597. None fails the first-record check, and none of
   the 11,769 WA cache files fails its header check. 10.75 % are missing: 30.6 % in 2007–2011, 14.5 % in 2012–2013 and
   3.1 % in 2014–2026. 806 of 4,859 station-months miss more than 3 %, every one for lack of a file on disk (none
   because a file is unreadable). Day-role files are written.
2. **Station health (W2): running** as a detached queue. W2 raw is the slow part (~3 h). Its section below is filled in
   automatically when the queue ends.
3. **Olympic Arm A dominance cap: not adopted. Production proceeds uncapped, and Olympic A gets stacks / group rates
   only.** Olympic A data/control at the ≤ 3-km group level was 1.00 with the 0.3 cap (needed ≥ 3; 1.00 uncapped) and
   1.42 when the dominant site is left out. NW precision did not stay "unchanged within 0.02": it **rose**, from 0.318 to
   0.365 with the cap and to 0.676 with the dominant site left out (but recall then falls from 0.124 to 0.018). By the
   fixed rule this is a fail. The NW gain is an open question for the main session or Merlin. I have not acted on it.
4. **Control clean-up: fails its first test. The two controls do not agree.** Contaminated share on PNSN-quiet days:
   NW 0.90, Olympic A 0.98, all regions 0.97 (needed ≈ 0; the chance level is 0.15). The mark is real, not a sign or
   timing mistake: with the shift sign flipped, or the times moved by ±97 s, the share drops to the chance level 0.15.
   But it marks the dominant site's own single-site burst, which also makes data candidates at the true time. It is not
   marking real LFEs. At ≥ 8.75 × MAD on the 21 core stations there is **no control excess on tremor days** for NW
   (7.5 per family-day on tremor days vs 8.4 on quiet days) or Olympic A (228 vs 282). So "share of the excess
   explained" cannot be computed there. Arm B (0.93) and SW (0.77) have an excess and pass the 70 % bar, but this means
   little while the first test fails. On ETS days the clean control and the time-reversed-template control do **not**
   agree within 1.5×: reversed / clean is NW 1.85, Olympic A 27.8 and all regions 9.2.

## Step 1 — W0 freeze, W1 completeness, day roles

- **W0** `w0/frozen_inputs.csv`: sha256, size and mtime of 139 files (1.39 GB). It covers the WA scripts, all
  templates, the four setups, the confidence and PNET models, the HF cuts, the node, tremor-day and episode tables, the
  PNSN catalogue, families v1.1/v1.2, truth_array_frac, the inventory, the plan, and the SVI/South scripts that
  production copies. The queue re-hashes them at its end (`w0/verify_*.csv`).
- **W1** (`scripts/pre_w01.py complete`; `w1/w1_station_days.parquet`, `w1_station_month.csv`, `w1_failing_months.csv`,
  `w1_station_summary.csv`, `w1_era_summary.csv`, `w1_raw_outside_epochs.csv`):
  - **Expected** = inside a metadata epoch of the station (the node tables, which follow the FDSN epochs).
  - **Raw**: a file in any raw root on disk, excluding the download and sealed roots.
  - **Readable**: the first miniSEED record parses and is an HH, BH or EH channel; the record length is also recorded.
  - **Cache**: the WA-chain 20 Hz file is present, its zip directory and `m`/`live` members are read, and the `x` header
    has shape (3, 1,728,000).
  - A full headonly read of every file ran at ~140 MB/s (~6 h for 2.9 TB), so I stopped it after 1 h and replaced it
    with this quick check.
  - The **full** read of every raw day file is done by W2 raw for the station-days that SVI R2 has not already read in
    full. The final readability that combines both reads is written to `w1/w1_readability_final.csv` when the queue
    ends; its line is in the W2 section.

| era | expected | readable on disk | missing | WA cache files (errors) |
|---|---|---|---|---|
| 2007–2011 | 34,633 | 24,024 | 30.6 % | 1,756 (0) |
| 2012–2013 | 15,351 | 13,128 | 14.5 % | 6,327 (0) |
| 2014–2026-09-15 | 97,461 | 94,445 | 3.1 % | 3,686 (0) |

  - Station-months missing more than 3 %: 474 in 2007–2011, 96 in 2012–2013, 236 in 2014–2026. All are "no raw file on
    disk".
  - The 2007–2011 gaps are mostly PGC, VGZ, B009–B011 and B926 (whole-station missing fractions 0.24 / 0.18 /
    0.32 / 0.28 / 0.32 / 0.29).
  - In 2014–2026 they are mostly B014 (33 months, incl. 2023–24), B004 (2019–20), LRIV (2015, 2018–19), VGZ (2022) and
    the 2026-06 → 09 tail of the 8 Olympic-set stations. Tier A of the running download is meant to fill these.
  - 270 readable station-days sit only in roots the WA cache builder does not read (`lfe_picker/traces`,
    `lfe_south/traces`, `data/waveforms`). Production should add those roots or link the files.
  - Stage guard (unreadable among present files > 3 %): 0 %, passes.
- **Day roles** (`w1/family_month_roles.csv`, `w1/family_month_roles_olympic.csv`, `w1/day_roles.csv`):
  - truth_array_frac per family-month is reused from `wa_gonogo/truth_array_frac.csv` (1,583 NW/SW family-months).
  - Roles:
    - **sealed one-look headline** = 2007-01 (211 families, 34,463 truth detections, median array fraction 0.71);
    - **sealed one-look, semi-independent** = 2009-05 (132 families, 0.56);
    - **excluded** = 2011-04 (0 array stations);
    - **fit** = 2007-10, 2008-05, 2009-02/04, 2010-03, 2010-08, 2011-08, 2011-12. The 2008-05 and 2010-08 build
      halves are in-sample; 2011-08 is in-sample for the stations whose templates come from tpl_step4b_a (D1).
  - Olympic A/B: 204 family-months in 2010-08/09, "fit out-of-fold only" (no sealed check exists for them).
  - Days 2007-01-01 → 2026-09-15: fit 305, sealed headline 31, semi-independent 31, excluded 30, other 6,801. Each truth
    day carries the truth runs that scanned it.

## Step 3 — Olympic Arm A dominance cap (V2010, 21 core stations, 20 usable)

- Method:
  - `scripts/pre_mf_cap.py` is a copy of the frozen `wa_mf.py`, changed only in how the site weights are set.
  - In frozen mode it reproduces the go/no-go run exactly (2 days, every candidate identical, max |Δnet_mad| = 0).
  - Variants: **cap30**, where no site's share of a family's site weight may exceed 0.3, enforced by water-filling (no
    family-day was infeasible); **lodo**, where the site with the largest weight share is left out.
  - Both variants use the same two scrambled controls (same seeds) and the same HF stage and scorer
    (`gonogo_score.py --wins V2010`).
- Definitions were written before the runs in `cap/PREREGISTERED.md`. Because no P model exists for these runs,
  "Olympic A nbr3 data/control" = the count of ≥ 8.75 × MAD candidates, data vs the mean of the 2 controls, merged at
  ±2 s within ≤ 3-km family groups, on the test half of the observable episode days.
- Results (`cap/cap_comparison.csv`, `cap/cap_verdict.json`):

| run | median Olympic A dominant-site share | Olympic A ≤3-km group data / control | Olympic A own / ≤3-km precision (control) | NW own precision (control) | NW own recall | NW family-windows ≥ 2× both controls (control vs control) |
|---|---|---|---|---|---|---|
| frozen | 0.54 | 63,363 / 63,469 = **1.00** | 0.048 / 0.080 (0.001 / 0.017) | **0.318** (0.009) | 0.124 | 81 (4) |
| cap30 | 0.30 | 38,673 / 38,733 = **1.00** | 0.082 / 0.141 (0.001 / 0.018) | **0.365** (0.007) | 0.107 | 111 (4) |
| lodo | 0.16 | 32 / 22.5 = **1.42** | 0.321 / 0.321 (0 / 0) | **0.676** (0.007) | 0.018 | 128 (19) |

- **Verdict: not adopted** (needed Olympic A ratio ≥ 3 AND NW precision within 0.02).
  - Cap30 fails (a): ratio 1.00.
  - It also fails (b) literally: NW precision is 0.048 higher, not "unchanged".
  - Production proceeds uncapped. Olympic A is followed by stacks and group rates only.
- Leaving B013+B943 out removes almost all of Olympic A: 95,675 → 32 candidates. Its single repeats therefore come
  almost entirely from those two boreholes.
- **Open, not decided here:** the cap improves NW: precision 0.318 → 0.365, family-windows 81 → 111, at the cost of
  recall 0.124 → 0.107, with control-vs-control unchanged at 4. Lodo makes NW much purer (0.676), but recall drops to
  0.018. The rule as written fails on improvement as well as on loss. Whether to cap NW only is a decision for the main
  session or Merlin, with a test on another episode before any use. This is a single 39-day episode.

## Step 4 — scrambled-control clean-up and the time-reversed-template control

- **Method** (`scripts/pre_ctlclean.py`):
  - For every control candidate, find the dominant site: the per-site z at its rolled position, against 1,000 random
    origins of that site-day.
  - Un-roll it: t_real = t − shift_site / 20 Hz.
    - The code builds the control as `roll(x, +shift)`, i.e. control[i] = data[i − shift]. So the real time is t − shift,
      not "t + shift" as the task text says.
    - The opposite sign was run as a check.
  - Mark the candidate contaminated if a data candidate of any family (net_mad ≥ 8) arrives within ±2 s **at that
    site**, using each family's template offset at that site.
  - Variant columns: origin-time match; coherent data candidates only (≥ 2 sites at z ≥ 3); other family only; chance
    (±97 s); flipped sign.
- **Runs:**
  - The go/no-go 21-core run (V2010 + A2011).
  - A new 21-core run on 60 trial-year days, 2010–2014 (`ctl/days60.csv`). Each day has ≥ 17 core stations cached.
    There are 20 PNSN-quiet days (0 windows in every WA table that day and the days either side), 20 moderate days
    (1–80 windows, NW < 100) and 20 ETS days (NW ≥ 100), picked with a fixed seed.
  - The 42-station T1 run with 4 controls.
  - The same 60 days with **time-reversed templates**: every channel reversed, same offsets and weights
    (`ctl/res/mf_rev60`).
- Bars:
  - T1: contaminated share on quiet days ≤ 0.05 (my reading of "≈ 0").
  - T2: explained = 1 − (clean_trem − clean_quiet) / (ctl_trem − ctl_quiet) ≥ 0.70.
  - T3: time-reversed / clean within 1.5× on ETS days.

Tables: `ctl/ctl_tests.csv`, `ctl/ctl_summary.csv`, `ctl/rev_compare.csv`.

| run (21 core unless noted) | region | control / fam-day: quiet → tremor | contaminated share: quiet / tremor (chance) | coherent-only share: quiet (chance) | T2 explained |
|---|---|---|---|---|---|
| go/no-go V2010+A2011 | NW | 8.40 → 7.46 | **0.90** / 0.93 (0.16) | 0.26 (0.035) | no excess |
| go/no-go | Olympic A | 281.6 → 227.6 | **0.98** / 0.99 (0.16) | 0.41 | no excess |
| go/no-go | Arm B | 2.55 → 5.87 | 0.82 / 0.88 | 0.39 | 0.93 |
| go/no-go | SW | 0.012 → 0.051 | 0.11 / 0.61 | 0.00 | 0.77 |
| 60 days | NW | 7.33 → 5.75 | 0.92 / 0.93 (0.16) | 0.28 | no excess |
| 42-station T1 | NW | 4.86 → 5.90 | 0.84 / 0.91 (0.13) | 0.36 | 1.22 (over-cleans) |

- **T1 fails everywhere except SW** (0.11, still above 0.05).
  - The marking is not chance: the shift-sign flip and the ±97 s shifts give ≈ 0.15.
  - But almost every scrambled-control candidate is carried by one site, B013+B943 in 99 % of them (smoke day
    2010-08-20). That site's own burst also makes data candidates at the true time, mostly in other families: the
    "other family" share is 0.78–0.97 on quiet days. So the method marks single-site bursts, not real LFE activity.
  - The coherent-only variant, which asks for a data candidate seen at ≥ 2 sites, still marks 0.26–0.43 on quiet days
    (chance 0.035). It explains only 0.17 of the excess where an excess exists.
- **T2** cannot be evaluated for NW or Olympic A at 8.75 × MAD on the core: the control is not higher on tremor days.
  Data and control are nearly equal on quiet days (NW 8.45 vs 8.40 per family-day, go/no-go run). This repeats the
  go/no-go finding that single repeats at this threshold are mostly single-site events. It holds for this network and
  threshold, not for the method in general (the P model separates them on 42 stations).
- **T3, 60 days** (per family-day, NW):

| day class | data | scrambled control | clean (spec) | time-reversed |
|---|---|---|---|---|
| quiet | 7.37 | 7.33 | 0.56 | 1.51 |
| moderate | 4.92 | 4.84 | 0.36 | 0.94 |
| ETS | 8.70 | 6.59 | 0.44 | 0.82 |

  - Reversed / clean on ETS days: NW 1.85, SW 0.82 (agrees), Olympic A 27.8, Arm B 0.22, all regions 9.2. So **no**
    agreement within 1.5× (SW only). The other-family variant happens to agree for NW (0.80), but that variant was not
    the specified method.
  - The time-reversed control fires ~5–8× less than the scrambled control (NW ETS 0.82 vs 6.59). It is lower on ETS
    days than on quiet days. So, unlike the scrambled control, it does not seem to pick up real activity. One possible
    reading, not tested here: the B013+B943 bursts that survive scrambling look like the templates, and reversing the
    templates removes that match.
- **Verdict:** the stored-shift clean-up does not meet its bars, so I would not use it to clean the control before the
  cut-off table or the rates.
  - Option: the time-reversed control as a second false-alarm reference.
  - Option: a control-only coherence rule (≥ 2 sites at z ≥ 3) that is not checked against data candidates.
  - Both options are for the main session or Merlin to decide. Neither was tuned here.

<!-- W2-AUTO-START -->
## Step 2 (W2) — station health of the 21 core stations, 2007 → 2026-09-15 (auto-written by scripts/pre_report.py)

- Raw metrics: 87112 station-days read by this run (72613 ok, 5 errors, 14494 expected but no file); the 13 SVI-shared stations 2014-01-01 → 2026-09-15 come from SVI R2's daily table (same files, same code). Files: `w2/r2_health_daily.parquet`, `w2/r2_health_stationday.parquet`.
- Stack test: 261 strongest-tremor days (2007–2026; sealed 2007-01 / 2009-05 and 2011-04 never used), cache {'ok': 3513, 'wa_cache': 1318, 'no_raw': 409, 'cached': 230, 'raw_not_prod_visible': 10, 'empty': 1}; 1111 station-months tested (≥ 3 families with ≥ 5 confident repeats, ≥ 30 repeats). Confident repeats for a station = detections whose network WITHOUT that station's site still reaches 8.75 × MAD with ≥ 2 other sites at site z ≥ 3. Scrambled control = the same station read on other selected days at the same time of day.
- Station-month status (all 4859 metadata-alive station-months): {'inherit_ok': 2508, 'inherit_control_level': 808, 'ok': 684, 'control_level': 426, 'nodata': 418, 'dead': 13, 'untested': 1, 'inherit_dead': 1}.
- Tested months at control level: 426 of 1111. Excluded station-months (health reason, not plain no-data): 1248 (`w2/r2_exclusions.csv`).

| station | tested months | at control level | median stack cc at placement | median control p95 |
|---|---|---|---|---|
| CN.PGC | 56 | 26 | 0.088 | 0.088 |
| CN.VGZ | 54 | 21 | 0.105 | 0.101 |
| PB.B001 | 59 | 9 | 0.110 | 0.092 |
| PB.B003 | 58 | 16 | 0.098 | 0.092 |
| PB.B004 | 57 | 25 | 0.088 | 0.086 |
| PB.B005 | 60 | 17 | 0.095 | 0.086 |
| PB.B006 | 60 | 19 | 0.091 | 0.083 |
| PB.B007 | 60 | 17 | 0.099 | 0.087 |
| PB.B009 | 53 | 19 | 0.085 | 0.082 |
| PB.B010 | 53 | 21 | 0.088 | 0.082 |
| PB.B011 | 53 | 21 | 0.087 | 0.083 |
| PB.B013 | 23 | 0 | 0.535 | 0.087 |
| PB.B014 | 53 | 26 | 0.086 | 0.086 |
| PB.B017 | 59 | 39 | 0.088 | 0.090 |
| PB.B018 | 58 | 33 | 0.084 | 0.086 |
| PB.B926 | 51 | 30 | 0.078 | 0.079 |
| PB.B943 | 23 | 0 | 0.591 | 0.088 |
| UW.DOSE | 50 | 24 | 0.100 | 0.092 |
| UW.GNW | 60 | 15 | 0.102 | 0.085 |
| UW.HDW | 60 | 25 | 0.155 | 0.146 |
| UW.LRIV | 51 | 23 | 0.091 | 0.088 |

Swap / epoch-change checks (`w2/r2_swap_checks.csv`; carry-over = the old node's template used on the new node):

| station | date | change | cc before / ctl p95 | cc after / ctl p95 | best channel transform after | verdict |
|---|---|---|---|---|---|---|
| CN.PGC | 2015-01-12 | node change (WA_PROD_PLAN W2 / node tables) | 0.136 / 0.099 | 0.074 / 0.088 | (0, 1, 2)(1, -1, -1);(1, 2, 0)(1, 1, 1) | AT CONTROL LEVEL |
| CN.PGC | 2017-08-10 | channels  -> HHE,HHN,HHZ; sensor  -> Nanometrics Trillium 120 Sec PH R | 0.088 / 0.088 | 0.097 / 0.082 | (1, 0, 2)(-1, -1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| CN.VGZ | 2007-02-28 | sensor Guralp CMG3ESP/CNSN GD1 digitizer -> Guralp CMG3T_30sec/CNSN GD | nan / nan | 0.376 / 0.092 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| CN.VGZ | 2009-07-31 | channels  -> HHE,HHN,HHZ; sensor  -> Guralp CMG3ESP_60sec/CNSN GD2_v21 | 0.270 / 0.093 | 0.314 / 0.094 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| CN.VGZ | 2017-03-03 | sensor Guralp CMG3ESP_60sec/CNSN GD2_v21 digitizer -> Nanometrics Tril | 0.135 / 0.109 | 0.108 / 0.119 | (0, 2, 1)(1, -1, -1);(0, 1, 2)(1, 1, 1) | AT CONTROL LEVEL |
| PB.B013 | 2007-01-06 | channels  -> EH1,EH2,EHZ; sensor  -> HS-1-LT/Quanterra 330 Linear Phas | nan / nan | 0.710 / 0.092 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| PB.B014 | 2008-03-05 | channels  -> EH1,EH2,EHZ; sensor  -> HS-1-LT/Quanterra 330 Linear Phas | nan / nan | 0.112 / 0.075 | (0, 1, 2)(1, 1, 1);(0, 2, 1)(1, 1, 1) | above control (carry-over OK) |
| PB.B017 | 2007-07-12 | channels  -> EH1,EH2,EHZ; sensor  -> HS-1-LT/Quanterra 330 Linear Phas | nan / nan | 0.136 / 0.089 | (0, 1, 2)(-1, -1, -1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| PB.B926 | 2007-11-13 | channels  -> EH1,EH2,EHZ; sensor  ->  | nan / nan | 0.236 / 0.080 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| PB.B943 | 2008-02-27 | channels  -> EH1,EH2,EHZ; sensor  -> HS-1-LT/Quanterra 330 Linear Phas | nan / nan | 0.721 / 0.085 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.DOSE | 2010-07-21 | channels  -> BHE,BHN,BHZ; sensor  -> TRILLIUM 120PA,Velocity Sensor,NA | nan / nan | 0.176 / 0.085 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.DOSE | 2019-05-07 | channels BHE,BHN,BHZ -> HHE,HHN,HHZ | 0.098 / 0.092 | 0.087 / 0.089 | (0, 1, 2)(-1, 1, 1);(0, 2, 1)(-1, 1, 1) | AT CONTROL LEVEL |
| UW.GNW | 2008-04-18 | gain BHE x0.533;BHN x0.533;BHZ x0.533 | 0.102 / 0.086 | 0.195 / 0.083 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.GNW | 2010-01-29 | gain BHE x1.878;BHN x1.878;BHZ x1.878 | 0.095 / 0.086 | 0.391 / 0.083 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.GNW | 2010-12-01 | channels  -> BHE,BHN,BHZ; sensor  -> Velocity Transducer | 0.284 / 0.084 | 0.193 / 0.084 | (0, 1, 2)(1, 1, 1);(2, 0, 1)(-1, 1, 1) | above control (carry-over OK) |
| UW.GNW | 2019-05-07 | channels BHE,BHN,BHZ -> HHE,HHN,HHZ | 0.143 / 0.079 | 0.123 / 0.079 | (0, 1, 2)(1, 1, 1);(1, 2, 0)(1, -1, 1) | above control (carry-over OK) |
| UW.HDW | 2008-04-01 | gain EHZ x0.505 | 0.183 / 0.149 | 0.216 / 0.147 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.HDW | 2008-10-01 | channels EHZ -> EH2; sensor L-4C VERTICAL,Velocity Sensor,SERCEL -> S- | 0.229 / 0.148 | 0.162 / nan | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | untested |
| UW.HDW | 2020-12-04 | channels EH2,EHZ -> EHZ; sensor L-4C VERTICAL,Velocity Sensor,SERCEL;S | 0.150 / 0.142 | 0.129 / 0.145 | (0, 1, 2)(-1, 1, 1);(0, 1, 2)(-1, 1, 1) | AT CONTROL LEVEL |
| UW.HDW | 2022-09-01 | channels EH2,EHZ -> HHE,HHN,HHZ; sensor L-4C VERTICAL,Velocity Sensor, | 0.135 / 0.146 | 0.165 / 0.132 | (0, 1, 2)(-1, 1, 1);(1, 0, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.LRIV | 2009-08-27 | channels  -> BHE,BHN,BHZ; sensor  -> Velocity Transducer | nan / nan | 0.292 / 0.096 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.LRIV | 2019-05-30 | channels BHE,BHN,BHZ -> HHE,HHN,HHZ | 0.111 / 0.109 | 0.082 / 0.079 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(-1, -1, -1) | above control (carry-over OK) |
| PB.B004 | 2020 | gap year check | nan / nan | 0.090 / 0.080 | nan | above control (carry-over OK) |
| PB.B014 | 2023 | gap year check | nan / nan | 0.080 / 0.081 | nan | AT CONTROL LEVEL |
| PB.B014 | 2024 | gap year check | nan / nan | nan / nan | nan | untested |
| PB.B926 | 2025 | gap year check | nan / nan | 0.068 / 0.074 | nan | AT CONTROL LEVEL |

- Positive control (LRIV after 2019-05-30, HH node, old BH template carried): tested months 27, at control level 14, median cc 0.085 vs control p95 0.087; before 2019-05: 24 tested, 9 at control, median cc 0.100. Flagged: True. (SVI R2: LRIV 2021-10 0.17 vs control p95 0.09 (weak, above control); A-2 -0.216 on its own event set)
- Effective station set per month: `w2/r2_effective_sets.csv` (median 14 of 21, min 2 in 2022-02).
- Bad-for-dv/v windows (KMR style, `w2/r2_bad_for_dvv.csv`): 16779 windows; by reason {'line': 12075, 'noisy': 1455, 'month': 1248, 'flat': 625, 'dc_jump': 526, 'gap': 501, 'dead': 171, '2-8': 165, 'meta_change': 13}.
- Read with care: the stack test has little power where a station's templates are weak for the families active on the test days (e.g. SW carriers B017/B018 on NW-dominated days); a 'control_level' month there means 'no evidence the station sees these families', not 'broken'. Months between test days inherit the nearest tested month of the same instrument epoch.
- W0 re-hash at the end: 139 frozen input files, changed 2 (`w0/verify_20261003T0252.csv`).
- W1 final readability (full reads by W2 raw + SVI R2): of 131597 present raw station-days, 131590 read in full without error, 7 errors, 0 not fully read by either (`w1/w1_readability_final.csv`).
<!-- W2-AUTO-END -->
