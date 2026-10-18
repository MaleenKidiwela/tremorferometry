# P4-2b — burst-robust thresholds, era thresholds, flag columns, completeness ladder

Prescription: `../PLAN.md`, section "P4-2b — Merlin review 2" (adopted 2026-10-11, user said "keep going").
Scripts `svi_science/scripts/bw_p4b_*.py` / `.sh`; outputs in this folder. Rules below written
**Sat 2026-10-10, 8:20 pm PT, before anything was computed**. Nothing here edits a forward script or a production output;
flag columns are sidecar tables. No Bostock events are read (step check 7d belongs to P6's one look).
"Scrambled control" = the same detector on the same day with each station site's data shifted to another time of that day.

## Part A — rules

### Copied from PLAN.md (P4-2b, adopted)

> Root cause of "2003–05 thresholds ~10× too strict": 25 PGC burst days (18 in Feb 2015, after the 2015-01-12 PGC
> change) carry 85 % of the CN4 scrambled-control candidates in the P4-2 reference days. The pooled-mean rule is
> hijacked by them. With burst-robust statistics every frozen set lands at 8.75 ± 0.15 (median-matched CN4 8.80, B1 8.70,
> B2 8.85, B3 8.80). Real 2003–05 CN4 control at 8.75 is 0.076–0.106, the same as S15. The 2005-11 "step" is the same artefact.
> Second flaw (structural): in the frozen detector the merged PGC site (PGC+B009/B010/B011) weighs ~2× any other site and
> is shifted as one unit in the control. So for CN4/B1/B2/B3 the scrambled control is a one-site-signal model, and
> data − control is biased low (more so on tremor days).
> Prescription:
> 1. P4-2b: same 509 days. (i) per-day median control rate; (ii) pooled mean after dropping the top-5 % S15fix control
>    days. Reference S15fix at 8.75. Pass: every frozen set within [8.50, 9.00] by both → T_core 8.75 all frozen sets,
>    T_ext 8.00. Statistics disagreeing > 0.5 → bracket.
> 2. Era detector: same statistics on its env-quiet days (quiet_env3_CN4), matched to the S15fix levels (expected ~8.25).
>    Cross-check frozen CN4 on the same days → 8.75 ± 0.25, else bracket. Use the control, never data FA.
> 3. Burst-day column (control ≥ 5× yearly median); flag, never drop.
> 4. Single-site column n_z2 ≥ 2. B2/B3 between rates = bracket [all, n_z2 ≥ 2]. Per-site-cap test logged as a decision
>    point after P6, on the 2005-11..2006-06 overlap.
> 5. Completeness ladder: old set run 60 d after each boundary (2005-07-01, 2005-11-02, 2007-11-13, 2009-08-27); k from
>    env-quiet ΣP-excess with bootstrap CI; k_CN4/era etc. from the overlap, ETS vs between separately; agree within ×2,
>    else bracket.
> 6. ΣP science threshold-free: keep all ≥ 8 MAD; rate = (ΣP_data − ΣP_control)/observable family-days; T only for the
>    headline count; never T before P.
> 7. Step checks before any 2003–13 vs 2014–26 claim: (a) 1–2 pass; (b) no step > 2σ after k in ±60 d; (c) slope of excess
>    on control ≈ 0 (env-quiet days); (d) Bostock ETS recall at P ≥ 0.75 agrees with ladder k within ×1.5; (e) 2010–13
>    between rates on env-quiet days (N ≤ 12).
> Keep the pooled-mean values beside the robust ones; this is a recorded protocol change.

### How the statistics are computed here (pinned before computing)

Common to steps 1 and 2 (the G1 control rule, as in P4-2): control candidates = `series == "c1"` rows of the feature files
with `net_mad ≥ T`, **no HF veto**; control family-days = `series == "c1"` rows of the `_fd` files with `observable`.
Per UTC day d: r_d(T) = n_c1(d, T) / m_c1(d). Days with m_c1(d) = 0 are left out. Data candidates are never used for a
threshold. T grid 8.00 … 14.00 step 0.01; matched T = the smallest grid T whose statistic is ≤ the target; at 8.00 it is
reported as "≤ 8.00 (floor)". Each matched T is also given rounded up to the 0.25 grid.

- **MED** (statistic i): median over days of r_d(T).
- **TRIM** (statistic ii): Σ n_c1(d, T) / Σ m_c1(d) over the days left after dropping the top-5 % days, i.e. the
  ⌈0.05 N⌉ days with the highest reference control rate r_d(8.75) (E0: the S15fix run; 509 days → 26 dropped).
  The same dropped days apply to every set (matching on identical days).
- **POOLED** (the P4-2 value, kept beside): Σ n_c1 / Σ m_c1 over all days (no trim).
- Core target = the reference's statistic at 8.75; ext target = at 8.25 (reported only; the rule adopts T_ext 8.00).

### Step 1 — P4-2b on the 509 E0 PNSN-quiet days
- Runs `p4/mf/feat/p4e0_{CN4,B1,B2,B3,S15fix}` (P4-2, unchanged), days `p4/tables/p4_e0_quiet_days.json`. Reference S15fix.
- **Pass** = for each of CN4, B1, B2, B3: MED-matched T_core and TRIM-matched T_core (0.01 grid) both in [8.50, 9.00].
  Pass → T_core 8.75 for all frozen sets (CN4, B1, B2, B3, S15), T_ext 8.00.
- **Bracket** = for any set where |T_MED − T_TRIM| > 0.5: that set's T_core is the bracket [min, max] of the two.
- A set that fails the window without a > 0.5 disagreement is not adopted: its MED, TRIM and POOLED values are carried
  and the choice goes back to the user / Merlin (logged in PLAN.md).
- Verification of the "25 burst days" claim (reported, not a pass rule): rank the E0 days by the CN4 control count ≥ 8.75;
  share of all CN4 control candidates ≥ 8.75 (and ≥ 8.00) on the top 25 days; how many fall in Feb 2015 and after
  2015-01-12; on those days the fraction of CN4 control candidates whose largest site z is the PGC site; overlap of the
  25 days with the 26 S15fix-trimmed days. Also the 2003–05 CN4 control rate at 8.75 from `p5f_2003..2005`
  (CN4-only days), pooled per year.

### Step 2 — era detector on its env-quiet days
- Days: `quiet_env3_CN4 == True` in `p4/envflag/tables/quiet_flag_env.csv`, 2003-02-01..2006-06-30, present in both the era
  run (`p5e_2003_kelb`, `p5e_2004`, `p5e_2005`, `p5e_2006`) and the frozen run (`p5f_2003..2006`).
- Targets: the S15fix E0 levels from step 1 (MED and TRIM at 8.75; ext at 8.25 reported).
- TRIM on these days: drop the ⌈5 %⌉ days with the highest **frozen-detector (p5f)** control rate r_d(8.75) — the
  analogue of trimming by the frozen S15fix at E0, and the same dropped days for era and frozen. Diagnostic beside it:
  trimmed by the era detector's own control rate.
- Era T_core (MED, TRIM, POOLED) on the full 2003-02..2006-06 span (adopted value) and per year (reported).
- **Cross-check:** the frozen detector on the env-quiet days of 2003-02-01..2005-06-30 where its stations ⊆ CN4
  ("frozen CN4"), same statistics, same targets, with the era detector on exactly those days beside it. Pass = frozen
  CN4 T_core within 8.75 ± 0.25 by the statistic. Pass → the era T (full span) stands for that statistic. Fail → bracket
  [era T, era T + (8.75 − frozen-CN4 T)] for that statistic. MED and TRIM era values differing by > 0.5 → bracket of the two.
- Needs ≥ 200 control family-days (always the case here; reported).

### Step 3 — flag columns (sidecar tables, the feature files are not rewritten)
- For every finished P5 run (its year marker `P5F_<Y>_DONE` / `P5E_<Y>_DONE` / `P5E_2003_KELB_DONE`; the superseded
  `p5e_2003` included) and for each ladder run when it is done. Runs not yet finished (p5f_2008, p5f_2009 window-only)
  get their sidecars from the queue when their year marker appears.
- Per run and UTC day: r_d = all stored control candidates (≥ 8.00 MAD) / observable control family-days.
  Reference = median of r_d over that run's days of the same calendar year with m_c1 > 0 (if that median is 0, the mean,
  recorded in `ref_stat`). **`burst_day` = r_d ≥ 5 × reference.** The flag applies to the data and control rows of that day.
- **`multi_site_z2` = n_z2 ≥ 2** (≥ 2 sites with z ≥ 2 at the candidate; `n_z2` copied).
- Key: `run, date, row` (0-based row of `<date>.parquet`), with `series, fam, t, net_mad` copied for checking.
  Files `p4b/sidecar/<run>_flags.parquet` (candidates) and `<run>_days.csv` (per day).

### Step 5 — completeness ladder (GPU queue; k waits for P6)
- Old station set, frozen detector, backward mode, exactly the p5f configuration (`bw_gate_mf.py --backward`, default
  site list, no prep, same caches) except `--sset eff --excl-file` = every S15 station outside the old set, each month
  touched (the P4-2 exclusion-file mechanism). One GPU process, one shard; scrambled control included (the detector always
  runs it).
- Boundaries (old set → new set; 60 days from the boundary, inclusive): CN4, 2005-07-01..2005-08-29; B1,
  2005-11-02..2005-12-31; B2, 2007-11-13..2008-01-11; B3, 2009-08-27..2009-10-25. Run names `lad_<OLD>_<YYYYMMDD>`.
- A day runs only when `P2_<YYYY>_DONE` exists for its year (2007 days now; 2008 / 2009 days wait). Markers
  `p4b/LAD_<OLD>_DONE` / `_FAILED`; queue `bw_p4b_q.sh`, PID in `p4b/Q_PID`.
- k (script `bw_p4b_kladder.py`, run after P6 gives P): k = excess rate of the new set (p5f, same days) / excess rate of the
  old set, excess rate = (ΣP_data − ΣP_c1) / observable data family-days, on env-quiet days (`quiet_env3_CN4`); ETS-window
  days and the rest reported separately; 2000-day bootstrap over days (paired), 68 % and 95 % CI. Overlap ks
  (era/frozen, 2003-02..2006-06, per frozen set) the same way, ETS vs between separately. Two k estimates of the same
  ratio agree if within ×2, else bracket. Until P exists the script only runs in `--smoke` mode (P ≡ 1, mechanics test;
  not a result).

## Part B — results (Sat 2026-10-10, 8:25 pm PT)

**Check progress (no waiting):** `cat p4b/logs/q.log`; `ls p4b/LAD_* p4b/Q_DONE`; queue PID in `p4b/Q_PID` (966687).

### Step 1 — E0, 509 PNSN-quiet days (`tables/p4b_e0_thresholds.csv`, `p4b_e0_decision.json`, curves `p4b_e0_curves.csv`)
Targets (S15fix at 8.75): MED 0.0348, TRIM 0.0562 (26 days dropped), POOLED 0.0909 control candidates per family-day.

| set | T_core MED | T_core TRIM | T_core POOLED (= P4-2) | MED vs TRIM | rule | diag: TRIM by the set's own top-5 % days |
|---|---|---|---|---|---|---|
| CN4 | **8.79** | 9.82 | 10.76 | 1.03 | TRIM outside [8.50, 9.00]; > 0.5 apart → **bracket [8.79, 9.82]** | 9.03 |
| B1 | **8.69** | 9.56 | 10.38 | 0.87 | same → **bracket [8.69, 9.56]** | 8.88 |
| B2 | **8.84** | **8.90** | 8.91 | 0.06 | both in window | 8.91 |
| B3 | **8.76** | **8.77** | 8.79 | 0.01 | both in window | 8.77 |
| S15fix | 8.72 | 8.75 | 8.75 | — | reference | 8.75 |

- **Step 1 FAILS as written**: MED puts every set at 8.69–8.84 (Merlin's 8.80 / 8.70 / 8.85 / 8.80 reproduced), but TRIM
  leaves CN4 and B1 high. T_ext (ext targets at 8.25): MED 8.16–8.35, TRIM 8.28–8.91.
- Why TRIM fails for CN4/B1: the 26 days dropped by the S15fix ranking contain only 15 of CN4's 25 burst days. The
  other 10 stay in and still dominate the pooled mean. Trimming each set by its own top-5 % days (diagnostic, not the
  registered rule) gives CN4 9.03 and B1 8.88.
- Burst-day claim **verified**: the top 25 CN4 control days hold 85.4 % of CN4's control candidates ≥ 8.75 (73.9 % of all
  ≥ 8.00). 18 of the 25 are in Feb 2015, and 23 are after 2015-01-12. On those days the PGC site has the largest z
  for 87 % of control candidates (71 % on other days), with median dom_top 0.91 (0.72 on other days).
- 2003–05 real CN4 control at 8.75 (p5f, CN4-only days): pooled 0.119 (2003), 0.086 (2004), 0.076 (2005 H1).
  Merlin's range was 0.076–0.106, and 2003 is a little higher. The medians are 0.027–0.036, against S15fix E0
  (0.091 pooled, 0.035 median). So on real data CN4 behaves like S15 at 8.75.

### Step 2 — era detector, env-quiet days (`tables/p4b_era_thresholds.csv`, `p4b_era_decision.json`)
- 475 `quiet_env3_CN4` days 2003-02..2006-06 with both runs; 24 dropped for TRIM (ranked by the frozen p5f control).
- **Era T_core: MED 8.25, TRIM 8.24 (→ 8.25), POOLED 8.07.** T_ext: 8.00 (floor) by both. Per year MED 8.20–8.31,
  TRIM 8.13–8.32. Diagnostic: TRIM ranked by the era detector's own control gives 8.15.
- **Cross-check PASS:** frozen CN4 on its 317 env-quiet days (2003-02..2005-06, stations ⊆ CN4) gives MED 8.57,
  TRIM 8.57, POOLED 8.58. That is within 8.75 ± 0.25 (−0.18), so the era T stands at **8.25** (no bracket).
  The era detector on the same 317 days gives MED 8.25, TRIM 8.20.

### Step 3 — sidecar flags (`sidecar/<run>_flags.parquet`, `_days.csv`, `_meta.json`)
- Written for p5f_2003–2007, p5f_2010–2013, p5e_2003, p5e_2003_kelb, p5e_2004–2006 and lad_CN4_20050701. The queue
  adds the other ladders, then p5f_2008/2009 when their year markers appear. Keys checked against the feature files.
- Burst days per run-year: 12–37, with a median control-rate ratio of 6–10× the yearly median. Between 8 and 19 of each
  run's burst days are envelope-active (tremor) days: on tremor days the scrambled control tracks the data, so the rule
  flags them too. It is a flag, never a drop.
- `multi_site_z2` holds for 72–84 % of frozen data candidates and 53–70 % of frozen control candidates. For the era
  detector the figures are 89–99 % (data) and 59–85 % (control).

### Step 5 — completeness ladder (queue `bw_p4b_q.sh`, PID 966687)
- `lad_CN4_20050701` done 8:21 pm PT (4 stations every day, R2-free exclusion file as pinned); `lad_B1_20051102`
  running (~5 s/day, ~8:27 pm); `lad_B2_20071113` 2007 part right after; its 2008 part waits for `P2_2008_DONE`
  (NRCan is at 2008-10, so roughly Sun 10-11); `lad_B3_20090827` waits for `P2_2009_DONE` (~Sun 10-11 late / Mon 10-12
  early PT, the P1 ETA). Each ladder takes ~5 min of GPU.
- `bw_p4b_kladder.py` is ready and waits for P6. Give it `--p-dir` with `<run>_P.parquet` (run, date, row, P). A smoke
  test (`--smoke`, P ≡ 1, `tables/k_smoke.csv`) ran end to end on the CN4 ladder and the three overlaps. It tests the
  mechanics only and is **not a result**.
