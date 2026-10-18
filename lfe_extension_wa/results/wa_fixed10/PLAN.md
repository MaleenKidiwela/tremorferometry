# WA fixed-10 re-run — plan (written 2026-10-07, ~10:25 pm PDT Oct 6, before launch)

**Why.** The production core network changes 12 → 37 → 19 stations across the eras E1 2007–11, E2 2012–14 and E3 2015–26.
So a family's rate can step at an era boundary with no physical cause (NW net per family-tremor-day 0.48 / 0.35 / 0.28).
This run uses only stations that are live in all three eras, so the rates can be compared across eras. Spec: Merlin step 4,
`results/episode_rates/MERLIN_DESIGN_2026-10-06.md`. The user approved the run.

**Stations.** The fixed 10 are B001, B003, B004, B005, B006, B007, B013, B017, B018 and GNW. This is exactly the production E1 ∩ E2 ∩ E3
(`wa_prod/setup/core_sets.json`), and it matches Merlin's list.

## 1. Design

| item | fixed-10 run | production |
|---|---|---|
| output root | `results/wa_fixed10/` (own setup, mf, hf, conf, post, logs, markers, lock `locks/wa_f10.lock`) | `results/wa_prod/` (read only here) |
| code | the production scripts, unchanged in detection logic, run with `WA_PROD_DIR=results/wa_fixed10` | same |
| caches / templates / raw | production `cache20hz_prod` (+ older caches), `results/templates/tpl_*.npz`, raw roots: **read only** | same |
| templates | the frozen W4a library restricted to the 10 stations: 2,760 of 19,622 templates (tags build, t2012new, step4b_t2012/13/14) | W4a; W4b only in the rich pass |
| detector | identical (20 Hz 1–8 Hz cc, guard, drop-dead 5 %, site merge B, ≥ 2 sites, 8 MAD, declustering) | same |
| networks in one pass | **`core` = F9** (the 9 PB boreholes) in every era; **`all` = F10** (every loaded station = the literal fixed 10) | `core` = era set (12/37/19), `all` |
| series | data, scrambled control c1, time-reversed rev (truth: + c2) | same |
| span | 2007-01-01 → 2026-09-15, 7,198 days (`setup/days_prod.json`) | same |
| truth run | `truth`: nets core,all; data,c1,c2,rev on the 397 truth days (fit months, V2008/V2010, A2011, 2011-04 excluded, sealed 2007-01 / 2009-05) | `truth` (E3,all) + `truth_e12` (E1,E2) |
| confidence | `wa_prod_fit.py fit` (same recipe: HGB + isotonic, 5 day-folds, ≤ 600k fit rows, V-window build halves in-sample, guards out, tpl_insample rows in-sample), refitted on the fixed-network truth rows (core + all pooled); `open_sealed` one look | same recipe on E1/E2/E3/all rows |
| PNET | **not used (WAVE-only)**, see D2 | WAVE+PNET headline (gate adopted) |
| post | `wa_prod_post_v2.py all` incl. the k fix (§29) and the zero-day fix (§30), run once per network → `post_v2_core9/`, `post_v2_all10/`; then `wa_f10_strict.py` (station liveness and strict days) | `post_v2/` |

**Data flow.** setup copy → matched filter (prod in the background, 2 shards; truth in the foreground) → HF ratio for every
candidate (CPU, follows the day files) → HF on truth → fit → sealed one look → apply to every production candidate
(4 processes, split by month) → post_v2 for core, then strict; post_v2 for all, then strict.

## 2. Deviations from production (each with its reason)

- **D1. The fixed network is F9 (`core`); F10 is computed too (`all`).** GNW is a fixed station but not a fixed template set:
  - Sensor epochs: BH.50 to 2010-09-09, then an 85-day outage, then BH.40 from 2010-12-04, then HH.100 from 2019-05-08 (HH.100 carries BH.40's templates, approved by W2).
  - Templates: 260 families have a GNW BH.50 ('build') template, but only 52 also have a BH.40 one (t2012new / step4b, built in 2012–14). Those 52 are 49 NW, 2 Olympic A and 1 Arm B.
  - So for 208 families, GNW contributes only until 2010-09. Using it would put a template-set step inside E1 for most families, which is exactly the kind of step this run exists to remove.
  - The 9 boreholes each have one sensor node and one frozen 'build' template set for 2007–2026 (258–274 families each).
  - Both networks come from one matched-filter pass, because the station cc is computed once. Both get the full chain: fit, apply, post, strict.
  - The `all` network keeps the literal Merlin/user spec.
- **D2. No PNET.** Production PNET features are site-level maxima of the picker's characteristic function, computed on the fly per day from the cache. They are not stored per station, so they cannot be restricted to the 10 stations without re-running the picker.
  - So the confidence is the WAVE model only.
  - The headline rules keep their production form with P_own_WAVE in place of P_own_WAVE_PNET: NW / Arm B P_own_WAVE ≥ 0.75 on every day class; SW P_own_WAVE ≥ 0.5 on tremor days; Olympic A groups P_nbr3_WAVE ≥ 0.75 (as before).
  - For comparison with production, use production's WAVE columns: label `own`, cut 0.75 (`rate075`), or the `wave_classrules` view. Do not use its PNET headline.
  - A PNET block could be added later by re-running `wa_prod_pnet.py` (truth + prod) on this run's det files: the picker on 10 stations, ~10 GPU-hours, plus a refit.
- **D3. One k per network for every era** (`WA_POST_KNET=core|all`). The truth days all lie in 2007–2011, but the network is identical in every era by construction. So k comes from that network's truth fit rows and applies to E1, E2 and E3. In production, E1/E2 had no PNET k.
- **D4. Zero-day rule** (`WA_POST_SCANNED=famday`): a unit counts as scanned on a day when the family is observable (≥ 2 sites) in the network that day (mf famday tables). Production instead required ≥ 1 live template node in the core run. With only 9–10 stations, '1 live station' days are common and cannot detect anything (E1).
- **D5. Feature z_<STA>** is for the 9 boreholes only (the `E3` key of this run's core_sets.json = F9). Production used the 19 E3 stations. The model is refit, so this is consistent within the run.
- **D6. Sealed months are no longer blind.** 2007-01 and 2009-05 were opened once for the production model (2026-10-04). Here they stay out of the fit and are scored once (`conf/SEALED_OPENED.json` written first), so they are out-of-sample but not blind.
- **D7. Code changes** (each default is byte-identical; backups `*.bak_pre_fixed10`):
  - `wa_prod_common.py`: `WA_PROD_LOG` sets the log path.
  - `wa_prod_post_v2.py`: `WA_POST_NET`, `WA_POST_V2_DIR`, `WA_POST_HEAD_LAB`, `WA_POST_KNET`, `WA_POST_SCANNED`; it also tolerates a missing PNET column.
  - `wa_prod_fit.py apply`: `--shard/--nshard` over months; the done marker gets a `_s<i>` suffix only when sharded.
  - Default check: `k_factors_v2()` with default environment reproduces production `post_v2/rates_k_factors_v2.csv` exactly (6,660 rows, max |Δ| 0).
  - New files: `scripts/wa_f10_setup.py`, `wa_f10_strict.py`, `wa_f10_queue.sh`. `wa_prod_mf.py` and `wa_prod_hf.py` are not changed.
- **D8. Not run:** rich pass (none of the 10 has a usable W4b template; GNW's 25 W4b candidates all failed QC), refit-risk, cut, delete list, report.

## 3. Station liveness and outages (decision)

Live = the station was loaded and passed the guards (mf run members). This is identical between production core and this run, because
it is decided per station. Share of days, from the production run files:

| station | E1 2007–11 | E2 2012–14 | E3 2015–26 |
|---|---|---|---|
| B001 | 0.923 | 0.994 | 0.990 |
| B003 | 0.878 | 0.976 | 0.927 |
| B004 | 0.925 | 0.996 | 0.929 |
| B005 | 0.853 | 1.000 | 0.973 |
| B006 | 0.917 | 0.996 | 0.975 |
| B007 | 0.863 | 0.993 | 0.991 |
| B013 | 0.868 | 0.975 | 0.980 |
| B017 | 0.804 | 0.993 | 0.986 |
| B018 | 0.790 | 0.995 | 0.983 |
| GNW | 0.784 | 0.974 | 0.977 |

| | E1 | E2 | E3 |
|---|---|---|---|
| all 9 boreholes live | 0.465 | 0.924 | 0.776 |
| ≥ 8 boreholes live | 0.748 | 0.995 | 0.977 |
| mean boreholes live | 7.8 | 8.9 | 8.7 |

Days with all 9 boreholes live, by year: 2007 0.11, 2008 0.33, 2009 0.56, 2010 0.75, 2011 0.58, 2012–19 0.85–0.94, **2020 0.36**,
2021 0.86, 2022 0.73, 2023 0.69, 2024 0.60, 2025 0.95, 2026 0.60.

Outages ≥ 30 days (`setup_check/outages_ge30d.csv`):

| station | outages |
|---|---|
| B017 | 2007-01-01..07-17 (198 d; node starts 2007-07-12); 2007-12..2008-01 (32 d); 2009-06..07 (37 d); 2011-07-20..09-10 (53 d) |
| B018 | 2007-09-04..11-28 (86 d); **2007-12-22..2008-05-13 (144 d, covers the V2008 ETS build window)**; 2011-07-20..09-10 (53 d) |
| B004 | **2020-03-12..10-14 (217 d)**; 2011-07/08 (35 d); 2011-11/12 (30 d) |
| B007 | 2007-11-06..2008-02-08 (95 d); 2011-02/03 (38 d) |
| B013 | 2007-02-23..05-19 (86 d); 2022-10/11 (30 d) |
| B005 | 2008-05-29..08-18 (82 d); 2007-11/12 (38 d); 2023-10/11 (31 d); 2024-02/04 (35 d) |
| B006 | 2007-11/12 (47 d); 2008-01/02 (34 d); 2022-12-22..2023-03-01 (70 d) |
| B003 | 2026-06-12..08-28 (78 d) |
| B001 | 2007-11/12 (37 d) |
| GNW | 2010-05/06 (34 d); 2010-09-10..12-03 (85 d, sensor change); 2023-10-17..12-07 (52 d) |

**Decision: report both.**
- (a) **All scanned unit-days** (post_v2 tables as they are). The network on a day is the live subset of F9 / F10, as in production.
- (b) **Strict unit-days only**: every station of the family's fixed template set contributes that day. For core, that is the family's borehole template stations. For all, it is the same plus GNW whenever GNW's sensor epoch on that date has a template for the family. An Olympic A group-day is strict when every observable member family is.
  - Outputs: `post_v2_<net>/unit_day_strict.parquet` (n_fixed, n_fixed_live, strict per unit-day) and `post_v2_<net>/strict/rates_summary_v2.csv` (same function and columns as post_v2's summary). Also `strict_share_by_region_era.csv`, `station_liveness_by_era.csv` and `station_outages_ge30d.csv`.
- Strict is the clean cross-era comparison, but E1 strict coverage is thin: 2007 ~11 %, 2008 ~33 % of days, and the 2008 ETS lacks B018. Mixing the two would let E1 rates depend on which boreholes are down.
- The per-episode step (later) should use strict days, or n_fixed_live as a column.

## 4. Smoke tests (done before launch; `smoke/`)

1. **Exact reproduction where the networks are identical.** These are days on which the production core network was a subset of the 10:
   - 2008-04-01 and 2008-07-01 (8 boreholes, core = F9 = production core);
   - 2008-04-01 (`all`) and 2011-01-04 (`all` = 8 boreholes + GNW BH.40 = production core).
   - Every data / scrambled-control / time-reversed row is identical: same positions, 0.0 difference over all 31 shared float columns (net_mad, net_cc, site features, z, HF position …). `smoke/compare_prod_f10.csv`.
2. **Plumbing reproduction on E3 days.** A scratch setup with the restricted-pool / restricted-library method, but with production's E3 set (19 stations):
   - 2016-03-01 and 2020-06-01 (GNW HH.100 carry-over) reproduce production core exactly (9,819 / 9,802 / 2,603 and 3,159 / 3,148 / 15,565 rows; 0.0 difference over 41 float columns). `smoke/compare_prod_e3.csv`.
3. **Full chain** (16 truth days incl. sealed, 7 production days across the eras; fit ≤ 200k rows): setup → mf → HF → fit → sealed → apply → post core / all → strict all ran.
   - Run twice (`smoke/chain`: 1 shard, single apply; `smoke/chain2`: 2 shards, 3 apply processes).
   - Identical P values and rate tables in both (max |Δ| 0).
   - k factors exist for every era from the one network.
4. **Default check:** post_v2 `k_factors_v2()` with default environment = production table exactly.

**Speed and memory.** The matched filter is **GPU-bound** (85–100 % utilisation):
- 5.5 s/day alone; 11–13 s/day each with 2 concurrent processes; ~20 GB GPU per process (95 GB card).
- CPU RSS ~1.9 GB per mf process.
- HF ratio ~25 worker-s per day.
- apply ~5 s per 35k rows per process.

Side note: on 2020-06-01, production itself has time-reversed control > data (core rev 15,565 vs data 3,159). The same happens
in this run (B004 out, 2020 outage period). This is not new; flagged for the write-up.

## 5. Runtime estimate (Pacific time; launch ~10:30 pm PDT Tue Oct 6)

| stage | estimate | finish (PDT) |
|---|---|---|
| matched filter, prod 7,198 d + truth 397 d (GPU-bound, ~5.6 s/day aggregate) | ~12 h | ~10:30 am Wed Oct 7 |
| HF prod (16 workers, follows the day files) | keeps pace | ~11 am |
| truth HF + fit + sealed (during the matched filter) | ~1.5 h | ~1 am |
| apply, ~250 M rows, 4 processes | ~2.5–3 h | ~2 pm |
| post_v2 + strict, core then all | ~3–4 h | **~5–7 pm Wed Oct 7** |

## 6. Run, check, resume

- Launch: `nohup setsid -f bash scripts/wa_f10_queue.sh > results/wa_fixed10/logs/queue.log 2>&1` (defaults: 2 matched-filter shards, 16 HF workers, 4 apply processes).
- Progress:
  - `results/wa_fixed10/WA_F10_HEARTBEAT` (stage, time, GPU);
  - `results/wa_fixed10/gpu_mem.log` (every 2 min: GPU, number of prod days done, HF days done);
  - `logs/queue.log`, `logs/wa_f10.log`, `logs/bg_prodmf.log`, `logs/prodmf_<i>.log`, `logs/bg_hfprod.log`, `logs/post_<net>.log`.
- Markers:
  - stages: SETUP_DONE, TRUTH_MF_DONE, HF_TRUTH_DONE, FIT_DONE, SEALED_DONE, PRODMF_DONE (+ mf/prod/MF_DONE), BG_hfprod_DONE, APPLY_PROD_DONE, POST_core_DONE, POST_all_DONE;
  - end: `WA_F10_DONE`; failure: `WA_F10_FAILED` (stage + rc).
- Resume: relaunch the same command. Finished stages are skipped. The matched filter skips days whose run_<date>.json exists, HF and apply skip existing day / month files, and the sealed look refuses a second run.
- Kill: by PID only (`WA_F10_PID` = queue shell; children via `pgrep -P`), never `pkill -f`.

## Run log
- 2026-10-07 12:38 am PT: queue FAILED at hf_truth: empty detection days (no observable family; wa_prod_mf writes det_<date>.parquet with only fid, t) — 1 truth day (2007-10-31) and 16 production days so far; production never had such days. Fix: wa_prod_hf.day_job writes an empty hf file and status "empty" (not counted by the guard) when the det file has no k column (backup wa_prod_hf.py.bak_pre_emptyday; tested on 2007-10-31 and a normal day). The other readers (wa_prod_fit load_run) already skip empty days. Relaunched ~12:55 am PT; finished stages resume.
