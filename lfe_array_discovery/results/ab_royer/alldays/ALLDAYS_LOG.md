# All-days detection pass -- running log

User requirement (2026-09-25): detections are NOT gated to tremor/ETS windows. Growth stays literal
on tremor days. A final pass scans every day on which a final template's stations recorded.
Tremor-day and ETS labels are columns and are never used as filters.

Scripts (none of the chain's files are edited; the chain's modules are imported read-only):
- `scripts/ab_alldays_fetch.py`: probe, plan, fetch and 40 sps cache. It uses `ab_fetch.fetch_one` and
  `ab_cache40.cache_one` unchanged.
- `scripts/ab_alldays.py`: the detection pass (`run`), labels (`labels`) and assembly.
- `scripts/ab_alldays_report.py`: writes `REPORT.md`.
- `scripts/ab_alldays_queue.sh`: the detached follow-on queue.

## Rules (documented choices)
- **Array lifetime.** NW 2006-07-21..2012-05-31 and SW 2006-07-09..2012-05-31 (CAFE XU lifetime).
  The SW lines end on 2008-12-31. After that, only the permanent pool stations and XU PL11/W060 remain,
  so read SW 2009-2012 detections with `n_sta` in mind. NVI 2005-06-01..2007-08-31: the PO
  deployment, found by a monthly probe. A later PO segment from 2009-10 is excluded.
- **Station-days.** Each station's window is intersected with its FDSN epoch and a monthly probe
  (10 min on the 15th; a month is kept if it or a neighbouring month has data). CN.MAYB returned no
  probe data (EHZ-only, NRCan), so MAYB is absent from the NVI all-days pool.
- **Days scanned.** Every day with >= 4 cached template stations. The detector skips any day where
  fewer than 4 stations have data, as growth does.
- **Templates.** `tpl_it<final_iter>.npz` from `<SUB>/grow/summary.json`. All templates in the file
  are scanned, with the same station subset. `n_det_growth` is a column.
- **Null.** `<SUB>/null<K>/grow` final templates go through the same pass over the same real data
  (AB_NULL forced to 0), so the comparison is between families built on real data and families built
  on noise, on identical days.
- **Tremor-day labels.** PNSN from 2009-08-06: >= 24 epicentres within 25 km of the pool (the
  `ab_tremor_days` rule). Before that, our WECC day labels (`tremor/wecc_days_<SUB>.csv`). NVI uses the
  paper's two episodes only. All other days are UNKNOWN (NaN), never treated as non-tremor.
- **ETS windows.** Tremor days separated by <= 3 non-tremor days are merged into runs, and a run
  with >= 5 tremor days is a window. Where no day labels exist, the PNSN published windows are added.
  NVI uses the paper episodes.

## Tests
### Detector identity (2026-09-25 20:08-20:10Z)
`ab_alldays.mf_day` is `ab_grow.matched_filter` with the zero-padded template built on the GPU.
With `--verify`, both ran on the same days with the smoke templates
(`smoke/grow/tpl_it1.npz`, NW, 135 templates, 14 stations). **Identical on 5/5 days**: same
count, same families, same times, cc within 1e-6. Speed is about the same (about 10 s/day for
135 templates), so the gain is memory only. Runtime scales with the number of final templates,
roughly 10 s/day per 135 templates: about 2800 templates x 2142 NW days would take about 5 days of
GPU time per job.

### Smoke numbers (NW smoke templates; `smoke/`, 19,893 detections)
| date | label | phase | stations | detections | families | median nmad | nmad >= 12 |
|---|---|---|---|---|---|---|---|
| 2010-08-15 | PNSN tremor | during (2010-08-08..27 window) | 14 | 8202 | 127 | 10.4 | 2643 |
| 2010-08-16 | PNSN tremor | during | 14 | 5532 | 100 | 10.5 | 1931 |
| 2010-09-20 | non-tremor | post, 23 d | 12 | 3356 | 53 | 9.5 | 481 |
| 2009-12-01 | non-tremor | post, 103 d | 10 | 1504 | 42 | 9.1 | 153 |
| 2010-06-15 | non-tremor | pre, 53 d | 13 | 1299 | 61 | 9.2 | 133 |

Reading: 8 x MAD per 24 h lets through roughly 10 detections per template per day on quiet days by
construction. The information is in the rate and in the nmad distribution compared with the null,
which is why the report gives rates per scanned day and separation curves over nmad. Tremor days
give 4-6x more detections and far more high-nmad hits, as expected. These are two-day smoke
templates, so no science is claimed.

NW ETS windows from the labels (`NW/ets_windows.csv`): 2007-01-13..31, 2007-10-11..15,
2008-04-24..05-18, 2009-04-28..05-17 (WECC); 2009-08-14..19, 2010-03-16..21, 2010-08-08..27,
2011-04-27..05-06, 2011-07-30..08-27 (PNSN). 2006 is unlabelled.

### Smoke prefetch
NW 2009-12-01, 2010-06-15 and 2010-09-20: 55 station-days. 18 fetched, and 49 cached (37 built,
12 z-only).

## Launch (2026-09-25 20:11Z)
- Queue `ab_alldays_queue.sh`, pid 1981520 (`results/AB_ALLDAYS_PID`), with lock
  `locks/ab_alldays.lock`. Its heartbeat subshell (`results/AB_ALLDAYS_HEARTBEAT`) closes the lock fd and
  is killed in the EXIT trap.
- Stage P (prefetch) started first: 75,416 station-days in the plan (NW 38,858 / SW 32,335 /
  NVI 14,422). 43,667 to fetch, 25,382 raw on disk to cache, 5,555 already cached, 812 known nodata.
  Rate is about 3 fetches/s, about 4-5 h in total.
- Then it waits on a `sleep 600` loop for `<SUB>/grow/summary.json` (real) and
  `<SUB>/null<K>/grow/summary.json`. It runs real NW, SW, NVI first, then null1-3, and rewrites
  REPORT.md after each job. `AB_ALLDAYS_DONE` is written when all 12 are done or skipped, and
  `AB_ALLDAYS_FAILED` on any failure.
- Logs: `logs/ab_alldays_{queue,fetch,report}.log`, `logs/ab_alldays_<kind>_<SUB>.log`.
