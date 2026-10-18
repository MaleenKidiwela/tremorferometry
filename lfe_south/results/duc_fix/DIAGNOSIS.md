# DUC_FIX diagnosis: why the 2008-2011 scan fails his cc x nchannel rule (2026-09-30)

## Plain-language summary

* **Cause: most of the station data never reached the detector.** When the continuous scan prepared each
  station-day, a GPU speed-up inside the response-removal step ran out of its 1.5 GB memory allowance. The
  step then failed without stopping the scan. **53,329 of 77,820 station-days (69 %) in 2008-01..2012-07 were
  lost this way.** Every 100 Hz and 50 Hz NC/PB station lost 80-95 % of its days; the 40 Hz BK broadbands
  (GASB, HOPS, JCC, WDC, YBH) mostly survived.
* **Why April 2008 worked.** The known-answer window (2008-04-01..05-10) did not go through this step. It reused
  the station-days already built by the earlier t6 chain (plain CPU). So the detector, templates,
  thresholds and variant ('nomute') were all correct. Only the new data built for the other months was broken.
* **How this breaks his rule.** A family normally sees about 15 channels. In the damaged months it saw about 3
  to 5. The mean correlation of a detection barely changed (median 0.129 in 2008-04, 0.140 in 2009-06). But his rule
  needs the cc summed over channels (cc x nchannel) to reach about 1.5. With 3 channels, detections summed to 0.40
  (median), so almost none passed. On 3 channels the hourly 8 x MAD also triggers on more noise. That is why the
  raw counts were high (2-3.5x his, for his 60 kept families) while the rule counts were 0.05-0.5x his.
* **Fix applied (a bug fix only; no change to his method).** `scripts/duc_fix_data.py` `_patch_fft`: the cuFFT
  plan cache is turned off, and any GPU failure falls back to numpy. The output is unchanged: the CPU path is
  bit-identical (max diff 0) to the known-answer cache, and the GPU path matches it to 1e-10..1e-9.
  The broken station-days are being rebuilt, followed by prep and matched filter for 2008-01..2012-07 (queue below).
* **Re-templating (0 kept, 26 dropped) is an artefact, not a result.** (1) The new station's own data were mostly
  lost (median 0 usable windows per family-channel). (2) The candidate LFEs came from the collapsed catalogue
  (median 3 per family). (3) The partial scan ended at 2012-07, and the queue then force-finalised 26 epochs whose
  24-month template window was still open. That gave them at most 9 months of detections. The other 17 new epochs
  start after 2012-07 and stayed pending. The force-finalise is now fixed in the queue (see §5).

## 1. Evidence

**Cache status, all planned station-days 2008-01..2012-07** (`data/plan_<m>.csv`, originals kept in `pre_diag/`):
error 53,329 / built 21,303 / reused 2,242 / no source 946. Of the errors, **99.98 %** read
`Out of memory allocating ... (limit set to: 1,61x,xxx,xxx bytes)`, which is the cupy pool limit set in
`_patch_fft`; the other 9 are "No matching response". Error share by station: B039 1434/1674, B040 1535, KCPB 1582,
KHBB 1548, L* SHZ stations ~1150-1300 each; GASB 45, HOPS 100, WDC 29, YBH 40, JCC 234.

**Reproduced in isolation** (fresh process, the unmodified function, B039 2009-06-01..12). The pool grew
277 -> 914 -> 1,470 MB as cuFFT plans accumulated (2 -> 11 plans, one per distinct FFT length).
`free_all_blocks()` cannot release plan work areas. The 2009-06-09 day (a different length) then failed with OOM.
With the fix, the same days all build (pool 0 MB after each call).

**Same data path, same numbers.** For 2008-04-15, B039 / KBN / KHBB / GASB were rebuilt through `duc_fix_data.cache_one`.
With the CPU path the result equals the t6 cache the known answer used (max |diff| = 0). With the GPU path it matches
to 6e-11..9e-11 (1.3e-9 for B039 2009-06-09). So the production path equals the known-answer path whenever it does not crash.

**Channels per detection (mf1 diag + parquet):**

| run | mean templates per family-day | median nchannel of detections (60 kept families) | median cc | median cc x nch | pass his rule |
|:--|--:|--:|--:|--:|--:|
| known answer 2008-04 (t6 cache) | 14.7 | 15 | 0.129 | 1.89 | 4,662 |
| production 2008-04 (same cache) | 14.67 | 15 | 0.129 | 1.89 | 4,662 |
| production 2009-06 (damaged)    | 4.6  | 3  | 0.140 | 0.40 | 26 |
| other damaged months (range)    | 4.0-12.0 | 3-4 | | | 0-0.1 per family-day |

His catalogue, for the detections that pass his rule, has median nchannel 17 (2008), 19 (2009), 18 (2010) and 16 (2011).
Across all his detections of the 60 kept families the median is 14 / 16 / 15 / 10. Ours was 4 / 3 / 3 / 3.

**Hypotheses ruled out.**
* (b) Per-period or per-network thresholds. His analysis scripts apply one `threshold_perm` per family to the whole
  2004-2011 catalogue (`cumulative_LFEs.py:42`: `cc * nchannel >= threshold_perm`). `threshold_cc.txt` has only
  FAME and perm columns. His `cc` is the channel MEAN, so cc x nch is the channel sum, exactly as `duc_fix_mf.py`
  computes it. The share of his detections that pass is steady at 0.24-0.52 per year.
* (c) Variant or configuration. The scan used `nomute` (`detector_choice.json`; the queue exports `DUC_MUTE=0`), and the
  same station list, template set and code as the known answer. The only difference is the data cache.
* (d) Time base. The joins are monthly counts, so no time matching is involved. The known answer matches one-to-one
  within +-2 s on his first-sample convention.
* (e) His 2011 raw count of 253,093. Of these, 204,087 are his ELIMINATED family 080326.07.048
  (threshold 0). From Dec 2010 to Jul 2011 it runs at nchannel = 3 with cc ~0.11, a 3-channel noise run (7 channels
  and cc 0.18 before and after). Another 22,113 are 080429.15.005, also eliminated. His raw count for his 60 kept
  families in 2011 is 24,309 (ours then 44,352: the same 3-channel noise excess as the other years).
  So the 0.21 raw ratio comes from his noise family, not from our detector.

## 2. Fix

`scripts/duc_fix_data.py::_patch_fft` (the only change to the data path):
`cp.fft.config.get_plan_cache().set_size(0)`, and every GPU rfft/irfft is wrapped with a numpy fallback on any exception.
Expected effect: every station-day that has source data builds again. Family-days should go back to about 13-15
templates, the same as the known answer, and the counts that pass his rule should return close to his
(the known answer gives ours/his = 4,734 / 4,717 under his rule).

## 3. Re-run (launch-and-exit queue)

`scripts/duc_diag_rerun.sh` (log `logs/duc_diag_rerun.log`, `logs/duc_diag_cache.log`; markers
`results/DUC_DIAG_{HEARTBEAT,DONE,FAILED}`; flocks `locks/duc_diag.lock` + `locks/duc_fix.lock`).
Nothing is deleted. The damaged `prep/<day>.npz` (1,634 days), `mf1/<m>` (54 months) and the re-templating outputs
were MOVED to `prep_pre_diag/`, `mf1_pre_diag/` and `pre_diag/`. 2008-04 was not damaged and is untouched.
The queue then: re-caches only the errored station-days (7 workers), re-preps those days (1 worker) and re-runs
the MF (GPU) per month, 2009-06 first as the test month. After that it re-runs re-templating `--upto 2012-07`
WITHOUT `--final`, plus a forced diagnostic in `diag_retpl/`; runs pass 2 if any template is kept; then
`duc_fix_report.py` and `duc_diag.py summary`, which appends the new joins below.
**Storage:** no downloads (all sources are already on disk). The rebuilt cache adds about 14 MB per station-day,
so roughly **+0.75 TB to `cache_duc/`** (366 GB now); the cache can be regenerated from the sources.

## 3b. Test month 2009-06 after the fix (the queue's first month, 23:28Z)

| 2009-06 | templates per family-day | raw 8 x MAD (all / 60 kept fams) | median nch (60 fams) | pass his rule |
|:--|--:|--:|--:|--:|
| before (OOM-damaged) | 4.60 | 11,311 / 10,728 | 3 | 26 |
| **after fix** | **16.06** | 6,346 / 2,824 | 16 | **1,723** |
| his catalogue | | 6,043 / 2,770 | 17 (rule-passing) | 1,715 |

One-to-one within 2 s, his rule on both sides: **recall 0.964 (1,653 / 1,715), precision 0.959**, which matches the April 2008
known answer (0.974 / 0.970). Ours/his under his rule: 0.015 before, **1.005 after**. The raw excess is gone too.
The remaining 53 months are running in the queue (about 10-11 min per month, so about 10 h). The full joins are appended below when it finishes.

## 4. Re-templating outcome explained

`templates_new.csv` (moved to `pre_diag/`) had 277 family-channels over 26 epochs. Median candidates per family
was 3 and median usable windows N was 0. Only 3 family-channels had N >= 20, and 1 of those also had >= 10 scrambled
controls. The one epoch with real numbers, GTC#2 (68 windows, family 080326.07.048), was built from detections of
an eliminated noise family. Its stack matched his SHZ template at cc 0.19 and did not beat its controls.
So "0 kept" measured missing data plus a collapsed candidate list, not the new instruments.

## 5. Second bug (queue logic), fixed

`duc_ext_fixed_queue.sh` ran `duc_fix_retemplate.py --final --upto $DUC_END` at the end of ANY run. With
`DUC_END=2012-07` that force-finalised epochs whose 24-month template window was still open. It now runs
`--final` only when `DUC_END` is the real end of the record (2026-09).

**Still to do (not applied):** the cache stage should fail loudly when more than a few % of the planned station-days
error. The queue currently logs the counts and carries on (`duc_ext_fixed_queue.sh` cache loop).

## Re-run results (auto-written by `scripts/duc_diag.py summary`)

_2026-10-01 14:40 UTC._ Per month: mean templates (channels) available per family-day, raw 8 x MAD detections, detections passing his cc x nchannel rule, median nchannel of the detections; before = the OOM-damaged scan (moved to `mf1_pre_diag/`), after = the re-run.

| month   |   ntpl_before |   raw_before |   rule_before |   nch_det_before |   ntpl_after |   raw_after |   rule_after |   nch_det_after |
|:--------|--------------:|-------------:|--------------:|-----------------:|-------------:|------------:|-------------:|----------------:|
| 2008-01 |          8.71 |         4155 |           229 |                6 |        14.27 |        2962 |          662 |              11 |
| 2008-02 |          5.97 |         4986 |            19 |                3 |        14.6  |        2039 |          180 |               9 |
| 2008-03 |          5.52 |         4525 |            18 |                3 |        14.67 |        2243 |          458 |               9 |
| 2008-04 |         14.67 |         7441 |          4662 |               14 |        14.67 |        7441 |         4662 |              14 |
| 2008-05 |          8.18 |         3956 |            85 |                4 |        14.67 |        3468 |          371 |               9 |
| 2008-06 |          6.07 |         5038 |           125 |                4 |        14.46 |        4169 |          200 |               9 |
| 2008-07 |          6.22 |         4727 |            23 |                3 |        14.48 |        5551 |          212 |               9 |
| 2008-08 |          6.62 |         3958 |            53 |                6 |        14.67 |        5167 |          574 |               9 |
| 2008-09 |          5.25 |         5539 |            32 |                5 |        14.68 |        3674 |          435 |               9 |
| 2008-10 |          5.14 |         4137 |            12 |                3 |        14.64 |        1884 |          524 |              11 |
| 2008-11 |          4.95 |         5603 |            45 |                3 |        16.08 |        2464 |          946 |              15 |
| 2008-12 |          5.21 |         3511 |           110 |                3 |        16.02 |        2830 |         1458 |              15 |
| 2009-01 |          3.98 |         5795 |            35 |                3 |        15.86 |        1658 |          326 |              11 |
| 2009-02 |          8.84 |         2056 |            41 |                4 |        16.18 |        1273 |          150 |               9 |
| 2009-03 |          4.72 |         4753 |            38 |                3 |        16.26 |        3149 |          868 |               9 |
| 2009-04 |          5.22 |         3405 |            11 |                3 |        16.26 |        2351 |          227 |               9 |
| 2009-05 |          4.71 |         5688 |            31 |                3 |        15.86 |        4444 |          722 |               9 |
| 2009-06 |          4.6  |        11311 |            26 |                3 |        16.06 |        6346 |         1723 |               9 |
| 2009-07 |          4.99 |         6423 |            42 |                3 |        16.26 |        5687 |         1627 |               9 |
| 2009-08 |          4.59 |         5705 |            26 |                3 |        16.21 |        4549 |          291 |               9 |
| 2009-09 |          8.89 |         4117 |            72 |                9 |        16.26 |        4447 |          674 |               9 |
| 2009-10 |          4.55 |         8520 |            48 |                3 |        16.26 |        3126 |         1482 |              14 |
| 2009-11 |          5.53 |         2557 |             8 |                3 |        16.26 |        1254 |          165 |               9 |
| 2009-12 |          5.72 |         4232 |            40 |                3 |        16.26 |        2545 |          394 |               9 |
| 2010-01 |          6.04 |         5615 |            33 |                3 |        15.51 |        1955 |          325 |               9 |
| 2010-02 |          6.02 |         4447 |            17 |                4 |        15.19 |        1757 |          171 |               9 |
| 2010-03 |          9.49 |         3524 |           159 |                6 |        15.42 |        3268 |          945 |              11 |
| 2010-04 |          5.72 |         4296 |            50 |                3 |        15.44 |        2503 |          403 |               9 |
| 2010-05 |          5.41 |         3833 |            11 |                3 |        15.44 |        3076 |          187 |               9 |
| 2010-06 |          4.59 |         4694 |            25 |                3 |        15.44 |        4248 |          386 |               9 |
| 2010-07 |          6.46 |         6457 |            78 |                6 |        15.03 |        7097 |          724 |               6 |
| 2010-08 |          6.3  |         5848 |            35 |                3 |        16.05 |        5186 |          195 |               9 |
| 2010-09 |          5.02 |         7151 |            50 |                3 |        16.15 |        4902 |         1548 |               9 |
| 2010-10 |          5.36 |         4998 |            29 |                3 |        16.09 |        2769 |          916 |              11 |
| 2010-11 |          6.9  |         1806 |            38 |                5 |        16.24 |        1873 |          639 |              11 |
| 2010-12 |          5.8  |         4402 |            11 |                4 |        15.91 |       10677 |          120 |               3 |
| 2011-01 |          4.48 |         6484 |            26 |                3 |        15.44 |        9108 |          531 |               3 |
| 2011-02 |          5.04 |         6144 |            17 |                3 |        15.44 |        6516 |          318 |               3 |
| 2011-03 |         12.03 |         3060 |           102 |                6 |        15.45 |       22077 |          186 |               3 |
| 2011-04 |          5.14 |         6185 |            44 |                3 |        15.29 |        4059 |         1036 |               9 |
| 2011-05 |          6.95 |         4976 |           158 |                3 |        15.46 |       48620 |         1112 |               3 |
| 2011-06 |          7.72 |         4119 |           104 |                9 |        15.4  |       31064 |          976 |               3 |
| 2011-07 |          5.41 |         4501 |            30 |                3 |        15.63 |       16689 |          226 |               3 |
| 2011-08 |          7.17 |         3619 |            37 |                3 |        15.88 |        5407 |          726 |               9 |
| 2011-09 |          4.35 |         5151 |            18 |                3 |        16.25 |        4344 |          602 |               9 |
| 2011-10 |          5.46 |         1756 |            70 |                4 |        16.1  |        3403 |         1396 |              11 |
| 2011-11 |          5.03 |         2524 |            39 |                3 |        14.3  |        2397 |          492 |               7 |
| 2011-12 |          4.32 |         3704 |            31 |                3 |        12.73 |        3059 |          340 |               6 |
| 2012-01 |          4.28 |         4711 |            12 |                3 |        11.89 |        2910 |          185 |               6 |
| 2012-02 |          5.72 |         4448 |            39 |                3 |        12.02 |        2772 |          428 |               9 |
| 2012-03 |          5.85 |         6231 |            25 |                3 |        12.27 |        3736 |          499 |               9 |
| 2012-04 |          5.82 |         5754 |            21 |                3 |        11.94 |        4640 |          225 |               9 |
| 2012-05 |          5.76 |         3154 |            41 |                3 |        11.63 |        5246 |          359 |               9 |
| 2012-06 |          4.23 |         7751 |            14 |                3 |        11.45 |        5857 |          527 |               9 |
| 2012-07 |          4.76 |         5087 |            21 |                3 |        11.56 |        5978 |          220 |               9 |

Per year:

|    y |   raw_before |   rule_before |   raw_after |   rule_after |
|-----:|-------------:|--------------:|------------:|-------------:|
| 2008 |        57576 |          5413 |       43892 |        10682 |
| 2009 |        64562 |           418 |       40829 |         8649 |
| 2010 |        57071 |           536 |       49311 |         6559 |
| 2011 |        52223 |           676 |      156743 |         7941 |
| 2012 |        37136 |           173 |       31139 |         2443 |

**New joins with his catalogue (his rule both sides; from `duc_fix_report.py`):**

|   year |   ours_rule |   his_rule |   ratio_rule |   ours_raw |   his_raw |   ratio_raw |   ours_raw_60fam |   his_raw_60fam |   ratio_raw_60fam |
|-------:|------------:|-----------:|-------------:|-----------:|----------:|------------:|-----------------:|----------------:|------------------:|
|   2008 |       10682 |      10740 |        0.995 |      43892 |     39961 |       1.098 |            23679 |           22809 |             1.038 |
|   2009 |        8649 |       8528 |        1.014 |      40829 |     37531 |       1.088 |            16869 |           16533 |             1.02  |
|   2010 |        6559 |       6610 |        0.992 |      49311 |     52622 |       0.937 |            16349 |           16244 |             1.006 |
|   2011 |        7922 |       5717 |        1.386 |     156703 |    253093 |       0.619 |            17595 |           24309 |             0.724 |

Forced re-templating diagnostic (all new epochs finalised at 2012-07, scratch `diag_retpl/`): 174/277 family-channels kept, 23/26 epochs; median N windows 43, median candidates 44.

Production re-templating state (no --final): 2 epochs finalised.
