# Royer & Bostock (2014) Washington -- THE PLAN, written before anything runs

Agent AB, 2026-09-25.  Code written and syntax-checked; **no stage of the detection chain has been
executed** -- no download, no GPU work, no 40 Hz cache, no detections.  The one exception, disclosed in
`results/WA_BUG_AUDIT.md`, is the read-only stage-0 quantifier `ab_gapscan.py`, which was launched
before the write-only instruction arrived and survived the `kill -TERM` as an orphan; its two CSVs are
on disk and stage 0 of the queue is guarded by them, so a launch will skip it.
Workspace `results/ab_royer/`, code `scripts/ab_{common,gapscan,fetch,cache40,auto,cluster,grow,locate,null,score}.py`,
queue `scripts/ab_royer_queue.sh`.  Bug audit: `results/WA_BUG_AUDIT.md`.

Source read directly: `notes/1-s2.0-S0012821X13004780-main.pdf` (sections 2, 3, 3.1, 3.2, 3.3, 5) and
`notes/royer2014_supp.pdf` (Table S1).  Every parameter below carries its section.

---

## 1. Their parameters, with our current values beside each

`OURS (before)` = what the Washington runs on record actually used, with file:line.
`OURS (ab_)` = what the new chain does.

| # | Parameter | Royer & Bostock 2014 (verbatim, with cite) | OURS (before) | OURS (`ab_`) |
|---|---|---|---|---|
| 1 | segment | "divided into 24 hour-long segments" (sec.3) | day files, hour blocks `k_auto.py:198` | day files, hour blocks |
| 2 | band | "band-pass filtered between 1 and 8 Hz" (sec.3) | 1-8 Hz `d_equiv.py:5-8` | **1-8 Hz** (match) |
| 3 | sample rate | "resampled to **40 sps**" (sec.3) | **20 Hz** `k_common.py:28` | **40 sps** `ab_common.FS` |
| 4 | autocorr window | "**15 s** windows" (sec.3.1) | **6.0 s** `k_auto.py:42` (from Bostock 2012, a different paper and region) | **15.0 s** `ab_common.WIN_S` |
| 5 | window lag | "lagged by **0.5 s**" (sec.3.1) | 0.5 s `k_auto.py:43` | 0.5 s (match) |
| 6 | windows per hour | "a total of **7170** individual windows per hour" (sec.3.1) | n/a (6 s window -> 7189) | 7171 on the 15 s/0.5 s grid; reported |
| 7 | what is correlated | "Each window was correlated with all other windows in the same 1-hour record for an individual (station-component) channel" (sec.3.1) | rows on the 0.5 s grid x **all sample lags** `k_auto.py:79-90` | **both**: `--cols window` (the literal 7171 x 7171 window-pair matrix, headline) and `--cols sample` (the old reading, reported) |
| 8 | stacking of the channels | "The resulting time-series from all stations and all components ... were stacked" (sec.3.1) | plain mean over channels, constant divisor `k_auto.py:90` | mean over the channels **valid at that pair**, per-pair divisor |
| 9 | autocorr threshold | "exceeded **8 times the median absolute deviation**" over the hour (sec.3.1) | 8 x MAD `k_auto.py:45`, pool includes masked zeros | 8 x MAD, **valid pairs only** in the pool |
| 10 | NW stations | "the **7 three-component stations B001, SQM, BS11, W020, W040, GNW and PL11**" (sec.3.1) | the 52-station probe-admitted roster `n_common.py:171-208`; **SQM, GNW, W020 absent**, GNW explicitly rejected `n_report.py:200` | exactly those 7 |
| 11 | SW stations | "the **7 three-component stations N050, S030, W070, N060, PL11, S040 and N070**" (sec.3.1) | 6 of the 7; **N070 absent** | exactly those 7 |
| 12 | NW episodes | "**January 17-30, 2007, May 01-23, 2008, August 07-24, 2010 and August 04-23, 2011**" (sec.3.1) | 2008-05 and 2010-08 only | all four (75 days) |
| 13 | SW episodes | "**January 14-February 01, 2007, May 01-25, 2008**" (sec.3.1) | not run as a subarray | both (44 days) |
| 14 | cull | "culled to exclude overlapping detections and retain only those with high SNR" (sec.3.1) -- **numerically SILENT** | `DECL_S = 4.0` `k_auto.py:46` (that number is Bostock 2015's post-hoc magnitude cull, not a detector rule); no SNR half at all | overlap cull at one window length (15 s), on both axes; **SNR measured per detection as a COLUMN, never a cut** |
| 15 | grouping | "combine all channels ... into a single 'super' trace and cross-correlated all such traces against one another ... populate a similarity matrix employed within a **hierarchical cluster analysis**" (sec.3.2, Rowe et al. 2002) | union-find **connected components** `k_group.py:51-65,96-99`, which chain | super traces -> lagged all-pairs cc -> similarity matrix -> `scipy` hierarchical clustering |
| 16 | linkage / cut | **SILENT** | n/a | swept: 4 linkages x 34 cuts; headline chosen by the **mean silhouette** (data-internal, never looks at 224/70); the cut nearest the published count is a separate reported column |
| 17 | min cluster size | **SILENT** | `--min-members 3` `k_group.py:72` | headline 2; 3/4/5 reported as columns |
| 18 | growth stations | "subsets of **7 to 10** high SNR, three-component stations available within a template" (sec.3.3) | 51-52 stations, each an independent vote | the subarray's 7 by default (the lower bound, literal); `--pool` widens to the 25-station regional pool of sec.2.1 as a stated, costed extension |
| 19 | growth scan | "scanning through **all available tremor dates**" (sec.3.3) | episode days | episode days (75 NW / 44 SW); a wider tremor-date scan is a reported extension |
| 20 | growth threshold | "the **summed** network cross-correlation coefficient exceeded **8 times the median absolute deviation**, now for a **24-hour period**" (sec.3.3) | 8 x MAD per group over the block `k_grow.py:202-204` | 8 x MAD over the 24-hour period; `--net mean` (mask-consistent) headline, `--net sum` reported |
| 21 | restack normalisation | "stacking waveforms (**normalized to unit maximum amplitude across three components**)" (sec.3.3) | **rms** normalisation `k_common.norm3` | `ab_grow.unit_max3` -- max\|x\| over the three components jointly |
| 22 | how many detections stack | "**all** newly registered detections" (sec.3.3) | `MAX_STACK = 200` `k_grow.py:39` | no cap (`--max-stack 0`); a cap only as a reported ablation |
| 23 | iterations | "After **several** iterations" (sec.3.3) -- count SILENT | fixed **3** `k_grow.py:234` | iterate to convergence, `--tol 0.01` on the detection count, `--max-iters 8` |
| 24 | template length | **SILENT** in sec.3.3 | 12 s `k_grow.py:36` | 15 s (internally consistent with sec.3.1's own window); 6 s reported as an ablation |
| 25 | location | "Hyp2000 (Klein, 2002)" then "**hypoDD** (Waldhauser, 2001)", Table S1: IDAT 3, 5 x NITER 8 = **40 iterations**, DAMP **50,50,40,40,40**, WTCCP/WTCCS 0.01,0.01,1,1,1, WTCTP/WTCTS 1,1,0.01,0.01,0.01 (sec.5) | Crosson 1-D grid search, VP/VS 1.78, home-grown damped-LSQR DD | the **rebuilt, calibrated** locator `L9_xcdiff` + `L12_comb` (`WA_LOCATION_REBUILD.md`: 4.0 km epicentral median against Bostock's published coordinates vs 27.6 km for the old one) with per-family station-bootstrap radii.  Neither Hyp2000 nor hypoDD is installed here; their numbers are the reference lines, and the double-difference pass is `uu_dd.py`, which was measured and **not adopted** in the north |
| 26 | depths | "**31 and 46 km** depth (NW) and **32 and 49 km** (SW)" (sec.5) | n/a | reported as the depth distribution and the fraction inside the band |
| 27 | location sigma | Hyp2000 mean nominal 2.0/1.3 km (NW), 3.7/3.5 km (SW) (sec.5) | n/a | our comparable number is the bootstrap radius, reported beside theirs (NOT the same statistic -- theirs is a linearised formal error, ours a station-resampling radius) |
| 28 | hypoDD retention | "**77 %** and **72 %** of the original events remain ... for the NW and SW datasets" (sec.5) | n/a | reported if a DD pass is run; it was not adopted in the north |

### What the paper does NOT give us
No catalogue, no event list, no template waveforms (Appendix A is Table S1 + Figures S1-S4 only), no
station list for the 25-station NW / SW pools beyond Figure 1D, no linkage, no cut, no minimum
cluster size, no iteration count, no numeric cull rule, no template length.  Every one of those is a
swept column here, never a guessed constant.

---

## 2. Their 14 stations, on disk or not

Networks and coordinates from `results/fdsn_stations_40-50N_2003-2014.txt`.  PL11 serves both
subarrays, so the 14 named slots are 13 distinct stations.  XU is the Flexible Array CAFE experiment;
the SW six (N050, S030, W070, N060, S040, N070) end **2008-12-31**, which is exactly why the paper
says only 2007 and 2008 were available for SW (sec.2.1).

| station | net | lat, lon | episodes needed | station-days needed | raw on disk | 20 Hz cache | to fetch |
|---|---|---|---|---|---|---|---|
| B001 | PB | 48.043, -123.131 | NW x4 | 75 | 55 | 74 | **20** (NW11) |
| SQM | UW | 48.074, -123.048 | NW x4 | 75 | **0** | **0** | **75** (all four) |
| BS11 | XU | 47.957, -122.925 | NW x4 | 75 | 55 | 75 | **20** (NW11) |
| W020 | XU | 47.873, -122.946 | NW x4 | 75 | 55 | 75 | **20** (NW11) |
| W040 | XU | 47.601, -123.121 | NW x4 | 75 | 55 | 75 | **20** (NW11) |
| GNW | UW | 47.564, -122.825 | NW x4 | 75 | 55 | 74 | **20** (NW11) |
| PL11 | XU | 47.454, -123.154 | NW x4 + SW x2 | 119 | 99 | 119 | **20** (NW11) |
| N050 | XU | 47.224, -123.396 | SW x2 | 44 | 44 | 44 | 0 |
| S030 | XU | 47.029, -123.297 | SW x2 | 44 | 44 | 44 | 0 |
| W070 | XU | 46.925, -123.228 | SW x2 | 44 | 44 | 44 | 0 |
| N060 | XU | 47.151, -123.164 | SW x2 | 44 | 44 | 44 | 0 |
| S040 | XU | 47.040, -123.027 | SW x2 | 44 | 44 | 44 | 0 |
| N070 | XU | 47.052, -122.715 | SW x2 | 44 | 44 | 44 | 0 |
| **total** | | | | **833** | **638** | 756 | **195** |

**The data gap is far smaller than assumed.** GNW, W020 and N070 are all on disk (the earlier claim
that they were missing was about the *probe-admitted roster*, not the archive).  Only two gaps exist:

1. **UW.SQM** -- absent entirely from both the raw traces and the 20 Hz cache.  75 station-days.
   (`traces_eq/UW.SQM` holds earthquake windows only, not continuous days.)
2. **the 2011-08-04..23 NW episode** -- present in the 20 Hz cache for six of the seven stations but
   its **raw miniSEED was freed**, and a 40 sps cache can only be built from raw.  120 station-days.

**Volume.** Measured mean on-disk size of these stations' day files is **19.7 MB** (range 14-36 MB;
XU at 50 Hz native, PB.B001 at 100 Hz, UW.GNW at 40/50 Hz).  195 station-days -> **~3.9 GB**
(budget 3.5-5 GB).  The 40 sps cache for all 833 station-days is **~37 GB**
(3 components x 40 sps x 86400 x 4 bytes + a 1-byte/sample mask).  `df` shows 69 TB free.
FDSN concurrency is capped at 8 in `ab_fetch.py` and the per-station-day result (bytes, channels,
native rate, coverage) is recorded in `results/ab_royer/fetch_manifest.csv`.

**One risk to check at fetch time, not assumed away:** the paper calls SQM a *three-component*
station.  If FDSN returns only a vertical channel for UW.SQM in 2007-2011, `ab_cache40.py` records
`status=error` with the channel list and the NW subarray runs with 6 of the 7, which is then a
reported column, not a silent substitution.

---

## 3. Stage-by-stage checkpoints, with their published yields

Every "ours" cell is a number the queue produces; every "theirs" cell is a REFERENCE LINE and never
a gate (`rank-dont-discard`).  The null receives the identical table.

| stage | script | inputs | outputs | checkpoint, theirs |
|---|---|---|---|---|
| 0 | `ab_gapscan.py` | the 20 Hz cache, Royer station-days | `gapscan_cache20.csv`, `gapscan_hours.csv` | **already done** (see the note above): 58 of 18,138 station-hours trip the degenerate guard, on GNW / PL11 / S040 / N060.  The stage will be skipped by its guard |
| 1 | `ab_fetch.py` | `fetch_plan.csv` | `traces/...`, `fetch_manifest.csv` | 195 station-days fetched, ~3.9 GB |
| 2 | `ab_cache40.py` | raw miniSEED | `cache40hz/...`, `cache40_{NW,SW}.csv` | 833 station-days at 40 sps, 1-8 Hz, with a mask |
| 3 | `ab_auto.py` | the 40 sps cache | `<SUB>/auto_window/*.parquet`, `_stats.csv` | windows/hour **7170**; **NW 4,915** / **SW 3,306** initial detections |
| 4 | `ab_cluster.py` | stage 3 | `detections.parquet`, `simmat.npz`, `sweep.csv`, `tpl_it0.npz` | initial templates **NW 224** / **SW 70** |
| 5 | `ab_grow.py` | `tpl_it0.npz` | `tpl_it<K>.npz`, `det_it<K>.parquet`, `iters.csv` | final templates **NW 122** / **SW 54**; "100's to 1000's of independent detections" per template |
| 6 | `ab_locate.py` | last `det_it`/`tpl_it` | `loc_L9_xcdiff.csv`, `loc_L12_comb.csv`, `depth_summary.csv` | depths **31-46 km (NW)**, **32-49 km (SW)**; their sigma 2.0/1.3 and 3.7/3.5 km; hypoDD keeps 77 %/72 % |
| 7 | `ab_null.py` | all of the above with `AB_NULL=k` | `<SUB>/null<k>/...` | the identical table; evidence is the real-vs-null separation curve |
| 8 | `ab_score.py` | everything | `score/stage_counts.csv`, `families_ranked.csv`, `separation_curves.csv`, `SUMMARY.md` | the scorecard |

**If a count comes in far below 122/54, that is a bug to find, not a finding**
(memory `reproduce-known-answer-first`).  The stage where the loss happens is visible in
`stage_counts.csv` because every stage's count sits in one row beside theirs.

---

## 4. The Stage-0 bug-class fix, carried by every `ab_` stage

`WA_BUG_AUDIT.md` is the full argument.  In short: the southern unbounded-cc bug is **not** in the
Washington engines (they set a rejected window's denominator to `+inf`, not to a 1e-15 floor, so every
cc is bounded by 1), but three weaker forms of the class are -- a guard relative to the median that
degenerates to `ew > 0` on a mostly-dead channel (`k_auto.py:69-70`, `k_grow.py:167-169`,
`o_auto.py:66-76`), masked pairs inside the median/MAD pool, and a constant channel divisor
(`k_auto.py:90`, `k_grow.py:199-200`).  Measured exposure: 24 of 3,867 station-days in the 2008 Puget
roster (0.62 %), 3 of 6,817 in the 2010 Olympic roster (0.04 %), and -- at the hour resolution the MAD
actually uses -- **58 of 18,138 station-hours (0.32 %) on Royer's own bed, on four stations the
day-level test calls clean (PL11, S040, N060, plus GNW)**, because a single dead component inside a
live three-component station passes every all-component zero test.  Two further defects: the 20 Hz
cache keeps **exact** 50,174-second zero runs, and the manifest column written to catch them
(`frac_nonzero`) reads 1.000000 on a file that is 41.9 % non-zero.  No recorded Washington number is
invalidated, and the bug cannot explain the 2008 zeros (it inflates counts; that result was zero).
The `ab_` chain does not inherit any of it: `ab_cache40.py` writes a per-sample validity mask,
`ab_common.window_stats` admits a window only when every sample is valid, and both detectors divide
by the **per-sample** count of contributing channels.

---

## 5. Scoring

Royer released no catalogue, so this is an **existence-and-count** reproduction, exactly as
`notes/MARGIN_CATALOG_PLAN_2026-09-24.md` prescribes for a region with no released event list:
do we produce located, repeating families above the whole-chain null, inside the published region and
depth band, and how do our per-stage counts stand beside theirs.  Rules: no gates anywhere; every
family keeps its row in `families_ranked.csv`; every criterion is a column (`n_det`, `n_sta`,
`nmad_med`, `cc_med`, `n_days`, `n_cluster_members`, `boot_km`, `dep`, `in_published_depth_band`,
`rail`); the null gets the identical table; the evidence is `separation_curves.csv` -- the real and
null counts as a function of each criterion -- plus a ranking.  Their 4,915/3,306, 224/70, 122/54,
31-46 km and 32-49 km are drawn as reference lines.

---

## 6. Cost -- wall time, GPU, disk

H100 NVL (95 GB, ~2 GB in use), 32-CPU cgroup cap, 69 TB free.  Two southern queues (pids 1509294,
1800120) share the GPU; peak here is ~8 GB in `--cols window` and ~20 GB in `--cols sample` and
stage 5, so they coexist.

| stage | wall time (estimate) | GPU | disk |
|---|---|---|---|
| 0 gapscan | done (8 min, 1 core) | - | 0.8 MB |
| 1 fetch | 30-90 min (FDSN latency-bound, 8 connections) | - | ~3.9 GB |
| 2 cache40 | 15-45 min, 24 workers | - | ~37 GB |
| 3 auto `window` (both subs, 2,856 hour-blocks) | **4-12 h** | ~8 GB, ~1.3 TFLOP/block | ~2-5 GB of parquet |
| 3b auto `sample` (the reported alternative reading) | **16-32 h** | ~20 GB, ~26 TFLOP/block | ~5-15 GB |
| 4 cluster (per subarray) | 0.5-1.5 h | ~10 GB for the 161-lag GEMM sweep | ~1 GB |
| 5 grow NW (224 templates x 75 days x up to 8 iterations) | **20-40 h** | ~20 GB | ~5-20 GB |
| 5 grow SW (70 x 44) | **5-10 h** | ~10 GB | ~2-6 GB |
| 6 locate (per subarray, boot 60) | 2-4 h | small | < 1 GB |
| **to the end of stage 6, both subarrays** | **~3-5 days** | | ~60-90 GB |
| 7 null, per realisation per subarray | ~1-2 days | as above | as above |
| 7 null, 3 realisations x 2 subarrays | **4-8 days** (the long tail) | | ~200-400 GB |
| 8 score | minutes | - | small |

The queue is ordered so the headline is reached before the optional work: the `window`
autocorrelation finishes and writes `AB_ROYER_AUTOCORR_DONE` **before** the `sample` variant runs,
and stage 8 is re-run after **each** null realisation so a partial null is still a usable null.

---

## 7. Exactly what the queue does when someone launches it

`scripts/ab_royer_queue.sh`, in order, each stage guarded by its own output file so a relaunch
resumes rather than redoing hours of work:

1. takes `locks/ab_royer.lock` with `flock -n` and exits quietly if another copy holds it; writes its
   pid to `results/AB_ROYER_PID`;
2. installs `TERM INT HUP EXIT` traps that write `results/AB_ROYER_FAILED` with the stage name and
   return code, and starts a background heartbeat that touches
   `results/AB_ROYER_HEARTBEAT` every 120 s;
3. **stage 0** `ab_gapscan.py` -> `results/AB_ROYER_AUDIT_DONE`;
4. **stage 1** `ab_fetch.py --workers 8` (the 195 station-days);
5. **stage 2** `ab_cache40.py --sub NW` then `--sub SW`, 24 workers;
6. **stage 3** `ab_auto.py --sub NW --cols window`, then SW -> `results/AB_ROYER_AUTOCORR_DONE`;
   then the `--cols sample` pair;
7. **stage 4** `ab_cluster.py` for each subarray;
8. **stage 5** `ab_grow.py --max-iters 8 --tol 0.01` for each subarray;
9. **stage 6** `ab_locate.py --boot 60` on the last growth iteration of each subarray;
10. **stage 7** `ab_null.py --only K` for K = 1, 2, 3 and each subarray, with **stage 8**
    `ab_score.py --sub NW,SW` re-run after every K;
11. writes `results/AB_ROYER_DONE` and exits 0.

Every stage runs as `env AB_STAGE=<n> python ...`.  **To stop it: `kill -TERM $(cat
results/AB_ROYER_PID)` -- never `pkill -f`.** It does not touch `lfe_south/` and does not read or
write any file of the two running southern queues.

### Decisions that are deliberately left open for review before launch
* **the growth station pool (row 18).** The default is the subarray's own 7, the literal lower bound
  of "subsets of 7 to 10".  Widening to the 25-station regional pool of sec.2.1 is the single change
  most likely to raise the final template count toward 122, and it costs roughly 40 extra stations x
  119 days ~ 4,760 station-days, ~95 GB of raw and ~200 GB of cache.  It is a `--pool` argument, so
  it needs no code change, but it should be an explicit decision, not a default.
* **`--cols window` versus `--cols sample`** (row 7): both run, `window` is the headline.
* **the number of null realisations**: three are queued one at a time; scoring after the first is
  already meaningful.

---

## 8. The launch command

```
cd /home/jovyan/tremorferometry/lfe_array_discovery && nohup setsid -f bash scripts/ab_royer_queue.sh >> logs/ab_royer_queue.log 2>&1 < /dev/null
```
