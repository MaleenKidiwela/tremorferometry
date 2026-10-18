# Royer & Bostock (2014) reproduction -- running log

A running log: tests, failures, fixes, stops and relaunches, and open items. It is appended after
every step. The method and audit are in `METHOD_FROM_PAPER.md`; the pools and scan dates are in
`POOL_RECONSTRUCTION.md`.

## v1 (2026-09-25 09:12Z, pid 1860580)
- **Smoke before launch** (NW, 2010-08-15/16, 7 stations; `smoke/`):
  - 213 autocorrelation pairs in 48 h;
  - 338 detections -> 136 templates;
  - growth, one iteration: 1237 -> 2314 detections;
  - L9 locator: 56/120 located, depth median 24 km.
- **Fixes before launch** (`METHOD_FROM_PAPER.md` table A1-A13, S1-S4). Three of them blocked the
  run:
  - the cut sweep stopped at distance 0.75, so every setting gave 0 templates;
  - the shift sign in the initial template stack was inverted;
  - `ab_locate.py` crashed on a missing `station_corrections`.
- **Full-run counts:**
  - autocorrelation window mode: NW 5025 (paper 4915), **SW 7981 (paper 3306)**;
  - sample mode: NW 39792, SW 50270;
  - initial templates: **2640 / 2192** (paper 224 / 70).
- **Stop:** 11:56Z, by the coordinator's instruction (the pool, scan dates and locator were not yet
  literal). Exact pids killed: 1882486 (ab_grow), 1860580 (queue), 1860588 (heartbeat).
  Outputs were renamed `*_v1_*`, not deleted.

## v2 (12:32Z, pid 1904278)
Changes from v1:
- 25-station pools from Fig.1D;
- scan days from the PNSN catalogue and ETS windows;
- Hyp2000 + hypoDD built in `/home/jovyan/tools` and wired as stage 6;
- strict time-overlap cull;
- cluster headline = max real-minus-null-1 count;
- the extension wired after DONE.

Tests:
- Hyp2000 reproduces its USGS test event exactly.
- hypoDD reproduces its example1 to <= 62 m.
- Synthetic recovery: Hyp2000 0.18 km, hypoDD 0.46 km.
- Wide-pool growth smoke: 3967 -> 13734 detections.
- Hyp2000/hypoDD on smoke templates: 28/53 located, depth median 34.3 km, hypoDD kept 3.

Full-run counts:
- strict cull: NW 3258 (-34 %), SW 2936;
- headline clustering: NW 2861, SW 2466.

**Bug found (15:04Z):** `ab_cluster.py` merged the null sweep on exact cut values. The real
run's median + 8 MAD cut (0.8942) is not in the null's grid, so `n_null` was NaN, filled with 0,
and that row "won" with excess = n_real. **The headline was wrong.**

**Stop:** 15:04Z, pids 1909262 (ab_grow NW, running 2.2 h), 1904278 (queue), 1904286
(heartbeat). `NW/SW clusters` -> `clusters_v2_nullmerge_bug`; `NW/grow` -> `grow_v2_stopped`.

## Round 3 (coordinator + user pushback), 15:05Z onward

### Item 4: NW autocorrelation -34 %, diagnosed
Hour and channel accounting (`NW/auto_window/_stats.csv`) found no processing loss:

| subarray / episode | hours | stations present | window pairs (v2) | strict pairs |
|---|---|---|---|---|
| NW 2007-01 | 336/336 (317 with 21 channels) | all 7 in 1747/1800 hours | 1212 | 742 |
| NW 2008-05 | 552/552 | | 1401 | 878 |
| NW 2010-08 | 432/432 | | 1083 | 754 |
| NW 2011-08 | 480/480 (SQM vertical-only) | | 1329 | 884 |
| SW 2007-01/02 | 456/456 | all 7 in 1053/1056 | 5312 | 1449 |
| SW 2008-05 | 600/600 | | 2669 | 1487 |

The per-hour rate was lower in NW under every cull reading, so the shortfall is not a data-loss bug.
It is the **reading of the MAD**. The paper says "the resulting time-series ... were stacked ...
When the network autocorrelation exceeded 8 times the MAD". Each window's correlation with all
other windows is its own time series (one row of the hour matrix). Taking the MAD **per row**
(`ab_auto.py --mad-pool row`), with the pair-overlap cull:

| | NW | SW |
|---|---|---|
| paper | 4915 | 3306 |
| hour MAD + pair suppression (v1) | 5025 | 7981 |
| hour MAD + strict cull (v2) | 3258 | 2936 |
| **row MAD + pair suppression (v3)** | **4603 (-6 %)** | **3217 (-3 %)** |
| row MAD + strict cull | 3243 | 2106 |

**Adopted: row MAD + pair suppression.** It is the most literal reading, and it is the only one that
matches both subarrays. Per episode (row MAD): NW 815 / 1406 / 1057 / 1325; SW 1610 + 109 / 1498.

### Item 2: velocity model, resolved
- Royer's PhD thesis (UBC 2014; full text via the UBC OC API, collections/24/items/1.0165914)
  repeats the paper and adds only "the 1-D, fixed VP/VS model used for LFE location".
- Bostock et al. (2015) use "the 1-D model used for routine earthquake location by the Geological
  Survey of Canada".
- Savard (2018, UBC PhD, Bostock group) locates LFEs with HYPOINVERSE and "the 1D GSC velocity
  model used for routine locations". Its Fig.4.7a draws the "GSC onshore model".
- I digitised that curve from the PDF's vector paths (exact coordinates):

  | depth (km) | Vp (km/s) |
  |---|---|
  | 0-1 | 5.01 |
  | 1-6 | 6.01 |
  | 6-30 | 6.71 |
  | 30-46 | 7.11 |
  | > 46 | 7.76 |

- GSC uses a uniform Vp/Vs of 1.73.
- **Headline = GSC onshore model, Vp/Vs 1.73.** Crosson (1976) and Vp/Vs 1.78/1.85 are columns
  (`ab_hypo.py MODELS`).

### Item 3: Rowe et al. 2002 clustering, resolved
- Rowe et al. (2002, PAGEOPH; technique in BSSA 92:1660) is paywalled. The best accessible detailed
  description is USGS PP1769 ch.5 (Augustine), which used Rowe's code. It describes "a
  dendrogram-based, hierarchical **pair-group** algorithm ... CC coefficients are used to fuse event
  pairs with high waveform similarity into clusters, and clustering is considered complete when a
  **user-defined CC coefficient cutoff threshold** is reached" (0.80 there). Minimum family size
  is a later, separate user choice (>= 10 at Augustine; Soultz families used cc > 0.85-0.9 and
  >= 6 stations).
- So: linkage = pair-group with arithmetic averaging (UPGMA = scipy `average`). The cut is a user
  cc threshold, which Royer does not state. No minimum size is implied.
- **Headline (A16):** average linkage, cut at the paper's own significance rule (similarity >
  median + 8 MAD), min size 2. The full linkage x cut x min-size table, null-1 counts, excess and
  null/real ratio are reported. The max-excess choice (A15) is dropped: it always favours pairs.
- Diagnostic: at the 0.88-0.92 cuts the site-shift null yields 30-100 % as many clusters as the
  real data. The 8 x MAD autocorrelation bar admits a fixed rate of noise pairs by construction.
  Separation from the null therefore has to come at growth, not at clustering.

### Item 1: pre-2009-08 tremor days, in progress
- No public machine-readable Wech list before 2009-08 was found:
  - the PNSN API starts on 2009-08-06;
  - `assets.pnsn.org/oldets` has only images;
  - the PNSN-hosted Wech PDFs did not download (error pages).
- Building our own WECC detector (`scripts/ab_wecc.py`, after Wech & Creager 2008): 1.5-5.5 Hz
  envelopes, 5-min windows, pairwise envelope cc, S-wave grid search at 40 km, in-region test,
  day = >= 12 tremor windows (1 h).
  - Calibration and validation against PNSN day labels, 2009-08-06..2011-12-31 (NW).
  - Then applied to 2007-01..2009-08-05 (NW) and 2007-2008 (SW).
- Raw fetch: 11,661 station-days, started 15:15Z (`logs/ab_wecc_fetch.log`).

### Item 1: WECC detector, built and tested (15:24Z)
- Smoke test on 2010-08-15 (a strong tremor day), NW 7 stations: 288 windows, coherence median
  0.43, max 0.89. In-region windows above 0.5/0.6/0.7: 82/39/12. The code path is verified.
- After CAFE XU paused (2010-10 to 2011-06), fewer than 4 of the 7 stations recorded. Those days
  are dropped by the coverage filter (>= 50 % of windows with >= 4 stations), so validation uses
  2009-08..2010-10 and 2011-06..12.
- The raw fetch (11,661 station-days) was run standalone first. I stopped it by exact pid
  (1918804, 1918802) at 2600/11661 so the queue owns it (it skips existing files and `.nodata`
  markers).

### Item 5: NVI added (15:21Z)
- Fetch: 351 station-days; data for 257, nodata for 94.
  - VW01-03, VI50 and TASB have data only in 2006.
  - VI03/07/09/54 have only 6 days each.
  - The CN EHZ stations have 1-3 days each.
- Cache: 250 built, 7 vertical-only.
- **Smoke** (2006-09-07 + 2007-06-15; `smoke/NVI/`):
  - 7-of-10 selection works (2006: VW02/VW03 used; 2007: VI08/VI06 used);
  - 361 autocorrelation pairs (about 180/day); 559 detection times;
  - 138 clusters at UPGMA + 8 MAD (nearest-published row: complete 0.96, min size 5 -> 54);
  - growth, one iteration: 4236 -> 9193 detections, 19-20 stations;
  - **Hyp2000 (GSC model, Vp/Vs 1.73, R 1.5): 44/69 templates located, depth median 36.9 km,
    91 % inside the paper's 33-39 km; ERH 7.1 / ERZ 7.2 km; hypoDD kept 55 %; bootstrap
    radius 2.5 km.** The Crosson model gives 36.4 km and 84 % in band.

### Queue v3
Changes: stage W (WECC fetch/env/detect/calibrate) -> scan days v3 -> fetch/cache v3 (NW, SW,
NVI) -> autocorrelation with row MAD -> null-1 sweep (reported) -> UPGMA + 8 MAD clusters ->
growth -> Hyp2000/hypoDD (GSC) + secondary L9 -> 3 nulls -> `AB_ROYER_DONE` (covers all three
subarrays) -> extension (NW, SW, NVI; `AB_ROYER_EXT_DONE/FAILED`).
Old v2 outputs are kept, renamed (`*_hourmad_v2`, `clusters_v2_nullmerge_bug`, `grow_v2_stopped`).

### Relaunch v3: 2026-09-25 15:25Z, queue pid 1928163
The first stage started (stageW1_wecc_fetch, pid 1928180).

Open items, to check in the full run:
- WECC calibration kappa and TPR/FPR (`tremor/wecc_calibration.csv`), and the resulting
  pre-2009 scan-day counts;
- NW/SW/NVI initial-template counts vs 224/70/54 (UPGMA + 8 MAD will be in the thousands, not
  tuned);
- growth runtime with thousands of templates (NW ~2500 x ~180+ days, several hours per iteration);
- hypoDD retention vs 77/72/91 %;
- whether the XU Big-Skidder BS elements dominate the extension ranking;
- the auto_sample variant is still hour-MAD (a reported column only).

### 2026-09-25 ~20:45Z  Stage-4 super-trace similarity audit (NW): NO BUG, queue left running
Suspicion: NW similarity median 0.048 / p99 0.083 and 2695 UPGMA+8MAD clusters (Royer 224)
looked inconsistent with detections registered at network autocorr cc ~0.13-0.15.
Tests on `NW/clusters/simmat.npz` + `NW/auto_window/*.parquet` (scratch scripts, no outputs written):
- 4603 autocorr pairs; 3168 have BOTH window times in the 7255-time catalogue (the other 1435 lost
  a partner to the per-time overlap cull in `detections()`).
- Known pairs: S median **0.136** (p10 0.119, p90 0.174) vs their autocorr cc median **0.134**;
  S - cc median 0.0000, p95 +0.019 (lag search only adds), S >= cc for 99.9 %; best lag 0, 21 channels.
- Random pairs: S median 0.048, p90 0.064.
- Direct recompute from the 40 sps cache (200 random pairs, autocorr definition: lag 0, mean over
  valid channels of per-channel normalised cc): identical to the autocorr cc (max |diff| 0.0000)
  and to the matrix entry. Window offset, 760-sample padding (15 s + 2x2 s lag), channel order,
  normalisation, filter: all consistent.
Why 2695 clusters: the pairs are high but NOT transitive. For a known pair (i,j), third detections k
with S(i,k) > 0.12 have median S(j,k) = **0.090** < the 8-MAD cut (similarity 0.1035). UPGMA
therefore fuses the pairs but rarely grows them: headline sizes = 2354 pairs, 145 triples, ...,
max 50; 593 singletons. Each detection has median 8 neighbours above the cut (p90 47), but they do
not form mutually similar groups. Average-linkage sweep (min size 2): cut 0.80 -> 171, 0.84 -> 593,
0.88 -> 2495, 0.8965 (headline) -> 2695, 0.92 -> 2177, 0.94 -> 745. Nearest-224 rows (reported
only): single 0.90 -> 228 (max size 6724, one giant chain), complete 0.88 min3 -> 218,
average 0.88 min3 -> 231. Interpretation: the 15 s, 21-channel super-trace cc of a real pair
(~0.13) is only ~12 MAD above background, so family structure is weak at this window length;
this is a property of the data/recipe, not a code defect. Queue pid 1968883 NOT stopped.

## 2026-09-25 -- stage 5 matched_filter GPU speed-up (v4 queue stopped, relaunched as v5)
`scripts/ab_grow.py::matched_filter` rewritten for speed, method unchanged (old code kept as
`matched_filter_ref`; pre-edit file `scripts/ab_grow.py.bak_v4`): overlap-save block correlation
(L = 8192, S = L - 599), every template channel's FFT computed once per call and resident on the GPU,
station-day spectra once per date, complex64 products / float32 inverse FFTs, 1/nrm folded into the
template spectra and 1/DEN (0 where invalid; DEN still inf-masked, never clamped) applied per channel,
per-offset vectorised shift-and-add (index_add), threaded CPU `window_stats` (unchanged arithmetic),
and an exact log-doubling sliding max replacing max_pool1d(k=1201) (verified torch.equal).
Equivalence test (real NW cluster stacks tpl_it0.npz, 200 random templates, days 2007-01-13 and
2009-05-11; network-statistic arrays dumped from both): max |d fn| = 4.0e-5 (max |fn| 0.32), NaN
(min_sta_valid) masks identical, detections 3456 vs 3456, all 3456 identical (template, sample).
Timing: 200 templates x 2 days 36.0 s -> 9.6 s (3.8x, dominated by fixed per-day cost); full 2695
templates one day: ~139 s/day (running v4 queue) -> 13.5 s (~10x), both measured with the GPU shared.

## 2026-09-25 late -- stage 5 fixes: GAP GUARD + duplicate-family MERGE (queue v4b stopped, relaunched v5)
Backups: scripts/ab_grow.py.bak_v5_pre_guard_merge, ab_common.py.bak_v5_pre_guard, ab_cluster.py.bak_v5_pre_guard.
**1. Gap guard.** v4b it0 had 6,658 detections with network cc > 1 (max 106) on 7 days / 482 templates.
Cause (2008-07-05): (a) window_stats' day-long cumsum loses precision after a huge glitch (FACB N |x| 7e6 ->
window energies read 0 or 3e-7); (b) dead channels carrying filter residue (SMW all comps ~1e-20, STW/BLN/HDW
horizontals) passed `e > 0`. New ab_common.window_stats_guarded / channel_refs / day_refs: exact block-local
energies; channel-window valid iff mask-valid, <=1 % zero/constant samples, e > 1e-3 x the channel's median over
its TRULY non-zero windows (not flat and > 1e-4 x the station's loudest component), channel dead if < 10 min of
such windows, station dead if < 1e-6 x the day's median station level. Used in matched_filter (the _ref keeps the
old guard); restack and ab_cluster zero guard-dead channels; ab_cluster super-traces use the segment guard.
Verify (tpl_it0 templates, 10 days, 19.5 s/day vs ~9-13 before): the 7 days max cc 0.806 (0 above 1), counts =
old minus the cc>1 rows within +-25 (e.g. 2008-07-05 26,246 -> 25,302; old cc<=1: 25,283); other days: 2007-01-13
667/667 and 2010-08-12 24,552/24,552 identical (grp,t); 2008-05-22 12,312 -> 12,395 (12,073 identical).
**2. Duplicate merge** (Thomas & Bostock 2015 rule, SOURCES.md S5: >10 % of the smaller family's detections within
5 s of the other's). After every iteration's detection, before restack: pair shared fractions (numba, per date);
GREEDY STAR merge -- families in descending detection count join the representative they share most with if
> 10 %, else become representatives (connected components were tried first: on v4b it0 single-linkage chaining
put 1,763/2,695 templates into ONE family, lags to 26 s -> rejected). Members' times shifted by the modal pair lag
(0.1 s bins) onto the representative, union deduplicated to one per 5 s (highest cc). Outputs: detraw_it<K>
(raw), det_it<K> (merged; used by restack/locate), families_it<K>.csv; iters.csv adds n_det_raw, n_dup_edges,
n_merge_groups, n_templates_merged_away, n_families, tpl_snr_median. Convergence unchanged (--tol on merged
detection count, --max-iters 8).
Merge checks: v4b it0 (3.56 M, incl. the bad days) 27 s: 2,695 -> 1,181 families, 857,949 detections; 10-day
guarded test: 2,454 -> 997 families, 189,713 -> 45,113 detections (2.3 s).
Stage 4 (ab_cluster.py) had the same flaw (super-trace channel valid iff sum x^2 > 0). Sampled 6 days/subarray:
new guard rejects 0.09 % (NW), 0.86 % (SW), 0.30 % (NVI) of old-valid channel-segments. Fixed identically; SW and
NVI clusters renamed clusters_v3_weak_e0_guard (rebuilt by the queue). NW clusters KEPT (0.09 %; growth restarts
from them as intended). null1 sweeps not rebuilt.
Relaunch: queue v5 (see results/AB_ROYER_PID).

## 2026-09-26 -- second location column: depth FIXED to the plate interface (user-approved)
Free-depth Hyp2000 headline stays the main result; this is a column beside it. `scripts/ab_hypo_slab.py`.
- Slab models in `data/slab_models/` (README with URLs): Slab2 cas_slab2_dep_02.24.18.grd (ScienceBase
  doi:10.5066/F7PV6JNV; headline of this column) and McCrory et al. 2012 = USGS DS633 v1.2
  Modeled_slab_XYZ_Grid (5', the smoothed ETS-suitable grid; Drafted grid kept too).
- Hyp2000 fix-depth test (USGS testone event, copy in scratch): terminator col 35 '-' + trial depth 6.00
  -> summary depth 6.00 exactly, remark '-', rms 0.06 -> 0.23; same solution from a trial epicentre
  3 km away; unmodified testone still reproduces testone.sum.
- Method: fixed-point iteration (not a grid search). Start = free-depth epicentre; depth = slab depth under
  the current epicentre, fixed; Hyp2000 re-solves epicentre + origin; undamped for 10 its, then relaxed
  x += 0.5 (T(x) - x); converged when the Hyp2000 epicentre moves < 0.1 km from its trial (max 60).
  Same picks.csv, GSC 1-D, Vp/Vs 1.73, R 1.5, MIN 4, station pool, control file (ab_hypo.write_hyp_inputs).
- NW (240 families; 48 with >= 4 phases, 9 of them exactly 4 = zero dof in free depth; 192 <= 3 phases,
  not locatable): slab2 39 converged / 9 not; McCrory 41 / 7. The non-converged are the ill-posed
  families (erh >= 17 km, most rms 1.3-8.6 s); they are kept with status=not_converged.
  Scored vs Royer Hyp2000 truth (ab_truth_score.score, 1000 rotated+shifted nulls):
  | set | n | med erh | frac 3/5/10 km | null mean 3/5/10 | p 3/5/10 |
  | free-depth headline | 48 | 5.7 | .123/.180/.254 | .072/.125/.202 | .03/.08/.17 |
  | slab2 fixed | 48 | 7.6 | .115/.189/.270 | .070/.121/.198 | .07/.05/.06 |
  | McCrory fixed | 48 | 6.3 | .123/.189/.254 | .071/.122/.199 | .05/.05/.14 |
  Median epicentral shift vs free depth 7.8 km (slab2) / 5.5 km (McCrory); median fixed depth 42.1 / 39.7 km
  vs free median 34.2 km. The 9 four-phase families: fixed-depth erh median 36 km (free 16.5), i.e. still
  unconstrained -- no extra usable families from the lower pick requirement at R 1.5.
  Verdict: fixing depth to the interface does NOT change the match to Royer's catalogue materially; both
  columns sit ~5-7 points above the rotated-positions expectation at every radius.
  Outputs: results/ab_royer/NW/hypo_slab/{loc_slab_slab2,loc_slab_mccrory,score_*,summary}.csv, summary.txt.
- SW: detached follow-up `scripts/ab_hypo_slab_sw_followup.sh` (setsid, lock locks/ab_hypo_slab_SW.lock,
  heartbeat killed in EXIT trap) waits for SW/hypo/loc_headline.csv then runs the same script;
  markers results/ab_royer/SW/hypo_slab/{PID,HEARTBEAT,DONE|FAILED}; log logs/ab_hypo_slab_SW.log.
  Launched 2026-09-26 08:11Z, pid 2245051.

## 2026-09-26 -- NW finish queue: converge growth + all-station final templates (queue nw_finish_queue.sh)
Main queue stays PAUSED. New detached queue `scripts/nw_finish_queue.sh` (lock locks/nw_finish.lock,
pid results/NW_FINISH_PID, markers results/NW_FINISH_{DONE,FAILED,HEARTBEAT}; log logs/nw_finish_queue.log).
**1. Growth continuation** `scripts/ab_grow_cont.py`: from det_it8, restack exactly as ab_grow would for
it9 (same restack arithmetic, run over 20 date-chunk processes; smoke: max |diff| vs ab_grow.restack
1.7e-16, identical keys), then cut_templates / matched_filter (GPU) / merge_duplicates unchanged.
Stop rule (Bostock 2012): |dN|/N < 0.01 AND SNR gain < 1 %, max 6 more (it9..it14). Outputs continue the
numbering in NW/grow (tpl_it9.. det_it9..), iters_cont.csv (it0-8 copied + new rows with rel_dN,
rel_dSNR, stop_rule_met), summary_cont.json, FINAL_ITER. iters.csv/summary.json untouched.
**2. All-station templates** `scripts/ab_allsta.py` -> NW/hypo_allsta/: literal NW 25-station pool only
(coordinator: no extra public stations; the ext stations within 75 km are anyway cached on only 2 of 188
scan days). Moveout: detection time t is the network alignment; every station is cut from the same absolute
[t-12, t+32] s, so each station's predicted moveout is already embedded (the "stack at predicted moveout"
step = the common-time stack; 44 s window holds S at 72 km). NOTE: ab_hypo.py already stacked all 25 pool
stations this way (picks.csv has 25 stations x 240 families) -- so the new ingredients are the converged
detections and the refinement test below.
**Per-detection x-corr refinement REJECTED (smoke, 5 largest families, 40 days):** leave-one-out station
stack, +-4 s window at its envelope peak, lag +-0.5 s, 2 passes. Non-growth-station SNR 2.24 -> 7.72, BUT the
same procedure on a NULL (every detection time +317 s) gives 2.23 -> 7.99; lags are uniform on +-0.5 s
(|lag| median 0.25 s, std 0.29 s = uniform's 0.289). It aligns noise onto the reference and fabricates an
arrival -> not used (code kept, --passes). Headline templates = common-time stacks (--passes 0); a +317 s
null stack set (hypo_allsta/null317/allsta_stack_stats.csv) gives a real-vs-null SNR column per
(family, station). Also noteworthy: at the median, non-growth-station common-time SNR (2.24) equals the
null's (2.23) for those 5 families.
**3.** ab_hypo.py unchanged on those stacks (GSC/CRO x 1.73/1.78/1.85 x R1.2/1.5/2.0, headline GSC 1.73 R1.5,
boot 60). **4.** `scripts/ab_allsta_report.py` -> hypo_allsta/allsta_comparison.csv (n with picks, located,
depth median, frac 31-46, erh/erz, hypoDD kept, ab_truth_score Hyp2000+hypoDD truth) vs the it8 hypo
(reproduces 48 located / 15 of 122 within 3 km), figure hypo_allsta/NW_catalogue_allsta.png.
Launched 2026-09-26T08:20Z; expected ~45 min/iteration.

### 2026-09-26 -- `scripts/ab_grow_fast.py`: faster growth, bit-identical detector output
ab_grow.py left untouched (NW QC run in flight uses it). Same CLI as ab_grow.py, plus --prep-workers (6,
<=8), --prep-threads (3), --restack-workers (8, <=8).
**Profile** (NW 2008-05-08, 2,695 tpl_it0 templates, 15 stations): ab_grow ~14 s/day on a free GPU
(20-25 s in the production log under contention); serial CPU prep ~5 s (day load 0.8, day_refs 1.6,
15x window_stats_guarded 2.3) + GPU ~9.5 s, which also stalled on host->device index copies.
**Changes (all exact):** (1) day prep (load, day_refs, window_stats_guarded, float32 DEN =
where(okc, sqrt(max(e,0)), inf) -- ab_grow's own GPU expression, done in float64 then cast) in a spawned
process pool, days prefetched ahead of the GPU, started before the template FFTs; okc = isfinite(DEN);
(2) irfft's 1/L folded into the template spectrum (L = 8192 = 2^13: scaling by 2^-13 is exact and
commutes with every multiply/add) + irfft norm="forward"; (3) per-(chunk, batch, offset) index tensors
built once on the GPU; no gather copy when all batch rows share one offset (always, cut="common");
(4) valid-channel / valid-station counts (small integers) by one exact float32 GEMM per chunk
(template live-mask x station OK; station indicator per live-component pattern) instead of per-batch
matmul + 2 index_add_ (ref path kept when any offset != 0); (5) restack = ab_grow_cont.restack_parallel.
**Equivalence** (`scripts/ab_grow_fast_equiv.py`, 300 templates = every 9th, 2008-05-08/12/16):
6,930 vs 6,930 detections, (grp, sample) sets IDENTICAL; network function max |diff| = 0.0 over all
samples, NaN pattern identical; cc / nmad / thr diff 0.0, n_sta / n_ch equal. All 2,695 templates x 12
days (2008-05-01..12): 310,330 = 310,330 identical. Restack (`--restack`, 4 dates, 18,098 det):
identical keys, max |diff| 6e-16 (float summation order), identical template-station selection.
**Speed** (all templates, free GPU): matched filter 13.9 -> ~6.6 s/day steady state (2.1x; ~3x vs the
20 s/day seen in production); one-off ~13 s pool start per iteration; restack 22 -> 8 s on 4 dates
(scales to ~8 workers on 188 dates). Remaining GPU time is the correlation arithmetic itself
(complex product, c2r + its input clone, weighting, acc index_add_, nanmedian, sliding max).

## 2026-09-26 -- BURST MUTING + line notching, then NW rerun from autocorrelation (queue nw_mute_queue.sh)
Problem: single-station bursts/ringing seed noise families (run 1: 8 noise catch-alls absorbed 2,137/2,695 starting
templates); dropping stations moved the noise. Analyst-style cleaning instead, as a MASK (never zero-fill without masking).
**Code.** `scripts/ab_mute.py` (energy / lines / mask modes) -> results/ab_royer/NW_mute/mute/{energy,lines,mask}/;
`ab_common.load_day` applies it when env `AB_MUTE_DIR` is set (`apply_mute`): muted samples -> m False and x 0,
exactly like a gap; every stage (window validity, gap guard, refs, cut, restack) already honours m. Without the env var
nothing changes. The cache mask has no channel axis, so a burst on any component mutes the station for that span.
**Rule** (2 s windows on a midnight grid, NW pool 25 stations, 188 scan days):
R = max over components of window energy / that channel's daily median; sites = stations < 3 km apart
(B007/B04A/LRIV, B013/B943/W030, DOSE/C04A, W060/SMW); Q = median over the 5 nearest OTHER sites of the
+-6 s max of R normalised to its daily median; QL = same with a 60 s running mean.
MUTE when R > K=15 AND Q < 2 AND QL < 1.3; pad +-3 s.
Why QL: with Q alone (K15 Q<2) validation on 2008-05-08 13h muted 17% of an hour of TREMOR at W020 (W020 x10, neighbours
BS11/C04A/B001/W040/HDW/PL11 x2-3 at the same time) incl. 17/32 real-LFE windows; tremor raises the neighbours for
minutes, which QL sees. After QL: 6% of that hour, 10/32.
**K from data** (scripts/ab_mute_qc.py -> mute/qc_sweep.csv, mute_ratio_distributions.png): per (detection, template
station) 15 s windows of run 1 (NW/grow det_it14; 33 eye-real vs 136 noise families, hypo_28/review/eye_class.csv):
| K | Q | QL | muted all | real-LFE windows muted | real det >=1 sta | noise windows muted | noise det >=1 sta |
| 10 | 2 | inf | 5.9% | 5.1% | 30% | 13.8% | 42% |
| 15 | 2 | inf | 4.2% | 2.6% | 16% | 11.4% | 36% |
| **15** | **2** | **1.3** | **2.9%** | **1.4%** | **8.3%** | **8.3%** | **28%** |
| 20 | 2 | 1.3 | 2.3% | 0.9% | 5.4% | 7.3% | 25% |
| 50 | 2 | 1.3 | 1.2% | 0.3% | 1.7% | 4.8% | 17% |
Noise windows are hit ~6x more often than real ones at the chosen point.
Muted fraction by station (K15/2/1.3): B013 11.3, FACB 10.7, B943 10.5, W030 10.4, N080 5.4, W060 3.8, SMW 2.8, BLN 2.4,
STW 2.1, N050 2.0, B04A 2.0, LRIV 1.8, OPC 1.7, PL11 1.6, B001 1.3, GNW 1.3, C04A 1.2, DOSE 0.9, B007 0.9, SQM 0.7,
N060 0.6, BS11 0.5, W040 0.5, HDW 0.3, **W020 0.3 %** (mute/muted_fraction_by_station.csv).
**Key caveat:** W020's noise-family detection windows are NOT energetic: median R 2.0 (real-LFE windows at W020: 7.0),
90 % < 5. The "W020 burst" packets in those stacks emerge from many low-amplitude detections, so energy muting removes
little of them; it mainly removes the B013/W030/B943 ringing (whole noisy hours, R~100) and W060/N080/FACB bursts.
**Validation** (scripts/ab_mute_validate.py -> mute/validate/hour_*.png, zoom_*.png, lfe_window_loss.csv; days
2008-05-08, 2008-05-15 (tremor), 2011-08-03 (most W020-noise detections)): real-LFE window loss 0.9 / 0.5 / 2.1 %
(7 / 3 / 8 % of detections lose >= 1 station). Zooms: W020 impulsive bursts, W060 bursts and B013 ringing packets
muted; the real-LFE windows are untouched.
**Spectral lines** (ab_mute.py lines): hourly Welch (0.01 Hz), prominence over the +-0.25 Hz running median, median over
hours >= 4 (present most of the day) and narrow (half-height within +-0.05 Hz), 1-8 Hz. 5,237 lines on 1,360 of 4,700
station-days: B007 3.49/6.9/7.13 Hz (1,397), N060 6.6-7.5, LRIV/OPC/B04A 3.49 (the same line at 4 stations = an
array-coherent line that muting cannot remove), B013/B943/W030 3.8-3.9, STW 1-5 Hz harmonics, SMW 2.0; 1.09 Hz x87,
5.47 x12, 3.91 x36. Notched in load_day (iirnotch, 0.06 Hz wide, zero-phase, each valid run >= 60 s, after muting).
Check B007 2011-08-10: line prominences 346 -> 0.44 etc.; 14 % of the energy was the lines.
**Rerun** `scripts/nw_mute_queue.sh` (lock locks/nw_mute.lock, pid results/NW_MUTE_PID, markers results/NW_MUTE_{DONE,
FAILED,HEARTBEAT}, log logs/nw_mute_queue.log): ab_auto (NW, --cols window) -> ab_cluster (--headline-from
NW/clusters/summary.json: average linkage, median+8MAD, min 2) -> ab_grow_fast (25-station pool, 14 iter, tol 0.01,
prep 6 / restack 8) -> ab_hypo (boot 60) -> ab_restack28 -> scripts/ab_review_sheets.py (new; layout of
hypo_qc28/review: 6 families/sheet by detections, stations by distance from the highest-SNR station, N and Z).
All under results/ab_royer/NW_mute/. Smoke (2008-05-08 13h): frac_valid 0.990, 3 detections (run 1: 6).
Launched 22:04:50Z, pid 2629125.

## 2026-09-27 ~06:00-06:40Z — NVI shortfall diagnosis (2,039 vs 4,267 first detections, -52 %)
Paper sec.3.1 (verbatim): "Network autocorrelation for the NVI dataset was applied to two tremor episodes (September
05-11, 2006, June 13-18, 2007) for different combinations of 7 of the 10 three-component stations VI10, VW03, VI11,
VW02, VI52, VI53, VI05, VI08, VI06, VI04." ... "After processing all data, we had selected 4915, 3306 and 4267 initial
detections". 4,267 over 13 days = 13.7 per hour (NW 2.7/h, SW 3.1/h).

**1. Data availability of the 10 named stations** (cache40_NVI.csv; frac_valid of the day):
| episode | VI10 | VW03 | VI11 | VW02 | VI52 | VI53 | VI05 | VI08 | VI06 | VI04 |
|---|---|---|---|---|---|---|---|---|---|---|
| 2006-09-05..11 | 1.0 | 1.0 | 1.0 | 0.93-1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| 2007-06-13..18 | 1.0 | none | 1.0 | none | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 0.25-0.74 |
All HH? 100 Hz, ZNE (no rotation needed). Autocorrelation is per channel with normalised cc, so gain and
orientation do not enter. **All 312 hours were processed with 7 three-component stations** (none dropped; the
<4-station floor never fired): 2006 VI10 VW03 VI11 VW02 VI52 VI53 VI05 (166 h; 2 h with VI08 for a VW02 gap);
2007 VI10 VI11 VI52 VI53 VI05 VI08 VI06 (144 h). VW01/02/03 2007: no data at EarthScope, NRCan, GEOFON, ORFEUS
(204); NRCan's PO holdings do not include the VI/VW stations. So the missing station-days do NOT cut the
autocorrelation.
**Missing pool station-days (growth/location only):** the 24 h NRCan requests failed with intermittent 404s.
`scripts/ab_refetch_nvi.py` (4 h pieces, all providers; then `--partial`: 1 h pieces, 6 tries, merged with the file on
disk) recovered CN.MAYB 9, CN.WOSB 9, CN.NCRB 7 station-days (EHZ, 17-100 % of the day; `refetch_NVI*.csv`).
Still absent everywhere: TASB 2007 (C8 epoch ends 2006-09-29), VI03/07/09/54 2006-09, VI56 4 days, VI50 2007, VW01-03 2007.

**2. The "7 of the 10" readings, two days** (`NVI/diag_7of10/per_hour_*.csv`, script `nvi_combo.py`; same
autocorr_hour, row MAD, pair suppression; the union of combinations gets the same +-15 s 2-D pair suppression):
| day | 3-C stations | 7-subsets | first 7 listed (current) | mean over 7-subsets | best 7-subset | union of all 7-subsets | all available (not "7") |
|---|---|---|---|---|---|---|---|
| 2006-09-08 | 10 | 120 | 16.7 /h | 15.9 /h | 31.2 /h | 197 /h | 20.0 /h |
| 2007-06-14 | 8 | 8 | 12.2 /h | 12.1 /h | 13.4 /h | 15.5 /h | 12.0 /h |
Published average 13.7 /h. Any single 7-station combination gives the same rate as the first-7 rule (15.9 vs 16.7,
12.1 vs 12.2), so WHICH 7 is not the cause. The union over all 120 subsets gives 197/h in 2006 (x12; mostly each
subset's own 8-MAD tail). The paper says "different combinations", not "all combinations", and its NW/SW
autocorrelation uses 7 stations; the wording is met by the station set changing with availability (VW02/VW03 in 2006,
absent in 2007). **Reading kept: per hour the first 7 of the listed 10 with 3-C data** (unchanged).
**Where the shortfall is:** by day, first-7 detections/h: 2006-09-05 4.0, 09-06 4.0, 09-07 7.6, 09-08 16.7, 09-09 9.1,
09-10 5.2, 09-11 1.1, 2007-06-13 3.6, 06-14 12.2, 06-15 7.5, 06-16 9.9, 06-17 2.2, 06-18 2.0. The two strongest days
match or exceed the published average; the deficit is in the weak days at the start and end of each episode. Nothing in
data availability, station choice, channels or the hour floor explains it; not resolved from the paper.

### NVI muted run launched + the shortfall resolved (06:44Z)
`scripts/nvi_mute_queue.sh` = copy of sw_mute_queue.sh with --sub NVI (outputs results/ab_royer/NVI_mute/, markers
results/NVI_MUTE_{DONE,FAILED,HEARTBEAT}, lock locks/nvi_mute.lock, log logs/nvi_mute_queue.log), plus two first steps:
wait for ab_refetch_nvi.py, then ab_cache40.py --sub NVI (builds only the re-fetched station-days). Mute K=15, qmax 2,
qlong 1.3, pad 3 s, lines notched; sites = stations < 3 km apart (ab_mute.sites: WOSB+VI09 is the only NVI pair).
Launched 06:29:49Z with `nohup setsid -f`, **queue pid 2824910**.
**Autocorrelation with the muted cache: 4,235 first detections (Royer 4,267, -0.7 %)**, same stations and 7-of-10
reading as the unmuted 2,039. Per day, unmuted -> muted: 09-05 95->143, 09-06 96->248, 09-07 182->305, 09-08 401->385,
09-09 218->428, 09-10 124->219, 09-11 27->177, 06-13 87->619, 06-14 292->543, 06-15 179->416, 06-16 237->341,
06-17 53->190, 06-18 48->221. Muted fraction of the 10 autocorrelation stations: VI06 14.2 %, VI52 11.0, VW02 9.8,
VI11 7.7, VI05 7.6, VI04 7.1, VI10 5.3, VI08 4.3, VW03 0.6, VI53 0.4 (NW typical ~1 %). **Cause of the -52 %:
single-station bursts in the POLARIS records** (not missing data, not the 7-of-10 rule): a burst on one of the 7
stations makes that station's channels uncorrelated for every pair in the span and dilutes the network average,
mostly on the weak-tremor days (the strongest day 09-08 is unchanged, 401 -> 385). Muting restores the published count.

## 2026-09-28 — Merlin step 1: no-merge growth test (176 seeds) → PASS; step 2 launched
- `scripts/ab_grow_fast.py --no-merge` added (skips merge_duplicates; default unchanged; .bak_nomerge kept).
- Seeds: 176 NW_mute it0 cluster stacks with n_sites >= 2 (nw_mute_seed_scores.csv) → test176/seeds176.npz; queue
  scripts/nw_nomerge_test176.sh (AB_MUTE_DIR=NW_mute/mute), 3 iterations, then ab_restack28 of det_it3.
- Detections (raw = kept, no merge): it0 133,507 · it1 454,808 · it2 711,498 · it3 973,974 (per-template median 212→1,778,
  max 11,905→31,925; template SNR median 2.04→8.13). Still rising ~1.4×/iteration at it3.
- Score (scripts/ab_nomerge_score.py = ab_famscore with two site groupings): it3 templates with >=3 impulsive sites =
  **65/176** (both NW-original B013/W030/B943 grouping and the <3 km grouping) — bar >=40, failure <15 → **PASS**
  (merge run: 8 of the 33 >=3-site seeds became real). By seed outcome: abs_into_lfe 35/36, lfe_root 15/15,
  abs_into_noise 15/113, noise_root 0/12 now >=3 sites. Of the 33 >=3-site seeds, 15 kept >=3 sites, 17 fell to 0-1.
- Caveat: many of the 65 are duplicates of each other. The final dedupe (scripts/ab_final_dedupe.py: shared fraction at the
  modal lag ±0.75 s within ±15 s, minus the same at +60 s shift, > 10 %; representative = higher site score, tie more
  detections) on it3: 591 duplicate pairs → 96 representatives, **30 distinct with >=3 sites** (test176/dedupe_it3/).
- Step 2 launched 04:29Z: scripts/nw_nomerge_queue.sh (nohup setsid -f), queue pid 2940177, grow pid 2940190,
  lock locks/nw_nomerge.lock (only the queue shell holds it), markers results/NW_NOMERGE_{DONE,FAILED,HEARTBEAT,RUNAWAY},
  log logs/nw_nomerge_queue.log, outputs results/ab_royer/NW_nomerge/full/. Chain: grow 14 it (--no-merge, tol 0.01) →
  restack28 (all) → score → final dedupe → ab_hypo + restack28 + review sheets + score on the representatives.
  Runaway guard: kills grow pid if any iteration's n_det_raw > 5 M (writes NW_NOMERGE_RUNAWAY, queue FAILED=runaway).

## 2026-09-28 ~14:45Z — FAIR per-iteration merge (user-approved) → NW_fair queue launched
- `scripts/ab_grow_fast.py --merge fair` (default `--merge most` = old merge_duplicates, unchanged; .bak_fairmerge kept).
  merge_fair, applied to each iteration's raw detections before the restack:
  (a) duplicate test = ab_final_dedupe logic: b-detections within ±15 s of each a-detection, 0.5 s bins, modal bin ±1 bin
      (±0.75 s), fraction = count / min(n_a, n_b); control = same with b shifted ±60 s (larger of the two signs);
      duplicate if excess > 10 % (pairs pre-pruned: total ±15 s count must exceed 10 % of min(n));
  (b) representative = higher site-count score (ab_famscore logic: impulsive = envelope max/median > 5 and half-max
      duration < 3 s; sites = stations < 3 km apart, ab_mute.sites) on the template's CURRENT stack band-passed 2–8 Hz
      (4th-order Butterworth, zero-phase, applied to the stack — no data re-read); tie → more detections;
  (c) greedy star in descending (score, n_det): join the representative with the largest excess, never chain;
  (d) a duplicate pair with one member ≥3 sites and the other <3 is never merged (either way) — kept separate and logged
      in grow/fair_straddle_it<N>.csv (count in iters.csv n_straddle);
  union: member times shifted onto the representative by the modal lag, then one detection per 15 s per family (strongest cc;
  = the matched filter's own tpl-len declustering; the old merge used 5 s). Scores per iteration: grow/fair_scores_it<N>.csv.
- Check of the fast score (stack band-passed) vs the per-segment 2–8 Hz restack score on the 299 converged NW_mute families:
  ≥3-site agreement 63/65, 2 fall below 3 (3→1, 4→2), no low family promoted. On the 2-member it0 seed stacks the fast score
  is looser: 119 of 3,080 seeds ≥3 sites (the nw_mute_seed_scores.csv reference had 33) — affects it0 representative order only.
- Cost (2 days 2011-08-21 + 2008-05-08, 3,080 templates): matched filter ~7 s/day (as before), score 9.4 s, merge 1.6 s;
  on the full it0 no-merge raw (1.30 M detections) score 9 s + merge 4 s → 2,362 duplicate edges, 638 straddle pairs kept
  separate, 3,080 → 2,635 families, 824 k detections. So ~15 s per iteration added to ~22 min of matched filter + restack.
- Queue `scripts/nw_fair_queue.sh` (nohup setsid -f; lock locks/nw_fair.lock held only by the queue shell; heartbeat subshell
  closes fd 9 and is killed in the EXIT trap; runaway guard: n_det_raw > 5 M in any iteration → kills grow + its children,
  writes results/NW_FAIR_RUNAWAY): grow (--merge fair, 14 it, tol 0.01, from NW_mute/clusters/tpl_it0.npz,
  AB_MUTE_DIR=NW_mute/mute) → ab_hypo (boot 60) → ab_restack28 → ab_review_sheets → ab_nomerge_score →
  NW_fair/hypo_28/review/auto_score.csv. Markers results/NW_FAIR_{DONE,FAILED,HEARTBEAT,RUNAWAY}; log logs/nw_fair_queue.log.
  Queue pid 3006027, grow pid 3006040 (results/NW_FAIR_PID, NW_FAIR_GROW_PID).
- Note: 8 orphaned restack workers of the stopped no-merge run (pids 2956299..2958758, ppid 1, ~32 GB RSS each, idle) are still
  alive; not killed (not asked).

## 2026-09-28 21:20Z — NW_lag15: clustering with the FULL lag window (Merlin step 3, user-approved)
- ab_cluster.py: new `--lag-mode {data,zeropad}` (default data = unchanged; `--lag-max` already existed, default 2 s) and `--simmat-only`.
  Backup scripts/ab_cluster.py.bak_lag2. Why a new mode: the old super traces are cut win+2*lag from the CONTINUOUS record and every
  shifted window is re-normalised, so at |lag| ~15 s the reference is compared with a different, adjacent 15 s window (not the detection).
  `zeropad` = literal linear (non-circular) normalised xcorr of the two finite 15 s windows: each channel window unit-normalised once, zero
  padded by lag_max each side, overlap sum per lag, channel-mean over shared channels; channels never slide into each other. Same LAG sign
  convention (verified on a synthetic pulse: lag +280 samples, cc 1.0). 1,201 lags on the GPU in ~40 s (8,107 x 21 ch).
- Step-2 check (NW_mute autocorr detections + masks; same headline rule average, median+8xMAD, min size 2):
  median sim 0.054 (2 s: 0.050), cut 0.8924 (0.8904) -> **3,010 clusters (2 s: 3,083); 3,006 stacked templates (3,081)**.
  Size distribution 2/3/4-5/6-10/11-20/21-50: 15 s 2669/171/93/42/21/14 (max 50) vs 2 s 2749/168/99/39/19/9 (max 36).
  Detections whose ancestry ends in the 66 real NW_mute families (1,332 dets; traced through families_it0..14 rep maps):
  in clusters >=3 48 % (2 s 44 %); noise-ancestry dets (5,866) 21 % (19 %). Real families' detections spread over median 4.5 lag-15 clusters
  (2 s: 4 initial templates) -> **no coalescence**: full lag barely changes the clustering. Script: scratchpad lag15_check.py (not in repo).
- Queue scripts/nw_lag15_queue.sh (pid in results/NW_LAG15_PID), launched 21:21Z: cluster (lag 15 zeropad) -> ab_grow_fast --merge most
  (original) 14 it, restack 4 -> ab_hypo (GSC, boot 60) -> ab_restack28 --workers 4 -> review sheets -> ab_nomerge_score (orig + 3 km sites).
  Runaway guard on post-merge n_detections > 5 M. Outputs results/ab_royer/NW_lag15/; markers results/NW_LAG15_{DONE,FAILED,HEARTBEAT,RUNAWAY}.
  Clustering done 21:42Z; growth started 21:42Z.

## 2026-09-28 22:20Z — OLY2010 (Chestler 2017 reproduction) relocated with USGS CVM v1.7 (USGS_NW) vs GSC, Hyp2000
- Script `scripts/ab_oly2010_cvm.py` (new; read-only on inputs), log `logs/ab_oly2010_cvm.log`, outputs `results/ab_royer/OLY2010_cvm/`
  {picks.csv, loc_USGS_NW.csv, loc_GSC.csv, compare.csv (+ *_null.csv), score_summary.csv, station_codes.json, run_*/}.
- Inputs: existing u_locate long stacks `results/u_locate/wa_{real,null}_stacks.npz` (1-8 Hz, 20 sps, 44 s, from frozen det_it0; 1,479 real /
  1,471 site-shift-null families, 22 sites = 10 beams + 12 singles, coords `results/o_oly2010/sites.csv`). No restack, no detection/growth.
  L9 = `results/u_locate/wa/wa_{real,null}_ranked.csv`. 280 rich = sf_top (dedupe_n_partners >= 20); null comparator = top 280 by partners.
- Picks ab_hypo.make_picks (A.FS set to 20); Hyp2000 via ab_hypo.write_hyp_inputs/run_hyp, R 1.5, free depth (ZTR 35), MIN 4;
  USGS_NW Vp/Vs 1.823 (LVZ merged, ab_hypo_cvm), GSC 1.73. Station codes mapped to <=5 chars (xPBB005->XPB05, xPBB013->XPB13, xBH->XBH ...).
- R>=1.5 is harsh on these stacks: of 32,517 real station stacks only 1,250 P / 833 S pass -> 500 families with >=2 picks -> Hyp2000 located
  110 (USGS_NW) / 109 (GSC) of 1,479 (null: 33 / 35 of 1,471). Median gap ~282-305 deg, num 5.
- Scores vs the 69 Chestler-derived coordinates (u_locate scorer verbatim + one-to-one Hungarian + rotated/shifted null; 1000 draws):

| locator | pop | n_loc | coords_hit 5/10/20 | fam<=10km | o2o-3D 3/5/10 (rot-null mean 10) | dep_med | frac 25-45 | erh/erz med | Royer NW frac 3/5/10 |
|---|---|---|---|---|---|---|---|---|---|
| L9 | real all | 788 | 16/61/69 | 65 | 6/8/28 (15.9) | 72.0 | 0.19 | boot 35.6 km | 0.66/0.80/0.97 (602 in region) |
| L9 | null all | 230 | 0/3/31 | 2 | 0/0/1 | 76.0 | 0.09 | — | 0.04/0.07/0.21 |
| GSC | real all | 109 | 4/47/69 | 17 | 1/2/5 (6.1) | 35.0 | 0.58 | 9.4/14.7 | 0.07/0.14/0.28 |
| GSC | null all | 35 | 1/9/52 | 3 | 0/1/3 | 35.0 | 0.57 | 81.6/77.1 | 0.05/0.08/0.11 |
| USGS_NW | real all | 110 | 22/56/69 | 24 | 4/6/12 (7.5), p(3,5 km)=0.01 | 34.4 | 0.65 | 8.9/12.2 | 0.12/0.20/0.34 |
| USGS_NW | null all | 33 | 1/3/34 | 2 | 0/0/1 | 35.0 | 0.61 | 85.4/86.8 | 0.05/0.07/0.07 |
| L9 | real 280 rich | 171 | 0/3/45 | 1 | 0/0/1 | 76.0 | 0.05 | | 0.03/0.09/0.19 |
| GSC | real 280 rich | 8 | 0/1/3 | 1 | 0/0/0 | 24.4 | 0.50 | 83.6/35.1 | 0.01/0.01/0.01 |
| USGS_NW | real 280 rich | 8 | 0/0/0 | 0 | 0/0/0 | 31.2 | 0.62 | 77.7/71.2 | 0.01/0.01/0.01 |

- Same-family control: L9 restricted to the 110 USGS_NW-located families: coords_hit 4/33/69, fam<=10 km 18, med nearest 3-D 19.2 km
  (USGS_NW on them: 22/56/69, 24, 22.0 km). USGS_NW vs L9 on 98 common families: median epi shift 24 km, dz -5.9 km; USGS_NW vs GSC
  (108 common): median epi 7 km, dz -0.45 km. Hyp2000 `edge` not defined -> u_rail column 0 for GSC/USGS_NW (not a rail measurement).
- Royer NW truth: L9's 0.97 at 10 km comes from 602 families in the 30 km region (the one-to-one fraction saturates with count; null 0.21);
  Hyp2000 sets have only ~90 families in-region, so their fractions are count-limited, not comparable to L9's.
- Reading (provisional): where Hyp2000 locates, USGS_NW beats GSC on every Chestler metric and pulls depths into 25-45 km (65 %) vs L9's
  railed 72 km median, and beats L9 on the same families; but the literal R 1.5 recipe locates only 7 % of families and 8 of the 280 rich
  (with erh ~80 km, gap >300) — the 280 remain effectively unlocated by Hyp2000.
