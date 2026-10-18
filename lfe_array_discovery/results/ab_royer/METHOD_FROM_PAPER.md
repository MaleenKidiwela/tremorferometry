# Royer & Bostock (2014, EPSL 402:247-256) -- the method, read independently from the paper

Written 2026-09-25 by the audit agent from `notes/1-s2.0-S0012821X13004780-main.pdf` (sec.2, 3,
3.1-3.3, 4, 5, 6) and `notes/royer2014_supp.pdf` (Table S1) BEFORE reading `AB_ROYER_PLAN.md` or the
`ab_*.py` code.  **SILENT** = the paper gives no value / no rule; those are swept and reported as
columns, never guessed constants and never gates.

## Recipe (Washington: NW and SW subarrays, processed independently -- sec.2)

1. **Data** (sec.2, 2.1). CAFE (XU) + PBO + TA + PNSN (UW). NW = 25 stations on the E flank of the
   Olympics; SW = 25 stations SW of Puget Sound (Fig.1D only; no station list for the 25). NW ETS
   recorded 2007, 2008, 2009, 2010, 2011; SW only 2007, 2008.
2. **Pre-processing** (sec.3). 24-hour segments; band-pass **1-8 Hz**; resample to **40 sps**.
   Filter type/order, response removal, gap handling: **SILENT**.
3. **Network autocorrelation** (sec.3.1).
   3a. Stations: NW = **B001, SQM, BS11, W020, W040, GNW, PL11** (7, 3-C); SW = **N050, S030, W070,
       N060, PL11, S040, N070** (7, 3-C).
   3b. Days: NW **2007-01-17..30, 2008-05-01..23, 2010-08-07..24, 2011-08-04..23**; SW
       **2007-01-14..02-01, 2008-05-01..25**.
   3c. Hour-by-hour. Each hour -> **15 s windows lagged by 0.5 s** (**7170** windows/hour).
   3d. Per single (station, component) channel: correlate **each window with all other windows of
       the same hour** (normalised cc of window pairs; the lag axis is the 0.5 s window grid).
       Whether the trivially overlapping near-diagonal pairs are excluded: **SILENT** (necessary,
       since a window correlates with itself / its 0.5 s neighbours at ~1).
   3e. Stack the per-channel matrices over all stations and components of that hour ->
       "network autocorrelation". Mean vs sum: **SILENT** (equivalent under a MAD-relative bar).
   3f. Detection: network autocorrelation **> 8 x MAD** -> the window PAIR is an LFE detection. The
       MAD population (whole hour matrix vs per row) and whether the bar is median + 8 MAD or 8 MAD:
       **SILENT**.
   3g. Cull: "exclude overlapping detections and retain only those with high SNR". Overlap rule and
       SNR definition/threshold: **SILENT**. Published yield: **4915 (NW), 3306 (SW)** detections
       (whether pairs or events: **SILENT**).
4. **Cluster analysis** (sec.3.2).
   4a. Per detection, concatenate all channels into one "super" trace (window length: **SILENT**,
       presumably the 15 s window).
   4b. Cross-correlate all super traces pairwise -> max cc and its lag (max lag: **SILENT**).
   4c. Max cc -> similarity matrix -> **hierarchical cluster analysis** (Rowe et al. 2002). Linkage,
       cut, minimum cluster size: **SILENT**.
   4d. Per cluster, shift (by the lags) and stack **all channels (stations/components) available**
       -> initial template. Stack normalisation: **SILENT**. Yield **224 (NW), 70 (SW)**.
5. **Network cross-correlation / growth** (sec.3.3).
   5a. Per template choose **7 to 10 high-SNR 3-C stations available within the template** (SNR
       definition: **SILENT**; the pool is the template's channels, i.e. potentially the 25-station
       array of sec.2.1).
   5b. Scan **all available tremor dates** (which dates beyond the sec.3.1 episodes, e.g. the 2009
       NW episode, and their day bounds: **SILENT**).
   5c. **Summed** network cc (per-channel normalised cc summed over channels, Gibbons & Ringdal 2006;
       Shelly 2006) **> 8 x MAD over a 24-hour period** -> detection. Template length, declustering:
       **SILENT**.
   5d. New template = stack of **all** newly registered detections, each normalised to **unit
       maximum amplitude across the three components**.
   5e. Repeat "**several**" iterations (count: **SILENT**). Yield **122 (NW), 54 (SW)** final
       templates, each with "100's to 1000's" of detections.
6. **Location** (sec.4, 5, 6-i/ii). Templates normalised to unit max and **90-degree phase rotated
   (Hilbert)**; time only unambiguous zero-phase P and S arrivals; locate with **Hyp2000** then
   **hypoDD** (catalogue picks + cc differential times). Table S1 (IDAT 3; 5 sets x NITER 8 = 40
   iterations): WTCCP/WTCCS 0.01,0.01,1,1,1; WRCC -9,-9,-9,5,5; WDCC -9,-9,2,2,0.5; WTCTP/WTCTS
   1,1,0.01,0.01,0.01; WRCT -9,-9,6,6,6; WDCT -9,8,8,8,8; DAMP 50,50,40,40,40. Velocity model,
   Vp/Vs, picking method: **SILENT**. Checks: depths **31-46 km (NW)**, **32-49 km (SW)**; Hyp2000
   sigma_h/z 2.0/1.3 km (NW), 3.7/3.5 km (SW); hypoDD keeps 77 % (NW), 72 % (SW).

## Audit -- my recipe vs `AB_ROYER_PLAN.md` and the `ab_*.py` code

Agreement: 1-8 Hz, 40 sps, 15 s / 0.5 s, hour blocks, the 7+7 stations, all episode dates, 8 x MAD
(hour for 3.1, day for 3.3), hierarchical clustering with linkage/cut swept, unit-max-3C restack of
all detections, iterate-to-convergence, published numbers as reference lines only.

Disagreements found and FIXED (the paper wins):

| # | where | what was wrong | fix |
|---|---|---|---|
| A1 | `ab_auto.py` pass 2 | the window-pair matrix is symmetric, so every pair was registered twice ((i,j) and (j,i)); detections doubled | keep only j > i (upper triangle); the MAD pool is unchanged (symmetric) |
| A2 | `ab_cluster.py stack_clusters` | shift sign inverted: `LAG[ref,k]=l` means member k aligns at `t_k + l/fs`, code cut at `t_k - l/fs` -> the initial templates were stacked MIS-aligned (2x the lag) | cut at `t_k + l/fs` |
| A3 | `ab_cluster.py silhouette` | pure-Python O(Nd x K) loop, run for 136 sweep configs -> hours-to-days at Nd ~ 10^4 | vectorised (sparse one-hot GEMM), identical definition |
| A4 | `ab_grow.py cut_templates` | invented gate `MIN_PEAK_MAD=3` (a station enters a template only if its stack peak > 3 MAD) plus `min_sta=7` scannability: with a 7-station pool ANY missing/weak station made the template unscannable. Paper: "choose 7-10 high-SNR stations available" = a RANKING | stations RANKED by stack SNR, top <=10 used; peak SNR kept as a column; no SNR cut |
| A5 | `ab_grow.py cut_templates` | per-station cut 2 s before the envelope peak (not in the paper; at 31-49 km S-P ~4-5 s, so P fell outside the template) | default `--cut common`: the template is the sec.3.1/3.2 15 s detection window itself, common to all stations (moveout preserved); envelope cut kept as `--cut envelope` ablation |
| A6 | `ab_grow.py matched_filter` | days with fewer than 7 stations were SKIPPED and samples needed all 7 valid (`min_sta=7`), so with the 7-station pool one gap/one missing station (e.g. SQM 2011, EHZ-only) removed the whole day; paper has no such rule | availability floor `--min-sta-valid 4` = the same floor stage 3 uses per hour; per-sample channel count is a column |
| A7 | `ab_grow.py matched_filter` | cc normalised jointly over the 3 components of a station and averaged over stations; paper's "summed network cross-correlation" (Gibbons & Ringdal; Shelly 2006) = per-CHANNEL normalised cc summed over channels | per-channel normalisation, divided by the per-sample valid channel count (`--net mean`) or summed (`--net sum`) |
| A8 | `ab_fetch.py` / `ab_cache40.py` | `BH?,HH?,EH?,SH?` returns >1 band on transition days (e.g. UW.SQM 2011-04..07 BH+EH, XU.BS11 2008 BH/HH) -> "mixed sampling rates" error | pick ONE 3-C band per station-day (preference BH > HH > EH > SH), recorded |
| A9 | `ab_null.py` | null ran `--grow-iters 4`, `--boot 30`; real runs `--max-iters 8 --tol 0.01 --boot 60` -> not the identical chain | null uses exactly the real settings |
| A10 | `ab_cluster.py CUTS` | cut sweep stopped at distance 0.75 (similarity 0.25); smoke similarity median 0.051, p99 0.104 -> EVERY configuration was all singletons, 0 templates | sweep spans distance 0.10-0.98 |
| A11 | `ab_cluster.py` headline | max-silhouette headline picked a degenerate config (1 cluster of 2) | headline = average linkage cut at the paper's own significance rule, similarity > median + 8 x MAD of the off-diagonal similarity matrix (sec.3.1/3.3 rule; linkage/cut SILENT). Max-silhouette and nearest-published configs stay reported (`summary.json`, `clusters.csv` columns) |
| A12 | `ab_common.null_remap` | per-station hash shifts collided (W020 and GNW got the same shift under k=1) -> those stations kept their TRUE mutual coherence inside the null | distinct shifts per station: seeded permutation of 1..n-1 per (realisation, episode) |
| A13 | `ab_cache40.py` | UW.SQM 2011-08 is EHZ-only (FDSN: 3-C BH? ends 2011-07-02) -> station-day errored and SQM dropped whole | vertical-only station-day cached with zero horizontals (`built_z_only`); every stage drops zero-energy channels, so SQM contributes its one real channel |
| S1 | `ab_locate.py RoyerMod` | `AttributeError: station_corrections` -- stage 6 crashed | method added, same convention as `l_common` (x_station_corrections.csv; missing -> 0, flagged) |
| S2 | `ab_score.py` | `to_markdown` needs `tabulate` (not installed) -> ImportError would kill the queue after null 1 | fixed-width fallback |
| S3 | `ab_grow.py restack` | template-outer loop reloaded each day once per template (~#templates x #days loads per iteration) | date-outer loop; identical sums (smoke: identical 1237/2314 detections) |
| S4 | `ab_auto.py` | empty hours wrote object-dtype parquet columns | typed empty frame |

### Round 2 (coordinator review, 2026-09-25): the former "recorded deviations" removed

Correction: an earlier draft said the 7-station growth pool was "by instruction". That was the
coordinator's default, NOT a user instruction, and it is now replaced by the paper's pool.

| # | where | what the paper says / what was wrong | fix |
|---|---|---|---|
| B1 | `ab_common.POOL`, `ab_grow.py`, `ab_cluster.py` | sec.3.3: 7-10 high-SNR stations "available within a template", i.e. drawn from the sec.2.1 **25-station array**. sec.3.2: stack "all channels available". The code used only the 7 autocorrelation stations | pools reconstructed from Fig.1D (`POOL_RECONSTRUCTION.md`: NW 25 labels, SW 28 labels vs the paper's 25, documented). Cluster stacks use the whole pool. Growth ranks the pool's 3-C stations by template SNR and keeps the top 10 (at least 7 when available) |
| B2 | `ab_tremor_days.py`, `ab_common.scan_days` | sec.3.3 "scanning through all available tremor dates", including the 2009 NW episode (sec.2.1). The code scanned only the sec.3.1 days | scan days = sec.3.1 days, plus PNSN in-region catalogue days (>= 1 h of tremor within 25 km of the pool) from 2009-08-06, plus the PNSN published ETS windows before that. NW 180 days, SW 56 |
| B3 | `ab_hypo.py` (new stage 6) | sec.5/6: Hilbert-rotated, unit-max templates; zero-phase P/S picks; **Hyp2000 then hypoDD, Table S1 exactly**. The code used the project's own locator | Hypoinverse-2000 v1.40 (USGS) and HYPODD 1.3 built from source under `/home/jovyan/tools/` (gfortran in `tools/gf`). Hyp2000 reproduces its USGS test event exactly; hypoDD reproduces example1 to <= 62 m; a synthetic 20-event test recovers Hyp2000 epicentres to 0.18 km and hypoDD to 0.46 km. Table S1 is used verbatim (IDAT 3, ISOLV 2 LSQR, 5 x NITER 8, all weights/cutoffs/DAMP). SILENT, swept as columns: Vp/Vs 1.73/[1.78]/1.85 with the Crosson (1976) PNSN model (the paper names no model); pick-unambiguity ratio 1.2/[1.5]/2.0; ph2dt at its documented defaults. Per family: Hyp2000 depth/ERH/ERZ, a 60-draw station bootstrap, and hypoDD retention. The calibrated L9/L12 locator remains only as a secondary `L9_*` column |
| B4 | `ab_cluster.detections` (A14) | sec.3.1 "culled to exclude overlapping detections". The 2-D pair suppression let many pairs share one window, which inflated SW to 7981 pairs vs 3306 | **strict time-overlap cull**: a pair is dropped if either of its windows overlaps a window of a stronger kept pair. Counts on the stopped run's autocorrelation: NW 5025 -> 3258 pairs (Royer 4915), SW 7981 -> 2936 (Royer 3306). Both counts are written to `summary.json`. The SNR half of the cull remains a column (no threshold given) |
| B5 | `ab_cluster.py` headline (A15, supersedes A11) | linkage, cut and minimum cluster size are SILENT (sec.3.2; Bostock et al. 2012 sec.2.1 likewise: 5775 detections -> 140 families, Rowe et al. 2002, no numbers). median + 8 MAD with min size 2 re-found the autocorrelation pairs themselves (2640 / 2192 "templates" vs 224 / 70) | full linkage x cut x min-size (2-30) table with a `nearest_published` column (`sweep_long.csv`). Headline = the configuration with the **largest real-minus-null template count** against the site-shift null-1 sweep of the same stage. This is a data-internal rule that never looks at 224/70. The null realisations copy the real headline (`--headline-from`) |
| B6 | `ab_common.null_remap` | extra scan days outside the sec.3.1 episodes were NOT shifted under the null | null shift groups = each (subarray, year)'s scan days. The 7 autocorrelation stations always get distinct shifts; pool stations collide only when a group has < 25 days |

Where the paper still cannot be followed exactly (stated, not substituted):
* Inter-ETS tremor days before 2009-08 cannot be identified: there is no public machine-readable
  catalogue. The ETS windows are used instead.
* SW pool: Fig.1D shows 28 labels for the paper's "25". All 28 are kept, and data availability
  decides.
* Velocity model, Vp/Vs, picking criterion, ph2dt settings: not in the paper. Swept or defaulted
  as above.
* hypoDD 1.3 is the public release. Which hypoDD version Royer used is not stated.

### Extension (user-approved, separate outputs under `results/ab_royer/ext/`)
The extension keeps the literal autocorrelation detections and initial templates unchanged. It
adds public stations (UW, PB, TA, XU, CN, US, CC, C8, GS, YW) within 75 / 100 / 150 km of the
subarray, ONLY at growth and location. Dense arrays XG/ZK/Z3 and OBS 7D/ZU are excluded, as
stated. Scripts: `ab_ext_pool.py`, `ab_grow.py --pool ext:R`, `ab_hypo.py --pool ext:R`,
`ab_ext_score.py`. The null gets the identical treatment. Markers: `AB_ROYER_EXT_DONE` /
`AB_ROYER_EXT_FAILED`.

## Smoke test (2026-09-25, `results/ab_royer/smoke/`)

NW subarray, 7 stations, 2010-08-15 and 2010-08-16 (B001 raw absent 08-16 -> 6 stations that day).

| stage | result |
|---|---|
| 1 fetch | UW.SQM 2010-08-15: BHE/BHN/BHZ 50 Hz, 14.8 MB, 2.6 s. UW.SQM 2011-08-10: EHZ only (100 Hz) |
| 2 cache40 | 7/7 built for 08-15 (frac_valid 0.99996-1.0; B001 EH1/EH2 rotated to N/E), 6/7 for 08-16; SQM 2011 `built_z_only` |
| 3 auto window | 48 h in ~10 s GPU; 7171 windows/hour; thr = med + 8 MAD = 0.114-0.154 (8 MAD = 5.4 sd); **213 pairs** after the overlap cull (08-15: 179 = 7.5/h; 08-16: 34 = 1.4/h; Royer's average 4915/1800 h = 2.7/h); cc 0.114-0.202, no NaN/inf, no pair with j<=i, no pair closer than 15 s |
| 3b auto sample | 1 h in ~5 s; 281 pairs -> 21 after the cull |
| 4 cluster | 338 culled detection times; super traces 20.4/21 channels; similarity median 0.051, MAD 0.0073 -> 8xMAD cut at distance 0.891; pair partners' similarity >= their autocorr cc (consistency check); **136 clusters >= 2 members** -> 136 initial templates, all finite |
| 5 grow | iter 0: 136 templates -> 1237 detections (7 sta: 837, 6 sta: 400); iter 1: 2314 detections, per template median 8, max 98; cc 0.108-0.83 (max = self-detection of a 2-member template), nmad >= 8.0, no NaN; ~10 s per day per 136 templates |
| 6 locate | 120 families stacked, 824 station stacks, P found 741; L9_xcdiff located 56 (rail 0.09), depth median 24 km (p05 7.5, p95 53.5), boot_km median 6.9; L12_comb 60 |
| 8 score | stage_counts / separation_curves / SUMMARY.md written (scratch copy) |

Smoke depths are NOT a check: two days of pair-templates dominated by single-day noise. WATCH in the
full run: if the final-template depths stay far shallower than 31-46 km, the 7-station
locator is the first suspect (our bug, not "irreproducible").

## Smoke test, round 2 (2026-09-25; NW 2010-08-15/16; `smoke/`)

| check | result |
|---|---|
| pool fetch/cache | 14 of 25 NW pool stations have data on these days. OPC, B04A, FACB and N050/N060/N080 were not operating; DOSE returned no data. STW, BLN, HDW and SMW are EHZ-only (`built_z_only`) |
| strict cull | 213 pairs -> 150 pairs -> 300 detection times |
| cluster (pool stacks) | 135 templates, stacked over the 25-station pool |
| growth, pool | every template uses 10 stations (BS11, W020, W030, W040, B001, B013, C04A, W060, B943, PL11, ...). Iteration 0: 3967 detections; iteration 1: 13734. cc 0.095-0.815, no NaN |
| Hyp2000 check | recovers 3 synthetic hypocentres to <= 0.05 km. With 20 synthetic events: Hyp2000 0.18 km, hypoDD (Table S1) 0.46 km epicentral |
| Hyp2000 + hypoDD, smoke templates | headline (Vp/Vs 1.78, R 1.5): 28 of 53 templates located, depth median 34.3 km (p05 4.2, p95 41.9), 71 % inside 31-46 km, hypoDD kept 3. These are two days of noisy pair templates, so they do not test depth |
| extension, R = 150 km | growth pool of 188 stations (55 operating on these days); iteration 0 identical to the literal run (3967), iteration 1: 16644. Hyp2000: 69 templates located, depth median 34.6 km, hypoDD kept 13. The XU Big-Skidder array elements (BS01-BS24) dominate the SNR ranking: this is a WATCH item |
| ext score | literal vs R150 table written, with the permanent-station carrier table (e.g. family 1: B943, B013) |

### Round 3 (2026-09-25 afternoon): readings resolved; details in `AB_ROYER_LOG.md`
| # | item | resolution |
|---|---|---|
| C1 | sec.3.1 MAD (NW had -34 %) | The MAD is taken over **each window's own network-autocorrelation record** (one row of the hour matrix; `ab_auto.py --mad-pool row`, now the default), and a pair is culled when both of its windows overlap a stronger pair (2-D pair suppression, `ab_cluster --cull pairnms`, default). Result: **NW 4603 vs 4915 (-6 %), SW 3217 vs 3306 (-3 %)**. The hour-MAD and strict-cull counts remain in the log. The `--cols sample` variant is still hour-MAD and is a reported column only |
| C2 | velocity model | **GSC routine-location onshore 1-D model** (Bostock et al. 2015; Savard 2018 Fig.4.7a, digitised from the vector paths): 0-1 km 5.01, 1-6 6.01, 6-30 6.71, 30-46 7.11, >46 km 7.76 km/s. **Vp/Vs 1.73**. Crosson (1976) and Vp/Vs 1.78/1.85 are columns |
| C3 | Rowe et al. 2002 clustering | "dendrogram-based hierarchical **pair-group** algorithm ... clustering complete when a user-defined CC cutoff threshold is reached" (USGS PP1769 ch.5, using Rowe's code). Headline = **UPGMA (average linkage), cut at the paper's median + 8 MAD significance, min size 2** (A16). Null-1 counts, excess and null/real ratio are columns in `sweep_long.csv` |
| C4 | pre-2009-08 tremor days | own **Wech & Creager (2008)-style envelope cross-correlation detector** (`ab_wecc.py`), calibrated against PNSN day labels 2009-08-06..2011-12-31, then applied to 2007-01..2009-08-05 (NW) and 2007-2008 (SW). The calibration table (accuracy, kappa, TPR, FPR per threshold) is `tremor/wecc_calibration.csv` |
| C5 | NVI added | sec.2.2/3.1. 27-station pool from Fig.1B (23 PO, TASB (C8), MAYB/WOSB/NCRB (CN, EHZ-only, NRCan FDSN)). Autocorrelation, literal "different combinations of 7 of the 10": per hour, the first 7 of the paper's listed 10 that recorded three-component data that hour (2006: VI10 VW03 VI11 VW02 VI52 VI53 VI05; 2007, no VW: VI10 VI11 VI52 VI53 VI05 VI08 VI06). Scan days = the two episodes (sec.2.2: only those two were available). Published: 4267 -> 54 -> 47, depths 33-39 km, sigma 1.8/2.16 km, hypoDD 91 % |
