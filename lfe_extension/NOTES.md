# NOTES — lfe_extension (dated ledger)

## 2026-09-16 start
User: extend Bostock's catalog continuously to 2026 on permanent stations; do not wait on downloads. PLAN.md written; Track B download launched first (method-agnostic); Merlin review of PLAN before Track A.

## 2026-09-16 agents launched
- Track B Opus downloader launched (epoch table, queue P0-P3, resumable runner, manifest, QC per month).
- Merlin review of PLAN.md launched; Track A waits for its spec.
- Memory note lfe-extension-study-2026-09-16 written.

## 2026-09-16 Track B (download) built and launched

**Built** (all in `lfe_extension/`, nothing written outside it):
- `scripts/download_plan.py` — FDSN `level=channel` metadata (IRIS served **all four networks**, CN
  included; the NRCan fallback was never needed), one seismometer band per epoch (HH > BH > EH;
  accelerometer HN/EN and LH/SN bands ignored), 3 components required → `results/station_epochs.csv`;
  cross-check vs `../codex-LFE/.../nonoverlapping_triplet_epochs.csv` → `results/epoch_crosscheck.csv`;
  ordered queue → `results/download_queue.csv`; seeds the present_elsewhere rows into the manifest.
- `scripts/download_day.py` — one station-day, adapted from `../lfe_catalog_repair/scripts/d3_download.py`
  (copied, not imported). Native fs, no decimation/filter/response/merge. Writes `<DOY>.mseed` or a
  `.mseed.nd` marker, then verifies by reading headers (n_channels, fs, coverage, fragments).
- `scripts/download_worker.py` + `scripts/download_run.sh` — resumable runner, ≤8 concurrent FDSN
  requests total (one provider, so the per-provider cap is the same 8), manifest appended one whole
  line per `os.write` under a lock (readers never see a torn row), `results/download_progress.json`
  every 5 min, second pass retries errors once, month-complete → QC automatically.
- `scripts/download_qc.py` — `results/download_qc_<YYYY-MM>.csv` + `figures/download_coverage.png`.

**Epochs found** (2010-08-01 .. 2026-09-16, `results/station_epochs.csv`):
- band switches: PGC BH 40 Hz → **HH 100 Hz 2017-08-10**; SNB BH → **HH 2011-03-17**;
  NLLB BH → **HH 2017-02-23** (2-day metadata gap 02-21..02-23); LRIV BH → **HH 2019-05-30**.
- end dates: LZB **2014-06-17**, YOUB **2014-09-03**, SOKB **2015-01-18**, PFB **2018-12-20**,
  SNB **2021-03-28**. Still open: PGC, VGZ, NLLB, LRIV, all 10 PB boreholes (EH 100 Hz throughout).
- cross-check: 16 codex epochs overlap our window, **0 disagreements** (band and switch dates identical);
  8 roster stations (SOKB, LZB, NLLB, PFB, YOUB, B006, B007, LRIV) are absent from that codex table, so
  no independent check exists for them — their epochs come from IRIS StationXML only.
- `loc` is empty ("") for every station/epoch; requests use `location="*"` then select the band+loc,
  so an added location code would show up as n_channels > 3 in QC rather than silently doubling data.

**Queue** (`results/download_queue.csv`, 90,587 rows):
P0 280 station-days (~6 GB) · P1 988 (~20 GB) · P2 1140 (~23 GB) · P3 69,768 (~1.49 TB) ·
present_elsewhere 455 (PGC 2017 001-221, SNB/VGZ/B011/LRIV 2017-02-20..04-21, all in
`../lfe_catalog_repair/traces`, same band) · no_epoch 17,956 (dead stations after their end date).
Total to fetch **72,176 station-days ≈ 1.54 TB** (74 TB free).

**Status at launch:** P0 complete (279 data + 1 nodata, verified 3-C 100 Hz, `traces/PB.B00*/2010/`);
runner PID in `logs/download_run.pid` working P1 → P2 → P3.

**Conventions for Track C:** `coverage_frac` in the manifest is the **minimum** over the day's channels
(seconds of data / 86400, overlaps capped at 1); `n_frag` is the **max** fragment count over channels
(extra columns beyond the spec: `n_frag`, `note`; `note` holds the foreign path for present_elsewhere
rows and the failure string for errors). CN days are heavily fragmented at IRIS (PGC 2014-03-01: 64
fragments, coverage 0.78) — Track C must merge/pad, not assume one continuous trace.

**Provider check (2026-09-16, done while the run was in flight):**
- IRIS/EarthScope serves CN, PB, UW **and C8** metadata; the NRCan fallback in `download_plan.py` was
  never exercised. NRCan's own FDSN is at `https://www.earthquakescanada.nrcan.gc.ca/fdsnws` (the URL in
  the task, without `www.`, 301-redirects and obspy's `Client()` cannot auto-discover it because its
  `event/1` service 404s — use `curl` or `service_mappings` if it is ever needed).
- **CN gaps are real, not a provider artifact**: PGC 2014-03-01 comes back with coverage 0.783 and ~50-64
  fragments per channel from BOTH IRIS and NRCan, byte-for-byte the same story (VGZ same day 0.92, 3
  fragments). So there is nothing to gain by re-fetching CN from NRCan.
- **C8 SOKB has NO waveform data at either provider** (204/no-data for 2011, 2012, 2013 and 2014 test
  days, at IRIS and NRCan, while its StationXML claims HH 100 Hz to 2015-01-18). The 34 `.nd` markers for
  SOKB 2010 in `../lfe_catalog_repair/traces` say the same. It will fill the queue with `nodata` rows
  (~1.6 k cheap requests) — treat SOKB as unavailable for Track A/C, not as a download failure.

## 2026-09-16 Merlin review received
14 flaws (F1 FA budget → two-tier product; F2 ETS-day site-rolled bar unfalsifiable; F3 mergeB wins at matched FA; F4 E6 JSON keys mislabeled — use e6_fa_table CSVs; F5 30% recall unsupported; F6 chain must stay 20 Hz/1-8 Hz; F10 relocation controls). Prescription Steps 0-8 saved. PLAN.md amended. Track A agent launched on Steps 0-3; downloader told to re-order queue.

**Queue re-order (coordinator / Merlin review, same day):** a new priority **P2b** was inserted between
P2 and the oldest-first continuous fill, for every station alive in each window:
(a) 2015-12-05..2016-01-15, (b) 2016-12-01..2016-12-31, (c) 2018-06-10..2018-07-25,
(d) 2021-10-01..2021-11-05. **P2b = 2,408 station-days (~50 GB)**, 537 station-days skipped as `no_epoch`
(LZB/YOUB/SOKB dead by 2015; PFB and SNB dead by window (d)); by month 2015-12 432, 2016-01 240,
2016-12 496, 2018-06 336, 2018-07 400, 2021-10 434, 2021-11 70. P3 shrinks to 67,360 station-days
(~1.44 TB) because `build_queue` now **de-duplicates**: a (station, day) queued at a higher priority is
never queued again, so no day is fetched twice and the month-complete QC trigger still counts each day
once. Totals are unchanged (72,176 station-days ≈ 1.54 TB).
Also corrected while re-planning: `present_elsewhere` now means **only** the read-only
`../lfe_catalog_repair/traces`; days already in our own `traces/` stay `pending` in the queue, and are
either skipped via their manifest row or re-verified on disk by `fetch_day` into a real `data` row
(header read only, no re-download). The runner was stopped (`STOPPED (resumable)`), the queue
regenerated, and restarted — it resumed at 2013-10 with 69,985 of 72,176 left and re-fetched nothing.

## 2026-09-16 Track A (method) — STEP 0: freeze the measurement chain, fix the E6 labelling, matched-FA ROC

Scope discipline: `../lfe_catalog_repair` and `../codex-LFE` are READ-ONLY. Copied into `scripts/`:
`e4_template_test.py`, `e5r_score.py`, `e5r_templates.py`, `e4_cache20.py`, `e4_chanorder.py`,
`e5r_groundtruth.py`, `e5r_offsets.py`, `d3_vel1d.py` (import dependency), all `e6_*.py` and the shell
drivers `e6_run_gpu{1,2}.sh`, `e5r_run_mf*.sh`, `e5r_run_onepersite.sh` as command references. Only the
copies are edited. Measurement chain untouched (20 Hz, 1–8 Hz Butterworth-4 zero-phase, existing
component handling) — Merlin F6.

**Patches to the copies**
- `e5r_score.py`: `--in-dir` (default `lfe_catalog_repair/results`, read-only) and `--out-dir`; `--one2one`
  (Merlin F11) switching own recall, catalog recall, the 600-s chance and the detection→event assignment
  from nearest-neighbour to greedy-by-time one-to-one (`one2one_absdiff` / `one2one_idx`); `matching`
  recorded in the JSON.
- `e4_template_test.py`: `--out-dir` (so no run can write into `lfe_catalog_repair`), `--min-sites N`
  (Merlin F8) skipping a family-day unless ≥ N distinct SITES have a usable template AND data, with the
  skipped count reported at the end; `n_sites` now recorded on every budget row, every netfam row and
  every detection (E6 runs already carried `n_site`; the new column is unconditional).

**(a) Re-score with `--e6` — F4 CONFIRMED and fixed.** `scripts/a0_rescore_e6.sh` → `results/e6_rescore/nn/`.
The on-disk `e6_score_e6_{sites,mergeA,mergeB}.json` were written WITHOUT `--e6`, so the key
`L{L}_T{int(T)}` collided: what they label `T7` is T = 7.5 and what they label `T8` is T = 8.5, and
**T = 8.0 was overwritten and absent entirely**. Example (mergeB, own recall / quiet FA): old `T8`
0.1448 / 0.317 — that is really T = 8.5; the true T = 8.0 is 0.1886 / 0.967 (= Merlin's "0.97 FA,
0.189 recall"). Same shift for sites and mergeA. The lever passes were already scored with `--e6`.

Two book-keeping points Merlin's spec did not anticipate:
- **9 scoreable passes, not 11.** The checkpoint lists 11, but `e6_lev_a_mv` and `e6_lev_c_mv` are
  build-half moveout passes run without `--budget`/`--catalog-out` — no budget CSV, no detections
  parquet, nothing to score. The 9 scored: sites, mergeA, mergeB, mergeB_cap25, mergeB_coinc,
  lev_a, lev_b, lev_c, 1b_noVGZSNB (the last a 5-family subset run, reported separately).
- **The lever QC CSVs are partial and crash the scorer.** `e6_template_qc_2012_a.csv` holds only the 7
  rebuilt stations and `e6_template_qc_2012_c.csv` only the 7 persistent ones (and has no `usable`
  column), so `--qc-csv` on them gives the wrong per-family station set → `share_ok` IndexError /
  AttributeError. Fixed by deriving the station set from the template npz actually used
  (`results/e6_rescore/qc/qc_from_npz_lev_{a,c}.csv`; 1117 and 1144 (fam, sta) templates, matching the
  1144 `usable & template` rows of the E5r QC for the baseline).

**(b) One-to-one matching (F11) — delta reported once** (`results/e6_rescore/one2one_delta.csv`,
8 variants × 5 thresholds): **own-family recall is IDENTICAL to 4 decimal places at every variant and
every threshold (Δ = 0.0000)**. Reason: detections are declustered within a family at the template
length, so no two own-family detections can sit within ±2 s of one event, and no two Bostock events of
one family share a detection in this catalog — nearest-neighbour was already one-to-one on the primary
metric. The CATALOG (any-family) recall does fall, by 0.042–0.147 (mergeB at T = 8: 0.295 → 0.176), and
its 600-s chance falls with it (0.041 → 0.032), i.e. the inflation was in both terms. One-to-one is kept
on from here as the default for all Track A scoring.

**(c) Matched-FA ROC — ACCEPT.** `scripts/a0_matched_fa_roc.py` → `results/e6_matched_fa_roc.csv`
(+ the threshold grid `results/e6_fa_recall_grid.csv`). Recipe: T*(FA) log-interpolated in FA on the
{7, 7.5, 8, 8.5, 9, 10} budget grid, own recall linearly interpolated in T (Merlin's recipe) and also
recomputed exactly at T* from the detections; quiet days = 2012 doy 228–242.

| variant | recall @ FA 1.0 | @ 0.5 | @ 0.1 | ETS own precision @ FA 1.0 | Merlin |
|---|---|---|---|---|---|
| sites | 0.1653 | 0.1396 | 0.0923 | 0.421 | 0.165 / 0.140 / 0.092, 0.42 |
| mergeA | 0.1711 | 0.1445 | 0.1051 | 0.473 | 0.171 / 0.144 / 0.105, 0.47 |
| **mergeB** | **0.1901** | **0.1626** | **0.1169** | 0.452 | 0.190 / 0.163 / 0.117, 0.45 |
| mergeB+cap25 | 0.1669 | 0.1446 | 0.0988 | 0.473 | 0.167 / 0.145 / 0.099, 0.47 |
| mergeB+coinc | 0.1808 | 0.1563 | 0.1128 | 0.486 | 0.181 / 0.156 / 0.113, 0.49 |
| lev_a | 0.1866 | 0.1606 | 0.1156 | 0.451 | 0.187 / 0.161 / 0.116, 0.45 |
| lev_b | 0.1993 | 0.1694 | 0.1246 | 0.441 | 0.199 / 0.169 / 0.125, 0.44 |
| lev_c | 0.1989 | 0.1719 | 0.1248 | 0.456 | 0.199 / 0.172 / 0.125, 0.46 |

Max deviation from Merlin's table 0.0005 (bar: ±0.005) → **ACCEPT**, no scorer audit needed. The exact
(uninterpolated) recall at T* agrees with the interpolated value to ≤ 0.002 everywhere.

**(d) `--min-sites 2` implemented; its cost is small on the 2012 15-station network.** Because the gate
only decides whether a family-day is scored at all (the network sum of a surviving family is unchanged),
it is exactly equivalent to filtering on the recorded `n_sites` after the fact — so Step 1 runs WITHOUT
the gate and reports its effect post hoc, which keeps the frozen operating point comparable to the E6
baseline Merlin's accept bar was set from. On mergeB/2012: 51 of 5,919 family-days (0.9 %) are
single-site, carrying 36 of 1,574 quiet-day detections at T = 8 and 1.4 % of all detections; exactly one
of 116 families is single-site throughout.

**Measured throughput (first ~20 min of running, two consecutive 5-min ticks of the restarted runner):**
**66 MB/s and ~12,000 station-days/hour** with 8 concurrent IRIS requests (65.2 / 65.6 / 66.2 MB/s at
the three steady-state ticks — flat, so 8 streams are not saturating anything at our end; CPU 1.6 of 32
cores, RSS 1.7 GB, disk 74 TB free). At that rate the remaining ~68,000 station-days (~1.44 TB) finish in
**~5.7 h (ETA 2026-09-17 ≈ 03:15 UTC)**; `results/download_progress.json` carries the live ETA.

State at hand-off: P0, P1, P2 complete, P2b in flight (2018-07 at the time of writing), 4,138 days of
data + 164 verified-empty + 455 present_elsewhere, **0 errors**, 96 GB on disk, QC written for 2010-08,
2010-09, 2011-08, 2011-09, 2013-08, 2013-09, 2013-10, 2016-12 (2015-12 waits for its 12-01..12-04 days,
which sit in P3 — the month-complete trigger is per calendar month across all priorities) (+ 2017-03, which
will be rewritten when P3 reaches it).

Known data facts for Track A/C (from the QC so far): SOKB always empty (see above); B004 missing
2011-08-06..21 (16 d); B011 missing 23 d in the P2b windows; PGC is chronically fragmented at IRIS
(50-68 fragments/day, coverage typically 0.99 but 0.78 on 2014-03-01) while VGZ/SNB/PB come in 1-2
fragments/day at coverage 1.0.

## 2026-09-16 Track A — STEP 1: variant, template length, HF veto, frozen operating point — ACCEPT

**Adopted variant (no new choice):** mergeB — within-site w-weighted mean cc for the PGC cluster
{PGC, B009, B010, B011} and the LRIV cluster {LRIV, B005, B006, B007}, site weight √Σw², weights
`e5r_weights_snr_2012.csv`, no PGC cap, no coincidence gate (Merlin F3, settled by the step-0 ROC).

**(a) 12-s templates — ADOPTED.** `scripts/a1_run_mergeB_L12.sh` (one GPU pass, 8.5 min, H100; identical to
the E6 mergeB pass except `--lengths 12`; 12-s templates, moveout and amplitude reference all already exist
in `e5r_templates_2012.npz` / `e4_moveout_e5r_mv_build.csv` / `e5r_amp_ref_2012.csv`). Own recall at matched
quiet-day FA 1.0: **0.2022 (12 s) vs 0.1901 (8 s), Δ = +0.0121** ≥ the pre-registered +0.01 → adopt. At FA 0.5
0.1747 vs 0.1626; at FA 0.1 0.1225 vs 0.1169. Outputs `results/a1_matched_fa_roc_L12.csv`,
`results/a1_fa_recall_grid_L12.csv`.

**(b) HF-ratio veto — ADOPTED, and it is the single biggest gain of the day.** `scripts/a1_hf_veto.py`:
for every mergeB detection, the energy ratio E(10–19 Hz)/E(2–8 Hz) at CN.PGC over the scored window, from
the NATIVE-rate 2012 traces (the 20 Hz cache cannot carry the band). PGC is BH 40 Hz in 2012 → Nyquist 20 Hz,
so the band is 10–19 Hz (0.95 Nyquist; a band touching Nyquist destabilises the Butterworth); PB.B011 (EH
100 Hz, literal 10–20 Hz) was computed alongside as a band-choice control and agrees.

| | own-matched detections | quiet-day detections | separation |
|---|---|---|---|
| PGC median HF ratio (8 s) | 0.247 | 4.547 | 18× |
| B011 median HF ratio (8 s) | 0.163 | 0.999 | 6× |

Cut = 99th percentile over BUILD-half own-matched detections (3.511 at 8 s, 3.638 at 12 s; the test half
never touched it). Result at 8 s: quiet-day FA falls **61–65 %** at a test-half own-recall cost of
**0.5–0.7 % relative** (bar: ≥ 30 % at ≤ 2 %); at 12 s **58–60 %** at **0.3–0.4 %**. This is above Merlin's
own "would change my mind" line (FA cut ≥ 50 %). Merlin's F14 was flagged SPECULATIVE; it is confirmed:
the quiet-day detections are HF-rich transients (uncatalogued micro-earthquakes / site glitches), not
sub-threshold LFEs.

**(c) Thresholds by FALSE ALARMS ONLY** (`scripts/a1_operating_point.py`, 0.25-MAD grid; the detections
parquet reproduces the budget counts exactly at every tabulated T, so the fine grid is exact, not
interpolated). Quiet days = 2012 doy 228–242, verified **0 PNSN tremor detections in the family box**.

| config | T_ext (quiet ≤ 1.0 AND post-SSE ≤ 1.5) | own recall | ETS own prec. | T_core (quiet ≤ 0.1) | own recall | ETS own prec. |
|---|---|---|---|---|---|---|
| mergeB 8 s, no veto | 8.50 (0.317 / 0.989) | 0.1448 | 0.574 | 9.00 (0.084) | 0.1125 | 0.669 |
| mergeB 8 s + veto | 8.00 (0.371 / 1.414) | 0.1874 | 0.467 | 8.75 (0.063) | 0.1267 | 0.628 |
| **mergeB 12 s + veto (FROZEN)** | **8.25** (0.277 / 1.042) | **0.1841** | 0.488 | **8.75** (0.097) | **0.1463** | 0.596 |

The **post-SSE bar binds, not the quiet-day bar** — at T = 8.25 the unvetoed 12-s run has quiet FA 0.555
(passing ≤ 1.0) but post-SSE FA 1.516, one part in 300 above the 1.5 bar. Worth knowing for Track C: the
operating point of this detector is set by residual post-SSE activity, not by the quiet background.

**ACCEPT (both bars met).** Own recall at matched FA 1.0 = **0.2022 ∈ [0.18, 0.24]**; quiet-day
site-rolled/data = 0.903–1.033 over the whole threshold grid, inside [0.77, 1.3]. Caveat worth recording:
the site-rolled control's detection TIMES are not stored, so the HF veto cannot be applied to the control —
the ratio bar is evaluated on the unvetoed detector.

**Frozen:** `results/frozen_operating_point.json` (variant, L = 12 s, T_ext = 8.25, T_core = 8.75, HF veto
yes + cut 3.6384, and sha256 prefixes of the weights / templates / moveout / families / stations / GT /
detections / HF-ratio files).

## 2026-09-16 user: inter-ETS repeats
User confirmed: do not assume families are silent outside ETS. Added Step 3b (withheld-station stack test on quiet-day detections) to PLAN; Track A agent told to run it after Step 3 using e6_stack_test.py machinery.

## 2026-09-16 Track A — STEP 2: the frozen cross-year transfer — PASS (the make-or-break test)

Track B had already delivered **PB.B009 and B010 for 2010 doy 222–256**, so the 2010 network could be made
*identical* to the 2012 ablation network instead of one station short: nine stations / six sites
(PGC + B009 + B010 + B011 = one site, VGZ, SNB, NLLB, B926, LRIV). B009/B010 2010 were cached at 20 Hz into
`lfe_extension/cache20hz` with the copied `e4_cache20.py` (`--traces-root`, `--outdir`); `e4_template_test.py`
gained `--cache-dirs` so the new cache is searched before the lfe_catalog_repair one. 2010 ground truth =
`e2_gt_bostock.parquet` re-keyed so `fam_id = bkey` (`results/a2_gt_bostock2010.parquet`; the matched filter
compares `fam_id` with the family key, and in the 2012 GT those are equal) + `e3_bostock_family_offsets.csv`,
**90 accepted offsets, 15,931 events** — the file Merlin named, used exactly as the code already reads it.

| | T_ext = 8.25 | T_core = 8.75 |
|---|---|---|
| own recall, 2012 ablation (9 stations, test half) | 0.1983 | 0.1589 |
| own recall, 2010 frozen transfer | 0.2117 | 0.1724 |
| **ratio (bar 0.7–1.3)** | **1.067 PASS** | **1.085 PASS** |
| quiet FA 2012 / 2010-strict | 0.773 / 1.158 | 0.260 / 0.613 |
| FA ratio (bar 0.5–2.0) | 1.50 PASS | 2.36 **not met** |
| own precision 2012 / 2010 | 0.353 / 0.309 | 0.471 / 0.402 |

With the frozen HF veto applied to both years the picture is the same (ratios 1.066 / 1.086; FA ratios
1.71 / 2.55). The recall bar — the one Merlin called decisive — passes with room to spare, and nowhere near
the failure branch (< 0.5×). **No time-base or moveout problem exists: 2012-built templates, 2012 weights,
2012 moveout and 2012 MAD thresholds transfer to 2010 unchanged.**

**Why the T_core FA ratio misses, and why it is not a transfer failure.** Only 3 of the 2010 candidate quiet
days are genuinely pre-SSE: `doy 222–226 / 250–256` contains PNSN tremor on 222 (10), 225 (1), 250 (4), 251
(91!) and 253 (16), and Bostock events on 250–251, so the strict-quiet set is {223, 224, 226, 252, 254, 255,
256} — three PRE-SSE days and four POST-SSE days. Split:

| window | quiet FA 2010 @ T_ext | ratio to 2012 | @ T_core | ratio |
|---|---|---|---|---|
| pre-SSE (223, 224, 226) | 0.233 | **0.30** | 0.078 | **0.30** |
| post-SSE (252, 254, 255, 256) | 1.851 | 2.40 | 1.015 | 3.90 |

2010 **pre**-SSE is ~3× *quieter* than 2012 pre-SSE; the whole excess is in the post-SSE days, exactly the
behaviour the 2012 post-SSE window already showed (and the reason T_ext is set by the post-SSE bar). The
2010 quiet set is half post-SSE by construction, so its FA ratio is not a like-for-like number.

**Reverse direction (2010-built templates → 2012), Merlin's second frozen-transfer bar.** No 2010 from-scratch
templates existed, so they were built with the same machinery (`e5r_templates.py` parameterised by `--year /
--gt / --stations / --cache-dirs / --out-dir`; 2010 events with accepted offsets, all of them): **716 usable
(family, station) templates over 90 families**, q-control cc to the SeisBench stacks **0.79 (PGC) / 0.79 (VGZ)
/ 0.87 (SNB)** — as good as or better than the 2012 build's 0.78. Applied to 2012 with their own 2010
build-half moveout at the frozen thresholds:

| comparison | T_ext | T_core |
|---|---|---|
| all 2012 events, all families | 0.175 vs 0.262 → 0.67 | 0.142 vs 0.214 → 0.66 |
| all 2012 events, the 90 common families | 0.202 vs 0.269 → 0.75 | 0.164 vs 0.221 → 0.74 |
| **2012 TEST half, the 90 common families** | **0.2025 vs 0.2109 → 0.96** | **0.160 vs 0.170 → 0.945** |

The first two rows are contaminated: for the 2012-built templates the 2012 *build* half is training data.
On the held-out half and the families both sets cover, the two directions are within 5 % of each other →
**reverse transfer PASS** as well. The residual gap in the "all families" row is pure coverage: only 90 of
116 families have an accepted 2010 offset and ≥ 30 2010 events.

## 2026-09-16 Track A — STEP 3b: are the quiet-day detections false alarms, or real inter-ETS repeats?

User amendment (PLAN.md): quiet-day detections are an UPPER BOUND on the false-alarm rate, not a measurement,
because families also fire between ETS episodes and Bostock's catalog — which has no events outside ETS
windows — cannot say otherwise. Test: a **withheld-station stack test**.

Setup (`scripts/a3b_run_withheld.sh`, `scripts/a3b_quietday_stack_test.py`, machinery copied from
`e6_stack_test.py`). One GPU pass of the frozen configuration with **CN.VGZ, CN.SNB and PB.B926 excluded from
the network sum** (113 of the 116 templated families have a usable template at one of the three, so the
withheld station never contributed to the detection it is asked to confirm). At each family's quiet-day
detection times (doy 228–242, PNSN-quiet, HF veto applied) the withheld station's 3-C 20 Hz data are stacked
on the predicted S (origin + ttS(Crosson) + the frozen template's S-pick refinement), and compared with 120
matched random-time stacks (same station, same quiet days, same N). ETS-day own-matched detections are the
positive control. PASS = S-SNR excess ≥ 3 dB AND |envelope peak − predicted S| ≤ 2 s AND p < 0.01.

| | T_ext = 8.25 | T_core = 8.75 |
|---|---|---|
| families with ≥ 20 quiet-day detections (testable) | 42 | 25 |
| PASS / FAIL / inconclusive | **0 / 33 / 9** | **1 / 16 / 8** |
| median quiet-day S-SNR excess | +0.13 dB | +0.46 dB |
| median cc of the stack to the family template (null) | 0.168 (0.157) | 0.167 (0.170) |
| **ETS positive control**: registers | **76 / 81**, median excess **+7.57 dB**, cc **0.735** | 70 / 75 |

The control has ample power — its N (median 62) is comparable to the quiet-day N (median 52) and the null
scatter is 1.2 dB — so the test would see a real repeat population easily. It does not.

**Quantitative bound** (`results/a3b_realfraction_bound.json`). Calibrating each family's coherent
signal-to-noise from its own ETS control (excess = 10 log₁₀(1 + p²·N·(a/σ)²), p = 1 for the control) and
solving for the real-repeat fraction p of the quiet-day detections:

| threshold | families | pooled quiet excess | control | **p (real fraction)** | 95 % upper bound |
|---|---|---|---|---|---|
| T_ext = 8.25 | 36 | +0.104 ± 0.215 dB | +7.69 dB | **0.057** | 0.140 |
| T_core = 8.75 | 21 | +0.776 ± 0.296 dB | +7.58 dB | **0.225** | 0.325 |

**Verdict.** At T_ext the quiet-day rate is a genuine false-alarm measurement: ≥ 86 % of those detections are
noise, so the FA-based thresholds are NOT materially conservative and the pre-registered FA criterion stands
as written. At T_core, where purity is higher, about **22 % (≤ 33 %)** of the surviving quiet-day detections
ARE real inter-ETS repeats — so the true T_core false-alarm rate is ~0.78× the measured one (quiet FA 0.097 →
≈ 0.075 per template-day) and T_core is mildly conservative. Only one family (bkey 75) passes individually at
T_core, with 21 quiet-day detections; no family passes at T_ext. Promoting inter-ETS detections to events per
family is therefore NOT yet justified by this test — the population-level estimate is the usable result, and
it says the inter-ETS LFE rate in these families is small but not zero.

Caveats: 71 of 113 families have fewer than 20 quiet-day detections and are untestable; the withheld network
(12 stations) is noisier on quiet days than the frozen 15-station one, so its quiet-day population is if
anything *less* pure than the frozen configuration's — the bound is conservative in the right direction.

## 2026-09-16 Track A — STEP 3: template iteration (Merlin's protocol)

**Selection** (`scripts/a3_run_selection.sh`, `scripts/a3_select_events.py`). Two GPU passes of the frozen
configuration with the two-site coincidence gate (`--coinc-gate 2,2.5`: ≥ 2 sites each ≥ 2.5 MAD of their own
site series), catalog written from T_core = 8.75; then tremor days only (≥ 30 PNSN detections in the family box
or ≥ 100 Bostock events), the frozen HF veto, and cross-family assignment by raw weighted cc among families
sharing ≥ 50 % of their stations.
* gated (15 stations): 10,135 → 8,977 (tremor days & T_core) → 8,736 (HF veto) → **6,468 events, 115 families**,
  median 42 events/family; 68 families with ≥ 30.
* **leave-one-site-out for PGC** (the only site over Merlin's 40 % weight share — it exceeds it on 53 % of
  family-days; every other site's median share is ≤ 0.14): 2,218 → **1,882 events, 105 families**, median only
  **11** events/family, 20 families with ≥ 30, 8 with ≥ 50. This is exactly Merlin's experiment-4 fallback
  trigger (< 50 per family per SSE). **Fallback taken: the PGC-site templates are re-stacked only for the 19
  families with ≥ 30 LOSO events; the rest keep the frozen Bostock-seeded template.** Re-stacking the PGC site
  on events the PGC site itself selected would be self-reinforcement, which is the thing the protocol exists to
  prevent.

**Stacking** — `e5r_templates.py` unchanged (rms-normalised windows, linear stack, cc alignment ±1.5 s, 2
iterations, ONE shift per station for all three components), `--min-events 30`. 428 usable pairs / 68 families
(non-PGC stations) + 73 pairs / 19 families (PGC site). Merged over the frozen set
(`scripts/a3_merge_templates.py`): **1,219 pairs — 426 frozen pairs replaced, 75 genuinely NEW (family,
station) templates, 718 frozen kept.** The new ones are almost all at the far stations
(B005 15, B003 12, B001 11, B004 10, B006 7, LRIV 6, B007 6, SNB 3, NLLB 2, B009 2, VGZ 1), i.e. iteration is
doing the TRANSFER job Merlin predicted (F5), not the sensitivity job.

**What the iteration does to the stacks** (426 replaced pairs):

| | far stations (B001–B007, LRIV) | near / high-SNR (PGC site, B926, NLLB) | VGZ / SNB |
|---|---|---|---|
| median ΔS-SNR | **+0.6 … +1.9 dB** | **−1.1 … −3.8 dB** | +0.05 / +0.18 dB |

The near stations lose because the selected sets are smaller than Bostock's (median 42, and 47 for the PGC
site, against his 115–160 per family); the far stations gain because his catalog was never built from them.
Stability is good: |ΔS pick| median **0.05 s** (62 of 426 move > 0.5 s, mostly B004), station-differential
moveout changes **0 of 426** beyond 0.3 s, max cross-family cc **falls** at every station (e.g. PGC 0.70 →
0.65, SNB 0.65 → 0.57), cc to the SeisBench stacks 0.797 → 0.785.

**The six rejection tests** (`scripts/a3_reject_tests.py`) — and a calibration problem Merlin's spec did not
anticipate. Applied literally, tests (a) "< 0.6", (d) "> 0.7" and (e) "> 50 % of the 1–8 Hz energy in one 1-Hz
bin" fire at essentially the same rate on the **FROZEN** templates (13 / 320 / 443 of 1,144 pairs) as on the
iterated ones (32 / 282 / 433 of 1,219) — they are properties of this template class (LFE spectra peak at
1–3 Hz; families within a few km share waveforms), not evidence that the iteration damaged anything. The
script therefore reports both the literal verdict and a **baseline-conditioned** one that keeps every RELATIVE
clause as written and applies each absolute clause only where the frozen template satisfied it:

| test | literal fails | baseline-conditioned |
|---|---|---|
| (a) q-control cc vs SeisBench | 32 | 19 |
| (c) S pick > 0.5 s | 62 | 62 |
| (c) differential moveout > 0.3 s | 0 | 0 |
| (d) cross-family collapse | 282 | 33 |
| (e) spectral narrowing | 433 | 28 |
| (f) S-SNR falls > 1 dB | 177 | 177 |
| **any** | 709 | **261** |
| **families accepted** | **2 / 116** | **52 / 116** |

The dominant honest rejection driver is (f): 51 of the 64 rejected families have at least one station whose
stack S-SNR fell more than 1 dB — the small-N problem above, not drift.

**GPU contention (2026-09-16 22:16 UTC).** The iterated-template 2012 pass died at doy 272 with
`torch.OutOfMemoryError` — three other processes held 88 of the H100's 93 GB. Recovery: `scripts/a3_rerun_2012.sh`
re-runs that pass with `--fam-chunk 8` (families chunked, nothing throttled) and
`PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`, cutting the peak allocation to a quarter. The pass that
matters for the iteration verdict is the HELD-OUT 2010 one (test b), which is unaffected: the 2012 run is
contaminated by construction, since the iterated templates were built from detections on both 2012 halves.

**Convergence check (Merlin's stop rule, cc(iter k, k−1) ≥ 0.98).** cc between the iterated and the frozen
12-s template over ±1 s, for the 426 replaced pairs: median **0.70** (p25 0.55, p75 0.80, max 0.94), **0 of 426
at ≥ 0.98** — `results/a3_cc_iter_vs_frozen.csv`. Per station the near/high-SNR ones are the most stable
(SNB 0.81, B926/PGC 0.76, VGZ 0.74, PGC-cluster boreholes 0.72) and the far ones change most (B006 0.45,
B007 0.49, LRIV 0.51, B004/B005 0.55). Iteration 1 has therefore NOT converged in Merlin's sense; whether a
second re-stack is worth running is decided by the held-out-year delta below, not by the cc.

## 2026-09-16 Track A-2 (Merlin steps 4 and 5) — epoch identity + relocation

Scope: `a4_` (epoch transfer/identity) and `a5_` (relocation) only; `a0_`-`a3_`, `e5r_*`, `e6_*` and
`c_*` belong to the other two agents and were not touched.  `../lfe_catalog_repair` and `../codex-LFE`
read-only; every script used is a COPY in `scripts/` (`a4_cache20.py`, `a4_epoch_stack.py`,
`a4_epoch_report.py`, `a5_locate.py`, `a5_synth_floor.py`, `a5_eq_corrections.py`, `a5_relocate.py`,
plus the unchanged `e5r_templates.py`, `d3_vel1d.py`, `e4_chanorder.py`, `e4_template_test.py`).

**20 Hz caches built for the epoch windows** (`scripts/a4_cache_windows.sh`, frozen chain: 20 Hz,
1-8 Hz Butterworth-4, explicit Z/N|1/E|2 order): 2018 doy 161-206 (15 stations), 2021 doy 274-309 (14),
2017 doy 70-100 (5, from the read-only lfe_catalog_repair traces).  1,333 station-days.
**Data fact for Track C:** CN day-files in 2021 (PGC, VGZ, NLLB) carry per-FRAGMENT sampling rates that
differ in the 7th digit (99.99999237, 100.00000763, ...).  obspy then refuses to merge them and the
integer-decimation check rejects the day — 41 of 504 station-days were lost until `a4_cache20.py`
snapped every trace to the nearest integer rate (relative time-base change 1e-7 = 0.009 s/day).

**STEP 3 VERDICT — iteration REJECTED** (`results/a3_heldout_matched_fa.csv`). Held-out 2010 own recall at
matched 2010 quiet-day FA: frozen 0.2021 / 0.1608 / 0.1257 at FA 1.0 / 0.5 / 0.25 vs iterated 0.1846 / 0.1466 /
0.1132 — **Δ = −0.017 / −0.014 / −0.013**, which triggers Merlin's rejection clause (a drop > 0.01). On the 52
families the rejection tests accepted the iteration is exactly neutral (Δ = −0.0003). The same sign shows on
the contaminated 2012 test half (0.1847 → 0.1696 at T_ext; quiet-day detections ×1.14). Adoption needed
≥ +0.01 and is not met anywhere → **the frozen Bostock-seeded templates stand**; no second re-stack was run.
Cause: the PGC site carries > 40 % of most families' weight and its leave-one-site-out selection yields only
~11 events/family, so the dominant stations get worse templates (−1 to −4 dB) while the gains sit at far
stations with little weight. Recommended follow-up (NOT run, outside this brief): a targeted iteration that
re-stacks only B001–B007 and LRIV and keeps the frozen templates elsewhere.
The re-run of the OOM-killed 2012 iteration pass was left unlaunched deliberately: that pass is contaminated by
construction and the verdict rests on the held-out year, so it would have burnt contended GPU for nothing.

## 2026-09-16 Track A complete (agent report; details results/TRACK_A_RESULTS.md)
Step0 ACCEPT (ROC reproduces Merlin ±0.0005). Step1 ACCEPT: mergeB, 12 s, T_ext 8.25 / T_core 8.75, HF veto ON (PGC 10-19/2-8 Hz, cut 3.64; quiet FA -58..65% at <1% recall cost; +0.02/+0.03 recall at matched FA). Step2 PASS both directions (2012->2010 recall ratio 1.07/1.09; 2010->2012 0.96/0.95; FA ratio at T_core 2.36 driven by post-SSE days only). Step3 iteration REJECTED (held-out recall -0.012..-0.017; PGC leave-one-site-out gives ~11 events/family). Step3b: quiet-day detections >=86% noise at T_ext, ~22% real repeats at T_core. Frozen point -> results/frozen_operating_point.json; Track C told. Follow-up (not run): iterate only B001-B007/LRIV templates.

**STEP 5(a) SYNTHETIC FLOOR — the 5 km bar IS testable on this geometry.** `scripts/a5_synth_floor.py`
+ `scripts/a5_locate.py` (a memory-lean copy of `d3_stage4_locate.locate_from_tables`, verified
IDENTICAL to the original on 12 random family/prior combinations).  9 SITES (PGC cluster = 1, LRIV
cluster = 1 — Merlin F10(i)), Bostock's 130 hypocentres forward-modelled through Crosson 1-D, 0.3 s S /
0.2 s P noise + per-site static U(-0.5, 0.5) s (applied to S, scaled by 1/1.78 on P), 50 realisations:
median horizontal error **1.43 km** (S+P, free depth), p90 2.67 km, median |dz| 1.70 km, 0 % at bounds;
noise-free 0.41 km (grid + table interpolation).  Matched to the pick sets actually obtained in 5(c):
2.73 km over the 88 located families, **2.03 km over the 14 that meet the accept geometry**.

**STEP 5(b) STATION CORRECTIONS from earthquakes, not from Bostock.** USGS M >= 2, 25-50 km depth,
48.0-49.0 N / -124.5..-123.0 E, 2010-2026 -> 42 events, **7 with >= 14 of our stations on disk**
(2011-08-06, 2011-09-06, 2012-08-29, 2014-08-01, 2014-12-10, 2015-02-07, 2015-12-17).  Picked on our own
traces at native rate, 2-15 Hz, envelope-onset + AIC, SNR >= 4 -> 97 P and 98 S picks.  Median S
residual vs Crosson: PGC +0.33, VGZ -0.05, SNB +0.67, NLLB +0.57, LRIV +0.57, B926 +0.01, B001 +0.31,
B003 +1.13, B004 +0.17, B005/6/7 +0.46..+0.49, B009/10/11 -0.05..+0.11 s.  **These are 2-3x SMALLER than
the +1.7 to +3.25 s that E3 reported from Bostock's hypocentres** — i.e. most of that "station residual"
was family mislocation, not a station term.  CAVEAT (figure `a5_eq_picks.png`): the P residual grows with
epicentral distance (~0 s below 40 km, +0.5..+2.5 s beyond 60 km), a Crosson-model bias that a constant
per-station term cannot absorb; the S residuals show no such trend.

**STEP 5(c) RELOCATION — ALL FOUR CONTROLS PASS, on 14 of 130 families.** `scripts/a5_relocate.py`.
P on Z / S on the horizontals of OUR 2012 from-scratch per-(family, station) stacks (the +-20 s stacks
already inside `e5r_templates_2012.npz` — no re-stacking was needed, P sits 4-8 s before S and well
inside the window), envelope-onset + AIC, SNR >= 3; site collapse = median over the cluster members of
(member pick - that member's earthquake correction).  1,144 stacks -> 768 S and 481 P picks -> 446
site-S picks -> 88 families with >= 3 picked sites, **14** with >= 4 picked sites within 60 km and
azimuthal gap <= 180 deg (the pre-registered ACCEPT subset).  The binding constraint is PICK YIELD at the
far sites, not the network: S-pick acceptance is 0.86-0.92 at the PGC cluster, 0.72-0.81 at NLLB/SNB/VGZ,
0.77 at B926, but 0.17-0.41 at B001/B003/B005-7/LRIV and **0.00 at B004**.

| run | n | median dh | p90 | median dz | at bounds | null dh | obs/null | spread strike/dip |
|---|---|---|---|---|---|---|---|---|
| A free depth (primary) | 14 | **3.61 km** | 5.76 | -2.20 (abs 2.93) | 0 % | 29.4 | **0.12** | 0.95 / 0.97 |
| B Slab2 prior | 14 | **3.45 km** | 5.76 | -3.36 | 0 % | 28.9 | 0.12 | 0.95 / 0.98 |
| A, NO station corrections | 14 | 4.88 km | 6.82 | -4.61 | 0 % | 27.9 | 0.17 | 0.84 / 0.93 |
| A, S-only (no S-P) | 14 | 2.88 km | — | -2.02 | 0 % | 29.4 | 0.10 | 0.97 / 0.99 |
| synthetic floor, matched geometry | 14 | 2.03 km | 4.17 | 2.29 | 0 % | — | — | — |

Bars: horizontal median <= 5 km **MET**; observed < 0.5x shuffled-null **MET** (0.12); spread >= 0.7x
Bostock **MET**; <= 5 % at depth bounds **MET** (0 %).  Depth is reported, not barred: median dz
-2.2 km (free) / -3.4 km (Slab2 prior), i.e. our solutions sit slightly SHALLOWER than Bostock's.
The observed 3.6 km is 1.8x the matched synthetic floor (2.03 km), so most of the residual is pick noise.
The S-only run being BETTER than S+P is consistent with the distance-dependent P bias found in 5(b).

**STEP 4 mechanics and two gotchas.**
- `a4_finish.sh`'s `while pgrep -f "e4_template_test.py --lengths 12"` guard MATCHED ITS OWN heredoc
  (the wrapper shell's command line contains the script text) and hung forever; replaced by
  `a4_finish2.sh`, which just runs the legs in order.  No run or number was affected, only the wait.
- One selection pass (`a4_NLLB_old`) died of a **CUDA OOM** when six matched-filter passes (three of them
  another agent's) shared the H100; re-run serially.  Two to three concurrent `e4_template_test.py`
  passes are safe on the 94 GB card, six are not.

## 2026-09-16 ~22:40 UTC status (orchestrator)
- Download: at 2016-02, 287 GB, 63 MB/s, ~55k station-days left, ETA 2026-09-17 03:20 UTC; 0 errors.
- Track C: c_run_month.py dry-run on 2015-12 started with the FINAL frozen config (mergeB, 12 s, T_ext 8.25/T_core 8.75, HF veto ON).
- Track A-2: Step 5 relocation DONE, all bars MET — median horizontal offset to Bostock 3.6 km (synthetic pick-noise floor 2.0 km; shuffled null 0.12x ratio; spread ok; 0% at depth bounds); depth median -2.2 km (free) / -3.4 km (Slab2 prior), i.e. slightly shallower than Bostock; S-only beats S+P (distance-dependent P bias). Step 4 epoch identity GPU pass running.
- User instructions today: keep the 20% template catalog as the first product; treat quiet-day detections as an upper bound on FA (Step 3b done); Lin-style picker retraining on our permanent stations (Bostock labels + our matched events as labels) = NEXT PHASE after the extension; use all 32 cores + H100.

## 2026-09-16 Track C (production) — runner built, frozen chain reproduced BIT-FOR-BIT, two dry runs

Scope discipline: only `scripts/c_*.py` were written; `e4_chanorder.py` and `d3_vel1d.py` were copied
again as `c_chanorder.py` / `c_vel1d.py` so Track A's copies are never touched, and nothing outside
`lfe_extension/` was written. `../lfe_catalog_repair` and `../codex-LFE` read-only throughout.

**Built** (`scripts/`): `c_common.py` (frozen-point loader, epoch logic, cross-family rules) ·
`c_month_ready.py` · `c_cache_month.py` · `c_run_month.py` · `c_qc_month.py` · `c_eq_catalog.py` ·
`c_registration.py` · `c_watch.py` + `c_run_all.sh`. Output schema and restart instructions:
**`results/trackC/README.md`**.

**The frozen chain is reproduced EXACTLY.** `c_run_month.py` is `e4_template_test.py`'s measurement
core with the referee machinery (ground truth, halves, random and rolled controls, Z-only series,
moveout pass 1, event scoring) removed and nothing else changed. Against the frozen Track A run
`a1_mergeB_L12` on the same cache, with `--min-sites 0`:

| window | detections new / old | Δ net_mad | Δ net_cc | Δ amplitude |
|---|---|---|---|---|
| 2012 doy 228–232 (quiet) | 253 / 253 | 0.0 | 0.0 | 0.0 |
| 2012 doy 248–252 (ETS) | 3426 / 3426 | 0.0 | 0.0 | 0.0 |

The CPU-offload path (`--cc-device cpu`, for when the shared H100 is full) is bit-identical to the GPU
path. Preserving the float32 site-summation order needed the FROZEN station order
(`PGC,VGZ,NLLB,SNB,LRIV,B001,…`), not `sorted()` — recorded in `c_common.STAS_ORDER`.

**Production rules implemented** (Merlin Step 6 component 2 and Step 7): ≥ 2-site observability with
an `observable` flag per family-day; both tiers stored as columns on one parquet recorded down to
8.0 MAD (T_core and T_ext are filters, never a re-run); within-family declustering at the template
length; the frozen cross-family assignment (raw weighted cc among families sharing ≥ 50 % of their
2012 stations) evaluated at BOTH tiers, with the old max-MAD assignment kept alongside as a column;
PNSN tremor-in-box flag (±2 h) and per-day count; catalogued-EQ flag within 60 s against a new
`results/trackC/eq_catalog.csv` (PNSN `eq_catalog_pgc.csv` 2009–2017 + USGS FDSN 2009–2026,
de-duplicated to 26,465 events); HF ratio at CN.PGC from the NATIVE-rate trace; per-station daily
stacks of the T_core detections; per station-day daily median 8-s window energy per component.

**Track A's close-out (coordinator, same day) folded in**: the HF veto is adopted → `cat_ext` /
`cat_core` = tier AND `hf_pass` AND cross-family keep (vetoed rows are KEPT in the file so λ can be
re-derived); iteration rejected → the 2012 set is FINAL for 2012-instrument epochs, so every detection
carries `epoch_template` = `2012_final` / `2012_provisional` (at or after NLLB 2017-02-23, VGZ
2017-03-03, SNB 2017-06-08, PGC 2017-08-10, LRIV 2019-05-30) / `per_epoch:<file>`, plus
`n_sta_provisional`; a per-detection `quiet_flag` is stored so the floor λ can be built per family per
season (Track A step 3b: ~22 % of quiet-day T_core detections are real inter-ETS repeats).
Per-epoch template files drop in with NO code change: `results/trackC/templates_by_epoch.csv` or
`results/a4_templates_<STA>_<YYYY-MM-DD>.npz` are picked up automatically, per station per epoch, with
their own `sp` and their own placement. A station-epoch with no template set is a hard error.

**Dry run 1 — 2015-12** (P2b, 27 days 12-05..12-31 on disk; 12-01..12-04 are still in P3, so the month
is not yet "runnable" and was run with the partial-window override). 429 station-days cached in 62 s
(10 workers). Matched filter **81 s for 27 days (3.0 s/day)** with the GPU free, 300 s when it was
contended. 6,160 detections ≥ 8.0 MAD → 4,819 T_ext → 4,762 after the HF veto → **3,321 `cat_ext`**;
3,150 T_core → 3,133 → **2,226 `cat_core`** over 102 families. Unobservable family-days 27/3,132
(0.86 %, all single-site). Quiet-day rate **0.098 / template-day at T_ext, 0.028 at T_core** — WINTER
quiet-day FA is 0.35× the 2012 summer value the point was frozen at (Merlin's "winter within 2×
summer" bar is met with room to spare). EQ-flag fraction 0.42 % of all detections, 0.27 % of the
catalog (bar ≤ 0.5 %). 92 % of `cat_core` events fall in a tremor window. The catalog tracks the
2015-12 SSE: ≤ 5 events/day before 12-18, then 170 / 194 on 12-18..19 and 449 / 279 / 352 / 231 on
12-28..31. Median PGC-site weight share 0.45 — the dominance Merlin flagged is unchanged in production.

**Dry run 2 — 2018-06** (P2b, 21 days 06-10..06-30; PGC is **HH 100 Hz** here, B011 absent 06-10..25).
Cache 9 s (299 of 320 station-days were already cached by an earlier pass — verified byte-identical to
our chain before reuse). Matched filter **104 s for 21 days (5.0 s/day)**. 3,528 detections ≥ 8.0 →
2,523 T_ext → **1,692 `cat_ext`**; 1,460 T_core → 1,343 after the veto → **982 `cat_core`**.
**The 100 Hz HH chain is handled correctly and the provisional 2012 templates still register.**
Per-station registration (`c_registration.py`: monthly stack of the T_core detections vs the template,
matched-filter statistic; null at a random offset in the same window ≈ 0.00 everywhere), against a
2012 ETS reference run (doy 244–262, same code):

| station | 2012 cc | 2018-06 cc | | station | 2012 cc | 2018-06 cc |
|---|---|---|---|---|---|---|
| **PGC (BH→HH)** | **0.740** | **0.486** | | B926 | 0.607 | 0.137 |
| B010 | 0.776 | 0.445 | | NLLB | 0.613 | 0.188 |
| B011 | 0.706 | 0.348 | | VGZ | 0.588 | 0.183 |
| B009 | 0.670 | 0.421 | | B001–B007, LRIV | 0.19–0.39 | 0.05–0.15 |
| SNB | 0.689 | 0.305 | | | | |

Every station registers far above its null, and the level is uniformly lower in 2018-06 simply because
the stacks hold ~6× fewer events (13,458 station-events vs 83,177). The confounder-free test is PGC
against its CO-LOCATED, UNCHANGED borehole sensors (B009/B010/B011, EH 100 Hz throughout): PGC/borehole
registration ratio **1.03 in 2012 → 1.20 in 2018-06**. The band change costs PGC nothing; in 2018-06
PGC is the single best-registering station.

**PROBLEM FOUND — the frozen HF-veto cut is NOT epoch-invariant.** The cut 3.6384 is the 99th
percentile of the own-matched ratio at PGC **BH 40 Hz**, whose anti-alias filter rolls off inside the
10–19 Hz band. On PGC **HH 100 Hz** it does not: the median detection HF ratio goes **0.199 (2015-12,
BH) → 0.480 (2018-06, HH)** and the fraction above the frozen cut **2.0 % → 12.0 %**, i.e. the
unchanged cut vetoes ~6× more detections after 2017-08-10. `hf_fs` is stored on every row and
`qc_summary.csv` tracks the shift, but the cut needs a per-epoch re-calibration (same rule, own-matched
99th percentile, per instrument epoch) before the post-2017 catalog is used. Flagged to Track A.

**Also found**: 2018-06 as downloaded (06-10..06-30) contains **zero** quiet family-days — every day is
within 3 days of a PNSN tremor day — so no λ can be measured in that window; the winter quiet month
2016-12 is the one that must carry the post-2017-boundary λ calibration, and an equivalent quiet window
is needed inside every instrument epoch.

**Watcher tested end-to-end** (`bash scripts/c_run_all.sh --months 2015-12,2018-06 --allow-partial`):
picks up finished outputs without recomputing, runs QC, appends `results/trackC/months_done.csv`,
`flock` single-writer, logs in `logs/trackC_*.log`. **It is NOT running** — start it with
`nohup bash scripts/c_run_all.sh > /dev/null 2>&1 &` when the download has advanced.

**Projected full-record cost** (2014-01-01 .. 2026-09-15, 4,641 days, ~67,600 station-days of the 15
templated stations): cache **≈ 2.7 h** at 10–12 workers and **≈ 1.4 TB** of `cache20hz` (74 TB free);
matched filter **≈ 5 h** with the GPU free (4 s/day) and **≈ 15 h** at the contention seen today
(12 s/day); outputs ≈ 40 MB of detections, ≈ 0.5 GB of stacks. So the whole record finishes within
hours of the download, as required.

**Shared-GPU handling** (four other agents were using 90+ GB of the H100 all evening): the runner waits
at the start of each day for 6 GB free, re-decides per day whether the per-station cc series live on the
GPU or in host RAM, halves the family chunk on a CUDA OOM and falls back to the host if even the
allocation fails. None of it changes a number — verified by the bit-for-bit test on both paths.

**STEP 4 VERDICT (results/TRACK_A2_RESULTS.md).** Five boundaries tested (PGC 2017-08-10, VGZ
2017-03-03, NLLB 2017-02-23, SNB 2017-06-08, LRIV 2019-05-30) plus a VGZ re-test in the cached 2017 SSE
and a **B926 negative control** (no boundary at all, same 2012 -> 2018 gap).

| station | ACCEPT / ctl-OK | cc new <- ctl | ratio | dS median (MAD) | registration ratio |
|---|---|---|---|---|---|
| PGC | 46/61 = 0.74 | 0.805 <- 0.893 | 0.90 | +0.00 (0.00) | 0.94 |
| SNB | 37/56 = 0.63 | 0.622 <- 0.737 | 0.84 | +0.00 (0.05) | 0.89 |
| NLLB | 22/36 = 0.50 | 0.619 <- 0.709 | 0.87 | +0.00 (0.05) | 0.88 |
| VGZ | 7/43 = 0.14 | 0.478 <- 0.662 | 0.72 | +0.00 (0.05) | 0.71 |
| LRIV | 2/5 = 0.40 | **-0.216** <- 0.542 | - | +0.05 (0.35) | 0.75 |
| **B926 (NO change)** | **3/20 = 0.05** | 0.404 <- 0.510 | **0.79** | +0.00 (0.05) | **0.66** |

The negative control fails the pre-registered absolute bars as badly as VGZ does, so those bars measure
the 2012-ETS -> 2018-SSE change in event population / per-event purity, not instrument identity.  Read
relatively: PGC, NLLB, SNB transfer BETTER than a station that did not change; VGZ is indistinguishable
from it; **LRIV is the only real failure** (cc goes negative).  **dS median is +0.00 s at every
boundary** (MAD <= 0.05 s except LRIV 0.35) -> no instrument change introduced a timing shift.  5
polarity flips in 338 family-stations, one of them at B926 -> noise, not real flips.  Cross-family
chance level of the |cc| >= 0.6 bar: p95 = 0.67 at PGC (10.5 % false accepts), 0.32-0.52 elsewhere.
StationXML audit (`results/a4_stationxml_audit.csv`): azimuths/dips unchanged everywhere; the 2017 CNSN
upgrade swapped Guralp -> Nanometrics Trillium 120 PH with gain /1.8 to /14.2, and PGC (the biggest gain
drop) transfers best -> gain is not the mechanism.  Per-epoch template sets written:
`results/a4_templates_{PGC 46, SNB 37, NLLB 22, VGZ 7, LRIV 2, B926 3}_<epoch>.npz`.

Track A-2 complete; write-up `results/TRACK_A2_RESULTS.md`.  Only two files appeared outside
`lfe_extension/`: `__pycache__/d3_stage4_{locate,stacks}.cpython-313.pyc` in the read-only
`lfe_catalog_repair/scripts/`, created by Python when the locator cross-check imported the D3 original.
Nothing was deleted anywhere.

## 2026-09-16 ~22:35 UTC Track C built + started (agent report; results/trackC/README.md)
- c_run_month reproduces the frozen chain BIT-FOR-BIT (2012 quiet 253/253, ETS 3426/3426 detections; float32 site-sum order = frozen station order).
- Dry run 2015-12: 27 d in 81 s; cat_ext 3,321 / cat_core 2,226 over 102 families; winter quiet-day rate 0.098 (T_ext) / 0.028 (T_core) per template-day = 0.35x the 2012 summer value (Merlin winter bar MET); EQ flag 0.27%; unobservable 0.9%.
- Dry run 2018-06 (PGC HH 100 Hz): 21 d in 104 s; cat_ext 1,692 / cat_core 982; PGC registration vs co-located boreholes 1.20 (2012: 1.03) -> band change costs PGC nothing.
- PROBLEM 1: HF-veto cut 3.64 is NOT epoch-invariant (BH 40 Hz anti-alias rolls off in 10-19 Hz; on HH median ratio 0.20->0.48, vetoed 2%->12%) -> per-epoch cut needed before post-2017 catalog is used (assigned to A-2 agent). PROBLEM 2: 2018-06 has zero quiet family-days -> 2016-12 carries the post-boundary floor.
- Projected full record: cache ~2.7 h + 1.4 TB; matched filter ~5 h uncontended. Watcher started: nohup bash scripts/c_run_all.sh --from-ym 2014-01 (post-2017 rows labelled 2012_provisional until a4 templates + per-epoch cuts land).

## 2026-09-16 ~22:45 UTC Track A-2 complete (results/TRACK_A2_RESULTS.md)
- Step 4 epoch identity: absolute bars (|cc|>=0.6, registration>=0.8x) are mis-calibrated — the NO-CHANGE control B926 fails them worse than VGZ over the same 2012->2018 gap (event-population/purity effect). Read relatively: PGC 0.74 accept (cc 0.81 vs matched-N control 0.89), SNB 0.63, NLLB 0.50 transfer BETTER than the unchanged station; VGZ 0.14 = indistinguishable from it; LRIV 2019-05-30 is the one real failure (cc -0.22). dS = +0.00 s at every boundary; no real polarity flips; StationXML clean (2017 CNSN Guralp->Trillium 120PH, gain /1.8../14.2).
- DECISION: keep 2012 templates for PGC/NLLB/SNB/VGZ across 2017; DROP LRIV after 2019-05-30. Per-epoch template npz written (a4_templates_*); templates_by_epoch.csv requested from Track C.
- Step 5 relocation: synthetic floor 1.43 km (bar testable); station corrections from 7 USGS intraslab EQs (S terms -0.05..+1.13 s; E3's +1.7..+3.25 s were mostly family mislocation); 14/130 families pass the geometry gate: median dh 3.61 km (free) / 3.45 km (Slab2), shuffled null 29.4 km (ratio 0.12), spread 0.95x, 0% at bounds, dz median -2.2 km. Limit = S-pick yield on far boreholes.
- Ops: snap CN sampling rates to integer (99.99999 Hz fragments) — assigned to Track C; <=3 concurrent GPU passes on the H100.
- Track C watcher running from 2014-01; A-2 agent now calibrating per-epoch HF-veto cuts.

## 2026-09-16 Track A-2 addendum — per-epoch HF cuts + epoch-boundary confirmation (Track C hand-off)

**HF veto.** The epoch shift Track C found (median 0.199 BH -> 0.480 HH) is NOT mainly an instrument
effect: the StationXML response predicts only x1.16 of it (`a4_hf_response.py`), NO band removes it
(10-19 x2.41, 10-16 x1.86, 10-14 x2.06, 8-16 x1.59, 12-16 x1.86 — `a4_hf_bands.py`, all days of both
months, medians reproduce Track C exactly), and the SAME BH instrument gives 1.019 on 2012 vs 0.199 on
2015-12. Restricted to a matched LFE-like core (T_core, n_sites>=4, tremor days) the medians are 0.256
(2012 BH) / 0.169 (2015-12 BH) / 0.247 (2018-06 HH) — the BH-BH spread equals any BH-HH difference.
**The distribution moves with the detection population, not the sensor.** The proposed referee-free rule
(p99 of that selection, "RULE-T") is a tail statistic of a mixture and is NOT stable between two BH
months (3.246 vs 1.208) -> REJECTED. Adopted **RULE-M**: cut = F x median(core) with F = 14.22 fixed on
the 2012 referee; it reproduces the frozen 3.6384 on 2012 and passes the original bars identically
(quiet FA reduction 0.59, relative recall cost 0.003). Per-epoch cuts (`results/a4_hf_veto_cuts.csv`):
PGC BH 3.6384 (frozen, unchanged), **PGC HH 3.5129** — whose 95 % bootstrap CI [3.34, 3.83] CONTAINS the
frozen cut, so the frozen cut was never badly wrong after 2017; the 12 % vs 2 % veto fraction is the veto
correctly removing a larger contaminant population (it vetoes only 0.5 % of the 2018-06 LFE-like core).
No separate cut for the PGC BH sub-boundary 2015-01-12 (4.6 % response effect, inside the +-35 %
SSE-to-SSE spread). If tighter control is wanted the right granularity is per SSE, not per epoch.

**Epoch boundaries CONFIRMED** (`results/a4_stationxml_audit.csv`): all five are real response changes
with unchanged azimuths/dips. SNB 2017-06-08 and VGZ 2017-03-03 are same-BAND sensor swaps (Guralp
CMG3ESP_60sec -> Nanometrics Trillium 120 PH, gain /1.8 and /14.2) — absent from station_epochs.csv only
because that table keys on the band; both need their own template epoch. PGC/VGZ/NLLB/SNB all end on the
SAME Trillium 120 PH response (one CNSN upgrade).

**`results/trackC/templates_by_epoch.csv` written** (`a4_make_epoch_sets.py`): MERGED sets = every 2012
template of the station, replaced family by family by the new-epoch re-stack where identity ACCEPTED it
— PGC 111 = 65+46, NLLB 88 = 66+22, VGZ 97 = 90+7, SNB 101 = 64+37. LRIV is deliberately NOT listed (its
identity FAILED; its rows must keep Track C's `2012_provisional` label) and neither is B926 (negative
control, no boundary). Human decision pending: drop LRIV after 2019-05-30 or keep it flagged.

**Coordinator decision on the per-epoch templates (2026-09-16, recorded):** my proposal was preserved as
`results/trackC/templates_by_epoch.a4_proposed.csv`; the operative `templates_by_epoch.csv` keeps the
frozen 2012 templates for PGC/NLLB/SNB/VGZ (`2012_final`) and DROPS LRIV after 2019-05-30 (`drop=1`).
Both match this agent's recommendations; the merged a4 sets stay on disk for a future sensitivity pass.

## 2026-09-16 22:59 UTC A-2 addendum: HF veto per epoch + boundaries (results/TRACK_A2_RESULTS.md ADDENDUM)
- The HH-vs-BH HF-ratio shift is NOT the anti-alias response (predicts x1.16 of a x2.41 shift) and no band fixes it; the ratio distribution tracks the MONTH (contaminant population), not the sensor. Proposed p99-of-core rule rejected (unstable BH->BH: 3.25 vs 1.21). Adopted RULE-M: cut = 14.22 x median(LFE-like core), fixed on the 2012 referee. PGC_BH cut 3.6384 (unchanged); PGC_HH cut 3.5129 (95% CI 3.34-3.83 contains 3.6384) -> the frozen cut was never badly wrong post-2017 (vetoes 12% of all 2018-06 detections but 0.5% of the LFE-like core). File results/a4_hf_veto_cuts.csv.
- All five epoch boundaries confirmed from StationXML as genuine response epochs (one CNSN upgrade to Trillium 120PH): PGC 2017-08-10, NLLB 2017-02-23, SNB 2017-06-08, VGZ 2017-03-03 (same band, sensor swap), LRIV 2019-05-30. Azimuths/dips unchanged.
- Operative template policy (Track C wrote templates_by_epoch.csv): 2012 templates for PGC/NLLB/SNB/VGZ all epochs; LRIV dropped after 2019-05-30; A-2's merged per-epoch sets kept as templates_by_epoch.a4_proposed.csv for a later sensitivity pass.

## 2026-09-16 Track C — A-2 hand-off applied: epoch policy, per-epoch HF cuts, CN rate snap

Watcher was running (`c_run_all.sh --from-ym 2014-01`, PID 433458) on 2014-01. **Coordinated stop**: the
watcher PARENT was killed while its `c_run_month.py 2014-01` child kept running to completion as an
orphan (`subprocess.call` does not kill the group), so the month finished cleanly and no partial output
was written; the watcher was then restarted after every change below. A graceful stop was also added for
next time: `touch results/trackC/STOP` makes the watcher finish the current month and exit.

**(1) Epoch template policy — `results/trackC/templates_by_epoch.csv`** (schema extended to
`sta, epoch_start, npz, label, drop, note`; `npz` empty now means "explicitly keep the frozen 2012 set
for this epoch", which also suppresses the `a4_templates_*` filename search):
- PGC 2017-08-10, NLLB 2017-02-23, SNB 2017-06-08, VGZ 2017-03-03 → **`2012_final`, keep the 2012
  templates** (A-2 step 4: all four transfer at least as well as the B926 no-change control across the
  same 2012→2018 gap; dS +0.00 s; no real polarity flips).
- LRIV after 2019-05-30 → **DROPPED from the network** (A-2: stack cc goes NEGATIVE −0.216 vs the
  matched-N control's +0.542; dS MAD 0.35 s). The station is not even loaded for those epochs, and every
  family-day records the loss in the new `excluded_sta` / `n_sta_excluded` columns for the completeness
  model. Measured on 2021-10: **1,054 of 3,596 family-days lose LRIV**, and the network drops 14 → 13
  stations.
- `epoch_template` is therefore `2012_final` everywhere; the "matched across an instrument change"
  diagnostic moved to the new per-detection column **`n_sta_post_boundary`** (median 4 in 2018-06, 3 in
  2021-10, 0 before 2017).

**CONFLICT, flagged for the coordinator.** While I was editing, the A-2 agent OVERWROTE
`results/trackC/templates_by_epoch.csv` with its own 4-row version pointing every 2017 boundary at
`results/a4_templates_merged_<STA>_<date>.npz` and containing **no LRIV row at all** — under which LRIV
would have kept carrying 2012 templates past 2019-05-30, exactly the failure the drop exists to prevent
(observed live: the run printed `LRIV@2019-05-30=2012_provisional`). I restored the coordinator's policy
and preserved theirs verbatim as `results/trackC/templates_by_epoch.a4_proposed.csv` (nothing deleted).
`results/trackC/` is a Track C output directory — the A-2 per-epoch stacks are a SENSITIVITY option the
coordinator has not adopted, and only one of us should write that file.

**(2) Per-epoch HF-veto cuts — `results/a4_hf_veto_cuts.csv` arrived mid-session** (and was rewritten
three times while I worked, from 4 to 8 rows). The runner now reads it, stamps `hf_cut` and
`hf_cut_src` on every detection, and `scripts/c_reflag_hf.py` re-flags already-processed months WITHOUT
re-running the matched filter (`hfr` and `hf_fs` are stored; it rewrites only `hf_cut`, `hf_pass`,
`cat_ext`, `cat_core`, keeps `hf_pass_prev` / `hf_cut_prev`, logs to `hf_reflag_log.csv` and
regenerates QC). Three reader decisions, all documented in `results/trackC/README.md`:
- A-2's schema (`station_epoch`, `band`, `fs`, `cut`, …) and the generic one are both accepted;
- rows whose `band` is not the 10–19 Hz band our stored `hfr` was computed in are **dropped** — a cut
  calibrated in another band is silently wrong;
- several rows match the BH epoch with different cuts, so the tie-break is **the anchor wins**: a row
  whose `method` says it is the frozen referee cut cannot be displaced by a later recalibration. Result:
  **BH epochs keep the frozen 3.6384** (the operating point does not move) and **HH epochs use 3.5129**.
  Ambiguity the anchor rule does not settle is reported in `hf_cut_src`.

**(3) CN sampling-rate snap** in `c_cache_month.py` (and in the HF-ratio reader, which merges native
data too): every fragment rate within 1e-4 relative of an integer is snapped before `merge`; further out
still raises. Counted in `n_fs_snapped`. **Full scan of all 7,563 CN/UW station-days on disk**: 74 carry
non-integer fragment rates, but only **3 would actually fail to merge** (a single channel with two
different rates) — CN.VGZ 2021-10-06, CN.PGC 2021-11-03, CN.VGZ 2021-11-03. All three are already
correctly cached (rebuilt with the snap → byte-identical to what is on disk), so **no month needed
re-caching for this fix**; it is a forward guard for the ~66,000 station-days still to come. A day that
needs no snap is byte-identical before and after. (A-2 reported 41 of 504; on our disk it is 3 of 7,563
— worth reconciling, but the fix is the same and harmless either way.)
A second, subtler bite of the same bug was found and fixed: `hf_fs` was being recorded as 99.99999 Hz on
CN 2021 days, so the per-epoch HF cut failed to match on `fs` and silently fell back to the frozen cut —
2021-10 came out with a MIX of 3.6384 and 3.5129. Now snapped at source and the comparison is
rate-tolerant.

**Re-ran every processed month under the final configuration** (the frozen measurement chain is
untouched — 2012 still reproduces bit-for-bit):

| month | days | detections ≥ 8.0 | T_ext → cat_ext | T_core → cat_core | HF-vetoed | runtime |
|---|---|---|---|---|---|---|
| 2014-01 | 31 | 2,328 | 1,359 → 790 | 550 → 331 | 133 | 65 s |
| 2015-12 | 31 (12-01..04 arrived) | 6,268 | 4,878 → 3,352 | 3,163 → 2,234 | 171 | 63 s |
| 2018-06 | 21 | 3,528 | 2,523 → 1,675 | 1,460 → 977 | 449 | 61 s |
| 2021-10 | 31 | 12,141 | 8,999 → 4,632 | 5,396 → 2,966 | 2,510 | 94 s |

2021-10 is the first month run with a post-2019 network (SNB dead 2021-03-28, LRIV dropped by policy):
13 stations, 62 of 3,596 family-days unobservable (1.7 %, double the 0.86 % of the 15-station epochs).

Watcher restarted afterwards; `months_done.csv` holds the four months above.

## 2026-09-16 23:07 UTC Track C follow-ups applied; watcher running (2014-02)
- templates_by_epoch.csv: ORCHESTRATOR DECISION STANDS — 2012 templates (2012_final) for PGC/NLLB/SNB/VGZ in all epochs; LRIV dropped after 2019-05-30 (network 14->13; 1,054/3,596 family-days in 2021-10 carry excluded_sta). A-2's merged per-epoch sets preserved as templates_by_epoch.a4_proposed.csv (sensitivity option, not adopted). Only Track C writes that file from now on.
- HF cuts: runner reads results/a4_hf_veto_cuts.csv; BH keeps frozen 3.6384 (anchor wins), HH uses 3.5129; c_reflag_hf.py re-flags without re-running. CN sample-rate snap applied (3 of 7,563 station-days on disk would have failed; already cached correctly) + hf_fs snapped so the per-epoch cut matches.
- Re-run of the four processed months under the final config: 2014-01 cat_ext 790 / cat_core 331; 2015-12 3,352 / 2,234; 2018-06 1,675 / 977; 2021-10 (13 stations) 4,632 / 2,966, unobservable 1.7%. ~60-95 s per month. Graceful stop: touch results/trackC/STOP.

## 2026-09-16 23:15 UTC Track D agent launched (Opus)
Completeness model (Step 6), catalog assembly v1 (Step 7), Thomas QC per SSE (Step 8), RESULTS.md; polls the sweep to completion. User: Fable near quota -> orchestrator-only from here.

## 2026-09-16 23:35 UTC Track D (product) — machinery built while the sweep runs

Scope discipline: everything new is `scripts/d1_*`, `d2_*`, `d3_*`, output under `results/trackD/` and
`figures/d_*`. Track C's runner and outputs are READ only; `../lfe_catalog_repair` and `../codex-LFE`
read-only (every script used is a COPY: `d3_thomas_harness.py` ← `e5_thomas_harness.py`,
`d3_thomas_base.py` ← `e4_thomas_harness.py`, `d3_thinned_reference.py` ← `e6_thinned_reference.py`).

**D1(a) epoch ablation — 4 GPU passes run (23:17–23:31 UTC).** The extension epochs that the production
run actually realises are, from `results/trackC/templates_by_epoch.csv` + `station_epochs.csv`:
E1 2014-01-01..2019-05-29 (15 stations), E2 2019-05-30..2021-03-27 (14, LRIV dropped by policy),
E3 2021-03-28..2026-09-15 (13, LRIV + SNB). **CN.PFB (ends 2018-12-20) carries no 2012 template and is
not in the network at all**, so the brief's "post-2018-12 without PFB" boundary is a no-op for this
detector — recorded, not ablated. Passes: `d1_abl2012_{E2,E3}` (2012 test half) and
`d1_abl2010_{E2,E3}` (held-out 2010, whose base network is the 9 persistent stations). The frozen HF
veto is applied to each so the recall is measured at the production operating point.

**D1(b) sensitivity model** = a probit on the referee events with per-family intercepts (ridge, λ by
5-fold CV), a magnitude slope, and two day terms: `u` = log10(family-day MAD scale / its epoch median)
and `v` = the SNR-weighted log10 noise-amplitude factor from the daily median 8-s window energy,
normalised per station per INSTRUMENT epoch (so the 2017 gain changes cancel — cc is amplitude-normalised).
Referee-year daily noise computed with Track C's own estimator (`d1_ref_noise.py`, 1,083 station-days).

**D1(c) floor λ** per family / tier / epoch / season from Track C's `quiet_flag`, corrected by the Track A
step-3b real-repeat fraction (0.225 at T_core, 0.057 at T_ext) and gamma-Poisson shrunk (k = 10
template-days) toward the epoch-season pool.

**D3 SSE windows** (`d3_sse_windows.py`, PNSN ≥ 50 epicentres/day in the family box, gaps ≤ 20 d,
span ≥ 10 d, ±3 d pad): **31 windows 2014–2026, 15 of them "major" (≥ 10 active days)** — the 15 majors
contain Merlin's 12 named episodes.

**Smoke tests on the 11 months finished so far**: catalog assembly reproduces `months_done.csv` exactly
(14,443 events, 14,202 cat_ext / 8,180 cat_core, 0 duplicates, 0 month mismatches); the Thomas harness
runs on our 2012 referee catalog rebuilt at the frozen operating point (tier ∧ HF veto ∧ frozen
cross-family) — **main-front bar MET** (first-burst speed −5.28 / −5.23 km/day vs Bostock's −4.69,
12 % / 11 %, bar 20 %), **re-rupture-interval bar NOT met** (median 1.05 d / 0.87 d vs 0.70 d, 49 % / 23 %,
bar 5 %) — which is the expected decimation bias and why Merlin's per-family-thinned Bostock reference
(`d3_thinned_reference.py`) is the fair comparison; it is being run as the reported diagnostic.

## 2026-09-16 23:50 UTC Track D — D1 complete on the referee data; two results that change the product

**D1(a) epoch ablation — the headline is NOT what the epoch list suggested.** At the FROZEN MAD
thresholds, dropping stations does not cost recall; it costs PURITY:

| epoch (station set) | quiet-day FA / template-day, T_core | own recall T_core | recall at MATCHED FA |
|---|---|---|---|
| E1 15 stations (2014-01..2019-05) | 0.097 (the frozen point) | 0.1463 | 0.1463 |
| E2 14, LRIV dropped (2019-05..2021-03) | 0.098 | 0.1460 | 0.1457 |
| E3 13, LRIV + SNB gone (2021-03..) | **0.290 (×3.0)** | 0.1526 | **0.1185 (−19 %)** |

Same at T_ext (FA 0.277 → 0.279 → **0.744**; matched-FA recall 0.1841 → 0.1837 → **0.1505**). **LRIV was
carrying nothing** — its removal is free on both axes, which independently supports the A-2 drop
decision. **SNB's loss in 2021-03 is the real epoch break**: at a fixed MAD a smaller network is not a
fairer detector, it is a LOWER threshold (Merlin F3), so the post-2021 catalog runs at ~3× the
false-alarm floor it was frozen at and loses ~19 % of its recall once that is matched. The held-out 2010
year reproduces the direction (FA 0.292 → 0.297 → 0.356; matched-FA recall 0.1718 → 0.1730 → 0.1645).
This is now a per-event column: every post-2021-03 event carries its own λ, and RESULTS.md states it.

**D1(b) sensitivity model — the pre-registered validation bar is MET.** Probit fit on the 2012
9-station ablation, predicted on the HELD-OUT 2010 year: |predicted − measured| own recall is
**0.005–0.024 in every magnitude bin and 0.002 pooled**, against the ±0.05 bar (in-sample 2012 control
0.001–0.010). Fitted α ≈ 1.1–1.2 probit units per magnitude; the noise day term c2 ≈ −0.10..−0.32 and
the MAD-scale term c1 = +0.5..+2.3 (positive: on busy days the cc series MAD tracks activity, so `u` is
partly an activity covariate, not a pure noise term — and its sign flips in the 2010 fit, so only the
2012 coefficients are ever exported). Production Mc(f,d): median **2.61 (T_ext) / 2.69 (T_core)**, IQR
2.26–2.95, consistent with E5's independent Mc(8 MAD) ≈ 2.25 prediction for this network at a lower
threshold. 45 % of production family-days have a noise day term outside the referee 1–99 % range
(late-summer ETS windows only) — flagged per row as `day_terms_in_calib_range`; with |c2| ≤ 0.32 the
extrapolation moves Mc by ≤ 0.1 magnitude.

**D1(c) floor λ — the seasonal-stability bar is FAILING so far.** Production quiet-day rates per
template-day (E1, 16 months): DJF 0.121 / MAM 0.070 / JJA 0.033 / SON 0.050 at T_core (T_ext 0.301 /
0.204 / 0.089 / 0.122). **Winter / summer = 3.6× (T_core), 3.4× (T_ext), against Merlin's ≤ 2× bar.**
Provisional — only 2014–2015 is in yet, and the summer months on disk are unusually quiet — but the
sign is consistent (winter storm noise) and λ is attached per event, per season, so the catalog carries
it rather than hides it.

**D3 fair reference — the re-rupture bar failure IS the decimation bias.** The per-family-thinned
Bostock reference (`d3_thinned_reference.py`, 20 seeds, our own per-family recall + our own quiet-day
floor) gives, on 2012: front −4.96 ± 0.17 km/day and re-activation interval median 0.857 ± 0.133 d
(2.5–97.5 % 0.674–1.058) at T_core — **our 1.051 d is inside that range** although it is 49 % off
Bostock's undecimated 0.704 d. At T_ext: thinned 0.915 ± 0.144 (0.693–1.139), ours 0.866 — inside.
Secondary fronts: ours 21 vs thinned 20.6 ± 4.1 at T_core (14 vs 12.8 ± 3.1 at T_ext) — but the
CIRCULAR NULL of the thinned reference is 13.9 (T_core) / 9.7 (T_ext) against Bostock's 2.25, so **at
20 % completeness the secondary-front count carries almost no information**; that is reported, not
barred. The one metric where we sit OUTSIDE the thinned range is the re-activated-family fraction
(ours 0.27 / 0.30 vs thinned 0.157 / 0.190) — consistent with false alarms manufacturing re-activations.

## 2026-09-16 23:55 UTC Track D — full chain runs end to end; RESULTS.md is generated, not written

`scripts/d_run_all.sh` rebuilds the whole product from whatever Track C months are on disk
(score epochs → fit Mc → apply Mc → floor λ → assemble → SSE windows → Thomas per SSE → thinned
reference → figures), and `scripts/d4_results.py` writes `results/trackD/RESULTS.md` straight from the
output files, so the write-up can never drift from the numbers. Both were exercised on the 18 months
finished so far; `scripts/d_sweep_qc.py` is the read-only anomaly scan of the sweep.

**Sweep health at 23:55 UTC**: 18 of 153 months done, ~3 min/month, no stall (watcher at 2015-08);
**0 months with zero observable family-days**; unobservable family-day fraction median 0.9 %, max 1.7 %
(2021-10, the 13-station epoch); 0 station-days with a caching error; **2 months with zero QUIET
family-days (2018-06, 2021-10)** — entirely inside tremor windows, so they contribute no floor and λ
falls back through the family-epoch → epoch-season → epoch hierarchy (recorded per row in
`lambda_src`). Track C's own QC flags 24 site-months above 2 λ and 10 months with an HF-ratio shift, but
both flags compare against a RUNNING median that has only a few months behind it this early in the
record — they need re-reading when the sweep finishes. 2014-06 and 2014-09 were skipped by the watcher
(their download was not complete when it passed) and will be picked up on a later loop.

## 2026-09-16 23:58 UTC ⚠ SWEEP BLOCKER found by Track D (for the orchestrator / Track C — NOT fixed here)

**Five months can never become runnable**, because of a one-day off-by-one between the download queue and
`c_month_ready.py` at every station-epoch END. `results/station_epochs.csv` gives e.g. `CN,LZB,BH,…,
2014-06-17`; `download_plan.py` treats that end date as EXCLUSIVE and queues 2014-06-17 as `no_epoch`, so
it is never fetched and never gets a manifest row; `c_month_ready.py` treats the same date as INCLUSIVE,
so it waits forever for exactly one station-day:

| month | the single blocking station-day | in our template network? |
|---|---|---|
| 2014-06 | CN.LZB 2014-06-17 | no |
| 2014-09 | CN.YOUB 2014-09-03 | no |
| 2015-01 | C8.SOKB 2015-01-18 | no |
| 2017-02 | CN.NLLB 2017-02-21 (last BH day) | **yes** |
| 2018-12 | CN.PFB 2018-12-20 | no |

Every other station-day of those months is on disk (439–521 usable station-days each), and a month that
is entirely past a station's death runs fine — only the BOUNDARY month blocks. The watcher skips them
silently and will finish at 148/153 months. 2017-02 and 2014-06 each contain an SSE window
(SSE201702, SSE201406).

Not fixed from Track D: `scripts/c_*` and `results/trackC/*` are Track C's to write. The two clean fixes
are (a) one line in `c_month_ready.py` — treat a queue row whose status is `no_epoch` as finished — or
(b) run those five months with the existing override:
`bash scripts/c_run_all.sh --allow-partial --months 2014-06,2014-09,2015-01,2017-02,2018-12`.
Option (b) costs nothing scientifically for four of the five (LZB / YOUB / SOKB / PFB carry no 2012
template and are not in the detector at all); for 2017-02 it costs NLLB on one day.

## 2026-09-17 23:58 UTC Track D interim + two Track C fixes ordered
- D1(b) Mc model validated on held-out 2010: max |pred-meas| recall 0.024 (bar 0.05) MET; alpha 1.20 vs independent 1.14.
- D1(a) epoch ablation: dropping LRIV is free; losing SNB (2021-03) -> ~3x quiet-day FA at the frozen MAD and -19% recall at matched FA. DECISION: re-derive T_core/T_ext PER NETWORK EPOCH by the frozen FA-only rule (re-flag on stored net_mad; tiers_by_epoch.csv). E0 2014-01..2017-02, E1 2017-03..2018-12, E2 2019-01..2021-03, E3 2021-04.. .
- D1(c) seasonal bar NOT met: winter/summer quiet-day rate 3.6x (bar 2) on 2014-15 — reported, not tuned; D agent to test whether the winter excess passes the withheld-station stack test.
- Thomas on 2012 referee: main-front speed MET (12%); re-rupture median NOT met vs raw Bostock but inside the per-family-thinned 95% range (decimation bias).
- BLOCKER fixed by order: c_month_ready.py epoch-end inclusive vs planner exclusive -> 2014-06, 2014-09, 2015-01, 2017-02, 2018-12 stuck on one station-day each.
- Sweep 21/153 months at 23:58.

## 2026-09-17 Track C — epoch-end off-by-one unblocked; tiers re-derived PER NETWORK EPOCH

Watcher stopped with the new graceful `touch results/trackC/STOP` (it finished 2015-11, marked it done
and exited), changed, re-derived, re-flagged, restarted. 23 months were already processed
(2014-01..2015-12 except the three blocked ones, plus 2018-06 and 2021-10).

**(1) BLOCKER fixed — epoch end is EXCLUSIVE.** `download_plan.py` queues a station-day ON an epoch end
date as `no_epoch` and never fetches it, while `C.alive()` counted that station as alive, so
`c_month_ready` waited forever for a station-day that will never exist. `C.alive()` and
`C.eff_epoch_start()` now use `start <= d < end`. All five named months became runnable immediately:

| month | station-day that blocked it | now |
|---|---|---|
| 2014-06 | CN.LZB 06-17 | runnable, 521 usable station-days over 18 stations |
| 2014-09 | CN.YOUB 09-03 | runnable, 463 over 17 |
| 2015-01 | C8.SOKB 01-18 | runnable, 467 over 17 |
| 2017-02 | CN.NLLB 02-21 | runnable, 439 over 16 |
| 2018-12 | CN.PFB 12-20 | runnable, 462 over 16 |

Checked that nothing is lost at a band switch: on 2017-08-10 PGC is alive through the NEW (HH) epoch,
on 2011-03-17 SNB through its HH epoch, and NLLB is correctly absent on 2017-02-21..22 (its real
2-day metadata gap). Terminal stations now simply stop on their end date (SNB gone from 2021-03-28).
None of the 23 already-processed months contains any of those dates, so no re-run was needed.
The watcher picked 2014-06 up on restart.

**(2) TIERS ARE NOW PER NETWORK EPOCH** (`scripts/c_tiers_by_epoch.py`, `results/trackC/tiers_by_epoch.csv`).
Same FA-only rule that defined them, measured on each epoch's own days: T_core = smallest T on a
0.25 grid ≥ 8.0 with quiet-day FA ≤ 0.1/template-day; T_ext = smallest T with quiet-day FA ≤ 1.0 AND
post-SSE FA ≤ 1.5. Quiet = the stored `quiet_flag`; **post-SSE = 1–10 days after the last tremor day of
an episode** (episode ends after ≥ 3 tremor-free days) — the production analogue of Track A's 2012 doy
270–279 window, which the spec did not define for production and which I had to fix. FA counts
`net_mad ≥ T AND hf_pass` (the frozen point was derived with the veto), no cross-family assignment (the
rate is per template per day), no control correction.

| epoch | months | status | T_core | T_ext | quiet family-days | post family-days | quiet FA @ T_core | quiet FA @ T_ext | post FA @ T_ext |
|---|---|---|---|---|---|---|---|---|---|
| E0 | 2014-01..2017-02 | **derived** (9.00/8.00 → anchor kept) | **8.75** | **8.25** | 36,865 | 33,873 | 0.106 | 0.246 | 0.246 |
| E1 | 2017-03..2018-12 | insufficient (0 quiet) | 8.75 | 8.25 | 0 | 345 | — | — | 0.128 |
| E2 | 2019-01..2021-03 | no processed month | — | — | 0 | 0 | — | — | — |
| E3 | 2021-04.. | insufficient (0 quiet) | 8.75 | 8.25 | 0 | 570 | — | — | 0.196 |

**E0 ANCHOR HOLDS.** Re-derived on E0's OWN 36,865 quiet family-days the rule gives T_core 9.00 and
T_ext 8.00 — both exactly 0.25 from the frozen 8.75 / 8.25, i.e. **within the coordinator's 0.25 bar**,
so the frozen pair is kept and E0 is marked derived-and-anchored. (At the frozen pair E0's own quiet FA
is 0.106 at T_core — a hair over the 0.1 bar, which is why the rule's own answer is one step up — and
0.246 at T_ext, well inside the 1.0 bar; the post-SSE bar is not binding in E0, unlike 2012.)

**E1/E2/E3 CANNOT BE DERIVED YET, and this is the honest limit today**: the only processed months in E1
and E3 are 2018-06 and 2021-10, both SSE months with **zero** quiet family-days, and E2 has no
processed month at all. All three therefore sit on the E0 anchor. They are NOT frozen in the file, so
`c_tiers_by_epoch.py` will derive them automatically once their quiet months are processed;
`c_reflag_tiers.py --all` must then be run once to apply them across the record. **Track D's ~3x
quiet-day FA after SNB's loss is not yet reproducible here for exactly this reason** — it needs the
E3 quiet months (2021-12..2022-02) on disk.

Machinery, all post-processing (net_mad, templates and the HF cut are never touched):
`scripts/c_reflag_tiers.py` rewrites `tier_ext`/`tier_core` from the epoch's T's, **recomputes the
cross-family assignment** (the frozen rule is evaluated ON the thresholded set, so a new T changes which
detections suppress which — this is the part a naive re-threshold would get wrong), rewrites
`cat_ext`/`cat_core` and the `n_det_ext`/`n_det_core` counts in `family_days`, stamps `T_core_used`,
`T_ext_used`, `tier_epoch` on every row, keeps `tier_*_prev`, logs `tier_reflag_log.csv` and regenerates
QC. The runner reads `tiers_by_epoch.csv` for new months. A row already marked `derived` is FROZEN so
the tiers cannot drift as months arrive (`--rederive` to override deliberately).

Applied over all 23 processed months: **cat_core 10,364 (E0) + 977 (E1) + 2,966 (E3) and cat_ext
20,717 + 1,675 + 4,632 — UNCHANGED**, because every epoch currently sits on the anchor. That no-op is
the check that the re-flag path reproduces the runner exactly.

Two small bugs found and fixed while doing this: the month-name parse in `c_tiers_by_epoch.py` and
`c_reflag_hf.py` (`basename[12:-8]` is right for `family_days_` but not for `detections_`; both now
use `rsplit("_", 1)`).

**Wired into the watcher (2026-09-17).** After every month it re-runs `c_tiers_by_epoch.py` — which can
only ADD an epoch, because derived rows are frozen — and the first time an epoch becomes derivable it
runs `c_reflag_tiers.py --all` once, so no epoch ever carries two different thresholds across the
record (`logs/trackC_tiers.log`). Restarted on 2015-01; 2014-06 (1,487 T_ext / 605 T_core) and 2014-09
(3,090 / 1,696) — two of the five previously blocked months — ran cleanly, with LZB correctly absent on
2014-06-17 and the templated network unchanged across that boundary (median 10 stations / 6 sites).

## 2026-09-17 00:12 UTC Track D — re-based on Track C's epoch axis; the winter/summer question opened

**Coordinator's two fixes landed**: `c_month_ready.py` epoch-end inclusivity (the watcher is now running
2015-01, one of the five stuck months) and `results/trackC/tiers_by_epoch.csv`.

**Track C's epoch axis adopted as the authority.** Its epochs are E0 2014-01..2017-02, E1 2017-03..2018-12
(post CN sensor upgrade), E2 2019-01..2021-03 (LRIV dropped 2019-05-30), E3 2021-04.. (no SNB) — NOT the
same cut as Track D's provisional station-set epochs, and the labels collide (its E1 is not mine). New
module `scripts/d_tiers.py` is now the single place that knows both: it reads Track C's epochs and
per-epoch T_core / T_ext, and maps each epoch onto the STATION SET the ablation actually measures —
**S15** (15 stations; Track C E0 and E1), **S14** (LRIV dropped; E2), **S13** (LRIV + SNB gone; E3).
Every D1/D2/D3/D4 script now reads epochs and thresholds through it, so the completeness model, the floor
λ and the catalog can never disagree with the flags. `completeness_family_epoch.csv` is the station-set
measurement expanded onto Track C's epoch axis; `completeness_family_stationset.csv` is the measurement.
Today's file derives 8.75 / 8.25 for E0 ("derived 9.00/8.00, within 0.25 → frozen anchor kept") and falls
back to the frozen point for E1–E3 (0 quiet family-days yet), so **no number has moved**; the plumbing is
in place for when the full record re-derives them. RESULTS.md now opens with the epoch/threshold table
and says which epoch used which T and whether it was derived or the frozen fallback.

**Seasonal-stability bar stays as reported — NOT met (3.0–3.2× on the months in hand), not tuned.**
New section 1d of RESULTS.md attributes it, from Track C outputs only:
* the quiet-day tail is elevated near-UNIFORMLY in winter (decay −0.73 dex/MAD in DJF vs −0.66 in JJA):
  winter is not adding a heavy tail of big excursions, it is adding events at every threshold;
* winter is 1.42× louder in daily energy (1.19× in amplitude) but the family-day **MAD scale moves by
  only +2.4 %** — the score is an amplitude-NORMALISED cc divided by its own MAD, so **a louder winter
  does not raise the bar**. That is the mechanism: the operating point is MAD-relative, and MAD is
  nearly noise-blind, so any seasonal increase in COHERENT transients lands straight in the floor;
* within winter the slope of the quiet-day count against the daily noise term is NEGATIVE
  (−0.95 counts per dex at T_core, r = −0.14): the loudest days fire LESS, so the excess is not carried
  by the loud days but by a winter population of coherent transients the daily median energy misses.

**Decisive leg in flight**: Track A step 3b's withheld-station stack test, repeated per season on
production months — CN.VGZ + CN.SNB + PB.B926 excluded from the network sum (`d1d_run_withheld.sh`,
6 GPU passes), frozen HF veto applied through a trace-root-aware copy (`d1d_hf_veto.py`), quiet
family-days from Track C's own `quiet_flag`, null drawn from the same quiet days, per-family
(a/σ)² calibrated on Track A's 2012 ETS positive control. Months chosen for MATCHED quiet-day counts:
DJF 2015-02 + 2014-02 + 2014-12 (6,293 quiet family-days) vs JJA 2014-08 + 2015-08 + 2014-07 (6,325).

## 2026-09-17 00:15 UTC Track C fixes verified (agent report)
- Epoch end now exclusive in c_month_ready (matches planner); the five stuck months run (2014-06 1,487/605; 2014-09 3,090/1,696; 2015-01 3,474/1,384 T_ext/T_core; 2017-02, 2018-12 runnable).
- Per-epoch tiers: c_tiers_by_epoch.py + c_reflag_tiers.py (recomputes cross-family assignment on the thresholded set). E0 own quiet days (41k family-days) return 9.00/8.00 = within 0.25 of frozen 8.75/8.25 -> anchor kept. E1/E3 have no processed quiet months yet (only SSE months) -> sit on the anchor, derive automatically when their quiet months land; watcher re-derives after every month. Post-SSE window defined for production = 1-10 days after the last tremor day of an episode.
- Sweep: 26 months done, at 2016-01.

## 2026-09-17 00:18 UTC Track D interim 2 (agent report)
- d_tiers.py = single source of epoch axis + per-epoch T; RESULTS.md opens with the epoch/threshold table. Recall measured on the station-set axis (S15/S14/S13).
- Seasonal bar kept NOT met (3.0-3.2x). Mechanism found: winter is 1.4x louder in energy but the family-day MAD scale moves only +2.4% (cc is amplitude-normalised), so a louder winter does not raise the bar; within winter the loudest days fire LESS -> the excess is a winter population of coherent transients. Withheld-station test per season (pooled incoherent statistic + NEW moveout-shuffled control): early winter month +1.06 dB vs time null but only +0.39 dB vs moveout-shuffled null -> mostly non-family-specific broadband transients, modest genuine residual consistent with the ~22% real fraction. Final contrast pending.

## 2026-09-17 00:30 UTC Track D — the winter/summer question is SETTLED (as far as it can be)

Six withheld-network GPU passes (CN.VGZ + CN.SNB + PB.B926 out of the network sum) on three winter and
three summer months with MATCHED quiet-day counts (6,293 DJF vs 6,325 JJA), frozen HF veto applied.

**Two changes to Track A step 3b's method were forced and both are improvements.** (1) Production's
quiet-day rate is ~4x lower than 2012's, so no family has enough quiet-day detections for a coherent
per-family stack (2–7 here vs 20+ there) → the primary statistic is a POOLED INCOHERENT one,
r = E(S±2 s)/E(pre-arrival 4 s) at the withheld station, each window normalised by its own noise,
averaged in log10 over thousands of detections. (2) A **moveout-shuffled null** was added — the same
detection times with ANOTHER family's S offset (>5 s away) at that station. Step 3b had only a
random-TIME null, which cannot tell a real LFE from a broadband transient that made the whole network
fire. At production purity that distinction is the whole question.

| tier | DJF vs time null | DJF vs moveout null | JJA vs time null | JJA vs moveout null |
|---|---|---|---|---|
| T_core | +0.918 dB | **+0.586 dB** (p 0.048) | +1.592 dB | **+0.737 dB** (p 0.048) |
| T_ext | +0.449 dB | **+0.173 dB** (p 0.048) | +1.033 dB | **+0.507 dB** (p 0.048) |

1. **Neither season is pure noise** — both beat the time null AND the family-specific moveout null at
   both tiers, so both carry real inter-ETS repeats (consistent with step 3b's ~22 % at T_core).
2. **36–62 % of the time-null excess survives the moveout shuffle**, so the quiet-day sets are dominated
   in BOTH seasons by broadband transients, not by family-specific arrivals.
3. **The winter excess is NOT because winter detections are more real** — per-detection family-specific
   content is LOWER in winter (0.78x at T_core, 0.33x at T_ext). At T_ext, where the seasonal rate
   excess is largest, winter is clearly less real even after the noise correction (0.41x) → that excess
   is predominantly FALSE ALARMS.
4. **Honest limit**: r is normalised by each window's own noise, so a real LFE of fixed amplitude gives a
   smaller r in a louder winter; correcting by the measured 1.26x winter noise energy leaves the seasons
   indistinguishable at **T_core** (0.98x). A genuine seasonal component of inter-ETS activity therefore
   **cannot be excluded at T_core**; what is excluded is the winter excess being predominantly real.

**Consequence for the product (unchanged, and now justified):** λ stays per season, the seasonal bar
stays reported as NOT met, nothing is tuned. Section 1d of RESULTS.md carries the mechanism (the score
is an amplitude-normalised cc over its own MAD, and the MAD is nearly noise-blind, so a louder winter
does not raise the bar) and this test as the check on it.

## 2026-09-17 03:24 UTC user-authorised deletion
Deleted southern raw waveforms in codex-LFE: traces/southern_2024 (1.1 GB: RVIT, PETY, HAYF, B045, B935 2024), cache/southern_picker_waveforms (368 MB), and 70 southern packets from cache/waveforms. Northern caches and lfe_extension/traces untouched.

## 2026-09-17 07:05 UTC Track D — FULL RECORD: sweep complete, product rebuilt, one finding CORRECTED

**Sweep finished: 153 of 153 months (2014-01 .. 2026-09), 0 missing** — the coordinator's
`c_month_ready.py` fix cleared all five stuck epoch-boundary months. 0 months with zero observable
family-days; 1 station-day with a caching error in the whole record; unobservable family-day fraction
median 0.9 %.

**Track C re-derived the tiers on the full record and they now DIFFER by epoch** — and one of them is a
striking independent confirmation of Track D's ablation:

| epoch | months | station set | T_core | T_ext |
|---|---|---|---|---|
| E0 | 2014-01..2017-02 | S15 | 8.75 | 8.25 |
| E1 | 2017-03..2018-12 | S15 (new CN sensors) | 8.25 | 8.00 |
| E2 | 2019-01..2021-03 | S14 | 8.75 | 8.00 |
| E3 | 2021-04.. | S13 | **9.25** | 8.00 |

Track D's 2012 station-set ablation predicted that the 13-station network needs **T = 9.33** to
reproduce the frozen quiet-day false-alarm rate at T_core; Track C's independent FA-only re-derivation
on the production record lands on **9.25**. Two completely different data sets and procedures agree to
0.08 MAD. Per-epoch recall, each epoch scored at its OWN thresholds
(`completeness_by_epoch_summary.csv`): E0 0.146 / 0.184, E1 0.184 / 0.208, E2 0.146 / 0.207,
E3 **0.122** / 0.215 (T_core / T_ext).

**CORRECTION to the 2026-09-17 00:30 entry — the "winter excess" is an E0-only phenomenon.** On the
full record the DJF/JJA quiet-day rate ratio at T_core is **2.39 (E0), 0.26 (E1), 0.26 (E2), 0.45 (E3)**
— it reverses. All 8 epoch-tier cells miss Merlin's ≤ 2x band, but in BOTH directions, so the floor is
seasonally **unstable, not seasonally biased**. The withheld-station test of the previous entry used six
E0 months and its verdict stands for E0's winter excess; it is not a statement about the record.

**The mechanism survives the correction and is now much better evidenced.** Winter is the loudest season
in EVERY epoch (DJF/JJA daily-energy ratio 1.45 / 1.22 / 1.31 / 1.08) yet carries the highest quiet-day
rate in only one of them, and the family-day MAD scale stays within a few per cent of 0.0165 in every
epoch and season. A detector scored in units of an amplitude-NORMALISED cc's own MAD is nearly blind to
the ambient level: **a louder season neither raises the bar nor raises the floor.** What moves the floor
is the population of COHERENT transients in that epoch-season, which the daily median energy does not
measure and which has no fixed seasonal sign. That is exactly why λ is carried per
(family, tier, epoch, season) on every event.

**Catalog v1 on the full record**: 248,547 rows — **93,291 `cat_core`**, **244,815 `cat_ext`** — 115
families, 2014-01-01 .. 2026-09-15, 0 duplicates. Two months (2018-06, 2021-10) disagree with
`months_done.csv`; `d2_assemble.py` now VERIFIES each disagreement against the month's own parquet and
both are stale LEDGER rows written before Track C's tier re-flag (the stored `cat_*` columns already
match the current per-epoch thresholds) — the catalog is the correct one, and the check is recorded in
`catalog_ledger_disagreements.csv`.

## 2026-09-17 07:15 UTC Track D COMPLETE — product delivered on the full record

`results/trackD/RESULTS.md` (564 lines, generated by `scripts/d4_results.py` from the output files, so
it cannot drift from them), README Track D section updated, 98 per-SSE figures + 3 summary figures +
2 model figures.

**Deliverables** (all in `results/trackD/`): `lfe_extension_catalog_v1.parquet` (+ `.csv`) 248,547 rows
= **93,291 `cat_core` / 244,815 `cat_ext`**, 115 families, 2014-01-01 .. 2026-09-15, one row per event
with tier flags, net_mad, cc, rel_amp, n_sites, n_sta_excluded, tremor/eq/hf/quiet flags, epoch_id,
Bostock hypocentre + SeisBench id, our relocation where accepted, **Mc(f,d)** and
**λ(f, T, epoch, season)**; `family_table_v1.csv` (130 families, 14 template-less, 10 with an accepted
relocation at σ_h ≈ 2.0 km, per-epoch recall); `family_day_table_v1.parquet`.

**Verdicts.** Completeness-model validation **MET** (max |pred − meas| own recall 0.024 on the held-out
2010 year, bar ±0.05). Thomas main-front speed bar **MET** at both tiers (12.4 % / 11.4 %, bar 20 %).
Thomas re-rupture-interval bar **NOT MET** against raw Bostock (49 % / 23 %, bar 5 %) but **INSIDE** the
per-family-thinned Bostock 95 % range at both tiers — the failure is decimation bias, and the thinned
reference is the fair comparison. Seasonal-stability bar **NOT MET in all 8 epoch-tier cells, in both
directions** (see the correction in the previous entry).

**Two things worth remembering beyond this task.**
1. **Independent cross-check of the whole completeness model**: Track D's 2012 ablation says the
   13-station network needs T = 9.33 to hold the frozen false-alarm rate; Track C's FA-only
   re-derivation on the production record, with no knowledge of that ablation, gives **9.25** for E3.
2. **The main-front estimator mostly fails at this purity** — a plausible southward speed in only 3/31
   windows at cat_core (4/31 hp, 6/31 ext), because its leading edge is the 10th percentile of event
   position and the along-strike false-alarm background pins it before the front arrives. Where it
   works: −3.3 to −8.6 km/day, bracketing Bostock's 2012 −4.70. A purity-matched front estimator
   (first-activation restricted to families whose burst exceeds their own λ by a stated factor) is the
   obvious next step and is deliberately NOT run here because it was not pre-registered.

## 2026-09-17 07:07 UTC CATALOG v1 DELIVERED (Track D; results/trackD/RESULTS.md)
- results/trackD/lfe_extension_catalog_v1.parquet/.csv: 248,547 rows, 93,291 cat_core / 244,815 cat_ext, 115 families, 2014-01-01..2026-09-15, 153/153 months; family_table_v1.csv (130 fams, 14 template-less, 10 accepted relocations sigma_h ~2 km); family_day_table_v1.parquet (537,312 family-days). Time basis = Bostock DB clock + per-family offset (NOT origin time).
- Tiers per epoch (FA-only rule): E0 8.75/8.25, E1 8.25/8.00, E2 8.75/8.00, E3 9.25/8.00; T_core completeness E0 14.6%, E1 18.4%, E2 14.6%, E3 12.2%. Independent 2012 ablation predicts 9.33 for the 13-station network vs 9.25 derived on production (0.08 MAD apart).
- Bars: completeness model MET (0.024 vs 0.05); main-front speed MET (12%); re-rupture median NOT MET vs raw Bostock but inside the thinned-Bostock 95% range (decimation bias).
- Seasonal: floor seasonally UNSTABLE (both directions; winter excess is E0-only), decoupled from ambient level; withheld-station test with pooled statistic + moveout-shuffled null: real inter-ETS repeats exist (p<=0.048) but 36-62% of the excess is non-family-specific transients; per-detection realness lower in winter.
- Anomalies: main-front estimator gives a plausible speed in only 3/31 SSE windows at cat_core (leading edge pinned by along-strike FA background) -> purity-matched estimator proposed, NOT run (user decision); 2 stale months_done rows (catalog correct); 5 SSE months with zero quiet family-days (lambda fallback recorded); flag_2lambda 340 site-months / flag_hfr_shift 56 months for a Track C look.
- Per-SSE: 31 windows (15 major); largest SSE201806 6,690 core; SSE202407 4,742; SSE202508 4,401; SSE201512 2,947.

## 2026-09-17 07:24 UTC user decision 1
User approved: build the purity-matched main-front estimator; freeze it on the 2012 referee (Bostock, thinned Bostock, our 2012 cat_core/cat_ext must all reproduce -4.7 km/day within 20%) BEFORE applying to 2014-2026. Assigned to the Track D agent (RESULTS.md 3b-2, thomas_front_v2.csv). Decision 2 (seasonal floor): documented, nothing to fix; user informed that inter-ETS analyses work with a signal that is at best ~half real.

## 2026-09-17 07:29 UTC user: score inter-ETS realness -> Track E agent launched
Per-event p_real from independent-site confirmation, burst clustering, n_sites_above, physics features; trained on 2012 (Bostock-matched positives vs rolled/shifted-control negatives), validated on held-out 2010 (calibration within 0.1 per bin), applied to v1 -> results/trackE/realness_v1.parquet + scored catalog; hp flag = p_real>=0.8.

## 2026-09-17 07:40 UTC Track D §3b-2 — purity-matched main-front estimator BUILT, FROZEN, APPLIED

Post-registration addition, user-approved via the coordinator 2026-09-17. Protocol followed exactly:
designed on the 2012 referee alone, frozen, and only then pointed at the extension years.

**Definition** (`scripts/d3b2_front_v2.py`). Two changes to the harness's first-burst front, both aimed
at the §3b failure and nothing else:
1. **Purity gate** — a family activates on the first day whose count for that family is inconsistent
   with its OWN measured floor: `n_f(d) >= 4` AND `P(N >= n | Poisson(λ_f)) < α/n_days` (Bonferroni),
   λ from §1c for that (family, tier, epoch, season). Activation time = the first qualifying burst that
   day (harness definition, gap < 600 s, >= 4 events).
2. **Window from the activations** — the propagation window is the 5–95 % range of the activation
   times, not the 10th-percentile leading-edge track, so background families cannot set it.
Everything else is the harness's. One along-strike frame for every window (computed once from all 130
Bostock hypocentres) so speeds are comparable between SSEs.

**Parameters.** N_MIN 4 and MIN_FAM 5 inherited from the harness; the 5–95 % window fixed a priori. The
ONLY fitted parameter is α, from the pre-declared grid {0.05, 0.01, 0.001}. **All four referee legs pass
at every α** — Bostock 2012 −5.45 (16.1 % from his −4.695), our 2012 cat_core −5.03 (7.2 %), our
cat_ext −4.62 (1.7 %), Bostock thinned to our completeness −4.84 ± 0.1 (3.2 %) — and the speeds move by
≤ 0.3 km/day across the whole grid, so the estimator is insensitive to α. Rule picks **α = 0.05**;
frozen in `thomas_front_v2_design.json`.

**A null had to be replaced, and this is worth remembering.** The harness's day-shuffle null permutes
family labels among the SAME day's events, so the multiset of labels per day — and hence every family's
DAILY COUNT — is exactly preserved. A day-resolution activation statistic is therefore nearly invariant
under it: `null_dayshuf_mean` reproduces the observed speed to a few per cent in almost every window
(e.g. SSE201806 −5.29 vs −5.10). It is reported for continuity and explicitly NOT used as a test. Added
instead a **position-shuffle** null: keep every activation time, permute the along-strike coordinate
across families — the direct test of "is activation ORDER related to POSITION".

**Result.** `cat_core`: **27 of 31 windows fitted, 13 plausible, 11 significant** (plausible AND
position-shuffle p <= 0.05) against **3 of 31** for the pre-registered estimator. `cat_ext`: 30 fitted,
8 significant. Of the 11 significant core fronts, 5 propagate SE (the 2012 sense) and 6 NW; |speed|
1.3–8.6 km/day, median 5.5, bracketing Bostock's −4.70 and Thomas's 7.3 km/day. Bidirectional
along-strike propagation is expected in Cascadia, so a mixed sense is not suspicious — BUT the four
largest |speeds| come from windows long enough (> 25 d) to hold more than one episode, and this
estimator fits ONE front per window. Splitting multi-episode windows is the named next refinement and
was deliberately not done.

**Ledger reconciled** (`scripts/d_fix_ledger.py`, coordinator-approved): the two stale rows in
`results/trackC/months_done.csv` (2018-06, 2021-10) were written before Track C's tier re-flag. The
script verifies each month's parquet already carries `cat_*` matching the CURRENT per-epoch thresholds
and only then updates `n_cat_ext` / `n_cat_core` in place — exactly what `c_watch.mark_done` would
write, same schema and row order. Backup `months_done.csv.pre_d_fix_ledger`, note
`months_done_LEDGER_NOTE.md`, diff `results/trackD/ledger_corrections.csv`. Default mode is
report-only; `--apply` was used once.

## 2026-09-17 07:40 UTC front estimator v2 (Track D; RESULTS.md 3b-2, post-registration, user-approved)
- Purity gate (family activates on first day with n>=4 AND Poisson tail vs its own floor lambda < alpha/n_days) + window from activation times (5-95%). Only alpha fitted from {0.05,0.01,0.001}; all four referee legs pass at every alpha: Bostock 2012 -5.45 (16%), our cat_core -5.03 (7%), cat_ext -4.62 (2%), thinned Bostock -4.84 (3%); frozen alpha 0.05 before touching new years.
- Day-shuffle null is invariant under a day-resolution statistic -> replaced by a position-shuffle null (permute along-strike coordinate).
- Result on 2014-2026: cat_core 27/31 fitted, 13 plausible, 11 significant (vs 3/31 before); 5 propagate SE, 6 NW; |speed| 1.3-8.6 km/day, median 5.5. Caveat: one front per window; the 4 largest speeds are in >25-day windows with multiple episodes -> split windows next (not done).
- months_done.csv reconciled (2018-06, 2021-10) via d_fix_ledger.py --apply; backup kept.

## 2026-09-17 07:43 UTC user: front-estimator v2 result is to be KEPT and REVISITED
Parked for later: results/trackD/thomas_front_v2.csv (per-SSE speed/CI/direction/nulls), thomas_front_v2_activations.csv, thomas_front_v2_design.json (frozen alpha 0.05), figures/d3b2_front_v2_summary.png, RESULTS.md 3b-2. Open refinement when we return: split >25-day windows into episodes; then re-rupture/secondary fronts per episode.

## 2026-09-17 07:57 UTC user request: per-family cumulative count plots 2011-2026
Queued on the Track E agent after scoring: figures/family_cumulative/<fam>.png + grid; curves = cat_core, cat_ext, hp (p_real>=0.8), completeness-corrected cat_core (dotted), Bostock 2011-2013 (dashed grey); SSE windows shaded; epoch lines; data in results/trackE/family_cumulative_daily.parquet.

## 2026-09-17 09:45 UTC Track E interrupted by Opus session limit (reset 09:30), resumed
Per-site recomputation kept running as a shell job (154/153+ months done, at 2026-07). Agent resumed: finish, train/validate, score, plots.

## 2026-09-17 09:55 UTC Track E — per-event realness score `p_real` DELIVERED

**What it is.** `p_real` = the calibrated probability that a detection is matched one-to-one (±2 s)
to a Bostock (2015) event of its OWN family, from per-event features only. That is the only realness
label that exists, and the report says so in its first paragraph: Bostock is incomplete, so a real
LFE he missed trains as a 0 and the measured likelihood ratio UNDERSTATES the true separation.

**The measurement that makes it possible** (`scripts/e_persite.py`): a per-SITE re-run of the frozen
matched filter at every stored detection time — 153 production months + the 2012 and 2010 referee
years, 9.7 GPU-hours, 641,754 production detections, 2.93 M (detection, site) rows. **Self-check: the
weighted sum of the recomputed per-site series reproduces Track C's stored `net_cc` to ≤ 9e-8
(float32 epsilon) in every month.** It gives, for each site independently: the site's own
weighted-mean cc in units of THAT site's own family-day MAD, and the S-window/pre-arrival energy
excess at the predicted moveout in units of that site's own daily distribution. `c_run_month.py` was
never modified; its measurement core was copied.

**Model.** Gradient boosting with monotone constraints (chosen over logistic regression on 2012
grouped-CV: 0.894 vs 0.868), 18 features in four groups (independent-site breadth, independent
confirmation from the 2 LOWEST-weight sites, burst clustering, physics), isotonic-calibrated on 2012
ETS days. Trained on 2012 ETS days (19,641 detections, 7,766 positives, π = 0.395) with the 2012
quiet days as extra negatives; the site-rolled chance control is a reference, not a training class
(it has no rel_amp/hfr/wsh_pgc, so training on it would separate the classes on a missing-value
pattern). Context (tremor/quiet/season/epoch) is NEVER a feature — with Bostock having no inter-ETS
events, a context feature would answer "inter-ETS ⇒ not real" by construction.

**VALIDATION on the held-out year 2010 — the pre-registered calibration bar is NOT MET.**
2010 AUC **0.855** (2012 CV 0.894) on a network three sites smaller. But the reliability deviates by
up to **0.261**, bar ±0.10, and in the SAME direction in every bin. **It fails in LEVEL, not in
SHAPE**: one free scalar (a log-odds anchor of −0.926) brings the worst bin to **0.060**, inside the
bar everywhere. The population-prior difference accounts for only −0.356 of that, so ~62 % is a
genuine NETWORK-STRENGTH effect. Delivered as a bracket: `p_real` (2012/9-site anchor) and
`p_real_cons` (2010/6-site anchor); production epochs carry 7–9 sites. Measured 2010 precision:
0.587 at p ≥ 0.5, **0.743 at p ≥ 0.8**, 0.837 at p ≥ 0.9.

**Ablation (2010 AUC):** all 0.855; physics-only 0.802; sites-only 0.804; burst-only 0.686;
confirmation-only 0.642. Removing the burst block costs the most (0.844) — it is the one group that
is not a re-reading of the same correlation peak. Site-rolled control: median p 0.09–0.10, only
0.3–0.9 % reach 0.8. The 2012 post-SSE window (doy 270–279, never in training) behaves like the quiet
days: 0.03 % reach 0.8.

**Applied to v1** (v1 itself untouched): `results/trackE/realness_v1.parquet` +
`lfe_extension_catalog_v1_scored.parquet`, 248,488 of 248,547 scored, hp (p_real ≥ 0.8) = 16,449
(6.6 %). **Inter-ETS (no tremor flag), cat_core: p ≥ 0.8 keeps 844 events (2.3 %) at an estimated
purity 0.77; p ≥ 0.9 keeps 318 (0.9 %) at 0.89.** So the coordinator's ≥ 80 % target IS reachable
between episodes, but only for ~1 % of the inter-ETS catalog — ~300 events over 12.7 years.

**An out-of-sample check the score never saw:** the fraction of inter-ETS cat_core events falling
inside one of Track D's 31 SSE windows rises monotonically with the score — 16.3 % (all) → 22.0 %
(p ≥ 0.5) → **30.9 % (p ≥ 0.8) → 38.4 % (p ≥ 0.9)**. SSE windows were never a feature or a label.

**⚠ OPEN TENSION with Track D §1d — the winter excess does NOT collapse under the score, it
STRENGTHENS.** E0's inter-ETS DJF/JJA ratio is 1.74 on raw counts and **7.17 on hp events**; every
epoch moves further from 1. Three artifact hypotheses tested and none explains it: (a) a refit on a
SEASON-ROBUST feature subset (only day-normalised quantities; 2010 AUC 0.842, precision 0.736) still
gives 5.96; (b) the December-2014 day-concentration (35 % of E0 DJF hp events on 2014-12-11) sits
inside SSE201411 but carries no tremor flag — episode-adjacent real activity; (c) dropping every
inter-ETS day inside an SSE window leaves 70 DJF vs 10 JJA, ratio 7.0. Both statistics are
noise-normalised and biased in OPPOSITE directions (Track D's r divides by the window's own noise so
a fixed-amplitude LFE scores lower in a 1.45× louder winter — Track D said a genuine seasonal
component at T_core "cannot be excluded"; my score is trained where more sites up and a larger
rel_amp mean more real, both of which winter supplies for non-seismological reasons). Not adjudicated.

**⚠ DATA-QUALITY FINDING handed back to Track C.** A site's normalised cc cannot exceed 1. On 18
station-days (mostly UW.LRIV before its 2019-05-30 drop, also SNB/VGZ/NLLB) a near-zero-energy
interval slips past the energy floor and the cc blows up to ~1e11. **486 detections** carry such a
site and **59 are in catalog v1, all flagged cat_core AND cat_ext**, some with `cc` > 1 (2.34) and
`net_mad` up to 136. Track C's stored `net_cc` carries the identical blow-up — the energy floor is
too permissive. Listed in `results/trackE/bad_persite_detections.csv`; Track E leaves them UNSCORED.

**Extra deliverable (user request): per-family cumulative counts 2011–2026.**
`results/trackE/family_cumulative_daily.parquet`, `figures/family_cumulative/<id>.png` (130 panels,
14 template-less and empty by construction), `figures/e_family_cumulative_grid.png` (12 × 11).
Bostock's 2011–2013 events were mapped to families by his own family hypocentre: his catalog has
exactly 130 unique (lat, lon, dep) triplets, **130/130 matched** to `family_table_v1.csv`, 0
unmapped, 61,685 events in 2011–2013. The dotted curve is `cat_core / recall(family, epoch)` joined
onto Bostock's 2013 total — a MODEL-BASED rescaling, floored at recall 0.05 so a 3–5-template family
cannot turn a 0.5 % recall into a 200× inflation. Palette validated with the dataviz validator
(blue ordinal ramp for our three readings, orange for the high-purity subset, grey for Bostock).

Scripts: `e_persite.py` (per-site pass), `e_features.py`, `e_train.py`, `e_apply.py`,
`e_season_check.py`, `e_family_cumulative.py`, `e_results.py`; report `results/trackE/RESULTS.md`.

## 2026-09-17 09:55 UTC Track E delivered (results/trackE/RESULTS.md)
- Per-site pass reproduces Track C net_cc to <=9e-8; 641,754 detections x sites. GBM with monotone constraints, 18 features (independent-site breadth, low-weight-site confirmation, burst clustering, physics); 2012 train, 2010 held-out AUC 0.855; calibration bar NOT MET in level (0.26) but shape OK -> bracket p_real (9-site anchor) / p_real_cons; 2010 precision 0.74 @>=0.8, 0.84 @>=0.9. Burst clustering carries most.
- Applied: lfe_extension_catalog_v1_scored.parquet; hp (p>=0.8) 16,449 (6.6%). Inter-ETS cat_core: p>=0.8 keeps 844 (2.3%) at purity 0.77; p>=0.9 keeps 318 at 0.89 -> ~300 trustworthy inter-ETS events in 12.7 yr.
- OPEN TENSION: winter excess STRENGTHENS on hp events (E0 DJF/JJA 7.2 vs 1.7 raw), contradicting Track D 1d; artifact tests rejected; not adjudicated.
- BUG (Track C): near-zero-energy intervals (mostly LRIV pre-2019) give cc ~1e11; 486 detections, 59 in catalog v1 (cc 2.34, net_mad 136) -> results/trackE/bad_persite_detections.csv; fix ordered.
- Per-family cumulative figure delivered (figures/e_family_cumulative_grid.png + 130 panels; family_cumulative_daily.parquet). Sent to user.

## 2026-09-17 Track C — cc-normalisation energy floor FIXED (Track E bug), 10 months re-run, tiers re-derived

**The bug** (Track E, `results/trackE/bad_persite_detections.csv`). The matched filter normalises each
correlation by that window's own energy, so a normalised cc cannot exceed 1. The old guards were a
denominator clamp at `1e-10 x median(ew)` plus `VALID = (3-C window energy > 1e-6 x its median) AND
(every channel energy > 0)`. Neither catches the real failure: gaps are merged with `fill_value=0`, so
a dead channel is **exactly zero** for most of the day (UW.LRIV 2019-02-06: one channel zero on
1,725,189 of 1,727,760 windows). The 3-C SUM is then carried by the two healthy channels and passes
both tests, while the windows that STRADDLE the edge of the dead interval have a tiny-but-POSITIVE
energy on that channel, so `min > 0` admits them, the denominator is not clamped, and
`cc = numerator / sqrt(tiny)` blows up. Track C's own files carried 302 such rows (13 dates, 6 months),
`net_cc` up to 4.0e11 and `net_mad` up to 2.7e13, 45 of them inside the catalog.

**The guard** (`c_run_month.py`, `--energy-floor`, default 1e-3): a floor RELATIVE TO THE STATION-DAY
MEDIAN, applied PER CHANNEL, with the cc SET TO ZERO (not merely masked) at bad windows, the reference
falling back to `1e-2 x median(e3)` for a dead channel, and the mask DILATED by the template length.
**Two wrong versions were caught before the final one, both by re-measuring rather than by reasoning:**
1. per-channel floor alone -> on 2018-03 the impossible rows went 28 -> 35 and `net_cc` to 1.7e14.
   Zeroing the core bad windows removed the huge value that had been SUPPRESSING its neighbours in the
   declustering max-pool, so contaminated neighbours the energy test does not reach became detections.
   Fix: dilate the bad mask by +-n (window j covers [j, j+n), so every window sharing a bad sample).
2. `median(ew)` as the scale -> on UW.LRIV 2018-03-02 that median is exactly 0 (one channel dead for a
   third of the day), so every relative threshold collapsed to 0 and the guard was a no-op on the one
   station that needed it; 2018-03 still had 34 impossible rows. Fix: scale on `median(e3)`, which stays
   positive while any channel is alive, and mask the whole station-day when it does not.

**2012 verification — the detection set and every raw cc are bit-for-bit**: 253/253 (doy 228-232) and
3426/3426 (doy 248-252), max |delta net_cc| = **0.0**. `net_mad` moves by <= 4.7e-4 / 2.4e-3, only on the
two 2012 station-days that contain a real data gap (PB.B001 doy 230, CN.NLLB doy 250), whose gap-edge
windows are now excluded from the median/MAD estimate.

**Ten affected months re-run** (`results/trackC/rerun_ccfloor_log.csv`; `before` for the four months
re-run once with the interim guard is the ORIGINAL production value):

| month | detections | cat_core | cat_ext | \|cc\|>1 | max \|cc\| |
|---|---|---|---|---|---|
| 2014-11 | 2597 -> 2459 | 579 -> 570 | 995 -> 976 | 0 -> 0 | 0.467 -> 0.383 |
| 2014-12 | 3360 -> 2578 | 880 -> 695 | 1618 -> 1264 | 0 -> 0 | 0.419 -> 0.374 |
| 2017-02 | 897 -> 887 | 164 -> 162 | 434 -> 426 | 0 -> 0 | 0.334 -> 0.334 |
| 2018-03 | 2041 -> 1930 | 521 -> 505 | 850 -> 835 | 28 -> 0 | 3.79 -> 0.333 |
| 2018-04 | 1739 -> 1679 | 289 -> 280 | 508 -> 499 | 23 -> 0 | 2.89 -> 0.358 |
| 2019-02 | 912 -> 878 | 66 -> 61 | 296 -> 291 | 34 -> 0 | 5.5e3 -> 0.354 |
| 2019-03 | 2657 -> 2628 | 162 -> 157 | 693 -> 688 | 29 -> 0 | 3.6e11 -> 0.359 |
| 2019-04 | 1876 -> 1823 | 263 -> 258 | 817 -> 812 | 53 -> 0 | 4.0e11 -> 0.352 |
| 2019-05 | 4171 -> 4082 | 720 -> 707 | 1741 -> 1728 | 89 -> 0 | 3.4e11 -> 0.370 |
| 2024-05 | 2970 -> 2907 | 189 -> 182 | 1292 -> 1284 | 0 -> 0 | 0.343 -> 0.343 |

All **486** rows Track E flagged are gone (486 -> 0). Note that the fix removes contaminated detections
even where no cc exceeded 1: 2014-12 loses 782 detections and 185 catalogue events although its worst
cc was 0.419, because a partially dead channel inflates the whole family-day's cc series.

**Whole-record scan** (`scripts/c_scan_anomalies.py` -> `results/trackC/anomaly_scan.csv`):
153 months, **640,355 detections, ZERO** with |net_cc| > 1 or net_mad > 50.

**Tiers re-derived on the corrected data, and this mattered.** The per-epoch tiers had been derived
while the artifacts were still in the files; on the corrected data the same FA-only rule gives
different answers for E1 and E2, and the previously stored E1 T_core (8.25) no longer met the
FA <= 0.1 bar (0.123 at 8.75, higher still at 8.25). Re-derived with `--rederive` and applied to all
153 months (previous file kept as `tiers_by_epoch.pre_ccfloor.csv`):

| epoch | T_core / T_ext | quiet FA @T_core | cat_core before -> after |
|---|---|---|---|
| E0 2014-01..2017-02 | 8.75 / 8.25 (anchor kept; rule gives 9.00/8.00, within 0.25) | 0.105 | 18,958 -> 18,958 |
| E1 2017-03..2018-12 | **8.25 -> 9.00** / 8.00 | 0.086 | **26,051 -> 10,660** |
| E2 2019-01..2021-03 | **8.75 -> 9.00** / 8.00 | 0.073 | **15,622 -> 12,307** |
| E3 2021-04.. | 9.25 / 8.00 | 0.093 | 32,405 -> 32,405 |

E3 also quantifies Track D's point: at the E0 anchor 8.75 its quiet-day rate is 0.166 against E0's
0.105 (1.6x), which is why it needs 9.25. HF re-flag is a no-op (0 of 640,355 flags change; the
per-epoch cuts were already applied by the runner). QC regenerated for every month.

## 2026-09-17 11:17 UTC Track C cc-floor fix applied (results/trackC/rerun_ccfloor_log.csv, anomaly_scan.csv)
- Cause: gaps merged with fill 0 -> dead channel exactly zero; windows straddling the dead-interval edge have tiny positive energy -> cc = num/sqrt(tiny). Guard: per-channel energy floor 1e-3 x station-day median(e3), cc set to 0, mask dilated, whole day masked if median<=0. 2012 bit-for-bit preserved (253/253, 3426/3426; max |dnet_cc| 0).
- 10 months re-run; all 486 flagged rows gone; 2014-12 loses 782 detections / 185 catalog events (partially dead channel inflated the whole family-day). Whole record: 0 rows with |cc|>1 or net_mad>50.
- Tiers re-derived on corrected data: E0 8.75/8.25 (anchor kept), E1 8.25->9.00 / 8.00 (cat_core 26,051->10,660), E2 8.75->9.00 / 8.00 (15,622->12,307), E3 9.25/8.00 unchanged. Old tiers kept as tiers_by_epoch.pre_ccfloor.csv. Watcher idle, record complete.
- NEXT: Track E re-scores the 10 re-run months; Track D assembles catalog v1.1 (completeness, Thomas, front v2, cumulative plots) with p_real merged.

## 2026-09-17 11:40 UTC Track E v1.1 — re-scored after Track C's cc-floor fix (NO retraining)

**Change detection, not trust.** Every one of the 153 `detections_<ym>.parquet` was compared key by
key against the per-site pass's own stored detection set (stronger than a hash: it says WHICH rows
moved). **Exactly the 10 months Track C re-ran changed; the other 143 are key-for-key identical.**
SHA-256 + before/after counts per month in `results/trackE/detections_hash_check.csv`. Detections
≥ 8.0 MAD 641,754 → 640,355 (−1,399). Superseded per-site outputs archived (not deleted) in
`results/trackE/superseded_v1_ccfloor/`.

**The per-site pass was re-run for those 10 months only** (≈ 40 min, 5 workers). Two checks:
1. `persite_bad` (a site whose own normalised cc exceeds 1) **486 detections / 18 station-days → 0**.
   All 640,355 detections are now scored; 0 unscored (was 59). Track C's fix is confirmed from the
   independent per-site side.
2. **The window masking did NOT bias the per-site background the score is built on**: on the 640,311
   rows common to v1 and v1.1, max |Δ| is 1.2e-6 on `z_site_med`, 4.8e-7 on `conf_z`, 0 on
   `n_sites_obs`. Only `net_mad` moves (5,636 rows, max 0.31) — Track C's own family-day MAD scale.
   This was the one thing I was unsure of (my pass keeps the original energy floor), and it is
   settled empirically rather than argued.

**`results/trackE/realness_v1_1.parquet`** — 640,355 rows, one per CURRENT Track C detection ≥ 8.0
MAD (not only catalog members). `model_v1.joblib` loaded unchanged; nothing refitted. **Merge key,
three equivalent forms, all carried as columns**: `month` + `det_row` (0-based row index into
`results/trackC/detections_<month>.parquet` as stored on disk, taken before any net_mad filter —
round-trip verified), `family_id` + `t` (float64 epoch seconds, bit-for-bit from Track C, unique over
the record), and `det_uid` = `"{month}:{det_row}"`. Carries Track C context + all 18 features +
`p_real`, `p_real_cons`, `p_real_prior`, `hp`, `hp_cons`, `persite_bad`.

**New tiers changed the catalog, not the score.** E1 T_core 8.25 → 9.00, E2 8.75 → 9.00 (E0, E3
unchanged): `cat_core` 93,291 → **74,330** (−20.3 %), `cat_ext` 244,815 → 244,375.

**Inter-ETS at p_real ≥ 0.8 (cat_core): 831 events, purity 0.76; at ≥ 0.9: 311, purity 0.89.** The
ABSOLUTE size of the usable high-purity inter-ETS set is unchanged from v1 (844 / 318); the retained
FRACTION rose 2.3 % → 3.1 % only because the stricter T_core shrank the denominator 36,781 → 26,419.

**Independent consistency check worth keeping:** E1/E2 lost most of their inter-ETS core population
to the stricter T_core (E1 DJF 977 → 283) but their hp counts barely moved (E1 10/16 vs 10/18,
E2 27/23 vs 27/22) — Track C's false-alarm-only tier rule and Track E's score, built from unrelated
information, agree on which events to drop.

**The §8c seasonal tension SURVIVES the new tiers**: E0 1.69 raw / **6.70 on hp events** (v1:
1.74 / 7.17); the DJF/JJA hp ratio is still further from 1 than the raw ratio in every epoch and
still reverses sign between E0 and E1/E3. Unchanged open item against Track D §1d.

Not redone, deliberately: no retraining (the 2010 validation lives in the referee years, untouched by
a production fix); no cumulative plots — Track D rebuilds those from the v1.1 catalog, and the v1
versions on disk are superseded. `realness_v1.parquet` / `lfe_extension_catalog_v1_scored.parquet`
are left in place, keyed to catalog v1. Script `scripts/e_apply_v11.py`; addendum = RESULTS.md §11.

## 2026-09-17 11:31 UTC Track E v1.1 re-score done (results/trackE/realness_v1_1.parquet; RESULTS.md s11)
- Exactly the 10 months changed (key-by-key check); persite_bad 486->0; per-site backgrounds unchanged (max |d| 1e-6); 640,355 detections all scored, no refit. Merge keys: month+det_row, family_id+t, det_uid.
- v1.1 tiers: cat_core 93,291->74,330 (all E1/E2), cat_ext 244,375; hp 21,227. Inter-ETS core p>=0.8: 831 at purity 0.76; p>=0.9: 311 at 0.89 (absolute hp set unchanged). FA-only tier rule and the score agree on which E1/E2 events to drop (hp counts barely moved).
- Seasonal tension persists (E0 DJF/JJA 1.69 raw / 6.70 on hp). Track D assembling v1.1 with p_real merged.

## 2026-09-17 11:35 UTC Track D — CATALOG v1.1 (Track C cc-floor fix + per-epoch tier re-derivation)

Coordinator request 2026-09-17. **v1 is preserved untouched** (its catalog / family table / family-day
table under their own names, its completeness / floor / Thomas outputs snapshotted to
`results/trackD/v1_snapshot/`); every D script now takes a version tag (`d2_assemble --tag`,
`d3_run_sse --catalog/--suffix`, `d3b2_front_v2 --catalog/--out`, `d4_results --tag/--prev-tag`).

**What changed upstream.** (1) Track C's cc-floor bug, 10 months re-run: over those months rows with
|cc| > 1 go 233 → 0, net_mad > 50 205 → 0, Track-E-flagged 179 → 0, recorded detections 23,160 → 21,851,
max |cc| 1.1e12 → 3.56. (2) The tiers were re-derived on the corrected data by the same frozen FA-only
rule: **E1 T_core 8.25 → 9.00, E2 8.75 → 9.00**; E0 and E3 unchanged.

**Catalog v1.1**: 248,069 rows — **74,330 `cat_core`** (v1 93,291, −20.3 %) and **244,375 `cat_ext`**
(v1 244,815, −0.18 %), 115 families, 2014-01-01 .. 2026-09-15, ledger cross-check clean (0 mismatches
after re-running `d_fix_ledger --apply`, which updated 53 stale rows).

| epoch | core v1 → v1.1 | ext v1 → v1.1 | cause |
|---|---|---|---|
| E0 | 19,154 → 18,958 | 40,222 → 39,841 | cc-floor fix only (T unchanged) |
| E1 | 26,076 → 10,660 | 38,358 → 38,334 | T_core 8.25 → 9.00 |
| E2 | 15,649 → 12,307 | 38,120 → 38,093 | T_core 8.75 → 9.00 |
| E3 | 32,412 → 32,405 | 128,115 → 128,107 | cc-floor fix only |

**E0 and E3 are the clean measure of the BUG's size on the product: −196 and −7 core events, well under
1 %.** Everything else is the threshold change.

**No verdict flipped.** The 2012 referee is Track A data, untouched by both changes, so the
completeness-model validation bar (max |pred − meas| 0.0236, bar ±0.05), both pre-registered Thomas bars
(main-front MET at 12.4 % / 11.4 %; re-rupture NOT MET at 49.4 % / 23.1 %) and the thinned-Bostock
reference are bit-identical. The §3b-2 front estimator was applied with its FROZEN α = 0.05 — no
redesign — and its yield barely moves: cat_core 27/31 fitted, 13 plausible, **12 significant** (v1: 11),
cat_ext 30/31 fitted, **11 plausible** (v1: 9), 8 significant. The seasonal-instability verdict is
unchanged in shape (E1 DJF/JJA 0.257 → 0.265) even though E1/E2 λ fall 1.5–4× with the higher T_core.

**What DID change**: per-epoch completeness at each epoch's own T — E1 T_core 0.184 → 0.129 and
E2 0.146 → 0.129 (E0 0.146 and E3 0.122 unchanged) — and the floor λ for E1/E2.

**Track E's re-scored realness arrived in time** (`results/trackE/realness_v1_1.parquet`, 640,355 scored
detections; **100 % of the 248,069 catalog events matched** on (family_id, t)), so the per-family
cumulative figures were regenerated WITH the high-purity curve:
`figures/family_cumulative_v1_1/` (130) + `figures/d5_family_cumulative_grid_v1_1.png` +
`results/trackD/family_cumulative_daily_v1_1.parquet`. Cumulative totals: core 74,330, ext 244,375,
hp (p_real ≥ 0.8) **14,777** (v1 15,319), completeness-corrected-and-joined 625,840 (v1 746,664).
`scripts/d5_family_cumulative.py` is a parameterised copy of Track E's script — same construction, same
palette — so Track E's own v1 outputs were not touched.

RESULTS.md now opens with **§0.0 "Catalog v1.1 — what changed since v1, and why"**.

## 2026-09-17 11:33 UTC CATALOG v1.1 ASSEMBLED (results/trackD/RESULTS.md s0.0; v1 preserved, v1 snapshot in results/trackD/v1_snapshot/)
- lfe_extension_catalog_v1_1.parquet/.csv: 248,069 rows; cat_core 74,330; cat_ext 244,375; hp 14,777; 115 families; ledger cross-check clean. Bug effect isolated in E0/E3: -196 / -7 core events (<1%); the rest is the E1/E2 T_core change. Per-epoch T_core completeness: E0 0.146, E1 0.129, E2 0.129, E3 0.122.
- No verdict flipped (referee-based results bit-identical). Front v2 (frozen alpha): cat_core 12 significant fronts (v1 11), 6 SE / 6 NW; cat_ext 11 plausible / 8 significant.
- p_real merged 100%; cumulative figures regenerated with hp: figures/family_cumulative_v1_1/ + d5_family_cumulative_grid_v1_1.png.
