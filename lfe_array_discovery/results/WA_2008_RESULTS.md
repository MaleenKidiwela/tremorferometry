# Washington Puget / S-Olympic 2008 — Bostock's published recipe on the last full-array Washington bed
Agent N, 2026-09-21.  **The same recipe and the same code** as the north reproduction and the N-Olympic 2010 run (`results/REPRO_BOSTOCK.md` §0 parameter table): 6-s windows, 1–8 Hz, per-component cc averaged over **all** admitted stations, threshold = median + 8 × MAD of the block's own distribution, |lag| ≥ 6 s, 4-s declustering, connected-component grouping at the 4-s node merge (≥3 members), Bungum–Husebye realignment before every stack, 12-s station-tapered templates, network matched filter at 8 × MAD, **no waveform-similarity merge**.  Engine files `scripts/k_auto.py`, `k_group.py`, `k_grow.py`, `k_locate.py` are agent K's, unedited; the bed is supplied by `scripts/n_common.py` and driven by `n_auto.py` / `n_run.py` / `n_locate.py` / `n_screen.py` / `n_relocate.py` / `n_refire.py`, which are the 2010/2011 `l_`/`m_` scripts with the bed adapter swapped.
**There is no Washington LFE catalogue and no tremor catalogue before 2010-01-06.**  Nothing below is a rediscovery count.  What replaces truth is §1 (a tremor-hour audit built from the waveforms), §4 (a five-clause screen whose every clause is one of the failure modes the 2010→2011 test exposed) and §5 (an out-of-sample re-fire on 2007-01), each with a whole-chain site-shift null run through byte-identical code.

## 1. The tremor-hour audit (there is no 2008 tremor catalogue)
The 2008 hour flag this project had (`results/hour_energy_flag_2008.csv`, network-median hourly energy over a fixed p95) selects 113 hours of which **86 % fall in 15–22 UTC** (`results/WA_V2_INPUTS.md` §2.2) — i.e. it is a cultural-noise detector, and the brief for this run required it to be replaced before anything else happened.  `scripts/n_hours.py` replaces it with three clauses (full statement in that file's docstring):
1. **the station's own diurnal cycle divided out** — E(s,d,h) = median over the hour of the 6-s-window mean square of the 2–8 Hz 3-C waveform, referred to that station's own median at the same hour-of-day over the literature-quiet days;
2. **simultaneity** — ≥ 60 % of the reporting near-field stations (admitted, within 60 km of an active cell, health-veto applied) above **their own** quiet-day p95 at once;
3. **spatial coherence** — the log-excess must fit an amplitude decay from ONE active cell better than a uniform field does, at p ≤ 0.05 against 200 station-label shuffles.

### 1.1 2008 bed
* near-field stations 49; literature-quiet days 19; hours judged 1368.
* clause 2 alone: **120 hours on 20 days**; clause 2 AND 3: **115 hours on 20 days**.
* **SELECTED: 115 hours on 20 days**; fraction in 15–22 UTC **0.009** (the old flag: 0.86); weekday fraction 0.774 (5/7 = 0.714 for a calendar-random set).
* **FALSE-POSITIVE CHECK — on the 456 hours of the literature-quiet days the same statistic fires 0 times on clause 2 and 0 times on clause 2+3.**
* selected days: 2008-05-05, 2008-05-06, 2008-05-07, 2008-05-08, 2008-05-09, 2008-05-10, 2008-05-11, 2008-05-12, 2008-05-14, 2008-05-15, 2008-05-16, 2008-05-17, 2008-05-18, 2008-05-19, 2008-05-20, 2008-05-21, 2008-05-22, 2008-05-23, 2008-05-24, 2008-05-26

UTC-hour histogram of the selected hours: `0:4 1:5 2:7 3:9 4:11 5:11 6:14 7:11 8:9 9:9 10:9 11:9 12:3 13:1 22:1 23:2`

The array's own diurnal NOISE FLOOR (network median of the quiet-day reference M(s,h) by UTC hour, normalised to its minimum): `0:3.7 1:2.4 2:2.0 3:1.5 4:1.5 5:1.5 6:1.2 7:1.1 8:1.0 9:1.0 10:1.2 11:1.7 12:3.1 13:5.5 14:7.2 15:8.8 16:9.3 17:9.3 18:9.8 19:10.8 20:9.4 21:10.3 22:6.7 23:5.1` — it is ~10× higher at 19 UTC than at 09 UTC.  A fixed tremor amplitude therefore produces a ~10× smaller ratio in the cultural afternoon, so **the clock distribution of the selected hours is detection-limited, not a claim that tremor is nocturnal**; the quiet-day false-positive rate by UTC hour (`0:0.04 1:0.042 2:0.042 3:0.042 4:0.042 5:0.049 6:0.06 7:0.06 8:0.086 9:0.078 10:0.048 11:0.07 12:0.079 13:0.06 14:0.054 15:0.052 16:0.055 17:0.053 18:0.042 19:0.027 20:0.041 21:0.034 22:0.051 23:0.047`) is flat, so the THRESHOLD is not what makes the set night-heavy — the SNR is.

### 1.2 2007 bed
* near-field stations 43; literature-quiet days 15; hours judged 1416.
* clause 2 alone: **17 hours on 9 days**; clause 2 AND 3: **15 hours on 8 days**.
* **SELECTED: 16 hours on 9 days**; fraction in 15–22 UTC **0.188** (the old flag: 0.86); weekday fraction 0.625 (5/7 = 0.714 for a calendar-random set).
* ⚠ **fallback applied**: only 15 hours pass clause 2 AND clause 3 at p<=0.05; the relaxed coherence bar p<=0.2 is used and RECORDED
* **FALSE-POSITIVE CHECK — on the 360 hours of the literature-quiet days the same statistic fires 2 times on clause 2 and 2 times on clause 2+3.**
* selected days: 2007-01-06, 2007-01-17, 2007-01-18, 2007-01-19, 2007-01-20, 2007-01-21, 2007-01-22, 2007-01-23, 2007-01-26

UTC-hour histogram of the selected hours: `7:1 8:1 9:1 10:4 11:1 12:1 13:3 14:1 15:1 20:1 21:1`

The array's own diurnal NOISE FLOOR (network median of the quiet-day reference M(s,h) by UTC hour, normalised to its minimum): `0:3.0 1:1.9 2:1.8 3:1.3 4:1.2 5:1.5 6:1.7 7:1.5 8:1.2 9:1.1 10:1.0 11:1.1 12:2.0 13:2.4 14:3.3 15:4.5 16:4.7 17:5.3 18:5.1 19:4.5 20:4.7 21:4.5 22:3.8 23:3.3` — it is ~10× higher at 19 UTC than at 09 UTC.  A fixed tremor amplitude therefore produces a ~10× smaller ratio in the cultural afternoon, so **the clock distribution of the selected hours is detection-limited, not a claim that tremor is nocturnal**; the quiet-day false-positive rate by UTC hour (`0:0.042 1:0.034 2:0.047 3:0.054 4:0.066 5:0.073 6:0.068 7:0.066 8:0.079 9:0.088 10:0.083 11:0.057 12:0.041 13:0.04 14:0.044 15:0.03 16:0.033 17:0.036 18:0.026 19:0.015 20:0.033 21:0.048 22:0.06 23:0.044`) is flat, so the THRESHOLD is not what makes the set night-heavy — the SNR is.

**Against the literature.**  The published Washington ETS of 2008 is the early-to-mid-May episode; the independent array-energy flag of `results/P_STEP2_SITES.md` §3.2 (a different statistic, fixed network p95) puts it on 2008-05-05..05-23 plus 05-29 and 06-10.  The audit here selects hours only on 20 days, all inside 2008-05-05 … 2008-05-26 — the two agree, and the audit reaches the answer with the diurnal cycle removed and a spatial-coherence test that the energy flag does not have.

## 2. The bed

**2008 discovery** — window 2008-04-15 … 2008-06-10 (tranche T4), sub-region Puget/S-Olympic (lat 46.2–47.8).
* **51 three-component stations** of the 52 probe-admitted ones (`results/site_probe_2008.csv`): {'PB': 10, 'XU': 40, 'YW': 1}.
* dropped by the > 30 % vetoed-station-hour bar: YW.FACB.
* literature-quiet days: 19 (2008-02-01, 2008-02-08, 2008-02-14, 2008-02-21, 2008-02-27, 2008-03-05, 2008-03-18, 2008-03-24, 2008-03-31, 2008-07-05, 2008-07-08, 2008-07-11, 2008-07-14, 2008-07-17, 2008-07-19, 2008-07-22, 2008-07-25, 2008-07-28, 2008-07-31).
* active cells in the sub-region: 264 (Slab2 28–45 km with ≥ 30 PNSN **2010–2011** epicentres within 8 km — there is no 2008 catalogue, so these are 2010–11 cells used for a 2008 bed, stated not hidden); 49 admitted stations within 60 km of one.
* intraslab-earthquake S station terms available for 3 stations; 48 carry 0.0 and are FLAGGED, never imputed.

**2007-01 confirmation** — window 2007-01-01 … 2007-02-28 (tranche T5), sub-region Puget/S-Olympic (lat 46.2–47.8).
* **45 three-component stations** of the 52 probe-admitted ones (`results/site_probe_2008.csv`): {'PB': 7, 'XU': 38}.
* absent from this window: PB.B201, PB.B941, PB.B943, XU.BS06, XU.C04A, YW.FACA, YW.FACB.
* literature-quiet days: 15 (2006-11-01, 2006-11-05, 2006-11-09, 2006-11-13, 2006-11-18, 2006-11-22, 2006-11-26, 2006-11-30, 2007-04-01, 2007-04-06, 2007-04-11, 2007-04-15, 2007-04-20, 2007-04-25, 2007-04-30).
* active cells in the sub-region: 264 (Slab2 28–45 km with ≥ 30 PNSN **2010–2011** epicentres within 8 km — there is no 2008 catalogue, so these are 2010–11 cells used for a 2008 bed, stated not hidden); 43 admitted stations within 60 km of one.
* intraslab-earthquake S station terms available for 2 stations; 43 carry 0.0 and are FLAGGED, never imputed.

**ZK is vertical-only in the archive** (1,604 of 1,642 delivered ZK station-days carry a single EHZ channel, `results/P_STEP2_SITES.md` finding 1) and contributes nothing; **UW.GNW fails the 2008 site probe** at 1.7 dB and is not in the sum.

## 3. Pair stage and grouping
| chain | blocks_tremor | det_tremor | per_hour_tremor | blocks_quiet | det_quiet | per_hour_quiet | tremor_quiet_ratio | nmad_in_sd | sta_per_block |
|---|---|---|---|---|---|---|---|---|---|
| real | 115 | 33694 | 292.99 | 38.00 | 12502.00 | 329.00 | 0.89 | 5.08 | 51.00 |
| site-shift null | 115 | 1867 | 16.23 | -- | -- | -- | -- | 5.05 | 50.00 |

nmad distribution of the tremor-hour detections: p50 9.5, p90 15.5, p99 24.2, max 40.9.  8 × MAD sits at 5.08 sd of the network function — the same measured bar as every other bed (2005: 5.15 sd, 2010: 4.90 sd).

| chain | components_ge3_members | grown | growth_budget |
|---|---|---|---|
| real | 160 | 160 | 1500 |
| site-shift null | 36 | 36 | 1500 |

## 4. The screen (every clause is a 2010→2011 failure mode)
Applied to **every** grown template before anything is counted.  The 2010 run's "50 candidates vs 17 null" died in 2011 because 46/50 sat on a face of the search box, the templates fired 1.89 per template-day on NON-tremor days against 0.66 on tremor days, and their locations moved 42–150 km when re-stacked on other data (`results/WA_2011_REFIRE_RESULTS.md`).  Clauses (a)–(e) are those three facts turned into tests; the wide search box 46.5–48.8 N / −125.2…−121.8 / 15–60 km is used **from the start** so that "interior" means something.

Network false-alarm rates measured directly: real chain **1.525** detections per template-day on the 19 literature-quiet days against **1.541** on the 20 audited tremor days; null chain 0.257 vs 0.311.  (The 2010 templates re-fired in 2011 at 1.89 quiet vs 0.66 tremor — the sign of this comparison is the whole test.)

| clause | real alone | real cumulative | null alone | null cumulative |
|---|---|---|---|---|
| (c) located on >= 5 stations | 156 | 156 | 36 | 36 |
| (c) shuffled-station control passes | 19 | 19 | 2 | 2 |
| (c) INTERIOR to the wide grid | 28 | 11 | 4 | 1 |
| (c) |dz| <= 8 km of Slab2 | 16 | 0 | 3 | 0 |
| (d) >= 3 distinct tremor days | 110 | 0 | 17 | 0 |
| (a) quiet-day FA <= 0.1/tpl-day and >= 3x lower than tremor | 46 | 0 | 6 | 0 |
| (b) tremor co-timing > same-day base rate | 149 | 0 | 36 | 0 |
| (e) stability: 2nd-half re-stack relocates within 10 km | 0 | 0 | 0 | 0 |

> **CANDIDATES after the whole screen: REAL 0 vs NULL 0.**

Clause (e) is the only one evaluated on the survivors of (a)–(d) rather than on all templates: a re-stack of every grown template is a second full pass over the waveforms. The intersection is unaffected by the order.

### 4.1 real chain — templates passing the whole screen (0)
_(empty)_

### 4.2 site-shift null chain — templates passing the whole screen (0)
_(empty)_

Full per-template screen tables: `results/n_puget2008/2008/screen_final.csv` (real), `results/n_puget2008/2008_null/screen_final.csv` (null); every clause is a column.

## 5. Confirmation on 2007-01 (tranche T5, the January 2007 ETS)
The rule was fixed before any 2007 number existed (`scripts/n_confirm.py` docstring): **≥ 3 distinct 2007 tremor days with detections, quiet-day false alarms ≤ 0.1 per template-day on the 2007 literature-quiet days, and an INTERIOR 2007 location within 10 km (3-D) of the 2008 one.**  The templates are read straight out of the 2008 run's `grow/tpl_it0.npz` — no re-stacking, no re-cutting, no re-thresholding.  The surviving site-shift-null templates go through byte-identical code on the identical 2007 data: **their confirmation count is the false-confirmation rate.**

Re-fired: 0 real and 0 null templates over 8 audited 2007 tremor days and 15 2007 literature-quiet days.

_(empty)_

## 6. Descriptive comparison with the Chestler-derived coordinates
`../codex-LFE/results/olympic_reference_discovery_v1/derived_coordinate_rows_not_a_catalog.csv` is, as its filename says, **not a catalogue**; it is used only to say where a confirmed family would sit relative to previously published Olympic LFE coordinates, never as truth.

* no confirmed family to compare.

## 7. Files
| what | where |
|---|---|
| hourly 2–8 Hz energy per station-day | results/n_puget2008/hour_energy28_{2008,2007}.parquet |
| hour audit (per hour, every clause) | results/n_puget2008/hour_audit_{2008,2007}.csv + .json |
| selected tremor hours / days | results/n_puget2008/tremor_hours_{2008,2007}.csv, tremor_dates_*.txt |
| pair stage (tremor / quiet) | results/n_puget2008/2008/auto/, auto_quiet/, pair_all.parquet, pair_pairstats.json |
| groups | results/n_puget2008/2008/gb_nodes.csv, gb_summary.csv |
| templates + tremor-day detections | results/n_puget2008/2008/grow/tpl_it0.npz, det_it0.parquet |
| quiet-day detections (clause a) | results/n_puget2008/2008/det_fa.parquet |
| locations (wide grid, edge-flagged) | results/n_puget2008/2008/loc_it0.csv |
| screen tables | results/n_puget2008/2008/screen.csv, screen_final.csv |
| stability re-stacks | results/n_puget2008/2008/reloc_h1.csv, reloc_h2.csv, reloc_all.csv |
| whole-chain site-shift null | results/n_puget2008/2008_null/ (offsets in results/n_puget2008/2008_null/null_offsets.json) |
| 2007 re-fire + confirmation | results/n_puget2008/2007/ |
