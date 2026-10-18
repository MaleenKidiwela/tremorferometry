# Washington Puget / S-Olympic 2008, QUIET-SUBSET bed — station quality as the lever

Agent P, 2026-09-21.  **The recipe is unchanged** (`results/REPRO_BOSTOCK.md` §0): 6-s windows, 1–8 Hz, per-component cc averaged over all admitted stations of the bed, median + 8 × MAD of the block's own distribution, |lag| ≥ 6 s, 4-s declustering, connected-component grouping at the 4-s node merge (≥3 members), Bungum–Husebye realignment, 12-s station-tapered templates, network matched filter at 8 × MAD, no waveform-similarity merge.  Agent K's `k_*.py` are unedited and agent N's `n_*.py` drive them; **the only thing that changed is the STATION SET** (`n_common.admitted()` honours an optional `N_STATION_SUBSET` whitelist — unset it and every agent-N number is reproduced byte-for-byte).  The audited tremor hours, the literature-quiet days, the health veto, the screen and the 2007 confirmation rule are agent N's, reused unmodified.

**Why.**  On the 51-station CAFE+PB bed the pair stage had a tremor/quiet ratio of **0.89** (`results/WA_2008_RESULTS.md` §3) — it fired as often on a literature-quiet day as inside the audited May-2008 ETS — and every grown template was a noise repeater.  Royer & Bostock (2014) built *composite arrays* from the stations that worked.  §1 below measures, station by station, who is manufacturing those pairs; §2–§5 re-run the unchanged recipe on the quiet subset the census defines.


## 1. The per-station noise-repeater census

`scripts/p_census.py` recomputes, at **exactly** the (i, j) pairs agent N's pair stage stored, each station's own 3-component mean cc, and compares it with that station's own median + 8 × MAD over 4,000 random pairs of the same block (the recipe's own threshold rule applied to one station instead of to the network).  Validation: the mean of the per-station values reproduces the stored network cc to 7 × 10⁻⁸.  Two numbers per station:

* **ratio** = contributions-above-its-own-8 × MAD per audited TREMOR hour ÷ the same per literature-QUIET hour (counts Laplace-smoothed by 0.5 on both sides; `ratio4`, the 10×-denser 4 × MAD version, is the stability check and is in the CSV);

* **top3_share_quiet** = the fraction of the QUIET-hour network detections in which this station is one of the top 3 contributors (chance = 3/N ≈ 0.059).

A **noise-repeater station** (the brief's definition) has ratio < 2 and top3_share_quiet above twice chance: **8 of 51 stations** qualify.


### 1.1 By network

| net | n | median_ratio | median_top3_share_quiet | median_per_hour_quiet | median_noise_db | n_noise_repeater |
|---|---|---|---|---|---|---|
| PB | 10 | 0.352 | 0.054 | 21.100 | 5.669 | 4 |
| XU | 40 | 2.678 | 0.007 | 0.609 | 19.682 | 4 |
| YW | 1 | 133.843 | 0.002 | 0.000 | 18.437 | 0 |


### 1.2 The eight worst (lowest ratio) and the eight best

**Worst (the pair stage's own noise engine):**

| netsta | ratio | ratio4 | per_hour_tremor | per_hour_quiet | top3_share_quiet | noise_db | noise_repeater |
|---|---|---|---|---|---|---|---|
| XU.BD01 | 0.013 | 0.126 | 2.877 | 222.692 | 0.625 | 21.205 | True |
| PB.B003 | 0.015 | 2.644 | 0.009 | 0.846 | 0.007 | 0.000 | False |
| XU.E070 | 0.016 | 0.146 | 0.000 | 0.250 | 0.003 | 19.608 | False |
| PB.B005 | 0.025 | 0.344 | 2.191 | 89.588 | 0.479 | 8.207 | True |
| XU.N040 | 0.025 | 1.353 | 0.000 | 0.156 | 0.006 | 15.161 | False |
| XU.N020 | 0.028 | 0.187 | 0.000 | 0.138 | 0.004 | 16.521 | False |
| PB.B006 | 0.035 | 0.551 | 1.878 | 54.241 | 0.147 | 5.898 | True |
| XU.BD11 | 0.037 | 0.201 | 3.157 | 85.844 | 0.519 | 21.659 | True |


**Best (tremor-locked contributors):**

| netsta | ratio | ratio4 | per_hour_tremor | per_hour_quiet | top3_share_quiet | noise_db | noise_repeater |
|---|---|---|---|---|---|---|---|
| XU.E030 | 1165.310 | 14.791 | 18.204 | 0.000 | 0.003 | 22.017 | False |
| XU.BS02 | 268.045 | 11.268 | 29.313 | 0.094 | 0.010 | 18.308 | False |
| XU.N050 | 174.470 | 8.778 | 2.722 | 0.000 | 0.002 | 15.709 | False |
| PB.B201 | 168.200 | 12.146 | 2.896 | 0.000 | 0.003 | 5.515 | False |
| YW.FACA | 133.843 | 7.760 | 2.087 | 0.000 | 0.002 | 18.437 | False |
| XU.W080 | 103.235 | 10.753 | 1.609 | 0.000 | 0.003 | 14.723 | False |
| PB.B001 | 76.522 | 10.346 | 1.191 | 0.000 | 0.002 | 3.500 | False |
| XU.BS06 | 50.312 | 6.165 | 16.504 | 0.312 | 0.024 | 12.313 | False |


`noise_db` is the P4-0 statistic — the station's median 2–8 Hz hourly power on the literature-quiet days, in dB above the QUIETEST PB borehole (0 dB = 1.81 in the cache's units).  Full table: `results/p_2008/station_noise_census.csv`; per-station-block rows in `census_tremor.parquet` / `census_quiet.parquet`.


## 2. The quiet subset

Rule, fixed in `scripts/p_subset.py` before any subset result existed: **the 10 PB boreholes in the bed, plus the CAFE (XU)/YW stations with ratio ≥ 3.0 and top3_share_quiet ≤ 0.118** (= 2 × chance), ranked by (ratio desc, top3 asc), floor 12 and cap 20 stations.

* CAFE/YW stations clearing both bars: **20**; taken: 10.

* **SUBSET (20 stations): PB.B001, PB.B003, PB.B005, PB.B006, PB.B007, PB.B013, PB.B018, PB.B201, PB.B941, PB.B943, XU.BS02, XU.BS04, XU.BS06, XU.E030, XU.N050, XU.N060, XU.PL04, XU.PL05, XU.PL11, XU.W060**

* boreholes-only control bed (10): PB.B001, PB.B003, PB.B005, PB.B006, PB.B007, PB.B013, PB.B018, PB.B201, PB.B941, PB.B943


**Coverage of the Puget/S-Olympic active cells** (Slab2 28–45 km with ≥30 PNSN **2010–2011** epicentres within 8 km — there is no 2008 catalogue, so these are 2010–11 cells used for a 2008 bed, stated not hidden), reported as `P_STEP2_SITES.md` §4 does:

| bed | n_cells | pct_cells_ge5_60km_hyp | pct_cells_ge5_60km_epi | median_sites_60km_hyp | median_sites_60km_epi |
|---|---|---|---|---|---|
| 51-station (agent N) | 264 | 90.5 | 96.6 | 11.0 | 19.0 |
| quiet subset (20) | 264 | 53.0 | 70.5 | 5.0 | 7.0 |
| boreholes only (10) | 264 | 7.2 | 15.5 | 1.0 | 2.0 |


⚠ **Read the screen with this number.**  Only 53.0 % of the active cells have ≥5 subset stations within 60 km, against 90.5 % on the 51-station bed.  Clause (c) of the screen requires a location on ≥5 stations, so a template from a real family in a poorly covered cell can fail this run for a GEOMETRIC reason.  The pair ratio in §3 is the measurement this bed was built to make; the screen counts are a lower bound, not a census of what is there.


Bed as actually run: **20 stations** {"PB": 10, "XU": 10}, 115 audited tremor hours on 20 days, 19 literature-quiet days.  2007 confirmation bed: 16 stations {"PB": 7, "XU": 9} (absent from 2007: PB.B201, PB.B941, PB.B943, XU.BS06).


## 3. The pair stage — the number this whole run was launched to move

| bed | stations | blocks_tremor | det_tremor | per_hour_tremor | blocks_quiet | det_quiet | per_hour_quiet | tremor_quiet_ratio | nmad_in_sd |
|---|---|---|---|---|---|---|---|---|---|
| 51-station bed (agent N) | 51.00 | 115 | 33694 | 292.99 | 38.00 | 12502.00 | 329.00 | 0.89 | 5.08 |
| QUIET SUBSET | 20.00 | 115 | 88702 | 771.32 | 38.00 | 33888.00 | 891.79 | 0.86 | 5.24 |
| quiet subset — site-shift null | 19.00 | 115 | 4593 | 39.94 | -- | -- | -- | -- | 5.24 |
| boreholes only | 10.00 | 115 | 324017 | 2817.54 | 35.00 | 220432.00 | 6298.06 | 0.45 | 5.31 |


Reference points: Vancouver Island POLARIS 2005 (the bed on which this recipe reproduces 58/96 of Bostock's families) **11×**; N-Olympic XG/PB 2010 **4.2×**; Puget 51-station 2008 **0.89×**.


nmad of the subset's tremor-hour detections: p50 8.5, p90 11.6, p99 18.0, max 30.3; 8 × MAD sits at 5.24 sd of the network function.


⚠ **Caveat that must be read with the ratio:** the quiet days are calendar-remote from the ETS window and the CAFE stations are less completely cached on them, so the quiet blocks run on fewer stations than the tremor blocks (`stations` column above is the tremor-block median; the quiet-block median is in `pair_pairstats.json` as `sta_per_block_median_quiet`).  A bed with fewer stations has a noisier network function and a higher absolute 8 × MAD bar, which **suppresses** quiet-hour detections; part of any ratio above 1 can be this and not tremor.  The template-level clauses of §4 — where the SAME template is fired on both day sets through identical code — are not exposed to it.


| chain | components_ge3 | grown |
|---|---|---|
| real | 175 | 175 |
| site-shift null | 133 | 133 |


## 4. The five-clause screen (agent N's, unchanged), real vs whole-chain null

Every clause is a 2010→2011 failure mode turned into a test (`scripts/n_screen.py` docstring): (a) per-template quiet-day false alarms ≤ 0.1/template-day AND ≥3× below the tremor-day rate; (b) tremor co-timing above the same-day base rate; (c) a ≥5-station location that passes the shuffled-station control, is INTERIOR to the wide grid (46.5–48.8 N / −125.2…−121.8 / 15–60 km, used from the start) and sits within 8 km of Slab2; (d) ≥3 distinct tremor days; (e) a second-half re-stack that relocates within 10 km.  The site-shift null is the same code on a network whose inter-station coherence has been destroyed by within-partition day shifts.


Network false-alarm rates measured directly: real chain **11.228** detections per template-day on the literature-quiet days against **12.526** on the audited tremor days; null chain 3.549 vs 2.318.  (The 51-station bed read 1.525 vs 1.541; the 2010 templates re-fired in 2011 at 1.89 quiet vs 0.66 tremor — the SIGN of this comparison is the whole test.)


| clause | real_alone | real_cumulative | null_alone | null_cumulative |
|---|---|---|---|---|
| (c) located on >= 5 stations | 109 | 109 | 61 | 61 |
| (c) shuffled-station control passes | 27 | 27 | 1 | 1 |
| (c) INTERIOR to the wide grid | 14 | 4 | 4 | 0 |
| (c) |dz| <= 8 km of Slab2 | 24 | 0 | 13 | 0 |
| (d) >= 3 distinct tremor days | 142 | 0 | 100 | 0 |
| (a) quiet-day FA <= 0.1/tpl-day and >= 3x lower than tremor | 31 | 0 | 15 | 0 |
| (b) tremor co-timing > same-day base rate | 127 | 0 | 128 | 0 |
| (e) stability: 2nd-half re-stack within 10 km | 0 | 0 | 0 | 0 |


> **CANDIDATES after the whole screen: REAL 0 vs NULL 0.**


### 4.1 real chain — templates passing the whole screen

_(empty)_


### 4.2 site-shift null — templates passing the whole screen

_(empty)_


Full per-template tables: `results/p_2008/2008/screen_final.csv` (real), `results/p_2008/2008_null/screen_final.csv` (null); every clause is a column.


## 5. Confirmation on 2007-01

The rule was fixed before any 2007 number existed (`scripts/n_confirm.py` docstring): **≥3 distinct 2007 tremor days with detections, quiet-day false alarms ≤ 0.1 per template-day on the 2007 literature-quiet days, and an INTERIOR 2007 location within 10 km (3-D) of the 2008 one.**  Templates are read straight out of the 2008 run's `grow/tpl_it0.npz` — no re-stacking, no re-cutting, no re-thresholding.  The surviving null templates go through byte-identical code on the identical 2007 data: their confirmation count IS the false-confirmation rate.


Re-fired: 0 real and 0 null templates over 8 audited 2007 tremor days and 15 2007 literature-quiet days.


_(nothing survived the 2008 screen, so there was nothing to re-fire)_


## 6. What the census says about Vancouver Island vs Puget

The network autocorrelation asks one question of every 6-second window: *does this waveform repeat somewhere else in the hour, coherently, at many stations at once?*  LFEs answer yes because the same patch of the plate interface ruptures the same way.  So does a pump, a mill, a road, a badly-clamped sensor and a wind-loaded tree root — and those repeat all day, every day, tremor or no tremor.  The census measures exactly that: the CAFE surface stations sit a median 14 dB above the quietest borehole in the 2–8 Hz band on quiet days, and their median tremor/quiet contribution ratio is 2.87 against 0.35 for the boreholes — i.e. on the surface sites the repeating-waveform population is very nearly the same whether or not there is an ETS underway.  Because the recipe SUMS the per-component cross-correlation over every station, adding more such stations does not average the noise repeaters away: it adds their coherent, always-on pairs to the sum, so the 8 × MAD bar rises to accommodate them and the LFE pairs — which are only present during tremor and only at the stations that see that patch — sit below it.  Bostock's southern Vancouver Island stations are POLARIS forest sites on rock with almost no cultural neighbours; there the same sum is dominated by whatever genuinely repeats, which is the LFEs, and the tremor/quiet ratio is 11×.  This is why 'more stations' made Washington worse (0.89× on 51 CAFE+PB sites) and why the published Washington work built *composite arrays* out of chosen stations rather than throwing the whole deployment into one sum.


What this run then tests is the direct consequence: if the obstruction is the noise repeaters, deleting them should restore the pair stage. §3 reports whether it did, and §4–§5 report whether a restored pair stage is enough to produce a template that survives a quiet-day false-alarm test, a location, a stability re-stack and an out-of-sample re-fire — which is the bar a family has to clear, not the pair ratio.


## 7. Files

| what | where |
|---|---|
| per-station-block census rows | results/p_2008/census_tremor.parquet, census_quiet.parquet |
| station census + by-network table | results/p_2008/station_noise_census.csv, census_by_network.csv |
| subset definition, coverage, rule | results/p_2008/subset_meta.json, subset_stations.txt |
| pair stage (tremor / quiet) | results/p_2008/2008/auto/, auto_quiet/, pair_all.parquet, pair_pairstats.json |
| groups | results/p_2008/2008/gb_nodes.csv, gb_summary.csv |
| templates + tremor-day detections | results/p_2008/2008/grow/tpl_it0.npz, det_it0.parquet |
| quiet-day detections (clause a) | results/p_2008/2008/det_fa.parquet |
| locations (wide grid, edge-flagged) | results/p_2008/2008/loc_it0.csv |
| screen tables | results/p_2008/2008/screen.csv, screen_final.csv |
| whole-chain site-shift null | results/p_2008/2008_null/ (offsets null_offsets.json) |
| boreholes-only pair stage | results/p_2008/bh/ |
| 2007 re-fire + confirmation | results/p_2008/2007/ |
