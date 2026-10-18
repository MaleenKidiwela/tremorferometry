# SVI backward extension 2003-01 → 2013-12, every day — scope

Written 2026-10-07, 2:07–3:00 pm PT. Scoping only: metadata queries and 10-minute FDSN samples, parsed in memory.
Nothing large was downloaded, nothing was deleted, and no production job was started.
Scripts: `svi_science/scripts/bw_scope_{inventory,sample,retry,truth,coverage,volume}.py`. Tables are in this folder.

"Scrambled control" = the same detector run on the same day, with each station site's data shifted to a different
time of that day.

## 0. In short

- **Availability.** The earlier finding (lfe_picker `P0O_RESULTS.md` §O1, "no reachable centre serves CN/C8
  before 2010") is **half right**:
  - Confirmed: EarthScope returns HTTP 204 for every CN and PO request before 2010-01. That is 307/307 CN
    requests, re-tried. From 2010 it delivers CN 4/4.
  - Wrong: **NRCan/CNDC FDSN is now reachable** at `https://www.earthquakescanada.nrcan.gc.ca/fdsnws` (the host
    without `www.` returns 301). It serves CN 2003–2009 and PO (POLARIS) 2003–2006.
  - The catch: any single NRCan request succeeds only **~1/3 of the time**. The same URL returns 404 "Not found",
    204 or 200 at random. Measured on dates where EarthScope has the data: NRCan first-pass success 0.33 (197
    requests). Five serial retries recovered 61 % of the CN misses, at the same rate before 2010 (59 %) as after
    (66 %). So the pre-2010 CN misses look like server flakiness, not gaps. Once a request succeeds, a 3-C day
    (16 MB) arrives in about 1 s.
  - C8 POLARIS-era sites (MGCB, TWBB, SHDB, SHVB, KHVB) are at EarthScope for 2003–2005. PB boreholes are at
    EarthScope from mid-2005, with 3–4 of 4 samples per year.
- **The pre-2010 template network is not "PB only".** It grows in steps:
  - CN4 (PGC, VGZ, NLLB, SNB) from 2003;
  - plus the Olympic PB stations from 2005-07;
  - plus B009–B011 from 2005-11;
  - plus B926 from 2007-11 (= the forward S14 set);
  - plus LRIV from 2009-08 (= S15).
- **Coverage, under the forward carrier rule (≥ 2 live carrier sites):**
  - CN4 era: 112/115 families have ≥ 2 sites and 100 have ≥ 3.
  - From 2005-11: 115/115 have ≥ 2 sites and 108–111 have ≥ 3, the same as the forward S13–S15 sets.
  - Coverage is not the bottleneck. **Sensitivity of the 4-station network** and **template transfer at VGZ/SNB**
    are, and both can be measured on data already on disk before any download.
- **Volume and time.**
  - Download: 0.48 TB for the 15 template stations (0.87 TB raw in total; 21,167 pre-2014 station-days are
    already on disk). With every non-template station within 80 km added: 2.0 TB.
  - 20 Hz cache: 0.99 TB. Matched filter (data + scrambled control): about 4.6 GPU-h for 4,018 days at the
    forward speed. All stages together: under 1 GPU-day.
  - Free disk: 58 TB of 100 TB.
- **Recommendation: go, in tiers, with gates.**
  - **2010–2013 is tier 1.** Same S15 network and template epoch as 2014–2019, and PNSN tremor is available.
  - **2005-11 → 2009 is tier 2.** Usable after the VGZ/SNB epoch-identity gate and a quiet-day rule that does not
    depend on PNSN.
  - **2003 → 2005-10 is tier 3.** CN4 only. Between-episode use only if a CN4 ablation on 2010–2013 truth days
    passes the control bars; otherwise ETS windows only.

## 1. Station availability 2003–2013

Method (`bw_scope_inventory.py`, `bw_scope_sample.py`, `bw_scope_retry.py`):
- Station metadata (level=channel) from both centres, for a box around the families. A station is kept if it is
  ≤ 60 km from the nearest family, or if it is a template station: 84 permanent stations and 218 temporary-array
  stations (XG, XU, YW, 4P).
- For every permanent station-year, four 10-minute requests (Feb/May/Aug/Nov 15, 10:00 UTC; HH?/BH?/EH?/SH?) at
  every centre that lists the station: 2,707 requests.
- Then every NRCan miss was re-tried up to 5 times, serially, with the exact Z channel (892 rows).
- Tables:
  - `availability_matrix.csv`: delivered / in-epoch samples per station-year; "-" = no metadata epoch.
  - `availability_station_year.csv`: centres, channels, rates, first-pass hits per centre.
  - `samples_raw.csv`, `samples_retry.csv`: one row per request.
  - `inventory_channels.csv`: every channel epoch, sensor and gain.

Template stations (delivered / in-epoch samples, any centre after retry):

| sta | net · cha · fs (2003–2013) | centre | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PGC | CN BH 40 Hz (CMG3T, GD1, 2000-06 → 2015-01, one epoch) | NRCan <2010, both ≥2010 | 4/4 | 3/4 | 3/4 | 2/4 | 4/4 | 2/4 | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 |
| VGZ | CN EH 100 (S-13) → BH 40 (2003-02-27) → HH 100 (2009-07-31) | as PGC | 3/4 | 3/4 | 3/4 | 3/4 | 3/4 | 4/4 | 3/4 | 4/4 | 3/4 | 4/4 | 4/4 |
| SNB | CN BH 40 → HH 100 (2011-03-17) | as PGC | 4/4 | 2/4 | 3/4 | 3/4 | 2/4 | 4/4 | 4/4 | 4/4 | 4/4 | 3/4 | 4/4 |
| NLLB | CN BH 40 | as PGC | **0/4** | 4/4 | 3/4 | 4/4 | 3/4 | 4/4 | **1/4** | 4/4 | 4/4 | 4/4 | 4/4 |
| B001/B003/B004/B005/B006/B007 | PB EH1/EH2/EHZ 100 (from 2005-06/07) | EarthScope | – | – | 1–2/2 | 3/4 | 2–4/4 | 3–4/4 | 3–4/4 | 3–4/4 | 4/4 | 4/4 | 4/4 |
| B009/B010/B011 | PB EH 100 (from 2005-11-02) | EarthScope | – | – | 1/1 | 4/4 | 3–4/4 | 4/4 | 4/4 | 4/4 | 4/4 | 3–4/4 | 4/4 |
| B926 | PB EH 100 (from 2007-11-13) | EarthScope | – | – | – | – | 0/1 | 3/4 | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 |
| LRIV | UW BH 40 (from 2009-08-27) | EarthScope | – | – | – | – | – | – | 0/2 | 4/4 | 4/4 | 4/4 | 4/4 |

For CN, a remaining miss is a lower bound on availability, because NRCan fails about 2/3 of single requests. NLLB
2003 (all four dates) and 2009 (three dates) failed ~8 times each and may be real gaps; a pilot month will tell.

Non-template stations within ~60 km that a template could be built for (`candidate_new_carriers.csv`). "Families
within r50": existing carriers sit at a median 47.5 km from their families.

| station | years live | median distance to families | families within r50 | note |
|---|---|---|---|---|
| CN LZB (BH 40) | 2003–2013 (ends 2014-06) | 26.6 km | 103 | **best single addition**; continuous; truth-year raw already on disk (2010/11/13) |
| CN PFB (HH 100) | 2003–2013 | 55.7 km | 42 | 2007-07 sensor change |
| CN YOUB, GOBB (EHZ only) | 2003–2013 | 60 km | – | GOBB is Z-only |
| C8 SHVB, MGCB, TWBB, SHDB, KHVB | 2003–2005/06 | 25–43 km | 73–106 | POLARIS (Bostock's array); EarthScope; MGCB has a flat-channel sample (0.70) |
| PO TWKB, TSJB, KELB/KLNB, SILB, SSIB, TWGB | 2003–2006 | 30–47 km | 59–99 | POLARIS; NRCan only; ends 2005–06 |
| C8 SOKB | 2008, 2013 partial | 24.2 km | 105 | patchy at both centres (matches the 2026-09 probe "0/8") |
| C8 PA01–PA12 | 2011–2013 | 38 km | ~80 | PGC mini-array, NRCan only, patchy |
| UW STW (EHZ only) | 2003–2013 | 42 km | 79 | single component, flat-sample fraction up to 0.79: excluded by the frozen 3-C rule |
| UW OPC (BH) / TA B04A (BH) | 2003–08 / 2005-10 → 2008-02 | 51–52 km | 46–49 | B04A sits at the LRIV site |

The temporary arrays are listed in `inventory_stations.csv`: XG (Array of Arrays, 2009–2011), XU (CAFE,
2006–2012), YW (FACES, 2005 / 2007–2010) and 4P (2007). They are short-lived and not useful for an every-day series.

## 2. Family coverage per year (forward carrier rule)

The carrier rule: a family counts a station if it has a usable 2012 template there (`s0/families.csv` carriers).
Sites are merged: PGC+B009–B011 is one site, LRIV+B005–B007 is another. The detector needs ≥ 2 live sites.
Table: `coverage_forward_reference.csv`. All stations are assumed live; delivery is the §1 matrix.

| network | years | stations | ≥ 1 | ≥ 2 | ≥ 3 | ≥ 2 sites | ≥ 3 sites | median carriers / sites |
|---|---|---|---|---|---|---|---|---|
| **forward S15** | 2014-01..2019-05 | 15 | 115 | 115 | 115 | 115 | 111 | 10 / 6 |
| **forward S14** | 2019-05..2021-03 | 14 | 115 | 115 | 115 | 115 | 111 | 9 / 6 |
| **forward S13** | 2021-04.. | 13 | 115 | 115 | 115 | 114 | 108 | 8 / 5 |
| B0: CN4 | 2003-01..2005-06 | 4 | 115 | 112 | 100 | 112 | 100 | 4 / 4 |
| B1: + Olympic PB | 2005-07..2005-10 | 10 | 115 | 114 | 107 | 114 | 107 | 6 / 5 |
| B2: + B009–B011 | 2005-11..2007-10 | 13 | 115 | 115 | 115 | 115 | 108 | 9 / 5 |
| B3: + B926 (= S14) | 2007-11..2009-08 | 14 | 115 | 115 | 115 | 115 | 111 | 9 / 6 |
| B4: + LRIV (= S15) | 2009-08..2013-12 | 15 | 115 | 115 | 115 | 115 | 111 | 10 / 6 |

Per-year coverage from the measured samples is in `coverage_by_year.csv` (median over the four sample dates;
CN under-counted because of the NRCan flakiness):
- 2006–2013: 112–115 families with ≥ 2 sites.
- 2003: 109; 2004: 97.5; 2005: 110.5.
- The "epoch_same" variant counts only stations on their 2012 instrument epoch. It drops 2003–2005 to 0–93,
  because VGZ (BH before 2009-07-31) and SNB (BH before 2011-03-17) then leave the network. §3 is about this.

What it means for detectability:
- The detector's threshold is set per network from quiet-day false alarms (forward E0–E3: T_core 8.75 → 9.25).
  Recall against Bostock was 0.12–0.21 for 13–15 stations, and 0.21 on 2010 with the 9-station S9 network.
- From 2005-11 the backward networks are the forward networks, minus at most B926/LRIV, so recall should be in
  the forward range.
- CN4 has no forward analogue. Its sensitivity is unknown until it is measured (gate G1).
- PGC carries > 40 % of most families' weight (Track A step 3), and PGC is present and on its 2012 instrument
  throughout 2003–2013. That is the main reason CN4 may still work.

## 3. Templates: can the 2012 templates be used on 2003–2009 data?

Instrument epochs, from FDSN metadata (`inventory_channels.csv`). The chain removes no response and correlates
1–8 Hz band-passed counts at 20 Hz, so gain and sample rate do not matter; sensor or digitiser phase in 1–8 Hz
could.

| station | 2012 build epoch | earlier epochs in 2003–2009 | families with a template | transfer status |
|---|---|---|---|---|
| PGC | BH CMG3T/GD1 (2000-06-28 → 2015-01-12) | **same** | 110 | no issue |
| NLLB | BH CMG3ESP-60/GD2 (2007-02-27 → 2017-02) | same sensor 2003-01 → 2007-01-27; CMG3T 2007-01-27 → 02-27 | 88 | 1 month differs |
| PB B001–B011, B926 | EH HS-1-LT/Q330 since install | **same** | 41–112 each | no issue (B926 is A-2's no-change control) |
| LRIV | BH (2009-08-27 →) | same | 34 | no issue |
| **SNB** | HH CMG3ESP/GD2_v21 (2011-03-17 →) | BH CMG3ESP/GD1 2008-12 → 2011-03; BH CMG3ESP/GD2_v12 2006-04 → 2008-12; BH CMG40T/GD1 → 2006-04 | **101** | GD1 epoch already used: the 2010 S9 truth run carried SNB-BH with 2012 templates and scored AUC 0.867 (dev table). Earlier two epochs untested |
| **VGZ** | HH CMG3ESP/GD2_v21 (2009-07-31 →) | BH CMG3T/GD2 2007-02 → 2009-07; BH CMG3ESP/GD1 2003-02-27 → 2007-02; EH S-13 short period → 2003-02-27 | **97** | untested; VGZ was the weakest transfer at the 2017 swap (cc ratio 0.72, the same as the no-change control 0.79) |

The test already exists. Track A-2 step 4 (`lfe_extension/results/TRACK_A2_RESULTS.md`):
- Use detections with the target station withheld.
- Stack the target's 3-C at the predicted S.
- Compare with the 2012 stack: cc ratio and registration ratio, against a no-change station across the same gap
  (B009–B011 / B926 / PB).
- It found 0.84–0.90 ratios at PGC, NLLB and SNB across the 2017 Trillium swap, dS = 0.00 s at every boundary, and
  no polarity flips.
- Backward it runs on the ETS windows: 2008-05 and 2009-05 for VGZ-BH-3T and SNB-GD2; 2004-07 and 2005-09 for
  VGZ-BH-3ESP and SNB-40T; 2007-01 is small. PB are the no-change control.

Where a transfer fails, the templates can be rebuilt from Bostock's catalogue on that epoch:
- Truth times plus the frozen 2012 per-family offsets, run through the `e5r_templates` code.
- Use the alternate 2-h build/test blocks, as in 2012, so recall is never scored on build blocks.
- Events per epoch: VGZ-3T has 2008 + 2009 = 35,156; VGZ/SNB pre-2006 have 2003–2005 = 141,738; SNB-GD2 has
  2007 + 2008 = 16,496.
- At most 2 stations and 97–101 families per epoch are affected. Track A-2 also has the per-epoch `--extra-tpl`
  mechanism ready.
- VGZ before 2003-02-27 (EH short-period, 57 days) is not worth a template: drop VGZ for those days.

Optional, a separate product: era templates on LZB (2003–2014, 103 families within 47 km) and POLARIS sites
(2003–2005).
- They raise CN4-era sensitivity. Bostock built his catalogue on exactly these stations.
- But they change the detector. Report them as an "era-augmented" variant beside the frozen detector, never
  spliced into the 23-year series.

## 4. Data volume, time, storage

`volume_summary.csv` and `volume_station_year.csv`:
- Station-days = metadata-epoch days × delivered fraction (lower) or all epoch days (upper).
- MB per day = median of our own past downloads per station.
- The download subtracts the 21,167 pre-2014 station-days already on disk (`ondisk_station_days_pre2014.csv`).
  Most are Olympic PB 2007–2013, VGZ/PFB/LRIV 2010–2013 and the truth windows. They are read-only inputs owned by
  other projects.

| plan | stations | station-days | raw total | to download | 20 Hz cache |
|---|---|---|---|---|---|
| A: frozen detector, 15 template stations | 15 | 42.9k – 47.4k | 0.80 – 0.87 TB | **0.42 – 0.48 TB** | **0.89 – 0.99 TB** |
| B: A + non-template stations within 80 km median distance | 60 | 108k – 130k | 1.95 – 2.40 TB | 1.55 – 1.99 TB | 2.24 – 2.71 TB |

- **Per station-year at native rate:**
  - CN BH 40 Hz: 9–15 MB/day, 3.3–5.5 GB/yr.
  - PB EH 100 Hz: 15–21 MB/day, ~7 GB/yr.
  - POLARIS HH 100 Hz: 16–54 MB/day, 6–20 GB/yr.
  - Cache: 20.7 MB per station-day (measured), 7.6 GB per station-year.
- **Download time:**
  - EarthScope: PB is 0.3 TB at the 63 MB/s seen in Track B, about 1.5 h.
  - NRCan: CN4 2003–2009 is ~10k station-days. At ~3 requests per success and ~1–2 s per request, with 4 lanes,
    that is about 2–4 h.
  - POLARIS for plan B: about 25k more NRCan/EarthScope station-days, ~0.5 day.
  - The download must retry until success (≥ 10 tries; see risk R1).
- **Cache (CPU):** ~70 s per station-month-set in the forward run, so about 2–3 h for 132 months.
- **GPU** (forward production: 4,641 days in 19,115 s on 3 shards = 4.1 s per network-day, data + scrambled
  control + features; `lfe_extension/results/refine/logs/prodmf_*.log`):
  - Production 2003–2013: 4,018 days, **~4.6 h**.
  - Station-set ablations on the 2010–2013 truth days already cached (B0/B1/B2, ~180 days each): ~0.6 h per set.
  - Epoch-identity stacks: < 1 h. Optional era-augmented run: ~5 h.
  - **Total: 8–12 GPU-h, under 1 day on the H100.**
- **Storage:** `df /home/jovyan` gives 100 TB, 43 TB used, **58 TB free**. Plan A peaks at about 1.5 TB (raw +
  cache) and plan B at about 5 TB. The cut-window policy is on hold: the user stopped the WA cut on 2026-10-07, so
  no cut stage is planned until that is re-decided.

## 5. Truth and calibration

Bostock's catalogue (`truth_by_year.csv`, `truth_episodes.csv`):
- 269,423 events; all 115 families map exactly.
- **Ten ETS windows only, 230 days**: 2003-02/03 (13 d, 31k events), 2004-07 (26 d, 55k), 2005-09 (23 d, 55k),
  2007-01/02 (11 d, 2.8k), 2008-05/06 (24 d, 13.6k), 2009-05/06 (29 d, 21.5k), 2010-08/09 (24 d, 17.7k),
  2011-08/09 (22 d, 18.0k), 2012-09 (27 d, 27.4k), 2013-09/10 (31 d, 14.1k).
- No events in 2006.

How the catalogue calibrates, per network:
1. **The forward models stay frozen.** WAVE/FULL per S15/S14/S13, sha256 686e4de45c423e11.
   - They apply directly to B3 (= S14, 2007-11..2009-08) and B4 (= S15, 2009-08..2013).
   - 2008 and 2009 truth are new, **fully out-of-sample** checks: never fitted, different years, and the
     VGZ/SNB epochs differ.
   - 2010 (S9) was reported, not fitted. 2011 and 2012 were fit years; 2013 was the sealed year, and its one look
     is spent.
2. **New networks B0/B1/B2** (CN4, CN4+Olympic, 13 stations) have no forward model. Fit them exactly like the
   S13/S14 reruns:
   - Re-run the 2011 and 2012 truth days with the station set ablated.
   - Fit out-of-fold by day on 2011 and the 2012 test blocks. No new look at 2013.
   - Then score them **once each** on their own era's truth: 2003 and 2004 for B0; 2005-09 for B1; 2007-01 for B2.
   - Those are the out-of-sample referees, about 145k events.
   - Bars per set, as forward: AUC; calibration slope 0.8–1.2; matched/ΣP 0.8–1.2; control P > 0.5 ≤ 1 %
     (forward S13–S15 "same": 1.8–2.6 %).
3. **Recall and thresholds per year:**
   - Recall vs Bostock on every ETS window, read against the 2010–2013 band (0.12–0.21).
   - The k = his/ΣP factor per family and set (forward median 4.4–4.7) is fitted on 2010–2012 and *checked* on
     2003–2009. This is the honest test that rates are comparable across eras.
4. **Scrambled control on every day**, as forward (one per day, `c1`), on every network.

Between episodes, Bostock's catalogue cannot validate anything: it has no events outside ETS windows. What
replaces it, in the forward era's own terms:
- **Scrambled-control rate and ΣP per day** (data − control bracket; forward inter-ETS mode
  "control_subtracted_bracket"). Report control/data per network era (forward E0–E3: 0.34–0.50).
- **Withheld-station S-energy test on quiet days.** This is the one inter-ETS check that passed forward
  (diff +1.8, Spearman 1.0), so repeat it per era.
- **Station-set invariance**: run 2010–2013 (and 2008–2009) with the B0/B1/B2 sets as well. The quiet-day rate
  ratio between networks on the same days tells whether CN4-era between-episode rates mean the same thing.
  Forward this failed (pooled ratio 1.09, monthly 1.24), so expect a bracket, not a pass.
- **Quiet-day stack test** (Track A step 3b): stack the quiet-day candidates against the family template.
- **Time-reversed control and burstiness/depth tests**, as in S7/A (`svi_science/results/s7`, `a_wa_between`).
- **The quiet-day flag needs a new source before 2010.** PNSN tremor starts in 2009-08 (the local catalogue
  starts in 2010-01). Options:
  - the array-energy hourly flag (`lfe_array_discovery w_energy_hours.py`, Merlin Q3 fallback);
  - days ≥ 30 d from any Bostock ETS window;
  - the Wech/WTD 2005–2014 tremor catalogue (not fetchable from here in 2026-09).
  - Decide this with Merlin before thresholds are derived; the forward T_core/T_ext come from quiet-day false
    alarms.

## 6. Recommended plan

Stages with gates. Each stage writes a DONE/FAILED marker, as in the shell queues. Nothing downloads before G1.

| # | stage | what | gate (go / no-go) |
|---|---|---|---|
| S0 | inventory (done) | this scope | – |
| S1 | **ablation on disk** (no download) | Re-run 2010–2013 truth days with B0 (CN4), B1, B2; fit B0–B2 models (2011 + 2012 test, OOF) | **G1:** CN4 AUC(same) ≥ 0.80 and control P > 0.5 ≤ 3 % → tier 3 goes ahead for between-episode use; otherwise 2003–2005-10 is ETS-window only. B2 recall ≥ 0.8 × S14 → tier 2 goes ahead |
| S2 | pilot download | 2008-05 (ETS) + 2008-08 (quiet) for all 15 stations; NRCan retry loop ≥ 10 tries | **G0:** ≥ 95 % of metadata station-days delivered for CN; the NRCan error mode is understood. If CN < 80 %, stop and ask about NRCan bulk / AutoDRM |
| S3 | epoch identity (A-2 protocol) | VGZ-BH-3T, SNB-GD2 on 2008-05/2009-05; VGZ-BH-3ESP, SNB-40T on 2004-07/2005-09; PB = no-change control | **G2:** cc ratio and registration ratio ≥ the no-change control → keep the 2012 templates; else rebuild that (station, epoch) from Bostock events (build blocks only) and re-run G2 on the test blocks |
| S4 | full download + cache | Plan A, 2003–2013, ≤ 4 NRCan lanes, > 3 % station-day error guard; R2-style station-health audit (flat / dead / swap) per station-month | **G3:** cache loss ≤ 3 % per station-month; health flags reviewed |
| S5 | quiet-day rule + thresholds | Pre-2010 quiet flag decided with Merlin; T_core/T_ext per backward network (B0–B4) by the forward false-alarm rule | **G4:** ≥ 200 quiet family-days per network; B4's thresholds reproduce E0's 8.75/8.25 within one grid step (anchor) |
| S6 | matched filter | Frozen 2012 template set (+ the G2 rebuilds), data + scrambled control every day, features as SVI_REF | – |
| S7 | confidence | B3/B4 = frozen S14/S15 models; B0–B2 = S1 models. One look each at 2003/04/05/07/08/09 truth | **G5:** per-era bars (§5); a failing era is reported as a ΣP bracket only |
| S8 | controls | Scrambled control, withheld S-energy, station-set invariance (B0 vs B4 on 2010–2013), quiet-day stack, time-reversed | **G6:** decides per era whether between-episode rates are "rates" or a "control-subtracted bracket" |
| S9 | catalogue | catalogue v1.4-backward: same columns as SVI_REF, plus `era`, `net_set`, `tpl_epoch`; rates from ΣP with step checks at every era/epoch boundary | – |

**Usable years for between-episode science (expected):**
- **2010-01 → 2013-12 (tier 1).** Same S15 network and same template epoch as 2014–2019 (except SNB-BH until
  2011-03, already exercised in 2010). PNSN tremor quiet flags exist. This extends the forward series by 4 years
  with the least risk; it is the first thing to run.
- **2005-11 → 2009-12 (tier 2).** Networks of 13–14 stations (B2/B3), PGC on the same instrument. Usable if G2
  passes for VGZ/SNB, or after the rebuilds, and G4 finds a pre-2010 quiet rule.
- **2003-01 → 2005-10 (tier 3).** CN4, or CN4 plus Olympic PB from 2005-07. Between-episode use only if G1 passes.
  ETS-window recall is worth having either way: 2003/2004/2005 are the three largest truth episodes.
  - The era-augmented variant (LZB + POLARIS templates) is the fallback for sensitivity, as a separate product.

**Main risks**
- **R1. NRCan flakiness.** ~1/3 success per request; 204 and 404 are both returned for data that exist.
  - A plain obspy fetch with tries = 3 gets ~70 %. That is what made the 2005 bed partial (PGC 35/44 days) and led
    to the old "CN not available" verdict.
  - Mitigation: explicit channel lists (the BH? wildcard failed more often than BHZ in tests), ≥ 10 retries with
    pauses, and a per-station-day completeness ledger.
  - NLLB 2003 and 2009 may be real gaps.
- **R2. Template transfer at VGZ/SNB before 2009/2011:** 97–101 families. Mitigated by G2 and the rebuilds.
- **R3. CN4 sensitivity 2003–2005.** No forward analogue; G1 measures it on disk before any download.
- **R4. No pre-2010 tremor catalogue for quiet days and thresholds.** This is a method decision for Merlin (S5).
- **R5. Station sets change 5 times in 2003–2013.** Every boundary is a potential rate step; step checks are
  mandatory, and a step next to an ETS is not separable from activity (forward R7 lesson).
- **R6. PB borehole start-up 2005–2006.** Partial samples (B004 1/2 in 2005; 2005 bed B004–B007 ~50 % of days)
  and possible orientation or noise issues. Mitigated by the R2-style health audit.
- **R7. Truth use.** 2013 is spent. 2011/2012 are fit years. 2003–2009 are untouched and must each be looked at
  once per network model. Freeze the day roles and model hashes before S7.

## 7. Files

- Scripts: `svi_science/scripts/bw_scope_inventory.py`, `bw_scope_sample.py`, `bw_scope_retry.py`,
  `bw_scope_truth.py`, `bw_scope_coverage.py`, `bw_scope_volume.py`.
- Tables (this folder):
  - `inventory_stations.csv`, `inventory_channels.csv`;
  - `samples_raw.csv`, `samples_retry.csv`;
  - `availability_station_year.csv`, `availability_matrix.csv`;
  - `coverage_by_date.csv`, `coverage_by_year.csv`, `coverage_forward_reference.csv`;
  - `candidate_new_carriers.csv`;
  - `truth_by_year.csv`, `truth_episodes.csv`, `truth_family_year.csv`;
  - `ondisk_station_days_pre2014.csv`, `volume_station_year.csv`, `volume_summary.csv`.
- Logs: `sample.log`, `retry.log`, `coverage.log`. The `*.partial.csv` files are intermediate checkpoints, kept.
