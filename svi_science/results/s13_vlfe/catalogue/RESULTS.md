# catalogue RESULTS — stage 1 (2011–12 tests on disk + far-reference metadata)

Written 2026-10-15 ~1:20 pm PT. Rules: `PLAN_CAT.md` (written ~1:10 pm PT before any number; one amendment ~1:15 pm PT
before any E3 number). Log: `NOTES.md`. "Scrambled control" = catalogue from data scrambled so no real event can line up
(A1 station-day shift, A2 polarity scramble); "noise template" = a random real 90-s network window used as a template.
Periods: ETS = H&G windows 2011-07-23 – 09-06 and 2012-08-30 – 10-11; between ETS = 2011-09-08 – 2012-08-28 (0.97 yr).

## Plain answer (so far)
- **e′ works as intended.** It brings back the 6 ETS near-misses (8 → 14 real-looking ETS events) and nothing between
  ETS (0/130); scrambled controls go from 0 to 2/650 (rule ≤ 2). Accepted.
- **E2: a tier separates real events from controls on the existing V4 run.** At 5.5 MAD with two-template agreement and
  e′: 18 ETS events, 1 between-ETS event, 0 in any of the 16 scrambled catalogues and 0 in 10 noise templates. Accepted.
- **Far references exist in every year 2003–2026.** At least 5 LHZ stations with ≥ 180 days of data per year in
  2003–2005 (US.NEW, US.MSO, US.HAWA, US.WVOR, IU.COR; BK and TA added from 2004–2005), 40–110 in 2006–2008 (TA).
  All lie east to south (azimuth 77–186°). A fixed 5-station set for 2003–2026 needs ~111 station-years, ~5.6 GB
  (5.2–8.9 GB). Not downloaded.
- **E3/E4: the empirical templates work on the full network but not on the four pre-2009 Canadian stations.**
  - At the 0.5-per-template-day noise threshold, H&G 9 re-finds 7/8 vetted events with templates that do not contain
    them (rule ≥ 7: met), all 13 stations 6/8 (just below the rule), 4 stations {PGC, VGZ, NLLB, SNB} 2/8 (rule ≥ 5:
    fail) with scrambled controls 1.16 per year (rule ≤ 1: fail).
  - So by the pre-set rule, 2003–2009 can be catalogued only at Mw ≥ 4.0 (+ Poisson bound); the uniform catalogue
    starts 2010 (pilot 2008 not run).
  - Full-network bank (all 13): 49 ETS events vs 1 between ETS vs 0.13 per year in scrambled controls (battery with e′).
- The single between-ETS E2 event (2012-02-13 05:44 PST) does not look like the ETS events in raw waveforms
  (`figures/e2_wave_2012-02-13_054440.png`): no pulse above the background, and the far stations carry a coherent
  long-period wavetrain of similar amplitude. Most likely noise that passed the battery.

## 1. e′ validation (Merlin step 3; `tables/eprime_validation.csv`, `eprime_changed_candidates.csv`)
Raw figure: `figures/eprime_ratio_vs_z.png` (ratio_far vs scan z for every row, candidates red, scrambled purple,
noise grey; black ring = passes the battery with e′).

| set | rows | battery with (e) | battery with e′ | in a head2 window | head2 of M ≥ 5.5 | head2 + ≥ 2 far refs excited |
|---|---|---|---|---|---|---|
| ETS candidates | 115 | 8 | **14** | 65 | 22 | 30 |
| between-ETS candidates | 130 | 0 | **0** | 94 | 43 | 62 |
| between-ETS scrambled control | 650 | 0 | **2** | 305 | 53 | 67 |
| between-ETS selection-matched noise | 1,300 | 1 | 1 | 611 | 109 | 205 |
| ETS scrambled control (column) | 575 | 2 | 5 | 247 | 51 | 67 |
| ETS noise (column) | 1,150 | 2 | 2 | 541 | 79 | 140 |

- The 6 events e′ adds are exactly the 6 near-misses: T2 detections on 2011-08-21 09:02:56, 2011-08-26 02:47:26,
  02:54:26 and 03:01:09, 2012-09-18 06:06:38, 2012-09-20 09:42:52 PDT. Each sat in a head2 window of a M 4.6–5.0
  teleseism 8,000–11,000 km away with no far reference excited; ratio_far 0.17–0.36; all with LFEs.
- Pass rule (≤ 14 ETS, ≤ 2/650 scrambled, 0 between ETS): met, at the edge for the scrambled controls (2 of 650).

## 2. E2 — the battery on every V4 detection (Merlin step 2; `tables/e2_separation.csv`)
Raw figures: `figures/e2_mf_vs_scanz.png` (MF level vs scan z for every T1–T3 detection, real / 16 scrambled
catalogues / 10 noise templates, ETS and between ETS, coloured by ratio_far, red ring = battery with e′),
`figures/e2_separation.png` (counts vs MF level). Counts are distinct events (detections within ±30 s merged).
Inputs: 115,209 detections ≥ 5.5 MAD of T1–T3, 2011-06-01 – 2013-01-01 (real 4,710 rows in ETS + between, scrambled
19,286, noise 62,547).

Detections before the battery (distinct events):

| MF level (MAD) | ETS real | ETS scrambled (mean) | between real | between scrambled (mean) | between noise template (mean) |
|---|---|---|---|---|---|
| 5.5 | 1,376 | 271 | 3,280 | 1,058 | 2,599 |
| 6 | 489 | 60 | 962 | 237 | 1,009 |
| 7 | 78 | 1.6 | 85 | 7.7 | 149 |
| 8 | 14 | 0 | 8 | 0.06 | 20 |

After the battery (selected tiers; all 40 tiers in the table):

| MF level | earthquake test | S5 | two templates | ETS real | between real | scrambled mean / max (between) | scrambled per yr | noise per template-yr | rule |
|---|---|---|---|---|---|---|---|---|---|
| 5.5 | e′ | no | no | 33 | 12 | 1.44 / 3 | 1.48 | 5.0 | fail |
| 5.5 | e′ | no | **yes** | **18** | **1** | **0 / 0** | **0** | **0** | **pass** |
| 5.5 | e′ | yes | yes | 15 | 1 | 0 / 0 | 0 | 0 | pass |
| 5.5 | e | no | yes | 12 | 0 | 0 / 0 | 0 | 0 | pass |
| 6 | e′ | yes | no | 21 | 2 | 0.13 / 1 | 0.13 | 0.62 | pass |
| 6 | e′ | no | yes | 10 | 1 | 0 / 0 | 0 | 0 | pass |
| 7 | e′ | no | no | 12 | 0 | 0.06 / 1 | 0.06 | 0.21 | pass |
| 8 | e′ | no | no | 6 | 0 | 0 / 0 | 0 | 0.1 | fail (ETS < 8) |

- **Chosen tier (rule: lowest MF level, then e′): 5.5 MAD, e′, two-template agreement.** ETS 18 vs between ETS 1 vs
  0 in all 26 scrambled/noise catalogues.
- The battery passes about 0.1–0.3 % of scrambled and noise detections and 2.4 % of real ETS detections
  (`tables/e2_pass_fractions.csv`).
- Single-template tiers at 5.5 MAD pass 12 real between-ETS events, against 1.4 per scrambled catalogue and 5 per noise
  template. The noise templates fire more than the scrambled controls (as in hg_repro): they keep the real
  long-period wavefield. Two-template agreement removes all of them.
- **The one between-ETS event of the chosen tier:** 2012-02-13 05:44:40 PST, T3 (agreement 6.05 MAD), scan z 8.3 at
  48.52 N, 123.58 W (SVI family 79, 2.4 km), ratio_far 0.40 (5 refs), S5 pass, **2 SVI LFEs within ±10 min**. It sits in
  a head2 window of a M 4.8 teleseism (so (e) rejects it; e′ keeps it, no far reference excited). A candidate, to be
  looked at in waveforms.
- The ETS events of the chosen tier include the 8 R3 events found from Table S1 plus events Table S1 does not list or
  our R3 matching did not reach (`tables/e2_real_battery_events.csv`).

## 3. Far references 2003–2026 (step 5, metadata only; `tables/farref_*`)
Sources: EarthScope FDSN station service (LHZ channel epochs, IU/US/UW/TA/BK/CC/CN/II, 250–800 km from 48.37 N,
123.61 W); MUSTANG percent_availability (daily, QC metadata; the fdsnws availability service is retired);
NRCan station service for CN (no LHZ channels; BHZ/HHZ broadband listed as a separate class). No waveforms requested.

| year | LHZ stations ≥ 180 d (MUSTANG) | in metadata ≥ 180 d | stations (≥ 180 d) |
|---|---|---|---|
| 2003 | 5 | 8 | NEW 480 km 89°, MSO 744 km 100°, HAWA 378 km 124°, WVOR 765 km 148°, COR 421 km 177° |
| 2004 | 8 | 8 | + BK.MOD 764 km, BK.YBH 741 km, BK.HUMO 643 km |
| 2005 | 11 | 12 | + US.BMO 620 km, TA.I05A, TA.M02C |
| 2006 | 40 | 41 | TA (Transportable Array, Oregon/Washington) |
| 2007 | 110 | 110 | TA |
| 2008 | 57 | 59 | TA |
| 2009 | 37 | 39 | UW + TA |
| 2010 | 32 | 35 | |
| 2011–2016 | 24–40 | 26–44 | |
| 2017–2026 | 5–10 | 14–25 | MUSTANG has no UW LHZ values after 2018: use the metadata count |

- **IU.COR pre-2009 is real data** (Merlin F4 "unverified"): 353–365 days per year 2003–2008 in MUSTANG.
- All references sit at azimuth 77–186° (east to south). There is no LHZ reference north or west (the ocean; CN has no
  LHZ). The CN broadband class adds 7–22 stations per year at 26–356°, mainly north-west (`farref_cn_broadband_per_year.csv`;
  BHZ, ~1 GB per station-year at 40 sps, candidate windows only).
- **Download estimate** (`farref_download_proposal.csv`): fixed set IU.COR, US.NEW, US.HAWA, US.MSO, US.WVOR, all on air
  2003–2026 (HAWA to 2024): ~111 station-years still needed (2011–12 COR/NEW/HAWA on disk) × 47–81 MB (measured on our
  2011 LHZ day files, median ~50 MB) ≈ **5.6 GB (5.2–8.9 GB)**. Awaiting the user's OK.

## 4. E3/E4 — empirical templates, four station sets (Merlin step 4; `tables/e3_summary.csv`, `e3_refind.csv`)
Raw figures: `figures/e3_mf_vs_scanz.png` (MF z vs scan z for every detection at the 1-per-day threshold, real / scrambled
/ noise templates, per set, ETS and between ETS; rings = battery with e′), `figures/e3_loo_refind.png` (best
leave-one-out MF z per vetted event and set).
Bank (`tables/e3_bank.csv`, `e3_events.csv`): 14 events in 6 complete-linkage clusters (C0: 5 events, C1: 4, C2: 2,
C3–C5: 1 each); 6 leave-one-out stacks; H&G T2/T4 re-cut on all stations; 25 noise templates. Recipe V4 (velocity,
0.02–0.05 Hz, 90-s window from origin − 9 s). Stations: 13 distinct (Merlin's "16" counted VGZ, NLLB, LRIV twice).
Frame 2011-06-01 – 2013-01-01. Data gaps: YOUB has 135 of 914 frame days, PFB 334, SNB 510, PGC 555.

Thresholds from the noise templates (MAD units; `tables/e3_thresholds.csv`):

| set | 2011, 0.25 / 0.5 / 1 per template-day | 2012 |
|---|---|---|
| all 13 | 7.43 / 6.98 / 6.57 | 7.01 / 6.60 / 6.22 |
| H&G 9 | 7.29 / 6.86 / 6.46 | 7.12 / 6.71 / 6.33 |
| PGC, VGZ, NLLB, SNB | 6.31 / 6.08 / 5.85 | 6.06 / 5.84 / 5.61 |
| + GNW | 6.56 / 6.29 / 6.02 | 6.30 / 6.05 / 5.79 |

Results at the 0.5-per-template-day threshold (battery with e′; distinct events; scrambled = mean of 16 catalogues):

| set | leave-one-out re-found (of 8) | ETS real | between real | scrambled between per yr (max) | noise templates per template-yr | rule |
|---|---|---|---|---|---|---|
| all 13 | **6** | 49 | 1 | 0.13 (1) | 0.12 | fail (≥ 7) |
| H&G 9 | **7** | 45 | 1 | 0.13 (1) | 0.25 | pass (column) |
| PGC, VGZ, NLLB, SNB | **2** | 21 | 8 | 1.16 (3) | 0.21 | **fail** (≥ 5 and ≤ 1/yr) |
| + GNW | **1** | 21 | 6 | 1.16 (3) | 0.16 | fail |

- At 1 per template-day: all 13 6/8, H&G 9 8/8, 4 stations 3/8 (scrambled 2.9/yr), + GNW 2/8. At 0.25: 6, 6, 2, 1.
- Misses on all 13: 2011-08-25 18:25:51 PDT (best 6.19 MAD, H&G T2) and 2012-09-18 05:02:57 PDT (6.19, T2) fall just
  under the 6.98/6.60 thresholds; H&G 9 finds the second with the C5 single-event template (8.3).
- On the 4 Canadian stations the best leave-one-out z is 4–6.3 MAD for 7 of 8 events (one, 2011-08-19 21:37:55 PDT, has no
  peak at all), against a threshold of 5.8–6.1. Events of Mw 3.1–3.5 are too small for 4 stations at 20–50 s.
- On 4 stations, 8 between-ETS events pass the battery at 0.5/day, above the scrambled catalogues (mean 1.1, max 3).
  The scrambled rate alone (1.16/yr) already fails the rule. The events are listed in `tables/e3_real_battery_events.csv`
  (2 "with LFEs"), and none has been checked by eye.
- **Decision by the pre-set rule:** the 4-station test fails → 2003–2009 = Mw ≥ 4.0 tier + Poisson bound; the uniform
  catalogue starts 2010 (pilot 2008 not run). Scope: the test events are Mw 3.1–3.5, below the Mw 3.6–3.8 target, so
  4-station completeness at Mw 3.6–4.0 is untested (the cross-template injection of PLAN_CAT §7 would measure it).
- The all-13 bank misses the ≥ 7/8 rule by one event while H&G 9 meets it. Adding the four CN stations lowers the summed
  CC of distant events (noise-dominated channels).

## 5. E2 between-ETS events in raw waveforms (`figures/e2_wave_*.png`)
- 12 between-ETS events pass the single-template e′ battery at 5.5 MAD; their raw records are plotted (near scan
  stations, H&G network with the LS-scaled template, far references).
- The two-template event 2012-02-13 05:44:40 PST shows no pulse above the background (60–150 nm at near stations),
  and the far references (OMAK, HAWA, COR, NEW, I02D) carry coherent 20–50 s wavetrains of 50–80 nm in the same window.
  It does not resemble the ETS events (largest at VGZ/LRIV, absent far away).

## 6. Far-reference download (approved ~1:20 pm PT)
Running: `../scripts/scan_fetch_ref_all.py`, IU.COR, US.NEW, US.HAWA, US.MSO, US.WVOR LHZ 2003-01-01 – 2026-10-15 →
`../data/raw_ref_all/`; log `../logs/ref_all_download.log`; marker `../REF_ALL_DONE`; stops at 10 GB. No year has < 2
LHZ far stations, so no other data class was fetched.

## 7. Early-era completeness, 4 stations {PGC, VGZ, NLLB, SNB} (Amendment P3; `tables/inj_early_completeness*.csv`)
Written ~2:55 pm PT. Raw figure first: `figures/inj_early_completeness.png` (fraction recovered vs injected Mw, summer and
winter; MF only, MF + battery′ with the 2011–12 scan, MF + battery′ with the early-era scan).
Method: the amplitude-calibrated stacks of clusters C0 (5 events) and C1 (4 events) injected into 2011–12 data at
1,000 random live times (50 per source × Mw × season); detected by E3 cn4 templates that contain none of the source's
events; threshold = 0.5 noise detections per template-day (5.84–6.08 MAD).

| Mw | season | MF ≥ threshold | MF + battery′ (2011–12 scan) | MF + battery′ (early-era scan, rule) | median best MF z | median scan z (early) |
|---|---|---|---|---|---|---|
| 3.4 | summer / winter | 0.02 / 0.04 | 0.01 / 0.03 | 0.00 / 0.00 | 4.6 / 4.7 | 10 / 9 |
| 3.6 | summer / winter | 0.09 / 0.06 | 0.04 / 0.04 | 0.01 / 0.00 | 4.9 / 4.8 | 18 / 18 |
| 3.8 | summer / winter | 0.21 / 0.15 | 0.14 / 0.13 | 0.00 / 0.01 | 5.2 / 5.1 | 39 / 33 |
| 4.0 | summer / winter | 0.24 / 0.26 | 0.13 / 0.21 | 0.00 / 0.01 | 5.4 / 5.4 | 72 / 64 |
| 4.2 | summer / winter | 0.29 / 0.30 | 0.09 / 0.13 | 0.00 / 0.00 | 5.4 / 5.3 | 151 / 123 |

- **Rule result: no tested Mw ≤ 4.2 reaches 80 %** → the 4-station early era (2003–2009) cannot be catalogued with this
  template bank at Mw ≤ 4.2 (PLAN_CAT P3: catalogue only above 4.2 as a bound, or not at all).
- **Why the matched filter stops at ~30 %:** the median best MF z rises only from 4.6 to 5.4 while the amplitude grows
  ×16. Once the injected event dominates the noise, the CC equals the waveform similarity between the injected cluster
  and the other templates. On 4 stations that similarity sums to ~5 MAD, just under the threshold. A larger event does not
  help; a template from the same place would.
- **Why the early-era battery fails even when the MF fires:** the scan without LRIV puts C1 injections on the grid's
  south edge (c passes 28–39 %), and C0 injections become single-station (z_loso / z_best 0.07–0.15, f passes ≤ 9 %). With
  LRIV (2011–12 scan) a, c, f pass 86–100 %.
- At high Mw the plane-wave test (b) also starts to fail (C0: 0.89 → 0.54 at Mw 4.2): the injected far synthetic plus the
  near signal look like one coherent wave across ≥ 250 km. This is an injection-model limit to keep in mind for large events.
- The scan statistic itself sees every injection strongly (median z 9–10 at Mw 3.4, > 30 at Mw ≥ 3.8). For a pre-2010
  catalogue the detector would have to be the scan at high z, not cross-templates on 4 stations. That is a design question
  for Merlin, not tested here.

## 8. Pilot 2014 (Merlin step 6; Amendment P2; `tables/p14_*`)
Written ~3:15 pm PT. Raw figures first: `figures/p14_mf_vs_scanz.png` (MF z vs scan z for every bank detection, real vs
16 scrambled catalogues, ETS / between / other episodes; rings = tier), `figures/p14_completeness.png`.
Tier: H&G 9 stations, V4, 5.5 MAD + two-template agreement + battery′. Bank: H&G T1–T5 (2011), C0–C5, H&G T6–T13 (2014).
2014 has many episodes: between-ETS live length is 0.56 yr.

**Verdict: NO-GO** (`tables/p14_verdict.json`).

| GO item | rule | 2014 | pass |
|---|---|---|---|
| truth T6–T13 re-found by a template not containing them | ≥ 5/8 | **0/8** (MF with two templates: 0/8) | no |
| between-ETS real vs scrambled | ≤ mean + 1 | **3** vs 0.44 (max 2) | no |
| scrambled rate | ≤ 1/yr | 0.79/yr | yes |
| completeness at Mw 3.8 (C0/C1 injections) | ≥ 80 % | **81.5 %** | yes |

Completeness (`tables/p14_completeness.csv`): Mw 3.4 / 3.6 / 3.8 / 4.0 / 4.2 → 0.45 / 0.73 / 0.82 / 0.85 / 0.80 (MF two
templates alone 0.61 / 0.85 / 0.97 / 0.99 / 1.00; the battery′ caps it at ~0.8–0.86).

**Why the truth test fails — the H&G 2014 events are not visible in our 2014 data** (`tables/p14_truth_battery_at_origin.csv`):
- At the eight Table-1 origins the scan statistic is 1.5–3.6 (2011 H&G T2/T4: 18–19). No battery test passes c or f.
- Raw 0.02–0.05 Hz vertical displacement in the 150 s after each origin: VGZ 24–51 nm and LRIV 17–74 nm, the same as their
  99th-percentile background in the surrounding ±1 h (25–63 nm). For 2011 T2/T4 the same stations show 103–201 nm pulses
  over a 57–79 nm background.
- The 2014 events lie at 48.65–49.10 N, 20–50 km NW of the scan core. SNB was excluded as noisy (December), YOUB was off
  after 2014-09-03, and LZB was excluded all year.
- **RETRACTED (Amendment Q0, ~3:37 pm PT):** "PGC background 99th percentile ~700–1,000 nm". PGC was at its normal level
  (median 10-min RMS 23–36 nm on the event days; month median Z 27 nm). The 99th percentiles came from transients inside the
  ±1-h windows (79–825 samples above 300 nm), and 5 of 8 origins sit in catalogued-earthquake windows. PGC shows 174 nm at
  T11 (median 22 nm) and 439 nm at T13 (inside a transient), so "no pulse anywhere" is not established for PGC; the scan
  statistic stays 1.5–3.6 at all eight origins.
- The 2014 templates do not find each other either (cross-template MF z 3–6.5 MAD; each self-detects at 20–22).
- So the truth item tests December-2014 station noise and geometry, not the method. The 2011–12 analogue (E3, events in
  the network core) re-found 6–7 of 8.
- **Between ETS, 3 events pass the tier** (all "quiet", no LFEs within ±10 min): 2014-01-23 10:59–11:00 PST (T3, C0 and T5
  agree; inside a M 4.6 head2 window), 2014-01-23 11:55 PST (T5 + C5), 2014-08-14 14:14 PDT (C1 + T13). Against 0.44 per
  scrambled catalogue they exceed the rule. Not yet looked at in waveforms.

## 9. Amendment Q0 fixes, 2014 re-summary and Q2 triage (written ~3:38 pm PT)
- Fixes: `hg_r3.lfe_lists()` window extended to 2027 (it clipped all LFE/tremor labels after 2012); US.MSO and US.WVOR
  prepped for 2014 and used as far references (5 far LHZ); validated bank only (T1–T5 + C0–C5); GO rule without the
  between-ETS gate; PGC statement in §8 retracted.
- **2014 re-summary** (`tables/p14_verdict.json`, `p14_counts.csv`; pre-Q0 versions kept as `*_v1_preQ0*`):
  - ETS 0, between ETS 3, other episodes 1 (Merlin expected 0 / 2 / 1). Scrambled 0.13 between ETS (0.22/yr).
  - Truth 0/8: H&G T6–T13 are not observable in our data (§8).
  - Completeness at Mw 3.4 / 3.6 / 3.8 / 4.0 / 4.2 = 0.42 / 0.70 / 0.80 / 0.81 / 0.82.
  - GO items: truth no, controls yes, completeness yes (exactly 0.80) → **NO-GO** (by the truth item only).
- **Q2 triage** of the 3 between-ETS events (`tables/q2_triage.csv`, figures `figures/q2_wave_*.png`, near stations + GNW
  + 5 far LHZ):

| time (PT) | templates | MF z | scan z | vr_thrust | far/near max | ratio_GNW | LFEs ±10 min | class |
|---|---|---|---|---|---|---|---|---|
| 2014-01-23 11:00:25 PST | T5 + C0 | 6.4 | 6.3 | 0.07 | 0.44 | 1.07 | 0 | likely noise |
| 2014-01-23 11:49:37 PST | C0 + C2, T3 | 6.2 | 7.1 | 0.11 | 0.57 | 0.81 | 0 | likely noise |
| 2014-01-23 11:55:37 PST | T5 + C5 | 6.0 | 5.7 | 0.07 | 0.30 | 0.93 | 0 | likely noise |

  All three fall within one hour on 2014-01-23, all fail vr_thrust ≥ 0.12 and the GNW test, and one also has far references
  at > 0.5 of the near amplitude. All are "likely noise", as Merlin expected. vr_thrust ≥ 0.12 is now a column in the pilots.
- **Sensor check against Merlin's density addendum:**
  - The CAFE stations BS11, W020, W040, PL11 and W060 carried 1-s CMG40T-1 sensors from mid-2011 to 2012-05 (response
    at 30 s = 0.16 % of 1 Hz). So the 2011 templates cannot be cut on them; they were broadband (CMG3T) in 2006–2010 only.
  - UW.SQM is broadband only until 2011-07-02, UW.HDW only from 2023, and UW.OTR only from 2020-10.
  - Q1′ sets were redefined by rule (dense2011 = 20 stations).

## 10. Q1 — 2011–12 leave-one-out re-test, plain vs weighted (`tables/q1_summary.csv`, `q1_refind.csv`, `q1_verdict.json`)
Raw figure: `figures/q1_loo.png` (best leave-one-out MF z minus the 0.5/day threshold, per vetted event and set).

| set (stations) | LOO plain | LOO weighted | ETS real (weighted) | scrambled between /yr (weighted) | rule |
|---|---|---|---|---|---|
| dense2011 (20) | 8/8 | **8/8** | 67 | 0.26 | **pass** (8/8, ≤ 1/yr) |
| hg9 (9) | 7/8 | **8/8** | 51 | 0.06 | pass (≥ 7/8) |
| era2016 (11) | 4/8 | 6/8 | 50 | 0.58 | fail (< 7/8) |
| era2010 (7) | 4/8 | 5/8 | 41 | 1.61 | (column) |
| era2007 (5) | 1/8 | 2/8 | 24 | 0.84 | **fail** (< 5/8) |
| cn4 (4) | 2/8 | 3/8 | 19 | 1.80 | control |

- **The weighted sum is ≥ the plain sum on every set**, so later stages use it. dense2011 lifts the best leave-one-out z to
  7.2–13.9 (weighted) against a ~7.0 threshold, against 3–7 on 4–5 stations: the √N gain Merlin predicted.
- era2007 (the 2011-template carriers that also exist in 2007) fails, so 2007–2010 do not enter with the 2011 bank. ETS
  counts here include the templates' own events.

## 11. Pilots under Amendment Q: 2013-09 NO-GO, 2008-05 non-informative (`tables/pilot2013_verdict_Q.json`, `pilot2008_verdict_Q.json`)
- **2013-09**, 14 stations, weighted:
  - T1 pass: 16 ETS tier events vs 0 in 16 scrambled catalogues and 0 noise pairs (Poisson p ≈ 0).
  - T2 fail: 94 % had LFEs, but the same-day LFE rate ratio was only 1.36.
  - T3 fail: the largest amplitude was at PGC, not at the 2 nearest stations.
  - Controls pass: 0.54/yr.
  - Completeness 0.42 / 0.69 / 0.78 at Mw 3.4 / 3.6 / 3.8 (n = 200), just under 0.4 / 0.7 / 0.8.
  - Between ETS: 4 events. Raw figure: `figures/pilot2013_mf_vs_scanz.png`.
- **2008-05** (era2007, 5 stations): 0 detections in every period and 0 % completeness. This is not informative:
  - The 2008 scan had only PGC N/Z active. VGZ and SNB were excluded in all 12 months by the noise rule, which is
    referenced to 2011 levels; PFB and YOUB have no 2008 data on disk.
  - The 5-station MF is shape-limited (median best z 4.3–5.2).

## 12. 2013 re-summarised under Amendment R (first pass on the n = 200 injections; `tables/pilot2013_R_*_n200`)
- T1′ pass: 16 vs λ from the between-ETS scrambled rate (p ≈ 0).
- **T2′ pass, reproducing Merlin:** T2a 0.94; T2b patch coincidence 48 LFEs vs 21.0 null mean (99th pct 41), p < 0.001.
- **T3′ 3/3, reproducing Merlin** (`tables/pilot2013_R_T3.csv`), top 3 events 2013-10-03 18:44, 18:48 and 2013-10-04 04:14 PDT:
  - near/mid median ratio 2.5–2.6;
  - far median ≤ 0.24 × the near maximum;
  - Spearman(R, amplitude) −0.53 to −0.84.
- Controls pass, 0.54/yr. Under R the 2013 items T1′, T2′, T3′ and controls all pass.
- Completeness (n = 200; n = 800 re-run queued) at Mw 3.8: ETS 0.84 (0.76–0.90), between 0.72 (0.63–0.80).
  - MF alone ≥ 0.96; battery′ given MF 0.86 (ETS) / 0.75 (between), so the battery is the ceiling.
  - Mw_c: ETS 3.8, between > 4.0.
- **Between ETS, 4 events, Tier A′ = 0** (`tables/pilot2013_R_between.csv`): all four fail Q2 (vr_thrust 0.04–0.12;
  ratio_GNW 0.42–1.21).
  - Three are on 2013-01-18 (cluster day; LFE status "candidates/tremor only" for two). One is on 2013-02-08.
  - My storm proxy (PGC Z 5–10 s, median 10-min RMS) puts 2013-01-18 at the 52nd percentile and 2013-02-08 at the 88th.
    This differs from Merlin's 88–94th for 01-18; the two proxies are not the same.

## 13. Production 2011 (Amendment S/T; scan v2, 2011 bank, weighted; `tables/prod2011_*_v2.*`) — checked ~9:30 pm PT
Network: 13 MF stations (PGC, VGZ, SNB, NLLB, PFB, LRIV, DOSE, GNW, SP2, A04D, D03D, WISH, NLWA).
- **Known answers: 7/7 of the 2011 vetted + near-miss events found at tier level** (best MF z 11–24; the other 7 of the 14
  are in 2012). Passes the T1 check, so the queue continues.
- 2011-08 episode (major): 41 tier events (38 with LFEs), summed scan M0 4.0e15 N m (Mw 4.3 equivalent). In H&G's
  2011-ETS window: 42, against 8 in E2's two-template e′ tier on the 9-station V4 network. The difference is the denser,
  weighted bank, as Q1 showed (dense2011 leave-one-out z 7–14 vs 6–11 on hg9).
- No tier events in the five minor 2011 episodes (scrambled 0–0.13 per episode).
- T1′ p ≈ 4e-82 (λ 0.17). T2′: T2a 0.93; T2b 121 vs 47 (99th pct 80). T3′ 3/3.
- Controls: scrambled 0.69/yr between ETS; noise pairs 0.
- Between ETS: 3 raw events (2011-09-23, 10-02, 12-08), all with vr_thrust ≈ 0.05 → **Tier A′ 0**.
- Completeness (n = 800 per Mw per period; Wilson 95 %):

| Mw | 3.4 | 3.6 | 3.8 | 4.0 |
|---|---|---|---|---|
| episodes | 0.52 (0.48–0.55) | 0.72 (0.69–0.75) | 0.79 (0.76–0.82) | 0.77 (0.74–0.80) |
| between | 0.42 (0.38–0.45) | 0.60 (0.56–0.63) | 0.71 (0.68–0.74) | 0.62 (0.58–0.65) |

  The MF alone recovers 0.97–0.98 at Mw ≥ 3.8; the battery′ given MF (0.63–0.82) is the ceiling. It falls at Mw 4.0
  between ETS, the large-event plane-wave effect seen in §7. Mw_c (≥ 0.8) is not reached (0.79 at Mw 3.8 in episodes).

## 14. Production 2010–2018 (scan v2, 2011 bank, weighted, n = 800 injections per Mw per period) — written 2026-10-16 ~3:20 am PT
**Partly superseded by §16 (2026-10-16 ~6:10 am PT):** the 2010 "with LFEs" counts were void (LFE-list clip, fixed;
2010-08b1 now 8/8 with LFEs), and 2016–2018 are replaced by the relay re-run (bank 2011 ∪ 2018 relay). The §14 numbers
below are kept as they were first reported.
Raw figure first: `figures/prod_episodes.png` (tier events and summed scan moment per SVI episode against time; large dots =
majors; colour = completeness at Mw 3.6; no error bars). Tables: `tables/prod_years.csv`, `prod_episodes_all.csv`,
`prod_ledger_all.csv`, `prod_between_all.csv`, per year `prod<Y>_*_v2.*`.

### Plain answer
- **Strong VLFEs come from the episodes that activate the C0/C1 LFE patches, and they come back in each such episode.**
  - Six episodes have ≥ 100 P ≥ 0.75 LFEs within 10 km of the patches: 2010-08b1, 2011-08, 2012-08, 2013-09, 2015-12 and
    2018-06. They hold 141 of the 163 tier events, and 126 of the 129 events that have LFEs within ±10 min.
  - The other 45 episodes give 22 events in total, 0–2 each and almost none with LFEs. That is about the scrambled-control
    level in the low-completeness years.
- **Between ETS, no real-looking VLFE in 9 years.**
  - 39 raw tier events between ETS, all failing the Q2 rule, so Tier A′ = 0 in every year.
  - Live between-ETS time is 5.5 yr, or 3.3 yr weighted by completeness at Mw 3.8.
  - 95 % Poisson bound: ≤ 0.9 per year (≤ 0.54 per year without the completeness weighting).
  - Scope: events like the C0/C1 patch VLFEs, Mw ≈ 3.8.
- **Known answers all found (`prod_ledger_all.csv`):** 2011–12 vetted 8/8, near-misses 6/6, 2013 pilot ETS events 18/18,
  2018 pilot ETS events 17/17.
- **Completeness is good in 2011–2014 and poor in 2015–2018**, where the matched filter (not the battery) is the limit.
  This is why relay templates are being built now.

### Per year (completeness = tier recovery of C0/C1 injections, Wilson 95 % in brackets; full MF/battery split in `prod_years.csv`)
| year | MF stations | episodes: Mw 3.4 / 3.6 / 3.8 / 4.0 | between: Mw 3.8 | MF-only / battery-given-MF at 3.8 (episodes) | scrambled /yr | between raw / Tier A′ | episode tier events |
|---|---|---|---|---|---|---|---|
| 2010 | 7 | 0.32 / 0.51 / 0.64 / 0.68 | 0.65 | 0.87 / 0.74 | 0.63 | 2 / 0 | 10 |
| 2011 | 13 | 0.52 / 0.72 / 0.79 / 0.77 | 0.71 | 0.97 / 0.82 | 0.69 | 3 / 0 | 41 |
| 2012 | 14 | 0.54 / 0.72 / 0.78 / 0.80 | 0.79 | 0.95 / 0.81 | 0.78 | 5 / 0 | 57 |
| 2013 | 14 | 0.52 / 0.74 / 0.83 / 0.85 | 0.83 | 0.97 / 0.86 | 0.72 | 5 / 0 | 18 |
| 2014 | 14 | 0.46 / 0.70 / 0.76 / 0.81 | 0.80 | 0.96 / 0.78 | **1.01** | 11 / 0 | 1 |
| 2015 | 13 | 0.23 / 0.31 / 0.41 / 0.44 | 0.64 | 0.69 / 0.60 | **1.01** | 7 / 0 | 4 |
| 2016 | 11 | 0.29 / 0.38 / 0.42 / 0.50 | 0.49 | 0.50 / 0.83 | **1.56** | 5 / 0 | 4 |
| 2017 | 11 | 0.16 / 0.20 / 0.25 / 0.35 | 0.19 | 0.34 / 0.72 | **1.21** | 1 / 0 | 10 |
| 2018 | 10 | 0.06 / 0.15 / 0.28 / 0.40 | 0.23 | 0.40 / 0.69 | 0.34 | 0 / 0 | 18 |

- Mw_c (≥ 0.8 in both periods) is reached only in 2013 (3.8); 2012 and 2014 reach it at Mw 4.0. 2015–2018 are flagged for
  relay templates.
- **Scrambled rate above 1/yr in 2014–2017** (1.01–1.56): the controls item fails there, so their few tier events outside
  the patch episodes should be read as noise-level.

### Per episode (`prod_episodes_all.csv`; episodes spanning New Year appear once per year frame)
| episode | major | tier events | with LFEs | summed M0 (Mw) | completeness Mw 3.4 / 3.8 | patch LFEs (all SVI P ≥ 0.75) | corrected count (p ≥ 0.3) |
|---|---|---|---|---|---|---|---|
| 2010-08b1 | yes | 8 | 0 | 1.6e15 (4.07) | 0.41 / 0.69 | 206 (1,474) | 20 |
| 2011-08 | yes | 41 | 38 | 4.0e15 (4.33) | 0.52 / 0.77 | 205 (1,445) | 78 |
| 2012-08 | yes | 54 | 51 | 5.4e15 (4.42) | 0.60 / 0.80 | 362 (2,864) | 90 |
| 2013-09 | yes | 18 | 17 | 2.1e15 (4.15) | 0.55 / 0.89 | 112 (1,359) | 33 |
| 2015-12 | yes | 3 | 3 | 5.2e14 (3.75) | 0.08 / 0.49 | 141 (1,138) | — |
| 2018-06 | yes | 17 | 17 | 2.5e15 (4.20) | 0.04 / 0.33 | 154 (1,498) | — |
| 2017-03 | yes | 8 | 3 | 1.2e15 (3.99) | 0.01 / 0.10 | 8 (42) | — |
| 2014-11 | yes | 0 | 0 | — | 0.34 / 0.80 | 11 (141) | 0 |
| other 43 | 9 majors | 0–2 each (14 total) | 0 | ≤ 3.6 | 0.00–0.78 | 0–115 | — |

- Across the 51 episodes, the tier count follows the patch LFE count (Spearman 0.57, p = 1e-5) and the total SVI LFE count
  (0.56).
- VLFE tier events per 100 patch LFEs: 2011-08 20, 2012-08 15, 2013-09 16, 2018-06 11 (at 4 % completeness for Mw 3.4).
  So each patch-active episode gives roughly one strong VLFE per 5–10 patch LFEs where the detector is complete.
- **2014-11 confirms Merlin's reading:** a major episode that did not reach the patches (11 patch LFEs) gives 0 VLFEs at
  0.80 completeness.
- Open points:
  - 2010-08b1: 8 events, none with an LFE within ±10 min, although the patches were active. Not yet looked at.
  - 2017-03: 8 events at ~1 % completeness and only 8 patch LFEs. The 2017 scrambled rate is 1.2/yr, so these are suspect.

### Scope
- One bank (2011 events, two patch clusters). Completeness is for events that look like C0/C1 at these stations.
- 2015–2018 completeness is low pending the relay banks. The relay queue is running: 2018 started 3:09 am PT; markers
  `RELAY_2018/2017/2016_DONE`.
- The scan moment is a whole-space thrust fit (×2 uncertain).

## 15. Relay templates (Amendment T2) and flagged episodes (written 2026-10-16 ~3:35 am PT)
### Relay step (`tables/prod<Y>_relay_*`)
| year | vetted relay events (ETS tier + LFEs + vr_thrust ≥ 0.12) | clusters | leave-one-out | scrambled /yr | verdict |
|---|---|---|---|---|---|
| 2018 | 17 (2018-07-05 – 07-09 PDT) | 8 (RC0 4, RC1 4, RC2–RC4 2 each, RC5–RC7 single) | **15/17 (0.88)** | 0.34 | **accepted** |
| 2017 | 3 | — | 0/3 | 0.22 | rejected |
| 2016 | 0 | — | — | — | no bank |

- With the relay bank, 2018 has 67 ETS tier events against 17 with the 2011 bank alone. This includes the 17 relay events'
  self-detections; the leave-one-out 15/17 is the honest benchmark.
- **Relay re-run running** (Amendment U): 2016–2018 with bank = 2011 ∪ 2018 relay, completeness from C0/C1 and from the
  relay-cluster injections. Markers `PRODR_2018_DONE`, `PRODR_2017_DONE`, `PRODR_2016_DONE`.

### Flagged episodes (`tables/flag_2017_2017-03.csv`, `flag_2010_2010-08b1.csv`; figures `figures/flag_<year>_<episode>_NN.png`)
- **2010-08b1, 8 events, all "quiet": this was a bug, not a property of the events.**
  - `hg_r3.lfe_lists()` also clipped the lower end at 2011-05-31, so every LFE/tremor label before June 2011 was void.
    Fixed; 2010 and 2011 are being re-scored (`RELABEL_DONE`).
  - The events themselves look like the 2011–13 VLFEs. They sit at the C0/C1 patches (48.07–48.34 N, 123.45–123.58 W) with
    vr_thrust 0.06–0.64 and ratio_far 0.14–0.34, and all pass a, b, c, e′, f.
  - The largest (2010-08-23 03:22:35 PDT, scan z 14.4) shows a clear pulse at LRIV (4 km) and VGZ, with the far references
    at noise level.
- **2017-03, 8 events at ~1 % injection completeness:** they look real, not noise.
  - They sit at the same patches (S_Vancouver_Is:3, Olympic_2010:907 / armB, families 18–28), with ratio_far 0.19–0.50 and
    scan z 5.4–16.6. 3 have LFEs within ±10 min and the others have LFE candidates or tremor.
  - The top event (2017-03-23 01:25:07 PDT, MF 8.5, scan z 13.3) shows a clean pulse at VGZ (26 km), PGC and SNB, decaying
    with distance, and nothing at the 5 far references.
  - Low injection completeness and real detections can coexist: the injected C0/C1 stacks resemble the 2017 records less
    than these events resemble the 2011 templates.
  - Most fail the GNW clause of Q2 (ratio_GNW 1.9–13.7) because GNW was noisy, so Tier A′ would reject them. The Q2 GNW
    clause needs a look for 2017.
- Station problems seen in these figures (affect the far/near column, not the battery): LRIV 2017 (7.2 µm peak) and LZB 2010
  (4.1 µm peak).

### Planned, not run (Amendment V)
- 2019–2026: forward relay chain from B(2018) = 2011 ∪ 2018 relay, with a relay bank accepted per year by the same rule.
- 2007–2009: backward relay via 2010-08b1 with CAFE + SQM broadbands, a scan extension, and pilots 2008-05 / 2009-05 / 2007-01.

## 16. Relay re-run 2016–2018, corrected 2010–11 LFE labels, updated 2010–2018 summary (written 2026-10-16 ~6:10 am PT)
Raw figure first: `figures/prod_episodes_final.png` (tier events and summed moment per SVI episode against time; colour =
completeness at Mw 3.6; large = majors). Tables: `tables/prod_years_final.csv`, `prod_episodes_all_final.csv`,
`prod_ledger_all_final.csv`, `prod_between_all_final.csv`, per year `prodr<Y>_*_v2.*` (relay re-run).

### Relay re-run (Amendment U; bank = 2011 ∪ 2018 relay RC0–RC7)
| year | stations | completeness at Mw 3.8, episodes (all injections / C0C1 / relay) | between (all / C0C1 / relay) | MF-only / battery-given-MF (episodes, all) | scrambled /yr | episode tier events | between raw / Tier A′ |
|---|---|---|---|---|---|---|---|
| 2016 | 11 | 0.38 / 0.61 / 0.15 | 0.40 / 0.64 / 0.17 | 0.84 / 0.46 | **3.72** | 11 | 11 / 1 |
| 2017 | 11 | 0.42 / 0.48 / 0.36 | 0.42 / 0.44 / 0.40 | 0.80 / 0.53 | **2.65** | 34 | 5 / 0 |
| 2018 | 10 | 0.64 / 0.56 / 0.72 | 0.56 / 0.50 / 0.62 | 0.84 / 0.76 | 0.67 | 66 | 3 / 0 |

(n = 800 per Mw per period per injection set; Wilson intervals ±0.02–0.03; full table `prodr<Y>_R_completeness_v2.csv`.)
- **2018: the relay bank works.**
  - C0/C1 completeness at Mw 3.8 doubles, 0.28 → 0.56 (episodes), and 0.23 → 0.50 between.
  - Relay-cluster injections are recovered at 0.72 / 0.62.
  - Controls 0.67/yr pass. Mw_c is not reached (> 4.0).
  - 2018-06: 65 tier events (58 with LFEs). 17 of them are the relay events themselves (self-detections) and 48 are new;
    all 17 events of the 2011-bank run are among them. The ETS-window scrambled mean is 3.25 (one draw 27) against 65.
- **2016–2017: the relay bank makes things worse.** The matched filter fires more (MF-only 0.80–0.84), but the scrambled
  rate rises to 2.6–3.7/yr and battery′-given-MF drops to 0.46–0.53. The relay templates (2018 windows) do not suit the
  2016–17 networks.
  - By the controls gate, 2016–2017 fail with either bank (2011 bank: 1.56 and 1.21/yr).
  - 2017-03: 31 tier events (6 with LFEs) at the patches vs 8 with the 2011 bank; the ETS-window scrambled mean is only 0.44.
  - 2015-12 (Jan 2016 part): 7 tier events, all with LFEs.
- Between ETS: one Tier A′ event, 2016-02-25 12:11:21 PST (vr_thrust 0.17, cluster day, no LFEs, far/near 0.006). It comes
  from a year whose scrambled rate is 3.7/yr; read it as noise until looked at.

### Corrected LFE labels (relabel after the second lfe_lists fix)
- 2010: T2a 0.0 → **0.80** (T2′ now passes; T2b was always 33 vs 8.2). 2010-08b1: "with LFEs" 0 → **8/8**.
- 2011: unchanged (41 tier events in 2011-08, 38 with LFEs). 2012–2018 used lists that were never clipped.

### Updated 2010–2018 summary (2010–2015 with the 2011 bank, 2016–2018 relay re-run)
- **242 tier events in 51 episodes (54 year-frame rows); 189 with LFEs within ±10 min.**
- **196 of them (182 with LFEs) are in the six patch-active episodes.** The other 45 episodes give 46 events (7 with LFEs);
  31 of those are 2017-03.

| patch episode | tier events (with LFEs) | summed M0 (Mw) | patch LFEs | VLFE per 100 patch LFEs | completeness Mw 3.4 / 3.8 |
|---|---|---|---|---|---|
| 2010-08b1 | 8 (8) | 1.6e15 (4.07) | 206 | 3.9 | 0.41 / 0.69 |
| 2011-08 | 41 (38) | 4.0e15 (4.33) | 205 | 20 | 0.52 / 0.77 |
| 2012-08 | 54 (51) | 5.4e15 (4.42) | 362 | 15 | 0.60 / 0.80 |
| 2013-09 | 18 (17) | 2.1e15 (4.15) | 112 | 16 | 0.55 / 0.89 |
| 2015-12 (both frames) | 10 (10) | 1.3e15 (4.02) | 141 | 7.1 | 0.14 / 0.49 |
| 2018-06 | 65 (58); 48 without the relay self-detections | 8.8e15 (4.56) | 154 | 42 (31 without self) | 0.42 / 0.65 |

- Tier count vs patch LFE count over all episodes: Spearman 0.55 (p = 2e-5). 2014-11 (11 patch LFEs, 0 events at 0.80
  completeness) stays the counter-example of a major episode that did not reach the patches.
- **Between ETS, 2010–2018:**
  - 1 Tier A′ event in 5.5 live years (2016, likely noise, see above); 0 if 2016 is excluded for its control failure.
  - 95 % Poisson upper bound with 1 event: ≤ 0.86/yr raw, ≤ 1.3/yr weighted by between-ETS completeness at Mw 3.8
    (3.6 effective years). With 0 events: ≤ 0.84/yr weighted.
  - Scope: C0/C1-like events; the Tier A′ completeness itself was not injected.
- Known answers 49/49 (unchanged).
- **What changed from §14:**
  - The total grew from 163 to 242: the relay re-run added 2016 +7, 2017 +24, 2018 +48.
  - Events with LFEs: 129 → 189.
  - 2010-08b1 now carries LFEs.
  - The 2018-06 moment rose from Mw 4.20 to 4.56.
  - The between-ETS Tier A′ count went from 0 to 1, from the 2016 relay run, where the controls fail.

## 17. 2016–2017 control failure and the 2016-02-25 candidate (Amendment W; written 2026-10-16 ~6:15 am PT)
### Diagnosis (`tables/diag_noise_<Y>.csv`)
- **UW.LRIV was broken through much of 2016–2017.** Its Z median 10-min velocity RMS was 420–480 nm/s (normal 2–5) in
  2016 Jan–Apr and Dec, and in 2017 Jan–Jun and Oct–Dec. LRIV is one of the most heavily weighted near stations in the bank.
- The scrambled-control tier detections (2011 bank) bunch into a few days. Most fall in LRIV-broken months: 8 of 10 events
  in 2017 and 6 of 13 in 2016 (2016-02-25 alone has 4).
- Polarity-draw leakage does not explain them: control counts per draw do not follow the draw's sign correlation with the
  true template.
- The other stations show only normal seasonal noise; GNW is normal at 20–50 s. The fetched/linked carriers are fine.
  PFB and SNB are present in 2016–17; YOUB is off (ended 2014-09).
- Year-own thresholds (6.3–7.2 MAD) do not separate 2016–17 from good years.
- Fix (W2b, defined after the diagnostics and disclosed as such): drop station-component-months whose median RMS is above
  3× that channel's 2011–2018 median, in all MF parts (data, controls, noise templates, injections). Re-runs of 2013
  (check), 2016 and 2017 (both banks) are queued: `scripts/w_queue.sh`, markers `W_prodq2013_DONE`, `W_prodrq2016_DONE`,
  `W_prodrq2017_DONE`, `W_prodq2016_DONE`, `W_prodq2017_DONE`, `W_DONE`.

### 2016-02-25 12:11:21 PST (`tables/flag_2016_2016-02-25_1211.csv`, `figures/flag_2016_2016-02-25_1211_00.png`)
- Detected by RC4 + RC5 (relay templates), MF 5.8 MAD, scan z 7.5. Best node 47.98 N, 123.18 W (south of 48.1 N, so the GNW
  clause does not apply; nearest family Royer_NW_WA:2388).
- vr_thrust 0.17, ratio_far 0.51 (5 refs), passes a, b, c, e′, f, S5.
- No SVI LFE or PNSN tremor within ±1 h. Storm percentile 83 (just under the 85 % storm line). 4 scrambled-control tier
  events on the same day.
- **Raw waveforms: noise.** A continuous coherent long-period wavefield of ~85–125 nm sits on every near and regional
  station, with 42–62 nm at the far references, and no pulse.
- By the Q2 rule it is mechanically "real-looking" only because the nearest station, LRIV, is broken (10.4 µm), which
  makes far/near 0.01. With the next-nearest working station (DOSE, 116 nm) far/near is 0.49, still just under 0.5.
- **Classified likely noise.** This also shows that Q2's far/near clause needs broken stations removed (W2b does this).

## 18. Why 2016–2017 still fail the scrambled-control rate with W2b (diagnostics only; written 2026-10-16 ~11:35 am PT)
Tables `tables/wdiag_<run>.csv` (`scripts/cat_wdiag.py`): every between-ETS scrambled-control tier detection, with the
per-station share of the weighted MF sum recomputed on the scrambled data.
- **W2b works as intended.** In the excluded LRIV months the LRIV CC is 0 and the channel is not live (2016: 42 % of the
  year, 2017: 75 %). Noise templates that fall in excluded months have no LRIV channel (30 and 51 cases). Injections read
  the same masked data.
- **No other broken channel.** In the months W2b keeps, LRIV's noise and transient tail match VGZ and DOSE
  (99th/50th percentile of 10-min RMS 59–252 vs 57–165). No other station-month is outside the 3× rule.
- **The control events sit on ordinary coherent-wavefield days and are led by the most heavily weighted near station,**
  just as in 2013:
  - LRIV leads 13 of 21 control rows in 2017 and 5 of 13 in 2016; in 2013, which passes, it leads 7 of 14.
  - The leading station carries a median 27–30 % of the sum in every year.
  - Weights are not unusually concentrated: N_eff 24–28 of 30 channels for both bank and relay templates; LRIV takes
    13–21 % of the 2016–17 weight vs 11–16 % in 2013.
- **The station set is smaller.** 2016–17 have 11 MF stations against 14 in 2013 (TA.A04D, TA.D03D and CN.LZB are
  gone). 2018 has 10 but passes, so station count alone does not decide it.
- **2011 bank: the "fail" is not significant.**
  - 2017: 1.21/yr means about 11 control events over 16 catalogues × 0.57 yr, against 9.1 expected at 1/yr
    (Poisson p ≈ 0.3).
  - 2016 passes (0.96/yr), but T1 fails there (3 ETS tier events).
- **Relay bank: the excess comes from adding templates.**
  - In 2017 the 2011-bank pairs give the same 21 control rows as without relay; the relay templates add another 21
    (12 relay+relay, 9 mixed).
  - In 2016: 27 rows involve a relay template, 26 do not.
  - The relay templates are not near-duplicates (pair CC ≤ 0.41) and their weights are spread out (top station VGZ
    13–16 %).
  - So the fixed 5.5-MAD, two-template rule simply lets more false pairs through as the bank grows from 11 to 19
    templates. This bites when the scan/battery also passes more often (2016–17). In 2018 the same relay bank stays at
    0.67/yr.
- Proposed fix (to Merlin, not applied): see the hand-back. It brings back Merlin's original "threshold from controls per
  template-year", with held-out controls for the check.

## 19. Merlin review 4 actions (Amendment Y; written 2026-10-16 ~11:40 am PT; Y1, Y3, Y4 running)
- **Y2: final 2016–2017 numbers are the 2011-bank W2b runs (prodq2016, prodq2017).**
  - 2016: T1 fails (3 ETS tier events); controls 0.96/yr.
  - 2017: T1 and T3 pass; controls 1.21/yr.
  - Both have Mw_c > 4.0. The relay-bank runs (scrambled 3.48 and 2.54/yr) are columns only; they fail the T2 relay rule
    in those years.
- **Y5: the 2017 T2 failure is an LFE-catalogue gap, so 2017 T2 is "untestable", not failed** (`tables/y5_tremor_patch.csv`).
  - During 2017-03, PNSN lists 2,506 tremor epicentres within 30 km of the C0/C1 centroids, the same level as the patch
    majors (2011-08 2,713; 2012-08 2,643; 2018-06 2,459). Yet the SVI catalogue has only 42 P ≥ 0.75 LFEs in the whole
    episode (8 at the patches).
  - Note against an earlier reading: 2014-11 also had 1,395 tremor epicentres within 30 km of the patches, although it
    gave 0 VLFEs and 11 patch LFEs. "2014-11 never reached the patches" therefore rests on the LFE catalogue alone.
- **Y6 (optional): LRIV's 2016–18 fault is a long-period sensor fault, not a gain error.**
  - Broken/normal amplitude ratio on Z: 435 at 200 s, 285 at 100 s, 203 at 50 s, 79 at 33 s, 10 at 20 s, 4–7 at 2.5–10 s.
  - N and E are only 1–7× higher.
  - So the excess is red and W2b is needed; a CC-harmless gain error would give a flat ratio.
- **Y1: 2018 re-run under W2b** (2011 bank, then a relay bank re-cut without LRIV, then 2011 ∪ new relay if
  LOO ≥ 85 % and scrambled ≤ 1/yr) is running: `scripts/y_queue.sh`, marker `Y_2018W2B_DONE`.
- **Forward chain: held** at 11:27 am PT. The 2019 outputs built with the old 2018 relay bank are kept as `*_preW2bBank*`
  (`data/prodf2019_preW2bBank/`, `tables/prodf2019_*_preW2bBank.*`, `PRODF_2019_DONE.preW2bBank`), to be redone from the
  new bank. The scan-year 2020 build continues (it does not depend on the bank).
- **Y3/Y4 (pooled controls, dispersion test, Tier A′ on controls)** run after Y1: `scripts/y_checks.sh`, marker
  `Y_CHECKS_DONE`.
  - A code dry run on 2013/2016/2017 alone (not the result): 27 control events, dispersion p = 0.53, 0 control Tier A′
    events.
- **Y7 (user decision, 10-16 ~12:05 pm PT):** the primary false-detection measure is now the noise-template control (H&G
  2019 style). The scrambled controls are columns only; relay acceptance = LOO ≥ 85 %, with the noise-template rate reported.
  - y_checks now reports the noise-template raw tier and Tier A′ rates per year and pooled (`tables/y7_noise_template_verdict.json`,
    `y3_controls.csv`), with the scrambled rates beside them. The between-ETS bound will be stated against the
    noise-template Tier A′ rate.
  - Dry run on 2013/2016/2017 (code check only): 2 noise-template tier events in 43.9 noise-template years (0.05/yr), 0 Tier A′.
- **Y1 result: 2018 under W2b (Y_2018W2B_DONE 1:06 pm PT).**
  - 2011 bank (prodq2018), LRIV excluded all year:
    - T1 20 ETS tier events (pass); T3′ 3/3.
    - T2a 0.90, but T2b 35 vs 17.5 is not above the 99th percentile, so T2′ does not pass.
    - Controls: scrambled 0.11/yr, noise-template 0.
    - Completeness at Mw 3.8: 0.34 (episodes) / 0.27 (between). Between ETS: 0 raw. Ledger 14/17.
  - New 2018 relay bank (cut without LRIV): 18 vetted events, LOO 18/18, scrambled 0.22/yr, noise-template tier rate 0
    → accepted (Y7 rule).
  - 2011 ∪ new relay (prodrq2018):
    - T1 83 ETS tier events (pass; the pre-W2b relay run had 66). T2′ passes (T2a 0.87; T2b 131 vs 53).
    - T3′ 2/3, still a pass (≥ 2 of 3).
    - Controls: scrambled 0.34/yr, noise-template 0. Between ETS: 2 raw, Tier A′ 0. Ledger 16/17.
    - Completeness at Mw 3.8: 0.68 (episodes) / 0.60 (between); C0/C1 injections 0.63 / 0.51; relay injections 0.74 / 0.70.
  - Accept rule met: T1 and T3 verdicts unchanged. The LRIV channel was not carrying the 2018 result.
  - Forward chain restarted 1:08 pm PT with bank 2011 ∪ W2b-2018 relay, capped at the latest accepted relay. It redoes 2019
    first.

## 20. DRAFT — ranked catalogue (Y8b/Y8c): why between-ETS real rows exceed the noise templates (written 2026-10-16 ~2:40 pm PT)
Scripts `scripts/cat_ranked2.py` (Y8b), `scripts/cat_ranked3.py` (Y8c + diagnostics). Times PT. Between-ETS live time
2010–2018 = 5.54 yr; noise templates 23–25 per year (136.3 template-years). Bank = 11 templates (T1–T5, C0–C5) every year.

### Y8c score fix
- Y8b was separable (26 positives vs ~70k noise-template rows): intercept −98.7, 1,088 real rows at p ≥ 0.9, many at 1.0.
- Y8c keeps the monotone logistic, picks L2 by 5-fold CV, and uses the **logit as the ordering score** (p is a column).
  CV still prefers the weakest L2 (held-out balanced log-loss 0.007; AUC 0.9999 at every L2), so the classes stay separable;
  the logit nevertheless spreads (100 distinct values in the top 100; 0.2 % of real rows at p ≥ 0.999).
- Calibration (`tables/ranked3_calibration.csv`): the probability is not a probability (all 26 positives in the top
  out-of-fold decile, mean p there 0.03, because of class balancing). The noise-template false rate is overdispersed between
  templates: at logit ≥ 3.6, even templates 0.54 and odd templates 1.20 per template-year; at ≥ 4.8, 0.14 vs 0.41.
- Known VLFEs: median rank percentile 99.89 %, lowest 98.98 % (in-sample); the 6 near-misses rank at 99.99 %.

### (a) Normalisation: same footing
- `noise_template_per_catalogue` was ONE noise template's catalogue, compared with a real catalogue of 11 templates.
  Worse, a noise template has one agreement partner (N_k with N_k+1) while each real template has 10, and the real templates
  resemble one another, so they agree far more easily.
- Raw MF ≥ 5.5 between ETS (`tables/r20_footing_raw.csv`, no battery): **single-template rates are the same** (real
  790–1,400 vs noise 880–1,490 detections per template-year; 2010 noise 4,420). **Agreement is not:** the real bank gives
  330–1,120 agreement events per year, a bank of 11 random noise templates 80–450 (20 random banks per year; 2010 and
  2012 are the exceptions, 820 and 650).
- Y8c score, between ETS, per bank-year (11 templates):

  | logit ≥ | real | noise template × 11 | scrambled A1 (day shift) | scrambled A2 (polarity) |
  |---|---|---|---|---|
  | 3.6 (top 428 real) | 77 | 9.4 | 10.7 | 21 |
  | 5.9 (top 69) | 12.5 | 1.0 | 1.1 | 2.9 |
  | 8.8 (top 16) | 2.9 | 0 | 0.4 | 0.7 |

  So on any footing there are 4–8× more real rows than controls at high score. The rest of this section shows what they are.

### (b) Why Y8b tierA gives 24 between-ETS events and production Tier A′ gives 0–1
- Production Tier A′ = Q2 applied to **tier** events, and tier requires battery′ (a far test, b not travelling, c located,
  e′ quake test, f z_loso). The final 2010–2018 runs have 45 between-ETS tier events and **0** Tier A′ (the single event,
  2016-02-25, was in a superseded run and was noise).
- Y8b `tierA` = Q2 on any agreement row, with **no battery′**. Of all 80 between-ETS tierA rows, 0 pass battery′ (a 22,
  b 65, c 30, e 25, e′ 43, f 60). Of the 24 at Y8b score ≥ 0.9: a 20, b 23, **c 8, f 8**, e′ 18; battery′ 0; storm day 1,
  cluster day 0; 14 on a grid-edge node; 20 with no LFEs (`tables/r20_condition_matrix.csv`). Agreement, MF threshold
  (5.5) and the storm/cluster rule are the same; the difference is battery′, mainly c (location) and f (leave-one-station-out).

### (c) The 24 between-ETS tierA rows (`tables/r20_tierA_between.csv`, all columns; figures `r20_wave_0..4.png`)
- The top 5 by Y8c logit all fail f (z_loso / z_best −0.10 to 0.09): the detection depends on one station.
- **Raw waveforms (±30 min, 0.02–0.05 Hz):** 4 of 5 are a **single-station transient on UW.LRIV** (77–298 nm in the window,
  while every other near, regional and far station sits at its 10–35 nm background): 2016-07-22 22:13, 2010-06-25 21:42,
  2016-07-04 20:36, 2015-06-21 21:44. LRIV is the nearest station (4–32 km), so far/near (Q2) passes mechanically.
  The fifth, 2015-01-28 16:28, is a **CN.PGC spike** 30 s before t0 (404 nm) on a coherent 30–70 nm wavefield seen equally at
  the far references (far/near 0.5). None shows a pulse that decays with distance. All 5 classified noise.
- Earthquakes beyond the mask: the largest event in the 6 h before each of the 24 is M4.5–6.4 at 900–18,700 km, 20–330 min
  earlier; one row has an M ≥ 6.5 within 24 h. No row is explained by a teleseism. Storm percentile (PGC 5–10 s) ≥ 0.7 for 3
  of 24.
- Across all high-score rows the single-station share is similar for real and controls (logit ≥ 3.6: f fails for 74 % of
  real, 92 % of noise-template and 80 % of scrambled rows), so the score's z_loso weight (0.28) is too weak to push these down.

### (d) The top days (PT): 2014-01-23 (13 rows at Y8b ≥ 0.9), 2013-01-18 (11), 2014-09-14 (7), 2011-10-04 (6)
- Noise-template rows on those days: 155, 160, 50, 135 (vs ~90 on a normal day), but only 1, 0, 0, 1 with agreement and
  **none** at logit ≥ 3.6; scrambled 0–0.5 per draw at logit ≥ 3.6.
- **Raw waveforms of the best f-passing rows on 2014-01-23 and 2013-01-18** (`figures/r20_topday_wave_0..3.png`): an
  energetic long-period wavefield on every station (hour RMS 50–105 nm near, 30–55 nm far, vs 5–20 nm on quiet days). The
  window peak equals the pre-window peak at each station, and far/near is ~0.4–0.5. No pulse: coherent wavefield, not a
  local source. CN.LZB is very noisy those days (hour RMS 230–330 nm).
- The storm metric misses this: 2014-01-23 has storm percentile 0.82 (PGC 5–10 s 434 nm), but 2013-01-18 only 0.52. The
  contaminating energy is at 20–50 s, not 5–10 s.
- No M ≥ 6 earthquake on any of the four days. Tremor/LFE activity: 2014-01-23 has ~220 PNSN tremor in ±3 days (none that
  day); 2011-10-04 has 146 SVI LFEs (c010 list) that day; 2013-01-18 13. So two of the top days sit inside small inter-ETS
  tremor/LFE bursts. That is worth a look, but the waveforms we checked are wavefield-dominated.

### Reading
- The between-ETS "excess" is real rows that the noise-template control cannot reproduce on its own footing: single-station
  transients (mostly LRIV, which the W2b month rule keeps) and days with a strong 20–50 s wavefield that the 11 similar real
  templates match together. Battery′ removes all of them (0 of 80 pass), which is why production has 0 Tier A′ between ETS.
- The ranked product is therefore only usable with battery′ items as strong features or as columns the reader filters on.
  Proposals for Merlin (not applied): (1) a fair noise control for agreement — noise "banks" of 11 templates with the same
  inter-template similarity as the real bank (or scrambled A2, which keeps the wavefield); (2) treat f (z_loso) and the
  max station share as missing-is-worst features with weight not learned from 26 positives; (3) a 20–50 s wavefield-level
  column (hour RMS at the near network) next to the 5–10 s storm percentile; (4) a short-transient (spike) column per
  station.

## 21. Amendment Y9 (Merlin review 5): transparent ladder, Q2 completeness, level-5 between-ETS rows by eye (written 2026-10-16 ~3:00 pm PT; step 4 added ~5:08 pm PT)
Scripts `scripts/cat_y9.py` (steps 1, 2, 3, 5) and `scripts/cat_y9nb.py` (step 4). Final runs 2010–2018 as in §19
(2013/2016–2018 = prodq). Between-ETS live time 5.54 yr. Tier A′/Q2 is retired as a gate and as the bound carrier (Y9);
its quantities stay as columns. `score_logit_insample` (Y8c features, loso fill now −1) is a column only.

### Step 1 — Q2 completeness from injections (`tables/y9_q2_completeness.csv`)
- Of the battery′-recovered injections, 93–100 % have vr_thrust ≥ 0.12 (all years, both periods; Mw 3.6 median 0.98,
  3.8 0.99, 4.0 1.00; lowest 2018, 0.93–0.94). So, by Merlin's reading, the vr clause is not what removes injections.
- But on the 26 known VLFEs (`tables/y9_q2_known.json`): Q2 passes 3; vr ≥ 0.12 16; far/near2 < 0.5 20; the GNW clause 11.
  The injected stack waveforms fit the thrust mechanism better than real weak events do, so injections overstate vr
  completeness. The GNW clause (north of 48.1 N) is the largest single loss. Either way Q2 removes most known VLFEs and
  is retired as a gate, as Y9 pre-registered.

### Step 2 — the ladder (`tables/catalogue_ladder_<Y>.csv`, `catalogue_ladder_2010_2026.csv`; checks `tables/y9_checks.json`)
- Level = battery′ checks passed among agreement rows (a far test, b not travelling, c located, e′ quake test, f
  leave-one-station-out); an unmeasurable check counts as not passed (`n_checks_missing`). Order: level, then MF z.
  Flag columns: eq_in_window, grid_edge, storm_pct_max_pm1d, plus every battery, Q2, emergence and T3 quantity.
- Mechanical checks:
  - Grid edge in the top 50 between ETS: 0 (pass).
  - **head2 earthquake windows in the top 50 between ETS: 24 (fail).** 21 of the 39 level-5 rows sit in head2 windows; they
    pass e′ because the event is below M 5.5 and fewer than 2 far references are excited. Ordering by (level, eq flag,
    MF z) does not change the top 50 (21 there), because only 18 level-5 rows are outside head2 windows.
  - **Known VLFEs at the top level allowed: 21 of 26 (fail).** None has a missing check; 4 fail c (located) and 1 fails b
    (all 2013-10-01/04 ETS events at MF 5.9–8.5).

### Step 3 — ladder tables (`tables/y9_ladder_table.csv`, completeness `tables/y9_ladder_completeness.csv`)
Pooled 2010–2018. Controls: noise template × 11 (one-template catalogue scaled to 11 templates; placeholder until step 4),
scrambled A1 and A2 per draw (2 draws each). Completeness = injections recovered with MF agreement and level ≥ L.

| level ≥ | between: real (days) | noise ×11 | A1 | A2 | compl. Mw 3.8 | ETS: real (days) | noise ×11 | A1 | A2 | compl. 3.8 |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 1,847 (896) | 167 | 259 | 483 | 0.75 | 1,499 (451) | 117 | 107 | 247 | 0.74 |
| 4 | 411 (312) | 24.5 | 36 | 82 | 0.73 | 497 (195) | 13.9 | 15.5 | 54 | 0.72 |
| 5 | **39 (29)** | **1.8** | **2.5** | **4.0** | 0.58 | **146 (51)** | 0.5 | 1.5 | 7.5 | 0.56 |

- Level 5 between ETS: 39 real (95 % Poisson 27.7–53.3) against ≤ 4 control, excess 37 (26–52) = 6.7 per year, or 11.7
  per year after dividing by completeness 0.58 at Mw 3.8. Per year: 2, 3, 4, 4, 10, 7, 5, 4, 0 (2010–2018); 10 of the
  2014 rows are 7 on 2014-01-23.
- Completeness at level 5 is 0.49–0.61 (Mw 3.6–4.0) in both periods; by year at Mw 3.8 it ranges 0.14 (2017) to 0.80 (2014).

### Step 5 — emergence and amplitude–distance (`data/ranked/y9_emerge_<Y>.parquet`, columns in the ladder tables)
- Per station: window peak (0–180 s) / median |Hilbert| envelope over ±30 min. Columns: median over network stations
  ≤ 150 km, number ≥ 2, max/second-max station, far-reference median, and the §12 T3 quantities.
- Calibration on the 26 known VLFEs: median near emergence 1.56–3.36; the threshold that keeps 85 % is **2.04**.
  On 2014-01-23, 4.5 % of all real rows and 14 % of its level-5 rows reach it (Merlin's acceptance: ≤ 20 %; pass).
  Known VLFEs: T3 passes 21/26, median Spearman −0.64, median max/second 1.10 (no single-station dominance).

### The 39 level-5 between-ETS rows by eye (`tables/y9_level5_between.csv`, figures `figures/y9_l5_00..38.png`)
Labelled with the pre-stated rule (spike / wavefield / candidate from emergence), then by eye, with a fourth class that
the rule lacked: **eq-wave** = a coherent arrival that is as large or larger on the far references (far emergence ≥ near),
usually with a catalogued teleseism whose surface-wave train spans the window (all are M 4.5–5.2 at 5,000–16,000 km,
below the M 5.5 e′ cut).
- **Wavefield 22** (7 on 2014-01-23, 2 on 2013-01-18, storm days 2015-12-10/15, earthquake codas, plain noise).
- **Earthquake wave 8**: 2010-03-03 05:14, 2010-03-13 11:02, 2011-09-23 14:49, 2012-06-28 13:26, 2014-07-26 20:41,
  2014-08-19 18:43, 2016-05-24 11:36, 2016-06-22 14:24.
- **Spike 4**: 2013-02-03 12:15 (PGC), 2013-05-13 14:19 (NLLB), 2015-01-22 12:46 (PGC + NLLB), 2015-03-18 05:23 (PGC spike at
  t0, then a regional wave; one of Merlin's three examples). Merlin's other two (2012-07-03 17:08, 2014-07-26 20:41) are
  wavefield and earthquake wave.
- **Candidates 5:**
  - **2016-06-13 07:41:24 PDT:** a clear pulse ~20 s after t0 on all 12 network stations (emergence 2.2–4.4, 44–52 nm),
    weaker at the far references (11–20 nm), Spearman −0.57, no catalogued earthquake. The best candidate.
  - **2015-04-13 06:06:11 PDT:** a sharp arrival on every regional station (emergence 6–10) with far references 60–150 s
    later (5.2), i.e. moveout from near the network; no catalogued event. Could be an uncatalogued regional earthquake;
    it needs a regional-catalogue / high-frequency check.
  - Weak / unclear: 2014-05-05 08:30, 2016-06-13 19:58, 2017-01-31 04:22.
- The rule alone over-calls candidates (20 vs 5 by eye), mainly because it lacks the far-versus-near emergence test.
  The eye labels are not blind (emergence values were on the figures), and control rows were not eye-labelled, so the 34
  explained rows cannot simply be subtracted from a control count.
- CN.PGC often shows isolated spikes on a very low background in 2013–2015 (and CN.LZB is very noisy; UW.LRIV 4,483 nm
  on 2017-01-31). A per-station spike and QC look at PGC and LZB in those years is worth doing.

### Step 4 — similarity-matched noise banks, 2011–2013 (`tables/y9nb_summary.csv`, `tables/y9nb_verdict.json`)
- Real bank mean inter-template CC (max over ±30 s, mean over channels and 55 pairs): 0.583 (2011), 0.579 (2012), 0.585
  (2013). 8 banks per year = base noise template + 10 siblings (mixing weight a 0.57–0.68); every bank matches the target
  CC exactly (`tables/y9nb_<Y>_banks.csv`).
- Pre-registered acceptance (between-ETS agreement rate within ×1.5 of the real bank's): **14 of 24 banks** pass (5 in
  2011, 4 in 2012, 5 in 2013). The agreement ratio spans 0.10–44.6 although every bank has the same mean CC: it follows the
  base window, not the similarity. The extreme bank (2012 N02, ×44.6) has a base that is itself a high-rate template
  (26,600 MF ≥ 5.5 detections per year vs ~1,300 for a typical base such as 2012 N03), and its siblings inherit that
  (6,700–30,400 each); its dominant period (26 s) is ordinary. Low banks (2011 N04 ×0.10, 2012 N01 ×0.27) have quiet bases.
- Level 5 between ETS, same footing (per between-ETS live year; 95 % Poisson):

  | | 2011–13 rate per yr | scaled to 2010–18 (5.54 yr) |
  |---|---|---|
  | real | **11** events → 6.1 (3.1–10.9) | **39** (27.7–53.3) |
  | accepted matched banks (14) | 22 in 8.51 bank-yr → **2.58 (1.62–3.91)** | **14.3 (9.0–21.7)** |
  | all 24 banks (column only) | 120 in 14.39 bank-yr → 8.3 (6.9–10.0) | 46.2 (38.3–55.2) |
  | noise template × 11 (step 3) | | 1.8 |

- **Outcome by the pre-registered rule: intermediate** (accepted-bank control 14.3 is neither ≥ 20 (A) nor ≤ 5 (B)).
  In 2011–13 the real count (11) exceeds the accepted-bank expectation (4.6) with Poisson p = 0.008; over 2010–18 the
  real 39 is 2.7× the scaled control, intervals not overlapping. So bank similarity explains roughly a third of the
  level-5 between-ETS rows, and the old noise-template control (1.8) understated the false rate about eightfold.
  The first verdict file used all 24 banks (mean of per-bank rates, 9.6/bank-yr → "A"); it is kept as
  `tables/y9nb_verdict_v1_allbanks.json` and superseded.
- The 11 real 2011–13 level-5 between-ETS rows by eye (`tables/y9_level5_between.csv`): 7 wavefield, 2 earthquake wave
  (2011-09-23 14:49, 2012-06-28 13:26), 2 spike (2013-02-03 12:15, 2013-05-13 14:19); **no candidate**. All 5 candidates of
  2010–18 are in 2014–2017, outside the years the matched banks cover.

### Reading
- Between ETS, 39 rows pass every battery′ check. Similarity-matched noise banks predict about 14 (9–22) over 2010–18, so
  the excess is about 25, not 37. By eye, 34 of the 39 are wavefield, earthquake-wave or spike rows, and 5 remain
  candidates (one clear: 2016-06-13 07:41 PDT; one regional-looking: 2015-04-13 06:06 PDT). In 2011–13, where the banks
  were built, the 11 real rows hold no candidate. The honest statement: N = 39 pass every check, matched controls predict
  ~14, 34 are explained by eye, 5 are candidates with the columns above. Extending the banks to 2014–2017 would test the
  candidates' years directly.
- Two things the ladder makes visible: the e′ quake test lets M < 5.5 teleseism windows through (21 of 39), and a
  far-versus-near emergence column separates earthquake waves from local sources. Both are reported as columns, not gates.

## 22. DRAFT — Y9b candidates and Y11 (check f′, stations) — interim (written 2026-10-16 ~5:18 pm PT)
### Correction (pre-registered in Y9b)
- `svi_c010_075` is the SVI LFE-candidate list (0.1 ≤ P < 0.75, all families), not patch LFEs. True patch LFEs (9 SVI
  families ≤ 10 km of the C0/C1 centroids, P ≥ 0.75; `tables/y9b_patch_lfe_2020_2021.csv`): 2020-11 111, **2021-01 0**,
  2021-05 0, 2021-07 0, 2021-08 0, **2021-10 55** (candidates 171 / 10 / 3 / 1 / 3 / 184). So 2021-10 did reach the patches,
  while 2021-01 did not; all-SVI P ≥ 0.75 counts are 617, 394, 38, 49, 0, 344.

### Y9b part 1 — the 5 between-ETS candidates (`tables/y9b_candidates.csv`, `tables/y9b_hf_<k>.csv`, `figures/y9b_cand_0..4.png`)
- No catalogued earthquake within ±5 min in the local lists (ComCat ≤ 300 km M ≥ 2.5, ≤ 1,500 km M ≥ 4, ≤ 2,500 km M ≥ 3,
  global M ≥ 4.5); no local list holds smaller events. No PNSN tremor within 30 km and ±1 h for any; SVI LFE detections within
  30 km and ±1 h: 0–1 candidates each.
- High frequency in the pre-registered window (t0 − 30 … t0 + 150 s): no near station reaches 3× (2–8 Hz medians 0.94–0.99,
  max 1.97), so the pre-registered rule labels all five "VLFE-like". **The figures show the window was too narrow for one:**
  2015-04-13 06:06 has a coherent 2–8 Hz burst on 7 stations (34–147 km) about 100 s before t0, then the long-period train —
  body waves then surface waves of an uncatalogued regional earthquake. Disclosed: the eye check overrides the rule here.
- Classes (rule → by eye):
  - 2015-04-13 06:06 PDT: **regular earthquake** (HF burst at −100 s on all near stations).
  - 2016-06-13 07:41 PDT: **VLFE-like** — no HF at any time in ±10 min, long-period pulse with moveout and decay (Spearman
    −0.57), 15 km from C1 and 1.2 km from SVI family 43, 1–2 LFE detections within 30 km in ±1 h. The one survivor.
  - 2014-05-05 08:30, 2016-06-13 19:58, 2017-01-31 04:22 PT: **unclear** — no HF, but no distinct long-period pulse above the
    surrounding hour; 2014-05-05 and 2016-06-13 19:58 are 56–72 km from the patches.

### Y11 step 2 — f′ re-scored from existing outputs (`tables/y11_rescore.json`, `y11_level5_f_vs_fprime.csv`,
`y11_completeness_by_month.csv`, `y11_completeness_by_year.csv`)
- Known VLFEs: f and f′ both pass 25 of 26; level 5 stays 21 (f is not what holds the other 5 back).
- Pooled 2010–2018 between-ETS level 5, f → f′: real 39 → 40; noise ×11 1.76 → 1.76; A1 2.5 → 2.5; A2 4.0 → 4.0; matched banks
  (accepted, 2011–13, scaled) 14.3 → 15.0. Every control rises by ≤ 1 (largest +0.65).
- C0 completeness at Mw 3.8: in the 38 year-months 2013–2021 where f passes < 0.5 of MF-agreeing injections, f′ adds 0.38 on
  average (79 % of those months ≥ 0.15). By year, C0 f → f′: 2013 0.76 → 0.82, 2015 0.41 → 0.64, 2017 0.14 → 0.25, 2018
  0.25 → 0.45, 2019 0.44 → 0.74, 2020 0.36 → 0.78, **2021 0.15 → 0.64** (2021-10: 0.01 → 0.47; Jul–Dec 0.45–0.76). C1 changes
  less (2021 0.57 → 0.77).
- The one real between-ETS row f′ promotes (2014-10-06 13:52 PDT, `figures/y11_prom_0.png`) is an **earthquake wave**: a
  coherent arrival ~60 s before t0 on every near station and on the far references, pre-window peaks larger than the window,
  in a head2 window (M4.6 at 5,900 km 109 min earlier).
- By Merlin's numeric rule f′ would be adopted (controls ≤ +1, C0 +0.38). Not applied yet: step 1 (which station the
  injections drop) is still running, and step 5 below meets Merlin's "stations before any rule change" condition.

### Y11 step 5 — stations (metadata only; `tables/y11_fdsn_stations_raw.csv`, `tables/y11_gobb_psd.csv`)
- Broadband HH stations (120 s sensors) within 100 km of 48.3 N 123.6 W, 2021–2026, nearest first: **UW.BUCKS 24 km (from
  2023-09)**, CN.VGZ 24, UW.LRIV 28, **CN.SYMB 32 km (since 2016)**, CN.PGC 40, **UW.ODUC 42 (from 2024-05)**, **UW.OSQM 45
  (from 2021-03)**, UW.ODEER 46 (2025-06), **UW.SLDQ 48 (2019, gaps 2021)**, **UW.RNWY 49 (from 2021-04)**, UW.OHOH 55, UW.LOPEZ
  55, UW.TURTL 56, UW.OTR 60, CN.PTRF 66, CN.CLRS 70, CN.GOBB 73, and ~25 more at 68–100 km.
- Already on local disk (shared trace archives): OSQM (2021 301 d, then full years), RNWY (251 d, then full), SLDQ (full),
  ODUC (2024 231 d, then full), OTR. **Not local: SYMB and BUCKS**, the two nearest healthy stations; they would need a waveform
  download (approval needed). The IRIS availability service is retired, so days per year come from local files only.
- CN.GOBB: 20–50 s vertical displacement noise is 1.4–12× CN.PGC (median ~3.6×) on 16 sampled days in 2021 and 2025, with one
  zero-signal day — unusable at long periods as is.

## 23. Y12 pilot 2021 — new stations (written 2026-10-16 ~7:45 pm PT)
Run `prodfy2021` = prodf2021 with the expanded v2 scan (matched filter reused). Tables `tables/y12_qc_2021.csv`,
`tables/y12_fvsfprime_2021.csv/json`, `tables/prodfy2021_*`; Y11 step 1 `tables/y11_step1_loso_dropped.csv`, `y11_step1.json`.

### Y11 step 1 (which station the 2021 injections drop; Mw 3.8)
- Jul–Dec: C0 drops UW.LRIV in 99.4 % of rows (521/524), C1 99.8 %, relay sources drop CN.PGC in 99.0 %. Merlin's acceptance
  holds: the f loss is geometry (LRIV 9 km from C0 carries the scan once VGZ-Z and SNB are gone), not a source bug.

### Station QC 2021 (daily 20–50 s Z level and largest-teleseism amplitude vs CN.PGC)
- All four candidates pass: CN.SYMB (level ×1.43, teleseism ×0.98, 365 d), UW.OSQM (×2.38, ×1.01, 295 d), UW.RNWY (×2.93,
  ×0.89, 247 d), UW.SLDQ (×2.68, ×1.01, 306 d); no flat days. (BUCKS starts 2023-09, ODUC 2024-05.)

### C0/C1 completeness at Mw 3.8 (injections; new-station sources = Green's functions at ×0.5 / ×1 / ×2)
| | pre-Y12 | Y12 ×0.5 | Y12 ×1 | Y12 ×2 |
|---|---|---|---|---|
| C0 2021-10, f | 0.006 | 0.068 | **0.362** | 0.492 |
| C0 2021-10, f′ | 0.475 | 0.480 | 0.463 | 0.497 |
| C0 Jul–Dec, f | 0.004 | 0.044 | 0.334 | 0.506 |
| C0 Jul–Dec, f′ | 0.557 | 0.586 | 0.574 | 0.571 |
| C0 year, f / f′ | 0.15 / 0.64 | 0.23 / 0.67 | 0.48 / 0.66 | 0.61 / 0.66 |
| C1 year, f / f′ | 0.57 / 0.77 | 0.60 / 0.83 | 0.73 / 0.81 | 0.80 / 0.80 |
- Monthly (C0, f): Jan–Jun 0.32–0.53 pre-Y12 → 0.74–0.81 at ×1; Jul–Dec 0.00–0.01 → 0.04–0.46 at ×1 (×2: 0.35–0.64).
- Under f the gain depends strongly on the assumed new-station amplitude (×0.5 gives almost nothing in Jul–Dec); under f′ it
  does not (0.57–0.59), because f′ already lets the LRIV-dominated detections through.

### 2021 on real data (pre-Y12 → Y12)
- T1: 2 → **5** ETS tier events (pass). T2: fail in both (2021 T2 is not testable well: few SVI LFEs near the patches).
  T3: 0/2 → **2/3** (pass). Noise-template rate 0 in both; scrambled rate 1.4 → 2.7 per year (column only). GO stays False
  on T2.
- ETS tier events (Y12): 2021-01 three (2021-02-08 03:34 and 2021-02-10 04:12 PST with LFEs; 2021-01-29 09:59 quiet),
  2021-05 one (06-05 13:27 PDT, quiet, vr 0.07), **2021-10 one** (10-27 05:03 PDT, candidates/tremor only, vr 0.20).
- Between ETS: 5 tier events in both runs (Tier A′ retired; ladder 2019–21 follows at Y11_LADDER_DONE).
- Relay vetting: 2 vetted events (ETS tier + LFEs ±10 min + vr ≥ 0.12), both 2021-01; **no 2021 relay bank** (needs ≥ 3).
  2021-10 gives 1 tier event without LFEs. The true patch LFE count in 2021-10 is 55 (§22).

### Decision (Y12 rule as amended)
- 2021-10 C0 completeness under f rises from 0.006 to 0.36 at the nominal source amplitude (0.07–0.49 across ×0.5–×2), and
  the real-data ETS count rises 2 → 5, independent of that choice. Judged "the stations help": the prodfy re-run for
  2022–2026 was started at ~7:40 pm PT (`scripts/prodfy_queue.sh`, CAT_NW = 6). The pre-Y12 chain was stopped cleanly after
  prodf2023 inj_list (to avoid CPU contention; its 2022/2023 matched-filter outputs are reused). All pre-Y12 outputs are kept.
- f′: by the Y12 rule it is not applied — with the new stations at ×1, C0 f passes 0.55 of MF-agreeing Jul–Dec injections
  (≥ 0.5). Caveat: at ×0.5 f still leaks (0.07), so whether f′ is needed rests on the Green's-function amplitude; a direct
  test comes once a relay bank with empirical new-station sources exists (2022+). f′ stays a column.

## 24. Y11 completed — ladder 2019–2021 and 2021 per episode (written 2026-10-16 ~8:08 pm PT)
Ladder extended to 2019–2021 (`tables/catalogue_ladder_2010_2026.csv`; 2019–2020 = prodf, 2021 = prodfy, the Y12 network).
- Level-5 between ETS (live 0.54 / 0.47 / 0.54 yr): real 0 / 9 / 5 → **14** (95 % 7.7–23.5); noise template ×11: 0 / 0 / 0;
  A1 0.5 / 0 / 1.0 per draw (1.5); A2 2.5 / 0 / 0.5 (3.0); matched-bank estimate (accepted banks, 2011–13 rate 2.58 per
  bank-year, a placeholder for this network) × 1.55 yr = **4.0 (2.5–6.1)**. Completeness at Mw 3.8: 0.65 / 0.55 / 0.42.
- The 14 by eye (`tables/y11_level5_between_2019_2021.csv`, `figures/y11_l5_00..13.png`): 9 wavefield (2020-12-02 ×3,
  2020-12-29, 2021-01-01 ×2, 2021-01-03, 2021-01-10 — storm days and an earthquake coda), 1 earthquake wave (2021-06-27
  08:03), 1 unclear (2020-09-21 16:21: clear near pulse but a 2–8 Hz excess of 3.5 on 2 stations in t0 − 300…+180 s), and
  **3 candidates**: 2020-06-11 13:42 PDT (27 km from C1), 2020-09-11 05:00 PDT (16 km from C0), 2020-10-25 19:25 PDT (16 km
  from C0; a week before the major 2020-11 ETS). For all three: no catalogued earthquake ±5 min, no 2–8 Hz excess ≥ 3 on ≥ 2
  near stations (pulse window and t0 − 300…+180 s), long-period pulse decaying with distance (Spearman −0.60 to −0.77),
  no PNSN tremor or P ≥ 0.75 LFE within 30 km and ±1 h (`tables/y11_candidates_2020.csv`).
- Mechanical check: 26 of the top 50 between-ETS rows are in head2 windows (as in §21); known VLFEs at the top level 21/26.
- 2021 per episode (Merlin step 7; `tables/y11_episodes_2021.csv`; Y12 network, f; C0+C1 completeness with Wilson 95 %;
  implied count = Poisson 95 % of the tier count / completeness at Mw 3.8; relay-set completeness "not comparable"):

  | episode | tier | level 5 | level 4 (f fails) | C0C1 Mw 3.6 | Mw 3.8 | Mw 4.0 | implied count | label |
  |---|---|---|---|---|---|---|---|---|
  | 2021-01 (major) | 3 | 2 | 2 | 0.81 (0.74–0.86) | 0.86 (0.80–0.90) | 0.88 | 0.7–10.2 | |
  | 2021-03 | 0 | 0 | 1 | 0.80 | 0.93 | 1.00 | 0–4.0 | |
  | 2021-05 (major) | 1 | 1 | 0 | 0.86 | 0.82 | 0.89 | 0–6.8 | |
  | 2021-07 | 0 | 0 | 3 | 0.46 | 0.59 | 0.78 | 0–6.2 | |
  | 2021-08 | 0 | 0 | 0 | 0.32 | 0.45 | 0.41 | 0–8.2 | |
  | 2021-10 (major) | 1 | 1 | 1 | 0.28 (0.23–0.32) | 0.42 (0.37–0.48) | 0.58 | 0.1–13.1 | |

  With the Y12 stations no 2021 episode is "bound only" (all C0C1 ≥ 0.2 at Mw 3.8); pre-Y12, 2021-07/08/10 were.

## 25. Y9b matched banks 2011–2017 — verdict corrected (written 2026-10-16 ~8:35 pm PT)
Files `tables/y9nb_summary.csv`, `y9nb_per_year.csv`, `y9nb_verdict.json` (as written by the queue) and
`tables/y9nb_verdict_2010_2018_corrected.json` (no re-run; recomputed from the same tables).
### 53 vs 39 — a bookkeeping error, not new rows
- The verdict's "real_2010_18 = 53" and "scaled_2010_18 = 22.0" read the ladder's pooled row, which since the Y11 extension
  (§24) covers **2010–2021**: 39 (2010–2018, unchanged year by year: 2, 3, 4, 4, 10, 7, 5, 4, 0) + 0 + 9 + 5 (2019–2021) = 53,
  over 7.08 live years instead of 5.54. Neither f′ (a column only) nor any dedup change added rows. The summary code now uses
  the 2010–2018 years explicitly.
- Corrected, 2010–2018: real **39** (95 % 27.7–53.3) vs accepted-bank control 3.11 (2.41–3.96) per bank-year × 5.54 yr =
  **17.2 (13.3–21.9)**.
### Outcome and ratio, both stated
- By the pre-registered rule (A if the scaled control ≥ 20, B if ≤ 5) the corrected control 17.2 gives **intermediate —
  neither A nor B**. The "A" in the queue's verdict file came from the wrong live time.
- Real exceeds the matched control by **2.8×** in 2011–2017 (37 vs 13.0 expected, Poisson p = 4 × 10⁻⁸) and 2.3× over
  2010–2018 (39 vs 17.2). Per year (real vs expected from that year's accepted banks): 2011 3 vs 1.8, 2012 4 vs 0.5, 2013
  4 vs 2.2, 2014 10 vs 2.3, 2015 7 vs 2.8, 2016 5 vs 1.5, 2017 4 vs 1.25. So similarity-matched noise explains roughly a
  third to a half of the level-5 between-ETS rows; an excess remains in every year, largest in 2014.
### The 2014–2017 rows by eye
- All 26 level-5 between-ETS rows of 2014–2017 were already labelled in §21: 15 wavefield, 4 earthquake wave, 2 spike,
  2 candidate, 3 weak candidate. After the Y9b checks (§22) the two candidates are one regular earthquake (2015-04-13 06:06)
  and one VLFE-like event (2016-06-13 07:41); the three weak ones are unclear. By year: 2014 7 wavefield, 2 eq-wave,
  1 unclear; 2015 4 wavefield, 2 spike, 1 earthquake; 2016 1 wavefield, 2 eq-wave, 1 VLFE-like, 1 unclear; 2017 3 wavefield,
  1 unclear.
- The excess over matched control is therefore made of rows that are mostly explained by eye (wavefield days, earthquake waves,
  spikes) — failure modes the noise banks, built from random noise windows, do not reproduce at the same rate. The controls
  were not eye-labelled, so the eye-explained rows cannot be subtracted from the control count.
### 2016-06-13 07:41 against the 2016 control
- 2016: 5 real level-5 between-ETS rows vs 1.5 expected from the matched banks. By eye 4 of the 5 are explained (2 earthquake
  waves, 1 wavefield, 1 unclear), leaving 2016-06-13 07:41 as the single VLFE-like row. One surviving row against an
  expectation of 1.5 control rows (not eye-filtered) is not an excess by itself; whether that event is real rests on its own
  evidence (HF-depleted pulse decaying with distance, 15 km from C1) and on the Y13 repeat test.

## 26. Amendment Y14 — f′ for all years (written 2026-10-16 ~11:00–11:15 pm PT)
- f′ adopted under Y11 (adopt rule met) and Y12 (f still leaks with the new stations: summer 2022 max station share
  0.84–0.90 and f pass 0.00–0.20 in Jul–Oct; summer 2020 likewise). Applied to every year by re-scoring only; f (pre-Y14) kept
  as columns and in `tables/preY14/`. Relay banks are not rebuilt and no extended-bank re-run is made.
- **2022 LFE-catalogue gap:** the March and August 2022 ETS had PNSN tremor within 30 km of the patches (139 and 218 events)
  but almost no P ≥ 0.75 patch LFEs (0 and 1), like 2017-03. T2 (LFE coincidence) is therefore "untestable" for 2022, not
  failed; February 2022 did not reach the patches (2 tremor events, 0 LFEs).
- Per year, f (pre-Y14) → f′ (`tables/y14_summary.csv`; ETS tier events from T1; between-ETS level 5 from the ladder; C0+C1
  completeness at Mw 3.8; vetted relay events = ETS tier + LFEs ±10 min + vr ≥ 0.12):

  | run | ETS events | between L5 | C0C1 Mw 3.8 | GO | vetted relay |
  |---|---|---|---|---|---|
  | prod2010 | 10 → 11 | 2 → 2 | 0.64 → 0.69 | T → T | 6 → 6 |
  | prod2011 | 41 → 41 | 3 → 3 | 0.75 → 0.82 | T → T | 33 → 33 |
  | prod2012 | 57 → 60 | 4 → 4 | 0.78 → 0.84 | T → T | 32 → 35 |
  | prodq2013 | 19 → 20 | 4 → 4 | 0.80 → 0.83 | T → T | 13 → 13 |
  | prod2014 | 1 → 1 | 10 → 11 | 0.78 → 0.81 | F → F | 0 → 0 |
  | prod2015 | 4 → 4 | 7 → 7 | 0.53 → 0.66 | F → F | 3 → 3 |
  | prodq2016 | 3 → 4 | 5 → 5 | 0.37 → 0.40 | F → F | 0 → 0 |
  | prodq2017 | 11 → 11 | 4 → 4 | 0.17 → 0.22 | F → F | 2 → 2 |
  | prodq2018 | 20 → 20 | 0 → 0 | 0.31 → 0.42 | F → F | 18 → 18 |
  | prodrq2018 | 83 → 85 | – | 0.57 → 0.74 | T → T | 70 → 72 |
  | prodf2019 | 36 → 38 | 0 → 0 | 0.63 → 0.81 | T → T | 29 → 31 |
  | prodf2020 | 28 → 28 | 9 → 9 | 0.63 → 0.84 | T → T | 22 → 22 |
  | prodfy2021 | 5 → 5 | 5 → 5 | 0.61 → 0.74 | F → F | 2 → 2 |
  | prodfy2022 (SYMB) | 1 → 2 | (ladder not yet extended) | 0.36 → 0.73 | F → F | 0 → 0 |
  | prodfy2023 | 57 → 63 | (ladder not yet extended) | 0.44 → 0.73 | T → T | 39 → 44 |

- No year's GO verdict changes. f′ mainly raises completeness (largest gains 2019–2023, +0.11 to +0.37) and adds 0–6 ETS
  events per year; between-ETS level 5 changes only in 2014 (10 → 11, the earthquake-wave row of §21's Y11 promotion). Pooled
  2010–2018 level 5: ETS 146 → 153, between 39 → 40; noise ×11 1.76 and 0.48, A1 2.5 / 1.5, A2 4.0 / 7.5 → 8.0 (ETS).
- Relay vetting: years with ≥ 3 vetted events and an existing LOO under f: 2019 (LOO 0.72), 2020 (0.36), 2023 (0.56) — none
  accepted. The LOO tests the matched filter re-finding held-out events and depends little on f, so f′ is unlikely to move
  them past 0.85 (2019 is the closest). Not re-run. Costs if wanted: relay re-vet (bank, relay MF, LOO, battery, verdict) ≈ 7 min
  per year at 6 workers (≈ 20 min at 2), keeping the f-era relay files as copies; an extended-bank re-run (prodfr/prodfyr:
  bank, MF, detect, battery, injection list, injections, summary, R) ≈ 45–60 min per year at 6 workers if a bank is accepted.

## 27. Y13 verification — do the candidates repeat? (written 2026-10-17 ~4:40 am PT)
Files `tables/y13_summary.csv` (Y13 as run), `tables/y13b_fair.csv`, `y13b_fair_by_year.csv`, `y13b_overlap.json`,
`y13b_g1_between_events.csv`, `y13b_top_dedup.csv` (labels), `y13b_timing.json`, figures `y13b_top_0..32.png`, `y13b_timing.png`.
### Y13 as pre-registered
- g1 (4 candidates, singletons; similarity with each other 0.56–0.64, with C0/C1 0.50–0.61): level-5 between ETS 122, 63, 130,
  107 and ETS 101, 57, 95, 90 over 17 years, against a noise-template expectation of 39.7 (between) and 29 (ETS) → outcome A by
  the rule. g2 (3 unclear): 38, 27, 24 between and 28, 16, 36 ETS → C.
### 1. Fair comparison — C0/C1 as single templates, same detector, years and battery (production rows, f′)
- C0: 46 between, **220 ETS**; C1: 17 between, 58 ETS. Noise templates: 39.7 between, 29 ETS per 17 template-years.
- g1 give 1.4–2.8× C0's between-ETS count but only 0.26–0.46× its ETS count. Their ETS/between ratio (0.8–0.9) is like the
  noise templates' (0.73), while the patch templates are strongly ETS-modulated (C0 4.8, C1 3.4). g1 templates are also
  high-rate templates (16,500–26,700 MF ≥ 5.5 detections between ETS vs C0 12,400, noise 13,900), like the over-agreeing noise
  bank in §21. So the g1 excess over noise is not patch-VLFE-like behaviour; the shared real-template background is part of it,
  and the g1 shapes add more.
### 2. Overlap
- 422 g1 level-5 between-ETS rows = 336 distinct events (±30 s): 269 seen by one template, 48 by two, 19 by three, none by all
  four. Only 7 coincide with a C0/C1 level-5 detection (8 with any C0/C1 MF detection).
### 3. Top-10 per g1 template, deduplicated (33 events), by eye with catalogue and 2–8 Hz checks
- **19 earthquake waves, 4 regular earthquakes (2–8 Hz excess ≥ 3 on ≥ 2 stations), 5 wavefield, 1 spike (network-wide
  simultaneous), 1 unclear, 2 re-detections of the source events themselves (2016-06-13 07:41, 2020-09-11 04:59), and 1 event
  11 min after the 2020-10-25 19:25 source** (the only repeat-like row).
- Many earthquake waves arrive at the far references (WVOR, COR) 150–250 s before the near network: surface waves of regional
  earthquakes to the south/southeast, several not in the local catalogue lists.
  - 2013-07-11 02:00:26 PDT (z_best 73.9): such a regional earthquake wave (near 200–650 nm, far references first, no near
    2–8 Hz excess) — not a VLFE.
  - 2023-08-25 01:40 (z_best 110.9): the same pattern.
  - 2021-06-18 02:27 (in_head2): earthquake wave (M4.7 at 6,100 km 28 min earlier; far references first).
  - 2020-09-21 16:21, the Y11 "unclear" row, is re-labelled an earthquake wave (far reference first, 2–8 Hz excess 3.2).
### 4. Timing
- The 336 between-ETS events do not cluster near ETS windows (17.6 % within 7 days vs 27 % for random between-ETS times) or with
  tremor near the patches (6.5 % with PNSN tremor ≤ 50 km within ±1 day vs 10 % random). The busiest days are 2014-01-23 (15)
  and 2020-12-02 (6), the wavefield/storm days of §20 and §24.
### Reading
- The Y13 outcome A comes from the pre-registered comparison with noise templates. With steps 1 and 3 done, it is not supported:
  measured against the patch templates C0/C1 on the same footing, the g1 templates behave like high-rate noise-like templates
  (no ETS modulation), their top detections are regional earthquake waves, earthquakes and wavefield, and their activity is not
  tied to tremor. Apart from re-finding their own source events, the candidates do not repeat in a VLFE-like way. The
  2016-06-13 07:41 and 2020 candidates stay isolated, unconfirmed events (effectively outcome C), and the "real-template
  background" lesson applies to single-template searches too.
