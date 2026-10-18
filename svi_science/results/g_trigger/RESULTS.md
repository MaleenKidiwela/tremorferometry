# G — Do large distant or regional earthquakes trigger deep LFEs under SVI and Washington? RESULTS

Written Sun 2026-10-11, ~6:55 am PT. Rules are in `PLAN.md`, written before any number was computed, with
Revisions R1–R7 at its end. Scripts are `svi_science/scripts/g_trigger_{eq,pgv,build,tests,posthoc,eye,figs}.py`.
Tables are in `tables/` and figures in `figures/`. "Scrambled control" = the same detector on the same day, with each
station site's data shifted to another time of that day. Times are Pacific.

## In plain words

- **No earthquake from 2003 to 2026 measurably triggered deep LFEs, under southern Vancouver Island (SVI) or under
  Washington (WA-NW).**
  - **No instant response** while the surface waves passed. This holds for all earthquakes, for the strong ones
    (≥ 1 kPa peak dynamic stress at the stations), for regional and for teleseismic events, for every network era,
    and both between and during ETS.
  - **No delayed response** in the 1–10 days after.
  - **No dose-response:** stronger shaking did not give more LFEs.
  - **ETS onsets do not follow big earthquakes** more often than chance.
  - No test came near the p < 0.01 pass bar. The smallest p of any headline or per-era stack was 0.05.
- **The biggest shakes did nothing.** Each of these gave 0 LFEs in its surface-wave window at SVI, and at WA where covered (2007+), against
  0.006–0.2 expected:
  - Haida Gwaii 2012 M7.8, 74 kPa measured at PGC;
  - Chignik 2021 M8.2, 41 kPa measured;
  - Ridgecrest 2019, 12 kPa;
  - Sumatra 2004, 11 kPa;
  - Tohoku 2011, 7.7 kPa;
  - Chile 2010, 6.5 kPa;
  - Kamchatka 2025 M8.8, 6 kPa.
- **The trap is real, and the confidence model catches it.**
  - In the first hour after a strong earthquake, the raw matched-filter candidates jump 5–8× at SVI: 3,434 vs 634
    expected in the strong windows. This is body waves, coda and aftershocks.
  - After the P ≥ 0.75 / counted-once rules almost none survive.
  - The scrambled control alone would **not** have caught it at SVI: its raw candidates did not rise (353 vs 412).
    Earthquake energy lines up across stations only in the real data.
  - One leak was found by eye:
    - In the 2018-10-21 off-Vancouver-Island M6.5 / 6.8 / 6.5 sequence, 3 SVI "LFEs" at P 0.75–0.89 sit on
      aftershock arrivals. The P ≥ 0.75 rows that follow (8 more over the next ~8 h) are by the same reasoning suspect.
    - Two WA "LFEs" during the 2021-12-07 Oregon offshore swarm sit at the same second on families 34 km apart.
    - Both were vetoed. Lesson: near big regional sequences, eye-check the record section. P and the scrambled control
      are not enough there.
- **One possible case, not proven: Kermadec Islands M7.3, Sat 2019-06-15 3:55 pm PDT** (about 0.2 kPa at PGC).
  Figure `figures/g_eye_WA_305.png`.
  - **What happened.** Near the end of the 2019 WA ETS there had been 5 days with no PNSN tremor anywhere from
    45.5 to 49.5 N, and no WA-NW LFEs. Then LFEs restarted at 4:58 pm, 64 min after the origin, in the late
    surface-wave coda (main Rayleigh waves at PGC about 40–60 min). PNSN tremor followed at 5:02 pm.
  - **Where.** At 47.9–48.2 N, 123.5 W, the same spot that was active on Jun 6–10. The burst lasted about 18 h:
    62 LFEs on 2 families, 146 tremor epicentres. It then stopped.
  - **Rule check.** The case met the pre-registered single-earthquake rule:
    - 26 LFEs in the window vs a ±24 h local mean of 1.5 and a maximum of 22;
    - all 26 with multi-site support;
    - 0 scrambled control;
    - the eye check shows real multi-station tremor bursts (LRIV, B007, B005, B003), not earthquake energy.
  - **Why it is not a detection.** It looks like a restart of the late-ETS patch, not tremor pulsing with the
    surface waves. A post-hoc chance check gives:
    - such "restart after ≥ 72 h quiet" bursts happen about 2 times a year in WA-NW;
    - the chance that one lands in any of our 309 earthquake windows is 12 %;
    - it is 3 % if only the 72 windows that fell during ETS are counted.
  - **Verdict:** a candidate worth keeping, not a clear trigger. The stress was also low (0.2 kPa). Several
    10–74 kPa earthquakes did nothing.
- **Upper bounds (99 %):**
  - **Between ETS, during the surface waves:** at most 0.043 (SVI) / 0.013 (WA) extra LFEs per strong earthquake.
    That is ≤ 2.0× / ≤ 0.5× the normal count in such a window; the background is tiny, about 0.02 LFEs per window.
  - **Per kPa:** ≤ 76 % (SVI) / ≤ 31 % (WA) extra LFEs per kPa.
  - **During ETS:** ≤ 57 % (SVI) / ≤ 34 % (WA) extra LFEs in the window, for strong earthquakes.
  - **Days 1–10 after, between ETS:** ≤ +61 % (SVI) / ≤ +74 % (WA) of the background rate, or ≤ 20 % / 27 % per kPa.
- **Compared with the literature (from memory, not checked here).**
  - Rubinstein et al. (2007, Nature): the 2002 Denali Mw 7.9 surface waves triggered tremor under Vancouver Island and
    Washington, with tremor pulses in step with the Love waves. Denali is before our catalogue.
  - Rubinstein et al. (2009, JGR, Vancouver Island 1997–2008): among the large teleseisms, Denali was the clear case
    and most others triggered nothing. A threshold of a few kPa was discussed.
  - Our negative result for 2003–2026 is consistent with "triggering here is rare". Reports for Haida Gwaii 2012,
    Chile 2010 or Tohoku 2011 in Cascadia are not reliably recalled, so I do not claim any.
  - The detection limit for small triggered bursts is set by our template families. Triggered tremor that does not
    repeat a known family would be missed. The PNSN tremor column partly covers this. It runs about 2× high in the
    strong windows, but never reaches p < 0.01 (best p 0.05), and its few "starts" are mostly earthquakes logged as
    tremor (§2).

## 1. Earthquakes and peak dynamic stress

- **The list.** ComCat gave 491 earthquakes (356 global M ≥ 7, 135 regional M ≥ 5.5 within 1,500 km, 2003-01 to
  2026-10). They form 441 trials after merging events within 6 h (`tables/eq_list.csv`).
- **Measured stress at PGC (0.02–2 Hz).** 460 events were measured on the CN.PGC broadband files on disk. The rest:
  - 27 with no file;
  - 2 with gaps;
  - 2 flat-topped.
- **The surface-wave relation works.**
  - Estimate vs measured: median log10 ratio +0.06, scatter 0.36 (×2.3). Teleseismic +0.12, regional −0.13.
    Figure `figures/g_sigma_calib.png`.
  - Used for: SVI 31 events; WA 142 (Δ < 1,000 km plus the unmeasured ones).
- **Strong events (≥ 1 kPa):** 167 at SVI, 175 at WA. ≥ 10 kPa: 15 and 16.
- **Trial status (SVI / WA-NW):**
  - between episodes: 297 / 223;
  - during episodes: 98 / 76;
  - in the ± 2-d pad, so excluded: 35 / 19.
  - WA 2007-01 to 2009-08 (only the WECC windows known) is a separate column: 45 between, 4 during.

## 2. Test 1 — instantaneous triggering (W_inst = P arrival to Rayleigh at 2.5 km/s + 1 h)

`tables/inst_stacks.csv`, `inst_trials_{SVI,WA,WA_pre}.csv`; figure `figures/g_superposed.png`.

| stack | n | LFEs obs / random mean | p_R | scr. control obs | p_S | multi-site p | raw candidates obs / random (data; scrambled) | PNSN tremor ≤100 km obs / random |
|---|---|---|---|---|---|---|---|---|
| SVI strong, between | 94 | 4 / 1.96 | 0.17 | 0 | 0.16 | 0.17 | 3,434 / 634; 353 / 412 | 43 / 23.3 (p 0.16) |
| SVI all, between | 297 | 11 / 6.99 | 0.22 | 1 | 0.24 | 0.22 | 7,626 / 2,188; 1,370 / 1,341 | 118 / 78 (p 0.15) |
| SVI strong, W_surf | 94 | 3 / 1.11 | 0.12 | 0 | 0.12 | 0.12 | 2,444 / 379 | 12 / 13.7 |
| SVI strong, 6 h | 94 | 9 / 6.5 | 0.24 | 0 | 0.23 | 0.24 | 5,933 / 2,375 | 106 / 84 |
| SVI regional strong | 36 | 4 / 0.64 | 0.054 | 0 | 0.054 | 0.054 | 2,602 / 200; 122 / 124 | 18 / 6.9 |
| SVI teleseismic strong | 58 | 0 / 1.29 | 1.0 | 0 | 0.99 | 1.0 | 832 / 437; 231 / 290 | 25 / 16.3 |
| SVI ≥ 10 kPa | 13 | 0 / 0.44 | 1.0 | 0 | 1.0 | 1.0 | 229 / 85 | 16 / 2.4 (p 0.058) |
| SVI strong, during ETS (vs ±24 h local) | 36 | 95 / 112 | 0.63 | 0 | 0.63 | 0.63 | 1,569 / 802; 149 / 287 | 131 / 103 |
| SVI all, during ETS | 98 | 246 / 261 | 0.57 | 0 | 0.57 | 0.57 | 3,814 / 2,367 | 601 / 558 |
| WA strong, between | 76 | 1 / 2.10 | 0.36 | 0 | 0.36 | 0.36 | 12,225 / 6,705; 11,642 / 6,721 | 43 / 14.6 (p 0.050) |
| WA all, between | 223 | 3 / 5.89 | 0.45 | 0 | 0.44 | 0.44 | 33,578 / 21,052; 31,134 / 20,908 | 81 / 46.8 (p 0.11) |
| WA regional strong | 32 | 0 / 0.88 | 1.0 | 0 | 1.0 | 1.0 | 5,340 / 1,824 | 7 / 5.2 |
| WA teleseismic strong | 44 | 1 / 1.24 | 0.24 | 0 | 0.24 | 0.24 | 6,885 / 4,871 | 36 / 9.5 (p 0.050) |
| WA ≥ 10 kPa | 9 | 0 / 0.17 | 1.0 | 0 | 1.0 | 1.0 | 2,212 / 452 | 0 / 1.9 |
| WA strong, during ETS | 26 | 5 / 14.9 | 0.82 | 0 | 0.82 | 0.82 | 4,266 / 2,575 | 64 / 81 |
| WA all, during ETS | 76 | 67 / 61.2 | 0.37 | 1 | 0.38 | 0.37 | 11,674 / 8,181 | 428 / 409 |
| WA 2007–09 strong, between | 16 | 1 / 1.52 | 0.28 | 0 | 0.28 | 0.28 | 1,714 / 1,125 | no PNSN |

- **Per network era: nothing passes.**
  - The smallest p is SVI S15f all, 6 vs 1.27 (p 0.085). Next: SVI FWD strong, 3 vs 0.73 (p 0.089), and these 3
    are the vetoed 2018-10-21 aftershock rows. Then WA E2 strong, 1 vs 1.1.
  - ERA, B2f, B3f, CN4f and WA E1 / E3 have 0 LFEs in every strong window.
- **PNSN tremor within 100 km** runs about 2× high in the strong windows, but no row reaches p < 0.01. A post-hoc
  look at the tremor "starts" in windows:
  - 2014-04-23 M6.5 off Vancouver Island: 9 epicentres beginning 5 min after the origin. This is very likely the
    earthquake or its aftershocks counted as tremor, with no LFE.
  - 2014-10-08 East Pacific Rise and 2017-07-05 Lincoln, Montana: small WA tremor starts 38–50 min after, with 0 LFEs.
  - Kermadec 2019 (below).
  - None is supported by LFEs, except Kermadec.
- **Individually named earthquakes (PLAN §3 rule a–d).** Only WA trial 305, Kermadec 2019-06-15. The eye check (e) is
  below and the post-hoc chance check is in §6.

## 3. Test 2 — delayed triggering (days +1…+10, between episodes, data − scrambled control per observable family-day)

`tables/delayed_stacks.csv`, `delayed_trials_*.csv`.

| stack | n | after (obs / random) | p | after − before p | multi-site p | scrambled alone p |
|---|---|---|---|---|---|---|
| SVI strong +1…+10 | 72 | 0.0018 / 0.0024 | 0.61 | 0.58 | 0.61 | 0.95 |
| SVI strong +1…+3 | 74 | 0.0020 / 0.0025 | 0.47 | 0.44 | 0.47 | 1.0 |
| SVI all +1…+10 | 238 | 0.0021 / 0.0026 | 0.72 | 0.77 | 0.72 | 0.85 |
| WA strong +1…+10 | 53 | 0.0030 / 0.0035 | 0.48 | 0.42 | 0.48 | 1.0 |
| WA strong +1…+3 | 60 | 0.0064 / 0.0039 | 0.24 | 0.22 | 0.24 | 1.0 |
| WA all +1…+10 | 172 | 0.0023 / 0.0034 | 0.78 | 0.69 | 0.78 | 0.56 |
| WA 2007–09 strong +1…+10 | 16 | 0.0075 / 0.0092 | 0.56 | 0.78 | 0.56 | 1.0 |

All FAIL. The delayed rate after strong earthquakes is, if anything, below the random-date level.

## 4. Test 3 — dose-response and ETS onsets

`tables/dose.csv`, `ets_onsets.csv`, `ets_onset_hits.csv`; figures `figures/g_dose.png`, `figures/g_ets_onsets.png`.

- **Instantaneous, between episodes (Poisson model N ~ μ(1 + βσ)):**
  - SVI: β̂ 0.19 / kPa, 99 % bound ≤ 0.76, ρ +0.04 (p 0.25);
  - WA: β̂ 0, ≤ 0.31, ρ +0.10 (p 0.080);
  - WA 2007–09: ≤ 0.77.
  - Does not scale.
- **Delayed (+1…+10 d):**
  - SVI ρ +0.08 (p 0.13), slope bound ≤ 20 % of background per kPa;
  - WA ρ +0.05 (p 0.25), ≤ 27 % / kPa.
- **ETS onsets with a strong earthquake in the 10 days before:**

  | region / earthquakes | onsets with one | expected (time shift / uniform) | p (shift / uniform) |
  |---|---|---|---|
  | SVI strong | 8 of 52 | 8.0 / 8.8 | 0.57 / 0.68 |
  | SVI ≥ 10 kPa | 2 of 52 | 1.0 | 0.26 |
  | WA strong | 15 of 66 | 10.8 / 11.9 | 0.105 / 0.185 |
  | WA all | 34 of 66 | 26.8 / 28.5 | 0.049 / 0.097 |

  - 30-day versions are all p ≥ 0.18.
  - Not supported. WA leans the right way at p 0.05–0.19, which is not a result.
  - The 2012-04-14 WA onset came 2–3 d after the Wharton Basin M8.6 + Oregon offshore M6.0 + Gulf of California M7.0
    group. It is listed in `ets_onset_hits.csv`, with no claim.

## 5. Eye check (2–8 Hz record sections; `figures/g_eye_<region>_<trial>.png`)

| figure | earthquake | status | LFEs in W | verdict |
|---|---|---|---|---|
| SVI_179 | Haida Gwaii 2012-10-27 M7.8, 74 kPa | between | 0 | earthquake coda + aftershocks only (simultaneous on all stations) |
| SVI_341 | Chignik 2021-07-28 M8.2, 41 kPa | between | 0 | nothing (one local event at +76 min) |
| SVI_149 | Tohoku 2011-03-10 M9.1, 7.7 kPa | between | 0 | teleseismic body phases only |
| SVI_125 | Maule 2010-02-26 M8.8, 6.5 kPa | between | 0 | nothing LFE-like; slow simultaneous envelope bumps 2.5–4 h later (aftershock / coda energy); 3 PNSN epicentres |
| SVI_157 | M6.4 off west Vancouver Island 2011-09-09 (49.5 N 126.9 W), 29 kPa | between | 0 | aftershocks; 177 raw candidates, all rejected by P = the trap, caught |
| SVI_293 | M6.5 / 6.8 / 6.5 sequence 2018-10-21, 6.3 kPa | between | 3 | **veto**: the P ≥ 0.75 detections sit on aftershock arrivals |
| SVI_144 | Argentina (deep) 2011-01-01 M7.0, 0.10 kPa | between | 4 | tremor burst already running before the earthquake (PNSN and LFEs from −35 min) = coincidence |
| SVI_22 | M5.9 off west Vancouver Island 2004-07-15 (49.7 N 126.9 W), 9.4 kPa | ETS 2004 | 53 (local 32) | ongoing ETS; brief detection gap during the earthquake coda; no change after |
| SVI_16 | New Caledonia 2003-12-27 M7.3 | ETS 2003-12 | 31 (local 13) | ongoing episode, no onset at the surface waves |
| WA_374 | Java 2023-04-14 M7.0, 0.06 kPa | ETS | 20 (local 10) | tremor and LFEs before and after = ongoing ETS |
| WA_349 | Oregon offshore swarm 2021-12-07 M5.8, ~1 kPa | between | 2 | **veto**: 2 detections at the same second on families 34 km apart, among swarm arrivals |
| WA_293 | 2018-10-21 off-Vancouver-Island sequence, 14.6 kPa | between | 0 | aftershocks only; 1,375 raw candidates, 0 kept |
| **WA_305** | **Kermadec 2019-06-15 M7.3, ~0.2 kPa** | ETS (late) | **26 (local 1.5)** | **LFE-like multi-station bursts, starting in the late surface-wave coda and lasting ~18 h; not pulsing with the surface waves** |

## 6. Kermadec 2019-06-15: post-hoc chance check

`tables/posthoc_onsets.csv`; the onset definition was fixed in the script docstring before computing.

- **Definition.** A "start" is the first item after ≥ 72 h with none, followed by ≥ 20 items in 12 h. Items are within
  100 km of the region point, 2010-01-06 to 2026-05-27.
- **WA-NW LFE starts.** These happen 2.1 times a year (35 in 16.4 yr).
  - Expected inside all 309 trial windows (541 h): 0.13; observed 1 (Kermadec), p 0.12.
  - Inside the 72 ETS-status windows: expected 0.030, p 0.030.
  - With a ≥ 10-LFE threshold: p 0.19 / 0.047.
- **PNSN tremor starts within 100 km.** WA: 11.5 a year; 3 in windows vs 0.71 expected (p 0.035). SVI: 3 vs 0.69
  (p 0.033). Most of these are the earthquake-as-tremor cases listed in §2.
- **SVI LFE starts in windows.** Only 2018-10-21, which is the vetoed aftershock case.
- **Reading.** Kermadec is the single best coincidence among ~300 earthquakes, at a 3–12 % chance level. It is not
  a clear trigger. If it were triggered, it would be a clock-advanced restart of a patch close to failure at the
  tail of an ETS, at 0.2 kPa. That is far below the few-kPa level where triggering is usually reported.

## 7. Caveats

- **Families.** LFEs are counted only where our template families exist: 115 SVI families and the WA-NW core set.
  Triggered tremor elsewhere, or triggered tremor that matches no family, is not seen. The PNSN column covers part of
  this, from 2010 only.
- **Low background.** Between-episode instantaneous windows are ~2 h and the background is tiny (~0.02 LFEs per
  window). The relative bounds between ETS are therefore loose (≤ 2× at SVI). The during-ETS bounds (≤ 33–58 %) and
  the delayed bounds (≤ 53–81 %) are the informative ones.
- **Episode lists.** The pre-2010 SVI list is the R23 fallback (Bostock windows plus unvalidated envelope minors). WA
  before 2009-08 knows only the 3 WECC windows, so it is a separate column.
- **Stress values.**
  - σ at WA is the PGC value for Δ ≥ 1,000 km and the estimate otherwise (e.g. Haida Gwaii 59 kPa at WA vs 74
    measured at PGC).
  - Three measured values disagree with the estimate by > 10× and are kept: Haiti 2010, the Illapel M7.0 aftershock,
    and 2018-07-24 Oregon offshore.
  - Two PGC data faults were found and fixed during the work (R5–R7): a sampling-rate jitter, and a midnight
    glitch / day-file overlap. The first clip rule was also wrong.
- **Guards are not perfect.** The P ≥ 0.75 filter leaks a few rows inside large regional aftershock sequences, and
  the scrambled control does not see that (§5). Any future triggered-LFE claim near a regional sequence needs the
  record section.
- **Rate statements.** Between-episode rates are a control-subtracted bracket (RESULTS_BACKWARD rules). Every
  comparison is within one network segment / era; nothing is trended across eras.
