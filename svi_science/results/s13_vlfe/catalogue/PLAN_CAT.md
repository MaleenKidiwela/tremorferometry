# PLAN_CAT — 2003–2026 strong-VLFE catalogue, southern Vancouver Island (stage 1)

Written 2026-10-15 ~1:10 pm PT (Thu), before any battery, e′, E3 or E4 number. Design: `MERLIN_CAT_DESIGN.md` (Merlin,
~1 pm PT). Inputs: `../hg_repro/` (R2 V4 matched filter, R3 per-event battery), `../scripts/scan_columns2.py` (scan
machinery). Disclosure: while sizing the compute I printed the raw V4 detection COUNTS per MF level for T1–T3, the 16
scrambled-control catalogues and 25 noise templates (2011-06-01 – 2013-01-01; no battery, no ETS split; logged in
NOTES.md). Nothing else was looked at.

Words: "scrambled control" = a catalogue made after scrambling the data so no real event can line up (station-day shift
A1: each station's record moved by a different whole number of days; polarity scramble A2: random sign per
station-component of the template). "Noise template" = a 90-s real network window at a random 2011-ETS time used as a
template (keeps the real station-to-station coherence of the long-period wavefield). PT = Pacific time.

## 1. Tiers (final catalogue; the stage-1 tests below decide whether they are reachable)
- **Tier A**: matched-filter detection ≥ the template-year threshold (§2) AND battery a, b, c, e′, f (§3) AND
  two-template agreement (§4) AND the scrambled-control rate of the same tier ≤ 0.3 per year between ETS.
- **Tier B**: detection ≥ threshold AND battery a, b, c, e′, f AND scrambled-control rate ≤ 1 per year between ETS.
- Everything else stays in the table with all columns (rank, don't discard), labelled "candidate".
- If e′ fails its pass rule (§6), (e) replaces e′ in both tiers.

## 2. Threshold rule (never from Table S1 counts)
- Per template and per calendar year: the summed-CC level at which the 25 noise templates, run under that template's
  channel configuration over the same year (peaks within ±90 s of a noise template's own window excluded), fire on
  average **0.5 per noise template per day**. Columns: 0.25 and 1 per day.
- In practice the level is found in MAD units: each statistic is expressed as (S − median)/MAD of its own year; the
  level z* where the mean noise-template rate equals 0.5/day is applied to the real template's own (S − median)/MAD.
  The scrambled controls of that template are cut at the same absolute summed CC as the real template.
- E2 (2011–12) uses the existing V4 run and reports fixed MF levels 5.5, 6, 7, 8, 10 MAD (real template's median/MAD
  over the 914-day run, `../hg_repro/tables/r2_thresholds.csv`), so that the battery can be judged at every level. Its
  controls and noise templates are cut at the same absolute summed CC as the real template.

## 3. Battery (per detection; each test a column; pass rules fixed here)
Machinery: `../hg_repro/scripts/hg_r3.py` (ev_row: scan best point within t0 ± 30 s, `scan_columns2.cols()`), scan years
2011 and 2012 only (stage 1). t0 = origin-aligned detection index.
- **(a) far/near: ratio_far** = median over live far LHZ references {IU.COR, US.NEW, US.HAWA, UW.OMAK, TA.I02D} of
  A_ref / A_near. **Pass a: n_far ≥ 2 AND ratio_far ≤ 0.551** (R1 threshold). n_far < 2 → untestable (not a pass).
  Columns: ratio6 (R3 rule incl. GNW), ratio_GNW (location-confounded, column only), ref_ratio.
- **(b) not a travelling plane wave** (R3 rule): travelling = pw2_n ≥ 8 AND pw2_vr ≥ 0.5 AND ≥ 2 kept stations ≥ 250 km
  from the best node. Pass b: not travelling.
- **(c) located**: z_best ≥ 5 AND best node not on the grid edge AND (inside the scan footprint polygon OR ≤ 20 km from
  the nearest LFE family, `data/catalogue/families_v1_3.csv`, R3 regions).
- **(d) LFE / tremor status** (column, never a rule): SVI P ≥ 0.75 + WA-NW selected_v2 within ±10 min → "with LFEs";
  else SVI 0.1 ≤ P < 0.75 or PNSN tremor ≤ 100 km → "candidates/tremor only"; else "quiet".
- **(e) earthquakes (R3 rule)**: no head2, bw5, reg or reg3 mask second in [t0 − 60, t0 + 120].
- **(e′) refined earthquake test (pre-registered alternative to e)**: FAIL if any bw5, reg or reg3 second in
  [t0 − 60, t0 + 120]; OR a head2 window there from an event with M ≥ 5.5; OR any head2 window there AND ≥ 2 of the five
  far LHZ references excited (ex_ref ≥ 3, the R1 "3 noise units" rule). Otherwise pass. (A head2 window of a distant
  M 4.5–5.4 teleseism with quiet far stations no longer vetoes.)
- **(f) not single-station** (R3 rule): z_loso ≥ 3 AND z_loso ≥ 0.5 × z_best.
- **Battery = a ∧ b ∧ c ∧ e ∧ f**; **battery′ = a ∧ b ∧ c ∧ e′ ∧ f**.
- **S5 column** (scan headline filter as validated, PLAN_SCAN S5/S8): n_sta ≥ 3, n_Z_agree ≥ 2, vr_thrust/vr_free ≥ 0.8,
  vr_thrust > vr_strikeslip and > vr_steepdip, isolation < 2, share_excess ≤ 0.25. Reported as "± S5" (battery with and
  without S5); not a tier rule unless stated.
- Compute shortcut (does not change any all-test count): the expensive columns (a, b, f, S5, far excitation) are computed
  only for detections that pass (c) and pass (e) or the time-mask part of (e′); others are battery failures by (c) or
  (e/e′) and keep NaN in a/b/f. Cheap columns (z_best, location, c, e, e′ time part, d) are computed for every detection.

## 4. Two-template-agreement tier
A detection of template K at level L agrees if another template K′ ≠ K (any of T1–T5 in E2; any template of the bank
in E3/E4) has a detection at level ≥ L within ±30 s. Scrambled controls: the same draw (A1_r / A2_r) of K′. Noise
templates: noise j and noise j+1 (mod n) under K's configuration.

## 5. Live-time rule
Rates per year use live time: seconds of the window in which (i) the MF template has ≥ 50 % of its channels live
(n_live), and (ii) the e′ time masks (bw5, reg, reg3, head2 of M ≥ 5.5) are off. Plain calendar rates are columns.
(For battery (e) counts the full head2|bw5|reg|reg3 mask is used for the live-time column.)

## 6. Pass rules (fixed now)
Periods: ETS = H&G windows 2011-07-23 – 09-06 and 2012-08-30 – 10-11 (UTC dates); between ETS = 2011-09-08 – 2012-08-28
(R3 window, 0.972 yr). "Control rate" = mean over the 16 scrambled-control catalogues (8 A1 + 8 A2) of battery passes
between ETS per year; "max(control)" = the largest single-catalogue count between ETS.
- **E2** (existing V4 T1/T2/T3 detections, 2011–12): accept a tier (MF level × (e or e′) × (± S5) × (single or
  two-template)) if **ETS real ≥ 8 AND between-ETS real ≤ max(control) + 2 AND control rate ≤ 1/yr**. Noise-template
  passes per template-year are a column. If several tiers pass, the lowest MF level wins (completeness), then e′ over e.
- **e′ validation (Merlin step 3)** on the R3 sets with ratio_far as (a): accept e′ if **ETS real-looking goes from 8 to
  ≤ 14, scrambled control ≤ 2/650, between-ETS stays 0/130.** Columns: selection-matched noise (1,300), ETS scrambled
  (575), ETS noise (1,150).
- **E3/E4** (empirical templates, 2011–12): accept the all-station bank if **≥ 7/8 vetted events are re-found by a
  template not containing them** (leave-one-out, detection ≥ threshold within ±30 s); accept the 4-station set
  {PGC, VGZ, NLLB, SNB} if **≥ 5/8 re-found AND its scrambled-control battery rate ≤ 1/yr**. The 4-station result decides
  whether 2003–2009 can be catalogued at ~Mw 3.6–3.8 (pass) or only Mw ≥ 4.0 (fail). {H&G 9} and {4 + GNW} are reported
  against the same rules as columns.
- **Pilot 2014** (later): GO if ≥ 5/8 of H&G 2014 T6–T13 template events re-found, between-ETS ≤ control + 1 with
  controls ≤ 1/yr, and ≥ 80 % completeness at Mw 3.8.
- **Pilot 2008** (later; only if the E3 4-station rule passed): same rules as 2014 with the 2011–12 vetted events as
  truth via the 4-station bank; else 2003–09 = Mw ≥ 4.0 tier + Poisson bound and the uniform catalogue starts 2010.

## 7. Completeness rule (cross-template injection; used from the 2014 pilot on, defined now)
Inject each vetted event's 90-s multi-station waveform (all live channels, scaled by 10^(1.5 ΔMw)) at 200 random times
per Mw step (3.4–4.2, step 0.2) per year, into the continuous data; detect with the bank built WITHOUT that event
(leave-one-out); recovered = detection ≥ the template-year threshold within ±30 s AND battery′ pass. Completeness(Mw) =
recovered / injected; Mw_c = the lowest Mw with ≥ 80 %. Injection times avoid ±10 min of catalogue detections.

## 8. E3/E4 design (Merlin step 4)
- Events: the 8 vetted ETS events + 6 near-misses (R3 ETS candidates passing a(ratio_far), b, c, f but failing (e) by a
  head2 window only): `../hg_repro/tables/r3_ets_events.csv`.
- Stations: all = H&G 9 (NLLB, VGZ, LRIV, DOSE, GNW, SP2, A04D, B05D, D03D) ∪ scan 7 (PGC, VGZ, SNB, NLLB, PFB, YOUB, LRIV)
  = 13 distinct (Merlin's "16" counts the 3 shared stations twice). Sets: {all 13}, {H&G 9}, {PGC, VGZ, NLLB, SNB},
  {PGC, VGZ, NLLB, SNB, GNW}.
- Recipe = V4: 1-sps displacement (v3 day files) → 0.02–0.05 Hz (hg_cont) → velocity (central difference); window 90 s
  starting at origin − 9 s; normalised CC per channel, plain sum over live channels, 90-s declustering.
- Templates: single-event windows clustered by best-node location (COMPLETE linkage, 15 km — amended ~1:15 pm PT before
  any E3 number: the scan grid is 10 km, so single linkage would chain adjacent nodes into one cluster); per cluster a stack of
  unit-normalised per-channel windows after ±5-s alignment to the cluster's first stack; plus H&G T2 and T4 re-cut on
  all stations. Leave-one-out: for each vetted event, its cluster's stack rebuilt without it.
- Controls per template and set: 8 A1 (station-day shift) + 8 A2 (polarity), cut at the real template's absolute
  threshold; 25 noise templates per set (the first 25 R1 noise times) → threshold (§2), per year 2011 and 2012.
- Battery on every detection (cached per unique time ±15 s) as §3; rates per set as §6.

## 9. Step 5 (metadata only)
FDSN station metadata (no waveforms) for LHZ 2003–2026 at 250–800 km from the SVI scan centre (48.37 N, 123.61 W),
networks IU, US, UW, TA (incl. 2006–08), BK, CC, CN, II; EarthScope availability extents (metadata). Table per year:
candidate references with distance and azimuth; count with ≥ 180 days; size estimate from the on-disk 2011 LHZ day
files (MB per station-year). No download until the user approves.

Compute: ≤ 8 CPU workers, OMP_NUM_THREADS=1, < 100 GB, GPU allowed. Stages > 20 min run under nohup with DONE markers in
`catalogue/`. Never delete files. Raw figures first, points coloured, no error bars.

## Amendment P (2026-10-15 ~1:45 pm PT, before any 2014 or injection number) — pilot 2014 (A) and early-era injection (B)
Coordinator (user-approved extension): run Merlin step 6 (pilot 2014) and a 4-station cross-template injection on 2011–12.
Facts used (no result numbers): E3 all-13 bank failed its rule (6/8), H&G 9 met it (7/8); E2 validated the V4 H&G-9 tier
5.5 MAD + e′ + two-template agreement. 2014 raw data are on disk for all 13 stations (TA.B05D only 37 days, CN.YOUB to
2014-09-03); far LHZ 2014 for IU.COR, US.NEW, US.HAWA from the approved download (OMAK, I02D not fetched, so n_far ≤ 3).

### P1. Common machinery
- **Source waveforms (cross-template)**: the two multi-event clusters of E3, C0 (5 events) and C1 (4 events). Per channel
  (13 MF stations, velocity, 0.02–0.05 Hz; 7 scan stations, displacement, scan band-pass), each member's window
  [origin − 60, origin + 121) (origin = R3 index + E3 alignment shift) is projected on the unit stack of the members;
  signal = mean projection × unit stack. Reference Mw = mean member Mw_tpl (R3 template scaling). Injected amplitude =
  signal × 10^(1.5 (Mw − Mw_ref)). Far references: the scan's thrust synthetic for the cluster's mean node scaled to the
  same M0 (a real VLFE is not silent there). Nothing is injected at other references.
- **Detection of an injection** (cluster c injected): bank templates that contain no member of c (all others). The MF
  statistic is recomputed locally on the injected data (exact: the CC is local), normalised by each template's
  uninjected median/MAD; detected = a template ≥ the tier threshold with its detection index within ±30 s of the
  injected origin.
- **Battery on injected data**: the scan machinery (zmax search over all grid nodes in t0 ± 30 s, cols, plane wave) run
  on the injected scan data (local recomputation of the scan statistic; checked against the stored zmax at uninjected
  times before use); e′ masks at the injection time.
- **Injection times**: uniform random in the period's live time (e′ time masks off, all MF channels of the set live at
  the injected window), ≥ 10 min apart, ≥ 10 min from any real detection of the tier, one injection evaluated at a time
  (others absent).

### P2. Pilot 2014 (A)
- Frame 2014-01-01 – 2015-01-01; scan year 2014 (same recipe; scan_common gains a 2014 branch, tables prefixed scan2014_).
- ETS = SVI 2014-11 episode ±2 d (2014-11-12 – 2014-12-16, `b_ets_record/episodes_b.csv`); between ETS = 2014 outside every
  SVI and NW episode ±2 d (live-year length reported); other episodes = a third period (column).
- **Primary tier = the E2-validated one** (all-13 bank failed E3): V4 recipe on the H&G 9 stations, threshold 5.5 MAD of
  each template's own 2014 statistic, two-template agreement, battery′ (a ratio_far n_far ≥ 2, b, c, e′, f).
  Column = the E3 H&G-9 configuration: threshold where 25 noise templates (random 2014 windows inside the ETS window with all
  9 stations live, ≥ 180 s from T6–T13 origins) fire 0.5 per template-day, single template + battery′.
- **Bank** = H&G T1–T5 (2011 windows, V4, as E2) + E3 cluster stacks C0–C5 (H&G-9 channels) + H&G T6–T13 (2014, window
  origin − 9 s). Controls: 8 day-shift + 8 polarity draws per template (same draws across templates).
- **Truth**: H&G T6–T13 origins (Table 1). Tk is re-found if bank \ {Tk} has a tier detection (≥ 5.5 MAD within ±30 s,
  agreement by a second template of bank \ {Tk}, battery′ pass). MF-only re-find is a column.
- **Completeness**: C0 and C1 injected at Mw 3.4, 3.6, 3.8, 4.0, 4.2, 100 times each per Mw per source (200 per Mw), into
  2014 between-ETS live time; recovered = primary tier on bank \ {templates containing the source's members}.
- **GO** if ≥ 5/8 truth re-found AND between-ETS real ≤ mean scrambled-control count + 1 (same period) AND scrambled
  rate ≤ 1/yr AND completeness ≥ 80 % at Mw 3.8. Otherwise NO-GO (report which item failed).

### P3. Early-era completeness (B) — decides the Mw cut for 2003–2009
- Data 2011–12 (E3 frame), MF stations {PGC, VGZ, NLLB, SNB}, bank = E3 bank minus templates containing the source's members
  (C0 injected → C1–C5, C1/C2 LOO stacks, HG_T2, HG_T4; C1 injected → C0, C2–C5, C0/C2 LOO stacks, HG_T2, HG_T4).
- Tier (the early-era tier): E3 cn4 threshold at 0.5 noise per template-day per year (5.84–6.08 MAD), single template,
  battery′ with an **early-era scan**: the 2011–12 scan without UW.LRIV (installed 2009-08) and with far references
  restricted to IU.COR, US.NEW, US.HAWA (UW.OMAK from 2008-11, TA.I02D from 2010). Columns: full 2011–12 battery,
  two-template agreement, MF only, thresholds 0.25 and 1 per day.
- Mw 3.6, 3.8, 4.0, 4.2 (3.4 column); seasons summer (Jun–Sep) and winter (Nov–Feb); 50 injections per source per Mw per
  season (100 per Mw per season, 200 per Mw).
- **Rule: Mw_c(early) = the lowest tested Mw with recovery ≥ 80 % in BOTH seasons; if none ≤ 4.2, the early era is
  catalogued only above 4.2 (bound) or not at all (stated).** Plus the early-era scrambled-control rate (E3 cn4: 1.16/yr
  at 0.5/day) must be ≤ 1/yr at the chosen tier; if not, the threshold column at 0.25/day is used for both (rule fixed now).
- Caveats fixed now: VGZ/SNB/NLLB sensors before 2007–09 differ (CMG40T/CMG3ESP); 2011–12 noise stands in for 2003–09.

## Amendment Q (2026-10-15 ~3:21 pm PT, before any number of its steps; the one-line lfe_lists edit was made ~3:18 pm PT just before it) — Merlin review 2 (`MERLIN_CAT_REVIEW2.md`)
### Q0. Bookkeeping fixes (before any re-summary)
- `hg_r3.lfe_lists()` clipped the LFE/tremor lists at 2013-01-01 → every LFE/tremor label after 2012 was void (all
  2014 "quiet" labels). Fixed to 2027-01-01 (pre-fix copy `data/hg_r3_before_lfe_fix.py`); 2014 labels recomputed.
- Far references: all five downloaded LHZ stations (IU.COR, US.NEW, US.HAWA, US.MSO, US.WVOR) are prepped for each
  pilot year and enter test (a) as far references (ratio_far = median over the live far set; n_far ≥ 2). For 2011–12
  the set stays COR/NEW/HAWA/OMAK/I02D (unchanged). MSO/WVOR get thrust reference templates in the pilot-year
  scan directories (the 2011 `templates_ref.npz` is not modified).
- **Validated bank = H&G T1–T5 + E3 C0–C5 (+ HG_T2/HG_T4 re-cuts where the recipe uses origin − 9 s).** Templates not
  validated by a leave-one-out test (e.g. H&G T6–T13) are columns only: they neither detect nor count for agreement.
- **GO rule rewritten (applies to every pilot from now on, and to the 2014 re-summary): the between-ETS count is a
  reported result, never a gate.** GO = truth items + scrambled-control rate ≤ 1/yr + completeness items.
- RESULTS §8 statement "PGC very noisy in Dec 2014 (99th pct 700–1,000 nm)" is checked against the daily 10-min RMS on
  the event days and kept or retracted as the numbers say.
- Expectation stated by Merlin (not a rule): 2014 re-summary ETS 0 / between 2 / other 1.
### Q1. Station-era table and 2011–12 leave-one-out re-test
- Table: per candidate near station (≤ ~150 km; H&G 9, scan 7, LZB, C8 SHVB/KHVB/MGCB/TWBB/SHDB, UW SQM/HDW/OTR/WISH,
  XU BS11/W040, US.NLWA, CC.OBSR, TA.D04D), per year 2003–2026: days on disk, BROADBAND days (BH/HH with a
  broadband sensor in the epoch table; EHZ/SH* = short period = not usable at 20–50 s), ≥ 180 broadband days flag,
  and data on each of the 14 template-event days (2011–12).
- Sets (Merlin): era2007 {PGC, VGZ, SNB, NLLB, GNW, SQM, HDW}, era2016 {PGC, VGZ, LRIV, SNB, PFB, NLLB, DOSE, GNW, SP2,
  HDW, SQM, OTR}, max2011 = every broadband station of the table with data on ≥ 7 of the 8 vetted-event days; plus hg9.
  In 2011–12 only broadband channels enter (a station without broadband data on an event day simply has no channel
  there); stations of a set that are short-period in 2011–12 are listed as absent (rule fixed now, applies to SQM/OTR
  if the table shows EHZ only).
- Statistic: plain sum (E3) vs **weighted sum S_w = Σ_c w_c cc_c, w_c = template RMS / RMS of the 10 min before the
  window (per channel; for stacks the mean over members of each member's ratio), normalised so Σ w_c = n_c**; z by the
  statistic's own median/MAD; thresholds from 25 noise templates per configuration (0.5/template-day per year), same
  weights for noise templates computed the same way at their own times; controls 8 day-shift + 8 polarity.
- Accept: weighted ≥ plain (re-found count, per set); hg9 ≥ 7/8; era2016 ≥ 7/8 with scrambled battery′ rate ≤ 1/yr;
  era2007 ≥ 5/8 with scrambled ≤ 1/yr (→ 2007–09 joins; 2008-05 pilot later). Statistic for later stages = weighted
  if it is ≥ plain on hg9 and era2016, else plain.
### Q2. Triage of the three 2014 between-ETS events (after Q0 re-summary; the events are whatever the Q0 tier leaves)
- Real-looking iff vr_thrust ≥ 0.12 AND no comparable far wavetrain (max |Z| over the 5 far refs in [t0, t0 + 180 s]
  < 0.5 × the near-station max in the same window) AND (ratio_GNW ≤ 0.551 if the best node is north of 48.1 N);
  else "likely noise". Figures: near scan stations, H&G network, GNW and all 5 far LHZ refs. vr_thrust ≥ 0.12 becomes a
  column for the next pilots.
### Q3/Q4. Pilots 2013-09 and 2018-06 (rules from Merlin review 2)
- Frames: calendar 2013 and 2018; scan years 2013 and 2018 built like 2014 (5 far refs). ETS = SVI 2013-09 / 2018-06
  episode ±2 d (`episodes_b.csv`); between = outside every SVI/NW episode ±2 d.
- Bank: validated templates on the pilot-year station set (2013: hg9 if E3 machinery; statistic from Q1; 2018: era2016).
- Truth T1: ETS tier count ≥ 5 AND Poisson p < 0.01 against the scrambled-control rate in the same live time.
  T2: ≥ 70 % of ETS tier events with LFEs within ±10 min AND LFE rate ≥ 3× random same-day times.
  T3: a near-station pulse (largest at the station nearest the best node, absent at the far refs) on the top event (figure).
- Controls: 8 + 8 scrambled and noise-template pairs, battery′ tier rate ≤ 1/yr between ETS (gate).
- Completeness (C0/C1 injections, ETS and between): ≥ 40 % at Mw 3.4, ≥ 70 % at 3.6, ≥ 80 % at 3.8 (gate).
- Between-ETS count reported only. A fail of T1/T2 with completeness passing = finding ("strong VLFEs did not repeat").

## Amendment Q1′ (2026-10-15 ~3:36 pm PT, before any Q1 number) — Merlin ADDENDUM (density-aware), with a sensor check
- Q0 re-summary of 2014 was run before this amendment (numbers in RESULTS §9); Q2 triage of its 3 events was run with the
  Q2 rule as written (RESULTS §9).
- **Sensor rule (Merlin): only sensors flat to ≥ 30 s.** Checked in the StationXML (velocity response at 30 s / 1 Hz):
  the CAFE stations XU.BS11, W020, W040, PL11, W060 are CMG3T (0.998) in 2006–2010 but **CMG40T-1 1-s sensors (0.0016)
  from 2011-06/07 to 2012-05** → excluded in 2011–12; the 2011 templates cannot be cut on them. UW.SQM is broadband only
  until 2011-07-02, UW.HDW broadband only from 2023, UW.OTR from 2020-10 (EHZ before) → absent in 2011–12.
  UW.ELW (CMG40T-60 s, 0.97), UW.LEBA (STS-2), UW.RATT, CC.OBSR (0.70), US.NLWA (0.98), UW.WISH (T120) pass.
- Sets on 2011–12 (rule-based, so they can be recomputed for any year):
  - dense2011 = every broadband ≤ 160 km with data on ≥ 7 of the 8 vetted-event days: PGC, VGZ, SNB, NLLB, YOUB, LZB,
    LRIV, DOSE, GNW, SP2, A04D, B05D, D03D, D04D, WISH, NLWA, OBSR, ELW, LEBA, RATT (20 stations; PFB has 3/8);
  - era2007 = stations broadband both in 2007 and on the 2011–12 event days: PGC, VGZ, SNB, NLLB, GNW (CAFE/SQM/HDW are
    not broadband in 2011–12, so this is the 2011-template subset of Merlin's era2007 list; a fail here does not rule out
    2007–2010 with era-local templates);
  - era2010 = era2007 + LRIV + PFB (W060 excluded, 1-s in 2011–12);
  - era2016 = broadband with ≥ 180 days in 2016 and data on ≥ 7/8 event days: PGC, VGZ, SNB, NLLB, LRIV, DOSE, GNW, SP2,
    WISH, NLWA (+ PFB, which has 3/8 event days but runs 2016; kept as listed by Merlin);
  - hg9 and cn4 = controls.
- Plain vs weighted (Q1 definition). Accept (ADDENDUM): weighted ≥ plain; hg9 ≥ 7/8; dense2011 = 8/8 with scrambled
  battery′ rate ≤ 1/yr; era2007 ≥ 5/8 with scrambled ≤ 1/yr (→ 2007–2010 enters).
- Pilot 4 = 2008-05 (CAFE era) instead of 2018-06 (backup 2010-08); pilot 2013-09 on the densest 2013 set. The scan/battery
  network is the fixed scan station list (7 + LZB) with its grid templates; extending it to all live broadbands needs new
  grid templates and is a separate step (not done here). For 2008 the scan gain table lacks the VGZ BH 2007–09, SNB BH
  2006–08 and YOUB 2003–09 epochs; the Q3/Q4 rule (g = 1, Z only, for epochs without a gain row) applies.
- Q3/Q4 implementation details (fixed ~3:45 pm PT, before any pilot number):
  - Pilot 2013 MF set = dense2011 stations with ≥ 180 broadband days in 2013: PGC, VGZ, SNB, NLLB, PFB, LZB, LRIV, DOSE, GNW,
    SP2, A04D, D03D, WISH, NLWA (14). Pilot 2008 MF set = era2007 = PGC, VGZ, SNB, NLLB, GNW (the stations that carry the
    2011 templates and exist in 2008; the CAFE broadbands of 2008 have no 2011 templates and the scan is not extended).
  - 2008-05 episode dates from the Bostock list used elsewhere in the project: 2008-05-10 – 06-02 (ETS window ±2 d);
    `episodes_b.csv` has no SVI episodes before 2010, so 2008 "between" = outside NW episodes and the ETS window.
  - Scan gain table: epochs without a gain row get g = 1, Z only (rows appended for PGC/VGZ/SNB/NLLB HH after 2017 and
    for VGZ BH 2007–09, SNB BH 2006–08, YOUB HH 2003–09).
  - Statistic for both pilots = Q1 verdict (`tables/q1_verdict.json`); injections C0/C1 at Mw 3.4–4.0, 50 per source × Mw ×
    period (ETS, between).

## Amendment R (2026-10-15 ~5:55 pm PT, before any 2018 number and before the 2013 re-summary) — Merlin review 3
Recorded first: **under Amendment Q the 2013-09 pilot was NO-GO** (T1 pass: 16 ETS vs scrambled 0, p ≈ 0; T2 fail: 94 % with
LFEs but LFE rate 1.36× same-day random < 3; T3 fail: largest amplitude at PGC, not among the 2 nearest; controls pass
0.54/yr; completeness 0.42/0.69/0.78 at Mw 3.4/3.6/3.8 vs 0.4/0.7/0.8). The 2008-05 pilot under Q was non-informative
(0 detections anywhere; scan had only PGC N/Z active: VGZ and SNB excluded in all 12 months by the 2011-referenced noise
rule, PFB/YOUB without 2008 data; MF on 5 stations shape-limited, median best z 4.3–5.2). These stay on record.
Merlin designed T2′/T3′ after seeing the 2013 numbers (disclosed in his review); they are applied unchanged to 2018.
- **T1′**: ETS tier events ≥ 5 AND Poisson p < 0.01 with λ = (between-ETS scrambled-control rate per year) × (ETS live
  years). The ETS-window scrambled mean is a column.
- **T2′ = T2a ∧ T2b.** T2a: ≥ 70 % of ETS tier events have ≥ 1 SVI P ≥ 0.75 or WA-NW LFE within ±10 min. T2b (patch
  coincidence): S = Σ over ETS tier events of the SVI P ≥ 0.75 LFEs from families ≤ 20 km of the event's best node within
  ±10 min; null = 1,000 sets in which every event time is replaced by a uniform random time on the same UTC day (same node);
  pass if S > 99th percentile of the null (p = fraction of null ≥ S).
- **T3′** on the top 3 ETS tier events by MF z: amplitude = max |Z| (0.02–0.05 Hz displacement) in [t0, t0 + 180 s] at every
  live scan and regional station (near network + GNW, DOSE, SP2, D03D, B05D, WISH, NLWA) and the 5 far LHZ; R = distance to
  the best node. Near = stations with R ≤ 40 km (if < 2, the 2 nearest), mid = 40 < R ≤ 150 km. An event passes if
  median(near)/median(mid) ≥ 1.5 AND median(far) ≤ 0.5 × max(near) AND Spearman(R, amplitude) ≤ −0.4 over near + mid + far.
  T3′ passes if ≥ 2 of 3 events pass.
- **Controls** (gate): between-ETS scrambled tier rate ≤ 1/yr and noise-pair tier rate ≤ 1/yr.
- **Completeness** (reported, not a gate): C0/C1 injections, 800 per Mw per period (ETS, between; 400 per source), Mw 3.4,
  3.6, 3.8, 4.0; recovery with Wilson 95 % intervals; split MF-only (two templates ≥ 5.5) vs battery′-given-MF; Mw_c = lowest
  Mw with recovery ≥ 0.8 (else "> 4.0").
- **Columns for between-ETS events**: storm day = the UTC day's median 10-min RMS of PGC Z in 5–10 s (from the 1-sps day
  files) ≥ 85th percentile of that year's days; cluster day = ≥ 2 between-ETS tier events on that UTC day; Q2 class
  (vr_thrust ≥ 0.12, far/near < 0.5, ratio_GNW ≤ 0.551 north of 48.1 N). **Tier A′** = tier ∧ Q2 "real-looking" ∧ not
  (storm day ∧ cluster day); reported next to the raw tier count.
- **Decision (Merlin (4))**: GO for production 2010–2018 if, under R, 2013 reproduces T2b p < 0.01 and T3′ ≥ 2/3, and the
  out-of-sample 2018-06 pilot passes T1′, T2′, T3′ and controls. Change of mind (reframe as "2011–13 VLFEs + bounds") if
  T2b p > 0.05, T3′ < 2/3, or 2018-06 gives 0 ETS tier events with completeness ≥ 0.8 at Mw 3.8 while the patches are
  active (P ≥ 0.75 LFEs within 10 km of C0/C1 in the episode).
- **2018-06 pilot**: stations = previous 2018 set + UW.WISH (PGC, VGZ, LRIV, SNB, PFB, NLLB, DOSE, GNW, SP2, WISH); 2011
  bank (T1–T5, C0–C5); weighted statistic (Q1: weighted ≥ plain on every set); injections n = 800/Mw/period.
- **Relay (era-matched) banks**: accepted only with leave-one-out ≥ 85 % of their vetted events and scrambled ≤ 1/yr.
- **Scan extension v2**: scan station list per year = every live broadband (sensor flat ≥ 30 s) ≤ 160 km with a bp year file;
  gains from `sc_gain.csv` where an epoch row exists, else g = 1 on Z only; same grid, recipe and noise rule. Accept v2 if,
  on the 2013 injections, battery′ acceptance given MF detection at Mw 3.8 is ≥ 0.84 (ETS) and ≥ 0.72 (between) AND the
  2013 between-ETS scrambled tier rate stays ≤ 1/yr and ≤ the v1 rate + 0.5/yr.
- Amendment R implementation detail (~5:58 pm PT, before any n = 800 number): with 400 injections per source × Mw × period
  the ETS window cannot hold all of them ≥ 10 min apart, so the ≥ 10-min spacing applies within each (source, Mw) set;
  sets may reuse times (each injection is evaluated alone, the others absent). The 2013 re-summary was first run on the
  existing n = 200 injections (reported as such), then re-run with n = 800.
- Scan v2 implementation (~6:00 pm PT, before any v2 number): the extra stations are the regional broadbands with a bp
  year file and ≥ ~180 live days (2013: from GNW, DOSE, SP2, D03D, B05D, WISH, NLWA). TA.A04D is left out because it has no
  entry in the scan coordinate table. Grid templates for the extra stations use the same synthetic recipe (scan_calib
  node_templates), appended to the v1 grid. Gain fallback g = 1 (Z only) via SCAN_GAIN_FALLBACK=1. The injected C0/C1
  signals at the extra stations are amplitude-calibrated stacks from their 2011–12 records (same rule as P1). Shared-code
  edits are env-gated with unchanged defaults (pre-edit copies in `data/*_before_v2.py`).

## Amendment S (2026-10-15 ~7:45 pm PT, before any of its numbers) — 2018 with scan v2, 2011 carriers, production 2010–2018
Recorded first: pilot 2018-06 under R with the v1 scan passed T1′ (18 vs λ 0.036), T2′ (T2a 0.94, T2b 27 vs 10.7, p = 0.004),
T3′ 3/3 and controls (0.34/yr); its completeness at Mw 3.8 was 0.27 (ETS) / 0.22 (between). Scan v2 on 2013 passed its
acceptance (coordinator's relay of `pilot2013_R_verdict_v2.json`).
- **(a) 2018 with scan v2**: same pilot, battery and injections (n = 800) re-run on the 2018 v2 scan (extra stations as the
  v2 rule). Reported: completeness at Mw 3.4–4.0 with Wilson CIs, MF/battery split, T1′–T3′, controls. Rule (Merlin): if
  completeness at Mw 3.8 stays < 0.5 (either period) with v2, 2016–2018 are still produced but flagged "relay templates
  needed".
- **(b) 2011 carriers usable**: (1) the hg_repro 1-sps displacement day files of UW.SP2, TA.A04D, TA.D03D (same v3 recipe,
  hg_prep) are linked into `data/disp1sps_v3` where a v3 day is missing (links, nothing overwritten); (2) the newly fetched
  raw days are prepped with `scan_prep_ref.py` (skips existing days); (3) the 2011/2012 reference bp arrays of these stations
  are rebuilt into new directories `bp2_ref_v2` (the old arrays, used by E2–E4, stay untouched) and production reads them
  via an env override; (4) the catalogue MF arrays for WISH/NLWA are rebuilt in a new frame directory (old ones kept).
- **(c) Production 2010–2018** (one queue, per-year markers `PROD_<year>_DONE`, FAIL lines in `logs/q_main.log`):
  - Frame = calendar year; MF stations = the 2011-bank carriers (dense2011 stations + PFB) with ≥ 180 broadband days that
    year (station-era table); bank T1–T5 + C0–C5 cut from 2011–12; weighted statistic; 25 noise templates per year (random
    times inside that year's SVI episodes, else the year's first 60 days if no episode); 8 + 8 scrambled controls.
  - Scan v2 per year (v1 scan year built first where missing: 2010, 2015, 2016, 2017).
  - Tier = MF ≥ 5.5 MAD + two-template agreement + battery′; Tier A′ for between-ETS events (R definition).
  - Periods: each SVI episode of the year ±2 d (`episodes_b.csv`; 2010 from the same file) = "episode <id>"; NW-only episodes
    = "other episode"; rest = between.
  - Per year: tier counts per episode and between (raw, Tier A′) with both controls (scrambled mean/max, noise pairs);
    summed moment per episode from the scan M0 (m0_cols); storm-day and cluster-day columns; completeness from C0/C1
    injections, n = 800 per Mw per period type (all episodes pooled vs between; Mw 3.4–4.0), Wilson CIs, MF/battery split,
    Mw_c; per-episode completeness from the injections falling inside each episode (with CIs).
  - Known-answer ledger: the 14 vetted + near-miss 2011–12 events, the 2013 pilot ETS events and the 2018 pilot ETS events —
    found / not found by production (±30 s), with tier, MF z and battery columns.
  - Order 2011, 2012, 2013, 2018, 2014, 2015, 2016, 2017, 2010. 2019–2026 and 2007–2009 held.

## Amendment T (2026-10-15 ~8:52 pm PT, before any production number and any relay number)
### T0. Production fixes (no result seen)
- First production attempt failed at 2011 "cont" (8:02 pm PT): the pilot CFG dict was indexed by year before the
  production override (KeyError 2011); leap years also had 365 days except 2008. Fixed: CFG is built per production year
  from `episodes_b.csv` (every SVI episode ±2 d, majors flagged) and `tables/prod_stations.csv` (2011-bank carriers with
  ≥ 180 v3 day files; carriers with ≥ 180 broadband days in the station-era table but < 180 prepped days are prepped first:
  2014 WISH/NLWA; 2010 and 2015–2017 after their scan years). Episodes are as pre-registered in S (all SVI episodes; the
  per-episode table carries a "major" column). Queue restarted 8:51 pm PT; finished stages skipped.
### T1. 2011 production end-to-end check (before later years are trusted)
- Known answers: the 8 vetted + 6 near-miss 2011–12 events present in the 2011 or 2012 production at tier level
  (found ±30 s); the 2011-08 episode tier count compared with the E2/E3 numbers for the same window (E2 chosen tier:
  2011 ETS part of 18 two-template e′ events; E3 dense, self-detections included). Reported as a ledger, not a gate; a
  shortfall of known events is investigated before later years are reported.
### T2. Relay templates for 2016–2018 (Merlin's rule; trigger: 2018 v2 completeness at Mw 3.8 = 0.29, MF-only 0.42)
- Relay source events per year Y ∈ {2016, 2017, 2018}: Y's production ETS tier events (2011 bank, scan v2) that also have
  LFEs within ±10 min and vr_thrust ≥ 0.12 ("vetted relay events").
- Relay bank: 90-s windows (origin − 9 s, V4 recipe, velocity) on every live broadband of Y with v3 data; complete-linkage
  clusters (15 km, best node); ±5-s aligned stacks; per-event leave-one-out stacks; weighted statistic.
- **Acceptance (Merlin): leave-one-out re-find ≥ 85 % of the vetted relay events (by a relay template not containing them, at
  0.5 noise detections per template-day from 25 noise templates of Y) AND scrambled (8 + 8) battery′ tier rate between ETS
  ≤ 1/yr.** Fewer than 3 vetted relay events → no relay bank for Y (stated).
- If accepted, Y is re-run with bank = 2011 bank + accepted relay bank (relay templates validated by LOO count as
  validated). Completeness reported twice: C0/C1 injections (comparable to other years) and relay-cluster injections
  (each relay cluster stack injected, detected by the bank without that cluster), n = 800 per Mw per period type.
- Queued after the whole 2010–2018 production queue (relay sources need the 2016–2017 production).

## Amendment U (2026-10-16 ~3:22 am PT, before any relay re-run number) — production with the 2018 relay bank (2016–2018)
Recorded first (relay step, T2 rule): 2018 accepted (LOO 15/17 = 0.88 ≥ 0.85; scrambled 0.34/yr), 2017 rejected (3 vetted,
LOO 0/3), 2016 no bank (0 vetted events).
- **Bank** for 2016, 2017 and 2018 = validated 2011 bank (T1–T5, C0–C5) ∪ the accepted 2018 relay primaries (RC0–RC7,
  their 2018 windows; channels of stations not live in 2018 are absent). All count as validated for agreement and tier.
- Everything else as Amendment S: each year's live broadbands (2011-bank carrier rule), scan v2, weighted sum, 25 noise templates
  of the year, 8 + 8 scrambled controls, tier = MF ≥ 5.5 MAD + two-template agreement + battery′, Tier A′ and Q2 triage
  for between-ETS events, per-episode tables, ledger. Outputs separate (`data/prodr<Y>/`, `tables/prodr<Y>_*`).
- **Completeness for the bank actually used**: (i) C0/C1 injections as before (400 per source per Mw per period), detected by
  the bank without the source's templates; (ii) relay-cluster injections: each 2018 relay cluster's amplitude-calibrated stack
  (MF channels from the 2018 frame, scan stations from the 2018 v2 bp arrays; Mw_ref from the members' scan M0, log-mean of
  M0), 800 per Mw per period split evenly over the clusters, detected by the bank without that cluster. Both reported
  (Wilson CIs, MF/battery split); Mw_c per group.
- In 2018 the relay templates were cut from 2018 events, so 2018 ETS counts include self-detections of the 17 relay events
  (stated next to the count; the LOO result is the honest benchmark).
- Markers `PRODR_<Y>_DONE`; FAIL lines in `logs/q_main.log`.
- Amendment U addendum (~3:30 am PT, before any relabelled number): `hg_r3.lfe_lists()` also clipped the LOWER end at
  2011-05-31, so every LFE/tremor label before 2011-06 (all of 2010, Jan–May 2011) was void ("quiet"). Fixed to 2003-01-01;
  2010 and 2011 summaries are re-scored (`scripts/relabel_queue.sh`, marker RELABEL_DONE; previous versions kept as
  `*_preRelabel.*`). T2b (which reads the SVI rows directly) was not affected; T2a and the "with LFEs" counts were.

## Amendment V (2026-10-16 ~3:32 am PT) — PLAN ONLY (not run): relay chain 2019–2026 and 2007–2009
### V1. 2019–2026 (forward relay chain)
- Start bank B(2018) = 2011 bank ∪ accepted 2018 relay bank. For each year Y = 2019, …, 2026 in order:
  1. scan year Y (v1) + scan v2 (all live broadbands with bp files; g = 1 Z-only fallback);
  2. production Y with B(Y−1) (Amendments S/U recipe; weighted; n = 800 per Mw per period; C0/C1 + all accepted relay
     cluster injections);
  3. vetted relay events of Y = ETS tier events with LFEs ±10 min and vr_thrust ≥ 0.12; if ≥ 3, build Y's relay bank on
     Y's live broadbands and accept it only with LOO ≥ 85 % and scrambled ≤ 1/yr;
  4. if accepted, B(Y) = B(Y−1) ∪ relay(Y) and Y is re-run with B(Y) (completeness for the bank used); else B(Y) = B(Y−1).
- Stop/flag rule: if 2 consecutive years have completeness < 0.5 at Mw 3.8 in both periods AND no accepted relay bank, the
  chain is broken there; later years are produced but labelled "bound only".
- Prerequisites to check before running: station epochs and gains after 2019 (HH replacements at PGC/VGZ/SNB/NLLB; the
  scan had 4 stations after 2018), sensor rule (≥ 30 s), far references 2019–2026 (on disk).
### V2. 2007–2009 (backward relay via 2010-08; Merlin review 3 §6)
- A) 2010 production exists (2011 bank; 7 stations; 2010-08b1: 8 tier events at the patches).
- B) Relay bank from 2010-08b1 events on the 2010 network extended with the CAFE broadbands (XU BS11/W020/W040/PL11/W060, CMG3T
  in 2006–2010) and SQM (broadband to 2011-07): vetted relay events = 2010-08b1 tier events with vr_thrust ≥ 0.12 (LFE
  labels after the Amendment U fix); accept with LOO ≥ 85 % and scrambled ≤ 1/yr, validated additionally on the
  11-station 2008 subset (stations live in both 2008 and 2010).
- C) Scan extension for 2007–2009: CAFE + SQM in the scan station list (grid templates), noise rule referenced to each
  station's own year (the 2011-referenced rule excluded VGZ/SNB in 2008), g = 1 Z-only fallback.
- D) Pilots 2008-05, then 2009-05 and 2007-01, with the Amendment R rules (T1′, T2′, T3′, controls gate; completeness reported).
  GO → production 2007–2009 with the relay bank; NO-GO → bound only for 2007–2009.
- No new downloads are needed for V1 or V2 (CAFE/SQM data on disk); the TA 2006–08 LH download stays deferred.

## Amendment W (2026-10-16 ~6:05 am PT, written before the diagnostics below) — 2016–2017 control failure, the 2016-02-25 candidate
### W1. Diagnostics (descriptive; no rule depends on them except through W2)
- For 2016, 2017 and a reference year with passing controls (2013):
  - per station-component, per month: median 10-min RMS of the MF velocity data, and glitch rate = fraction of 10-min
    windows with RMS > 10× that station-component's yearly median;
  - v3 day coverage per station (PFB, SNB, YOUB, LZB, NLWA status) and the MF channel weights;
  - for every scrambled-control tier detection between ETS: template, draw, time, the scan battery columns, and the
    per-station share of the summed MF statistic (local recomputation on the scrambled data);
  - GNW monthly noise.
### W2. Fix-selection rule (fixed now, before any re-run)
- A fix is adopted only if it is a data-quality rule that (i) is defined without reference to control or real detection
  counts, (ii) applies to every production year, and (iii) changes the passing reference year 2013 by ≤ 10 % in tier events
  and ≤ 0.05 in completeness at Mw 3.8. Candidates, in order:
  1. **Station-month QC for the MF**: a station-component-month enters the MF only if its median 10-min RMS ≤ 3× the
     station-component's median over all its production months AND its glitch rate ≤ 2 %.
  2. **Year-own threshold**: tier MF threshold = max(5.5 MAD, the level at which that year's 25 noise templates fire
     0.5 per template-day).
- The re-run (2016–2017 with each bank, plus 2013 as the check) is done only for a candidate that the W1 diagnostics show
  to act on the dominant cause; the outcome numbers are reported whatever they are. If neither candidate addresses the
  cause, 2016–2017 stay flagged and no re-run is made.
### W3. 2016-02-25 12:11:21 PST candidate
- Waveform figure (near, regional, 5 far LHZ), full battery columns, storm percentile, LFE/tremor ±10 min / ±1 h,
  location; Q2 class as written (vr_thrust ≥ 0.12, far/near < 0.5 with near = the 2 nearest live scan stations to avoid a
  single glitchy station setting the near maximum — this near definition is new and is stated here before the figure).
### W4. Forward relay chain 2019–2026 (Amendment V1) starts only after W1–W2 are resolved; 2007–2009 held.
### W1 result and W2b (2026-10-16 ~6:08 am PT, after the W1 diagnostics, before any re-run; disclosed as post-diagnostic)
- W1 found (tables `diag_noise_2013/2016/2017.csv`): **UW.LRIV is broken in 2016–2017** — Z median 10-min RMS 420–480 nm/s
  (normal 2–5) in 2016 Jan–Apr, Dec and 2017 Jan–Jun, Oct–Dec; LRIV is a heavily weighted near station. Most 2017 scrambled
  tier detections (8 of 10 events) and half of 2016's (6 of 13) fall in LRIV-broken months; they cluster on a few days
  (2016-02-25: 4 control events; 2017-01-18 …). Polarity-draw leakage (sign correlation with the true template) does not
  explain the counts. Other stations show only the usual seasonal noise.
- The W2 candidate 1 as written cannot pass W2(iii): its 2 % glitch clause (> 10× yearly median) is exceeded by 2–11 % of
  windows at every station in every year, including 2013. Candidate 2 (year-own threshold) raises the threshold in every
  year alike (6.3–7.2 MAD), so it does not act on the cause. Neither is adopted.
- **W2b (new candidate, defined from data quality only): station-component-month QC on RMS alone** — a station-component-
  month enters the MF (real, controls, noise templates and injections alike) only if its median 10-min velocity RMS
  ≤ 3× that station-component's median of monthly medians over 2011–2018. Computed from `diag_noise_<Y>.csv` for all
  production years. W2(ii)–(iii) apply: it is applied to 2016 and 2017 (both banks) and to 2013 as the check (≤ 10 %
  change in tier events, ≤ 0.05 in completeness at Mw 3.8); the results are reported whatever they are, and W2b is
  adopted for production only if the 2013 check passes.
- W2b implementation fix (2026-10-16 ~7:45 am PT, after the 2013 check had run, before any 2016/2017 QC number): in QC mode
  the noise-template window search required every "needed" channel to be live; with LRIV/NLLB/NLWA months excluded no 2016
  window could qualify, so the search never ended (bank step ran 60 min; killed by PID 374447). Now a needed channel that
  is QC-excluded at that time is not required, the per-channel masks are cached, and the search stops with an error after
  200,000 tries. This rule only differs where a needed channel is excluded (2013's check bank was built with the old rule;
  2013 exclusions: 3 PGC, 1 LZB component-months).
- **W2(iii) 2013 check result: tier events 18 → 19 (+5.6 %, ≤ 10 %); completeness at Mw 3.8 0.831 → 0.831 (episodes),
  0.826 → 0.779 (between, −0.047, ≤ 0.05); controls unchanged (0.725/yr). → W2b passes the check and is adopted.**

## Amendment X (2026-10-16 ~9:16 am PT, before any 2019–2026 number) — forward relay chain 2019–2026 (user: extend to 2026 now)
Implements Amendment V1 with the adopted W2b noise filter. Runs in parallel with the W re-runs; CPU split 4 + 4 workers.
- Bank B = 2011 bank ∪ accepted relay banks, starting from {2018} (`data/fwd_bank_years.txt`). Relay templates of year y
  keep their own-year windows; names are prefixed `R<y>` (2018 unprefixed) and members `<y>:` so that leave-one-out and
  injection exclusions stay per year.
- Per year Y = 2019, …, 2026: scan year Y (v1 recipe; prep of the scan, regional, far and new broadband stations
  UW.OTR [broadband from 2020-10] and UW.HDW [from 2023]; ComCat regional list); MF stations = 2011-bank carriers plus OTR
  and HDW when they have ≥ 180 prepped days (they carry relay templates only); W2b exclusions for Y (reference = 2011–2018
  channel median, or the year's own median for channels without 2011–2018 data); scan v2; production `prodf` with B; tier,
  Tier A′/Q2, controls, per-episode tables, injections n = 800 per Mw per period (C0/C1 400 each per period + all accepted
  relay clusters together 800).
- Relay: vetted relay events = Y's ETS tier events with LFEs ±10 min and vr_thrust ≥ 0.12 (≥ 3 needed); relay bank on Y's
  stations; accept with LOO ≥ 85 % and scrambled ≤ 1/yr; if accepted, B ← B ∪ {Y} and Y is re-run (`prodfr`) with the
  extended bank (the year's final numbers are the prodfr ones).
- Failed year = completeness at Mw 3.8 < 0.5 in both periods with the final bank AND no accepted relay bank for that year.
  **Stop after 2 consecutive failed years**; later years are not run.
- 2026 is the year to 2026-10-15 (data on disk). 2007–2009 held.
- Markers PRODF_<Y>_DONE, RELAYF_<Y>_DONE, FWD_DONE; FAIL lines in `logs/q_main.log`.

## Amendment Y (2026-10-16 ~11:28 am PT, before any number of its steps) — Merlin review 4
Merlin's verdict adopted: no per-year control-tuned threshold anywhere.
### Y1. 2018 under W2b (LRIV excluded all 12 months)
- Re-run 2018 with the 2011 bank under W2b (`prodq2018`). Then the vetted relay events of prodq2018 (ETS tier, LFEs ±10 min,
  vr_thrust ≥ 0.12) give a new 2018 relay bank cut under W2b (no LRIV channel), accepted only with LOO ≥ 85 % and
  scrambled ≤ 1/yr (T2 rule). If accepted, 2018 is re-run with 2011 ∪ new relay (`prodrq2018`; relay-cluster injection
  sources rebuilt from the new bank, the previous RC sources kept as a copy).
- Accept: 2018 T1 and T3 verdicts unchanged from the pre-W2b 2018 runs (T1 pass, T3′ 3/3); the new control rate is
  reported whatever it is. Failure signature (Merlin): T1/T3 lost → the LRIV channel was carrying the 2018 result.
- Forward chain: held since 11:27 am PT (scan-year 2020 build allowed to finish). The 2019 outputs made with the pre-W2b 2018
  relay bank are kept as `*_preW2bBank*` and redone. The chain restarts with bank = 2011 ∪ the new 2018 relay bank if it
  is accepted, else with the 2011 bank only until a new relay year passes. The forward-chain bank is capped at
  2011 ∪ the latest accepted relay bank (Merlin flaw 5), not the union of all.
### Y2. Final 2016–17 = the 2011-bank W2b runs (prodq2016, prodq2017); the relay runs are columns (they fail the T2
  relay rule's scrambled ≤ 1/yr in those years).
### Y3. Controls: pooled rate replaces the per-year gate
- Final 2011-bank runs per year (2010–2012, 2014–2015 as produced; 2013, 2016, 2017, 2018 under W2b): scrambled-control tier
  events between ETS, deduplicated per draw (±30 s across templates). Pooled rate = Σ events / Σ (16 × live between-ETS
  years) with a Poisson 95 % interval; A1/A2 split columns; per-year rates as columns.
- Dispersion test: χ² of per-year counts against the common rate × (16 × live years), df = years − 1. p > 0.05 → per-year
  differences are noise and 2016–17 are plain low-completeness years (Mw_c > 4.0, stated); p ≤ 0.05 → 2016–17 flagged
  "elevated noise" and still reported.
### Y4. Tier A′ on the controls
- Q2 is applied to every scrambled-control tier event, with the scan columns at the control time (real data, as in the
  battery). far/near = max |Z| over the 5 far LHZ / max |Z| over the 2 nearest live scan stations in [t0, t0 + 180 s];
  vr_thrust ≥ 0.12; ratio_GNW ≤ 0.551 north of 48.1 N; storm day and cluster day (≥ 2 events of the same draw on a UTC day)
  as for real events. Tier A′ control rate per year and pooled. Accept: pooled ≤ 0.3/yr; otherwise the between-ETS bound
  is stated at the raw-tier false rate.
### Y5. PNSN tremor within 30 km of the C0/C1 centroids during 2017-03 (2017-03-16 – 04-07), vs other patch majors; if
  abundant, 2017 T2 = "untestable" (LFE-catalogue gap), not failed.
### Y6 (optional, < 10 min): LRIV PSD shape, broken vs normal months (flat ratio → gain/metadata; red → sensor fault).
### Y7. User decision (2026-10-16 ~12:05 pm PT): follow the published VLFE practice
- The primary false-detection measure is the **noise-template control** (H&G 2019 style: templates cut at random quiet
  times; our noise-pair tier: noise template j with agreement by noise template j+1, same battery′). The scrambled controls
  (polarity, station-day shift) are columns only and no longer a pass/fail gate in any year, including relay acceptance.
- Relay acceptance = LOO ≥ 85 % (with ≥ 3 vetted events); the year's noise-template false rate is reported next to it; the
  scrambled rate is reported, not gated. (The running 2018 re-run is not restarted; its relay verdict step reads this rule.)
- Y3/Y4 report, per year and pooled: noise-template raw tier rate (per noise-template year), noise-template Tier A′ rate
  (Q2/Tier A′ applied to noise-template tier events), and the scrambled rates as columns. **The between-ETS bound is
  stated against the noise-template Tier A′ rate.**
### Y8. User instruction (2026-10-16 ~12:10 pm PT): rank, do not discard — catalogue product
- Per year `tables/catalogue_ranked_<year>.csv` and a combined `catalogue_ranked_2010_2026.csv` (years appended as they
  finish): every REAL detection of the validated bank above the MF MAD threshold (5.5 MAD, the only hard cut), deduplicated
  to one row per event (±30 s across templates; the row with the highest MF z represents it; n_templates and agree list as
  columns). Every check is a column: zmad, n_templates, agree, z_best, ratio_far, n_far, ratio_GNW, vr_thrust, z_loso/z_best,
  e and e′ masks (head2 magnitude), footprint/location, plane-wave flag, nearest family and distance, LFE status and count
  ±10 min, period/episode; the "strong" tier is a boolean label only. Columns that the battery computes only for times
  passing (c) and (e or e′-time) are empty elsewhere (stated).
- Combined score: logistic regression of known VLFEs (2011–12 vetted + near-miss, 2013 pilot ETS tier events; positives)
  against noise-template detections 2011–2013 (negatives) on standardised check columns (missing values filled with the
  failing value plus a missing indicator), class-balanced, L2 (C = 1). Coefficients and cross-validated AUC reported. The
  score is applied unchanged to every year.
- Curves: count of real events with score ≥ s, against the noise-template control (events per noise template, i.e. one
  catalogue's worth) as the primary reference and the scrambled controls (mean per draw) as a column, for ETS and between
  ETS, per year and pooled. The between-ETS statement is read off these curves.
- Existing tables are not changed.
### Y8b (2026-10-16 ~1:33 pm PT, before any new ranked number) — the Y8 score was wrong; fix
- Fault: the battery fills ratio_far / n_far / vr_thrust / z_loso only for times passing (c) and (e or e′-time); 90 % of real
  rows lacked them, the fill values plus an unconstrained logistic fit (32 positives) let missing values RAISE the score
  (grid-edge rows with no far ratio scored 1.0). Y8 outputs are kept but superseded.
- Fix: (1) the expensive columns are computed for EVERY row with two-template agreement, identically for real rows,
  noise-template rows and scrambled rows (scrambled: draws A1_0, A1_1, A2_0, A2_1 only, counts scaled per draw; all 25 noise
  templates), with the same scan-v2 battery as the year's run; rows without agreement keep missing values (they cannot rank
  high because agreement is a feature). (2) Score = logistic regression with every feature oriented "higher = more
  VLFE-like" and coefficients constrained ≥ 0 (L-BFGS-B, L2, class-balanced); missing values are set to the worst value of
  their feature, so a missing check can only lower the score. Features: agreement (0/1), MF z − 5.5, log10 scan z, far test
  pass (n_far ≥ 2 and ratio_far ≤ 0.551), −ratio_far (missing → −3), vr_thrust (missing → 0), z_loso/z_best (missing → 0),
  e′ pass, location pass (c), not travelling (missing → 0), not on grid edge. Columns n_checks_missing, grid_edge, strong_tier
  (raw tier), tierA (Q2 rule as Y4: vr_thrust ≥ 0.12, far/near (2 nearest) < 0.5, GNW clause north of 48.1 N, not storm+cluster).
  (3) Sanity: rank of the 32 known 2011–13 VLFEs; between-ETS real vs noise-template at high score; if there is an
  excess, the top 10 with all columns. Outputs `catalogue_ranked2_<Y>.csv`, `catalogue_ranked2_2010_2026.csv`,
  `ranked2_curves_*.csv`, `figures/ranked2_curves_<Y>.png`; marker Y_RANKED2_DONE.
- Y8b detail (1:40 pm PT, written at launch, before any number): positives = the 8 vetted 2011–12 events (e3_events role
  "vetted") + the 18 2013 ETS events (real rows within 30 s), negatives = all 2011–13 noise-template rows; the 6 "near-miss"
  events are not positives and their rank is reported separately. Because positives are in the fit, their rank is in-sample;
  5-fold CV AUC is reported alongside. Cluster day for tierA = ≥ 2 Q2 events (30-s declustered) on the same day in the
  same catalogue (real; one noise template; one scrambled draw). Queue scripts/y_ranked2.sh (CAT_NW = 2).

### Y8c (2026-10-16 ~2:45 pm PT, after the Y8b numbers; method change after seeing data, flagged for Merlin) — score fix 2
- Fault in Y8b: the monotone logistic is separable (26 positives vs ~70k noise-template rows; intercept −98.7, coefficients
  up to 5), so 1,088 real rows sit at score ≥ 0.9, many at 1.0, and the score cannot order the top.
- Fix: (1) same features and monotone constraint; the L2 strength is chosen by 5-fold CV (grid 1e-3…1e3 per row,
  minimum held-out class-balanced log-loss). (2) The ordering score is the logit (no sigmoid saturation); the class-balanced
  probability is a column only. (3) Calibration: out-of-fold probability by decile (positives vs predicted), and the
  noise-template false rate per template-year at fixed logit thresholds estimated on even templates and checked on odd
  templates. (4) Per-row false-rate columns: noise-template per template-year, the same × 11 (bank of 11 templates), and
  scrambled per bank-year. (5) Variant score (column only): negatives = noise-template rows + an equal number of scrambled
  rows, because the noise-template control has only one agreement partner per template while the real bank has 10 similar
  partners (RESULTS §20 (a)). (6) Curves: real, noise template per catalogue and × 11, scrambled A1 and A2 per draw.
- Outputs `catalogue_ranked3_*`, `ranked3_curves_2010_2026.csv`, `ranked3_calibration.csv`, `ranked3_cv_l2.csv`,
  `ranked3_sanity.json`, `figures/ranked3_curves_pooled.png`; diagnostics `r20_*`. Y8b outputs are kept.

### Y9 (2026-10-16 ~2:45 pm PT, before any Y9 number) — following Merlin review 5 (`MERLIN_CAT_REVIEW5.md`) and the user's no-invented-gates rule
- Tier A′/Q2 is retired as a gate and as the carrier of the between-ETS bound (it passes 3 of the 26 known VLFEs); vr_thrust,
  far_over_near2 and ratio_GNW stay as columns. The learned score is not the ordering; `score_logit` stays as a column labelled
  "in-sample" (Y8c features, loso fill changed from 0 to −1, the clip floor).
- Step 1 (Q2 completeness): from `data/<run><Y>/inj_v2.parquet`, the fraction of battery′-recovered injections (MF best and
  second ≥ 5.5 and battery′) with vr_thrust ≥ 0.12, per year, period and Mw 3.6/3.8/4.0. < 85 % → Tier A′ retired as
  bound-carrier (columns kept); ≥ 85 % everywhere → re-examine the GNW clause.
- Step 2 (ladder): among agreement rows, level = number of battery′ checks passed (a far test, b not travelling, c located,
  e′ quake test, f leave-one-station-out); a check that cannot be measured (n_far < 2, travelling / z_loso / z_best missing)
  counts as not passed and is counted in `n_checks_missing`. Order: level, then zmad. Flag columns: eq_in_window (fails the
  strict e: head2 / bw5 / reg / reg3 window), grid_edge, storm_pct_max_pm1d (max of the PGC 5–10 s percentile over the day
  and ±1 day). Mechanical checks: no in_head2 or grid-edge row in the top 50 between ETS; all 26 known VLFEs at the top level
  their n_checks_missing allows.
- Step 3 (ladder tables): per year and pooled, ETS and between, for each level ≥ L: real events, distinct days, noise-template
  control (per catalogue × 11, placeholder until step 4), A1 and A2 per draw, completeness at level ≥ L from inj_v2
  (Mw 3.6/3.8/4.0), and the excess (real − control) with a 95 % Poisson (Garwood) interval on the real count.
- Step 4 (similarity-matched noise banks, 2011–2013): the real bank's mean inter-template CC (11 templates, per channel,
  max over ±30 s lag, mean over channels and the 55 pairs). For each of 8 base noise templates per year, 10 siblings
  = normalise(a·base + √(1−a²)·Q_j), Q_j another quiet window drawn like the noise templates, a set so the bank's mean
  pairwise CC equals the real bank's (achieved value reported). Each bank (base + 10 siblings) runs the same MF (base noise
  weights, same channels, ±90 s around every source window excluded), the same agreement (≥ 2 bank members within 30 s,
  MF ≥ 5.5) and the same battery′ (scan v2 env of the year). Accept the banks if their between-ETS agreement rate is within
  ×1.5 of the real bank's. Outcome A: level-5 between-ETS control per bank ≥ 20 (pooled 2011–13 scaled to 2010–18) → the
  excess is bank similarity; Outcome B: ≤ 5 → the excess stands.
- Step 5 (emergence + amplitude–distance): for every agreement row, per station window peak (0–180 s, 0.02–0.05 Hz Z) /
  median envelope (|Hilbert|) over ±30 min; columns: median near emergence, n near stations with emergence ≥ 2, max/second
  station emergence ratio (spike), and the §12 T3 quantities (near/mid ratio, far/near-max, Spearman of amplitude vs
  distance). Calibrate the emergence threshold on the 26 known VLFEs (≥ 85 % at or above it) and report the share of
  2014-01-23 rows above it. Figures (r20 style) for every level-5 between-ETS real row, each labelled wavefield / spike /
  candidate by eye in a table, with the labelling rule stated before looking: spike = one station carries the window
  (max/second emergence ≥ 3 or one station > 3× all others); wavefield = window peak not above the surrounding ±30 min
  (median near emergence < threshold) on most stations; candidate = neither.
- Markers Y9_STEP1_DONE … Y9_STEP5_DONE, Y9_DONE; FAIL lines → logs/q_main.log. CPU ≤ 2 workers (forward chain has 6).

### Y10 (2026-10-16 ~3:20 pm PT; user's rule, relayed by the coordinator) — no stop rule in the forward chain
- The Amendment X rule "stop after 2 consecutive failed years" is removed. It fired on 2019 and 2020 (FWD_DONE at 3:10 pm PT)
  although both years pass T1/T2/T3 with GO = True under Y7; "failed" was the combined completeness at Mw 3.8 < 0.5, which
  averages in the 2018-relay injection sources (2020: 0.0; 2019: 0.27) while C0/C1 alone give 0.58/0.69 (2019) and
  0.65/0.61 (2020).
- Every year is catalogued with its completeness stated per injection source (C0/C1 and each relay set); a low-completeness
  year is flagged in the log and the tables, never skipped. Bank for 2021+: 2011 ∪ 2018-W2b relay unless a later relay is
  accepted (Y7). W2b, Y7 and Y9 rules unchanged. The old FWD_DONE is kept as FWD_DONE.stopRule2020.

### Y10b (2026-10-16 ~3:40 pm PT; bug fix approved by the coordinator — W2 already says W2b applies to every MF part, injections included)
- The relay injection sources (`data/inj/sources.npz`, RC*) are rebuilt with W2b masks: a member window in a W2b-excluded
  station-component-month is not used for that channel (scan-side: per component). `cat_inj.run_relay_sources` now does this for
  every later relay too. For the 2018 relay set this removes all UW.LRIV keys (32 keys; LRIV excluded all of July 2018);
  every other key is unchanged. Old file kept as `data/inj/sources_before_Y10b.npz`.
- Re-run inj + summary + R + psum only, for prodrq2018, prodf2019, prodf2020 (`scripts/y10b_queue.sh`, after Y9_STEP4_DONE,
  2 workers; old outputs in `tables/preY10b/` and `data/<run>/inj_v2_preY10b.parquet`; marker Y10B_DONE). The forward chain
  reads the fixed sources from 2021 on (its first inj step is hours away). Relay templates kept.
- UW.WISH: July 2018 passes the 3× rule (MF-band median RMS 0.56–0.73× its 2011–2018 median on all components), so it is not
  masked and no rule is added. Its large relay-source amplitudes (RC7 N 2.8 µm, E 1.7 µm) are horizontal only (Z 32 nm, normal;
  scan-band Z at the relay events 21 nm, emergence 3.7); the scan-band horizontals at the relay events are 0.38–0.50 µm,
  i.e. long-period horizontal noise stacked into sources built from 2–4 members. Reported, not changed.

### Y9b (2026-10-16 ~5:06 pm PT, before any Y9b number) — the 5 between-ETS candidates; matched noise banks 2014–2017
- Candidates (Y9 eye labels): 2014-05-05 08:30, 2015-04-13 06:06, 2016-06-13 07:41, 2016-06-13 19:58, 2017-01-31 04:22 (PT).
  For each: (i) catalogued earthquakes ±5 min, any magnitude, from the local copies (PNSN/ANSS ComCat regional lists in the
  scan eq tables: reg300 M ≥ 2.5, reg1500 M ≥ 4, reg2500 M ≥ 3, global M ≥ 4.5); if no local list has smaller events, say so
  (no new download without approval); (ii) high-frequency energy at the near stations from the raw traces: RMS in 2–8 Hz and
  > 1 Hz in the pulse window (t0 − 30 s … t0 + 150 s) over the median of 3-min windows in ±30 min (an ordinary earthquake
  shows an HF excess at the pulse; a VLFE is HF-depleted); (iii) PNSN tremor within 30 km and ±1 h, and SVI LFE (P ≥ 0.75)
  and LFE-candidate (0.1 ≤ P < 0.75) detections within ±1 h, with family distance ≤ 30 km where family locations exist;
  (iv) the scan grid location and its offset from the C0/C1 patches. Class: regular earthquake (catalogued, or HF excess ≥ 3
  on ≥ 2 near stations at the pulse), VLFE-like (no catalogue event, no HF excess, LP pulse decaying with distance), else
  unclear. One figure per event.
- Matched noise banks for 2014 (prod), 2015 (prod), 2016 (prodq), 2017 (prodq) with the Y9 step-4 method unchanged
  (8 banks per year, base + 10 siblings, mean CC matched to the year's real bank, same MF/agreement/battery′) and the same
  acceptance rule (agreement within ×1.5 of the real bank). Report per year and pooled 2011–2017: real level-5 between-ETS
  count and rate vs accepted-bank control rate (95 % Poisson), and whether the years with candidates sit above control.
  Queue scripts/y9b_nb.sh, ≤ 2 workers, starts after Y10B_DONE; marker Y9B_DONE.
- Correction noted before numbers: the list `svi_c010_075` is LFE candidates (0.1 ≤ P < 0.75, all SVI families), not
  "C0/C1-patch LFEs"; the 2021 counts reported at ~5:05 pm PT were mislabelled and are redone here as true patch counts.

### Y11 (2026-10-16 ~5:12 pm PT, before any Y11 number) — following Merlin review 6; f′ is a post-diagnostic rule change, disclosed (as W2b was)
- Step 1: injection battery outputs now carry `loso_dropped`, `max_share`, `max_share_ch` (added to cat_inj.battery_local /
  cat_pilot.inj_row; pure additions). Mw 3.8 injections of prodf2021 re-run once with these columns
  (`data/prodf2021/inj_y11_m38.parquet`); tabulated by month and source. Accept: ≥ 80 % of C0 Mw 3.8 rows in Jul–Dec drop
  UW.LRIV and relay rows drop CN.PGC; if a regional station dominates C0, it is a source bug and f′ waits.
- Step 2: f′ = (z_loso ≥ 3 AND z_loso ≥ 0.5·z_best) OR z_loso ≥ 5 ("the network without its strongest station alone reaches
  the (c) location threshold"). Re-scored from existing outputs only (bat_v2, r2bat_<Y>, inj_v2, y9nb bat): (i) the 26 known
  VLFEs; (ii) pooled 2010–2018 between-ETS level-5 counts, real vs noise template ×11, matched noise banks (accepted, scaled),
  A1 and A2, under f and f′; (iii) C0 and C1 completeness at Mw 3.8 by year and month, 2013–2021, under f and f′; (iv) every
  real between-ETS row that f′ promotes to level 5, figured and labelled with the Y9 emergence rule (plus eye check).
  **Adopt f′ if each control's level-5 count rises by ≤ 1 (pooled 5.54 yr) AND C0 completeness at Mw 3.8 rises by ≥ 0.15 in
  the months where f < 0.5; otherwise reject** (then 2021+ is reported bound-only and the station route is the fix).
- Step 3 (only if adopted): re-derive the 2019–2021 tier, episode and relay-vetting tables under f′ (summary/R/relay, no MF
  re-run); 2022+ chain outputs re-scored the same way afterwards. Accept: 2021-10 C0 completeness ≥ 0.4 at Mw 3.8.
- Step 4: ladder extended to 2019–2021 (then each new year). Step 5: metadata-only FDSN station query (broadband, corner
  ≥ 30 s, ≤ 100 km of 48.3 N 123.6 W, ≥ 180 d/yr in 2021–2026; CN, UW, CC, NV) and a CN.GOBB PSD/response check against
  CN.PGC from local data. No waveform downloads. Horizontals pilot deferred (coordinator).
- Step 7 reporting 2021–2026 per episode: tier count; level ≥ 4 count with f status; C0/C1 completeness at Mw 3.6/3.8/4.0 with
  Wilson CI; implied-rate interval (count + Poisson 95 %)/completeness; label "bound only" when C0/C1 completeness at Mw 3.8
  < 0.2 (label only); relay-set completeness as a column marked "not comparable (sources lack LRIV)". Marker Y11_DONE.

### Y12 (2026-10-16 ~5:35 pm PT, before any Y12 number) — stations first (Merlin review 6 condition; user approved the download ~6 pm PT relay)
- Fetch (`scripts/fetch_y12.py`, from fetch_2011carriers.py): CN.SYMB 2021-01-01…2026-10-16 and UW.BUCKS 2023-09-17…2026-10-16,
  HH? from EarthScope FDSN dataselect (serves CN too), ≤ 2 parallel, polite retries, default user agent, no personal data;
  decimated in memory to 1 sps (100 sps: 10 × 10), merged per sampling rate, STEIM2 in `data/raw` (a TRACE_ROOT); stop at 3 GB
  stored; log `logs/fetch_y12.log`; marker FETCH_Y12_DONE. StationXML responses for SYMB, BUCKS, RNWY merged into
  `data/resp/{CN,UW}_bb.xml` (originals kept as *.preY12).
- Candidate stations: CN.SYMB, UW.BUCKS (fetched), UW.OSQM, UW.RNWY, UW.SLDQ, UW.ODUC (local). QC per station-year, from the
  prepped 1-sps displacement: daily median |Z| at 20–50 s over CN.PGC's on the same day — pass if the median ratio is within
  0.33–3 and ≤ 10 % of live days are zero/flat; response sanity: peak 20–50 s Z displacement of the year's largest teleseism
  within ×2 of PGC's. Failing station-years are reported, not used. W2b applies as for every station.
- Use: QC-passing stations join the 2021–2026 scan as near verticals (scan v2 station list via `Y12_EXTRA`, same synthetic-GF
  scan templates as every v2 extra station). In the MF they enter only through relay templates (Merlin review 6 step 6); the
  2011 C0–C5/T1–T5 channel lists are unchanged (those windows predate the stations), so MF scores of the 2011 bank are reused.
- Injections: the new scan stations get scan-side sources from the same thrust Green's functions the far references already
  use (`far_synth` precedent): m0(Mw_ref) × scan template at the source node, for C0, C1 and every relay source, in a separate
  file `data/inj/sources_y12.npz` (`INJ_SOURCES` env; the shared sources.npz is untouched). Check before numbers: on the
  existing near stations this recipe matches the empirical C0/C1 Z peaks within ×0.5–2 (cc 0.73–0.92). Disclosed as a method
  choice; flagged for Merlin.
- Pilot 2021 first, as run `prodfy2021` (bank/MF/detect/inj_list reused from prodf2021; battery, inj, summary, R re-run on the
  expanded scan). Pre-Y12 scan files kept (`templates_grid_v2.preY12.npz`, `runs/real_v2_preY12`). Then f vs f′ on the expanded
  network by Merlin's rule; f′ is applied only if f still leaks (C0 f pass < 0.5 of MF-agreeing injections in Jul–Dec 2021).
- 2022–2026: re-run order proposed after the pilot; chain outputs so far are labelled pre-Y12. Marker Y12_DONE.
- Y12 addendum (2026-10-16 ~5:45 pm PT, coordinator): the Green's-function injection sources on the new stations are accepted
  for completeness only, with a sensitivity column: the 2021 pilot reports C0/C1 completeness with the new-station sources at
  ×0.5, ×1 and ×2 (`inj_v2_x0p5`, `inj_v2`, `inj_v2_x2` in data/prodfy2021; only those keys scaled). The 2022–26 prodfy re-run
  starts only if Y12_PILOT2021_DONE shows 2021-10 C0 completeness at Mw 3.8 rising clearly under f (before f′); otherwise the
  coordinator is told first.

### Y13 (2026-10-16 ~8:15 pm PT, before any Y13 number) — do the between-ETS candidates repeat?
- Events: group 1 (candidates) 2016-06-13 07:41, 2020-06-11 13:42, 2020-09-11 05:00, 2020-10-25 19:25 PDT; group 2 (unclear,
  kept separate) 2014-05-05 08:30, 2016-06-13 19:58, 2017-01-31 04:22 PT.
- Templates: 90-s windows at the event's MF origin − 9 s, cut from the event year's MF continuous data (W2b masks), unit norm
  per channel, weights = window RMS / pre-noise RMS (as every bank template). Within a group, events cluster if their
  similarity (mean over common channels of the max normalised CC within ±30 s, the Y9 step-4 metric) is ≥ 0.7 (single
  linkage); clusters of ≥ 2 are stacked per channel (zero lag, unit-norm mean); singletons are used alone. Each template's
  similarity with the 2011 C0 and C1 templates on common channels is reported (≥ 0.7 = same source as the patch VLFEs).
- Detection 2010–2026: the year's MF stations and data (W2b), weighted CC sum on the channels the template has, peaks ≥ 5.5 MAD
  (template-year MAD), ±90 s around the template's own source windows excluded. A single template has no agreement partner,
  so detections are single-template rows; the battery′ ladder (a, b, c, e′, f; unmeasurable = not passed) is computed with
  the year's scan-v2 battery (gated as production: expensive columns only where c and e/e′-time pass, which cannot change
  level-5 counts).
- Controls on the same footing: three noise templates per year (N00–N02 of that year's run, their existing MF outputs) put
  through the same single-template detection and gated battery; control = level-5 detections per noise-template-year, ETS and
  between separately. (Matched banks control agreement and do not apply to single templates.) Clusters also get the
  relay-style LOO (template from the other members detects the held-out one ≥ 5.5 MAD within ±30 s).
- Report per template: level-5 detections ETS vs between (and per live year), Poisson 95 %, noise-template expectation and
  p-value, top-5 between-ETS level-5 detections figured and eye-labelled (Y9 classes), similarity with C0/C1.
- Outcomes per template: (A) between-ETS level-5 count above control (Poisson p < 0.01, counting only detections other than
  the source events) → a between-ETS VLFE source; (B) above control only in ETS → patch-like source firing occasionally
  between ETS; (C) neither → isolated, unconfirmed.
- Queue `scripts/y13_queue.sh`, ≤ 2 workers, starts after Y9B_DONE; years 2022–2026 wait for the prodfy run's cont/qcx steps.
  Marker Y13_DONE.

### Y12b (2026-10-16 ~8:20 pm PT; post-hoc change under the user's standing no-invented-gates rule, disclosed)
- The Y12 station QC's teleseism-ratio limit (×0.5–2 of PGC) is removed as a gate; the ratio stays a QC column. Stations are
  excluded only by the noise-level rule (median daily 20–50 s Z level 0.33–3× PGC) and the flat-day rule (≤ 10 %). Trigger:
  CN.SYMB 2022 failed only on a 2.25 ratio for one M7.6 (2022-09-19, 3,812 km) while its noise level is ×0.99 PGC — plausibly
  site amplification. `cat_y12.qc` changed before the 2023 QC ran, so 2023–2026 use Y12b as the re-run reaches them.
- prodfy2022 is re-run with SYMB (QC, bp_ref, v2 setup/detect, sources, battery, inj, summary, R, relay vetting; MF reused);
  the first run is kept as the no-SYMB variant (`data/prodfy2022_noSYMB`, `tables/noSYMB_2022/`, scan files
  `templates_grid_v2.y12noSYMB.npz`, `runs/real_v2_y12noSYMB`). 2022 is reported both ways. Lane: after Y9B_DONE, ≤ 2 workers,
  before Y13 (`scripts/y12b_2022.sh`, marker Y12B_2022_DONE).

### Y14 (2026-10-16 ~10:58 pm PT, before any Y14 number) — f′ adopted for all years (cites Y11 and Y12)
- Conditions met: Y11 — f′ passed Merlin's adopt rule (every control's level-5 count +≤ 1, largest +0.65; C0 completeness
  +0.38 in the months where f < 0.5); Y12 — with the new stations f still leaks (summer 2022: max station share 0.84–0.90,
  f pass 0.00–0.20 Jul–Oct; summer 2020 likewise). Coordinator decision 2026-10-16 ~10:50 pm PT.
- f′ = (z_loso ≥ 3 AND z_loso ≥ 0.5·z_best) OR z_loso ≥ 5, in `cat_common.apply_rules` and, for injections, `cat_common.rescore_inj`
  (from stored z_loso/z_best). f (pre-Y14) is kept as columns `pass_f_preY14`, `battery_e2_preY14`.
- Applied to ALL years 2010–2026 by re-scoring only: summary, R, per-year summary and ladder re-derived for prod2010–2015,
  prodq2013/2016/2017/2018, prodrq2018, prodf2019/2020, prodfy2021–2023 (`scripts/y14_rescore.sh`, marker Y14_RESCORE_DONE);
  f outputs kept in `tables/preY14/` and `rows_scored_v2_preY14.parquet`. Relay vetting is counted from the re-scored event
  tables; relay banks are NOT rebuilt (the accepted 2018 bank stays as is) and no extended-bank re-run is made — the cost of
  each is listed instead. The running prodfy chain (2024–26) and Y13 pick up f′ in every step started after the change.
