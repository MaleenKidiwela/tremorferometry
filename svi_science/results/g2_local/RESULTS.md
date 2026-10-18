# G2 — Local earthquakes and deep LFEs under SVI and Washington, both directions. RESULTS

Written Sun 2026-10-11, ~7:10 am PT. Rules are in `PLAN.md`, written before any number was computed, with Revisions
R1–R5 at its end (R4 = post-hoc checks, labelled below). Scripts `svi_science/scripts/g2_local_{dl,build,tests,
posthoc,eye,figs}.py`. Tables in `tables/`, figures in `figures/`. "Scrambled control" = the same detector on the
same day, with each station site's data shifted to another time of that day. Times are Pacific.

## In plain words

- **Local earthquakes do not set off deep LFEs.** We looked at ~4,000 local earthquakes per region (M ≥ 1.5, within
  150 km of the LFE families, 2003–2026 for SVI, 2009–2026 for Washington). In the 1, 6 and 24 hours after them, and on
  days 1–3 and 1–10 after, the LFE count is what random dates give. This holds for crustal, intraslab and
  near-interface earthquakes, for big and small ones, for near and far ones, between ETS and during ETS. No headline
  test came near the p < 0.01 bar (smallest headline p 0.10).
- **The largest local shaking did nothing.** The 2015-12-29 M4.7 intraslab earthquake under Sidney / Victoria
  (11:39 pm PST, 58 km deep, 14 km from an SVI family, ~120 kPa estimated dynamic stress) struck during the
  2015-12-19 … 2016-01-11 big SVI ETS: 30 LFEs in the next 6 h vs 30.4 in the same-length windows around it. Figure `g2_eye_SVI_2559.png`.
- **Deep slow activity does not change the local earthquake rate.** During ETS episodes, on days with many
  multi-site LFEs, on PNSN-tremor days, and during the small between-ETS deep LFE bursts, the declustered rate of
  M ≥ 1.5 crustal, intraslab and near-interface earthquakes is the same as at other times (rate ratios 0.80–1.10,
  none above the shuffled-episode-date control at p < 0.01). Tremor does not hide small earthquakes at M ≥ 1.5:
  completeness is the same during ETS (Mc 1.4 vs 1.3 SVI, 1.2 vs 1.2 WA).
- **Where:** nothing is concentrated at the up-dip edge, the down-dip edge, or in the slab right under the
  families. No location bin passes.
- **The trap again, smaller.** In the first hour after a local earthquake the raw matched-filter candidates rise 18 %
  at SVI (17,373 vs 14,700; p 0.0001) while the SVI scrambled control does not (8,949 vs 8,815). The P rules remove
  it: LFE counts in that hour are 38 vs 44 expected. Same lesson as study G: the scrambled control alone would not
  have caught it at SVI.
- **Two named cases worth remembering (not proofs).** Both are moderate deep earthquakes that happened *inside* an
  ongoing between-ETS deep LFE burst right next to them (figure `g2_cases_bursts.png`):
  - **Washington, Mon 2018-11-19 3:09 am PST, M4.1 intraslab, 39 km deep, 10 km from a family** (38 km NE of Amanda
    Park). The deep burst (with PNSN tremor) had started ~31 h earlier and peaked in the 3 h before the
    earthquake. LFEs paused for ~85 min after it, then went on for another ~15 h. Borehole B003/B004 tremor bursts
    ran straight through. Figure `g2_eye_WA_2227.png`.
  - **SVI, Sun 2006-01-15 4:29 am PST, M3.6 near-interface, 44 km deep, 7 km from a family** ("near Victoria"). It
    fell inside the deep Jan 14–16 2006 burst (R23); the burst started ~10 h before and was strongest 7–12 h after.
  - Post-hoc chance: of 13 such earthquakes (M ≥ 3.5, intraslab or near-interface, ≤ 15 km from a family, between
    ETS) 2 fell in a burst window vs 0.21 expected (p 0.018). Below the bar, post-hoc, and the timing runs both ways
    (burst first in both). Reading: a deep slow burst and a nearby slab earthquake can coincide; we cannot say one
    caused the other.
- **Upper bounds (99 %).**
  - Q1 (earthquakes → LFEs), between ETS, 6 h after: ≤ 0.03 (SVI) / ≤ 0.08 (WA) extra P-weighted LFEs per earthquake,
    i.e. ≤ +43 % / ≤ +80 % of normal; per class SVI crustal ≤ +87 %, intraslab ≤ +40 %, near-interface ≤ +70 %; WA
    crustal ≤ +117 %, intraslab ≤ +209 % (held up by the 2018-11-19 case), near-interface ≤ +32 %.
  - During ETS, 6 h after vs the ±24 h around: ≤ +8.5 % (SVI), ≤ +18 % (WA).
  - Days 1–3 after, between ETS: ≤ +22 % (SVI) / ≤ +26 % (WA) of the background rate.
  - Per kPa of estimated dynamic stress: ≤ +20 % (SVI) / ≤ +16 % (WA) extra LFEs.
  - Q2 (deep slow → earthquakes), rate during ETS / outside: SVI ≤ 1.06 (all), crustal ≤ 1.10, intraslab ≤ 1.21,
    near-interface ≤ 1.06; WA ≤ 1.05, crustal ≤ 1.07, intraslab ≤ 1.23, near-interface ≤ 1.42. During between-ETS
    bursts: SVI ≤ 1.14, WA ≤ 1.20 (all classes).
- **Compared with the literature (from memory, not checked here).**
  - Cascadia: I recall that ETS episodes in northern Cascadia do not obviously change regional seismicity at catalogue
    level, and one array study (Vidale et al. 2011, GRL, "tiny intraplate earthquakes triggered by nearby ETS") reported
    a small number of very small slab earthquakes near the tremor during an episode. Our M ≥ 1.5 network catalogue
    cannot see such tiny events; at M ≥ 1.5 we find no change, within ~5–20 % for all classes together.
  - Mexico (Guerrero): large slow slip events have been linked to intraslab / interface seismicity, and the 2014
    Papanoa M7.3 struck during an ongoing slow slip event (Radiguet et al. 2016, Nature Geoscience). There the slow
    slip is shallower, larger and adjacent to the seismogenic zone.
  - New Zealand (Hikurangi): shallow slow slip events off Gisborne come with swarms of small earthquakes
    (Delahaye et al. 2009 and later work), and the 2016 Kaikōura M7.8 triggered widespread slow slip. Those slow slip
    events are at < 15 km depth, next to brittle rock.
  - Our deep (30–45 km) ETS sits well below the locked zone and below most crustal seismicity, so a null link at
    catalogue level is what one would expect; the two burst-coincident slab earthquakes are the only hint, and are
    at the chance level.

## 0. Data

- **LFEs**: study G's tables (`g_trigger/tables/`), i.e. SVI 2003-01 … 2026-09 (backward 2003–13 under the
  RESULTS_BACKWARD usage rules + forward 2014+, counted once, P rules), WA-NW 2007-01 … 2026-09 (P ≥ 0.75, counted
  once); scrambled control rows through the same steps; n_z2 multi-site flag; episode lists SVI 145 (pre-2010
  envelope-flag minors included, unvalidated) and WA-NW 118; PNSN tremor. Headline counts are Σ P (P-weighted);
  plain counts are columns and agree.
- **Earthquakes** (`tables/eq_local.parquet`, `dl/`): ComCat 33,540 rows (26,976 earthquakes, 6,560 explosions …)
  + NRCan 12,655 rows (8,264 earthquakes, 4,390 mining explosions, 1 induced); 4,988 cross-catalogue duplicates;
  9,942 non-earthquakes removed (57 more vetoed because the other catalogue called them blasts); volcano circles
  removed. Neither catalogue holds LFE / tremor types. 18,363 events within 150 km of the families (M ≥ 1.0),
  13,612 Gardner–Knopoff mainshocks.
- **Classes** (Slab2; McCrory 2012 as a column): of the analysis mainshocks within 150 km, crustal 3,155,
  overriding-deep (≥ 20 km, > 5 km above the slab) 1,347, near-interface (± 5 km of the slab top) 507, intraslab 609
  (map `figures/g2_map.png`).
- **Completeness** (`tables/completeness.csv`): maximum curvature + 0.2 per 3-year block gives 1.2–1.5 at SVI and
  1.2–1.3 at WA → Mc_R = 1.5 both (the floor). Trials: SVI 4,748 (between 3,435, ETS 995, pad 318), WA 4,035
  (2,677 / 1,043 / 315), WA 2007–09 687.

## 1. Q1 — do local earthquakes trigger LFEs?

### Hours (`tables/q1_hours_stacks.csv`, `q1_trials_*.parquet`; figure `g2_superposed.png`)

Headline = 6-h window after origin, between ETS. Σ P observed vs random-date mean; p_R (data), p_S (data minus
scrambled), p_mz (multi-site); 99 % bound on extra Σ P per earthquake.

| region / class | n | Σ P obs / random | p_R | p_S | p_mz | after − before p | bound per EQ (relative) |
|---|---|---|---|---|---|---|---|
| SVI all | 3,435 | 211.5 / 235.9 | 0.60 | 0.62 | 0.60 | 0.20 | 0.030 (+43 %) |
| SVI crustal | 1,641 | 138.3 / 116.5 | 0.27 | 0.27 | 0.27 | 0.13 | 0.062 (+87 %) |
| SVI intraslab | 405 | 12.1 / 26.2 | 0.71 | 0.73 | 0.71 | 0.28 | 0.026 (+40 %) |
| SVI near-interface | 373 | 17.1 / 23.3 | 0.53 | 0.57 | 0.53 | 0.73 | 0.044 (+70 %) |
| WA all | 2,677 | 329.5 / 267.8 | 0.22 | 0.22 | 0.21 | 0.28 | 0.080 (+80 %) |
| WA crustal | 1,414 | 205.7 / 141.5 | 0.14 | 0.15 | 0.14 | 0.021 | 0.117 (+117 %) |
| WA intraslab | 278 | 61.8 / 29.3 | 0.097 | 0.10 | 0.097 | 0.98 | 0.220 (+209 %) |
| WA near-interface | 254 | 8.3 / 25.8 | 0.72 | 0.72 | 0.72 | 0.74 | 0.033 (+32 %) |
| WA 2007–09 all (column) | 611 | 118.5 / 411.2 | 1.0 | 1.0 | 1.0 | 0.52 | — |

- **All FAIL.** The 1-h and 24-h windows give the same (1 h: SVI 33.8 vs 39.4, p 0.58; WA 33.9 vs 44.1, p 0.63).
- **Clustering check** (segment-wise circular shift, keeps aftershock-like clustering, `q1_shift_check.csv`): p 0.15–0.63
  for every headline class.
- **Columns** (magnitude, distance, dynamic / static stress, McCrory classes, segments, location bins; ~120 stacks):
  one stack passes the three legs, **WA M ≥ 3.5** (31 earthquakes, Σ P 43.7 vs 2.5, p 0.007). It is the single
  2018-11-19 M4.1 trial; its "before" window holds more LFEs than its "after" window (after − before p 0.998) and
  the raw-candidate gauge is also up (p 0.008). Eye check: the earthquake fell inside an ongoing deep burst (§4), so
  it is not a triggered response → **vetoed**. Next smallest: SVI "beneath, intraslab" (12 earthquakes, 8.0 vs 0.6,
  p 0.021; one M1.6 on 2013-04-21) and WA "beneath, intraslab" (84, p 0.020; the same 2018-11-19 trial).
- **ETS column** (6 h after vs the ±24 h around): SVI 4,634 vs 4,946 (p 0.84; bound ≤ +8.5 %), WA 2,993 vs 2,939
  (p 0.40; ≤ +18 %). 1 h: SVI 807 vs 817 (p 0.53); WA 596 vs 493 (p 0.047; raw candidates also up, p 0.006 = the
  earthquake's own energy at the very start).
- **Artefact gauge (the trap).** First hour, between ETS: SVI raw candidates 17,373 vs 14,700 (p 0.0001), scrambled
  8,949 vs 8,815 (p 0.37); WA raw 166,758 vs 145,316 and scrambled 162,408 vs 144,176 (both p 0.0001; the WA scrambled
  series has far more candidates). After 6 h the gauge is back to normal (SVI p 0.29, WA p 0.08). P removes it.
- **Dose-response** (`q1_dose.csv`): Spearman ρ −0.001 (SVI, p 0.54), +0.027 (WA, p 0.089) against estimated
  dynamic stress. β ≤ 0.20 / kPa (SVI), ≤ 0.16 / kPa (WA). WA's β lower bound is above 0 only because of the
  2018-11-19 trial (75 kPa); ρ fails, so "does not scale".
- **Named earthquakes (rule a–d).** 28 at SVI (27 during ETS) and 19 at WA (18 during ETS). **Post-hoc** (R4,
  `posthoc_named_rate.csv`): random pseudo-earthquakes at the same status meet the rule at the same rate —
  expected 20.4 + 2.7 (SVI) and 17.3 + 2.1 (WA); p 0.09 / 0.93 / 0.46 / 0.88. During ETS the ±24 h background has
  only 8 tiles, so any ETS burst that switches on near an earthquake qualifies. The ones drawn (§4) are bursts that
  were already running, or that started 2–3 h later at negligible stress.

### Days (`tables/q1_delayed_stacks.csv`; G's test 2, Σ P per observable family-day, data − scrambled)

| region / class | n | days +1…+3 obs / random | p (after) | p (after − before) | p (multi-site) | bound |
|---|---|---|---|---|---|---|
| SVI all | 2,849 | 0.0022 / 0.0024 | 0.73 | 0.79 | 0.72 | ≤ +22 % |
| SVI crustal | 1,372 | 0.0018 / 0.0024 | 0.91 | 0.96 | 0.91 | ≤ +17 % |
| SVI intraslab | 321 | 0.0017 / 0.0023 | 0.66 | 0.58 | 0.66 | ≤ +53 % |
| SVI near-interface | 306 | 0.0033 / 0.0022 | 0.15 | 0.10 | 0.15 | ≤ +121 % |
| WA all | 2,180 | 0.0023 / 0.0027 | 0.75 | 0.92 | 0.75 | ≤ +26 % |
| WA crustal | 1,148 | 0.0022 / 0.0027 | 0.76 | 0.96 | 0.75 | ≤ +33 % |
| WA intraslab | 221 | 0.0021 / 0.0027 | 0.58 | 0.47 | 0.58 | ≤ +65 % |
| WA near-interface | 204 | 0.0035 / 0.0027 | 0.27 | 0.17 | 0.27 | ≤ +119 % |

All FAIL. Days +1…+10 and counts (not P-weighted) the same. Smallest column: SVI within 30 km, +1…+10, p 0.025
(after − before p 0.004), with the scrambled-control-alone leg also at p 0.022 → fails the rule.

## 2. Q2 — does deep slow activity trigger local earthquakes?

`tables/q2_rates.csv`, `q2_detection.csv`, `q2_bursts_*.csv`; figure `g2_q2_rates.png`. RR = rate of declustered
M ≥ 1.5 mainshocks in the state / rate outside (outside every episode ± 2 d). Exact 99 % range from the
conditional binomial; p from the two controls (circular shift of the episode list; random placement).

| region | state | all | crustal | intraslab | near-interface | overriding-deep |
|---|---|---|---|---|---|---|
| SVI | ETS days (1,755 d) | 0.97 (0.89–1.06) p 0.79/0.79 | 0.98 (0.86–1.10) | 0.94 (0.72–1.21) | 0.80 (0.60–1.06) | 1.10 (0.93–1.29) p 0.08 |
| SVI | major ETS only | 0.98 | 0.97 | 0.92 | 0.85 | 1.13 p 0.07 |
| SVI | LFE-active days (555) | 0.98 (0.85–1.12) | 0.93 | 1.02 | 0.99 | 1.11 |
| SVI | between-ETS bursts (65 bursts, 180 d) | 0.89 (0.69–1.14) | 1.04 | **0.16** (low, p 0.001) | 0.91 | 0.94 |
| WA | ETS days (1,572 d) | 0.96 (0.88–1.05) p 0.81/0.81 | 0.95 (0.84–1.07) | 0.94 (0.71–1.23) | 1.08 (0.82–1.42) | 0.95 |
| WA | major ETS only | 0.96 | 0.92 | 1.08 | 1.11 | 0.95 |
| WA | LFE-active days (430) | 0.98 (0.85–1.14) | 0.89 | 0.84 | 1.21 | 1.10 |
| WA | between-ETS bursts (47 bursts, 114 d) | 0.90 (0.67–1.20) | 0.88 | 0.64 | 1.12 | 0.95 |

- **No pass.** Every headline p ≥ 0.07 against both controls. The multi-site LFE leg (LFE-active days) agrees.
- **PNSN-tremor days** (≥ 10 epicentres ≤ 100 km, 2010+): RR 0.99 (SVI), 1.00 (WA).
- **10 days after / before an episode**: all within 0.8–1.25 except one column, SVI near-interface in the 10 days
  *before* a major ETS, RR 1.56 (46 events, shift-control p 0.004). One column of ~600, placement control not run,
  and the same class is *lower* during ETS (0.80) → a lead only.
- **Lower, not higher.** SVI intraslab earthquakes were rarer during between-ETS bursts (2 vs ~12 expected; low-tail
  p 0.001). Opposite to the question; not seen in WA (5 vs 7.9). Not interpreted.
- **Detection check (veto, passed):** Mc during ETS 1.4 (90 % 1.2–1.5) vs 1.3 (1.2–1.4) at SVI, 1.2 vs 1.2 at WA;
  masking index (M 1.0–1.5 per M ≥ 1.5) 1.42 vs 1.39 (SVI), 1.57 vs 1.53 (WA) → tremor does not hide M ≥ 1.5 events.
- **Post-hoc** (R4, `posthoc_burst_after_eq.csv`): burst starts with a local earthquake in the 6 h before: SVI 18 %
  vs 13 % random (p 0.15), WA 13 % vs 16 % (p 0.77).

## 3. Q3 — where?

- **Q1 by position** (between ETS, 6 h): up-dip, LFE-zone, down-dip and outside-strike bins all p ≥ 0.05, except
  single-case bins: "beneath the families, intraslab" (SVI p 0.021, WA p 0.020) and WA LFE-zone intraslab (p 0.017)
  are the 2018-11-19 / 2013-04-21 trials; WA down-dip crustal (634 earthquakes, p 0.042) is the 2017-11-10 Tanner
  M1.7 trial (G019 burst, tremor already running). Days +1…+10: SVI LFE-zone p 0.051, the rest
  ≥ 0.09.
- **Q2 by position**: no bin passes. The smallest p values are low-tail (fewer earthquakes during ETS or bursts):
  SVI up-dip during ETS RR 0.78 (low p 0.02), WA down-dip in the 10 days before ETS RR 0.79 (low p 0.001,
  shift only). "Beneath the families, intraslab" is small (SVI 16 events, WA 112): RR 1.18 (SVI), 0.83 (WA).
- **Answer:** no link is concentrated anywhere. The only spatial hint is the post-hoc pair of slab earthquakes
  within 10 km of a family inside a deep burst (§4).

## 4. Eye check (2–8 Hz record sections, `figures/g2_eye_<region>_<trial>.png`)

| figure | earthquake (PT) | status | LFEs in 6 h (local / random mean) | verdict |
|---|---|---|---|---|
| WA_2227 | 2018-11-19 3:09 am M4.1 intraslab, 10 km | between | 50 (0.05) | earthquake inside an ongoing deep burst (LFEs and B003/B004 tremor from ≥ 1 h before; LFE pause +5…+88 min, then resume) → not a trigger; veto |
| SVI_448 + `g2_cases_bursts.png` | 2006-01-15 4:29 am M3.6 near-interface, 7 km | between | 3 (in 6 h; 68 in 24 h) | inside the Jan 14–16 2006 deep burst, which started ~10 h before and peaked 7–12 h after |
| WA_2002 | 2017-11-10 9:54 pm M1.7 crustal, 44 km | between | 153 (0.07) | PNSN tremor already running 30–55 min before origin; ~0.02 kPa → coincidence with the G019 burst |
| SVI_491 | 2006-05-25 4:46 pm M1.9 crustal, 80 km | between | 39 (0.29) | LFEs from 40 min before → ongoing burst |
| SVI_1306 | 2010-01-22 1:09 pm M1.8 crustal, 102 km | between | 16 (0.02) | LFEs start 150 min later, stress 5 × 10⁻⁵ kPa → coincidence |
| WA_2339 | 2019-06-15 1:52 pm M1.8 near-interface, 13 km | ETS (late) | 47 (1.9) | LFEs start ~3 h later, ~1 h after the Kermadec M7.3 surface waves (study G's candidate); neither is a clear trigger |
| WA_4015 | 2026-04-25 11:37 pm M1.5 deep crustal, 12 km | ETS | 88 (0.25) | patch switches on 17–28 min later during ETS (tremor first); scrambled control also 4 → ETS onset coincidence |
| SVI_2559 | 2015-12-29 11:39 pm M4.7 intraslab, 14 km, ~120 kPa | ETS | 30 (30.4) | earthquake coda only; ETS rate unchanged |
| SVI_41, SVI_1060, SVI_3339, WA_2348, WA_3305, WA_3659 | largest between-ETS local events (2003-04-25 M4.8, 2009-01-30 M4.7, 2020-01-24 M4.5, 2019-07-12 M4.6, 2023-10-08 M4.3, 2025-03-03 M4.5) | between | 0 each | nothing |
| SVI_502, 2806, 2045, 2132, 905, 262, 45, 401; WA_1034, 128, 1035, 283, 1581, 1053, 549, 1565 | remaining §5 list | — | — | drawn, not individually reviewed (no test outcome depends on them; R5) |

## 5. Caveats

- **Catalogue mixing.** NRCan ML (Canadian side) and PNSN Md / Ml (US side) are merged at a crude border; Mc = 1.5 is
  the floor of the rule and maximum-curvature estimates tend to run low. The rate tests compare states within one
  merged catalogue, so a stationary mixture does not bias them. Pre-2010 NRCan depths are often fixed → those events
  are "unknown depth" and only in "all".
- **Stress values are order-of-magnitude** (point-source static, Brune far-field dynamic with Q = 300). They rank
  events; they are not calibrated to measured strain.
- **Families only.** LFEs are counted on our template families; triggered tremor elsewhere is not seen (PNSN columns
  partly cover it from 2010).
- **Between-episode rates** are a control-subtracted bracket (RESULTS_BACKWARD rules); every Q1 comparison is within
  one network segment / era.
- **WA 2007–09** knows only the 3 WECC episode windows, so its "between" random dates include unlisted ETS days
  (random mean 411 vs 118 observed). Column only.
- **Many columns.** ~120 Q1 stacks and ~600 Q2 rows; a few p ≈ 0.01–0.05 are expected by chance. Only the 16 Q1 and
  ~24 Q2 headline rows carry the pass rule.
