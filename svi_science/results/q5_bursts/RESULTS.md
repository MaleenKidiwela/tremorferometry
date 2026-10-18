# Q5 results — are the deep between-ETS LFE bursts small slow slips?

Run 2026-10-07 ~2:45 pm PT. Rules frozen beforehand in `PLAN.md` (written 2:25 pm PT). Script
`svi_science/scripts/q5_bursts.py` (CPU, ~1 min). All numbers are in `q5_summary.json` and the log is `q5.log`
(log clock is UTC). No Merlin review; provisional.

## Plain answer
**Partly.** Most between-ETS LFE activity is the LFE side of small tremor bursts. Those bursts are what Wech & Creager
(2011) mapped as inter-ETS tremor, and they are usually read as small slow slips. So in that sense most of these bursts
are "mini-ETS". But we have **no direct evidence of slip** for any burst, and the strictly tremor-free deep bursts on
Vancouver Island show no tremor link at all.

**Evidence for:**
- **Tremor within 30 km on the same day.** This uses the broad set (between episodes, any local tremor).
  - SVI: 53 % of bursts vs 4 % on random days (p = 0.001).
  - WA: 33 % vs 4 % (p = 0.001).
  - Counted by events: 89 % (SVI) and 75 % (WA). Most between-ETS LFEs happen during a small local tremor burst.
- **WA's two biggest "quiet" bursts are really small tremor bursts.** These are 2017-11-10 and 2018-11-18 PT, both
  Hood Canal G019 (+G005).
  - They have 143 and 63 PNSN epicentres within 100 km during the burst ± 6 h (112 within ± 24 h for 2018-11).
  - The tremor spans 36–58 km along strike.
  - In 2018-11 the tremor moves along strike over about 40 h (`q5_along_strike_time.png`).
  - The 2018-11-18 burst also has the largest single GNSS offset of any burst: +0.84 mm, z = 2.3. That is suggestive
    only.
- **The same deep patches keep firing.**
  - SVI: G10 (44 km deep) is in 9 of 40 bursts.
  - WA: G006 (42 km) is in 22 of 62 bursts over 11 years. G032 is in 10 bursts. G019 has 94 events in 5 bursts.
  - The same sources lead in 2009–14.

**Evidence against, or missing:**
- **SVI tremor-free bursts have no tremor link at any distance.** They have no tremor within 30 km by definition.
  - Tremor 30–60 km away: 5 % of bursts vs 5 % on random days.
  - Tremor 60–100 km away: 8 % vs 10 %.
  - The biggest SVI burst (2015-04-12 PT, G10, 11 LFEs over 21 h) had no PNSN epicentre within 100 km.
  - This is not a PNSN coverage gap at ETS scale: 46 % of SVI major-ETS family-days have an epicentre within 30 km,
    against 37 % in WA. Tiny tremor may still go undetected there.
- **No migration we can resolve.**
  - Bursts migrate in 0 of 4 eligible SVI bursts and 0 of 6 WA bursts (cross-family pairs; plain Theil–Sen: 1 and 0).
    The bar was > 2.
  - Most bursts are 1–3 families at 1–2 positions. LFEs alone cannot show movement.
- **Short and small.**
  - Big bursts last a median of 8 h (SVI) and 1.1 h (WA). They involve 2–3 families about 10–13 km apart.
  - The smallest catalogued minor episodes last 5–20 days, with tremor spread over 60–110 km.
- **No regular repeat time.**
  - The gaps between bursts are as irregular as random timing: CV 1.2 vs 1.0 for random days in both regions, and the
    same for G10 and G019 alone.
- **GNSS sees nothing, and could not.**
  - The stack bounds each burst at **≤ 0.45 mm, i.e. Mw ≤ 6.0**. Both regions combined give Mw ≤ 5.9.
  - The slow-earthquake scaling predicts Mw 5.0–5.2 for these durations. That would give an offset of only 0.01–0.03 mm,
    15–40× below the bound.

**Bottom line:** the deep between-ETS bursts behave like the small tremor bursts that sit down-dip of the ETS zone. In WA
the biggest ones are such bursts. Whether each one is a slip event cannot be decided with our GNSS (too insensitive by
~1.5 magnitude units) or our family positions (too few). It is consistent with hidden slow slip, but not shown.

## 1. Tremor with the bursts, by distance (`test1_tremor_by_distance.csv`, `figures/q5_tremor_by_distance.png`)
- **Unit:** the burst (events chained by gaps < 24 h). A burst counts as "with tremor" if any of its events has a PNSN
  epicentre in that distance bin on the same UTC day.
- **Control:** random days. Each burst is moved to a random eligible day of its first family, 1000 draws.
- **Pass rule:** above the 99th percentile of the control.

| region | set | bin | bursts with tremor | random days, median (99th pct) | p | events with tremor (random) | verdict |
|---|---|---|---|---|---|---|---|
| SVI | BE (76 ev, 40 bursts) | 30–60 km | 5.0 % | 5.0 % (12.5 %) | 0.50 | 3.9 % (2.6 %) | not tied |
| SVI | BE | 60–100 km | 7.5 % | 10.0 % (22.5 %) | 0.79 | 15.8 % (9.2 %) | not tied |
| SVI | broad (688 ev, 79 bursts) | 0–30 km | **53 %** | 3.8 % (8.9 %) | 0.001 | 89 % (3.2 %) | **tied** |
| WA | BE (220 ev, 62 bursts) | 30–60 km | **13 %** | 3.2 % (8.1 %) | 0.001 | 52 % (3.2 %) | **tied** |
| WA | BE | 60–100 km | **16 %** | 6.5 % (14.5 %) | 0.007 | 28 % (5.9 %) | **tied** |
| WA | broad (934 ev, 84 bursts) | 0–30 km | **33 %** | 3.6 % (9.5 %) | 0.001 | 75 % (2.8 %) | **tied** |
| WA 2009–14 (secondary) | BE (172 ev, 41 bursts) | 30–60 / 60–100 km | 12 % / 12 % | 2.4 % / 4.9 % | 0.014 / 0.054 | 16 % / 15 % | same direction |

- **SVI does not repeat the WA result.** S7's suggestive 17 % vs 9 % (events, 30–100 km) comes from one burst,
  2018-10-21 PT (6 epicentres at 60–100 km). At burst level there is nothing.
- **Where the WA link comes from.**
  - It comes from G019/G005 bursts: 2017-11-10, 2018-11-18, 2020-05-13 and 2020-09-13 PT.
  - Measured from the source's primary position, these bursts even have tremor within 30 km. That is 42 % of WA BE
    events, nearly all from fid 75.
  - The "quiet" label for fid 75 was set at its own detected position, 40 km east of its G019 primary. In plain words:
    these WA "between-episode" bursts sit inside small Hood Canal tremor bursts.
- **±6 h instead of same day:** almost the same shares (SVI 60–100 km: 17 % of events; WA 30–60 km: 52 %).

## 2. Migration inside bursts (`bursts.csv`)
Eligible bursts have ≥ 4 events and families ≥ 5 km apart. Theil–Sen slope of position vs time, compared with event
times shuffled inside the burst (5000 draws). The bar is p < 0.025 per axis, and the region passes if more bursts pass
than Binomial(n, 0.05) allows at 99 %.

| set | eligible | migrating (plain / cross-family) | bar | verdict |
|---|---|---|---|---|
| SVI 2014–26 | 4 of 40 | 1 / 0 | > 2 | no |
| WA 2015–26 | 6 of 62 | 0 / 0 | > 2 | no |
| WA 2009–14 (secondary) | 7 of 41 | 1 / 1 | > 2 | no |

- **SVI 2015-04-12 PT (G10).** The plain slope says up-dip at 0.43 km/h (~10 km/day, ETS-front-like; p = 0.006). The
  cross-family version does not hold (p = 0.11), so it is not robust: three families, 8.6 km apart.
- **SVI 2018-10-21 PT** spans 88 km and 8 families: +3.7 km/h along strike, p = 0.065. Not significant.
- **Tremor column** (PNSN epicentres within 100 km, burst ± 6 h, when there are ≥ 5). WA has 2 of 10 bursts at
  p < 0.025, which does not exceed the bar.
  - 2019-11-08 PT: 5 km/h along strike and down-dip, p < 0.001. That is a rapid-tremor-reversal-like speed.
  - 2010-05-30 PT: 0.6 km/h (~15 km/day, ETS-front-like), p = 0.002.
  - 2017-11-10 PT: 0.53 km/h along strike, p = 0.029.
  - All descriptive.
- **Deviation:** the cross-family columns were added after the first run. One family often holds most events, and
  same-family pairs have zero position change, which forces the plain slope to 0. The same pass rule was applied to both.

## 3. Duration and size (`bursts.csv`, `minor_episodes.csv`)
| | n | duration | families | max family separation | PNSN tremor within 100 km |
|---|---|---|---|---|---|
| SVI big bursts (≥ 3 events) | 7 | median 8.3 h (0.0–48 h) | median 3 | median 13 km (max 88) | 2 of 7 bursts, 1–6 epicentres |
| WA big bursts | 14 | median 1.1 h (0.1–25 h) | median 2 | median 10 km (max 81) | 5 of 14, 9–143 epicentres, 24–58 km span |
| WA 2009–14 big bursts | 13 | median 3.0 h | median 2 | median 11 km | 4 of 13, 56–138 epicentres |
| SVI minor episodes | 40 | median 6 d (5–20 d) | 55 % have no P ≥ 0.75 LFE | — | median 180 epicentres, 62 km span |
| WA minor episodes | 40 | median 7 d (5–13 d) | 68 % have no LFE | — | median 177 epicentres, 114 km span |

- The bursts are 10–100× shorter than the smallest catalogued episodes, and smaller in extent.
- They sit below the episode-catalogue floor (≥ 5 days). They fill the gap between single LFEs and minor episodes.
- This matches Wech & Creager's picture of a continuum of ever smaller, more frequent tremor bursts down-dip.
- One oddity: most minor (tremor-defined) episodes have no confident LFE at all. The bursts and the minor episodes are
  largely different populations.

## 4. Same patches, regular repeats? (`q5_summary.json` → test4)
- **The same deep patches, yes.**
  - SVI: G10 (44 km) is in 9 of 40 bursts over 4 years, and G1 (42 km) in 8 bursts over 6 years.
  - WA 2015–26: G006 (42 km) is in 22 of 62 bursts over 11 different years, G032 (45 km) in 10, and G019 (46 km) holds
    94 events in only 5 bursts.
  - WA 2009–14: the same three (G006 14 of 41, G032 8, G019 3).
  - These are persistent deep sources, active for a decade.
- **Regular repeat interval: no.** The control moves each burst to a random eligible day of its first family.

  | set | n | CV of intervals | random days, median (1–99 %) | median interval | verdict |
  |---|---|---|---|---|---|
  | SVI all bursts | 40 | 1.20 | 0.96 (0.71–1.41) | 74 d | random |
  | SVI G10 bursts | 9 | 1.23 | 0.83 (0.41–1.50) | 116 d | random |
  | SVI big bursts | 7 | 0.35 | 0.78 (0.30–1.57) | 680 d | random (p_regular = 0.017, misses 0.01) |
  | WA all bursts | 62 | 1.17 | 1.02 (0.81–1.35) | 45 d | random |
  | WA G019 bursts | 5 | 0.50 | 0.70 (0.18–1.44) | 377 d | random |
  | WA 2009–14 all | 41 | 0.81 | 1.01 (0.74–1.45) | 46 d | random |

## 5. GNSS on the biggest bursts (`figures/q5_gnss_bursts.png`)
- **Method:** the S7/A estimator, a 5-day median after minus before, projected on each station's own ETS direction and
  averaged over SVI-14 / WA-12. The control is season-matched random days.
- **Excluded bursts:** those within 8 days of an episode (as S7/A): SVI 2015-04-12; WA 5 of 14.
- **From mm to Mw:**
  - Pure thrust slip on the Slab2 mesh patches within 10 km of each burst. The patches sit at 33–54 km depth.
  - Cutde half-space Green's functions, μ = 40 GPa.
  - This gives 0.031 mm (SVI) and 0.034 mm (WA) of mean station offset per 10¹⁷ N·m. The value barely changes for a
    5 or 20 km patch.

| set | bursts | offset (mm) | random days mean ± σ | z | bound (mm) | bound Mw |
|---|---|---|---|---|---|---|
| SVI | 6 | −0.37 | −0.08 ± 0.15 | −1.9 | 0.45 | **≤ 6.0** |
| WA | 9 | −0.02 | −0.08 ± 0.13 | +0.4 | 0.46 | **≤ 6.0** |
| WA 2009–26 (secondary) | 20 | −0.10 | −0.07 ± 0.08 | −0.4 | 0.24 | ≤ 5.8 |
| both regions combined | 15 | — | — | −0.9 | — | **≤ 5.9** |
| single: SVI 2018-10-21 PT | 1 | −0.47 | σ 0.36 | −0.9 | 1.08 | ≤ 6.3 |
| single: WA 2018-11-18 PT | 1 | +0.84 | σ 0.40 | **+2.3** | 2.12 | ≤ 6.6 |
| single: WA 2017-11-10 PT | 1 | −0.01 | σ 0.41 | +0.1 | 1.27 | ≤ 6.3 |

- **What a slow slip of this duration would be.** The slow-earthquake scaling is M0 ≈ T × 10^12.5±0.5 N·m/s
  (Ide et al. 2007).
  - The median big burst would be Mw 5.2 (4.9–5.6) in SVI and Mw 5.0 (4.6–5.3) in WA.
  - That predicts an offset of 0.03 mm (SVI) and 0.01 mm (WA), 15–40× below our bound.
  - With a single-burst σ of ~0.4 mm, a 3σ detection would need roughly 1,000–10,000 bursts.
- **The stack cannot test the slow-slip idea.** A null here is expected even if every burst is a slip event.
  - Weighting stations by their predicted response would gain about 2–3×, about 0.3 magnitude units. That is not
    enough.
- **2018-11-18 PT.** This WA burst gives +0.84 mm (z = 2.3), in the slip direction, with 63 tremor epicentres within 100 km (± 6 h; 112 within ± 24 h) and a
  migrating tremor front. It is the best single candidate, but it is below the bar. The M 4.1 earthquake 74 km away at
  03:08 PT is far too small to move GNSS by that much.

## Comparison with earlier work
- **Wech & Creager (2011, northern Cascadia).**
  - They found tremor between ETS in many small swarms, more frequent and smaller down-dip of the ETS zone, and read it
    as a depth continuum of slip.
  - Our between-ETS LFEs fit that picture. They are deeper than the ETS LFEs (S7/A), and most sit inside small local
    tremor bursts (53 % / 33 % of bursts, 89 % / 75 % of events).
  - Bursts recur on the same deep patches for over a decade.
  - New here: LFE family identity shows that one or two deep sources carry most of it (G10; G006/G032/G019). We also
    find tremor-free deep bursts on SVI that have no catalogued tremor anywhere nearby.
- **Frank et al. (2015, GRL; Frank 2016, GRL; Guerrero and Cascadia).**
  - They showed that LFE or tremor bursts between large slow slip events carry hidden slow slip. Stacking GPS on very
    many burst times brought out a sub-millimetre slip signal.
  - We **cannot** reproduce that here. We have 6–20 bursts, not hundreds, and our burst list is restricted to days with
    no catalogued episode.
  - Our bound (Mw ≤ 5.9–6.0 per burst) is far above what their result implies for burst-sized slips.
  - Their approach needs every LFE-active day, not only the "quiet" ones. With our 2014–26 catalogue that is the
    obvious next test: stack GNSS on all between-episode LFE days, ranked by LFE count, as Frank (2016) did with
    tremor.
- **Honest difference.** In Guerrero the inter-SSE LFE bursts recur quasi-periodically, and many of them do. Ours are
  random in time (CV ~1.2) and few. We see neither the regularity nor the geodetic signal. That is not evidence that
  slip is absent.

## Failed / not passed (as pre-registered)
- Test 1, SVI BE at 30–60 and 60–100 km: **not tied**. The WA result does not repeat in SVI.
- Test 2, migration: **not passed** in any set.
- Test 4, regular recurrence: **not passed**. All sets are consistent with random timing.
- Test 5, GNSS: **no slip-like signal; bound only**, and the bound is above the expected size.
- Passed: test 1 in WA (both BE bins) and the 0–30 km broad-set link in **both** regions.

## Deviations and caveats
- The cross-family Theil–Sen columns were added after the first run (see §2). No rule was changed.
- **WA positions:** family copies use the source's primary position (as in A). Their "quiet day" label was set at the
  copy's own position. For fid 75 these differ by 40 km, so WA BE bursts include days with tremor within 30 km of the
  primary position.
- The first run died on a day-index unit bug (it was killed at 120 GB RSS) and then on a date-parse error. Both are
  fixed. The numbers above are from the clean run.
- The PNSN catalogue ends 2026-05-27. Bursts after that are not used.

## Files
- `PLAN.md`, `q5_summary.json`, `q5.log`, `run.out`
- `test1_tremor_by_distance.csv`, `events_{SVI,WA}_tremor_bins.csv`
- `bursts.csv` (every burst: PT times, families, duration, separation, tremor counts/span, migration slopes and p)
- `minor_episodes.csv`, `gnss_ctrl_{SVI,WA}.npy`
- Figures:
  - `svi_science/figures/q5_map_bursts.png`: map of the 4 biggest bursts per region, with tremor ≤ 100 km, ± 12 h.
  - `q5_along_strike_time.png`: along-strike position vs time for the 2 biggest bursts per region, with tremor.
  - `q5_tremor_by_distance.png`: burst–tremor coincidence by distance against random days, both regions.
  - `q5_gnss_bursts.png`: GNSS stack against the random-day distribution.
