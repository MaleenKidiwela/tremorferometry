# Step 1 — north positive control of the array-discovery engine (v2: 1a, 1b, 1c + diagnostics D0–D7)

> **This file replaces the 2026-09-20 15:46Z Step-1a write-up**, which is kept verbatim as
> Appendix A (also at `results/STEP1_RESULTS_v1_superseded.md`).  Sections 0 and 1 record the
> FROZEN runs and their pre-registered consequence; section 2 is **DIAGNOSTIC** and changes
> no verdict, retunes no constant and proposes no new bar.  Written by agent G after the
> three runs, 2026-09-20.

## 0. Frozen verdicts

All three runs were frozen before launch (`results/freeze_step1{,b,c}.json`: sha256 of every
`w_` script and of the imported `lfe_picker` scripts).  **No engine constant was retuned in
or after any of them.**  Scoring is `w_score_step1.py` against the 10 *isolated* hidden
families of `lfe_picker/results/holdouts.json`; accept bars, pre-registered: ≥3/10
rediscovered, spurious-split ≤0.10, site-shift null ≤1 gated family per realisation.

| | **1a** permanent bed | **1b** + dense CN/C8 | **1c** 2005 POLARIS bed |
|---|---|---|---|
| window | 2012-08-15..10-05 | 2012-08-15..10-05 | 2005-08-29..09-21 |
| bed | `lfe_catalog_repair/cache20hz` | + `lfe_array_discovery/cache20hz` | `lfe_array_discovery/cache20hz` |
| detection stations / sites | 12 / 6 | 16 / 5 | 12 / 12 |
| site distance to the families | 46–91 km (median hypocentral) | 28–91 km (LZB 28 km) | **6 sites at 25.7–35.6 km** median epicentral (GLBC 25.7, SHVB 26.3, TWKB 29.5, MGCB 30.8, TSJB 30.9, TWBB 35.6) |
| withheld | VGZ SNB B926 | VGZ SNB B926 YOUB | VGZ SNB LZB JRBC |
| station-days | 618 | 740 | 880 |
| quiet hours → **T** | 45 → **7.25** (FA rule) | 45 → **7.25** (FA rule) | 40 → **12.50** (FA rule) |
| tremor hours | 563 (PNSN epicentre rule) | 563 | 299 (array-energy flag) |
| pairs recorded | 199,579 (≥5.0) | 163,973 (≥5.0) | 143,524 (≥5.0) |
| pairs ≥ T | 5,089 | 3,993 | **25** |
| connected components | 126 | 182 | 4 |
| **seed clusters** | 3 | 2 | **0** |
| pass-0 templates ≥4 sites | 3 | 2 | 0 |
| survivors after growth + scramble null | 0 | 0 | 0 |
| deduped / located / gated | 0 | 0 | 0 |
| **isolated hidden rediscovered** | **0/10** | **0/10** | **0/10** |
| all-30 hidden rediscovered | 0/30 | 0/30 | 0/30 |
| spurious-split rate | 0.00 | 0.00 | 0.00 |
| site-shift null, gated families ×3 | not run | 0 / 0 / 0 | 0 / 0 / 0 |
| **verdict** | `ENGINE_INADEQUATE` | `ENGINE_INADEQUATE` | `ENGINE_INADEQUATE` |

Ceiling for the same scorer: the frozen *template* detector — which is given the answers —
scores 5/10 isolated and 9/30 absolute.  So 0/10 is a floor result, not a near miss.

**Pre-registered consequence (PLAN §10 / Merlin Step 1): `ENGINE_INADEQUATE` → *fix on the
north only; WA is not touched*.**  Accordingly `results/WA_GO` has NOT been created, and the
WA Step-5 discovery queue (`scripts/x_step5_queue.sh`, pid 1122190) remains parked at its GPU
gate.  Merlin's stop rule — "Step 1 failing at <2/10 after one repair cycle" — has not yet
been reached: one repair cycle is still allowed, and these diagnostics precede it.

## 1. Beds and deviations

**1a — the north permanent bed.**  The 2012 south-Vancouver-Island cache on disk: PB
boreholes + PGC/NLLB/LRIV, 6 sites at 46–91 km median hypocentral distance from the 130
Bostock families.  Deviation on record: 45 quiet *hours* were taken as 5 hours on each of 9
quiet days, because the 52-day bed has only 9 days ≥3 days from a PNSN tremor day — Merlin's
binding clause (≥40 quiet hours on the same sites) is met, the ≥20 quiet *days* is not.

**1b — the same window with the dense stations that exist in it.**  POLARIS-BC ran 2002–2007,
so 2012 offers only CN.LZB (28 km) plus the C8 PA02–05/PA12 cluster and BPCB.  Sites were
added by increasing median epicentral distance while the median S-moveout spread over the
families stayed ≤ `AC_WIN_S` = 10 s (9.7 s achieved), which dropped the far PB boreholes.
Constants unchanged; YOUB added to the withheld set.

**1c — the 2005-09 ETS on the POLARIS array: the first bed with the geometry the WA arrays
offer.**  12 sites, six of them at 21–36 km, moveout spread 9.81 s.  Three deviations, all
frozen before launch in `results/freeze_step1c.json`:

1. **Window from the literature.**  No tremor catalogue exists for 2005 (the PNSN product on
   disk starts 2010-01-06) and the Slow Earthquake Database's WTD catalogue could not be
   retrieved.  The 2005-09 southern-Vancouver-Island ETS dates come from the published
   14-month sequence.  No Bostock event time is used to choose the window, the hours or the
   quiet days — only inside the scorer.
2. **Hours from the array-energy flag** (`w_energy_hours.py`), which Merlin's Q3 sanctions
   for hour selection only: hourly network-median band energy over each site's own
   quiet-day median, hot above the quiet p95 (2.26) → 299 candidate hours on 19 days.
   Section 2 D7 measures what that flag actually caught.
3. **The 2012 template reference is transferred unchanged.**  The relative clauses of the
   scorer (recall ≥ 0.5×, purity ≥ 0.8× the frozen template detector's own per-family
   figures) have no 2005 reference — the frozen templates are built on stations that do not
   exist in 2005 — so the 2012 bar was transferred, fixed before scoring, with the absolute
   recall/purity written alongside in `*_isolated_detail.csv`.

Two engine defects were found and fixed *before* the runs that section 0 reports, and are
recorded here because they invalidate earlier numbers: (i) hour cuts were indexed from each
cache file's own `t0` with exact-equality matching across stations, and the frozen cache
stores midnight + 0.0000 s (CN/C8) / +0.0084 s (PB) / +0.0150 s (LRIV), so whole stations
were dropped from every hour — the first 1a ran on 5 of 6 sites and the first 1b on 3, below
`MIN_SITES`; (ii) `v_engine.DayData` reuses one station's rFFT length for every station.

## 2. Diagnostics D0–D7 (DIAGNOSTIC — no verdict is changed)

### D0 — the crashed diagnostic stage, fixed, and run for 1b and 1c

`w_diag.py` died at the very end of both the 1b and the 1c queue with
`KeyError: 'BPCB'` / `KeyError: 'GLBC'` inside `s_common.fam_tt`: that function defaults to
the `lfe_picker` P3D roster station table, which contains none of the Step-1b/1c bed
stations.  Fixed by passing `w_common.station_table()` — the same roster *extended* with
those stations' FDSN coordinates that every other engine stage already uses.  Two further
defects surfaced once it ran: D1 cut its self-correlation window at the family's median
event over the whole bed, which on a bed with thousands of events per family lies outside
the loaded hour (D1 came back empty); and the family's busiest hour could fall outside the
bed's own days.  Both fixed; the fix is confined to `scripts/w_diag.py`, which is a
diagnostic and is not called by any discovery stage.

| | 1a (2012 permanents) | 1b (2012 + LZB) | **1c (2005 POLARIS)** |
|---|---|---|---|
| D1 median z_self (perfect-repeat ceiling) | 14.4 | 14.7 | **14.2** |
| D2 median true-pair z (own-S cut) | 2.13 | 2.15 | **2.11** |
| D2 median jitter-null z | 1.98 | 2.00 | **1.98** |
| D2 fraction of true pairs above their own null max | 0.16 | 0.15 | **0.17** |
| D3 pair-invariant matches / pairs tested | 0 / 2000 | 3 / 2000 | 544 / 2000 |
| D3 lag-jitter null matches | [0, 0, 0] | [0, 1, 2] | [98, 81, 94] |

D1 is the ceiling of the whole statistic: a window correlated with *itself* reaches z ≈ 14,
so a real pair needs single-station cc ≳ 0.35 to reach z = 5.  D2 is what a real pair
actually reaches at the most favourable cut that exists (each station cut at that family's
own predicted S).  **The 2005 dense bed does not move either number**: the true-pair z at
21–36 km is the same ≈2.1 against a ≈2.0 jitter null as at 46–91 km.

**D3 separates the two beds, and it is the one encouraging number in this whole file.**
On the 2012 permanent bed none of the 2,000 strongest reported pairs had a lag equal to
a kept family's inter-event interval at a common offset (null 0/2/0).  On the 2005 dense
bed **544 of 2000 do, against a lag-jitter null of 91** — a
**6.0× enrichment**.  So the pairs the engine records above its 5.0 engineering
floor on this bed are substantially *real LFE pairs*: the pairing stage works on array
geometry.  What removes them is the threshold (D4), not the pairing.

### D4 — why T = 12.50 on the 2005 bed

T is derived on the quiet hours by the FA-only rule (≤1 false pair per quiet hour AND zero
quiet-hour components of size ≥5), floor 5.0.  On this bed the component clause is satisfied
already at T = 6.50; **the pair-count clause is what pushes T to 12.50, and it is driven by a
single quiet hour**.

Leave-one-quiet-hour-out (the frozen T is **not** changed by this; it is a diagnostic):

| dropped_date | dropped_hour | n_pairs_dropped | T_without | binds |
|---|---|---|---|---|
| 2005-07-09 | 23 | 1669 | 7.5 | fa_rule |
| 2005-06-25 | 0 | 105 | 12.5 | fa_rule |
| 2005-06-26 | 0 | 297 | 12.5 | fa_rule |
| 2005-06-25 | 23 | 278 | 12.5 | fa_rule |
| 2005-06-27 | 0 | 194 | 12.5 | fa_rule |

Dropping **2005-07-09 23** alone takes T from
12.50 to 7.50; dropping any other
quiet hour leaves T at 12.50.  Above net_mad 8 that one hour supplies
**100 % of all quiet-hour pairs**, and every one of them lies inside a single ~70-second
burst (2005-07-09 23:48:03–23:49:10), on 9 of 9 observable sites, at |lag| 13–56 s — the
autocorrelation of one transient against its own coda.  It is a regional/teleseismic
earthquake in a "quiet" hour, which the literature-quiet day rule cannot exclude because no
2005 tremor or local-earthquake screen was applied to the quiet days.

Pooled pairs and the tremor-vs-quiet excess against T:

| T | quiet_pairs | quiet_per_hour | tremor_pairs | tremor_per_hour | ratio | quiet_hours_with_pairs | quiet_max_component |
|---|---|---|---|---|---|---|---|
| 5 | 3910 | 97.75 | 143524 | 480 | 4.911 | 28 | 110 |
| 5.5 | 1385 | 34.62 | 55596 | 185.9 | 5.37 | 20 | 17 |
| 6 | 697 | 17.43 | 24057 | 80.46 | 4.617 | 8 | 7 |
| 6.5 | 502 | 12.55 | 12101 | 40.47 | 3.225 | 5 | 4 |
| 7 | 419 | 10.47 | 6622 | 22.15 | 2.114 | 3 | 2 |
| 7.5 | 361 | 9.025 | 3833 | 12.82 | 1.42 | 3 | 2 |
| 8 | 301 | 7.525 | 2292 | 7.666 | 1.019 | 2 | 3 |
| 8.5 | 268 | 6.7 | 1385 | 4.632 | 0.6914 | 1 | 2 |
| 9 | 231 | 5.775 | 875 | 2.926 | 0.5067 | 1 | 2 |
| 9.5 | 179 | 4.475 | 577 | 1.93 | 0.4312 | 1 | 2 |
| 10 | 143 | 3.575 | 358 | 1.197 | 0.3349 | 1 | 3 |
| 10.5 | 116 | 2.9 | 234 | 0.7826 | 0.2699 | 1 | 4 |
| 11 | 80 | 2 | 141 | 0.4716 | 0.2358 | 1 | 4 |
| 11.5 | 61 | 1.525 | 89 | 0.2977 | 0.1952 | 1 | 4 |
| 12 | 50 | 1.25 | 47 | 0.1572 | 0.1258 | 1 | 3 |
| 12.5 | 38 | 0.95 | 25 | 0.08361 | 0.08801 | 1 | 3 |
| 13 | 15 | 0.375 | 8 | 0.02676 | 0.07135 | 1 | 3 |
| 13.5 | 4 | 0.1 | 1 | 0.003344 | 0.03344 | 1 | 3 |

The tremor-hour rate exceeds the quiet-hour rate by ~5× at the engineering floor, and the
ratio crosses 1 at T ≈ 8.25 — *because the quiet tail is
that single earthquake, not because tremor hours go quiet*.  At the frozen T = 12.50 only 25
tremor pairs survive, spread over 4 components, and the ≥5-member seed-cluster rule then
returns nothing.

Who makes the pairs (per-site z of the pairs the engine itself reported, recomputed on the
GPU at the stored pair times):

| label | n_pairs | hours | med_net_mad | med_n_sites_z3 | med_abs_lag_s | frac_lag_lt_20s | top_site | top_site_z | top_site_share_of_sum |
|---|---|---|---|---|---|---|---|---|---|
| quiet_busy | 120 | 2 | 6.181 | 3 | 1089 | 0.025 | MGCB | 4.964 | 0.3228 |
| quiet_high_netmad | 120 | 2 | 9.854 | 4.5 | 17.2 | 0.6833 | GLBC | 5.307 | 0.2696 |
| tremor_busy | 240 | 4 | 7.983 | 3 | 43.3 | 0.2833 | TWBB | 6.463 | 0.2797 |

No single site or station dominates: the strongest site carries only ~27–32 % of the clipped
site sum in every class, 8–11 sites are observable, and 3–4.5 sites clear z ≥ 3.  So the
quiet-hour pairs are **network-coherent transients, not site-coherent noise and not one bad
station** — which is why the Step-3 stop rule ("≥5-member components at every T ≤ 20") does
not fire.  The high-net_mad quiet pairs sit at |lag| ≈ 17 s with 68 % below 20 s, i.e. just
above the 10-s exclusion, as a single impulsive arrival correlating with its own coda must.

### D4b — an engine defect found by the diagnostics: hours voided by a single dead channel

`w_autocorr` normalises each station's cross-correlation by the per-window energy
`E = ew.clamp_min(1e-10 * ew.median()).sqrt()`.  When more than half of a station-hour is
flat — a gap filled with zeros, a dead component — `ew.median()` is 0, the clamp is a no-op,
`E` is 0, the cc is ±inf and the station z is NaN.  `sbuf += z * vm` then makes
`NaN × 0 = NaN`, so the *valid-sample mask cannot rescue it*: the site mean is NaN, the
equal-weight site sum is NaN at **every** lag, `hit = (nm >= T) & …` is False everywhere,
and the **whole hour silently records zero pairs**.  A single degenerate member of a
multi-station site is enough.

| | tremor hours | quiet hours |
|---|---|---|
| hours on the 2005 bed | 299 | 40 |
| hours with a degenerate station at some site | 26 (8.7%) | 0 |
| hours that recorded zero pairs | 26 | 1 |
| both (the defect firing) | **26** | 0 |
| zero pairs without the defect | 0 | 1 |
| voided but still recorded pairs | 0 | 0 |

Offending sites (2005 bed): `B006q` × 16, `B006q,SHVB` × 8, `B003,B006q` × 2.

The same scan on the Step-1b bed: 0 of 563
tremor hours voided by the same mechanism — so this is a property of the statistic as
implemented, not of the 2005 data.

The correspondence is exact — every hour that recorded zero pairs is a voided hour and
every voided hour recorded zero — so the mechanism is confirmed, not inferred.  **This is a
code defect, and it costs the 2005 run
8.7% of its tremor hours.**  It is reported here as a
defect to repair, not as an explanation of the result: D5 (below) shows that the hours which
*did* run contained catalogued pairs the statistic could not recover either.

### D5 — can a *real* pair be recovered on this bed at all?

For the 99 KEPT families with ≥50 Bostock events
in the window, 934 same-family event pairs inside one hour were
measured **exactly as `w_autocorr.autocorr_hour` measures a pair** — per-station 3-C cc of
the 10-s window against the whole hour, lower-half-MAD z over valid lags, per-site mean
clipped at 8, equal-weight site sum, per-row median/lower-half-MAD renormalisation — at the
true lag.  Two cuts are reported: the frozen **common absolute window**, and an **oracle
moveout-aware cut** in which every station is cut at that family's own predicted S and the
site z are summed on the common inter-event-lag axis (the statistic Merlin's 1a audit asked
about).  The ±0.5 s columns take the maximum within ±0.5 s of the true lag, against a
same-width null drawn at 20 random lags in the same hour.

| statistic at the true lag | frozen common window | oracle moveout-aware cut |
|---|---|---|
| median net_mad (exact lag) | 1.00 | 0.06 |
| median net_mad (max ±0.5 s) | 2.11 | 2.22 |
| median of the matched jitter null | 1.78 | 1.79 |
| fraction above their own null maximum | 0.152 | 0.184 |
| **fraction reaching T = 12.50** | **0.0000** | **0.0000** |
| fraction reaching 7.75 (the 1a/1b T) | 0.0011 | 0.0021 |
| fraction reaching 5.00 (the engineering floor) | 0.0064 | 0.0128 |
| fraction with ≥2 sites at z ≥ 3 | 0.0064 | 0.0075 |
| maximum over all pairs measured | 7.50 | 7.37 |

Control, same code path: the *self* lag of the same windows reads net_mad
25.5 — the statistic has ample dynamic range on this
bed; it is the real pairs that are absent from it.

Per site, against distance:

| site | dist_km | n | z_med | z_p90 | frac_z_ge3 | zmo_med | zmo_p90 | frac_zmo_ge3 |
|---|---|---|---|---|---|---|---|---|
| GLBC | 28.08 | 934 | 0.2666 | 1.61 | 0.006424 | -0.01727 | 1.358 | 0.001071 |
| SHVB | 28.18 | 934 | 0.09746 | 1.551 | 0.0182 | 0.003961 | 1.435 | 0.01392 |
| TWKB | 30.02 | 773 | 0.3949 | 1.764 | 0.01423 | 0.003294 | 1.496 | 0.009056 |
| MGCB | 30.66 | 934 | 0.3467 | 1.761 | 0.008565 | -0.02431 | 1.359 | 0.006424 |
| TSJB | 31.35 | 934 | 0.3667 | 1.815 | 0.006424 | 0.03934 | 1.399 | 0.005353 |
| TWBB | 36.57 | 934 | 0.3494 | 1.634 | 0.01071 | 0.07902 | 1.464 | 0.005353 |
| PGC | 38.15 | 880 | 0.4043 | 1.76 | 0.01364 | 0.09134 | 1.58 | 0.01364 |
| SSIB | 45.04 | 806 | 0.3929 | 1.719 | 0.007444 | -0.06936 | 1.374 | 0.004963 |
| SILB | 47.67 | 751 | 0.3314 | 1.613 | 0.01332 | -0.01394 | 1.445 | 0.01332 |
| LCBC | 48.79 | 870 | 0.3127 | 1.643 | 0.003448 | -0.00239 | 1.366 | 0.002299 |
| B003 | 57.72 | 934 | 0.2886 | 1.612 | 0.005353 | -0.03184 | 1.402 | 0.008565 |
| B006q | 59.01 | 934 | 0.1586 | 1.232 | 0.003212 | 0.03924 | 1.134 | 0.001071 |

There is **no distance trend inside the bed**: the nearest POLARIS sites return the same
z ≈ 0 median and z ≈ 1.5 ninetieth percentile as PGC and as the boreholes three to four times
further away.  Two 10-s windows containing two occurrences of the same LFE, at 25 km, in the
frozen 20-Hz 1–8 Hz chain, do not correlate above the hour's own noise.

### D6 — sensitivity ablations (diagnostic only; none of these constants is retuned)

On the same true pairs:

| variant | n | med | p90 | mx | frac_ge5 | frac_ge7p75 | frac_ge12p5 |
|---|---|---|---|---|---|---|---|
| dist_weighted | 75 | 1.043 | 2.555 | 5.165 | 0.01333 | 0 | 0 |
| frozen | 75 | 1.042 | 2.383 | 3.371 | 0 | 0 | 0 |
| mad_full | 75 | 1.041 | 2.4 | 3.354 | 0 | 0 | 0 |
| noclip | 75 | 1.042 | 2.383 | 3.371 | 0 | 0 | 0 |
| win15 | 75 | 0.8725 | 2.468 | 3.143 | 0 | 0 | 0 |
| win6 | 75 | 1.37 | 2.516 | 3.444 | 0 | 0 | 0 |
| win8 | 75 | 1.181 | 2.45 | 3.908 | 0 | 0 | 0 |

On four quiet hours, the same variants, counting windows whose best admissible lag exceeds T:

| variant | quiet_windows>=T5 | quiet_windows>=T7.75 | quiet_windows>=T12.5 |
|---|---|---|---|
| dist_weighted | 10.75 | 0 | 0 |
| frozen | 10.75 | 0 | 0 |
| mad_full | 10.5 | 0 | 0 |
| noclip | 10.75 | 0 | 0 |
| win15 | 11.5 | 0 | 0 |
| win6 | 4.5 | 0 | 0 |
| win8 | 6.5 | 0 | 0 |

Band and sample rate.  The frozen cache is 20 Hz, Butterworth-4 zero-phase 1–8 Hz; Bostock's
own network autocorrelation used higher-rate data in a 2–8 Hz band.  One station-hour was
recomputed from the raw 100 Hz mseed already on disk (no download):

_(not available)_

### D7 — did the array-energy hour flag starve the run?

| subset | n_events | frac_in_selected_hours | n_in_selected_hours | frac_in_hours_actually_run | n_in_hours_actually_run |
|---|---|---|---|---|---|
| all_2005 | 57504 | 0.7841 | 45088 | 0.7841 | 45088 |
| all_2005_kept | 42890 | 0.7809 | 33493 | 0.7809 | 33493 |
| window | 56181 | 0.8025 | 45088 | 0.8025 | 45088 |
| window_kept | 41808 | 0.8011 | 33493 | 0.8011 | 33493 |

**No.**  80.1% of the kept-family Bostock events
inside the window fall in the 299 hours the flag selected, and all
299 of those hours were actually processed.  The flag is also selective in
the right direction: the median selected hour holds
92.0 catalogued events against
30.5 in an unselected one.  The 60 % clause of the
diagnostic plan is comfortably met, so hour selection is not the explanation for the null
result.  (57,504 Bostock events exist in the 2005 window, so same-family pairs within an
hour are abundant — D5 measured 1,200 of them.)

## 3. What this does and does not show

**Shown.**

- The three runs are a genuine engine test, not a plumbing failure: every stage executed, the
  FA rule bound above the engineering floor on all three beds, the site-shift nulls returned
  0/0/0 on 1b and 1c, and the self-lag control reads net_mad ≈ 26, so the statistic has its
  full dynamic range on the 2005 bed.
- The 2005 threshold T = 12.50 is set by **one quiet hour containing one earthquake**, not by
  array noise and not by any site: leave-one-hour-out moves T to 7.50 and no other hour moves
  it at all.  The quiet-day rule for 2005 was literature-only and screens tremor, not
  transients.
- Catalogued same-family pairs are **not recoverable by this statistic on this bed at any
  threshold that the quiet hours would tolerate**, with or without the moveout-aware cut, and
  with no improvement from the 21–36 km sites over the 46–91 km ones.  The single-pair SNR,
  not the geometry and not the threshold, is what is missing.
- The array-energy hour flag caught ~80 % of the catalogued events; hour selection is not the
  cause.
- **The pairing stage works on this bed.**  544 of the 2,000 strongest pairs the engine
  reported carry a kept family's inter-event interval, against a lag-jitter null of
  91 — a 6.0× enrichment, where the 2012 permanent bed gave 0
  against 0.  The information the chain needs is present in the pair list at the
  engineering floor and is thrown away by the threshold.

**Not shown.**

- **This is not evidence that array discovery is impossible.**  Bostock et al. found these
  families on this geometry.  What is measured here is one statistic (single-pair network
  autocorrelation of 10-s windows) on one data chain (20 Hz, 1–8 Hz, 3-C, frozen cache) on
  one bed.  A negative for that triple is not a negative for array discovery, for the WA
  arrays, or for Arm B.
- Nothing here licenses a retune of the frozen constants.  The numbers above are
  measurements, not proposed bars; T, the window, the hop, the clip and the seed rule stay
  exactly as frozen until a repair is designed and audited.
- The 1c withheld set lost LZB (no 2005 data), so 1c's Gate-1 power rested on JRBC plus
  partial VGZ/SNB.  Nothing downstream ran, so this did not bind, but it would have.

- The D4b defect is a genuine code bug, but at 8.7% of tremor hours it cannot
  account for a 0/10, and repairing it does not make any of the D5 pairs recoverable.


**The single change the diagnostics point to.**  The quiet-hour set — and therefore T — is
the one place where a single, local, testable defect was found: the 2005 quiet days were
selected by literature distance from the ETS only, with no screen for regional or
teleseismic earthquakes, and one such event in one hour multiplies T by 1.7 and removes the
entire tremor-hour pair population.  A quiet-hour *robustness* rule (an hour-level outlier
screen, or a per-hour FA quantile instead of a pooled mean) is the change the evidence names.
It is stated here as the finding, **not as an authorised retune**: with T = 7.50 the 1c run
would have had 3,833 tremor pairs to cluster instead of 25, and D3 says a large minority of
the pairs at that level really are LFE pairs — but D5 says the *individual catalogued* pairs
still do not reach 7.75, or even 5.00, at their own true lag.  D3 and D5 are both right and
point at the same place: the pair list is enriched in the aggregate while no single pair is
individually significant, which is what a threshold at the edge of the noise produces.
Repairing T is therefore necessary to make the 1c run *informative* and is not by itself
enough to make it *pass*; a statistic that accumulates over many pairs, rather than one that
thresholds each pair, is the direction this evidence points, and that is a redesign, not a
retune.  Any repair cycle needs Merlin's audit first.

## 4. Files

Frozen runs: `results/freeze_step1{,b,c}.json`, `results/step1{,b,c}_score.json`,
`results/step1{,b,c}_score_{rediscovery,splits,isolated_detail}.csv`,
`results/step1/`, `results/step1b/`, `results/step1c/` (T.json, T_grid.csv, ac*.parquet,
cl_seeds.csv, bank*.npz, p*_survival.csv, dedupe.csv, lg_*, null1-3/),
`results/STEP1_DONE`, `results/STEP1B_FAILED`, `results/STEP1C_FAILED` (the diag-stage
tracebacks), `logs/step1{,a_refresh,b,c}_queue.log`.

Diagnostics (this file): `scripts/w_diag.py` (fixed), `scripts/y_step1diag.py`,
`scripts/y_write_step1_results.py`, `scripts/y_step1diag_queue.sh`,
`results/step1{,b,c}_diag.json` + `_D1_dynamic_range.csv` + `_D2_true_pair_z.csv`,
`results/step1diag_D4{.json,_per_hour.csv,_ratio_vs_T.csv,_leave_one_hour_out_T.csv,
_site_attribution.csv}`, `results/step1diag_D4b{.json,_voided_hours.csv}`,
`results/step1bdiag_D4b*`, `results/step1diag_D5{.json,_true_pairs.csv,
_site_z_vs_distance.csv,_by_site.csv}`, `results/step1diag_D6{.json,_true_pair_variants.csv,
_true_pair_summary.csv,_quiet_variants.csv,_native_rate.csv,_native_summary.csv}`,
`results/step1diag_D6native.json`, `results/step1diag_D7{.json,_hour_coverage.csv,
_events_per_hour.csv}`, `logs/step1diag_queue.log`.

---

## Appendix A — the superseded 2026-09-20 15:46Z Step-1a write-up, verbatim

# Step 1a — north PERMANENT bed (6 sites, 46–80 km): the autocorrelation stage has no SNR margin

> **SUPERSEDED NUMBERS (added 2026-09-20 15:xxZ by the orchestrator).** The run described below used 9 of 12 stations / 5 of 6 sites because of the hour-cut t0 bug (cache midnight offsets +0.0084 s PB / +0.0150 s LRIV; fixed 2026-09-20). The corrected re-run (`w_step1a_refresh.sh`, summary `results/STEP1_DONE`, 15:06Z) on the full 12 stations / 6 sites gives the same verdict: T = 7.25 (FA rule), 563 tremor hours → 199,579 pairs ≥5.0 → 3 seed clusters → 0 pass-0 survivors → 0/10 isolated hidden families (ceiling 5/10); D2 true-pair z 2.13 vs null 1.98. Diagnostics D1–D3 below stand. Full Step 1a/1b/1c write-up follows when 1c completes.

### *uninformative about the engine; informative about distance*

Run 2026-09-20, agent E. Frozen before launch in `results/freeze_step1.json` (sha256 of 12 `w_`
scripts + 13 imported lfe_picker scripts; nothing in lfe_picker edited). Engine per PLAN §10 /
Merlin Q1(a). **No constant was retuned during or after the run.**

## 1. What was run

`lfe_catalog_repair/cache20hz`, 2012-08-15..10-05, 52 day-files, **618 station-days**,
12 detection stations merged by the data-driven 1-km rule into **6 sites**; VGZ/SNB/B926 withheld
(Gates suite). Chain: network autocorrelation → clusters → align_stack handoff → matched-filter
growth → survival + template-scramble null → 2 iterations → dedupe → locate → Gates 0–2.
**The three site-shift null realisations were NOT run** (coordinator decision: the outcome is
already determined at the autocorrelation stage, and a null on a chain that produces nothing
measures nothing).

| site | members | median hypocentral distance to the 130 Bostock families |
|---|---|---|
| PGCq | B009+B010+B011+PGC | 52.6 km (40–76) |
| LRIVq | B005+B006+B007+LRIV | 64.9 km (33–117) |
| B003 | B003 | 70.1 km (40–108) |
| B004 | B004 | 74.2 km (46–111) |
| B001 | B001 | 77.4 km (44–133) |
| NLLB | NLLB | 91.4 km (55–137) |

## 2. Engine prerequisites (Merlin F3) — all satisfied, and the FA rule BINDS

| prerequisite | state |
|---|---|
| P1 no `T = 8.0` floor | **T = 7.75 derived by the FA rule**, engineering floor 5.0 never binding; growth-stage `T_seed` 6.25 / 7.75 with `floor_binds = False` for both templates |
| P2 data-driven `SITE_MEMBERS` | 1-km single linkage reproduces the hand-written merge exactly (6 sites) |
| P3 `station_z` over valid samples only | implemented in the `w_` subclass |
| P4 `MIN_SITES` 4, availability never intersected | 6/6 sites usable on all 52 days |

**T derivation** (45 quiet hours on 9 quiet days — the bed has only 9 days ≥3 days from any PNSN
tremor day within 100 km, so Merlin's "≥20 quiet days × 2 h" was met as "≥40 quiet hours" with 5
hours per quiet day; deviation recorded in the freeze):

| T | pairs | false pairs / quiet hour | quiet components ≥5 |
|---|---|---|---|
| 6.50 | 853 | 18.96 | 6 |
| 6.75 | 436 | 9.69 | 2 |
| 7.00 | 216 | 4.80 | 0 |
| 7.50 | 61 | 1.36 | 0 |
| **7.75** | **31** | **0.69** | **0** |

## 3. Result

| stage | count |
|---|---|
| candidate tremor hours (≥5 PNSN epicentres ≤50 km of an active cell) | 563 over 32 days |
| autocorrelation, median cost | 0.81 s / station-hour (H100, 3.8 GB), 11 min total |
| pairs recorded at net_mad ≥ 5.0 | 248 397 |
| pairs at T = 7.75 | 3 740 → **6.64 / tremor hour vs 0.69 / quiet hour (9.6× excess)** |
| connected components → **seed clusters** (≥5 in an hour, ≥2 ten-minute bins) | **2** |
| pass-0 templates with ≥4 sites | 2 |
| growth hits over 52 days | 16 (12 + 4), 9 + 4 independent |
| **survivors** (≥20 independent, ≥3 days, ≤20 % near seed, scramble ≥3×) | **0** |
| iterations / dedupe / locate / Gates 0–2 | 0 templates throughout |
| **isolated hidden families rediscovered** | **0 / 10** (ceiling: the frozen template detector itself 5/10) |
| spurious-split rate | 0.00 |

Verdict from `w_score_step1.py`: `ENGINE_INADEQUATE` **on this bed** — i.e. Merlin's
"<2/10 with null ≤1" failure signature. Accept bars were ≥3/10, split ≤0.1, null ≤1.

There *is* a real tremor-related excess of coherent pairs (9.6× more at T than in quiet hours), so
the statistic is not inert; those pairs simply do not organise into repeating families.

## 4. Why — three diagnostics (`results/step1_diag.json`, `w_diag.py`, KEPT families only)

**D1 — dynamic range.** A window correlated with itself (a perfect repeat) scores
**z_self = 14.4**, because the lower-half MAD of an hour's cc distribution is ≈0.07. So a pair of
occurrences needs **cc ≳ 0.35** to reach z = 5 at one station. That is the ceiling of the whole
statistic on this bed.

**D2 — true-pair z.** For the five largest kept families, at the TRUE lag of catalogued same-family
event pairs, cut moveout-aware at the family's own predicted S (the most favourable cut that
exists): **median z = 2.13** against a same-width jitter null of **1.98**; only **16 %** of true
pairs exceed their own null maximum. Measured at window lengths 3, 4, 6, 8, 10, 15 and 20 s the
site-median z is flat at 1.6–2.4 — window length is not the free parameter that rescues this.
Real LFE pairs here sit at cc ≈ 0.10–0.15, a factor 2.5–3.5 below what z = 5 requires.

**D3 — pair invariant.** Of the 2 000 strongest pairs the engine reported, **0** have a lag equal to
an inter-event interval of a kept family at a common offset; the lag-jitter null gives 0, 2, 0. The
pairs that pass T are not LFE pairs.

**Geometry.** The bed's S-moveout spread across the 6 sites is **median 13.8 s** (9.4–24.3 s;
93 % of families > 10 s). Merlin sized the 10-s common window for the **45-km-aperture WA array**
(5–6 s spread). On a 46–137 km permanent network a single common window cannot contain S at all six
sites for 93 % of the family field, so four of six sites contribute noise to every sum — an
additional √3 loss on top of the cc deficit.

## 5. What this does and does not establish

- It does **not** establish that the Arm-A engine cannot discover families. Every stage ran, the FA
  rule bound above the engineering floor, the handoff/growth/iterate/dedupe/locate/gates path was
  exercised end to end (smoke tests and this run), and the failure is located precisely at the
  first stage, for a measurable physical reason.
- It **does** establish that **single-pair network autocorrelation needs stations at LFE distances
  of order 10–40 km**, not 46–137 km. This is the P4-0b starvation lesson in a new place: the north
  PERMANENT network is the wrong bed for a positive control, because Bostock's families were
  discovered on temporary stations that this bed does not contain.
- Consequence: Step 1 is re-run as **Step 1b** on the same 2012 window with the POLARIS/CN dense
  stations added (`results/freeze_step1b.json`). Constants are carried over unchanged — the window,
  hop and T are NOT retuned before that result is seen.
- **For Merlin to audit**: (i) is the 10-s common window still the right choice once the bed has
  10–40 km stations (the WA arrays' own aperture is 40 km, spread ≈5–6 s, so probably yes); (ii)
  should the statistic fall back to a moveout-aware, cell-scanned pairing on wide-aperture beds;
  (iii) does the 9.6× tremor excess of coherent pairs carry usable information that the
  components/seed-cluster rule discards.

## 6. Defect found in the inherited engine

`v_engine.DayData` stores a single `self.NF` (the last station's rFFT length) and uses it for every
station in `station_z`. It is harmless only while every station-day has an identical sample count;
a short or partial day raises `RuntimeError: The size of tensor a ... must match ...`. Fixed in the
`w_` subclass (per-station `NFs`); **`v_engine.py` on disk still carries it** — for the P4-0b
close-out.

## 7. Files

`results/freeze_step1.json`, `results/step1_score.json`, `results/step1_diag.json`,
`results/step1_diag_D1_dynamic_range.csv`, `results/step1_diag_D2_true_pair_z.csv`,
`results/step1_geometry.json`, `results/step1/` (T.json, T_grid.csv, ac*.parquet, cl_seeds.csv,
bank*.npz, p*_survival.csv, dedupe.csv, lg_*), `logs/step1_queue.log`, `logs/step1a_tail.log`.

