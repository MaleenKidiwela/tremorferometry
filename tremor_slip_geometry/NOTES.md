# NOTES — tremor_slip_geometry (append-only ledger)

## 2026-09-10 setup
Folder created; PLAN.md moved here from notes/. User directives: separate folder, document everything, use Opus sub-agents for
coding/simple tasks, no dv/v in this study yet. Steps S1-S6 per PLAN.md.

## 2026-09-10 orchestration — agents launched
Three Opus sub-agents launched in parallel from specs written by the main session (specs mirror PLAN.md sections 1-5):
- S1 co-registration (scripts/s1_coregistration.py): per confirmed north/central patch Q1-Q5 in the local strike/dip frame,
  mirror + time-shifted-tremor controls, Theil-Sen trends with permutation p.
- S2 depth occupancy (s2_depth_occupancy.py) + S3 efficiency maps (s3_efficiency_maps.py): yearly tremor/slip depth percentiles
  and edges with completeness proxy; per-triangle slip-per-detection by epoch (raw + network-normalised), change maps, profiles.
- S6 model-only checks (s6_model_checks.py): inter-ETS slip fraction and burst confirmation vs tremor (shift null); southern
  +-5 mm oscillation: sign symmetry, spectrum, seasonal amplitude, band coherence, dominant vbICA components, stationarity.
S4 (fronts) and S5 (site recurrence) to be specified after reviewing S1 (needs the confirmed-event list and the frame checks).

## 2026-09-10 S1 co-registration
Ran `scripts/s1_coregistration.py` (OMP_NUM_THREADS=4, seed 20260910, ~85 s end to end; two identical runs, the second
only added the group-level z test). Per-event Q1-Q5 for the 27 CONFIRMED thrust patches (catalog=='thrust' & confirmed &
seg in {north, central} & Mw>=6.0 & not broad_pulse; 21 north / 6 central) and for the 72 MIRROR controls
(catalog=='mirror' & Mw>=6.0 & seg in {north, central}). Slip from `-(Ld @ A_raw.T)` (mirror: `+`), accumulated
`ts-30 .. te` and clipped at >=0; footprint F = triangles labelled with the patch id on any day of the event; tremor =
PNSN detections inside the event window and within 40 km of an F centroid; depths = Slab2 at the epicentre (finite for
93.6 % of the 748,002 detections); local down-dip/strike frame from the +-0.05 deg finite-difference slab gradient at the
slip centroid. Sanity check before the run: `slab_depth_interpolator` at the 742 triangle centroids reproduces the npz
`patch_depth_km` to 1.8e-15 km, so the two depth definitions are the same object.

Printed output (verbatim):

    --- CONFIRMED (thrust) (n=27) ---
    group         n        Q1        Q2        Q3       Q4a       Q4b   log10Q5   fQ1>0    medN    medPot
    combined     27      4.68     23.45      0.94     11.03     -5.25     -4.37    0.81    3817    0.1272
    north        21      4.68     23.45      1.05     11.38     -5.25     -4.32    0.86    5085    0.2039
    central       6      4.67    -11.26      0.43      9.67     -3.50     -4.43    0.67    2751    0.1000

    --- MIRROR control (all) (n=72) ---     [Q-medians skip the 32 rows with N=0]
    combined     72     18.61     18.95      5.29     13.01     -5.68     -3.06    0.82      10    0.0572
    north        64     20.53     -0.42      5.29     11.83     -5.08     -3.06    0.85       5    0.0544
    central       8     13.94     77.34      5.13     24.07     -5.68     -3.10    0.67     222    0.2047

    --- spread comparison (IQR of the offsets; S1 decision point) ---
    Q1_downdip_km      thrust med    4.68 IQR   8.98 (n=27)   mirror med   18.61 IQR  28.58 (n=40)
    Q2_strike_km       thrust med   23.45 IQR  52.19 (n=27)   mirror med   18.95 IQR 125.24 (n=40)
    Q3_depth_km        thrust med    0.94 IQR   1.66 (n=27)   mirror med    5.29 IQR   6.18 (n=40)

    --- time-shift null z-scores (confirmed) ---
    z1: n=27 median -0.52 mean -0.48 (SE 0.12)  |z|>2 in 0/27   group t p=0.00061; sign test 22/27 negative p=0.0015
    z3: n=27 median -1.07 mean -0.92 (SE 0.12)  |z|>2 in 0/27   group t p=3.4e-08; sign test 24/27 negative p=4.9e-05
    thin nulls (fewer than 50 of 100 shifts gave any tremor): 2 events -> [1434, 1655]

    --- TRENDS (Theil-Sen per year, permutation p over 1000 date shuffles) ---
    set        seg       quantity             n   slope/yr               95% CI  p_perm    rho    p_sp
    confirmed  combined  Q1_downdip_km       27    -0.1388 [ -0.8250,  0.5747]   0.644  -0.14   0.479
    confirmed  combined  Q3_depth_km         27    -0.0388 [ -0.1596,  0.0898]   0.579  -0.19   0.351
    confirmed  combined  Q4a_updip_km        27     0.1258 [ -0.3765,  0.6553]   0.560   0.13   0.532
    confirmed  combined  Q4b_downdip_km      27     0.0236 [ -0.4589,  0.5022]   0.888   0.02   0.921
    confirmed  combined  log10Q5             27    -0.0006 [ -0.0251,  0.0217]   0.909  -0.04   0.835
    confirmed  north     Q1_downdip_km       21    -0.1326 [ -0.9056,  0.8661]   0.708  -0.13   0.563
    confirmed  north     Q3_depth_km         21    -0.0350 [ -0.1541,  0.1444]   0.680  -0.19   0.414
    confirmed  north     Q4a_updip_km        21    -0.0205 [ -0.5716,  0.5762]   0.967  -0.03   0.898
    confirmed  north     Q4b_downdip_km      21     0.0232 [ -0.6075,  0.4913]   0.875   0.04   0.871
    confirmed  north     log10Q5             21     0.0003 [ -0.0271,  0.0198]   0.960   0.01   0.955
    confirmed  central   (n=6, all p>=0.14; Q3 -0.206/yr p=0.188, Q4a +0.784/yr p=0.144 -- not interpretable at n=6)
    mirror     combined  Q1_downdip_km       40     2.8477 [  0.8645,  4.8558]   0.013   0.37   0.018
    mirror     combined  Q3_depth_km         40     0.4997 [  0.1149,  0.8416]   0.014   0.34   0.034
    mirror     combined  Q4a/Q4b/log10Q5     40     p_perm 0.167 / 0.327 / 0.697
    (full table incl. mirror-north/central in results/s1_trends.csv)

### Reading
1. **H1 holds and the offset is TIGHT.** Tremor sits a median **+4.7 km down-dip** of the slip centroid (81 % of events
   positive; 86 % in the north), median depth offset **+0.94 km deeper** than the slip-weighted depth. Both are small
   compared with the ~28 km triangle size -- i.e. tremor and modeled slip are co-located to well inside one mesh cell,
   with a small consistent down-dip bias.
2. **PLAN's S1 DECISION POINT: the offsets are NOT mesh geometry.** The mirror (back-slip) offsets are both larger and
   ~3x more scattered: Q1 med 18.6 km IQR 28.6 vs thrust 4.7 / 9.0; Q3 med 5.3 km IQR 6.2 vs 0.94 / 1.66. So the
   geometry line survives -- continue to S2.
3. **Time-shift null: the co-registration is real but modest per event.** No single event exceeds |z|=2, yet the group
   is decisively displaced: z1 mean -0.48 (t p=6e-4, 22/27 negative), z3 mean -0.92 (t p=3e-8, 24/27 negative). Sign
   convention: during its own event, tremor sits LESS far down-dip and LESS deep relative to the slip than the same
   40-km patch neighbourhood does at random other times. Interpretation: the event-time tremor pulls up-dip toward the
   slip, i.e. the association is genuine; the background tremor around these patches is on average deeper.
4. **H2 (secular drift): nothing, in the confirmed set.** Every combined and north trend has p_perm >= 0.56 and a CI
   straddling zero; the tightest bound is Q3 at -0.039 km/yr [-0.16, +0.09], i.e. **no deepening or shallowing of the
   tremor-slip depth offset beyond ~0.16 km/yr (~2.6 km over 16 yr)**. Q1 is bounded at |slope| < ~0.83 km/yr.
5. **H3 (efficiency): flat.** log10 Q5 slope -0.0006/yr [-0.025, +0.022], p=0.91 -- slip per detection is constant to
   within ~6 %/yr over 2010-2026 for confirmed events. Note the confirmed log10Q5 (median -4.37) is 1.3 orders BELOW
   the mirror (-3.06): mirror "events" have far fewer detections, as they should.
6. **The mirror set DOES trend (Q1 +2.85 km/yr p=0.013, Q3 +0.50 km/yr p=0.014).** This is the important negative
   control result: an apparent secular drift of exactly the kind H2 predicts arises in the NON-event catalog. Any future
   drift claim must beat this. Most plausible cause is PNSN network growth (more, and more widely distributed,
   detections after ~2019) interacting with near-empty mirror windows -- 32/72 mirror rows have N=0 and are excluded
   from every median and trend, and those that survive are dominated by late-period events.

### Issues / caveats
- **Small n.** 27 confirmed events, only 6 central. Central-segment trends (n=6) are reported but are not interpretable;
  a Theil-Sen slope on 6 points with a 1000-permutation p has no power.
- **Chained footprints (PLAN 6).** Q2 is large and erratic (thrust north med +23 km, mirror central +77 km, IQR 52-125 km):
  footprints that chain along strike put the tremor and slip centroids far apart along strike even when they co-register
  across strike. Q1/Q3 use slip-WEIGHTED centroids so they are much better behaved, but **Q2 should not be read as an
  along-strike lead/lag** -- that is what S4 (fronts, per-day tracks) is for.
- **Network growth is unmodelled here.** N ranges 75 to 10,880 across events; Q5 = potency/N is therefore a detection-rate
  quantity as much as a physical one. The flatness of log10Q5 is reassuring, but S3 must do the per-band-year
  normalisation before any efficiency claim.
- **Mirror windows are often empty** (32/72 with N=0), so the mirror comparison rests on 40 rows, biased late.
- **Two thin nulls** (ids 1434, 1655): fewer than 50 of the 100 shifted windows contained any tremor in the patch
  neighbourhood, so their z1/z3 are noisy. 1434 is also the smallest event (N=75) and the strongest Q1/Q3 outlier
  (-10.4 km, -2.8 km) -- treat it as unreliable.
- **6.4 % of detections have no Slab2 depth** (outside the mesh convex hull) and are dropped from the depth
  percentiles but kept in N and in the map centroid. This is a small inconsistency between Q1/Q2 and Q3/Q4.
- **Slip edges are a "within the latitude band" definition**, not "within F": they can be set by triangles outside the
  event footprint that happen to share its latitudes and exceed 3 mm. That makes the edges systematically wide
  (median up-dip edge ~22 km, down-dip ~49 km), which is why Q4a is strongly positive (+11 km) and Q4b strongly
  negative (-5 km): the tremor band sits well INSIDE the slipping depth range on both sides. Correct per spec, but the
  edge quantities describe the depth extent of the slipping band at that latitude, not of the patch.
- The absolute log10 Q5 level depends on the 40-km tremor radius and on the accumulation window; only relative
  comparisons (thrust vs mirror, and vs time) are used.

Files written: `results/s1_events.csv` (27 rows), `results/s1_controls.csv` (72 rows), `results/s1_trends.csv` (30 rows),
`figures/s1_offsets_vs_time.png`. Nothing outside `tremor_slip_geometry/` was created or modified.

## 2026-09-10 S6 model-only checks

Ran `scripts/s6_model_checks.py` (PLAN.md step S6). Model + PNSN tremor only, no dv/v. Window 2010-01-06..2026-05-27
(5986 d, the tremor span); 10 one-degree bands by triangle centroid latitude ([40,41)..[49,50); 2 triangles at lat<40
dropped). Tremor assigned to the nearest mesh-triangle centroid in a local km frame, kept if <=40 km: 741,176 / 748,002
detections (99.1%). Daily thrust rate = d/dt of `-1e3*(Ld @ A_raw.T)` with a 7-d centred mean; band series are
area-weighted means. CONFIRMED patches = thrust & confirmed & seg in (north,central) & Mw>=6.0 & ~broad_pulse -> 27
(21 north, 6 central). EVENT DAYS per band = [start-30, end+30] of any confirmed patch whose label-footprint touches the
band. Full definitions in the script docstring.

### PART A — inter-ETS slip
Event-day coverage per band ranges 0 d (40-41, 41-42: no confirmed patch reaches them) to 1864 d of 5986 (47-48).

| band | frac of POSITIVE slip on inter-event days | event/inter-event mean positive rate | frac of tremor inter-event | slip bursts (conf) | tremor bursts (conf) |
|---|---|---|---|---|---|
| 40-41 S | 1.000 | – | 1.000 | 8 (0.50) | 32 (0.09) |
| 41-42 S | 1.000 | – | 1.000 | 20 (0.55) | 48 (0.10) |
| 42-43 S | 0.974 | 0.81 | 0.917 | 25 (0.56) | 31 (0.35) |
| 43-44 S | 0.936 | 1.15 | 0.815 | 8 (1.00) | 24 (0.21) |
| 44-45 C | 0.867 | 1.05 | 0.461 | 1 (0.00) | 36 (0.00) |
| 45-46 C | 0.759 | 1.19 | 0.384 | 2 (1.00) | 41 (0.02) |
| 46-47 N | 0.633 | 1.36 | 0.387 | 2 (1.00) | 40 (0.03) |
| 47-48 N | 0.566 | 1.69 | 0.310 | 1 (1.00) | 28 (0.04) |
| 48-49 N | 0.615 | 1.69 | 0.393 | 3 (1.00) | 29 (0.07) |
| 49-50 N | 0.784 | 1.22 | 0.744 | 2 (1.00) | 40 (0.00) |

Aggregates: inter-event share of positive slip 0.704 (north/central bands) vs 0.978 (south); inter-event share of tremor
0.446 vs 0.933. The modeled daily positive rate on ETS days is only ~1.4x the inter-event rate in the north — the model
puts most of its positive slip outside the confirmed-ETS windows.

72 inter-event SLIP bursts (band 30-d slip >3 mm, runs >=5 d): 11 north/central, 61 south; 47 tremor-confirmed
(30-d count at the burst peak >= 3x the band's inter-event median AND >=10 detections). Rigid tremor-shift null, 200 draws:

| group | burst kind | n | observed | null mean | null p95 | p |
|---|---|---|---|---|---|---|
| north/central | slip burst | 11 | 0.909 | 0.299 | 0.545 | 0.005 |
| north/central | tremor burst | 214 | 0.023 | 0.009 | 0.021 | 0.015 |
| south | slip burst | 61 | 0.607 | 0.209 | 0.311 | 0.005 |
| south | tremor burst | 135 | 0.178 | 0.062 | 0.090 | 0.005 |

IMPORTANT caveat found while reading the burst list: the 10 tremor-confirmed north/central inter-event slip bursts are
2010-03-11, 2014-12-10/13, 2018-01-09/18, 2018-01-28, 2018-02-12, 2019-12-16, 2021-02-25/26 — i.e. recognisable ETS
episodes that the strict CONFIRMED filter (Mw>=6.0 & confirmed & ~broad_pulse) removed, not aseismic inter-ETS creep.
Southern bands have zero event days by construction (no confirmed patch is southern), so their "inter-event" bursts are
just the southern ETS episodes. The converse rate is tiny (2-18%) because a 1-degree BAND-MEAN >3 mm is a strict bar for
a localised patch — the asymmetry is a definition effect, not evidence that tremor happens without slip.

### PART B — the southern +-5 mm oscillation
(1) Sign symmetry, all years, triangle-days: south (<44N) frac(inc30>+5 mm) = 0.0301, frac(<-5 mm) = 0.0255, ratio 1.18;
north (46-49.5N) 0.0249 / 0.0095, ratio 2.62. The southern excursions are nearly sign-symmetric — an oscillation, not
slip episodes; the north is 2.6:1 thrust-positive as real ETS should be.

(2) Spectrum (Welch, nperseg = 1461 d): peak period in the 60 d - 3 yr range is 122 d (bands 40-41, 41-42, 42-43) and
86 d (43-44 and ALL six bands 44-50). Identical in north and south -> it is the model/30-d-boxcar timescale, not a
southern feature. Annual power ratio (330-400 d vs 180-300 & 430-900 d) = 0.48/1.91/1.69/0.98 south, 0.70-1.49 north —
no consistent annual line. Seasonal amplitude (month-of-year max-min) 0.50-0.73 mm south (null mean 0.51-0.78) with
p = 0.25 / 0.97 / 0.55 / 0.76; north 0.58-1.09 mm, p = 0.11-0.59. NOT seasonal in any band.

(3) Coherence: mean Pearson r between band-mean inc30 series = 0.721 among southern bands, 0.515 among northern bands,
0.418 south x north. The south is MORE internally coherent than the north. Component ranking: vbICA components 7 and 3
carry 36.3% and 33.1% of the southern inc30 variance (top-3 with comp 11 = 77.0%), while carrying only 10.6% and 12.1%
in the north; the north's dominant component 1 carries 43.9% north but 1.8% south. So the southern signal is essentially
a two-component projection. Correlation of the band series with the sum of the 13 raw amplitudes' 30-d means is -0.07 to
+0.06 (uninformative: the ICA sign of each component is arbitrary, so their raw sum has no physical meaning).

(4) Stationarity: southern frac(|inc30|>5 mm) = 0.0552 (2010-2014) / 0.0520 (2015-2019) / 0.0588 (2020-2026), std
2.43 / 2.41 / 2.51 mm; north 0.0348 / 0.0358 / 0.0330. Flat — no GPS-network-growth step.

### Interpretation (honest)
The southern +-5 mm oscillation is NOT seasonal (no annual power excess, month-of-year amplitude indistinguishable from a
circular-shift null in all four southern bands), NOT per-triangle noise (southern bands correlate at r = 0.72 with each
other, higher than the northern bands do), and NOT an artefact of network growth (its amplitude is stationary to ~5%
across the three epochs). It is a spatially coherent, sign-symmetric, few-component mode: two of the 13 vbICA components
carry two-thirds of the southern variance, and those same components are minor in the north. The most economical reading
is that the vbICA decomposition assigns southern Cascadia a pair of long-wavelength patterns whose amplitudes ring at the
model's own smoothing timescale (~86-122 d) — i.e. common-mode GPS noise projected onto a poorly-constrained part of the
mesh, filtered by the 30-d window. It behaves like model noise with spatial structure, and the near-1:1 positive/negative
ratio means a southern "5 mm patch" carries no more evidence of slip than its mirror does; the south should keep being
treated as below the model's noise floor. On inter-ETS slip: the model does claim most of its positive slip between the
confirmed events (70% in the north/central bands, 98% in the south), and inter-event slip bursts ARE tremor-confirmed far
above the shift null (91% vs 30% north/central, p=0.005; 61% vs 21% south, p=0.005) — but on inspection those confirmed
bursts are ETS episodes that the strict CONFIRMED filter discarded (Dec 2014, Jan-Feb 2018, Dec 2019, Feb-Mar 2021, and
the whole southern ETS sequence), not aseismic creep. So the honest statement is: the high inter-event slip fraction is
mostly a catalogue-completeness effect plus the model's temporal smearing of ETS slip across the +-30 d window, and this
test does NOT establish tremor-confirmed aseismic slip between ETS episodes. A proper version would define event days
from ALL slip patches (or from a tremor-based ETS calendar) rather than from the 27 Mw>=6.0 confirmed ones.

Files: `scripts/s6_model_checks.py`; `results/s6_interets_bands.csv`, `results/s6_interets_bursts.csv` (421 rows),
`results/s6_south_oscillation_bands.csv`, `results/s6_sign_symmetry_year.csv`, `results/s6_south_components.csv`,
`results/s6_band_corr.csv`; `figures/s6_model_checks.png` (A1 inter-event slip/tremor fractions, A2 burst confirmation vs
null, B2a spectra, B2b month-of-year cycle, B3 band correlation matrix, B4 stationarity).

## 2026-09-10 orchestrator review — S1 and S6 (first pass)
S1 (agent report reviewed): H1 holds — tremor sits a median 4.7 km down-dip and 0.9 km deeper than the slip centroid (27
confirmed N/C events; 81 % Q1>0), i.e. INSIDE one mesh cell: at this model's resolution tremor and slip are co-located.
Decision point PASSED: mirror offsets are larger (18.6 km, 5.3 km) and ~3x more scattered -> the offsets are not mesh geometry.
Time-shift null: no single event significant, but the group z is displaced (z3 mean -0.92, 24/27 negative): tremor during an
event is shallower than the background tremor of the same neighbourhood -> inter-ETS tremor is deeper (consistent with the
literature's down-dip inter-ETS tremor). H2/H3 (drift of offsets or efficiency): NOTHING in the confirmed set (all p >= 0.56;
Q3 slope -0.04 km/yr [-0.16, +0.09]). CRITICAL CONTROL: the mirror set DOES drift (Q1 +2.9 km/yr p 0.013, Q3 +0.5 km/yr p 0.014)
— network growth acting on near-empty windows produces exactly the kind of drift H2 predicts. Any future drift claim must beat
that. Q2 (along-strike) is unusable from union footprints (IQR 50-125 km) -> S4 daily tracks launched.
S6 first pass: Part B is informative — the southern +-5 mm oscillation is sign-symmetric (ratio 1.2 vs 2.6 north), NOT seasonal
(p 0.25-0.97), spatially coherent (r 0.72 between southern bands), carried by vbICA components 7 and 3 (77 % of southern
variance incl. a third; those components carry ~11 % in the north), stationary across epochs, peaking at the 86-122 d
timescale of the model's own smoothing -> a common-mode GNSS mode projected onto a poorly constrained part of the mesh. Part A
(inter-ETS slip) was invalid as first defined (CONFIRMED set leaves the south with no event days and drops many real ETS in the
north) — agent asked to redefine event days from all substantive patches (A1) and a tremor ETS calendar (A2) and rerun.
S4 (fronts) and S5 (site recurrence) agents launched with specs.

## 2026-09-10 orchestrator — two S1 observations from the figure
1. log10 Q5 (potency per tremor detection) is remarkably STABLE for confirmed events: -4.1..-4.7 over 16 years (a factor ~3
   range) with slope -0.001/yr — tremor count is a near-constant-ratio proxy for slip potency (~4e-5 km^3 per detection) in
   the north and centre, and the mirror control scatters over 3 decades. Candidate result: tremor productivity per unit slip is
   stationary 2010-2026 at the model's resolution (to be normalised for network completeness in S3 before claiming).
2. The three northern events starting in 2022 have Q1 +20..+29 km and Q3 +3.5..+5.4 km (all other years 0-10 km / 0-2 km):
   tremor sits far down-dip of the modeled slip in 2022 only. Follow-up S1b launched (slip side vs tremor side, component mix,
   footprint depth, mirror in 2022).

## 2026-09-10 S2 depth occupancy
Run: `OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s2_depth_occupancy.py` (script docstring holds the exact
definitions). Outputs: results/s2_yearly.csv (51 rows = 17 years x 3 segs), results/s2_trends.csv (42 rows =
7 quantities x 3 segs x 2 year-sets), results/s2_half_record.csv, figures/s2_depth_occupancy.png.

Setup numbers
- mesh 742 triangles, Slab2 depth band 11.2-57.4 km; segments north=392, central=144, south=206 triangles.
- tremor 748,002 detections -> 741,176 (99.1 %) assigned to a triangle centroid within 40 km; Slab2 depth finite for
  94.1 % of those (the interpolator is only defined inside the 10-60 km mesh band, so the shallowest/deepest 6 % of
  epicentres get no depth and are counted in N but not in the percentiles).
- confirmed patches (thrust & confirmed & seg in {north,central} & Mw>=6.0 & ~broad_pulse): 27 (21 north, 6 central,
  0 south, so all DUR-only columns are NaN for south and sparse for central).
- 128,994 detections fall inside a same-segment confirmed DUR window.
- slip: pos_rate 742 x 7168 days 2007-01-03..2026-08-18, margin total positive slip 386,640 mm-triangle.

Levels (16-yr picture)
- north:   tremor p10/p50/p90 ~ 34.8 / 39.4 / 45.3 km; slip-weighted depth ~ 35.1 km, slip 10-90 % weight ~ 19-50.7 km.
  So the tremor MEDIAN sits ~4 km DOWN-DIP of the slip-weighted depth, but the slip weight is spread over the whole
  mesh depth range (the 13 fixed vbICA patterns are broad) -> the "edges" are a property of the patterns, not a
  resolved slip edge.
- central: tremor p10/p50/p90 ~ 31.8 / 34.3 / 36.8 km; slip-weighted depth ~ 32.9 km (tremor ~1.4 km down-dip).
- south:   tremor p50 ~ 47-49 km vs slip-weighted depth ~ 32.6 km -> a ~15 km offset. This is the known south problem
  (PLAN section 1: the south is at the model's noise floor) plus the FALLBACK mesh geometry; do not read it as physics.

Trends (Theil-Sen vs year, PRIMARY set 2010-2025, permutation p over 2000 year-shuffles)
    seg      quantity      slope km/yr   p_perm   spearman rho (p)
    north    p10            +0.006       0.870    +0.056 (0.84)
    north    p50            +0.003       0.918    +0.021 (0.94)
    north    p90            +0.023       0.458    +0.162 (0.55)
    north    dur_p50        -0.008       0.896    +0.003 (0.99)
    north    slip_depth     -0.008       0.653    -0.035 (0.90)
    north    slip_edge10    +0.012       0.589    +0.159 (0.56)
    north    slip_edge90    -0.022       0.178    -0.347 (0.19)
    central  p10            -0.112       0.104    -0.424 (0.10)
    central  p50            -0.030       0.399    -0.224 (0.41)
    central  p90            +0.005       0.853    -0.012 (0.97)
    central  dur_p50        -0.062       0.355    -0.500 (0.39)   [only 5 years have DUR tremor]
    central  slip_depth     -0.013       0.271    -0.291 (0.27)
    central  slip_edge10    -0.003       0.550    -0.106 (0.70)
    central  slip_edge90    -0.032       0.077    -0.374 (0.15)
    south    p10            -0.255       0.092    -0.482 (0.06)
    south    p50            -0.237       0.039    -0.618 (0.011)
    south    p90            -0.091       0.068    -0.321 (0.23)
    south    slip_depth     -0.020       0.075    -0.385 (0.14)
    south    slip_edge10    -0.009       0.064    -0.435 (0.09)
    south    slip_edge90    -0.003       0.087    -0.400 (0.12)
  21 tests per year-set -> Bonferroni FWER 0.05 needs p < 0.0024. NOTHING passes. The only nominally significant
  entries are in the SOUTH (p50 p_perm 0.039; with the partial 2026 added, p90 p_perm 0.020) and the south is exactly
  where completeness grew most and where the slip model is at its noise floor.
- Adding the partial 2026 changes no verdict (results/s2_trends.csv year_set 2010-2026).

Half-record tremor depth, 2010-2017 vs 2018-2026 (equal-size subsample, 200 draws, median Mann-Whitney p)
    north    n 124,872 vs 151,699  p50 39.4 -> 39.5 km  (dp50 +0.04)  p_full 0.12   p_sub 0.16     -> no change
    central  n  19,727 vs  31,856  p50 34.4 -> 34.4 km  (dp50 +0.03)  p_full 0.015  p_sub 0.026    -> no usable change
    south    n 133,382 vs 235,651  p50 48.6 -> 46.9 km  (dp50 -1.66)  p_full <1e-16 p_sub <1e-16   -> SHALLOWING, but
      it tracks the completeness growth (0.63 in 2010 -> 1.5 in 2021/2024), i.e. new southern stations detecting
      shallower/updip tremor is the leading explanation.

Verdict for S2: the tremor depth band and the modelled slip depth are STATIONARY in the north and central segments over
2010-2026 within the year-to-year scatter (|slope| < 0.03 km/yr north, < 0.11 km/yr central, none surviving FWER).
H2 (secular drift of the depth edges) is NOT supported in the segments where the slip model is trustworthy. The one
apparent southern shallowing (~-0.24 km/yr, ~-1.7 km between halves) is confounded by network growth and by the
southern noise floor -> report as a detection-completeness effect unless a completeness-matched test survives.

Caveats
- Completeness: the "completeness proxy" is only N/segment-mean; it does not correct the depth distribution. Central
  swings 0.12 (2015) to 2.11 (2025) because the segment alternates ETS/no-ETS years, so central yearly percentiles rest
  on as few as 376 detections. Any yearly trend in a low-N segment is fragile.
- Partial 2026: tremor ends 2026-05-27, slip 2026-08-18. 2026 is flagged (`partial` column), greyed in the figure and
  excluded from the PRIMARY trend fit.
- Mesh resolution: 0.25 deg (~28 km) triangles from a FALLBACK Slab2 mesh whose own header warns it "does not reproduce
  Gualandi's inversion"; slip is 13 FIXED patterns, so slip_depth can only move by re-weighting patterns and the
  slip_edge10/90 are pattern properties, not resolved slip edges. Do not read a few-hundred-metre depth change.
- Tremor depth = Slab2 depth at the epicentre (on-interface assumption); PNSN epicentral error is 5-10 km, which in the
  north maps to ~1-3 km of depth, far larger than any fitted trend.

## 2026-09-10 S3 efficiency maps
Run: `OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s3_efficiency_maps.py`. Outputs:
results/s3_efficiency_triangles.csv (742 rows), results/s3_profiles.csv (21 lat bins x 3 epochs),
results/s3_tests.csv, figures/s3_efficiency_maps.png.

Epoch totals and normalisation
    E1 2010-01-01..2014-12-31 (5.00 yr)  N(all mesh) 196,648  N(lat>=44N) 105,724  slip 100,069 mm-tri  1826 d
    E2 2015-01-01..2019-12-31 (5.00 yr)  N(all mesh) 207,822  N(lat>=44N)  89,699  slip  99,645 mm-tri  1826 d
    E3 2020-01-01..2026-05-27 (6.40 yr)  N(all mesh) 336,706  N(lat>=44N) 146,164  slip 124,114 mm-tri  2339 d
    f[e] = tot[e]/mean(tot):  E1 0.9285, E2 0.7878, E3 1.2837   (mean tot 113,862)
    triangles passing N>=50 & S_mm>=5: E1 275, E2 301, E3 323 of 742; 265 have an E3/E1 ratio (186 at lat>=44N).

Margin-wide (lat>=44N) medians over defined triangles
    E1  eff 0.354  eff_norm 0.329 mm/detection
    E2  eff 0.432  eff_norm 0.340
    E3  eff 0.368  eff_norm 0.473

Tests on log10(eff_E3/eff_E1)
- north (lat>=46, n 144, median -0.0195 = x0.96) vs central (44-46, n 42, median -0.1913 = x0.64):
  Mann-Whitney U 4360, p 1.36e-05. Central efficiency drops ~2.6x more than northern.
- margin-wide vs 0 (Wilcoxon, n 186): raw median -0.0710 (x0.849), W 5841, p 1.03e-04;
  normalised median +0.0696 (x1.174), W 5642, p 3.28e-05. Both are "significant" and they point in OPPOSITE
  directions -- the normalisation is a constant log-offset of +0.1406 = log10(f_E3/f_E1), which flips the sign of the
  median. (The one-sample signed-rank test is not shift-invariant, so the two p's differ slightly; the Mann-Whitney
  above IS shift-invariant and is therefore the same for raw and normalised.)
- E2/E1 (n 181): raw median +0.0342 (x1.08), p 0.021; normalised -0.0372 (x0.92), p 0.015 -- same sign flip.

Along-strike profile (median mm slip per detection, raw / normalised; full table in results/s3_profiles.csv)
    lat bin   eff E1 / E2 / E3        eff_norm E1 / E2 / E3
    42.75     1.255 / 0.402 / 0.263   1.166 / 0.317 / 0.338
    43.25     1.135 / 0.440 / 0.214   1.054 / 0.347 / 0.275
    44.25     1.113 / 0.880 / 0.533   1.034 / 0.693 / 0.684
    45.25     0.367 / 0.336 / 0.355   0.341 / 0.265 / 0.456
    46.25     1.027 / 0.920 / 0.735   0.954 / 0.724 / 0.944
    47.25     0.337 / 0.468 / 0.331   0.313 / 0.369 / 0.425
    48.25     0.219 / 0.252 / 0.230   0.203 / 0.198 / 0.295
    49.25     0.216 / 0.399 / 0.278   0.200 / 0.315 / 0.357
  The spatial PATTERN is stable across epochs (peaks near 42.8-44.3 and 46.3, trough 47.5-49) and the 90 % bootstrap
  bands of the three epochs overlap in nearly every bin; the epoch differences are a near-uniform vertical shift.

Verdict for S3 (hits the PLAN S3 decision point)
The margin-wide efficiency "change" IS the network-growth artefact the plan warned about: the raw ratio says efficiency
fell 15 % from E1 to E3, the network-normalised ratio says it rose 17 %, from the same data, both with p < 1e-3 --
the sign of the headline number is decided entirely by the normalisation choice, so no margin-wide change in slip per
tremor detection is claimable. What survives the constant normalisation is the north-vs-central CONTRAST (p 1.4e-5):
central Oregon lost ~36 % of its slip-per-detection between E1 and E3 while the north lost ~4 %. That contrast is
still a candidate detection artefact -- the largest PNSN station additions since 2015 are in Oregon and the central
detection count is the most volatile (see S2: central completeness 0.12-2.11) -- so it needs the completeness-matched
control (subsample detections per triangle to the E1 rate and recompute) before it means anything physical. Not run
here; flagged as the S3 follow-up.

Caveats
- Completeness: eff = S/N is inversely proportional to detection sensitivity. A triangle that gained a nearby station
  gets more N and therefore LOWER efficiency with no change in the fault. The f[e] normalisation removes only the
  margin-average rate change, never a spatially varying one; the north-vs-central result is exactly the case where the
  spatial part matters.
- Epoch length: E3 is 6.4 yr vs 5.0 yr for E1/E2. The ratio S/N is per-detection so length largely cancels, but the
  N>=50 & S>=5 mm gate is easier to pass in E3 (323 vs 275 triangles), so the E3 triangle SET is slightly more
  permissive; the E3/E1 ratio uses only the 265 triangles defined in both.
- Mesh resolution: 0.25 deg (~28 km) triangles, tremor assigned to the nearest centroid within 40 km, so a triangle's N
  mixes detections up to ~20 km away; slip is 13 fixed broad patterns, so per-triangle S_mm is not independent between
  neighbours. Triangle-level bootstraps and the Wilcoxon therefore OVERSTATE the effective sample size.
- Partial 2026: E3 is truncated at the tremor catalogue end 2026-05-27 for both N and S_mm, so the two are consistent,
  but E3 contains only ~5 months of 2026.
- Sign/scale: only positive (thrust-sense) daily rate is accumulated, so S_mm is an upper bound that includes model
  noise rectified into positive slip; this matters most where slip is small (south, inter-ETS periods).

## 2026-09-10 S6 Part A rerun (event-day definitions A1, A2)

Rewrote Part A of `scripts/s6_model_checks.py` after the first run showed that every "tremor-confirmed inter-event slip
burst" was a real ETS episode the CONFIRMED filter had failed to exclude. Part B is unchanged (its RNG stream is
reproduced bit-for-bit by burning the 2 x 200 shift draws the old Part A consumed, so all Part B numbers above still
stand verbatim). Part A now runs under THREE event-day definitions, all kept in the CSVs under a `definition` column:

- **A0** (old, traceability): CONFIRMED patches — thrust & confirmed & seg in (north,central) & Mw>=6.0 & ~broad_pulse,
  n = 27. Event day = within [start-30, end+30] of any such patch touching the band.
- **A1**: SUBSTANTIVE thrust patches — `slip_patches_5mm.csv` with ~minor & ~broad_pulse, **any Mw, any segment**,
  n = **284** (Mw 4.98-6.79, median duration 5 d). Same +-30 d padding. This is the model's own full ETS calendar.
- **A2**: model-independent TREMOR ETS calendar — per band, runs of >=5 consecutive days with the band's 30-d tremor
  count above that band's 90th percentile, padded +-30 d. Band p90 = [1672, 1282, 1224, 520, 136, 290, 319, 1034, 1705,
  1051] detections/30 d; 15-22 ETS runs per band.

Mean event-day coverage over the 10 bands: **A0 15%, A1 51%, A2 28%**.

### Per-band inter-event shares and burst counts

| band | frac positive slip inter-event (A0 / A1 / A2) | frac tremor inter-event (A0 / A1 / A2) | event/inter-event rate ratio (A0 / A1 / A2) | slip bursts (A0 / A1 / A2) |
|---|---|---|---|---|
| 40-41 S | 1.000 / 0.542 / 0.639 | 1.000 / 0.536 / 0.438 | – / 1.36 / 1.25 | 8 / 1 / 3 |
| 41-42 S | 1.000 / 0.349 / 0.634 | 1.000 / 0.343 / 0.381 | – / 1.44 / 1.29 | 20 / 2 / 8 |
| 42-43 S | 0.974 / 0.194 / 0.643 | 0.917 / 0.201 / 0.438 | 0.81 / 1.44 / 1.22 | 25 / 2 / 8 |
| 43-44 S | 0.936 / 0.254 / 0.675 | 0.815 / 0.187 / 0.323 | 1.15 / 1.36 / 1.21 | 8 / 0 / 0 |
| 44-45 C | 0.867 / 0.457 / 0.720 | 0.461 / 0.222 / 0.112 | 1.05 / 1.17 / 1.04 | 1 / 0 / 1 |
| 45-46 C | 0.759 / 0.485 / 0.724 | 0.384 / 0.237 / 0.151 | 1.19 / 1.34 / 1.25 | 2 / 0 / 0 |
| 46-47 N | 0.633 / 0.370 / 0.645 | 0.387 / 0.169 / 0.207 | 1.36 / 1.40 / 1.57 | 2 / 0 / 0 |
| 47-48 N | 0.566 / 0.318 / 0.576 | 0.310 / 0.169 / 0.269 | 1.69 / 1.71 / 1.94 | 1 / 0 / 0 |
| 48-49 N | 0.615 / 0.475 / 0.590 | 0.393 / 0.262 / 0.312 | 1.69 / 1.70 / 1.88 | 3 / 0 / 0 |
| 49-50 N | 0.784 / 0.658 / 0.656 | 0.744 / 0.615 / 0.349 | 1.22 / 1.35 / 1.23 | 2 / 0 / 0 |

Group means — inter-event share of positive slip: A0 0.704 (N/C) / 0.978 (S); **A1 0.461 / 0.335**; **A2 0.652 / 0.648**.
Inter-event share of tremor: A0 0.446 / 0.933; A1 0.279 / 0.317; A2 0.233 / 0.395. Event/inter-event mean positive-rate
ratio stays modest under every definition: A0 1.37 / 0.98, A1 1.45 / 1.40, A2 1.49 / 1.24 — the model never concentrates
its daily slip rate in ETS windows by more than ~1.5x, which is a genuine (definition-independent) property of the vbICA
reconstruction: it smears slip broadly in time.

### Burst confirmation vs the rigid tremor-shift null (200 draws)

| def | group | burst kind | n | observed | null mean | null p95 | p |
|---|---|---|---|---|---|---|---|
| A0 | north/central | slip burst | 11 | 0.909 | 0.299 | 0.545 | 0.005 |
| A0 | north/central | tremor burst | 214 | 0.023 | 0.009 | 0.021 | 0.015 |
| A0 | south | slip burst | 61 | 0.607 | 0.215 | 0.328 | 0.005 |
| A0 | south | tremor burst | 135 | 0.178 | 0.066 | 0.097 | 0.005 |
| **A1** | north/central | slip burst | **0** | – | – | – | – |
| **A1** | north/central | tremor burst | 162 | 0.000 | 0.001 | 0.006 | 1.000 |
| **A1** | south | slip burst | 5 | 0.200 | 0.225 | 0.600 | 0.741 |
| **A1** | south | tremor burst | 48 | 0.021 | 0.017 | 0.037 | 0.418 |
| **A2** | north/central | slip burst | 1 | 0.000 | 0.280 | 1.000 | 1.000 |
| **A2** | north/central | tremor burst | 205 | 0.000 | 0.004 | 0.010 | 1.000 |
| **A2** | south | slip burst | 19 | 0.211 | 0.217 | 0.368 | 0.597 |
| **A2** | south | tremor burst | 94 | 0.043 | 0.039 | 0.069 | 0.393 |

(The A0 southern null shifted very slightly from the first run — 0.209 -> 0.215 and 0.062 -> 0.066 — because the null now
draws ONE shift per iteration shared by all bands, groups and definitions, i.e. it is paired across definitions rather
than drawn independently per group. Observed values and every conclusion are unchanged.)

### Verdict: no tremor-confirmed inter-ETS slip once real ETS days are excluded
Under A0 the picture looked strong — 91% of north/central and 61% of southern inter-event slip bursts were tremor-
confirmed, far above a rigid-shift null (p=0.005 both). It was an artefact of an incomplete ETS calendar. Widening the
calendar to all 284 substantive thrust patches (A1) leaves **zero** north/central inter-event slip bursts at all — every
band-mean >3 mm/30 d episode north of 44N falls inside some substantive patch's +-30 d window — and only 5 in the south,
confirmed at 20% against a null of 22.5% (p=0.74). The model-independent tremor calendar (A2) agrees: 1 north/central
burst (unconfirmed) and 19 southern bursts confirmed at 21% vs a 22% null (p=0.60). The converse test collapses the same
way: inter-event tremor bursts coincide with band slip >3 mm at 0.0-4.3%, indistinguishable from their nulls (p=0.39-1.0)
under both A1 and A2, whereas under A0 they were weakly above null. So the honest statement is that **the model shows no
tremor-confirmed slip between ETS episodes at 1-degree band scale**; what remains between events is a broad, low-level
positive rate (35-65% of the total positive slip, mean rate only ~1.2-1.5x lower than on ETS days) that is spread evenly
rather than concentrated in confirmable bursts, and is exactly what a temporally smooth 13-component reconstruction
would produce whether or not real inter-ETS creep exists. Two caveats keep this from being a clean null result: the
band-mean >3 mm bar is strict for a localised patch, so small inter-ETS transients are not testable this way; and A2 is
partly circular by construction (inter-event days are defined as tremor-quiet days, which caps how tremor-confirmable an
inter-event burst can be) — A1, which is defined purely on the slip side, is the definition to quote, and it gives zero
north/central candidates to begin with.

Files rewritten: `scripts/s6_model_checks.py`; `results/s6_interets_bands.csv` (30 rows, `definition` column),
`results/s6_interets_bursts.csv` (955 rows), NEW `results/s6_interets_conf.csv` (12 rows);
`figures/s6_model_checks.png` now 9 panels (A1a/A1b/A1c per-band fractions and rate ratio for the three definitions,
A2a/A2b burst confirmation vs null, B2a/B2b/B3/B4 unchanged). Part B CSVs unchanged.

## 2026-09-10 orchestrator review — S2 and S3
S2: no secular depth drift of tremor or slip where the model is trustworthy (north |slope| <= 0.02 km/yr, p >= 0.18; central
<= 0.11 km/yr, p >= 0.08; none pass Bonferroni over 21 tests). Levels: north tremor p50 39.4 km vs slip-weighted 35.1 km
(tremor ~4 km down-dip, same as S1's Q1/Q3); central 34.3 vs 32.9. The south's -1.7 km "shallowing" tracks completeness
(0.63 -> 1.5) -> artefact. H2 not supported.
S3: the margin-wide efficiency change flips sign with the network normalisation (raw x0.85, normalised x1.17, both "p<1e-4")
-> PLAN decision point fired: no claimable margin-wide change. Shift-invariant contrast: central Oregon slip-per-detection fell
x0.64 (E3/E1) vs x0.96 north (MWU p 1e-5) — but Oregon received most new PNSN stations, so likely detection sensitivity.
Follow-up S3b launched: efficiency per TREMOR-DAY (days with >= 1 / >= 3 detections) instead of per detection, plus the N/D
(detections per tremor-day) completeness diagnostic per segment. Along-strike efficiency profile is stable across epochs (peaks
42.8-44.3N and 46.3N, trough 47.5-49N) — a persistent spatial pattern of slip-per-tremor, worth keeping.

## 2026-09-10 orchestrator review — S6 Part A rerun
With event days defined from ALL substantive patches (A1) or a tremor ETS calendar (A2), the "tremor-confirmed inter-ETS slip"
of the first pass vanishes completely: zero north/central inter-event slip bursts under A1 (every band-mean > 3 mm/30 d episode
north of 44N lies inside some patch's +-30 d window); southern bursts confirm at the null rate; converse test at null.
VERDICT (quote A1): no tremor-confirmed inter-ETS slip at the 1-deg band scale. Definition-independent: 35-65 % of the model's
positive slip falls outside event windows with an event/inter-event rate ratio <= 1.5 — the vbICA reconstruction smears slip
in time (a model property), but the smeared slip is featureless, not bursts. This bounds what the geodesy resolves between
events: nothing tremor can confirm.

## 2026-09-10 S3b tremor-days efficiency
Follow-up requested after S3, to test whether the central-vs-north efficiency contrast is a detection-sensitivity
artefact. Run: `OMP_NUM_THREADS=4 python tremor_slip_geometry/scripts/s3b_efficiency_tremordays.py`.
Outputs: results/s3b_efficiency_triangles.csv (742), s3b_profiles.csv (63), s3b_tests.csv, s3b_completeness.csv,
s3b_decomposition.csv; figures/s3b_efficiency_tremordays.png.

Idea: replace the detection COUNT N by a count of tremor-ACTIVE DAYS D, which saturates (a day is one day whether the
network reports 3 detections or 30). eff_days = S_mm / D, gate D >= 20 & S >= 5 mm. Four variants:
D_own (>=1 detection in the triangle), D3_own (>=3), D_nb (>=1 in the triangle or a vertex-sharing neighbour;
neighbourhood size median 13 triangles), D3_nb (>=3 in the neighbourhood).

Why NO normalisation here: the S3 factor f[e] = tot[e]/mean(tot) is built from detection counts, the very quantity a
tremor-day is designed to be insensitive to; multiplying eff_days by f[e] would re-inject the growth the definition
removes, and S3 already showed that this constant log-offset alone decides the sign of the margin-wide median. The
residual network effect on tremor-days is a THRESHOLD effect (a weak day becomes detectable), which no scalar can fix;
D3 and the N/D diagnostic are the checks for it instead.

Counts
- 4,605 distinct tremor days on the mesh, 2010-01-06..2026-05-27.
- tremor days on the mesh: E1 1,273 of 1,826 calendar days; E2 1,397 of 1,826; E3 1,935 of 2,339.
- median per-triangle D_own (lat>=44N): E1 3 -> E2 5 -> E3 11. The tremor-day count is itself NOT stationary.
- triangles passing D>=20 & S>=5 mm: D_own 252/298/345, D3_own 194/214/267, D_nb 499/515/550, D3_nb 464/478/503.
- margin-wide (lat>=44N) median eff_days (D_own): E1 3.43 -> E2 3.34 -> E3 2.67 mm per tremor-day.

Tests (raw, no normalisation)
    variant  north (>=46)          central (44-46)        MWU p        margin-wide (lat>=44N) E3/E1  Wilcoxon p
    D_own    -0.1328 (x0.737)      -0.2926 (x0.510)       5.4e-09      -0.1489 (x0.710)              8.4e-24
    D3_own   -0.0609 (x0.869)      -0.1641 (x0.685)       2.9e-04      -0.0783 (x0.835)              1.2e-09
    D_nb     -0.0635 (x0.864)      -0.2885 (x0.515)       8.5e-34      -0.1142 (x0.769)              3.5e-36
    D3_nb    -0.0624 (x0.866)      -0.2669 (x0.541)       7.2e-24      -0.1152 (x0.767)              9.3e-31
  E2/E1 margin-wide is a much smaller drop in every variant (x0.92, x0.96, x0.92, x0.92).
  So in tremor-days the central-vs-north contrast is not only preserved, it is STRONGER and more significant than in
  S3's detection counts (S3: central x0.64 vs north x0.96, p 1.4e-5). Taken at face value that would be a big result.
  It is not one -- see the decomposition.

DECOMPOSITION (the decisive test): log10(eff_days E3/E1) = log10(S E3/E1) - log10(D E3/E1), D_own, triangles defined
in BOTH epochs. Note the epoch-length factor is log10(2339/1826) = +0.1076.
    seg      n     dlog10 S            dlog10 D            -> dlog10 eff_days
    north    137   +0.0944 (x1.243)    +0.2201 (x1.660)    -0.1328 (x0.737)
    central   31   +0.1081 (x1.283)    +0.3917 (x2.464)    -0.2926 (x0.510)
    south     78   +0.0889 (x1.227)    +0.4223 (x2.644)    -0.3444 (x0.452)
    north vs central: dlog10 S  U 1551, p 0.019 |  dlog10 D  U 619, p 7.8e-10 |  dlog10 eff_days U 3551, p 5.4e-09
  The modelled SLIP term is +0.094 to +0.108 in all three segments, i.e. indistinguishable from the pure epoch-length
  factor +0.1076 -> the modelled slip RATE is flat (within a few per cent) everywhere, in every epoch. 100 % of the
  efficiency change, margin-wide AND in the north/central contrast, lives in the tremor-DAY count: per unit time the
  tremor-day rate grew x1.30 in the north, x1.92 in the centre and x2.06 in the south.

Completeness diagnostic (detections per tremor-day, N / D_own; triangles with D>=20 in both E1 and E3)
    north    n 137   N/D 10.53 -> 7.56   median dlog10 -0.1008 (x0.793)
    central  n  31   N/D 10.83 -> 7.36   median dlog10 -0.1590 (x0.693)
    south    n  78   N/D  7.65 -> 8.36   median dlog10 +0.0092 (x1.021)
    north vs central dlog10(N/D): MWU U 2118, p 0.98  -> the two segments' detections-per-day change is the SAME.
    Spearman dlog10(N/D) vs S3 detection-based log10(eff_E3/eff_E1), lat>=44N, n 167:  rho -0.744, p 9.8e-31
    Spearman dlog10(N/D) vs S3b tremor-day log10(eff_E3/eff_E1),      lat>=44N, n 168:  rho -0.130, p 0.093
  Two things follow. (i) The S3 per-triangle pattern WAS largely detection sensitivity: 55 % of its variance is
  explained by a pure catalogue statistic. Tremor-days remove that channel (rho -0.74 -> -0.13). (ii) But N/D FALLING
  (x0.79 north, x0.69 central) while D rises steeply (x1.66, x2.46) is the signature of a LOWERED DETECTION THRESHOLD:
  the network is adding many marginal, sparsely-detected tremor days, not making existing days denser. Consistent with
  this, hardening the day threshold shrinks the effect monotonically -- margin-wide x0.710 (D>=1) -> x0.835 (D>=3),
  and the central/north gap 0.510/0.737 = 0.69 -> 0.685/0.869 = 0.79.

VERDICT: NO -- the central Oregon efficiency drop does not survive as a physical result.
It survives NUMERICALLY in tremor-days (indeed it strengthens: x0.51 central vs x0.74 north, p 5e-9, in all four D
variants), but the decomposition shows the modelled thrust slip per unit time is IDENTICAL in the two segments to
within a few per cent, so nothing about the slip changed; the whole contrast is that central Oregon's tremor-day rate
grew x1.92 while the north's grew x1.30. That differential is itself threshold-limited (falling N/D, and the gap
closing as the day threshold hardens), and central Oregon is where PNSN added the most stations. What S3 + S3b
actually measure is the growth of the PNSN catalogue in Oregon relative to Washington -- a catalogue result, not a
fault result. Efficiency (slip per tremor) is therefore NOT usable as a change detector with this catalogue, in either
detection-count or tremor-day form; drop the H3 "efficiency changed after ~2015" line unless one of these is done
first: (a) recompute on a FIXED station subset operating continuously 2010-2026, (b) count independent ETS EPISODES
per triangle (declustered) instead of days, or (c) use an amplitude-based tremor measure with a fixed detection
threshold. Note also that the same decomposition kills the margin-wide number: the slip term is pure epoch length.

Caveats
- Completeness is the whole story here and is NOT fixed by tremor-days; it is only pushed from the count channel to
  the day-threshold channel (quantified above). D3 is the sensitivity check, not a cure.
- Effective sample size: the neighbourhood variants smooth over a median of 13 triangles, so D_nb/D3_nb triangles are
  strongly correlated and their n (352, 334) hugely overstates the independent sample; the same is true to a lesser
  degree for D_own (28-km triangles, detections assigned up to 40 km away). Read the p-values as ordering, not as
  calibrated evidence.
- Epoch length: E3 is 6.4 yr vs 5.0 yr. S/D is a per-day ratio so length largely cancels, but the D>=20 gate is easier
  to pass in E3 (345 vs 252 triangles for D_own); the E3/E1 ratios use only the 246 triangles defined in both.
- The gate D >= 20 over a 5-yr epoch is ~1 tremor-day per quarter, a low bar; central has only 31 triangles defined in
  both epochs, so the central medians rest on a small, spatially clustered set.
- Slip: only the positive (thrust-sense) daily rate is accumulated, so S_mm rectifies model noise upward; and slip is
  13 fixed broad vbICA patterns on a 0.25-deg FALLBACK mesh, so neighbouring triangles are not independent.
- Partial 2026: E3 is truncated at the tremor catalogue end 2026-05-27 for both tremor and slip.

## 2026-09-10 S4 fronts

Ran `scripts/s4_fronts.py` (PLAN.md step S4, quantity Q6; OMP_NUM_THREADS=4, seed 20260910, ~50 s end to end).
Per-event DAILY tracks — never a footprint average — for the 27 CONFIRMED north/central patches of
`results/s1_events.csv` (21 north / 6 central), because S1 showed the event-averaged along-strike offset Q2 is
erratic under chained union footprints (Q2 IQR 52-125 km).

Definitions actually used (full statement in the script docstring). Window W = [start-45, end+15] d.
Footprint F = triangles labelled with the patch id on any day in [ts, te]; neighbourhood NB = triangles within 60 km
of any F centroid (NB is fixed over W). Slip rate = diff of `-1e3*(Ld @ A_raw.T)` -> 7-d centred mean -> positive
part. SLIP TRACK y_s(t) = area*rate-weighted mean of triangle y over NB, defined only where the NB rate sum >= 20 %
of its W maximum; z_s(t) = the same-weight Slab2 depth. TREMOR TRACK from PNSN detections within 60 km of any F
centroid: daily mean y on days with >= 5 detections (median Slab2 depth for z_t), then 3-d centred median smoothing.
L = lead/lag in -30..+30 d maximising Pearson r between y_s(t) and y_t(t+lag) over days where both are defined
(>= 10 overlapping days required per lag); **L = -argmax so that POSITIVE L means the SLIP track LAGS the tremor
track**. Speeds = Theil-Sen slope of each track vs day index (>= 5 defined days). Depth change = median of the last
third minus the first third of each track's defined days. Mixture diagnostic: per vbICA component the positive
area-weighted rate summed over NB and W, normalised by the sum over the 13 components; n_components = # with >= 10 %.
CONTROL: 50 rigid shifts (>= 365 d, wrapped in 2010-01-01..2026-05-27) of the same neighbourhood's tremor dates,
same track rules, L recomputed against the unchanged slip track.

Scale of the objects measured: F = 23-310 triangles (lat span 1.08-5.75 deg), NB = 102-546 of the 742 triangles,
NB along-strike extent 222-722 km; slip track defined on 20-85 of the ~62-115 window days, tremor track on 19-82.

### Headline numbers (median [q25, q75], combined n=27)

| quantity | combined | north (21) | central (6) |
|---|---|---|---|
| L (d, + = slip lags tremor) | **0.0 [-3.5, +3.0]** | 0.0 [-4.0, +3.0] | +1.0 [-0.8, +2.8] |
| onset lag (d, + = slip starts later) | **0.0 [-2.5, +3.0]** | 0.0 [-2.0, +3.0] | -2.0 [-2.8, +1.8] |
| corr at best lag | 0.690 [0.628, 0.767] | 0.707 | 0.598 |
| corr at lag 0 | 0.493 [0.312, 0.688] | 0.494 | 0.428 |
| v_t (km/d) | +0.13 [-1.33, +2.04] | +0.17 | -0.29 |
| v_s (km/d) | +0.13 [-0.35, +0.54] | +0.13 | -0.12 |
| v_s / v_t | **0.176 [0.070, 0.371]** | 0.176 | 0.310 |
| excursion tremor (km) | 281 [236, 362] (p10-p90 192) | 277 | 323 |
| excursion slip (km) | 139 [110, 165] (p10-p90 81) | 139 | 135 |
| dz_t (km, + = deepening) | -1.55 [-5.28, +0.75] | -1.89 | +1.76 |
| dz_s (km) | -0.28 [-1.22, +0.20] | -0.28 | -0.28 |
| n_components (>= 10 % share) | 3 [2.5, 4] | 3 | 3.5 |

Direction (sign of v) agreement tremor vs slip: **78 % (21/27)**; north 76 %, central 83 %.
`|v_s| < |v_t|` in 24/27 events. Onset lag: 12 positive / 2 zero / 13 negative, median |onset| 3 d.

### CONTROL — two different answers for the two halves of the cross-correlation
- **The peak-correlation VALUE is meaningless.** Observed |r at best lag| median 0.690 vs pooled shifted-tremor
  null 0.652 (n=1104 draws); mean per-event p_null 0.42 (0.5 = no signal, Wilcoxon vs 0.5 p=0.10); only
  **2/27 events (7.4 %)** have p_null < 0.05, i.e. the chance rate. Scanning 61 lags over two smooth,
  strongly autocorrelated tracks with >= 10 overlapping days buys r ~ 0.65 for free.
- **The LAG itself is decisively real.** Observed L is sharply peaked at zero — **63 % of events have |L| <= 5 d
  vs 18.5 % of null draws (binomial p = 4.3e-07)**; median |L| 3 d observed vs 15 d null. North 57 % vs 19 %
  (p = 1.3e-04); central 83 % vs 15 % (p = 3.9e-04). Under the shifted-tremor null L is essentially uniform on
  -30..+30 (figure panel A).

**Reading (H4).** The modeled slip front and the tremor front are locked in time to within a few days:
median lead/lag **0 d**, IQR ±3 d, and 63 % of events inside ±5 d against a flat null. There is **no measurable
lag of slip behind tremor** — PLAN's H4 ("slip front lags tremor front by a constant lag") is not supported in
the direction it predicted; what survives is synchrony, bounded at |L| <~ 3-5 d. And no CHANGE in the lag:
Theil-Sen L vs date is **-0.00 d/yr [-0.55, +0.56], p_perm = 1.00** (north -0.00 [-0.86, +0.67], p = 1.00), so
the "changing lag would be new" branch of H4 is closed at roughly ±0.6 d/yr (±9 d over 16 yr).

**Speeds.** The tremor front migrates at ~0.1-4 km/d (both signs, as expected for bilateral ETS), the slip
centroid at ~1/6 of that (median ratio 0.18), but with the **same sign in 78 % of events** — the model does
reproduce the propagation DIRECTION and not the speed. That is the expected signature of a 13-pattern mixture:
the centroid can only move by re-weighting fixed patterns, so it is heavily damped.

**No trends, no Mw dependence that survives multiplicity.** All 12 date trends (3 groups x L, v_t, v_s, ratio)
have p_perm >= 0.20; combined L -0.00 d/yr p=1.00, v_t -0.042 km/d/yr p=0.75, v_s -0.026 p=0.38, ratio -0.002
p=0.90. Spearman vs Mw: the only nominal hits are north L rho=+0.48 p=0.028 and central ratio rho=-0.83 p=0.042
(n=6); with 12 tests in the block neither survives FWER (Bonferroni alpha 0.0042).

### Caveats
- **MIXTURE (the big one).** Only 2-5 vbICA components (median 3) carry >= 10 % of the neighbourhood's positive
  slip in a window, and the top component carries a median 25 %. The slip track is therefore a re-weighting of
  ~3 fixed spatial patterns, which is why v_s is ~6x slower than v_t and why excursion_s (139 km) is half
  excursion_t (281 km). The **direction agreement and the synchrony are the defensible results; the slip SPEED
  is not a physical propagation speed** and the v_s/v_t ratio is a model-damping number, not a rheology number.
- **The neighbourhood is not local.** F is a union over [ts, te] and chains along strike, so NB spans 222-722 km
  (546/742 triangles for event 1317). "Along-strike track" for the biggest events is a margin-wide centroid, not a
  front position; that is why the raw excursions reach 400-710 km — single sparse days at the far end of the
  footprint. The p10-p90 excursions (192 km tremor / 81 km slip) are the robust version. Because both tracks live
  on the same wide neighbourhood, part of the r ~ 0.65-0.7 is shared margin-wide structure — exactly what the
  shifted-tremor null shows.
- **Days with tracks.** The slip track exists on 20-85 window days and the tremor track on 19-82; the median
  overlap at the best lag is 24 d (minimum 10, the imposed floor). Event 1434 (2023-11-11, N=75 in S1, 20 track
  days) has a negative best correlation (-0.11), p_null = 1.0, and only 9/50 null draws admissible — treat its
  L = -2 d as noise. Events 1390 (19/50) and 1655 (28/50) also have thin nulls.
- **Network completeness.** The >= 5-detections/day rule is an absolute threshold on a catalogue whose detection
  rate grows with the PNSN network (S1: N per event 75-10,880). Early, sparse events get fewer defined tremor days
  and hence noisier tracks; this cannot bias L (the null is computed per event with the same rule), but it does
  mean the *number* of track days is not comparable across epochs, and it is one more reason the date trends
  above should be read as bounds, not as measurements.
- **Along-strike = y (north).** Per spec, no strike rotation. Fine inside one event; do not compare v across the
  margin (the strike rotates from ~N in Oregon to ~NW off Vancouver Island).
- **Sign of L** is a convention: L = -argmax_lag corr(y_s(t), y_t(t+lag)), so +L = slip lags tremor. Both L_days
  and the raw lag_argmax are in the CSV.
- Event 1625 has v_t = -3e-7 km/d (a genuinely flat tremor track); it counts as a direction DISAGREEMENT and its
  v_s/v_t is left NaN.
- The 6.4 % of detections with no Slab2 depth are dropped from z_t only (they still count toward the >= 5/day rule).

Files written: `scripts/s4_fronts.py`, `results/s4_fronts_events.csv` (27 rows x 53 cols),
`results/s4_fronts_summary.csv` (102 rows), `figures/s4_fronts.png` (8 largest-Mw events' tracks + panel A
L histogram vs null, panel B v_s vs v_t with 1:1, panel C L vs date). Nothing outside `tremor_slip_geometry/`
was created or modified.

## 2026-09-10 orchestrator review — S3b
Decomposition log(eff) = log(S) - log(D) settles H3: the modeled slip term is the pure epoch-length factor in EVERY segment (slip
rate flat 2010-2026), while tremor-days per unit time grew x1.30 north, x1.92 central, x2.06 south, with the growth threshold-
limited (N/D fell as D rose; hardening the threshold shrinks the effect monotonically). The "central Oregon efficiency drop" is
PNSN catalogue growth in Oregon relative to Washington. H3 CLOSED at this data level; a fixed continuously-operating station
subset, declustered episodes per triangle, or an amplitude-thresholded measure would be needed to reopen it. What stands: the
along-strike efficiency PATTERN (persistent peaks 42.8-44.3N, 46.3N; trough 47.5-49N) and the stationarity of potency per
detection in confirmed events (S1) — both relative statements insensitive to the margin-wide completeness factor... the S1
statement only if the confirmed events' detection rate is itself stable (north grew x1.3 in tremor-days) — treat as "stable to
within a factor ~1.3".

## 2026-09-10 S1b 2022 outliers
Coordinator follow-up on S1: the three northern 2022 events (id 2015 2022-03-31, 1583 2022-05-26, 1905 2022-12-02) have
Q1 = +21.1 / +19.6 / +28.9 km and Q3 = +3.5 / +4.6 / +5.4 km against ~0-10 / 0-2 km everywhere else.
Ran `scripts/s1b_2022_outliers.py` (OMP_NUM_THREADS=4, ~43 s). It recomputes S1's Q1/Q3 from scratch for all 21 northern
confirmed events (asserts they reproduce s1_events.csv exactly) and adds: an area-only footprint centroid/depth, the
slip-inside-footprint offset, the tremor aperture sweep (20/40/60 km), the vbICA component decomposition of the
accumulated slip, a per-event margin control, a window-phase control, and a chained-footprint check.

### (1) Slip side or tremor side? -- group medians (north, n=21)
    group        n     Q1     Q3  Q3_fp  slip_depth  fp_depth  tremor_depth  slip-fp  tremor-fp  n_tri     N
    2010-2019   12   6.98   1.25   1.22       37.67     37.75         38.62     2.10      11.14    147  5597
    2020-2021    2  -0.17   0.50  -0.05       38.02     38.57         38.52     2.51       2.34     97  6220
    2022         3  21.06   4.61   4.08       36.16     39.37         40.76    -0.08      24.59     60  1388
    2023-2025    4   1.99   0.61   0.25       37.62     38.01         38.26     3.88       5.36    141  7881
Per event: 2015 slip 34.8 / fp 34.2 / tremor 38.3 km (fp is the SHALLOWEST of the 21; tremor normal);
1583 slip 36.2 / fp 39.4 / tremor 40.8 (slip sits 12.5 km UP-DIP inside its own footprint -- the most negative
slip-minus-footprint of the 21); 1905 slip 39.8 / fp 39.8 / tremor 45.2 (footprint normal, TREMOR very deep).
So the three do NOT share one mechanism: 2015 is footprint/slip-side, 1583 is slip-concentration-side, 1905 is
tremor-depth-side.

### (2) Component mix -- did the model switch to a different mix in 2022?
Fixed component depths are all within 30-35 km (k0..k12 = 32.4 34.0 33.2 34.3 31.8 30.7 34.0 34.0 34.9 30.2 32.1 32.4
33.9), so `mix_depth_km` barely varies (33.5-34.3) and is a WEAK diagnostic -- the components are margin-wide patterns,
not depth-specific ones. The mix IDENTITY does change, though. 18 of 21 events are k1/k2-dominated; the exceptions are
    2015  k8:+0.62 k12:+0.26   (the only k8-dominated event)
    1583  k7:+0.47 k1:+0.29 k6:+0.28
    1935  k8:+0.42 k2:+0.25 ...
and 1905 is ordinary (k1:+1.14, with compensating negatives). Test: frac_k1k2 = frac_1+frac_2.
    unusual mix (k1+k2 < 0.5): n=3, ids [1583, 1935, 2015] -> Q1 med 19.6 km, Q3 med 3.5 km
    usual   mix (k1+k2 >= 0.5): n=18                      -> Q1 med  4.1 km, Q3 med 1.0 km
    Spearman frac_k1k2 vs Q1 rho=-0.47 p=0.031; vs Q3 rho=-0.29 p=0.20
So an unusual component mix does go with a large Q1 -- but it covers only 2 of the 3 (1935 has Q1=8.9, and 1905 has an
ordinary mix), so it is a contributing factor, not the explanation.

### (3) Footprint up-dip test
2022 fp area-weighted depth 39.4 km vs tremor 40.8 -> Q3_fp = +4.08 km, still the largest group (others 1.22 / -0.05 /
0.25). Measuring against the FOOTPRINT instead of the slip centroid removes only ~0.5 km of the 4.6 km Q3, so the
anomaly is not simply "the slip centroid moved inside a normal footprint". For 2015 specifically the footprint IS the
problem (fp 34.2 km, the shallowest of the record, Q3_fp +4.1 vs Q3 +3.5).

### (4) Mirror control in 2022 -- the reverse of an artefact
    mirror north, N>0, by year: 2019 Q1 22.1 / 2020 31.2 / 2021 46.3 / **2022 3.0** / 2024 17.2 / 2025 40.6
    mirror both segments, N>0: 2022 n=4 Q1 med 3.04 Q3 med 1.28  vs  other years n=36 Q1 med 18.61 Q3 med 5.43
2022 is the most NORMAL mirror year on record. Whatever produces the thrust-side 2022 outliers is not a model-wide
2022 pathology -- it does not touch the back-slip catalog.

### (4b) The dominant control: the DEPTH OF ALL NORTHERN TREMOR IN THE WINDOW
    Spearman (n=21):  north-wide tremor depth in window vs Q3  rho=+0.72  p=0.0002   <-- strongest predictor found
                      north-wide tremor depth in window vs Q1  rho=+0.57  p=0.0071
                      northern tremor rate/day in window  vs Q1 rho=-0.62  p=0.0027
Per-event margin control (median Slab2 depth of ALL lat>=46N tremor in the SAME window, no spatial restriction):
    2015  near-patch 38.31  north-wide 38.63  delta -0.32 km
    1583  near-patch 40.76  north-wide 41.23  delta -0.47 km
    1905  near-patch 45.23  north-wide 45.17  delta +0.07 km
group medians of delta: 2010-2019 -0.01, 2020-2021 +0.17, 2022 -0.32, 2023-2025 +0.00.
**The tremor near each 2022 patch is at exactly the depth ALL northern tremor had during that window.** Nothing local
happened. 1905's window (2022-12-02..12-05, only 4 days) is simply the deepest northern tremor window in the record
(45.2 km vs a typical 38-39 km); 1583's is 41.2 km. Since the modeled slip depth is pinned near 36-39 km by the 13 fixed
patterns, Q3 largely measures "how deep was the northern tremor band on the days the 5 mm/30 d threshold selected".

### (4c) Not chained footprints
2022 footprint along-strike extents are 130 / 231 / 120 km -- among the SHORTEST. The long chained footprints
(1317 638 km, 1588 361 km, 1605 342 km) all have small Q1; Spearman extent vs Q1 rho=-0.59 p=0.005. Chaining averages
Q1 DOWN, it does not inflate it. (Hypothesis raised from the 1583 map and rejected.)

### (5) Aperture and size controls
Aperture: Q1 is flat from 20 to 60 km for all three (2015 21.02/21.06/21.06; 1583 19.27/19.57/19.65;
1905 27.57/28.90/28.97). NOT a measurement-aperture artefact.
Size: Q1 vs n_tri_F rho=-0.59 p=0.005; vs N rho=-0.69 p=0.0005; vs potency rho=-0.54; vs Mw rho=-0.53. The 2022 three are
among the 6 smallest northern events (37-88 triangles, N 571-1899). But size is not sufficient: among the 10 small
(<=90 triangle) northern events, 2022 still has Q1 med 21.1 km vs 10.2 km for the others (Q3 4.6 vs 1.9).
Year-wide tremor: northern (lat>=46N) tremor depth p50 = 39.72 (2019-21) / 40.00 (2022) / 39.09 (2023-25); p10/p90
35.0/45.9 in 2022 vs 35.1/45.2 and 34.6/45.4. Mann-Whitney p=3e-14 and 1e-80 purely from n~1.6e4-6.5e4 -- the effect is
0.3-0.9 km, an order below the 3.5-5.4 km Q3 of the outliers. **2022 as a YEAR is not deeper in any useful sense.**
(Also: the 40-km band around the union of the 21 northern footprints contains 100.0 % of all lat>=46N detections, so the
"restricted to the event band" subset is identical to "all northern" -- the restriction is vacuous, keep that in mind
for S2.)

### VERDICT
**Both, but neither is a genuine tremor-vs-slip geometry change.** Ranked:
1. **Window-phase / event-definition effect (dominant, tremor-side but NOT local).** Q3 is mostly set by how deep the
   entire northern tremor band was during the days the 5 mm/30 d threshold happened to select (rho +0.72). The
   near-patch tremor is within 0.5 km of the north-wide tremor for all three events, so there is no local down-dip
   migration. 1905 is the extreme case: a 4-day December window that caught the deepest northern tremor of the record.
2. **Model/footprint effect (secondary, slip-side).** 2015's footprint is the shallowest of the 21 (34.2 km) and is the
   only k8-dominated event; 1583's slip concentrates 12.5 km up-dip inside its own footprint and is one of only three
   events not built from k1/k2. An unusual component mix goes with Q1 med 19.6 vs 4.1 km.
3. **Small-event amplifier.** All three are small; Q1/Q3 rise systematically as events shrink (rho -0.59 to -0.69).
4. **Ruled out:** measurement aperture, a year-wide 2022 tremor deepening, a model-wide 2022 artefact (the mirror
   catalog's most normal year), and chained footprints.
No change to the S1 conclusions: these are three separate flukes of the event definition, not evidence of a 2022
down-dip shift, and S1's trend tests already found nothing (Q1 p=0.64, Q3 p=0.58 -- with these three included).

### Consequence for S2 (act on this)
Depth quantities must be **referenced to the concurrent north-wide tremor depth** (or restricted to matched window
phase), otherwise Q3 and the Q4 edges are largely reporting which days the slip catalog picked, not tremor-vs-slip
geometry. Recommend carrying `delta_vs_northwide` (already in the CSV) as the S2 depth quantity alongside raw Q3.

Files: `results/s1b_2022_outliers.csv` (21 event rows + tremor_epoch/tremor_year/tremor_mw/mirror_year blocks, `block`
column), `figures/s1b_2022_outliers.png` (3 event maps with footprint slip, DUR tremor, three centroids and Slab2
contours; plus Q1 vs footprint size, northern tremor depth by year, mirror Q1/Q3 by year). Nothing outside
`tremor_slip_geometry/` was created or modified.

## 2026-09-10 orchestrator review — S4
Daily slip-rate track vs daily tremor track (27 confirmed events): SYNCHRONOUS — |lag| <= 5 d in 63 % of events vs 18.5 % under
the shifted-tremor null (binomial p 4e-7; north p 1e-4, central p 4e-4), median |L| 3 d vs 15 d null; no lag drift with time
(-0.00 d/yr, bound +-0.6 d/yr). The peak-correlation VALUE is meaningless (null r 0.65 from scanning 61 lags over autocorrelated
tracks) — only the lag distribution carries information; good catch by the agent. So H4 (slip lags tremor) is NOT supported:
the vbICA daily amplitudes move with the tremor front to within a few days; the 35-d "tremor leads slip" of the earlier catalog
comparison was the trailing-window threshold crossing, not model latency. Propagation direction agrees in 78 % (21/27). The slip
centroid speed is ~1/6 of the tremor speed and its excursion half — the 13-pattern mixture cannot follow a front (median 3
components carry >= 10 % of an event's slip): v_s is not physical. Neighbourhoods span 220-720 km because union footprints
chain along strike — S4's statistics are margin-scale synchrony, not local front tracking.

## 2026-09-10 orchestrator review — S1b (2022 outliers)
Not a geometry change. (1) Window selection dominates: Q3 tracks the depth of ALL northern tremor during the event window
(rho 0.72, p 2e-4); relative to concurrent north-wide tremor the three 2022 events are -0.3/-0.5/+0.1 km, i.e. normal. The
5 mm/30 d threshold picked short windows (one is 4 days) that coincided with the deepest northern tremor days in the record,
while the modeled slip depth is pinned at 36-39 km by the fixed patterns. (2) Secondary slip-side effect: the 2022 events are
the only ones not dominated by components k1/k2 (k8 for #2015); events with frac_k1k2 < 0.5 have Q1 19.6 km vs 4.1 km (p 0.03),
and #1583 concentrates slip 12.5 km up-dip inside its footprint. (3) Small-event amplifier (Q1 vs N rho -0.69). Ruled out:
aperture, a year-wide 2022 tremor deepening (0.3-0.9 km only), a model-wide 2022 artefact (mirror 2022 is its most normal year),
chaining. METHOD LESSON: Q3 and the Q4 edges must be referenced to the concurrent north-wide tremor depth; the 40-km band
around the union of northern footprints contains 100 % of northern detections, so "near the patch" is not a local filter in the
north. S1's null trends are unchanged.

## 2026-09-10 S5 site recurrence

Ran `scripts/s5_site_recurrence.py` (PLAN.md step S5) end to end, ~2 min 14 s, OMP_NUM_THREADS=4.
SITES = 20 bands of 0.5 deg latitude, 40.0-50.0; window 2010-01-06..2026-05-27 (N = 5986 d).
Tremor -> band by the detection's OWN epicentre latitude, keeping only epicentres within 40 km of a
mesh-triangle centroid (cKDTree, km frame x=(lon+124)*111*cos47, y=(lat-47)*111): 728,548 of 748,002
detections kept (97.4 %). TREMOR EPISODE = local maximum of the 30-d trailing count >= the band's
90th percentile, greedy highest-first with >= 90 d separation; onset = first day in [peak-45, peak]
with count > 25 % of the peak. Peaks only sought from day 30 of the grid (partial trailing window).
SLIP = substantive patches (~minor & ~broad_pulse, any Mw), start inside the window: 239 of 284
(126 with Mw >= 6.0); footprint = triangles with label == id at any time; a patch enters every band
its footprint touches -> 1041 (band, patch) instances, 292 of them from CONFIRMED thrust patches.
MATCH: episode -> the touching patch whose start is closest to the onset within [onset-30, onset+90];
patch -> the episode whose onset is closest to the start within [start-90, start+30].

### Headline numbers
- 359 tremor episodes (N 139 / C 68 / S 152); 501 episode-patch pairs, 283 of them mutual.
- Episodes matched by a slip patch: **0.791** (284/359). Patch instances matched by tremor: **0.480**
  (500/1041); Mw >= 6.0 instances **0.532** (355/667); confirmed instances 214/292 = 0.733.
- Rigid tremor-shift null (1000 draws, circular shift 365..5621 d, episodes recomputed each draw):
  episode-matched 0.791 vs null mean 0.609 [p5 0.561, p95 0.659], **p = 0.001**; patch-matched 0.480
  vs 0.357 [0.318, 0.393], **p = 0.001**. Null episode count 359.5/draw vs 359 observed, so the shift
  preserves burstiness. Per band, 10/20 beat the null at raw p < 0.05, and **3 survive Holm**:
  46.0 (1.000 vs 0.676), 47.0 (1.000 vs 0.693), 47.5 (0.941 vs 0.649) — all north.
  So the association is real but MODEST: chance matching alone already reaches ~0.61 because the
  patch catalogue is dense and the window is 120 d wide.
- RECURRENCE, tremor: median interval 203-444 d across bands (segment medians N 349 d, C 298 d,
  S 294 d) — the familiar ~10-14-month ETS clock, longest in the north.
  RECURRENCE, slip patches: much shorter, segment medians N 111 d, C 125 d, S 71 d
  (Mw >= 6.0: 171 / 228 / 101 d). In 19 of 20 bands the model emits patches FASTER than tremor
  recurs; the excess is largest in the south, i.e. the model's noise-floor region (consistent with
  the S6 finding that the southern +-5 mm signal is a coherent model mode, not slip).
- PHASE (lag = slip start - tremor onset), 284 matched pairs: median **+5 d**, IQR -3..+14 d.
  Segment medians N +1 d, C +9.5 d, S +8 d; confirmed-only median +2 d. Slip start is at or just
  after tremor onset, as H4 expects — but see the caveat on the trailing frame below.
- PHASE DRIFT: pooled Theil-Sen slope +0.0000 d/yr; tie-robust Spearman rho = +0.031, p = 0.607.
  Per band, no lag trend reaches raw p < 0.05 by either statistic. **No detectable phase drift.**

### Trends and the network-growth control
Raw p < 0.05 (Holm over the 20 bands, per family): tremor recurrence 41.0 (-14 d/yr, p=0.039),
42.5 (-41, p=0.005), 43.0 (-31, p=0.001, Holm 0.020 SURVIVES), 44.5 (-32, p=0.002, Holm 0.038
SURVIVES); slip recurrence 42.5 (+1.4, p=0.023) and 43.5 (+3.7, p=0.039), neither surviving Holm.
All four raw-significant tremor trends are NEGATIVE (intervals shortening) and sit in the bands with
the fastest detection growth (+2 to +13 %/yr), which is the signature of the artefact: the 90th-pct
peak gate is applied to the WHOLE record, so a band whose rate grows puts most of its peaks late.
Pre-registered control (PLAN.md sec. 4): divide each band's daily count by a 365-d centred running
mean of itself (floored at 10 % of the band's median rate), redo episodes and the trend.
**All four die**: 41.0 -14 -> +3.3 (p=0.52); 42.5 -41 -> +9.9 (p=0.030, SIGN FLIPPED); 43.0 -31 ->
+4.7 (p=0.47); 44.5 -32 -> -2.7 (p=0.68). No band survives Holm in the rate-normalised family
(min Holm p = 0.50). **VERDICT: no site's recurrence interval has a demonstrable trend over
2010-2026; the apparent shortening is detection growth.** Panel (d) is annotated to say so.

### Yearly
`results/s5_yearly_match.csv` carries per (year, segment) the matched fractions both ways, the
counts, the yearly kept detections (context) and partial = True for 2026. Tremor->slip matched
fraction has no trend in the north (Spearman rho -0.11, p=0.70) or centre (+0.04, p=0.87); the
south declines (rho -0.55, p=0.03 raw; p ~ 0.09 Bonferroni over the 3 segments) from 1.00 in
2010-2013 to 0.69-0.88 after 2020 — most plausibly because network growth adds southern episodes
faster than the model adds patches, not a physical change. 2026 (3 N, 2 C, 4 S episodes) is partial
and its 0.00 northern/central fractions are an end-of-record truncation, not a signal: a patch that
would match a 2026 episode may start after 2026-05-27.

### Cross-check against the existing table
`res_catalog_sse/tremor_episodes.csv`: 348 episodes; S5 finds 359. Same band and peak within 20 d:
**335 of 359 S5 episodes (93.3 %) have an old counterpart, 335 of 348 old ones (96.3 %) have an S5
counterpart.** The two detectors agree on the ETS calendar; the small excess is S5's peaks-from-day-30
rule plus the greedy 90-d de-clustering.

### Caveats
1. TRAILING FRAME. Slip start/end are 30-d trailing threshold crossings, so every slip time is late
   by 0-30 d relative to true slip onset. The +5 d median lag is therefore NOT evidence that slip
   follows tremor; corrected for the frame it is consistent with slip leading tremor by up to ~25 d.
   Only the ABSENCE of a drift in that lag is frame-independent (the bias is constant in time).
2. CHANCE FLOOR. With 239 patches spread over 1041 band-instances and a 120-d matching window, the
   shift null already matches 61 % of episodes. Read 0.79 against 0.61, not against 0.
3. UNIT. Yearly and unmatched fractions count (band, item) INSTANCES, so a patch touching 5 bands
   counts 5 times. Distinct-patch counts are in the same CSV (n_patch_uniq).
4. TIE ATOM. lag_d is a tied integer, so the Theil-Sen pairwise-slope set has a large atom at 0;
   when the median lands there the slope is exactly 0 and its permutation p is 1 BY CONSTRUCTION.
   Spearman is reported alongside for every lag test and is the statistic to read.
5. Band 49.5 has only 2 substantive patches (1 interval, 3070 d) — its slip recurrence is not a
   statistic. Band 49.5 also has the worst episode-matched fraction (0.176), consistent with the
   model mesh thinning out at the north end.
6. The south's slip "recurrence" of ~71 d is dominated by the model's noise-floor patches (S6);
   it should not be read as a physical recurrence.

Files: `scripts/s5_site_recurrence.py`; `results/s5_tremor_episodes.csv` (359),
`results/s5_matches.csv` (501), `results/s5_site_recurrence.csv` (20 bands x 42 cols),
`results/s5_yearly_match.csv` (51 year-segment rows); `figures/s5_recurrence.png`
(a: per-site median recurrence tremor vs slip, log x; b: lag vs date by segment with Theil-Sen +
Spearman; c: yearly matched fractions with the shift-null chance level; d: north/central recurrence
intervals vs date, annotated with the artefact verdict).

## 2026-09-10 orchestrator review — S5
Association tremor episode <-> slip patch per 0.5-deg band is real but modest against chance: 79 % of episodes matched vs 61 %
under the rigid-shift null (p 0.001); 48 % of patch instances matched vs 36 % (confirmed patches 73 %); only 3 bands (46.0,
47.0, 47.5N) survive Holm. Slip patches recur 2-4x faster than tremor episodes in 19/20 bands (N 111 vs 349 d; S 71 vs 294 d)
— catalog fragmentation + southern noise, not physics. NO phase drift of slip start vs tremor onset (0.00 d/yr, Spearman 0.03,
p 0.61). The apparent shortening of recurrence in bands 43.0/44.5 (-31 d/yr, Holm-significant) disappears under the pre-
registered rate-normalisation control (p 0.47/0.68) -> detection growth again. The +5 d median lag is frame-limited (trailing
window) and not evidence of slip following tremor; S4 gives the frame-free answer (synchronous). STUDY COMPLETE for S1-S6.

## 2026-09-10 round 2 launched (user: "don't be lazy" — answer Q5, Q7, Q2 properly; add mechanism tests Q13, Q14, Q17, Q18)
Agents S7, S8, S9 launched with specs (PLAN.md round 2).

## 2026-09-10 S7 occupancy + terminations

Ran `scripts/s7_occupancy_terminations.py` (PLAN.md round 2, S7; quantities Q5 and Q17) end to end,
`OMP_NUM_THREADS=4`, seed 20260910, ~37 s. Window 2010-01-06..2026-05-27; epochs E1 2010-01-06..2014-12-31
(4.99 yr), E2 2015..2019 (5.00 yr), E3 2020..2026-05-27 (6.40 yr). 742 triangles, 741,176 of 748,002 PNSN
detections assigned to a triangle within 40 km, 4,605 distinct tremor days on the mesh.
Outputs: `results/s7_occupancy_triangles.csv` (742), `s7_occupancy_summary.csv` (87),
`s7_quiet_slip_candidates.csv` (247), `s7_terminations_events.csv` (27), `s7_termination_bins.csv` (32);
`figures/s7_occupancy_maps.png`, `figures/s7_terminations.png`.

### The completeness trap and how S7 dodges it
S3b's decomposition is the binding constraint: per unit time the modelled slip rate is FLAT everywhere
(dlog10 S = +0.089..+0.108 vs the pure epoch-length factor +0.1076 in all three segments) while the
tremor-day rate grew x1.30 north / x1.92 central / x2.06 south, and detections-per-tremor-day FELL
(x0.79 north, x0.69 central) — a lowered detection threshold. So no absolute tremor level is comparable
across epochs. S7 therefore never compares levels: within each epoch AND each segment, S (mm/yr of
accumulated positive thrust slip rate) and D (tremor-days/yr, the D3_nb variant = >=3 detections in the
triangle or a vertex-sharing neighbour, S3b's least sensitivity-dependent one) are converted to
PERCENTILE RANKS. A segment-wide monotone sensitivity change leaves ranks untouched by construction; only
a change in the SPATIAL pattern can move a class. Classes: SLIP+TREMOR (rS>=50 & rD>=50), SLIP-ONLY
(rS>=50 & rD<25), TREMOR-ONLY (rD>=50 & rS<25), QUIET (both <25), MIXED otherwise.

### Part A numbers
- Class composition is nearly fixed by construction (it is a rank rule), so the content is the SPATIAL map.
  North: SLIP+TREMOR 39.8/39.8/40.8 %, SLIP-ONLY 0.0/0.0/0.0 %, TREMOR-ONLY 2.8/1.8/2.3 %, QUIET
  17.3/17.1/17.6 %, MIXED 40.1/41.3/39.3 % (E1/E2/E3; area fractions within 0.2 pt of the count fractions
  everywhere, so no class is an artefact of triangle size). Central: SLIP-ONLY 6.2/3.5/2.1 %, TREMOR-ONLY
  7.6/6.2/6.2 %, QUIET 11.8/13.2/15.3 %. South: SLIP-ONLY 3.9/3.9/4.4 %, TREMOR-ONLY 0.5/1.0/0.5 %.
- CO-LOCATION (the reason SLIP-ONLY is essentially empty in the north): Spearman(rank_S, rank_D) per epoch
  is +0.71/+0.77/+0.76 north, +0.70/+0.73/+0.72 south, but only +0.40/+0.46/+0.44 CENTRAL. Modelled slip
  and tremor occupancy are strongly co-located in Washington/Vancouver Island and markedly less so in
  central Oregon.
- REPRODUCIBILITY of the underlying ranks E1 vs E3: rank_S +0.99 (all segments), rank_D +0.94/+0.93/+0.89.
  This is the honest reading of the persistence test below — the ranks barely move.
- PERSISTENCE (same class in all three epochs) vs 5,000 independent within-segment label shuffles:
  north 86.5 % vs null 13.5 % [10.7, 16.1], p 0.0002; central 71.5 % vs 11.8 %, p 0.0002; south 86.4 % vs
  13.5 %, p 0.0002; lat>=44N pooled 82.5 % vs 13.0 %, p 0.0002. Per class (lat>=44N): SLIP+TREMOR 186 tri
  (34.7 %, depth 39.3 +- 6.3 km), TREMOR-ONLY 15 tri (2.8 %, 38.0 +- 6.2 km), QUIET 70 tri (13.1 %,
  16.8 +- 3.3 km), MIXED 171 tri, SLIP-ONLY 0 tri. All p 0.0002 except SLIP-ONLY (p 1.0).
  CAVEAT: this null is very weak — it destroys the spatial pairing of labels that are drawn from
  quantities correlating at rho 0.89-0.99 between epochs, so a large excess over the shuffle was
  essentially guaranteed. Read the persistence as descriptive bookkeeping of a stable map, not as
  evidence of anything surprising.
- ALONG-DIP structure (Slab2 depth of each class, pooled over epochs, mean [median]):
    north    SLIP+TREMOR 39.9 [38.7] | TREMOR-ONLY 45.3 [45.4] | QUIET 18.4 [17.0] | MIXED 30.3 [27.1]
    central  SLIP+TREMOR 35.8 [35.3] | SLIP-ONLY 37.9 [44.5] | TREMOR-ONLY 33.6 [35.4] | QUIET 17.3 [16.1]
    south    SLIP+TREMOR 39.1 [39.1] | SLIP-ONLY 17.3 [18.6] | TREMOR-ONLY 27.2 [27.2] | QUIET 15.8 [15.6]
  In the north the persistent TREMOR-ONLY patches sit ~4-5 km DOWN-DIP of the SLIP+TREMOR band; QUIET is
  the shallow up-dip locked zone at 15-18 km in every segment.
- PERSISTENT SLIP-ONLY: 7 triangles, ALL in the far south (40.33-41.42N, -124.4..-124.1, depth 13-20 km,
  up-dip of that segment's tremor band), 5 of the 7 on the mesh BOUNDARY, median 0.83 deg from the mesh
  latitude limit. Zero in the north, zero in the central segment. This is exactly the Mendocino / model
  noise-floor corner S6 flagged, at the shallow mesh edge -> NOT credible silent slip.
- PERSISTENT TREMOR-ONLY: 15 triangles in two clusters.
    (a) central Oregon 44.33-45.42N, -123.9..-123.1, depth 27-38 km, 9 triangles, NONE on the mesh
        boundary, 4.5-4.8 deg from the mesh latitude limits (tri 234,235,236,248,249,254,255,286,304).
    (b) north Vancouver Island 49.33-49.67N, -125.2..-124.4, depth 42-47 km, 6 triangles, 4 of them ON the
        mesh boundary, 0.25-0.58 deg from the mesh north limit (tri 698,699,724,725,727,729).
  Cluster (a) is the only one that is not an edge effect.
- ABSOLUTE map: the literal thresholds in the spec are DEGENERATE. Because only the POSITIVE daily rate is
  accumulated, every triangle rectifies model noise upward: S_tot spans 123..1284 mm (median 407 mm) over
  16.4 yr, so "S_tot >= 20 mm" passes 742/742 and "S_tot < 5 mm" passes 0/742. The literal quiet-slip
  criterion therefore reduces to a pure tremor criterion (231 triangles = simply D3_nb_tot < 30 d: north
  122, central 48, south 61; 219 of them up-dip of their segment's tremor band, 12 down-dip) and the
  literal converse is the empty set. Percentile versions relative to each segment's own S_tot are the ones
  to read (segment p25/p75 in mm: north 258/522, central 292/469, south 363/701):
    quiet-slip_p75 (S_tot >= segment p75 AND D3_nb_tot < 30 d): 9 triangles — 6 in central Oregon at
      44.08-44.67N, -122.7..-122.6, depth 44-54 km, all DOWN-DIP of the central tremor band (p90 37.0 km);
      1 north (46.33N, -122.33, 54.1 km, down-dip); 2 south (40.33/40.58N, 14-19 km, up-dip).
    tremor-no-slip_p25 (D3_nb_tot >= 300 d AND S_tot <= segment p25): 7 triangles — 3 north at
      49.33-49.58N (45 km), 2 central at 44.33/44.42N (37-38 km), 2 south at 41.17/43.67N (26-28 km).
  Segment tremor Slab2-depth bands (p10/p50/p90 over the detections): north 34.8/39.4/45.4 km,
  central 31.5/34.4/37.0 km, south 33.0/47.7/52.0 km.

### Part B numbers (27 confirmed patches, 21 north / 6 central — the S1/S4 set)
- Extents: tremor p5-p95 latitude span median 1.69 deg [0.37, 5.05]; modelled-slip span (>=3 mm accumulated
  over start-30..end, within F's latitude span +-1 deg) median 2.58 deg [1.50, 5.58]. Slip is WIDER than
  tremor in 89 % of events.
- Terminations, tremor minus slip: NORTH end median -11.7 km (median |diff| 48.1 km); SOUTH end median
  +70.7 km (median |diff| 71.7 km); pooled median |diff| 58.6 km. The modelled slip runs systematically
  FURTHER SOUTH than the tremor, and matches the tremor much better at the north end.
- Propagation direction: 13 northward, 14 southward, 0 unclear (23 from S4's `v_t_km_d`, 4 re-derived as
  the Theil-Sen slope of the daily tremor-centroid latitude because |v_t| <= 0.2 km/d: ids 1625, 1878,
  1124, 1583).
- BARRIER TEST, 54 tremor terminations in 0.25-deg bins (32 bins, 42.25-50.25N) against 5,000 draws that
  keep each event's extent LENGTH and reposition it uniformly inside its segment domain:
    widened domain (north 43.93-50.08, central 42.30-46.84; mean freedom 3.97 deg; all 32 bins testable):
      NOTHING survives Bonferroni. One bin has raw p < 0.05: 48.88N, obs 6 vs null 1.42, p 0.0028,
      Bonferroni p 0.090. Next strongest: 47.62N (obs 4, p 0.157), 47.38N (obs 4, p 0.176).
    strict segment bounds (north 46-49.92, central 44-46; mean freedom 1.80 deg; only 24 of 32 bins
      testable — 5 observed terminations fall in bins the null can never reach and were EXCLUDED rather
      than being scored as barriers): again nothing survives Bonferroni; raw p < 0.05 at 45.12N (obs 3 vs
      0.67, p 0.0086), 45.62N (obs 3 vs 0.81, p 0.020) and 48.88N (obs 6 vs 2.57, p 0.036).
- EPOCH PERSISTENCE of the 48.88N candidate: 2 of 16 terminations in E1 (null 0.44, p 0.069), 2 of 14 in E2
  (null 0.38, p 0.047), 2 of 24 in E3 (null 0.63, p 0.125). The same bin is used in all three epochs by the
  same number of events, but 2 terminations per epoch is not a testable quantity on its own.
- SLIP vs TREMOR terminations: 6 of 54 tremor terminations (11 %) sit in the 48.88N bin, but only 2 of 54
  slip terminations (4 %). The most-used SLIP termination bins are elsewhere: 45.62N (4), 46.12N (4),
  46.62N (4), 47.38N (4), 49.12N (6).
- DIRECTION: north-end latitude N-propagating median 48.67 (n 13) vs S-propagating 48.06 (n 14), MWU
  p 0.396; south-end 46.72 vs 46.29, MWU p 0.369. Terminations in the 48.88N bin: 3 of 26 for
  N-propagating (12 %) and 3 of 28 for S-propagating (11 %). No difference in either test.

### Caveats
1. Rank normalisation kills a segment-WIDE monotone sensitivity change, not a spatially heterogeneous one
   WITHIN a segment — and S3b measured exactly that heterogeneity. So any class CHANGE across epochs is
   still exposed; only class PERSISTENCE is robust, and persistence is what is reported.
2. The persistence null (shuffled labels) is far too weak to be evidence (rank autocorrelation E1-E3 is
   0.89-0.99). It confirms the map is stable; it does not test any hypothesis.
3. Ties: 10-19 % of triangles per segment-epoch have D3_nb = 0 and share an average rank, so a large tie
   block can land above the 25th percentile and be scored MIXED rather than QUIET. QUIET is an undercount.
4. Only the positive daily rate is accumulated, so S rectifies model noise upward — this is what makes the
   literal absolute thresholds meaningless and is why the percentile criteria were added.
5. Slip is 13 fixed broad vbICA patterns on a 0.25-deg fallback mesh (~28 km); adjacent triangles are not
   independent, so triangle counts hugely overstate the independent sample. The south is at the model noise
   floor (S6). Mesh-boundary triangles are flagged (`on_mesh_edge`) and most of the SLIP-ONLY set is there.
6. Tremor depth is Slab2 interface depth at the epicentre, not a hypocentre.
7. Part B: footprints are large (F latitude span 1.08-5.75 deg) and chained, so an "event" extent is partly
   a catalogue-segmentation choice; p5/p95 of detection latitude is robust but blunt. Event 1434 has only
   75 detections in window; 1655 and 1905 have 511 and 572.
8. Part B power: 27 events, extents that fill much of their segment, and Bonferroni over 32 bins. The
   strict-domain variant has essentially no power by construction. Bins the null cannot reach are excluded
   (they would otherwise score the minimum possible p from a single termination).

### Plain-language answers
(i) PERSISTENT SLIP-WITHOUT-TREMOR OR TREMOR-WITHOUT-SLIP PATCHES?
    Slip-without-tremor: NO, not anywhere the model can be trusted. There is not a single triangle in the
    north or the central segment that is high-slip-rank and low-tremor-rank in all three epochs — indeed
    SLIP-ONLY is 0.0 % of the north in EVERY epoch, because modelled slip and tremor occupancy are strongly
    co-located there (Spearman 0.71-0.77). The only persistent SLIP-ONLY triangles (7) are in the far south
    at 40.3-41.4N and 13-20 km depth, up-dip of the tremor band, 5 of 7 on the mesh boundary, in the very
    corner S6 identified as the model's noise floor — so they are a model artefact, not silent slip. The
    honest statement is a NULL RESULT with a resolution caveat: at 28-km triangles and with slip described
    by 13 broad fixed patterns, a genuinely silent patch smaller than a triangle or narrower than a vbICA
    pattern could not be seen.
    Tremor-without-slip: YES, one credible cluster. Nine triangles in central Oregon, 44.33-45.42N,
    -123.9..-123.1, Slab2 depth 27-38 km, are TREMOR-ONLY in all three epochs; none is on the mesh
    boundary and they sit 4.5-4.8 deg from the mesh limits, so this is not an edge effect. They straddle
    the central tremor band (p10 31.5, p90 37.0 km): 3 in-band, 3 up-dip, 3 down-dip. A second cluster of
    six at 49.3-49.7N and 42-47 km (down-dip end of the northern band) IS at the mesh edge (4 of 6 boundary
    triangles, 0.25-0.58 deg from the north limit) and is most likely the model thinning out, consistent
    with S5's finding that band 49.5N has the worst episode-patch matching. In DEPTH, then, the robust
    along-dip picture is: shallow 15-18 km = QUIET (locked, no tremor, no transient slip); ~35-40 km =
    SLIP+TREMOR; and the persistent TREMOR-ONLY patches sit at the DOWN-DIP edge in the north (45.3 km
    mean, ~5 km deeper than SLIP+TREMOR) but essentially AT the band depth in central Oregon (33.6 km).
    The complementary percentile map says the same thing from the other side: the 9 "high slip, no tremor"
    triangles are almost all DOWN-DIP of the tremor band (7 of 9), i.e. the model puts slip below where
    tremor stops.
(ii) PERSISTENT ALONG-STRIKE BARRIERS, AND DOES SLIP STOP WHERE TREMOR STOPS?
    Barriers: NOT DEMONSTRATED. No 0.25-deg bin survives Bonferroni correction under either null. One
    candidate stands out — 48.88N, where 6 of the 54 tremor terminations fall against a null expectation of
    1.4 (raw p 0.0028, Bonferroni p 0.090) — and it is used by exactly 2 terminations in each of E1, E2 and
    E3, so it is at least not an artefact of one epoch or of the network's growth. Under the strict segment
    bounds two Oregon bins (45.12N and 45.62N) also reach raw p < 0.05. With 27 events and extents that
    already fill much of their segment, the test simply lacks the power to promote any of these; they are
    candidates for a longer catalogue, not results. Nothing distinguishes northward- from
    southward-propagating events: their north-end and south-end termination latitudes are statistically
    indistinguishable (MWU p 0.40 and 0.37) and each group puts ~11-12 % of its terminations in the 48.88N
    bin, so the candidate barrier, such as it is, acts symmetrically.
    Does modelled slip stop where tremor stops? Only loosely, and asymmetrically. At the NORTH end the two
    agree well (median difference -11.7 km, median |difference| 48.1 km); at the SOUTH end the modelled slip
    runs systematically ~71 km FURTHER SOUTH than the tremor. Slip extents are wider than tremor extents in
    89 % of events (2.58 vs 1.69 deg median), and the slip terminations pile up in different bins from the
    tremor ones (only 4 % of slip terminations vs 11 % of tremor terminations in the 48.88N bin). Some of
    that is physical (the vbICA slip is smooth and a 3-mm contour is generous) and some is the well-known
    southern noise floor, so the safe reading is: the model does not independently confirm a tremor
    termination latitude anywhere, and its southern edges in particular should not be used as slip
    boundaries.

## 2026-09-10 orchestrator review — S7
Q5 answered. (i) Persistent SLIP-WITHOUT-TREMOR: none in the north or centre (the only 7 are far-south mesh-edge noise). Persistent
TREMOR-WITHOUT-SLIP: a credible cluster in central Oregon 44.3-45.4N at 27-38 km (9 triangles, interior of the mesh) and a
mesh-edge cluster off northern Vancouver Island (discount). Central Oregon is also where slip-tremor rank co-location is weakest
(rho 0.40-0.46 vs 0.70-0.77 north and south) in every epoch — the model's slip and the tremor decouple there. Two readings:
GNSS coverage/resolution in the Oregon Coast Range, or genuinely tremor-productive, slip-poor fault (Siletzia hypothesis region)
— cannot be separated with this model; flag for the dv/v stage and for a higher-resolution slip model. (ii) Along-dip: in the
north, tremor-only triangles sit ~5 km down-dip of slip+tremor, and 7 of 9 high-slip/low-tremor triangles are ALSO down-dip
of the tremor band -> modeled slip extends below the tremor band while tremor stays in a ~40 km band; tremor-free deep slip is
model-resolved but its amplitude is smoothing-limited. Quiet class = the 15-18 km locked zone everywhere (sanity). Absolute
"quiet-slip" thresholds are degenerate because positive-rate accumulation rectifies noise (every triangle accumulates 120-1280 mm
in 16 yr) — read the percentile versions only.
Q17 answered: NO persistent along-strike barrier survives Bonferroni (candidate 48.88N, 6/54 terminations vs 1.4 null, p_Bonf
0.09, used in every epoch). Modeled slip extents are wider than tremor extents in 89 % of events and run ~70 km further south
(smoothing/mixture), so modeled slip does not stop where tremor stops. Propagation direction does not change termination sites.
Persistence "86 % vs 13 % null" is near-vacuous (rank autocorrelation 0.9-1.0) — agent flagged it correctly.

## 2026-09-10 S9 within-event mechanics

Ran `scripts/s9_within_event.py` (PLAN.md Round 2, S9; OMP_NUM_THREADS=4, seed 20260910, **70 s** end to
end, deterministic — two runs produced byte-identical numbers). Loaders copied verbatim from
`scripts/s4_fronts.py` (Ld/A_raw -> `-1e3*(Ld @ A_raw.T)`, diff -> 7-d centred mean -> positive part;
`slip_patches_5mm.npz` labels; PNSN `pnsn_tremor_cascadia_full.csv`; Slab2 depth per epicentre from
`sse_strain.megathrust.slab_depth_interpolator`; x/y km frame). The 27 CONFIRMED north/central events of
`results/s1_events.csv` (21 north, 6 central). Full definitions in the script docstring.

Common setup: window WA = [start-10, end+10]; footprint F = triangles labelled with the patch id on any
day of DUR; neighbourhood NB = triangles within 60 km of an F centroid (102-546 of 742); near tremor =
detections within 60 km of an F centroid (64k-366k per event, all-time). A day is DEFINED at >= 5 near
detections; an event needs >= 8 defined days. Nulls = 200 rigid shifts of the near-tremor DATES by
>= 365 d, wrapped over 2010-01-01..2026-05-27 (the S4 control), with the SLIP series left unchanged.

Because S1b showed raw depth quantities mostly report *which days the 5 mm/30 d threshold selected*
(Spearman(north-wide window depth, Q3) = +0.72), **every depth quantity here is an anomaly**: Part A
subtracts each event's own overall median (only the within-event SHAPE survives), Part C subtracts the
CONCURRENT same-segment margin-wide percentile, exactly as S1b's "Consequence for S2" demanded.

### PART A — tremor depth migration within events (Q14)

Stacked fifth-binned anomaly (median over 27 events; + = deeper; bootstrap 90 % CI over events):

    fifth        1        2        3        4        5      first->last   null mean+-sd    p
    tremor    +0.44    -0.41    -1.08    +1.33    +2.70      **+2.26**    -0.14 +- 0.89   0.020
      90%CI  [-0.06, [-0.47, [-1.69, [-0.22, [+1.24,
              +1.82]  -0.12]  -0.59]  +2.71]  +5.04]
    p10       +0.21    -0.11    -0.65    +1.45    +3.15      **+2.94**    -0.11 +- 0.92  <0.005
    slip      -0.11    +0.20    +0.06    -0.39    +0.77        +0.85   boot [-1.15,+1.61]  0.70
    log10 N   +0.04    +0.14    -0.14    -0.69    -0.66        -0.87 dex  (x0.13 of fifth 2)

- **Tremor DEEPENS through an event**: flat/slightly shallow for the first three fifths, then +2.3 km into
  the last fifth (200-shift null p = 0.020; the null band is [-0.9, +1.5] km in every bin).
- **The up-dip-most tremor deepens MORE** (+2.9 km, 0/200 null draws as large): the daily p10 front
  **retreats down-dip during the event — it does not advance up-dip**.
- The modeled slip depth is **flat** (+0.85 km, CI straddles zero, 15/27 positive, p = 0.70). Treat this
  as a null of the MODEL, not of the fault: the 13 fixed patterns pin the slip-weighted depth at 36-39 km
  (S1b: component depths all 30-35 km) and cannot express a 2-3 km migration.
- **Strength is modest per event.** Per-event z against its own 200 shifts: median +0.51 [-0.23, +1.24],
  18/27 positive (sign p = 0.12), Wilcoxon p = 0.039, and only **2/27 events individually p < 0.05**.
  This is a group tendency, not a per-event rule.
- **Robustness to the +-10 d pads:** recomputing the whole profile on DUR only gives +1.86 [-0.78, +4.71]
  (n = 24; the 2-4 d events 1655/1905/1434 lose the profile) vs +2.56 [-1.29, +5.21] padded, Spearman
  between the two versions rho = +0.51, p = 0.010. The pads are not making the result.
- **COUNT-BIAS CONTROL (important).** Detections/day collapse by 0.87 dex (a factor ~7) from the first to
  the last fifth. A daily p10 from fewer detections is pulled toward the daily median, i.e. biased DEEP,
  so part of the +2.9 km p10 change may be an estimator artefact; across events Spearman(count change,
  depth change) = -0.40 (p = 0.037) for the median and -0.49 (p = 0.009) for the p10. The daily MEDIAN is
  count-unbiased as an estimator, so **the +2.26 km median result stands; the p10 number is the
  secondary, partly count-biased one**. Physically the two facts are the same fact: the intense phase is
  shallow, the waning tail is deep.
- **Depth vs intensity within events**: per-event Spearman(daily median depth, daily count) median
  **-0.107** [-0.32, +0.09]; **19/27 negative** (sign p = 0.052), Wilcoxon p = 0.049; 7/27 nominally
  significant, 6 of them negative. So tremor is **shallower when it is most intense**.
- **No dependence survives.** segment: central +3.70 [+2.67, +10.91] (n=6) vs north +1.84 [-3.19, +4.75]
  (n=21), Mann-Whitney p = 0.065; epoch 2010-14 +3.89 / 2015-19 +0.22 / 2020-26 +2.03, Kruskal p = 0.34;
  Mw rho = +0.03 p = 0.88; date rho = -0.21 p = 0.29; propagation direction (sign of v_t) northward +2.95
  (n=14) vs southward +0.22 (n=13), p = 0.64.

### PART B — tremor rate vs slip rate scaling (Q13)

b = Theil-Sen slope of log10(C+1) on log10(R); C = daily near-footprint tremor count, R = area-weighted
positive modeled slip rate over NB (mm/d km^2, a potency rate); days with C >= 1 and R > 0 (8-72 per event).

    group                 b  median [IQR]          |  pooled, event offsets removed (OLS / Theil-Sen, n days)
    combined  (27)   +0.401 [+0.181, +0.478]       |  +0.328 / +0.397   (947 d, r = +0.36)
    north     (21)   +0.432 [+0.193, +0.503]       |  +0.357 / +0.413   (750 d)
    central    (6)   +0.356 [+0.078, +0.407]       |  +0.220 / +0.288   (197 d)
    2010-2014  (8)   +0.410 [+0.355, +0.455]       |  +0.313 / +0.373   (283 d)
    2015-2019  (7)   +0.409 [+0.299, +0.461]       |  +0.327 / +0.427   (290 d)
    2020-2026 (12)   +0.297 [+0.148, +0.513]       |  +0.342 / +0.379   (374 d)
    NULL (shifted C, unchanged R): per-event mean -0.004 [-0.016, +0.008], per-event sd 0.278 [0.24, 0.31]

- **The scaling is real and strongly sub-linear: C ~ R^0.4.** Per-event z against its own null: median
  +1.18 [+0.55, +1.90], **26/27 positive (sign p = 4.2e-07)**; 10/27 events individually p_null < 0.05.
  A 10x rise in modeled slip rate buys only ~2.5x more detections.
- **The exponent is CONSTANT across segments and epochs.** Every group median is 0.30-0.43 and every
  pooled within-event fit 0.22-0.43, with IQRs that overlap completely; nothing distinguishes north from
  central or any epoch from another.
- **Hysteresis: a weak, non-significant tendency for a steeper decaying limb.** rising +0.258
  [+0.10, +0.44] vs decaying +0.413 [+0.17, +0.55]; rise - decay -0.119 [-0.28, +0.10], Wilcoxon
  p = 0.075, 9/23 positive. Not claimable at n = 23.
- CAVEATS (stated in the script output): **R is 7-d smoothed and mixture-limited** (13 fixed patterns), so
  its daily dynamic range is compressed and b is biased HIGH relative to a true count-vs-slip-rate
  exponent — read 0.4 as an upper bound on the exponent, i.e. the true response is at least this
  sub-linear. **C is completeness-limited** (PNSN growth), which moves the per-event LEVEL far more than
  the slope; the per-event offsets/dummies remove exactly that level. Thin events: 1905 (8 days, b = 1.00)
  and 1434 (8 days, b = -0.05) carry the tails of the distribution.

### PART C — edges referenced to concurrent tremor (Q2 redo)

Reference percentiles: north all-time tremor p10/p50/p90 = 34.77 / 39.43 / 45.36 km (n = 276,633);
central 31.53 / 34.38 / 36.97 km (n = 51,574).

    quantity (km)                      combined                north (21)             central
    rel_p10  (near - concurrent seg)   -0.018 [-0.28,+0.04]    -0.005 [-0.09,+0.01]   -0.33 [-0.83,+0.10] (5)
    rel_p90                            +0.003 [-0.08,+0.19]    -0.000 [-0.13,+0.01]   +1.02 [+0.45,+8.05] (5)
    slip_updip  - near p10 (concurr.) -11.034 [-15.15,-8.16]  -11.38 [-15.18,-9.76]   -9.59 [-13.35,-7.66]
    slip_downdip- near p90 (concurr.)  +5.246 [+1.56,+8.96]    +5.25 [+2.83,+11.29]   +3.29 [-3.05,+7.68]
    slip_updip  - seg all-time p10     -9.964 [-17.30,-7.97]  -11.02 [-17.44,-7.88]   -9.82 [-9.96,-8.46]
    slip_downdip- seg all-time p90     +4.791 [+1.44,+8.30]    +5.41 [+1.76,+8.59]    +3.76 [+1.35,+7.20]

- **The modeled slipping band strictly overhangs the tremor band on BOTH sides.** The slip up-dip edge is
  **~11 km shallower than the tremor up-dip edge in 27/27 events** (median -11.0 km vs the concurrent
  near-patch p10; -10.0 km vs the segment all-time p10), and the slip down-dip edge is ~5 km deeper than
  the tremor p90 in 22/27.
- **Neither edge drifts.** 21 Theil-Sen trends (7 quantities x combined/north/central, 1000 permutations):
  **min p = 0.143** (central rel_p90, n = 5); Bonferroni needs 0.0024. Tightest bounds: north rel_p10
  -0.0005 km/yr [-0.014, +0.031] p = 0.92; combined slip_updip_rel -0.000 km/yr [-0.58, +0.40] p = 1.00;
  combined slip_updip - near_p10 -0.126 km/yr [-0.65, +0.37] p = 0.58.
- **The concurrent reference is near-vacuous in the north** — rel_p10 and rel_p90 are 0.00 +- 0.1 km — and
  that IS the finding: during its own event, a patch's tremor depth distribution is **indistinguishable
  from the whole northern band on the same days**, exactly as S1b's "the 40-km band around the union of
  northern footprints contains 100 % of lat>=46N detections" predicted. Central is a real filter (n = 5)
  and there the patch tremor reaches ~1-8 km DEEPER than the concurrent central-wide p90 — n = 5, do not
  read it as more than a hint.

### Plain answers

**(i) Does tremor migrate up-dip or down-dip during an event, and does slip?**
Tremor migrates **DOWN-DIP**: the stacked daily median Slab2 depth is flat for the first ~60 % of an event
and then deepens **+2.3 km into the last fifth** (200-shift null p = 0.020); the up-dip-most tremor of each
day deepens **+2.9 km**, i.e. **the up-dip front retreats, it does not advance**. This is a group tendency
(18/27 events positive, only 2/27 individually significant) and it is the same fact as "tremor is shallower
when most intense" (within-event Spearman(depth, count) median -0.11, 19/27 negative): the deep tremor is
the sparse waning tail (detection rate falls ~7x over the same interval), which also means part of the p10
number is a small-sample bias — the count-unbiased median result is the one to quote. The modeled slip depth
shows **no** migration (+0.85 km, CI straddles zero), but the 13-fixed-pattern mixture cannot express a
2-3 km depth change, so that is a null of the model, not of the fault.

**(ii) How does tremor rate scale with slip rate, and is the exponent constant?**
**Sub-linearly, C ~ R^0.40** (median exponent 0.401 [0.18, 0.48]; pooled with per-event offsets removed,
OLS 0.33 / Theil-Sen 0.40 over 947 event-days). Decisive against a time-shifted null (null exponent
-0.00 +- 0.28; 26/27 events positive, p = 4e-07). **It is constant**: north 0.43 vs central 0.36, epochs
0.41 / 0.41 / 0.30, pooled 0.22-0.43 — no segment or epoch is distinguishable. The rising and decaying
limbs give 0.26 vs 0.41 (p = 0.075): a hint of hysteresis toward a steeper decay limb, not claimable.
Because R is 7-d smoothed and mixture-compressed, 0.4 is an **upper bound** on the exponent — the real
tremor-rate response to slip rate is at least this sub-linear.

**(iii) Where is the slip up-dip edge relative to the tremor up-dip edge, and does either move?**
The slip up-dip edge sits **~11 km SHALLOWER (up-dip) of the tremor up-dip edge, in 27 of 27 events**
(median -11.0 km against the concurrent near-patch tremor p10, -10.0 km against the segment's 2010-2026
p10), and the slip down-dip edge ~5 km deeper than the tremor p90 (22/27): the modeled slipping band
brackets the tremor band. **Neither edge drifts over 2010-2026** — 21 trend tests, min p = 0.143, with the
tightest bounds ~0.03 km/yr (north tremor edges) and ~0.6 km/yr (slip edges). Caveat that governs the
whole answer: S1/S2 already established that the slip "edges" are properties of the 13 broad fixed
patterns, not resolved slip edges, so **11 km is an upper bound on the up-dip overhang, not a measurement
of it**; and referencing the tremor to concurrent margin-wide tremor removes essentially nothing in the
north, because the patch tremor IS the northern band.

### Issues / caveats specific to S9
- **n = 27, of which 6 central.** Every central number (n = 5-6) is a hint at best.
- **Three very short events** (1655 3 d, 1905 4 d, 1434 2 d) fail the DUR-only robustness profile and carry
  the extreme Part B exponents; 1434 was already flagged unreliable in S1 (N = 75) and S4 (negative peak r).
- **Event 1434 has NO concurrent same-segment reference**: its footprint spans 44.9-46.9N (labelled central
  by its slip-weighted centroid), and during its 2-day DUR there was zero tremor in 44-46N, so rel_p10/p50/p90
  are NaN for it (Part C n = 26 for those three columns). All 75 of its near detections are north of 46N.
- **The null shifts the tremor, never the slip**, so it tests "is this specific event window special for the
  tremor", not "is the slip model special". A slip-side null would need the mirror catalogue (S1) and is not
  part of S9.
- **Phase normalisation uses the first and last >= 5-detection day inside WA.** In the north tremor is
  almost always above that threshold, so the phase axis is close to the window itself; the +-10 d pads were
  kept for symmetry and the DUR-only control shows they do not drive the result.
- **6.4 % of detections have no Slab2 depth** and are dropped from every depth percentile but still count
  toward the >= 5/day rule and toward C in Part B (same convention as S1/S4).
- **Multiplicity.** Part A reports one headline with its own null; Part C's 21 trend tests are reported with
  the Bonferroni threshold and nothing approaches it. Part B's segment/epoch comparisons are descriptive
  (no test claimed significant).

Files written: `scripts/s9_within_event.py`, `results/s9_depth_migration.csv` (27 rows x 52 cols),
`results/s9_depth_migration_stack.csv` (20 rows = 4 profiles x 5 fifths),
`results/s9_rate_scaling.csv` (27 event rows + 6 pooled rows), `results/s9_edges.csv` (27 rows),
`results/s9_edges_trends.csv` (21 rows), `figures/s9_within_event.png` (9 panels: stacked depth vs phase
with null band and the log-count overlay; p10 front vs phase; rho(depth, count) histogram; exponent vs
null; exponent by segment/epoch; rising vs decaying exponents; tremor edges vs date; slip edges vs date;
summary). Nothing outside `tremor_slip_geometry/` was created or modified.

## 2026-09-10 orchestrator review — S9
(i) Q14: tremor deepens through an event: first->last fifth +2.3 km (p 0.02 vs shifted windows; 18/27 events positive), and
tremor is shallower when most intense (rho -0.11, 19/27, p 0.05) -> the late-stage tremor "tail" sits down-dip of the main
front, the up-dip front does not advance. No dependence on segment, Mw, epoch, direction. Modeled slip depth is flat (pinned by
the fixed patterns; not a test). Count-bias caveat handled by the agent (quote the median-depth result, not p10).
(ii) Q13: daily tremor count vs modeled potency rate: C ~ R^0.40 [0.18, 0.48], 26/27 events positive (null 0.00 +- 0.28,
p 4e-7), the same in north/central and in all epochs. Two OPPOSITE biases bound it: 7-d smoothing + mixture compresses R
(pushes the exponent UP) while PNSN detections are 5-min windows, so daily counts saturate near 288 (pushes it DOWN at high
rates) — see the saturation check below. Read as: tremor rate is a monotone, sub-linear (or saturating) function of slip rate,
with an exponent that has not changed 2010-2026. Hysteresis (rising 0.26 vs decaying 0.41, p 0.075) is a hint only.
(iii) Q2 answered: the modeled slipping band brackets the tremor band — slip up-dip edge 11 km shallower than the concurrent
tremor p10 in 27/27 events, down-dip edge 5 km deeper in 22/27; both are upper bounds (broad fixed patterns). No drift of any
edge (21 trends, min p 0.14; tremor-edge bounds ~0.03 km/yr). In the north a patch's tremor is indistinguishable from the whole
northern band on the same days (rel_p10/p90 = 0.00 +- 0.1 km) — the north behaves as one tremor band, patches are not local.
Saturation check (orchestrator): detection times sit on 2.5-min multiples and days reach up to 1011 detections (p99 695), so the
daily count is NOT capped at 288 — several epicentres per window are catalogued. The downward saturation bias on the exponent is
weaker than feared; the smoothing/mixture upward bias remains. Exponent 0.4 = upper bound, sub-linear scaling stands.

## 2026-09-10 round 3 launched (publishability): S10 literature, S11 scaling robustness (energy/duration 2017+), S12 representability
Catalog check (orchestrator): energy/duration/num_stas/magnitude populated from 2017 (64.5 % of detections). Repo resolution_test()
is a field-roughness test, not a network-resolution test -> S12 builds a representability map instead.

## 2026-09-10 S12 representability

Ran `scripts/s12_representability.py` (PLAN.md Round 3, S12) end to end, `OMP_NUM_THREADS=4`,
~55 s, deterministic (no random numbers anywhere). Outputs:
`results/s12_representability_triangles.csv` (742), `s12_summary.csv` (17),
`s12_pattern_energy.csv` (17), `s12_stations.csv` (213), `s12_misfit_by_latitude.csv` (12);
`figures/s12_representability.png` (9 panels). Nothing outside `tremor_slip_geometry/` was
created or modified.

### What the solution files actually contain (checked, not assumed)
`/home/jovyan/cascadia-sse-strain/data/2026-08-18/` — the six files of one dated Gualandi
near-real-time solution:
- `X.mat` -> `Xd`: `llh` (639 x 3), `name` (639 strings `ALBHe/ALBHn/ALBHu`...), `timeline`
  (7169 decimal years), `ts`, `var_ts`, `unit` = **mm**, `data_mask`, `centering_offsets`.
  **Station coordinates ARE recoverable**: 639 rows = **213 stations x (e, n, u)**;
  `Solution.station_llh` decimates by 3 and the names give the 4-character IDs.
- `ICA.mat` -> `ICA_essential`: `U` (639 x 26 spatial responses), `var_U` (639 x 26 data
  variances, the Cd of the published inversion), `S` (26), `V` (7169 x 26), `var_V`.
- `ind_comps.mat`: the **13** retained components, 1-based 5, 10, 11, 15, 16, 17, 20, 21, 22,
  23, 24, 25, 26.
- `misfit_comps.mat`: a **(9 sigma) x (13 component) cross-validation grid**, NOT a per-station
  misfit. Its per-column argmin is what recovers the per-component smoothing hyperparameter
  (`ind_sigma0_comps` is absent from the public download); here sigma0 = 0.00316 (1 comp),
  0.01 (7), 0.0316 (5).
- `options.mat`: `inversion.lambda_strike = lambda_dip = lambda0 = 21.0895 km`, exponential
  prior, `sigma` grid 1e-4..1e-1 (9 values), `m0 = 0`, `rake_pos`, `fault.nu`; plus `scen`
  (27 selection circles, outlier settings) and `tm` (trajectory-model) blocks. No geometry.
- `fitresult.mat`: trajectory-model fit; not used.
Reproduced inversion: G = 639 x 1484 (cutde), L = 1484 x 13, per-component fractional misfit
0.52-0.94 (median 0.80). Thrust patterns P = -Ld, 742 x 13, **rank 13**. Mesh depth confirmed
= `-centroid_llh[:,2]` and identical to `patch_depth_km` in `slip_patches_5mm.npz`.

### A1 REPRESENTABILITY — can the fixed basis express a compact patch there?
Unit thrust Gaussian test patch u_j (peak 1 at triangle j, 3-D on-fault centroid distance),
projected onto span(P) by thin QR. `r` = fraction of the patch's L2 energy in the span
(a random 742-vector scores 13/742 = 0.018); `pot_frac` = fraction of its POTENCY reproduced.
                                        r30     r50    pot_frac30  pot_frac50
  north SLIP+TREMOR (145 tri)          0.486   0.710     0.951       0.925
  **S7 central-Oregon TREMOR-ONLY (9)  0.295   0.432     0.763       0.854**
  north, ETS band 25-45 km (191)       0.395   0.626     0.867       0.861
  central, ETS band 25-45 km (88)      0.330   0.511     0.913       0.900
  south, ETS band 25-45 km (94)        0.507   0.703     0.956       0.925
  north TREMOR-ONLY, mesh edge (6)     0.114   0.307     0.392       0.521
  QUIET / locked 15-18 km (91)         0.134   0.249     0.638       0.716
So the basis is **worse but not blind** at the central-Oregon cluster: it recovers **76 % of the
potency of a sigma=30 km patch there against 95 % in the north (0.80x), and 85 % vs 93 % at
sigma=50 km (0.92x)**; shape fidelity r30 is 0.61x the north's. Per 1-degree band the median r30
over 25-45 km triangles is 0.598 (40N), 0.417, 0.522, 0.494, **0.311 (44N), 0.370 (45N)**, 0.507
(46N), 0.594, 0.531, 0.240 (49N, mesh edge) — 44-45N is the worst-represented interior band, but
by a factor <2, and the mesh-edge cluster S7 already discounted scores 3x worse still (0.114 /
0.392), which is the correct behaviour of the metric.

### A2 PATTERN ENERGY — where is the model's slip capacity vs where is the tremor?
Area-weighted |thrust|, aggregated over the 13 components weighted by the s.d. of each amplitude
time function, as a share of the margin:
                       north    central 44-46N   south
  slip capacity        0.507      **0.163**      0.330
  mesh area            0.503        0.196        0.301
  tremor-days D3nb     0.429        0.106        0.465
  tremor detections    0.391        0.070        0.539
Central Oregon holds 16.3 % of the model's slip capacity against 19.6 % of the mesh area — a
**0.83x deficit, not a hole**; per component the central share runs 0.10-0.21 with no component
near zero there. Crucially the tremor share of central Oregon (10.6 % of tremor-days, 7.0 % of
detections) is **LOWER than its slip-capacity share**, so on a margin-wide budget there is no
"lots of tremor, no capacity" mismatch; the S7 result is about the local RANK pattern, not about
a margin-scale imbalance.

### A3 GNSS COVERAGE — recoverable, and central Oregon is NOT under-covered
213 stations; 161 of them lie within 150 km of the deep (30-45 km) ETS band. Per 1-degree band
(stations / within 150 km of the deep band / median nearest-station distance for deep triangles):
40N 38/32/14.0 km, 41N 14/7/22.1, 42N 15/11/16.9, 43N 12/8/22.6, **44N 15/14/16.2, 45N
21/15/17.9**, 46N 28/18/14.5, 47N 38/29/14.7, 48N 23/20/15.9, 49N 7/6/28.2. Median
nearest-station distance in the ETS band: **central 17.2 km, north 18.0 km, south 17.8 km**.
For the 9-triangle cluster: nearest station 17.1 km, 17 stations within 100 km (north
SLIP+TREMOR: 15.1 km, 24 stations). Modestly fewer neighbours, essentially the same proximity.

### A4 MISFIT — flat with latitude
`misfit_comps.mat` is a CV grid, so per-station misfit was computed from the reproduced
inversion: fraction of each station's 3-component ICA spatial response, summed over the 13
retained components, that the slip model fails to explain. Median by segment: **central 0.816,
north 0.793, south 0.826**. By 1-degree band: 0.828, 0.922, 0.715, 0.737, **0.776 (44N), 0.855
(45N)**, 0.796, 0.808, 0.688, 0.825, 0.940. 45N is mildly the worst of the interior bands but
41N and 50N are worse; there is no central-Oregon misfit anomaly.

### A5 SYNTHETIC RECOVERY — the decisive test (done, not skipped)
Station coordinates and the displacement Green's functions are both available, so the full test
was run: the sigma=30 km thrust patch at each triangle was forward-modelled to the 213 stations
(`D = G M`), then recovered with the **real** Tarantola-Valette machinery (per-component sigma0,
lambda 21.09 km exponential prior, Cd = diag(var_U) — the published inversion's own ICA-level
noise), `M_rec = GC' A^-1 D`, **without the 13-pattern restriction**. Because the estimator is
linear, the signal recovery is deterministic and the noise is fully described by the posterior
covariance, so no Monte Carlo is needed; `det_threshold_rel = 2 sqrt(Cm_post_jj)/res_peak_j`
normalised per component to the median over the north SLIP+TREMOR triangles.
                                   res_pot  res_peak  det_rel   max stn disp   stns >0.5 mm
                                                                (mm per 1 cm slip)
  north SLIP+TREMOR (145)           0.919     0.517    1.000        1.93            19
  **central TREMOR-ONLY (9)         0.986     0.556    0.909        2.02            15**
  north TREMOR-ONLY, mesh edge (6)  0.635     0.155    3.767        1.59             5
  QUIET, shallow locked (91)        0.903     0.393    1.379        1.93             6
**The network and the prior recover 98.6 % of the potency of such a patch under the central-Oregon
cluster — marginally better than in the north — and it is 0.91x as hard to detect there as a
typical northern ETS patch.** The same test flags the mesh-edge northern TREMOR-ONLY cluster as
3.8x harder to detect with only 64 % potency recovery, i.e. it reproduces S7's edge-effect
judgement from first principles and therefore has the discriminating power claimed for it.

### A6 CLASS ROBUSTNESS — undoing the basis handicap does not dissolve the cluster
Re-ranking each central triangle by `S_tot_mm / pot_frac30` (i.e. crediting every triangle with
the slip the basis could not express) moves the 9 from a median within-segment slip rank of
13.2 to 9.7 (5 rise, 4 fall — the four with pot_frac30 > 0.85 lose ground to neighbours that gain
more). **The worst case is rank 23.6: all 9 stay below the 25th percentile, so the TREMOR-ONLY
classification survives a full first-order correction for representability.**

### Caveats
1. **The representability metric is not fully independent of the answer.** P is estimated FROM
   the GNSS data; if little slip occurred in central Oregon the vbICA will not have built
   patterns there, so a low r at 44-45N is a *consequence* as much as a *cause*. This is why the
   verdict rests on A3/A4/A5, which use only the station geometry, the Green's functions and the
   published noise model, and never touch P.
2. A5 tests the linear step (ICA spatial response -> slip), not the vbICA step (positions ->
   components). The surviving loophole is slip that is temporally independent of all 13 time
   functions: it would need its own component and would not get one. Central-Oregon tremor is
   largely synchronous with the adjacent episodes, so most of it should be carried by existing
   components — but a genuinely asynchronous central-Oregon transient is outside what S12 tests.
3. A1 is generous to the model: in the real inversion a component's amplitude drives its
   strike-slip pattern too, so expressing a thrust-only patch is at least as hard as stated.
4. `pot_frac` exceeds 1 for some triangles (up to 1.29 in central Oregon): the basis overshoots
   the potency by adding slip elsewhere. `leak30` (the share of |u_hat| beyond 2 sigma) is
   0.74-0.91 everywhere, i.e. every reconstruction is heavily smeared — the 13-pattern model is
   a broad-scale description everywhere on the margin, not only in Oregon.
5. Test-patch size is a choice: sigma 30 and 50 km were run. The central/north contrast shrinks
   as sigma grows (0.80x -> 0.92x in potency), as it must — the basis is a smooth one.
6. The reproduced inversion is not bit-identical to Gualandi's (sigma0 recovered from the CV
   argmin, mesh is the local `cascadia_mesh.txt`), and per-component fractional misfits are high
   (0.52-0.94) because the ICA spatial responses contain non-tectonic signal. Relative
   comparisons between latitudes are unaffected; absolute misfit levels should not be quoted.
7. The 9 triangles are adjacent on a 28 km mesh and the prior correlation length is 21 km, so
   they are ~2-3 independent resolution elements, not 9.

### Plain verdict
**Could the model have represented slip under the central-Oregon tremor cluster? YES — with a
modest handicap.** The 13 fixed patterns reproduce **76 % of the potency of a 30-km patch there
against 95 % in the north (a 20 % relative shortfall; 8 % at 50 km)**, and the model's
amplitude-weighted slip capacity in 44-46N is 0.83x its share of the mesh area. That is a
handicap of order 20 %, not a blind spot — and a 20 % potency shortfall cannot turn a
top-quartile slip triangle into a bottom-quartile one, which A6 confirms directly: crediting
every triangle with the slip the basis cannot express leaves all 9 below the 25th percentile
(worst case 23.6).
**Is it a coverage limitation? NO.** 15 and 21 stations sit in the 44N and 45N bands, 14 and 15
of them within 150 km of the deep ETS band; the median nearest-station distance for deep
triangles is 16-18 km, the same as the north; the per-station data misfit is flat with latitude
(central 0.816 vs north 0.793); and a synthetic sigma=30 km thrust patch under the cluster is
recovered by the real inversion with **98.6 % of its potency and a detection threshold 0.91x**
that of a typical northern ETS patch. The same machinery correctly flags the northern
Vancouver Island TREMOR-ONLY cluster as unresolvable (0.64 potency, 3.8x threshold), which is
what gives the negative result its force.
**Therefore: the S7 central-Oregon decoupling is a CANDIDATE PHYSICAL RESULT, not a model or a
coverage artefact.** The residual model caveat is a ~20 % potency handicap of the fixed basis
plus the untestable loophole that a central-Oregon transient temporally independent of all 13
components would have no component to live in (caveat 2). The honest phrasing for a paper is
that central Oregon hosts persistent tremor with modelled aseismic slip in the bottom quartile
of its own segment, at a place where the geodetic network can resolve a 30-km slip patch as well
as anywhere in Cascadia — an apparent tremor-rich / slip-poor stretch of interface, to be
confirmed with a slip model that does not impose a fixed low-rank spatial basis.

## 2026-09-10 orchestrator review — S12
The central-Oregon tremor-without-slip cluster is NOT a coverage or representability limitation: 213 stations recovered from the
solution files; nearest-station distance (17 km) and misfit (0.82 vs 0.79) as good as the north; a synthetic 30-km thrust patch
under the cluster is recovered at 98.6 % potency with the real prior and data covariance (detection threshold 0.91x the north),
while the same test correctly condemns the mesh-edge Vancouver Island cluster (0.64, 3.8x). Fixed-basis handicap only ~20 %
(pot_frac30 0.76 vs 0.95), far too small to move a top-quartile slip triangle to the bottom quartile. Candidate PHYSICAL result:
persistent tremor over a part of the interface where the resolvable slow slip is bottom-quartile, in the up-dip part (27-38 km)
of the central Oregon tremor band. One loophole left (temporal: a transient with its own time history has no component) ->
S12b launched (component time functions vs central tremor; cluster slip-rate vs tremor covariance; data-space check at the
nearest GNSS stations).

## 2026-09-10 S11 scaling robustness

Ran `scripts/s11_scaling_robustness.py` (PLAN.md Round 3, S11; OMP_NUM_THREADS=4, seed 20260910,
**57 s** end to end). Loaders copied verbatim from `scripts/s9_within_event.py`; same 27 CONFIRMED
north/central events of `results/s1_events.csv`, same window [start-10, end+10], same 60-km
footprint neighbourhood. **Cross-check passed exactly**: b(C vs R7) for event 1601 = 0.542195 in
S11 vs 0.542195 in S9, and the whole combined median reproduces S9's +0.401 [+0.181, +0.478]. The
fast Theil-Sen used everywhere agrees with `scipy.stats.theilslopes` to 10 decimals.

Question: is S9's C ~ R^0.40 real, or (a) count saturation, (b) network sensitivity, (c) an
artefact of how the slip rate is smoothed?

### Setup — measures and slip-rate variants
Tremor measures per event-day (all from the same near-footprint set): **C** detections;
**U** distinct 2.5-min windows with >= 1 detection (100.0 % of PNSN times are exact multiples of
150 s -> hard cap 576/day; a duration-like, saturating measure); **Uloc** distinct (window,
0.1-deg cell) pairs; **CperU** = C/U multiplicity; and, for the 18 events from 2017 on (energy /
duration / num_stas / magnitude are 0/NaN before 2017-01-01 and populated for every row after),
**E** daily summed energy, **Dur** summed duration, **Mmax** daily max magnitude, **EperDet** =
E/C, **C_ns8** count restricted to num_stas >= 8.
Slip-rate variants, all area-weighted sums over the 60-km neighbourhood (mm/d km^2): **R7** =
7-d centred mean of the positive daily potency rate (S9's), **R3** 3-d, **R1** raw daily,
**Rcum** = daily increment of the signed cumulative neighbourhood potency, clipped at 0.
Nulls: 200 rigid >= 365-d wrapped shifts of the near-tremor dates with the slip series unchanged
(the S4/S9 control), for **every** measure x variant x subset; the 2017+ measures use a separate
2017-01-01..2026-05-27 wrap so a shifted window never lands where those columns are structurally
zero.

### 1. Exponent matrix (median over events [IQR], 90 % bootstrap CI of the median)

| measure | vs R7 | vs R3 | vs R1 | vs Rcum | n |
|---|---|---|---|---|---|
| C       | **+0.401** [+0.181,+0.478] CI [+0.31,+0.45] | +0.198 | +0.108 | +0.128 | 27 |
| U       | +0.377 [+0.174,+0.447] CI [+0.30,+0.42] | +0.175 | +0.089 | +0.088 | 27 |
| Uloc    | +0.400 | +0.198 | +0.108 | +0.128 | 27 |
| CperU   | +0.011 | +0.007 | +0.001 | +0.000 | 27 |
| E       | **+0.644** [+0.166,+0.840] CI [+0.33,+0.76] | +0.392 | +0.232 | +0.250 | 18 |
| Dur     | +0.368 | +0.190 | +0.116 | +0.065 | 18 |
| Mmax    | +0.147 mag/decade | +0.111 | +0.072 | +0.057 | 18 |
| EperDet | +0.231 | +0.143 | +0.093 | +0.061 | 18 |
| C_ns8   | +0.382 | +0.275 | +0.133 | +0.091 | 18 |

Against the shifted-tremor null, C/U/Uloc vs R7 are 26/27 positive (sign p = 4.2e-07), vs R1
24/27 (p = 4.9e-05); E vs R7 15/18 (p = 7.5e-03), E vs R1 16/18 (p = 1.3e-03). Rcum is the one
variant that fails (17/25, p = 0.11): it is zero on 48 % of event-days (the neighbourhood is net
back-slipping) and its log range is only 2.9 dex vs 4.6-5.1 dex for the others.
**b < 1 in 27/27 events for C and 14/18 for E**; the bootstrap CI of the median excludes 1 for
every measure in the table.
**Uloc is identical to C to the last digit** — within a 2.5-min window no two catalogued
epicentres ever share a 0.1-deg cell, so the space-time occupancy measure adds nothing.

### 2. Saturation — NOT the explanation
- Max U over all 956 event-days = **546 of the 576 cap**; 1.46 % of days have U >= 500 (6 events
  ever reach it), 9.2 % have U >= 400. Multiplicity C/U has median 1.056 and never exceeds 1.66.
- **b(C vs U) = +1.028 [+1.017, +1.043]**, range 1.000-1.105 across events: the count tracks the
  occupancy essentially one-for-one, so C is not being carried by multiple epicentres per window.
  Decomposition checks: median b_U +0.377 + median b_CperU +0.011 = +0.387 vs median b_C +0.401.
- Dropping the days with U >= 400 changes nothing: b = +0.385 vs +0.401, paired difference
  **+0.000 [-0.025, +0.000], Wilcoxon p = 0.41**.
- **Upper vs lower half of log10 R**: a saturating relation must flatten at high rate. It does not
  — C vs R7 lower +0.254, upper +0.248, upper - lower **+0.099 [-0.256, +0.623], p = 0.35**,
  14/25 positive; C vs R1 +0.069 p = 0.20; U vs R7 +0.098 p = 0.51. If anything the exponent is
  slightly LARGER at high rate.

### 3. Network sensitivity — not the explanation either
- b by epoch: 2010-2014 +0.410 / 2015-2019 +0.409 / 2020-2026 +0.297 (Kruskal p = 0.94; Spearman
  vs date rho = -0.096, p = 0.63). Same for U (p = 0.91) and R1 (p = 0.26).
- Spearman(b, event mean num_stas) over the 18 2017+ events: +0.317 (p = 0.20) for C, +0.278
  (p = 0.26) for E, +0.313 (p = 0.21) for U — nothing.
- Restricting to detections seen by **>= 8 stations** RAISES the exponent slightly: b_C +0.370 ->
  b_C_ns8 +0.382, paired difference +0.074 [+0.003, +0.109], Wilcoxon p = 0.012, 14/18 positive.
  A detection-threshold artefact would push the other way.

### 4. Energy — sub-linear too, and counts DO under-represent high-rate tremor
On the same 18 events: **b_E = +0.644 [+0.166, +0.840]** vs **b_C = +0.370 [+0.166, +0.494]**;
b_E - b_C = **+0.222 [+0.021, +0.387], Wilcoxon p = 0.0034**, 14/18 positive (against R1: b_E
+0.232 vs b_C +0.127, difference +0.082, p = 0.0028). Equivalently **energy per detection rises
with slip rate: b(E/C vs R7) = +0.231 [+0.004, +0.296], Wilcoxon p = 0.0028**, and the daily
maximum magnitude rises **+0.147 magnitude units per decade of R** (CI [+0.03, +0.23]).
So the count under-represents the high-rate end by ~0.23 in the exponent — but the amplitude-
weighted measure is **still sub-linear**, with the bootstrap CI of the median [+0.33, +0.76]
excluding 1 and 14/18 events below 1. Duration is not an independent check: the column takes only
the two values 150 s and 300 s, so Dur is ~150*C by construction (b_Dur +0.368 ~ b_C).

### 5. Smoothing — the bias runs OPPOSITE to what S9 assumed, and the ladder settles which rate drives the tremor
S9 warned that 7-d smoothing biases b HIGH. The simulation says the reverse. Counts generated as
Poisson(k * (Rtruth/max)^b_true) per event and refitted exactly as in Test 1, for every
generating rate:

    generating rate    b(fit R1)/b(fit R7)   b(fit R3)/b(fit R7)
      OBSERVED DATA          0.269                 0.494
      R1  (b_true 0.4)       4.597                 2.182
      R3  (b_true 0.4)       1.250                 2.732
      R7  (b_true 0.4)       0.190                 0.367

The fit is faithful whenever the fitted variant IS the generating variant (recovered/true =
0.985-0.988 in every case), and it collapses otherwise. The observed ladder — the exponent is
**largest** against the smoothed rate and smallest against the raw one — is reproduced only by
"the counts follow the SMOOTH rate and R1 is R7 plus differentiation noise" (|miss| 0.21 vs 3.2
for R3 and 6.0 for R1). Independent confirmation:
- **Errors-in-variables identity.** Reliability ratio lambda1 = Theil-Sen slope of log10 R7 on
  log10 R1 = **+0.170 [+0.101, +0.319]**; predicted b(C,R1) = b(C,R7)*lambda1 = +0.063 [+0.031,
  +0.113] vs observed +0.108 (observed - predicted +0.044, p = 0.12 — consistent). For R3,
  lambda3 = +0.382, predicted +0.126 vs observed +0.198 (difference +0.082, p = 0.006 — the 3-d
  rate carries a little real information beyond pure attenuation).
- **Rate character.** Lag-1 autocorrelation of log10 R: **R1 +0.061**, R3 +0.391, **R7 +0.661**;
  corr(log10 R1, log10 R7) = +0.208; peak amplitudes 11.79 (R1) / 5.19 (R3) / 3.03 (R7) mm/d, a
  factor 3.9 between R1 and R7. A_raw comes from TRAILING-30-DAY GNSS displacements, so the model
  cannot support a daily slip rate at all: R1 is essentially white differentiation noise.
Consequence: **b(C,R7) is an essentially unbiased estimate** of the exponent with respect to the
only slip rate the model can resolve (de-biased +0.401/0.988 = +0.406), and the R1/R3/Rcum
exponents are attenuated versions of it, not competing measurements.

### VERDICT — the sub-linearity is real, not a measurement artefact
Every artefact channel was tested and none of them produces it:
1. **Count saturation: ruled out.** The 576-window/day ceiling is approached on 1.5 % of days,
   removing those days moves b by 0.00 (p = 0.41), the upper half of the rate range gives the
   same-or-larger exponent (+0.10, p = 0.35), and multiplicity is negligible (b(C vs U) = 1.03).
2. **Network sensitivity: ruled out.** No epoch dependence (p = 0.94), no dependence on num_stas
   (p = 0.20), and the >= 8-station subset gives a slightly LARGER exponent (p = 0.012).
3. **Duration/occupancy measures agree with the count**: b_U +0.377 vs b_C +0.401, b_Uloc = b_C.
4. **The amplitude-weighted measure agrees and is still sub-linear**: b_E = +0.644, CI [+0.33,
   +0.76], 14/18 events below 1; Mmax +0.15 mag/decade.
5. **Smoothing: settled, and it does not rescue linearity.** The observed dependence on the
   smoothing choice is errors-in-variables attenuation of a relation to the SMOOTH rate, not a
   destroyed relation to the raw rate; fitting against R7 recovers the true exponent to 1 %.
6. **The one bias that survives strengthens the result.** If the 13-pattern vbICA model COMPRESSES
   the log dynamic range of the true slip rate (log R_model = alpha * log R_true with alpha < 1),
   then b_obs = b_true/alpha > b_true — model compression makes the true exponent even smaller.
   Only an EXPANDED model dynamic range could inflate a linear relation into an apparent 0.4, and
   a smooth 13-pattern mixture of trailing-30-day amplitudes does not expand dynamic range.

**Best estimate.** Tremor response to modeled slip rate is **sub-linear with exponent b ~ 0.4 in
detection counts and ~0.65 in radiated energy**. Quoting the energy exponent as the physical one
(counts miss the +0.23 of the response that goes into larger detections) the best estimate is
**b = 0.65, honest range 0.3-0.8** (90 % bootstrap CI of the median +0.33 to +0.76; per-event IQR
+0.17 to +0.84; the count-based number +0.40 [+0.31, +0.45] is a lower bound on it, and the
pooled event-dummy fits give +0.40 for C and +0.62 for E). **b = 1 is excluded** for both measures
(27/27 and 14/18 events below 1; CIs exclude 1). A 10x rise in modeled slip rate buys ~2.5x more
detections and ~4.5x more tremor energy — a factor 2 short of proportionality.

### Caveats specific to S11
- **n = 18 for every energy/magnitude number** (the catalogue's amplitude columns start
  2017-01-01), and 6 of those events are the short/thin ones already flagged in S1/S4/S9 (1905
  8 days b = 1.43, 1434 8 days b = -0.04, 2015 and 1127 with b ~ 0). They carry the tails.
- **The exponent is measured against a MODEL slip rate**, not against slip. Everything above
  bounds the tremor-side artefacts; it cannot bound an error in the vbICA rate itself beyond the
  compression argument in point 6.
- **Rcum is not usable** as a slip-rate variant here (48 % zeros, 2.9 dex range, fails the null);
  it is reported for completeness only.
- **The 2017+ null uses its own wrap window** (3434 d instead of 5991 d), so those nulls have
  fewer independent shifts available than the pre-2017 ones.
- The lower/upper-half tests use >= 6 days per half, so 2 of the 27 events drop out of them.
- `EperDet` and `CperU` are ratios of two measured quantities on the same days; their exponents
  are exactly b_E - b_C and b_C - b_U to within the eps convention, and are reported as
  decompositions, not as independent evidence.

Files written: `scripts/s11_scaling_robustness.py`, `results/s11_exponents.csv` (144 rows =
9 measures x 4 slip variants x 4 day-subsets), `results/s11_saturation.csv` (27 rows),
`results/s11_energy_events.csv` (18 rows), `results/s11_smoothing_sim.csv` (36 rows = 27
simulation cells + 9 ladder rows), `figures/s11_scaling_robustness.png` (A exponent heatmap,
B C vs U with the 576 cap, C lower- vs upper-half exponents, D energy vs count against R7,
E the ladder discriminant, F observed exponents vs the shifted-tremor null). Nothing outside
`tremor_slip_geometry/` was created or modified.

## 2026-09-10 orchestrator review — S11
Sub-linear tremor-slip scaling SURVIVES every artefact test: no window saturation (max U 546/576, 1.5 % of days >= 500; removing
them changes b by 0.000), counts track window occupancy 1:1 (b(C,U) 1.03), upper half of the slip-rate range gives the same or
larger exponent (opposite of saturation), no epoch or station-count dependence (>= 8-station subset RAISES b). Energy (2017+, 18
events): b_E = 0.64 vs b_C = 0.37 on the same events (p 0.003) — energy per detection rises with slip rate, so counts under-
represent strong tremor by ~0.23 in exponent, yet energy is still sub-linear; b = 1 excluded (27/27 events below 1 in counts,
14/18 in energy). Smoothing caveat CORRECTED (S9 had it backwards): the observed exponent ladder is reproduced only if tremor is
generated by the SMOOTH rate; b vs R7 is unbiased to 1 %, and the low exponents vs unsmoothed rates are attenuation. Best estimate
b ~ 0.65 (energy), honest range 0.3-0.8; counts 0.40 [0.31, 0.45] = lower bound. Remaining caveat: the slip rate is a MODEL
rate (any dynamic-range compression by the 13-pattern model would push b_obs above b_true, i.e. strengthens sub-linearity).
Physical reading: a 10x higher slip rate yields ~4.5x more tremor energy — tremor productivity per unit slip FALLS as slip
accelerates, the same everywhere and in every epoch. This is now the strongest candidate result of the study.

## 2026-09-10 S8 big events + predictability

Ran `scripts/s8_bigevents_predictability.py` (PLAN.md round 2, S8 = Q7 + Q18) end to end,
OMP_NUM_THREADS=4, seed 20260910, ~365 s. Model + PNSN tremor only; no dv/v. Nothing outside
`tremor_slip_geometry/` was created or modified.
Files: `results/s8_bigevents_before_after.csv` (56 rows), `results/s8_bigevents_events.csv` (127),
`results/s8_predictability.csv` (72), `results/s8_predictability_controls.csv`,
`results/s8_growth_check.csv`; `figures/s8_bigevents.png`, `figures/s8_predictability.png`.

### PART A — before vs after the three largest northern events (Q7)
Reference events (the three largest confirmed N/C patches of `results/s1_events.csv`):
1317 2015-12-28..2016-03-21 Mw 6.79 (footprint 310 triangles), 1863 2017-02-24..05-04 Mw 6.63 (177),
1605 2024-08-14..10-01 Mw 6.58 (166). SITE = the reference footprint; an event joins when
|F_E ∩ F_R| / min(|F_E|,|F_R|) >= 0.30 (0.60 also run as a sensitivity). BEFORE = start < R.start,
AFTER = start > R.end (R itself falls in neither). NULL = 200 random reference dates uniform
2012-2024, one per site per draw, same sites and the same BEFORE/AFTER gap; only the split date moves.

**FIRST FINDING, before any test: the 30 % overlap rule makes all three "sites" essentially the whole
northern margin.** The three sites contain 27 / 18 / 21 of the 27 confirmed events (BEFORE/AFTER
8/18, 5/12, 18/2). The small event footprints are nested inside the big chained ones, so the overlap
criterion is nearly vacuous and the three "sites" are largely the same event list. Ref 1605 (2024)
has only 2 AFTER events. This is a weak test by construction; read it as a bound, not a search.

Pooled (stacked: each event once per site it belongs to; 31 BEFORE / 32 AFTER instances), overlap >= 0.30:

| quantity | med BEFORE | med AFTER | diff | p_MWU | null sd | z | p_null |
|---|---|---|---|---|---|---|---|
| Q1 down-dip offset (km)        |   5.065 |   4.112 |  -0.953 | 0.545 | 0.735 | -1.30 | 0.612 |
| Q3 raw depth offset (km)       |   1.190 |   0.806 |  -0.384 | 0.221 | 0.090 | -4.24 | 0.318 |
| tremor depth - margin-wide (km)|   0.000 |  -0.005 |  -0.005 | 0.434 | 0.006 | -0.78 | 0.652 |
| tremor-days(>=3)/potency (d/km3)| 120.27 |  115.54 |  -4.73  | 0.773 | 16.16 | -0.29 | 0.935 |
| potency (km3)                  |   0.204 |   0.127 |  -0.077 | 0.191 | 0.058 | -1.33 | 0.667 |
| duration (d)                   |  42.0   |  34.5   |  -7.5   | 0.036 |  5.01 | -1.50 | 0.458 |
| recurrence to previous (d)     | 268.0   | 183.0   | -85.0   | 0.269 | 53.76 | -1.58 | 0.617 |

Per reference the only p_MWU < 0.05 anywhere are ref 1863 duration 52 -> 33 d (p 0.0051, p_null 0.153)
and ref 1863 potency 0.299 -> 0.166 km3 (p 0.048, p_null 0.199), plus the pooled duration row above.
28 tests are run per overlap threshold; none survives Holm and none survives the random-reference-date
null. The 0.60-overlap sensitivity gives the same picture (pooled duration p 0.086).

Two methodological points that matter for reading the CSV:
1. **Read `p_null`, not `z`.** The null distribution of a DIFFERENCE OF MEDIANS over 5-18 events is
   discrete and heavy-tailed (the median hops between data values), so its sd understates the tail:
   Q3_depth_km_raw shows |z| ~ 4.2-4.4 at p_null 0.12-0.32. `p_null` is the rank-based p against the
   same 200 draws and is the statistic to quote.
2. **The margin-wide depth reference is nearly degenerate here.** New diagnostic column
   `frac_marginwide_in_near` = the fraction of the concurrent margin-wide detections that are ALREADY
   inside the 40-km near-patch set: median 0.976 (north) and 0.999 (central). So
   `Q3_rel_marginwide_km` is ~0 BY CONSTRUCTION for almost every event — which is itself the S1b
   conclusion made quantitative: **there is no LOCAL tremor-depth anomaly at any confirmed patch; the
   entire Q3 signal of S1/S1b is the depth of the whole margin's tremor during the days the 5 mm/30 d
   threshold selected.** Only three events have a genuinely independent reference
   (frac < 0.8): 1448 (0.43, +0.41 km), 1430 (0.61, -3.67 km), 1873 (0.78, +0.48 km). Event 1434
   (central, 2-day window) has NO 44-46N detection with a finite Slab2 depth in its window — its
   40-km neighbourhood actually reaches north of 46N — so its Q3_rel is NaN.
Integrity check inside the script: the recomputed near-patch tremor median depth reproduces
`s1_events.csv`'s `tremor_depth_med` to 7e-15 km for all 27 events.

### PART B — time- vs slip-predictable recurrence (Q18)
SITES = the 12 S5 bands with 44.0 <= band < 50.0 (all have >= 6 episodes; 13-20 each), 207 episodes,
195 T_prev and 195 T_next pairs, peak-to-peak as in S5. Size proxies: `count30` raw; `count30`
divided by a 5-yr running median of the band's own count30; and the S5-matched slip patch's Mw and
potency (Mw is a monotone function of potency, so their Spearman rhos are identical by construction).
Pooled = within-site fractional ranks; p from 2000 within-site permutations.

Network growth: **no band's episode size trends with date** (12 bands, all p >= 0.22; pooled per-band
z-scores rho +0.006, p 0.94). The 5-yr detrend therefore changes little.

Pooled (n = 195 for count30, 141/144 for the patch proxies):

| proxy | rho(size, T_prev) | p | rho(size, T_next) | p |
|---|---|---|---|---|
| count30 raw        | +0.160 | 0.030 | **+0.330** | **0.0005** |
| count30 detrended  | +0.136 | 0.071 | **+0.329** | **0.0005** |
| matched patch Mw   | +0.031 | 0.733 | +0.075 | 0.394 |
| matched patch potency | +0.031 | 0.735 | +0.075 | 0.392 |

Segments: north rho_next +0.337 raw / +0.356 detrended (p 0.0005 / 0.0005), rho_prev +0.209 / +0.168
(p 0.016 / 0.058); central rho_next +0.316 / +0.268 (p 0.012 / 0.040), rho_prev +0.066 / +0.072 (ns).
Per site: nothing survives Holm over the 96 site x proxy x direction tests (best raw p: band 45.5
T_next +0.665 p 0.006, band 46.5 T_next detr +0.616 p 0.008, band 48.0 T_next detr +0.580 p 0.016).

**Artefact controls on the T_next result** (added because a positive size~T_next can be manufactured
by the episode DETECTOR: S5 picks peaks greedily highest-first with a >= 90 d exclusion, so a big peak
wins its slot and mechanically stretches the interval after it):
- C1 interval clustering, rho(T_prev, T_next) pooled = **+0.118, p 0.12** (only band 49.0 significant,
  +0.516) — short intervals do not strongly cluster.
- C2 rank-partial rho(size, T_next | T_prev) = **+0.305 raw / +0.307 detrended, p < 0.001** — not
  routed through T_prev.
- C3 using only intervals >= 150 d (n 179): **+0.272 / +0.279, p 0.0005**.
- C4 rebuild the episode catalogue from the band trailing counts at PK_SEP = 45 / 90 / 150 d
  (the 90 d rebuild reproduces all 207 S5 episode keys exactly, asserted in code):
  rho_next detrended = **+0.344 / +0.329 / +0.234**. It DECREASES as the exclusion radius grows —
  the opposite of what the greedy-exclusion artefact predicts. Cleared.
- C5 bootstrap over the 12 bands, rho_next - rho_prev = **+0.193, 95 % CI [-0.042, +0.426]** (raw:
  +0.171 [-0.015, +0.361]) — the two directions are NOT separable.

**Power (the interpretation guard).** With ~14-month recurrence and ~16 yr of catalogue a band has
13-20 episodes, i.e. 11-19 pairs: at n = 16 only **|rho| >= 0.65** is detectable at 80 % power
(alpha 0.05 two-sided). Per-site tests are therefore essentially uninformative — no per-site claim,
positive or negative, should be made. Only the pooled test has power, and even it can only see
|rho| >= 0.20 (n 195) / 0.23 (n 144, patch proxies).

### Plain answers
**(i) Did the 2015-16 / 2017 / 2024 big northern events change subsequent events at their sites?
NO — nothing measurable, on any of the seven quantities.** Down-dip offset, raw depth offset,
margin-referenced tremor depth, tremor-days per unit potency, potency, duration and site recurrence
are all indistinguishable BEFORE vs AFTER, both per reference and pooled, and every nominal p < 0.05
(ref 1863 duration 52 -> 33 d, p 0.005) is matched by the random-reference-date null (p_null 0.15),
i.e. it is the record-wide drift toward shorter/smaller confirmed events, not an effect of the big
event. Honest strength: this is a **bound**, and a loose one — the 30 % overlap rule makes all three
sites nearly the whole northern margin, ref 1605 has only 2 AFTER events, and the quantities are
event-averaged over ~20 events. The clearest positive statement to come out of Part A is the
degeneracy diagnostic: once tremor depth is referenced to the concurrent margin-wide tremor depth,
the depth offset at every confirmed patch is ~0, so there is no local tremor-depth signal to change.

**(ii) Is ETS recurrence here time-predictable, slip-predictable, or neither? NEITHER cleanly; the
only signal is a weak, margin-wide, SLIP-predictable TENDENCY in the tremor-count proxy.** Pooled
across 12 bands, episode size correlates with the FOLLOWING interval at rho +0.33 (p 0.0005, ~11 % of
rank variance) — bigger ETS, longer wait to the next one — and this survives detrending, partialling
out T_prev, dropping intervals < 150 d, and the de-clustering-radius sweep. The time-predictable
direction (size ~ preceding interval) is weaker and marginal (+0.14 to +0.16, p 0.03-0.07). But three
things stop this from being a result: (a) the independent MODEL-SLIP size proxy (matched patch
Mw/potency) shows **nothing** at comparable power (rho_next +0.075, p 0.39, detectable |rho| 0.23), so
the relation lives entirely in the tremor count; (b) the bootstrap cannot separate rho_next from
rho_prev (CI [-0.04, +0.43]); (c) no individual band is significant and no band CAN be at n ~ 16.
So: **not per-site predictable in any usable sense; at the margin scale a weak lean toward
slip-predictable, unconfirmed by the slip model.**

### Caveats
- **n.** Part A: 27 confirmed events total, 31/32 pooled instances built from them by counting an
  event once per reference site (so the pooled rows are NOT independent observations — the same
  events recur across the three sites; treat the pooled p as descriptive). Part B: 195 pairs but
  only 11-19 per site.
- **Frame.** Slip start/end are 30-d trailing threshold crossings, so every Part A date and every
  `recur_prev_d` is late by 0-30 d. The bias is constant so the BEFORE/AFTER split is unaffected, but
  individual recurrence intervals carry +-30 d of frame noise. Part B uses tremor peaks only and is
  frame-free.
- **Completeness.** `count30` is a raw detection count; handled by the 5-yr running-median detrend
  and checked by the growth test (no band trends, pooled rho +0.006). `tremordays3_per_potency` uses
  tremor-DAYS with >= 3 detections precisely to blunt the S3b detection-density channel, but S3b
  showed tremor-days still carry a threshold effect, so it is hardened, not clean.
- **Truncation.** The tremor catalogue ends 2026-05-27, so the last episode of every band has no
  T_next and the 2025 events have few successors; ref 1605's AFTER group is 2 events for this reason.
- **Multiple testing.** Part A runs 7 quantities x 4 reference sets x 2 thresholds; Part B runs
  12 sites x 4 proxies x 2 directions. Nothing in either block survives Holm within its family; only
  the POOLED Part B T_next rows (p 0.0005 = the permutation floor at 2000 draws) would.

## 2026-09-10 S10 literature check
Ran as a read-only literature sub-agent (WebSearch/WebFetch: publisher sites, PMC, Semantic Scholar API,
author/institutional PDFs). Full report with tables, DOIs and URLs: `results/s10_literature_check.md`.
Only that file and this section were written. RESULTS.md was updated with S11/S12 mid-task; the verdicts
for candidates 1 and 2 and the ranking were revised against the new RESULTS.md.

Papers read IN FULL (text extracted): Bartlow et al. 2011 GRL; Brudzinski & Allen 2007 Geology;
Wells et al. 2017 Geology; Wech & Creager 2011 Nat. Geosci.; Gualandi 2025 GJI. Abstract/summary only
(AGU and Springer returned 403 throughout): Meng et al. 2026 GRL; Costantino et al. 2026 GRL;
Hall et al. 2018 and 2019 G-cubed; Bartlow 2020 GRL; Wech & Bartlow 2014 GRL; Schmidt & Gao 2010 JGR;
Aguiar et al. 2009 JGR; Michel et al. 2019a PAGEOPH; Sweet et al. 2019; Bombardier et al. 2024;
Zhang & Houston 2025; Wech 2021. Read via PMC: Frank & Brodsky 2019 Sci. Adv.; Passarelli et al. 2021
Sci. Adv.; Bletery & Nocquet 2020 Nat. Commun.; Luo & Liu 2021 Nat. Commun.

**Verdicts.**
- **C1 / R9 sub-linear tremor-vs-slip-rate exponent — B (quantitative extension).** No published exponent
  exists for daily tremor count (or energy) vs slip rate within events; the field's default is
  proportionality, stated as an assumption (quoted as such by Passarelli et al. 2021). Closest priors:
  Frank & Brodsky 2019 Sci. Adv. eaaw9386 Guerrero, Mdot_geo = alpha*A^(3.17+-0.29) with A = median LFE
  displacement amplitude, i.e. A ~ Rate^0.315; Passarelli et al. 2021 Sci. Adv. abg9718,
  log M_geod = 15.3 + 0.20(0.17-0.24) log M_seis (0.29 deep-only) between EVENT-integrated moments;
  Bletery & Nocquet 2020, daily tremor count corr 0.63 with moment rate, 0.76 with daily max slip rate
  (no exponent); Hawthorne & Rubin 2013, strain rate ~ tremor envelope amplitude, LINEAR but only over a
  1.5-2x amplitude range; Aguiar et al. 2009, "highly LINEAR" GPS moment vs total tremor duration over 23
  events. KEY NUMERICAL POINT for the writeup: our energy exponent 0.64 implies summed amplitude
  ~ Rate^0.32, numerically identical to Frank & Brodsky's 1/3.17 = 0.315 — but their A is a PER-LFE median,
  and our per-detection energy exponent is 0.64-0.40 = 0.24 (per-detection amplitude ~ Rate^0.12). The two
  readings disagree by ~2.6x; match the aggregation before putting 0.32 anywhere near an abstract.
  S11 already retired the smoothing/saturation/network blockers. Remaining blocker: repeat the exponent
  against a slip model that is NOT the 13 fixed vbICA patterns (Bartlow 2020 or Costantino et al. 2026).
- **C2 / R11 central-Oregon tremor without slip — B, upgraded from C by S12.** The resolution test passing
  (98.6 % recovery, threshold 0.91x north) answers Schmidt & Gao 2010's explicit caveat ("limited station
  coverage south of Portland (45.5N) restricted the analysis to events on the northern half"). No prior
  report of tremor-WITHOUT-slip as a persistent localised Cascadia feature; the published asymmetry runs
  the other way (Wech & Bartlow 2014; Hall et al. 2018; Wells et al. 2017 explicitly). BUT two published
  results must be answered before this is claimable: (i) **Bartlow 2020 GRL** resolves a time-averaged ETS
  slip rate over the whole margin INCLUDING Oregon and says slip correlates well with tremor locations —
  now the highest-priority PDF to obtain; (ii) **Brudzinski & Allen 2007** report GPS slow-slip episodes at
  central-Oregon stations (COR, Siletzia zone, 19+-3 mo). Also: DROP the word "tremor-productive" — Wells
  et al. 2017 ("intensity is lowest in the central Oregon forearc", Corvallis/Gales Creek low-tremor
  transverse zones, Siletzia) and Wech 2021 ("tremor in central Cascadia is systematically weaker and less
  frequent") establish the opposite. Correct framing: tremor-PRESENT, slip-ABSENT, and a test of Wells
  et al.'s fluid-escape model.
- **C3 / R8 within-event deepening — A for the statistic, but it CONFLICTS with a Nature Geoscience paper,
  and the intensity half is probably C.** Wech & Creager 2011 Nat. Geosci. 4, 624-628 (doi:10.1038/ngeo1215)
  plot "median daily downdip distance" through four ETS episodes — essentially our statistic — and conclude
  events initiate ~80 km downdip and MIGRATE UPDIP; they also give recurrence decreasing exponentially with
  depth (15 mo updip -> ~2 mo downdip) and minor swarms <20 h concentrated downdip. Luo & Liu 2021 and
  Peng et al. 2015 (2011 ETS: initiated at the downdip edge, migrated updip over 5 days) say the same.
  Most plausible reconciliation: their result is the onset ramp in a restricted N-Washington wedge, ours
  integrates the whole window margin-wide and is dominated by the tail — SPLIT THE PHASE AXIS into onset
  ramp vs tail and report both. "Shallower when most intense" is very likely Wech 2021's SPATIAL
  depth-amplitude gradient ("a clear shift to higher magnitudes at shallower portions … higher proportion
  of MeL > 1.25 events updip") re-expressed in time; drop it unless it survives a null that shuffles daily
  depth within the same day's footprint. Also needs: fixed station subset + MeL amplitudes (deep tremor is
  weaker and more continuous, so an improving detector deepens the late median for free), and a control for
  Slab2 depth varying along strike as the event propagates.
- **C4 / R3-R5-R7 stationarity — A (null), publishable only if aimed at a named claim.** The only verified
  published secular-evolution claim for the Cascadia ETS zone is **Meng, Liu, Adhikari & Luo 2026 GRL
  e2025GL120690**: "a long-term change in the locations where ETS events nucleate and their propagation
  directions" plus an incongruent tremor seasonality. ABSTRACT-ONLY (AGU 403) — get the full text; our
  R3/R7 do NOT currently test nucleation latitude or propagation direction, which is what they claim.
  No published claim exists of secular drift in the tremor-slip RELATIVE GEOMETRY. Michel et al. 2019b
  (M0 ~ T^3) and Costantino et al. 2026 ("scaling laws are strongly influenced by the chosen detection
  threshold and the SNR" — reinforces our X5) are scaling, not secular. Needs a power analysis (inject
  synthetic drifts, report the smallest detectable slope) and the bounds expressed physically (<2 km over
  16 yr = inside one mesh cell and inside the PNSN depth error).
- **C5 / Gualandi 2025 — leaves essentially everything we measure unmeasured.** Read in full. Tremor
  appears only as black dots on a latitude-time map (Fig. 7) plus one number: "The correlation with the
  maximum slip potency is typically of the order of 0.7 over the all Cascadia subduction. The time lags
  observed are typically between 0 and a few days, with positive values indicating that the slip potency
  rate precedes the tremors." NO co-location offsets in km, NO depth comparison, NO scaling, NO per-event
  analysis, NO resolution/checkerboard test, NO Oregon discussion (all verified by full-text search).
  He disclaims novelty on the correlation (cites Wech et al. 2009, Bartlow et al. 2011, Michel et al.
  2019a, Bartlow 2020) and flags his own lag as possibly a filter artefact: "further analysis on the
  effect of the adopted filter are needed"; also "more aggressive filters … may anticipate the rising of
  the slip rate for past events". That caveat bears directly on our R4 and on S11's smoothing test.
  Costantino et al. 2026 GRL e2025GL117446 (Costantino, Radiguet, El Yousfi, Socquet): 2007-2023 denoised
  GNSS, ST-DBSCAN catalogues at 4 thresholds, upper bound on slow-slip moment rate, 2-D expanding crack
  propagating laterally at the along-dip limits of the transition zone. No tremor comparison in the
  abstract; full text 403, so not settled (their 2023 Comm. Earth Environ. paper DID validate against
  tremor).

**Unsolicited finding — R10 is largely already published.** Hall, Houston & Schmidt 2018 G-cubed 19,
2706-2718 (doi:10.1029/2018GC007694), 6 ETS 2007-2016, 80 GPS: "slow slip of up to 2 cm extends updip of
tremor by about 15 km beneath the Olympic Peninsula … Slip updip of tremor is a persistent feature of all
six ETS events", with along-strike variation tied to overlying crustal lithology. Peng & Rubin 2017 restate
it as 10-20 km. Our 11 km over 27 events / 16 yr is an EXTENSION, not a discovery — cite Hall et al. 2018.
The new part of R10 is the no-drift bound on the edges (0.03-0.6 km/yr), which belongs with C4.

**Also worth knowing (bears on several results).** Hall, Schmidt & Houston 2019 G-cubed
(doi:10.1029/2019GC008510) find PEAK TREMOR RATES LEAD PEAK SLIP RATES in the 2010 and 2012 Mw 6.8 events,
independent of propagation direction, and "ongoing slip continuing several days after tremor has passed" —
the opposite sign to Gualandi's "slip precedes tremor". Our R4 sits between the two and should cite both.
PNSN catalogue is non-stationary by construction (algorithm updated Feb 2018; Wech 2021 reports +55 %
epicentres with the new method) — S11 controlled for this, keep the citation.

**Ranked publishability (post-S11/S12):** 1) C1 exponent — GRL stand-alone or JGR-SE headline, closest to
done, one blocker (model independence). 2) C2 central Oregon — GRL if Bartlow 2020 is answered. 3) C4
stationarity — GRL only if aimed at Meng et al. 2026, else a JGR-SE section. 4) C3 deepening — GRL if it
survives, else cut; least controlled of the four. 5) C5 — the gap statement for the introduction, not a
result. Recommended package: one JGR-SE paper with C1 leading, C2 second, C4 as the stationarity section,
C3 held back; or a fast GRL on C1 alone.

**Could not verify online (must be obtained before writing):** Meng et al. 2026 full text; Costantino et al.
2026 full text; **Bartlow 2020 full text and slip map (highest priority — decisive for C2)**; Hall et al.
2018/2019 full texts; Wech & Bartlow 2014 full text (no slip-rate threshold value obtained); Aguiar et al.
2009 full text (no regression slope/R^2 obtained); Michel et al. 2019a full text (the "10-20 km halo"
checkerboard statement is second-hand).

## 2026-09-10 orchestrator review — S8
Q7: no detectable change in any of seven quantities after the 2016 / 2017 / 2024 events (pooled 31 before / 32 after; all
p_null >= 0.3; nominal duration/potency hits for 2017 match the record-wide drift toward shorter confirmed events). The 30 %-
overlap "site" rule is vacuous because small footprints nest inside the big chained ones -> the test is margin-wide, not local;
that is the best this catalog allows. Degeneracy quantified: the concurrent margin-wide tremor reference shares 98-100 % of its
detections with the near-patch set -> there is NO local tremor-depth anomaly at any confirmed patch (S1b made rigorous).
Q18: recurrence is neither cleanly time- nor slip-predictable. Weak slip-predictable TENDENCY in the tremor proxy: episode
tremor count vs time-to-next rho +0.33 (p 0.0005, 195 pairs, 12 bands), vs time-since-previous +0.16 (p 0.03); survives
detrending, partialling, interval floors and re-detection radii. But the model-slip size proxy shows nothing (rho +0.08, power
to detect 0.23), per-site tests are powerless (n ~ 16), and bootstrap over bands cannot separate rho_next from rho_prev. Keep as
a candidate to test with an independent size measure (tremor energy 2017+, or a NIF catalog), not a claim.

## 2026-09-10 orchestrator review — S10 literature check; round 4 launched
Ranking from the literature agent: (1) sub-linear tremor-vs-slip-rate exponent — no published exponent for daily count/energy
vs slip rate within events; the field assumes proportionality; closest: Frank & Brodsky 2019 (per-LFE amplitude ~ rate^0.32,
Guerrero) — our energy exponent 0.64 implies summed amplitude ~ rate^0.32 (coincidence to be checked against their per-event
aggregation); Gualandi 2025 itself has no co-location, no scaling, no per-event analysis (tremor = black dots + one correlation
sentence). Blocker: model dependence -> S11b. (2) central Oregon tremor-without-resolvable-slip — upgraded to B by S12; must
answer Bartlow 2020 (margin-wide average ETS slip incl. Oregon) and Brudzinski & Allen 2007 (GPS episodes at central-Oregon
stations); DROP "tremor-productive" (Wells 2017, Wech 2021: central Oregon is the least tremor-productive segment). (3) the
within-event deepening CONFLICTS with Wech & Creager 2011 (up-dip migration through four ETS) -> S9b reconciliation; the
"shallower when intense" half is probably Wech 2021's spatial gradient. (4) stationarity bounds publishable only if aimed at Meng
et al. 2026 GRL (claimed long-term change of nucleation locations and propagation directions) -> S13 tests their quantities +
power analysis. R10 (slip up-dip of tremor by ~15 km) is already Hall, Houston & Schmidt 2018 — cite; our new part is the no-
drift bound. Hall et al. 2019: tremor peak rate LEADS slip peak rate (opposite sign to Gualandi) -> cite both for R4.

## 2026-09-10 S12b temporal loophole

Ran `scripts/s12b_temporal_loophole.py` end to end, `OMP_NUM_THREADS=4`, seed 20260910, ~4 min.
Follow-up to S12's caveat 2: S12 tested the linear step (ICA spatial response -> slip) but not
the vbICA step, leaving open the possibility that a central-Oregon transient with its own time
history has no component to live in. Outputs `results/s12b_{component_tremor_corr,
slipseries_corr,matched_filter,data_vs_model,stations}.csv`,
`figures/s12b_temporal_loophole.png`. Nothing outside `tremor_slip_geometry/` touched.

Setup: solution 2026-08-18 (P = -Ld and `component_amplitudes()` from the SAME solution — the
shipped `components_at_cells_trailing30.npz` `A_raw` comes from an earlier snapshot and
correlates only 0.87-0.99 with it, so it was NOT used). Timeline 7169 d; analysis window
2010-05-14..2026-03-28 (5798 d). Tremor per set = detections within 40 km of a set centroid,
30-d trailing sum, detrended by a 365-d centred rolling median (kills the PNSN growth ramp).
Slip rate = 30-d trailing increment of each amplitude. Sets: CLUSTER (the 9 S7 triangles,
27,275 detections), CENTRAL_ST (41), CENTRAL_ALL (144), NORTH_ST (145, 207,605 detections).
Spearman(cluster tremor, north tremor) = +0.100, so the two ETS series are largely independent
and the north/central comparison is a real contrast, not the same signal twice.

### T1 — the 13 amplitude functions DO carry central-Oregon timing
Best-lag |Pearson| of each component's 30-d amplitude rate against each tremor series:
**9 of 13 components correlate better with CENTRAL than with NORTH tremor.** Sign-corrected by
the component's own thrust sense at the cluster (+ = adds thrust when tremor is high):
  comp 11  r_signed +0.573 @ lag -2 d   (|r| north 0.536)
  comp 26  r_signed +0.358 @ lag  0 d   (|r| north 0.089)
  comp 24  r_signed +0.234 @ lag -4 d   (|r| north 0.116)
  comp 20  r_signed +0.203 @ lag -1 d   (|r| north 0.156)
  comp 10  r_signed +0.135 @ lag -25 d  (|r| north 0.502 — this is the northern ETS component)
So the basis is not temporally blind to central Oregon: several components fire in phase with
central-Oregon tremor and deliver thrust to the cluster while doing it (0.06-0.86 mm of thrust
per 1 sd of amplitude). The loophole is closed at the component level.

### T2 — the modeled slip ON the cluster tracks the cluster's own tremor
Spearman(modeled thrust potency rate, own 30-d tremor), 200 rigid tremor shifts >= 365 d as null:
  CLUSTER      +0.186  null +0.003 +- 0.045  p 0.0050   best-lag Pearson +0.508 @ -3 d
  CENTRAL_ST   +0.169  null +0.001 +- 0.044  p 0.0050   best-lag Pearson +0.488 @ -3 d
  CENTRAL_ALL  +0.101  null +0.003 +- 0.037  p 0.0149   best-lag Pearson +0.307 @ -3 d
  NORTH_ST     +0.380  null -0.004 +- 0.056  p 0.0050   best-lag Pearson +0.653 @ -2 d
The cluster's modeled slip follows its tremor at least as well as its SLIP+TREMOR neighbours'
does. **But the model has almost no freedom to do otherwise**: corr(CLUSTER, CENTRAL_ST modeled
potency-rate series) = **+0.965** (vs +0.521 against NORTH_ST). With 13 fixed patterns the
cluster and its neighbours are given essentially the same time history by construction, so T2
alone is weak evidence — its real content is the AMPLITUDE: potency-rate sd per km^2 is
**1.03e-3 (CLUSTER) vs 2.05e-3 (CENTRAL_ST) and 2.25e-3 (NORTH_ST)** — the cluster gets exactly
HALF the slip-rate amplitude per unit area of its slip+tremor neighbours, at the same timing.
**The S7 decoupling is spatial/amplitude, not temporal.**

### T3 — data-space matched filter, no patterns involved
sigma=30 km unit thrust patches on tri 252 (44.67N, -123.42, 34.5 km; 5 nearest stations P375
3 km, CORV 13, P374 34, P376 39, P378 42) and tri 524 (47.92N, -123.42, 36.5 km; P435 17 km,
PTAA 23, SC03 24, P436 25, CHCM 49). Green's-function cross-talk cos(d_cluster, d_north) =
**-0.003**, so the two patches are essentially orthogonal over the 213-station network and the
JOINT fit equals the single-patch fit (the cross-talk worry was unfounded, but the joint fit is
reported anyway). Channel noise on the detrended 30-d increment: raw GNSS 1.54 mm, ICA-13comp
0.73 mm, model-predicted 0.46 mm.
Spearman(matched-filter amplitude, own tremor), raw GNSS, JOINT: **CLUSTER +0.219, NORTH +0.379**
(both p 0.005 against the shifted null, null sd 0.05).
Composite = mean estimated 30-d peak slip in the top tremor decile minus the bottom decile,
se from the number of independent 30-d windows. The SAME filter is applied to (a) the raw GNSS,
(b) the 13-component ICA reconstruction, (c) the slip model's own predicted displacement
`G L a(t)` — so data and model suffer identical 30-km smearing and the difference is meaningful:
                       raw GNSS      13-pattern model   data - model        model captures
  CLUSTER  (JOINT)   10.2 +- 2.7 mm   6.4 +- 1.7 mm   **+3.7 +- 3.2 (1.2 sig)**    63 %
  NORTH_ST (JOINT)   27.6 +- 4.8 mm  21.4 +- 3.9 mm   **+6.2 +- 6.2 (1.0 sig)**    77 %
  CLUSTER  (5 stn)    6.3 +- 2.4 mm   3.8 +- 1.6 mm     +2.6 +- 2.9 (0.9 sig)      60 %
  NORTH_ST (5 stn)   23.6 +- 4.7 mm  16.6 +- 3.0 mm     +6.9 +- 5.6 (1.2 sig)      71 %
The raw GNSS at the five central-Oregon stations DOES see a tremor-synchronous transient
(3.8 sigma), the model reproduces ~63 % of it, and the residual is **1.2 sigma — not
significant, and no worse than the 1.0 sigma residual at the northern control**. Formally:
**2-se bound on unmodelled tremor-synchronous slip under central Oregon = 6.3 mm of peak slip
per 30-day window on a 30-km patch (equivalent Mw 5.95)**, against an observed total of 10.2 mm.

### Caveats
1. **T3's resolution is 30 km, not 9 triangles.** The patch centred on tri 252 overlaps central
   SLIP+TREMOR triangles as well as the cluster (S12: peak retention 0.56). T3 therefore says
   "central Oregon near 44.7N has ETS slip and the model has most of it", NOT "these 9 triangles
   slip". The 9-triangle statement remains a rank statement inside the model.
2. T2's +0.965 cross-correlation means the timing test has little discriminating power by
   construction; the informative half of T2 is the 2x amplitude deficit.
3. The composite uses tremor deciles on a 30-d trailing count; the top decile is dominated by
   the big central-Oregon episodes, so the estimate is an episode-average, not a per-event value.
4. Detrending removes anything slower than ~1 yr, so a secular or multi-year central-Oregon
   slip deficit is outside this test (S12/S7 cover the long-term amplitude).
5. The station-level east/north correlations are negative (median -0.14 central, -0.21 north),
   the expected seaward ETS reversal; the matched filter carries the sign through the Green's
   functions, so the positive composite is the physically correct sense.
6. 200 null draws give a p-value floor of 0.005; every reported p sits at or near that floor.

### Verdict on the loophole
**CLOSED — the loophole is not open.** (i) The 13 amplitude functions are not northern-only:
nine of them track central-Oregon tremor better than northern tremor, and four deliver thrust to
the cluster in phase with its tremor at lags of -4 to 0 d, the strongest at r = +0.57. (ii) The
modeled slip on the cluster co-varies with the cluster's own tremor (rho +0.186, p 0.005) as
well as its SLIP+TREMOR neighbours do (+0.169), so the model is not mis-TIMING the cluster; what
it withholds is amplitude — half the slip-rate per km^2 of its neighbours. (iii) Bypassing the
patterns entirely, the raw GNSS at the five nearest stations shows a tremor-synchronous transient
of 10.2 +- 2.7 mm equivalent peak slip, of which the 13-pattern model already accounts for 6.4 mm;
the 3.7 +- 3.2 mm residual is 1.2 sigma, statistically indistinguishable from the 1.0 sigma
residual at the northern control, and bounds any unmodelled central-Oregon tremor-synchronous
slip at 6.3 mm peak per 30-d window (Mw 5.95).
**Consequence for S12/S7.** The last model-side escape route for the central-Oregon decoupling is
shut: the deficit is not a missing time function and not slip hiding in the data. Combined with
S12 (network resolves the patch, misfit flat, basis handicap ~20 %), the S7 central-Oregon
TREMOR-ONLY cluster stands as a **candidate physical result** — a stretch of interface that is
tremor-productive while its transient slip amplitude per unit area is about half that of the
slip+tremor patches on either side, at the same ETS timing. The honest residual uncertainty is
resolution, not model architecture: at 30-km smearing neither the model nor the geodesy can say
whether the deficit is confined to the 9 triangles or spread over the surrounding 44-46N band.

## 2026-09-10 orchestrator review — S12b
Temporal loophole closed, but the central-Oregon candidate is DOWNGRADED by the same test: 9 of 13 components correlate better with
central than northern tremor (comp 11 r 0.57); the cluster's modeled potency rate co-varies with its own tremor (rho 0.19, p
0.005) as well as the neighbours' do; and in DATA space the 5 nearest stations see a tremor-synchronous transient (10.2 +- 2.7 mm
composite, 3.8 sigma) of which the 13-pattern model already captures 63 % (north control 77 %), residual +3.7 +- 3.2 mm (1.2
sigma) — indistinguishable from the north. Bound on unmodelled tremor-synchronous slip: 6.3 mm per 30-d window on a 30-km patch.
So central Oregon near 44.7N does slip with its tremor; what the model withholds on the 9 triangles is AMPLITUDE (potency rate per
km^2 exactly half the adjacent slip+tremor patches), a rank statement at the resolution limit (the synthetic test resolves 30 km,
the cluster is smaller than that in the along-dip direction). Honest reading: "tremor with about half the slip per unit area of
its neighbours", not "tremor without slip". Together with the literature (central Oregon = least tremor-productive segment) this
is a modest amplitude heterogeneity, not a headline. Input note from the agent: components_at_cells_trailing30.npz A_raw is from an
EARLIER solution snapshot (corr 0.87-0.99 with 2026-08-18 amplitudes) — the whole SSE-vs-dv/v chain used that snapshot; flag for
the main project (re-derive A_raw from the current solution before any final dv/v numbers).

## 2026-09-10 S10b (literature check, follow-ups: candidate 6 + Bartlow data + Oregon quotes)
Report section appended to `results/s10_literature_check.md` under "# 2026-09-10 S10b". New file
downloaded: `results/external/results_for_paper.mat` + `results/external/PROVENANCE.md`.

**C6 recurrence predictability — VERDICT C (already published), and THE LABEL IS INVERTED.**
Shimazaki & Nakata (1980) convention, quoted verbatim from Gao et al. (2012, p. 357): "For the time
predictable model, the time to the next event is equal to the static stress drop released by the most
recent earthquake divided by the stressing rate, which implies that a larger event requires a longer
recovery time until the next one. For the slip predictable model, the fault slip released for an event is
proportional to the time interval since the last event." So size -> NEXT interval = **TIME-predictable**;
size <- PRECEDING interval = slip-predictable. Our +0.33-with-next / +0.16-with-previous is **weakly TIME
predictable**. Calling it "slip-predictable" will be caught in review.

**Prior publication (near-exact): Gao, Schmidt & Weldon II (2012), BSSA 102(1), 352-360,
doi:10.1785/0120110096.** 11 Cascadia SSEs near the Olympic Peninsula from the Schmidt & Gao (2010)
catalog, intervals from successive transient start days at GPS station ALBH, size = geodetic moment.
p.357 verbatim: "we find that the seismic moment is anticorrelated to the preevent recurrence interval for
each event with a R2 of 0.48, but has a weak positive dependence on the postevent recurrence interval with
an R2 of 0.16. This suggests that SSEs in northern Cascadia are moderately antislip predictable and weakly
time predictable." Abstract: "...weakly time predictable and moderately antislip predictable ... which may
indicate that healing between events strengthens the fault with time." Figs. 6a (slip-pred) and 6b
(time-pred). They close: "we hope to explore this more as the catalog of events grows."
=> SAME REGION, SAME SIGN, SAME MAGNITUDE, SAME NAME as our headline. C.

What is actually new in ours, and what the write-up must be built on:
(a) 195 pairs / 12 bands / 44-50N / 2010-2026 vs their 11 events at ONE station, 1998-2008 — they invited
    exactly this;
(b) a tremor-count size proxy instead of geodetic moment, with a p-value where they gave only R2;
(c) **THE INTERESTING PART — a SIGN REVERSAL on the other relation.** They found size vs PRE-event interval
    moderately ANTI-correlated (R2 0.48, r ~ -0.69) and built their whole fault-healing interpretation on
    it; we find it weakly POSITIVE (+0.16). If that reversal survives, the paper is "the anti-slip-
    predictable behaviour behind Gao et al.'s fault-healing interpretation does not survive an
    order-of-magnitude larger sample" — much better than restating the time-predictable half;
(d) an internal discrepancy to explain first: Gao used geodetic moment and DID get the weak positive
    signal; we get +0.33 from tremor count but NOTHING from the model-slip proxy.

Blockers for C6: (1) relabel to time-predictable, cite Gao 2012; (2) reproduce Gao exactly on the
overlapping region/epoch before generalising; (3) **INDEPENDENCE — biggest threat to p=0.0005**: our own
R10/S1b says "the north is ONE tremor band", so the 46-50N bands are near-duplicates and the effective n is
far below 195 — needs a block bootstrap over episodes or per-band-then-Fisher/Stouffer; (4) rule out the
episode-merge artefact (a merged large episode swallows the next one and lengthens the apparent next
interval); (5) third size proxy via Wech 2021 MeL 2017+; (6) PNSN detection-rate drift again; (7) mechanism
= Wech & Creager (2011) exponential recurrence-vs-depth under a fixed failure threshold, which PREDICTS
time-predictability.

Other checks, all negative for prior art on our specific statistic: Gualandi, Avouac, Michel & Faranda 2020
Sci Adv eaaz5548 (Cascadia 2007-2017 is low-dimensional deterministic chaos, horizon days-weeks — system
predictability, not size-vs-interval); Michel et al. 2019b Nature (M0~T^3); Michel et al. 2019a PAGEOPH;
Wech & Creager 2011 (recurrence decreases exponentially with depth, 15 mo updip -> ~2 mo downdip — the
mechanism, not the test); Brudzinski & Allen 2007; Bartlow 2020; Nuyen & Schmidt 2021; Hawthorne & Bartlow
2018 JGR 123, 4243-4265 (moment-rate SPECTRUM: flattens 1-10 d then f^-nm, nm 1.1-1.4 tremor / 0.4-1.5
strain); Beeler et al. 2013 JGR 118, 5976-5990 (LFE rate as strain-rate proxy for rheology at 400-600C);
Obara & Kato 2016 Science. Elsewhere: **Boso peninsula is reported as SLIP-predictable** (opposite pole;
Ozawa 2014 GRL, Fukuda 2018, Gardonio 2018, Kobayashi 2019 EPS); Japan renewal-process study EPS 2026
(CoV 0.25-0.54, quasi-periodic); Truttmann et al. 2024 GRL (Hikurangi chaos).

**Bartlow (2020) supporting data — DOWNLOADED OK.**
`results/external/results_for_paper.mat`, 169,711 bytes, sha256 093e48cd51bbfc9db16ab8a17b5f98b4ce09259f
cde58c4b3a45a42f199294d4 — matches the Mendeley manifest. Dataset doi:10.17632/mc49zmg7n7.1 (CC-BY-4.0).
URL: https://data.mendeley.com/public-files/datasets/mc49zmg7n7/files/66392424-e44d-4d66-b949-44eaa5b1799a/file_downloaded
Not downloaded (not requested, large): PBO_stations.mat 787 MB, UNR_stations_new.mat 1.51 GB.
Contents: Lats/Lons/ID_G (443 GPS stations), Ev/Nv/Uv_SSE + sigmas, el (3607x3 triangles), slip_rate
(3607,), strike/rake/ds_fact/ss_fact, Pole.
**BLOCKER: the mesh NODE coordinates are NOT in the deposit** — they live in homog_GFs_Nov2019.mat /
Heter_GFs_May2019, which the README says were "not computed specifically for this study" and which are
neither on Mendeley nor in github.com/noelbartlow/Long-Term-Cascadia-ETS-velocities (4 .m files only). So
slip_rate cannot be georeferenced from public data; to read Bartlow's Oregon slip we must digitise her
figure, email the author, or rebuild the mesh.
**What the deposit DOES settle:** her station distribution. 443 stations, 330 west of -120 (her ind_GPS
cut). Per-degree west of -120: 39-40:37, 40-41:40, 41-42:21, 42-43:21, 43-44:21, 44-45:26, 45-46:31,
46-47:38, 47-48:35, 48-49:29, 49-50:15, 50-51:4. **Central Oregon 44-46N = 57, north 46-48N = 73, ratio
0.78** — a mild deficit, not a hole, and consistent with S12's 0.91x detection threshold. Supports C2.

**Schmidt & Gao (2010) JGR 115 B00A18 — PDF obtained (faculty.washington.edu/dasc/PAPERS/2008JB006045.pdf),
read in full. THE CAVEAT IS TIME-BOUNDED AND DOES NOT APPLY TO US.**
Abstract p.1: "Limited station coverage south of Portland (45.5N latitude) restricts our analysis to events
on the northern half of the subduction zone."
Para [15], PDF p.5: "Addition of GPS stations on Vancouver Island helps to improve resolution in 2004.
Fault plane resolution greatly improves in Washington beginning in 2005 with the densification of GPS
benchmarks in Puget Sound (47.5N latitude). **In Oregon, resolution is poor PRIOR TO 2007 from inadequate
station coverage.** Therefore, we focus our analysis on events that occur north of Portland (>45.5N)
because the GPS stations are too widely distributed in the south for proper inversions. Likewise, we also
limit our analysis to regions south of 49.5N on Vancouver Island."
=> Our window is 2010-2026, entirely post-PBO. **Correcting my S10 first pass: citing Schmidt & Gao as a
general Oregon-resolution objection to R11 is a misreading of the source.**
Also PDF p.7: "We attributed the patchiness of the slip distribution to the fault resolution because slip
tended to be a MAXIMUM in the low-resolution regions between stations." => their low-resolution failure
mode makes spurious slip MAXIMA, not holes; a persistent slip HOLE is the harder artefact. Favours R11.
No formal resolution test in the paper ("resolution test"/"checkerboard" do not occur).
PDF p.11 (useful for C3's geometry control): interface dip ~7 deg in N Washington vs "~15 deg in central
Oregon" and ~20 deg on Vancouver Island.

**Brudzinski & Allen (2007) Geology 35, 907-910 — central Oregon quotes.**
p.908 c.1: "Previous observations from Oregon are sparse; ... data obtained using ONE GPS STATION in
central Oregon showed slow slip episodes with a longer recurrence interval (Szeliga et al., 2004)."
p.908 c.1: "nonvolcanic tremors that correlate with slow slip episodes are apparent in several new
locations along the subduction zone, particularly along central Cascadia, despite more limited
observatories (Figs. 2D-2F)." And: "30 stations reporting slow slip episodes and 55 stations reporting
nonvolcanic tremor (Fig. 1), of the >300 stations that have been investigated."
p.908 c.2: "in the Siletzia zone ... station COR in Corvallis, Oregon ... (Fig. 2E), giving a robust
estimate of 19 +- 3 months ... The average interval across the Siletzia zone (19 +- 4 months) is longer
than those observed on Vancouver Island to the north (14 +- 2 months), and is nearly twice as long as that
from California to the south (10 +- 2 months)."
Fig. 2E "North Oregon 19+-3 mo" / 2F "South Oregon 19+-3 mo"; GPS codes listed: albh, nano, cpxf, **newp**,
**ddsn**, ybhb, ntka. **newp is almost certainly Newport, Oregon (~44.58N) — inside the 44.3-45.4N band of
R11.** CAVEAT: panel-to-station mapping could not be verified from the PDF text layer; definitive check is
Fig. DR4 in GSA Data Repository item 2007224 (not retrieved).
**Method caveat that limits its weight against R11** (p.908 c.1): their GPS detection is a single-station
hyperbolic-tangent fit with a 99% f-test, and they say "These methods are NOT replacements for network
solutions, but are simply used as surrogates". Not a slip inversion; cannot localise slip to 44.3-45.4N at
27-38 km. No Oregon slip-amplitude estimate anywhere in the paper.

**NET EFFECT ON C2: strengthened again** — (i) Schmidt & Gao's caveat is pre-2007; (ii) that class of
inversion fails by making slip maxima, not holes; (iii) Bartlow's own network gives central Oregon 0.78x
the northern density, matching S12's 0.91x; (iv) Brudzinski & Allen's Oregon detections are single-station
tanh fits the authors decline to call network solutions. Only decisive test left: Bartlow's slip map.

**REVISED RANKING:** 1) C1 exponent (GRL/JGR-SE headline; one blocker, non-vbICA slip model). 2) C2 central
Oregon (GRL once Bartlow's slip map is obtained). 3) C4 stationarity (GRL only if aimed at Meng et al.
2026). 4) C3 deepening (GRL if it survives; least controlled). 5) **C6 recurrence — C, a section not a
paper**; BSSA note or short GRL only if the pre-event sign reversal holds and survives a block bootstrap.
6) C5 = the introduction's gap statement.

## 2026-09-10 orchestrator review — S10b
Predictability (S8 candidate): LABEL INVERTED — size predicting the NEXT interval is TIME-predictable (Shimazaki & Nakata 1980);
and Gao, Schmidt & Weldon 2012 (BSSA) already reported northern Cascadia SSEs "weakly time predictable" with the same sign; the
only new element is our weakly POSITIVE size-vs-previous-interval (+0.16) where they found anti-correlation (R^2 0.48, n 11) —
but our effective n is far below 195 (the north is one tremor band) and the model-slip proxy shows nothing. Verdict: a section at
most, not a result. RESULTS R13 corrected. Bartlow 2020 deposit downloaded (results/external/, SHA verified) but its 3607-triangle
slip_rate has no coordinates in the deposit -> cannot be georeferenced without the author's mesh. Coverage: Bartlow's network has
57 stations in 44-46N vs 73 in 46-48N (ratio 0.78), matching our S12 threshold ratio 0.91. Schmidt & Gao 2010's Oregon-resolution
caveat is explicitly PRE-2007 (irrelevant to 2010-2026), and poor resolution produces spurious slip MAXIMA between stations, not
holes -> the central-Oregon amplitude deficit is not a coverage artefact by that mechanism either. Brudzinski & Allen 2007: one
central-Oregon GPS station (Newport ~44.6N) with single-station fits. Literature ranking now: C1 scaling > C2 central Oregon
(as an amplitude heterogeneity) > C4 stationarity (aimed at Meng 2026) > C3 deepening (pending S9b) > C6 predictability.

## 2026-09-10 S11b model-free scaling

Ran `scripts/s11b_scaling_gnss.py` (OMP_NUM_THREADS=4, seed 20260910, **40 s**). Follow-up to S11:
the S10 literature check left one blocker on the sub-linear exponent — the regressor was the
vbICA slip rate, i.e. a 13-fixed-pattern inversion. S11b re-measures the same exponent against a
slip-rate proxy built **straight from GNSS displacements, with no slip patterns, no Green's
functions and no inversion**. Events, window, footprint, tremor measures, nulls and fitting are
identical to S9/S11; only the regressor changes. The vbICA R7 is rebuilt inside the same script
as the reference column and reproduces S11 exactly (+0.401 for C, +0.644 for E).

### The GNSS proxies
`X.mat` (read via `sse_strain.load_solution`, the reader `scripts/s12_representability.py` uses):
639 channels = **213 stations x [e,n,u]**, 7169 daily epochs 2007-01-02..2026-08-18 — the SAME day
grid as the slip model, so nothing is re-gridded. Per channel: gap mask applied (91.5 % real
data), then one least-squares removal of constant + linear trend + 365.25-d and 182.6-d sinusoids
over the whole record, then a NaN-aware centred moving mean, then a first difference (mm/d).
Per event, the stations within 120 km of the footprint (**74 [65, 84]** of 213) give
  * PROJECTED rate = (station-mean v_e, station-mean v_n) . u, where u = the unit vector of the
    station-mean horizontal displacement CHANGE over [start, end]; u is **west-south-west in
    21/27 events** (global median azimuth **247.3 deg** = the trench-ward ETS reversal; the 6
    off-azimuth events are the lowest-SNR ones, incl. 1878, 2015, 1905);
  * SUM-|E| rate = sum over stations of |v_e|.
Two sources: **RAW** = the masked displacements (fully independent of the slip model) and
**ICA13** = the 13 retained ICA components reconstructed as U S V^T (denoised, still no slip
patterns, but NOT independent of the vbICA decomposition — reported as an intermediate).
Amplitudes: data rms 4.71 mm, residual after detrend/deseason 3.61 mm, **13-component
reconstruction only 1.07 mm** — the ETS signal is ~1/4 of the daily scatter, which is the whole
difficulty.

**Bandwidth matters and had to be matched.** The vbICA `A_raw` amplitudes are TRAILING-30-DAY
quantities, so the modeled R7 is effectively a ~30-d rate: its lag-1 autocorrelation of log rate
is **+0.664**, against **+0.231** for a 7-d GNSS rate. Comparing a 7-d GNSS rate with R7 is
therefore not like for like. A **centred 31-d** GNSS proxy (same bandwidth, no trailing-window
lag) has lag-1 **+0.410** and is the one to quote. Per-event SNR (mean rate over DUR / sd of the
same proxy over [start-70, start-10]): raw 7-d **0.43**, raw 15-d **0.67**, raw 31-d **1.00**,
ICA13 31-d 1.29, vbICA R7 1.75.

### Exponents (median [IQR] over events; 90 % bootstrap CI of the median; pooled = event dummies)

| measure | proxy | n | b median [IQR] | boot90 CI | pooled OLS/TS | b<1 | sign vs null |
|---|---|---|---|---|---|---|---|
| C | Graw 7-d      | 26 | +0.128 [+0.00,+0.26] | [+0.04,+0.23] | +0.154/+0.173 | 26/26 | 19/26 p=0.029 |
| C | Graw 7-d, fixed dir | 26 | +0.181 [+0.10,+0.29] | [+0.13,+0.26] | +0.209/+0.241 | 25/26 | 21/26 p=0.0025 |
| C | Graw 15-d     | 26 | +0.235 [-0.01,+0.38] | [+0.14,+0.31] | +0.296/+0.311 | 26/26 | 20/26 p=0.0094 |
| C | **Graw 31-d (matched)** | 26 | **+0.278** [+0.05,+0.65] | **[+0.095,+0.514]** | +0.358/+0.383 | 22/26 | 21/26 p=0.0025 |
| C | Gica 31-d     | 25 | +0.322 [+0.04,+0.64] | [+0.09,+0.59] | +0.369/+0.383 | 23/25 | 20/25 p=0.0041 |
| C | Graw sum-\|E\| | 27 | +0.806 [-0.03,+1.69] | [+0.32,+1.13] | +1.046/+1.091 | 16/27 | 20/27 p=0.019 |
| C | R7 vbICA (ref)| 27 | +0.401 [+0.18,+0.48] | [+0.31,+0.45] | +0.328/+0.397 | 27/27 | 26/27 p=4.2e-07 |
| E | Graw 7-d      | 17 | +0.247 [+0.14,+0.64] | [+0.15,+0.38] | +0.278/+0.337 | 16/17 | 13/17 p=0.049 |
| E | **Graw 31-d (matched)** | 17 | **+0.305** [+0.02,+1.13] | **[+0.043,+0.712]** | +0.558/+0.628 | 12/17 | 13/17 p=0.049 |
| E | Gica 31-d     | 16 | +0.592 [+0.07,+1.17] | [+0.11,+1.10] | +0.540/+0.645 | 10/16 | 13/16 p=0.021 |
| E | Graw sum-\|E\|| 18 | +1.001 [+0.12,+2.37] | [+0.44,+1.95] | +1.590/+1.608 | 9/18 | 14/18 p=0.031 |
| E | R7 vbICA (ref)| 18 | +0.644 [+0.17,+0.84] | [+0.33,+0.76] | +0.500/+0.615 | 14/18 | 15/18 p=0.0075 |

(Event 1905 — 8 days, already flagged in S1/S4/S9/S11 — is the one event with no usable GNSS
proxy: its 31-d rate is negative through the window, so n = 26 rather than 27.)

### The three answers asked for
**(1) Is the exponent still < 1 with a raw-GNSS regressor?  YES.** For the bandwidth-matched raw
proxy, b(C) = **+0.278**, 90 % CI of the median **[+0.095, +0.514]**, **22/26 events below 1**,
21/26 positive against the shifted-tremor null (p = 0.0025); b(E) = +0.305, CI [+0.04, +0.71],
12/17 below 1. Every raw-GNSS variant (7-d, 15-d, 31-d, fixed direction) has an upper CI bound
below 1 for the count. The only variant that cannot exclude 1 is the crude **sum-|v_e|** proxy,
whose CI is [+0.32, +1.13] — that measure has no directional discrimination and its apparent
SNR ~5 is just the positive mean of a sum of absolute noise, so it is uninformative rather than
contradictory.

**(2) Value with an honest range.** Model-free count exponent **b ~ 0.3, range 0.1-0.5**
(bootstrap CI of the median [+0.095, +0.514]; pooled event-dummy fit +0.38); model-free energy
exponent **b ~ 0.3-0.6, range 0.04-0.7** (median +0.305, pooled +0.63, ICA13 version +0.59).
Combining with S11, the defensible statement is **b = 0.3-0.7, centred near 0.4 for counts and
0.6 for energy, with b = 1 excluded**.

**(3) Does it agree with the vbICA 0.40 / 0.64?  YES, once the bandwidth is matched.** Paired
per event:
- b(C, Graw31) - b(C, R7) = **+0.023 [-0.334, +0.321], Wilcoxon p = 0.92**; Spearman rho = +0.52,
  p = 0.0063.
- b(E, Graw31) - b(E, R7) = **-0.156 [-0.439, +0.748], p = 0.78**; rho = +0.46, p = 0.066.
- ICA13 31-d: differences +0.029 (p = 0.98) for C and -0.095 (p = 0.98) for E.
With the MISMATCHED 7-d proxies the difference is significantly negative (C: -0.201,
p = 0.0036; E: -0.187, p = 0.064) — i.e. the apparent disagreement is entirely a
bandwidth/SNR effect, not a disagreement about the physics.

### Why the raw exponent starts low: the same attenuation ladder S11 found
b(C) rises monotonically with the proxy's SNR — 0.128 (SNR 0.43) -> 0.235 (0.67) -> 0.278 (1.00)
-> 0.401 (1.75, the model) — exactly the errors-in-variables behaviour S11 demonstrated by
simulation for R1/R3/R7. Independent check: the reliability ratio lambda = Theil-Sen slope of
log R7 on log Graw(7-d) is **+0.436 [+0.273, +0.638]**, and the de-attenuated exponents are
b(C) = **+0.324 [+0.06, +0.75]** and b(E) = **+0.648 [+0.18, +1.16]** — i.e. the raw 7-d numbers
are consistent with the model-based 0.40 / 0.64 once the regressor noise is accounted for. (The
same correction applied to the 31-d proxy gives +0.66 / +1.05, but there lambda = 0.26 mixes
residual bandwidth difference with noise and over-corrects; treat those as the upper bracket.)

### Robustness inside S11b
- **Direction estimate is not doing the work.** Replacing the per-event ETS direction with one
  global fixed direction (azimuth 247.3 deg) gives b(C) = +0.181 [CI +0.13, +0.26], the same
  answer as the per-event version (+0.128) within its CI.
- **Saturation, again nothing.** Upper minus lower half of log10 proxy: Graw31 -0.127
  (p = 0.47), Gica31 +0.033 (p = 0.92), Graw7 +0.305 (p = 0.37), R7 +0.099 (p = 0.35). No
  flattening at high rate on any regressor.
- **The requested "rate > 2 x pre-event noise" restriction is not usable at this SNR** and is
  reported as such: it leaves >= 8 usable days in only 16/27 events for Graw31 and 4/27 for
  Graw7, giving b(C, Graw31, hisnr) = +0.220 [-0.34, +0.70] (10/14 positive, p = 0.18) — and
  truncating on a noisy regressor is itself a biased selection. Its one clear message is that a
  DAILY raw-GNSS slip rate is at best a ~1-sigma measurement, which is the point of the
  attenuation analysis above.

### Caveats
- The events, windows and footprints still come from the slip-patch catalogue; S11b removes the
  slip model from the **regressor**, not from the event DEFINITION. A fully model-free event
  definition (tremor-episode based) is a separate piece of work.
- `Gica13` shares the 13 time functions with the vbICA inversion, so it is not independent
  evidence; only the RAW rows are.
- The GNSS proxy is a surface-velocity scalar, not a potency rate: it mixes the slip amplitude
  with the source-station geometry, so its ABSOLUTE level is not comparable across events. The
  per-event demeaning / per-event Theil-Sen removes exactly that, which is why only the exponent
  is reported.
- corr(log Graw31, log R7) is only **+0.25 [+0.02, +0.37]** day by day: the two regressors agree
  on the EXPONENT but track each other loosely at daily resolution, consistent with both being
  noisy views of the same slow transient.
- n = 26 (count) and n = 17 (energy) usable events; the energy column starts 2017-01-01.

Files written: `scripts/s11b_scaling_gnss.py`, `results/s11b_exponents.csv` (72 rows = 2 measures
x 9 proxies x 4 subsets), `results/s11b_events.csv` (27 rows), `results/s11b_proxy_check.csv`
(27 rows), `figures/s11b_scaling_gnss.png` (A example event with all four rate proxies,
B exponent by proxy for count and energy, C per-event agreement with the model-based exponent,
D pooled event-dummy fit on the model-free regressor, E exponent vs proxy SNR, F observed
exponents vs the shifted-tremor null). Nothing outside `tremor_slip_geometry/` was created or
modified.

## 2026-09-10 orchestrator review — S11b
The sub-linear scaling survives removal of the slip model from the regressor: with a raw-GNSS slip-rate proxy (213 stations from
X.mat, detrended/deseasoned, projected on the per-event ETS direction, azimuth 247 deg = trench-ward in 21/27 events; 31-d
bandwidth-matched), b(count) = 0.28 [0.10, 0.51] (22/26 events below 1; vs shifted-tremor null p 0.0025), b(energy) = 0.31
[0.04, 0.71]; paired difference to the vbICA-based exponents +0.02 (p 0.92) at matched bandwidth. Exponent rises monotonically
with proxy SNR (0.13 -> 0.40 from SNR 0.4 -> 1.75), i.e. the raw values are attenuated; de-attenuation gives 0.32 (counts) /
0.65 (energy) = the vbICA numbers. b = 1 excluded in every count variant. Remaining dependence on the slip model: only the event
windows/footprints (selection), not the regressor. THE SCALING RESULT IS NOW CLAIMABLE: tremor count ~ slip-rate^0.3-0.4,
tremor energy ~ slip-rate^0.6-0.65, constant across 27 events, two segments and 16 years, in both a vbICA model and raw GNSS.

## 2026-09-10 S14 launched (LFE-based scaling exponent) — see PLAN.md round 4

## 2026-09-10 S9b depth-phase reconciliation

Coordinator follow-up on S9 Part A. Literature check: Wech & Creager (2011, Nat. Geosci. 4, 624) plotted
median daily down-dip distance through four ETS and reported initiation far down-dip with UP-DIP migration
during the episode -- apparently opposite to S9's +2.26 km first-to-last-fifth DEEPENING. Ran
`scripts/s9b_depth_phase_reconcile.py` (OMP_NUM_THREADS=4, seed 20260910, **146 s**, loaders identical to
S9/S4). Four parts, all on the same 27 confirmed events.

### Two measurement diagnostics found first (both matter)

1. **Down-dip distance and depth are the same axis only WITHIN a latitude band.** Built the mesh's up-dip
   edge (shallowest triangle centroid per 0.1 deg latitude, 82 points, lat 39.83-49.92, edge depth
   11.2-50.6 km) and set ddist = horizontal km from it. Spearman(depth, ddist) **pooled over the margin is
   +0.080**, but **within 0.5 deg latitude bands the median is +0.970** [+0.284, +0.999] over 20 bands. The
   slab dip changes along strike (the up-dip edge runs lon -124.4 at 40N to -126.8 at 49N while tremor depth
   goes the other way, ~48 km south vs ~39 km north), so a detection can be far down-dip AND shallow.
2. **Therefore a confound S9 did not control**: a neighbourhood spans 220-720 km along strike (S4), so
   tremor merely MIGRATING ALONG STRIKE into a deeper part of the margin raises the daily median depth with
   no across-strike motion. Fixed by latitude-detrending every detection against the margin-wide median at
   its own 0.25 deg latitude (41/44 bins usable; **3 all-NaN bins outside the mesh hull had to be dropped
   before interpolating -- left in, `np.interp` silently NaNs a whole latitude range**). Residual depth
   sd 4.60 km vs raw 6.93 km. Everything below is reported raw AND latitude-detrended.

### (1) Finer phase axis and the onset / plateau / tail split

10-bin stacked anomaly (median over 27 events, + = deeper/down-dip; null = the same 200 rigid >= 365 d
tremor-date shifts as S9, 200/200 draws kept):

    bin              1      2      3      4      5      6      7      8      9     10
    depth        +0.57  +0.87  -0.19  -0.17  -0.07  -1.34  +0.49  +1.71  +0.98  +3.91   (null band ~ +-1.5)
    depth_resid  +0.86  +0.00  -0.74  -0.09  +0.10  -0.05  +0.03  +1.75  +0.47  +3.24
    ddist_resid  +5.09  +0.00  -4.46  -1.30  +0.18  -0.87  -0.11  +6.72  +4.43 +14.69

    quantity      bin1 -> minimum          minimum -> bin10        net bin1 -> bin10
    depth        -1.91 km  (p = 0.28)      +5.25 km  (p < 0.005)   +3.34 km  (p < 0.005)
    depth_resid  -1.60 km  (p = 0.22)      +3.99 km  (p < 0.005)   +2.39 km  (p = 0.005)
    ddist_resid  -9.55 km                  +19.15 km               +9.60 km   (no null; see caveat)

- **The profile IS non-monotonic in exactly the direction that would reconcile the two results**: it
  shallows (UP-DIP) over the first half and deepens (DOWN-DIP) over the last two tenths, with the deepest
  bin last (+3.9 km raw, +3.2 km across-strike, both clear of the null band).
- **But only the LATE leg beats the null.** The early up-dip leg (-1.6 to -1.9 km) is *inside* the
  shifted-window null (null mean -1.4 +- 0.8 raw, -1.1 +- 0.7 detrended; p = 0.22-0.28), because taking the
  minimum of a noisy 10-bin profile is itself a selection statistic that the null reproduces. **So our data
  are CONSISTENT with an up-dip early phase but cannot independently confirm it**; only the tail deepening
  is a measurement.
- **The deepening is ACROSS-STRIKE, not along-strike**: it survives latitude detrending at +2.39 km,
  p = 0.005. Tremor really does move down-dip relative to the local tremor band, it is not just sliding
  along the margin into deeper geometry. In Wech & Creager's own coordinate the last tenth sits
  **+14.7 km [+11.3, +20.5] further down-dip** than the event's own median.

ONSET / PLATEAU / TAIL (peak = smoothed daily-count maximum; plateau = the run of days with >= 50 % of the
peak; onset/tail = the defined days before/after it):

    onset -> peak   stacked median -0.911 km   6/20 positive   null median -1.102 sd 0.940 -> **p = 0.58**
    peak  -> end    stacked median +0.345 km  14/27 positive   null median +0.380 sd 0.930 -> **p = 0.74**
    per-event z vs its own null: -0.23 [-0.41, +0.74] and -0.13 [-0.44, +0.39]

- **This decomposition finds NOTHING, and that is the single most informative result in S9b.** The null
  reproduces both numbers almost exactly (-1.10 vs -0.91; +0.38 vs +0.35). The split is conditioned on the
  COUNT series, so the null inherits the generic depth-intensity coupling; once you compare high-count days
  with high-count days, the onset-versus-tail depth contrast disappears entirely. The phase-binned test
  survives only because a random window has no ETS ramp aligned to its phase axis.
- Mechanistic consequence, stated plainly: **the defensible statement is "the sparse waning tail of an
  episode is deeper", not "the tremor source migrates down-dip"**. S9 already showed the detection rate
  collapses ~7x over the same interval; S9b shows that conditioning on that rate removes the effect.

### (2) Wech & Creager's own plot on our events (no phase normalisation)

Latitude-detrended daily medians vs days from event start, stacked over the 27 events:

    series        anchor              days 0-9 -> 10-24        days 10-24 -> 25-60
    ddist_resid   catalog slip start   -0.64 km (n=27)          +10.56 km (n=23)
    ddist_resid   first tremor day     -3.81 km (n=27)           +2.54 km (n=24)
    depth_resid   catalog slip start   -0.00 km (n=27)           +3.27 km (n=23)
    depth_resid   first tremor day     -0.64 km (n=27)           +0.69 km (n=24)

- Same sign pattern (early up-dip, late down-dip) but **much weaker and much noisier than the
  phase-normalised version**, and the two anchorings disagree by a factor 4. The reason is structural: our
  durations span **2 to 85 days**, so a fixed day-lag axis mixes the early phase of a long event with the
  tail of a short one; by lag +45 only 2-10 events remain. Wech & Creager's four similar-length ETS do not
  have this problem -- **their un-normalised axis is appropriate for their sample and not for ours**, which
  is why S9 phase-normalised.
- **Overlap with their episodes**: the PNSN file starts 2010-01-06, so at most our two pre-2012 northern
  events can overlap their 2010-2011 Washington ETS -- 1601 (2010-08-14: onset->peak -2.92 km,
  peak->end +0.78 km, net fifth +1.84 km) and 1585 (2011-07-24: -0.53, +1.94, +3.14 km). Both show the
  up-dip-then-down-dip shape. **Their exact four-episode list could not be verified from this environment
  (no literature access)**, and any episode of theirs before 2010 is outside our catalogue entirely.

### (3) Does the sign depend on the start definition or the radius?

    anchor    R(km)  n_ev   depth      depth_detr   null mean+-sd     p        raw ddist
    tremor      40    27    +2.430      +2.432     -0.377+-0.826    0.040       -0.84
    catalog     40    24    +2.217      +1.439     -0.456+-1.015    0.023       +0.88
    tremor      60    27    +2.256      +2.160     -0.117+-0.936    0.010       -2.80
    catalog     60    24    +1.766      +1.312     -0.250+-0.975    0.087       +0.93
    tremor     100    27    +1.534      +2.223     -0.046+-1.094    0.160       -4.67
    catalog    100    24    +2.142      +2.317     -0.224+-1.105    0.084       +5.12

- **Positive (deepening) in 6/6 cells**, raw +1.53 to +2.43 km, detrended +1.31 to +2.43 km; p between
  0.010 and 0.160. The sign does not depend on anchoring the phase to the first tremor day or to the
  catalog slip start, nor on the aperture.
- **CROSS-CHECK: the (tremor, 60 km) cell reproduces S9 exactly, +2.256 km** (S9 p = 0.020 with 200 draws,
  0.010 here with 100) -- the sensitivity table uses S9's statistic (difference of stacked bin medians),
  so the two blocks are directly comparable.
- The RAW down-dip-distance column flips sign with radius (-0.84 / -2.80 / -4.67 tremor-anchored vs
  +0.88 / +0.93 / +5.12 catalog-anchored): the un-detrended distance is dominated by along-strike geometry,
  a third demonstration that only the latitude-detrended version is interpretable.

### (4) Is "shallower when most intense" temporal, or the along-strike depth-amplitude gradient?

    raw     Spearman(daily depth, daily count)        -0.107 [-0.321, +0.092]  19/27 neg  Wilcoxon p = 0.049
    partial Spearman(depth, count | along-strike y)   -0.123 [-0.323, +0.017]  19/27 neg  Wilcoxon p = 0.017
    partial - raw                                     -0.001 [-0.075, +0.113]             Wilcoxon p = 0.86
    controlled-for gradient: Spearman(depth, y) +0.027 [-0.594, +0.666];  Spearman(count, y) +0.219 [+0.04, +0.35]

- **The intensity-depth relation SURVIVES the spatial control completely unchanged** (median -0.107 ->
  -0.123; the per-event shift is -0.001, p = 0.86; the same 19/27 events negative, and the Wilcoxon
  actually strengthens). So S9's "tremor is shallower when most intense" is a genuine TEMPORAL relation
  within events, not the along-strike depth-amplitude gradient.
- The reason the control changes nothing is visible in the gradient itself: **within a single event there is
  no consistent along-strike depth gradient** -- Spearman(depth, y) has median +0.027 with IQR
  [-0.59, +0.67], i.e. it is strongly event-specific and averages out. Count does rise weakly northward
  (+0.219 [+0.04, +0.35]).

### Answers to the two questions asked

**Is there up-dip migration during onset and down-dip migration in the tail (reconciling both), or a
genuine disagreement?** **A partial reconciliation, with one leg measured and one leg merely consistent.**
The 10-bin profile has exactly the reconciling shape -- shallowing to a minimum around mid-event
(-1.6 km across-strike) then deepening steeply into the last tenth (+4.0 km across-strike) -- and the
deepening is a genuine across-strike migration (survives latitude detrending, p = 0.005; +14.7 km in
Wech & Creager's down-dip-distance coordinate) that is robust in 6/6 anchoring x radius cells. But the
early up-dip leg is inside the shifted-window null (p = 0.22-0.28), so **we cannot claim to have
reproduced their up-dip migration; we can only say our data do not contradict it**. There is therefore no
demonstrated disagreement of sign: the two statements concern different parts of the episode and only the
late leg is resolved here. The strongest caution is the onset/plateau/tail test, which conditions on
intensity and finds nothing at all (p = 0.58 / 0.74) -- so the honest mechanism is **"the waning, sparse
tail of an episode is systematically deeper"**, and neither our data nor this design can separate that
from a true migration of the active source.

**Does the intensity-depth relation survive the spatial control?** **Yes, unchanged.** Removing the daily
along-strike position from both variables moves the stacked partial rho from -0.107 to -0.123 (per-event
shift -0.001, p = 0.86) and leaves 19/27 events negative. The Wech (2021) spatial depth-amplitude gradient
is not the explanation, because within a single event's neighbourhood there is no consistent along-strike
depth gradient to begin with (Spearman(depth, y) median +0.03, IQR -0.59 to +0.67).

### Caveats specific to S9b
- **The up-dip leg is not measured, only "not excluded".** Do not quote -1.6/-1.9 km as a result.
- **The onset/plateau/tail null kills the like-for-like version of the effect.** Any mechanistic reading of
  the tail deepening must state that it is inseparable from the ~7x intensity collapse (S9's count control).
- **The `ddist` series has no shifted-window null** (only the two depth series do); its numbers are point
  estimates with bootstrap CIs, quoted to give the result in Wech & Creager's units, not as tests.
- The raw (un-detrended) down-dip distance is confounded by along-strike slab geometry and should not be
  used at all; it is reported only to demonstrate the confound.
- **Wech & Creager's episode list is unverified here** (no literature access); the two overlapping events
  named are simply every northern confirmed event of ours before 2012.
- Same small-n and short-event limits as S9: 7 events have an onset segment shorter than 3 defined days and
  get no onset->peak value (n = 20 of 27), and 1434/1905/1655 remain thin throughout.

Files written: `scripts/s9b_depth_phase_reconcile.py`, `results/s9b_phase10.csv` (40 rows = 4 series x 10
bins), `results/s9b_segments.csv` (27 event rows + 2 stack rows), `results/s9b_downdip_stack.csv` (86 lags),
`results/s9b_sensitivity.csv` (6 rows), `results/s9b_partial_rho.csv` (27 rows),
`figures/s9b_depth_phase.png` (6 panels: 10-bin raw and latitude-detrended profiles with nulls;
onset->peak / peak->end boxes; the day-aligned Wech & Creager view; sign stability over anchoring x radius;
raw vs partial rho). Nothing outside `tremor_slip_geometry/` was created or modified.

## 2026-09-10 orchestrator review — S9b
Reconciliation with Wech & Creager 2011 is partial: the 10-bin profile shallows to a mid-event minimum then deepens steeply in the
last tenth (late leg +4.0 km across-strike, p < 0.005; early up-dip leg -1.6 km, p 0.22 = not confirmed, not contradicted).
Latitude-detrending (slab dip changes along strike, so along-strike migration could masquerade as depth change) KEEPS the result:
+2.4 km, p 0.005 — genuinely across-strike. But the count-conditioned onset/plateau/tail split finds nothing, so the defensible
statement is "the sparse waning tail of an ETS is deeper", not "the source migrates down-dip". The intensity-depth relation
survives the spatial control (partial rho -0.12, unchanged). Verdict: a supporting observation with caveats (detection geometry
in the sparse tail cannot be excluded), not a headline. RESULTS R8 reworded.

## 2026-09-10 S13 nucleation, direction, drift power

Ran `scripts/s13_nucleation_direction_power.py` end to end (OMP_NUM_THREADS=4, seed 20260910, 791 s;
one earlier aborted run + one `S13_FAST=1` smoke run, neither wrote published output). Motivated by
S10 candidate C4: the only VERIFIED published secular-evolution claim for the Cascadia ETS zone is
Meng, Liu, Adhikari & Luo 2026 GRL e2025GL120690 — "a long-term change in the locations where ETS
events nucleate and their propagation directions" — and S1/S4/S5/S7/S9 never measured either
quantity. A null aimed at a named claim also needs a power analysis. S13 supplies both.
**PART A is MODEL-FREE**: PNSN catalogue only, no vbICA slip, no mesh, so it cannot inherit any
basis artefact.

### PART A.1 — margin-wide tremor events straight from the catalogue
Detector (deterministic, no eps tuning; ST-DBSCAN was rejected because its eps is a DENSITY
calibration and the PNSN density triples across the record — the very non-stationarity under test):
latitude bins 0.25 deg over 39.5–51 N; `C[b,d]` = detections on day d within a 200-km moving
latitude window centred on bin b; `Active = C >= 50`; 1-day holes closed; EVENT = an **8-connected
component** of Active; a detection joins the event owning its own (bin, day) cell; keep events with
**>= 500 detections and >= 5 days**.

    191 tremor events, 547,288 assigned detections = 73.2 % of the 748,002-detection catalogue
    duration  d : min 5   q25 7    med 10   q75 18   max 51
    n_det       : min 504 q25 841  med 1667 q75 3590 max 16325
    extent   km : min 65  q25 216  med 315  q75 401  max 985
    by NUCLEATION segment: south 118 / north 63 / central 10
    events per year 2010..2026: 8 11 10 13 12 15 7 10 10 16 9 14 12 9 14 16 5(partial)
    migration slope defined for 180/191; median |v| 6.0 km/d
    "bilateral/stationary" (|v| < 0.1 deg/d = 11 km/d, the spec's bar): 75.9 %
    northward fraction: 0.492 (all events) | 0.370 (|v| >= 11 km/d, n=46) | 0.535 (|v| >= 2 km/d, n=159)

Catalogue variants for the completeness controls (threshold and size floor rescaled by the variant's
own rate ratio r): `full` r=1.00 TH=50 → 191 ev; `bright` (days with >= 100 detections margin-wide)
r=0.87 TH=43 → 173 ev; `full2017` (2017+) TH=50 → 115 ev; `sta8_2017` (num_stas >= 8, 2017+;
num_stas is populated only from 2017-01-01) r=0.485 TH=24 → 100 ev.

### PART A.2 — do the 27 confirmed slip events sit inside tremor events?
Matched against tremor events overlapping `[slip_start - 30 d, slip_end]` (slip times are 30-d
TRAILING crossings, hence the 30-d back-shift); a "time+space" match also requires the tremor event's
latitude range (+-0.5 deg) to contain `slip_lat`.

    >= 1 time-overlapping tremor event   : 26/27      exactly one : 1/27
    >= 1 time+space matching tremor event: 25/27      exactly one : 16/27
    n time-overlaps  {0:1, 1:1, 2:7, 3:9, 4:8, 5:1}
    n time+space     {0:2, 1:16, 2:7, 3:2}     23 distinct tremor events used as best matches

So YES in the sense that matters (25/27 slip events have a tremor event at the right place and time),
but "contained in ONE" is only true for 16/27: the slip windows are 3–84 d long and chain along
strike (S1/S4 caveat), so they routinely span 2–4 of our tremor events. The 2 unmatched are 1905
(2022-12-02, 3 d) and 1434 (2023-11-11, 2 d, already flagged unreliable in S1/S4).
Cross-check: 198/359 (55.2 %) of the S5 per-band episodes fall inside an S13 event — S5 counts every
per-band peak including ones below our 500-detection floor, so the two catalogues are different
objects by construction.

### PART A.3(i) — nucleation-latitude trends (Theil-Sen, 1000 date permutations)
`p_yb` = year-block permutation p (whole calendar years exchangeable, so within-year clustering
survives the null). `p_H` = Holm over the whole 16-test family (4 variants x 4 segments).

    variant    group      n  km/yr   95% CI (km/yr)   p_perm   p_yb    p_H    med2010-17 -> 2018-26  MW p
    full       combined 191  +4.45  [ -2.40, +12.11]  0.210   0.266   1.00    42.575 -> 42.626      0.683
    full       north     63  +3.00  [ -1.99,  +9.06]  0.287   0.279   1.00    48.380 -> 48.355      0.550
    full       central   10  -6.02  [-21.22,  +8.64]  0.359   0.282   1.00    45.183 -> 45.311      0.833
    full       south    118  +4.07  [ -0.23,  +9.11]  0.064   0.103   0.70    41.210 -> 41.580      0.131
    bright     combined 173  +3.60                    0.320   0.272   1.00
    bright     north     57  +2.29                    0.408   0.375   1.00
    bright     south    104  +2.91                    0.212   0.286   1.00
    full2017   combined 115 +21.33                    0.008   0.017   0.112
    full2017   south     71 +18.33                    0.001   0.003   0.016
    full2017   north     36 +10.70                    0.140   0.027   1.00
    sta8_2017  combined 100 +24.19                    0.021   0.003   0.273
    sta8_2017  south     57 +18.28                    0.002   0.033   0.030
    sta8_2017  north     36 +14.49                    0.026   0.013   0.312
The `nucleation_lat_p10` definition (median of the first 10 % of detections) gives the same answer
(full/combined +3.33 km/yr p 0.310).

**Over the full 2010-2026 record: nothing.** Every full-catalogue and bright-catalogue trend has
p >= 0.06 and a CI straddling zero, and no median split reaches p < 0.13.

### PART A.3(ii) — northward fraction per epoch (exact binomial CIs, 3x2 chi2, Fisher first-vs-last)

    full / all events   combined  0.490 (51) | 0.491 (55) | 0.568 (74)   chi2 p 0.596   Fisher p 0.466
                        north     0.500 (20) | 0.389 (18) | 0.417 (24)   chi2 p 0.766   Fisher p 0.762
                        south     0.467 (30) | 0.576 (33) | 0.652 (46)   chi2 p 0.278   Fisher p 0.154
                        central   n = 1/4/4, not testable
    full / |v|>=2 km/d  combined  0.479 (48) | 0.574 (47) | 0.547 (64)   chi2 p 0.628
    full / |v|>=11 km/d combined  0.286 ( 7) | 0.727 (11) | 0.412 (17)   chi2 p 0.131
    bright / all        north     0.750 (16) | 0.231 (13) | 0.333 (24)   chi2 p 0.0079  Fisher p 0.0225

**No change in propagation direction anywhere in the full catalogue.** The single nominal hit is
`bright/north` (p_chi2 0.0079): it is ABSENT in the full catalogue for the same segment (p 0.766),
rests on 16 events in epoch 1, and with 8 variant x segment direction tests its Holm p is ~0.06.
Note also that the spec's "bilateral" bar (11 km/d) sits ABOVE the observed median |v| of 6 km/d, so
it classifies most genuine ETS migrations as bilateral; that is why the 2 km/d subset is reported too
(it keeps 159/180 events and says the same thing).

### PART A.3(iii) — nucleation latitude vs event size, and the size drift itself
Spearman rho of nucleation latitude against `n_det` / `peak_daily`: combined +0.130 (p 0.072) /
+0.100 (p 0.171); north -0.221 (p 0.081) / +0.014; south +0.089 / +0.144; central n=10. Nothing.
But event SIZE itself drifts: **log10 peak_daily +0.0086/yr [+0.0033, +0.0138] p_perm 0.002**
(x1.4 over 16 yr) and log10 n_det +0.0101/yr p 0.100. That is the network, and it is the reason (iii)
has to be reported alongside it.

### PART A.4 — control (b), the 2018-02 PNSN algorithm change
    daily detection rate, 365-d means: after/before at 2018-02-01 = 0.834  (112.3 vs 134.6 det/d)
    null at 200 random dates: mean 1.055 sd 0.298 [p5 0.656, p95 1.572]
    -> observed sits at the 20.5th percentile; two-sided p on |log ratio| = 0.488   NO STEP
    tremor EVENTS per 365 d: after/before 1.000 (null 1.089, p 1.000)               NO STEP
    before vs after 2018-02-01 (full catalogue):
      combined  nucleation med 42.590 (87) -> 42.620 (104)  +3 km   MW p 0.747 | frac N 0.512 -> 0.531 Fisher 0.881
      north               48.335 (32) -> 48.428 (31)       +10 km   MW p 0.306 | 0.500 -> 0.367 Fisher 0.317
      south               41.210 (52) -> 41.580 (66)       +41 km   MW p 0.131 | 0.500 -> 0.639 Fisher 0.173
Wech 2021 reports +55 % epicentres with the new method; in the Cascadia-full catalogue as we hold it
that does NOT appear as a 365-d rate step, and nothing about nucleation or direction changes across
that date.

### PART A.4b — the short-window control (added after the first run found the 2017+ result)
The 2017+ window is 3433 d = **9.40 yr**. Sliding that window in 30-d steps and rerunning the same
test inside each (>= 20 events required, 300 permutations):

    combined  86 windows | slope range [-0.086, +0.184] deg/yr |  5 with p<0.05 ( 6 %), 5 N / 0 S
              window nearest 2017-01 (2016-12-25): +0.184 deg/yr p 0.027 = 100th percentile
    north     86 windows | [-0.084, +0.151]                    | 11 with p<0.05 (13 %), 11 N / 0 S
              nearest 2017-01: +0.119 p 0.033 = 90th percentile
    south     86 windows | [-0.086, +0.165]                    | 24 with p<0.05 (28 %), 24 N / 0 S
              nearest 2017-01: +0.165 p 0.003 = 100th percentile
The window trend rises MONOTONICALLY from about -7 to -10 km/yr (windows starting 2010) to +18 to
+20 km/yr (window starting 2017) — the signature of one low-frequency swing with a minimum around
2014-2016, not of a secular drift. Because the date-permutation null does not preserve that
structure, its p is anti-conservative inside a short window: 28 % of southern 9.4-yr windows reach
p<0.05, all of them northward.

Southern sub-zones (the south spans 39.5-44 N and contains two distinct belts; split at 42 N):

    epoch fraction of southern events nucleating in 42-44 N: 0.303 (2010-14) | 0.171 (15-19) | 0.460 (20-26)
      -> NON-MONOTONIC; median nucleation lat 41.29 | 40.95 | 41.65
    WITHIN 39.5-42 N: full record +0.88 km/yr p 0.487 | 2017+ +4.82 km/yr p 0.172
    WITHIN 42-44 N  : full record -3.38 km/yr p 0.060 | 2017+ -1.77 km/yr p 0.651
    zone MIX pre-2017 0.277 (n=47) vs 2017+ 0.366 (n=71), Fisher p 0.327
Neither belt moves. The apparent southern "drift" is a change in the relative productivity of two
fixed belts.

### PART A.5 — robustness to the detector parameters (18 sets: TH 30/50/80 x WIN 150/200/300 km x DLAT 0.25/0.5)
    n_events 157-222, median duration 10-13 d
    nucleation trend margin-wide: slope [+0.011, +0.040] deg/yr, min p_perm 0.169, 0/18 with p<0.05
    nucleation trend north      : slope [-0.005, +0.047] deg/yr, min p_perm 0.086, 0/18 with p<0.05
    epoch chi2 on northward frac: min p 0.202, 0/18 with p<0.05
The full-record answer does not depend on any detector choice.

### PART B — power of the stationarity bounds
Injected drift `slope*(date - 2018)` added to bootstrap-resampled Theil-Sen residuals at the observed
event dates (observed scatter and date sampling preserved), then the identical Theil-Sen permutation
trend test. 1000 realisations x 1000 permutations for the n=27 quantities; 300x300 (n>=100) or
500x500 for the nucleation-latitude groups. The slope-0 row is the empirical false-positive rate and
comes out 0.027-0.059 as it must.

    quantity              n   resid sd   MDD80 (80 % power)        over 16 yr      vs 28-km cell   vs 5-10 km PNSN loc err
    Q1 down-dip offset   27   8.95 km    0.970 km/yr               15.5 km         0.55 cells      1.6-3.1 x
    Q3 depth offset      27   1.79 km    0.198 km/yr                3.2 km         0.11 cells      0.32-0.63 x
    L tremor-slip lag    27   9.54 d     0.882 d/yr                14.1 d          --              --
    nucleation lat N     63   0.95 deg   0.068 deg/yr = 7.6 km/yr  121 km          4.3 cells       12-24 x
    nucleation lat S    118   1.01 deg   0.059 deg/yr = 6.5 km/yr  104 km          3.7 cells       10-21 x
    nucleation lat all  191   3.31 deg   0.102 deg/yr = 11.3 km/yr 181 km          6.5 cells       18-36 x
    nucleation lat C     10   0.60 deg   0.169 deg/yr = 18.8 km/yr 300 km         10.7 cells
    northward fraction   51 vs 74 events: smallest detectable epoch1->epoch3 change = 0.249 (absolute)

The published bounds restated as a total 2010-2026 drift (95 % CI half-width x 16 yr):
Q1 13.2 km (0.47 cells, 1.3-2.6 x the location error), **Q3 2.55 km (0.09 cells, 0.26-0.51 x)**,
Q4a 10.5 km, Q4b 8.0 km, L 9.0 d; nucleation latitude combined 194 km, north 145 km, south 146 km.

### PLAIN ANSWERS
**(i) Did nucleation locations or propagation directions change over 2010-2026?  NO — not detectably,
and the completeness controls neither create nor hide a change.** Margin-wide nucleation latitude
+4.5 km/yr [-2.4, +12.1], p 0.21; north +3.0 km/yr p 0.29; south +4.1 km/yr p 0.064 (the closest to
significant); 2010-17 vs 2018-26 medians differ by at most 0.37 deg (south, MW p 0.13). Northward
fraction 0.49/0.49/0.57 across the three epochs, chi2 p 0.60, and no segment or subset reaches
significance. This survives 18 detector parameter settings (0/18 significant in either family), the
bright-only catalogue, the >= 8-station 2017+ catalogue, and Holm correction. The 2018-02 PNSN
algorithm change is not even visible as a rate step here (p 0.49) and nothing changes across it.
The one apparent positive — a northward nucleation drift of +18 km/yr in the SOUTH confined to
2017-2026 (p 0.001-0.003, Holm 0.016/0.030) — is **not claimable**, for four independent reasons:
(a) it is absent from the full 16-yr record (+4.1 km/yr, p 0.064); (b) 28 % of southern 9.4-yr
windows of this record reach p<0.05 by the same test and ALL of them are northward, with the window
trend sweeping smoothly from -10 to +18 km/yr as the window slides — one low-frequency swing, and
the date-permutation p is anti-conservative against it; (c) it is a sub-zone MIX effect (neither the
39.5-42 N nor the 42-44 N belt moves; their epoch mix goes 0.30 -> 0.17 -> 0.46, non-monotonic);
(d) the two Holm-surviving rows share 57 of 71 events, so they are one test, not two. Honest wording:
*no secular change; a suggestive northward shift of southern-Cascadia nucleation confined to the last
~9 years that our record is too short to separate from a low-frequency oscillation.*

**(ii) What is the smallest drift our bounds could have detected?**
- **Q3 (tremor-slip DEPTH offset) is the strong bound: 0.20 km/yr, i.e. 3.2 km over 16 yr** — one
  ninth of a mesh triangle and *inside* the PNSN epicentre location error (5-10 km). It cannot be
  usefully tightened with this catalogue: we are already testing at the resolution of the data.
- Q1 (down-dip offset): 0.97 km/yr = 15.5 km over 16 yr = 0.55 mesh cells — a moderate bound only.
- L (tremor-slip lag): 0.88 d/yr = 14 d over 16 yr.
- **Nucleation latitude: 6.5-11.3 km/yr, i.e. 104-181 km over 16 yr — a WEAK bound.** Event-to-event
  scatter in nucleation latitude is +-1 deg (+-111 km) within a segment, so ~21x more events (or a
  far less noisy nucleation estimator) would be needed to see a 30-km drift.
- Northward fraction: only a 0.25 absolute change between the first and last epoch is detectable.

**Consequence for S10/C4 (act on this).** Our stationarity result does **not** refute Meng et al.
2026 on their own quantities: at 104-181 km / 0.25-in-fraction sensitivity we could not have seen a
change of the size they are likely claiming. What we CAN put against them is the tremor-slip
RELATIVE geometry, headed by Q3 (no depth-offset drift larger than 3.2 km over 16 yr, below the
location error), then Q1 (15.5 km) and L (14 d) — quantities nobody has published a drift for. C4
should therefore be re-framed as "the tremor-slip relative geometry is stationary to within the
location error, while the nucleation/direction quantities that ARE claimed to change cannot be
tested to better than ~100 km with this catalogue" — the second half is itself worth saying, because
it means any nucleation-migration claim (theirs included) needs an explicit power statement and a
short-window control of the kind in PART A.4b.

### Caveats specific to S13
- **Our tremor-event definition is not theirs.** Meng et al.'s "ETS events" are defined by an unknown
  procedure; comparing verdicts requires their event list or their exact rule. Our detector was swept
  over 18 parameter settings, which bounds sensitivity to OUR choices, not to theirs.
- **Central Cascadia is untested.** The 500-detection floor leaves only 10 events nucleating in
  44-46 N, so every central number (including MDD80 = 300 km) is uninterpretable.
- The "south" bin is 4.5 deg wide (39.5-44 N) and contains two belts; segment-level southern
  statistics are mix statistics unless the sub-zone split is carried along.
- Segment membership is assigned by NUCLEATION latitude, so an event that nucleates at 45.9 N and
  runs north counts as central. Events are also not independent of the segment mix: a margin-wide
  nucleation trend can be driven purely by a change in how many events each segment produces.
- The 8-connected grid can merge a northern and a southern episode when the intervening latitudes are
  briefly active; that is bounded by the WIN_KM/DLAT sweep but not eliminated.
- `num_stas` is zero before 2017-01-01, so the >= 8-station control can only test 2017+ and therefore
  cannot address the pre-2017 half of the record at all.
- The residual-bootstrap power analysis assumes exchangeable residuals; if the per-event errors are
  correlated in time (they may be, via shared network state) the true MDD80 is LARGER than quoted,
  i.e. the bounds are if anything weaker than reported.
- MDD80 is read off a 6-12 point power curve by linear interpolation; +-0.02-0.05 in the last digit.
- All PART A quantities are latitude-only (no strike rotation), consistent with S4.

Files written: `scripts/s13_nucleation_direction_power.py`; `results/s13_tremor_events.csv` (191),
`results/s13_nucleation_trends.csv` (20), `results/s13_direction_epochs.csv` (192),
`results/s13_completeness_controls.csv` (333 rows, `block` column names the analysis),
`results/s13_drift_power.csv` (90); `figures/s13_nucleation_direction.png` (4x3: a nucleation vs date
by segment with trends; b migration speed vs date; c northward fraction per epoch; d the same trend
under all four catalogue variants; e daily rate with the 2018-02 marker; f the step-ratio null;
g/h/i power curves for Q1+Q3, L and nucleation latitude; j the sliding-window control; k the event
latitude-time map; l power for the northward fraction). Nothing outside `tremor_slip_geometry/` was
created or modified.

## 2026-09-10 orchestrator review — S13
Model-free tremor-event catalog (191 events from the PNSN grid): nucleation latitude and propagation direction show NO detectable
change 2010-2026 (margin-wide nucleation +4.5 km/yr [-2.4, +12.1], p 0.21; northward fraction 0.49/0.49/0.57, p 0.60; robust to
18 detector settings, bright-only and >= 8-station catalogs, Holm; the Feb-2018 algorithm change is not even a rate step). The one
nominal positive (southern nucleation +18 km/yr in 2017-26) is a low-frequency swing + sub-zone mix effect — killed. POWER: the
stationarity bounds are only tight for the RELATIVE geometry: depth offset MDD80 = 0.2 km/yr = 3.2 km over 16 yr (inside PNSN
location error; un-tightenable with this catalog); down-dip offset 15 km; lag 14 d; nucleation latitude 104-121 km (!); direction
fraction 0.25 absolute. So our null does NOT refute Meng et al. 2026 on their quantities (we could not have seen a plausible
nucleation migration); it does establish that tremor-slip relative geometry is fixed to within one location error. Event SIZE
(peak daily count) grew x1.4 over 16 yr — network growth, flagged. C4 reframed accordingly in RESULTS.

## 2026-09-10 S14 LFE scaling

Ran `scripts/s14_lfe_scaling.py` (OMP_NUM_THREADS=4, seed 20260910; **12.3 min** end to end on the
cold pass, of which **8.3 min streaming 118.6 GB**; 4.3 min on later passes with the per-station
count matrices cached in the session scratchpad). Question, from the publishability memo TO DO (a):
repeat S9/S11's `C ~ R^0.40` with the project's OWN LFE detections as a third, independent tremor
measure (Frank & Brodsky 2019: per-LFE amplitude ~ rate^0.32 in Guerrero). Events, window
[start-10, end+10], footprint, 60-km neighbourhood, R7, Theil-Sen fitting, shifted-tremor nulls and
event dummies are IDENTICAL to S9/S11; loaders copied verbatim from `s11_scaling_robustness.py`.
**Machinery cross-check passed exactly**: the script rebuilds the PNSN measure and returns
b(C, R7) = **+0.401 [+0.181, +0.478]**, digit-for-digit S9/S11.

### What `n_stack` turned out to be — NOT an LFE count
`data/daily_dvv_<STA>_Z_2to4.csv` column `n_stack` is, from `scripts/dvv_roll30cal.py`,
`counts = idx + 1 - searchsorted(dates, dates - 30 d, 'right')`: the **number of DAILY STACKS
(days that produced a stack, i.e. days with >= 20 detections) inside the trailing 30 calendar
days**. It is an ACTIVE-DAY OCCUPANCY count with a **hard ceiling of 30** and a floor of 5
(`--min-stacks`). Verified by recomputing it from each CSV's own dates: reproduces `n_stack`
EXACTLY on 89.3 % and is `<=` it on 100 % of 987k rows in 400 random families — exactly what the
definition predicts, because rows with counts < 5 are never written, so a recount from the CSV can
only undercount. Over all 111.7 M rows: min 5, **median 30**, max 30, 62.5 % at the ceiling; inside
the 27 event windows 59.8 % of family-days are at the ceiling.
=> **n_stack cannot measure a 30-d LFE rate.** Tests 1a/1b were run as specified and return
b = **+0.000 [0.000, 0.000]** in every group (n = 229,435 family-events; pooled OLS +0.004,
r +0.035; 200-shift null mean -0.008 +- 0.086), and the cross-event version (peak n_stack vs peak
event slip, 24,915 families with >= 4 events) also gives **+0.000** (Spearman +0.137). That is
saturation, not physics, and it is reported as such.

### What was streamed (Level B, `results/s14_coverage.csv`)
12 stations chosen for family coverage of the confirmed footprints x epoch, north + central:
WISH 2009-19 (20.5 GB), RADR 2011-19 (14.0), B943 2010-26 (25.8), LEBA 2013-23 (7.2), PR01 2017-26
(1.7), CARB 2019-26 (7.3), ANGUS 2021-26 (4.7) | HEBO 2009-19 (7.5), YACT 2013-26 (19.0), SSO2
2019-26 (3.2), ALSE 2019-26 (6.1), WLOO 2016-26 (1.6).
**118.6 GB of `mf_*.csv` apparent bytes (budget 150 GB), 123 station-year files, 1,584,478,589
detections, 3,342 templates**, reduced on the fly with pyarrow (~210-380 MB/s, IO-bound) to dense
per-(template, day) counts in five cc bins (>= 0.80/0.85/0.90/0.95/0.98) plus summed cc.
DEVIATION FROM THE SPEC, deliberate: **all days were kept, not only event-window days** — the file
has to be parsed either way, the product is tiny, and the full series is what makes the >= 365-d
shifted null and the 30-d rolling counts possible. All 27 events are covered; 145 (station, event)
pairs have the whole window inside the station's recording span, 135 of them yield an exponent.

### The exponents
| measure | b median [IQR] | 90 % CI of median | n |
|---|---|---|---|
| n_stack vs inc30 (30-d, saturated) | **+0.000** [0.000, 0.000] | [0.000, 0.000] | 229,435 fam-events |
| TRUE 30-d LFE count vs inc30 (T1c) | **+0.021** [-0.020, +0.092] | [+0.009, +0.039] | 142 |
| daily LFE count vs R7 (T2) | **+0.032** [+0.000, +0.105] | [+0.021, +0.049] | 135 |
| daily LFE count, floor removed (T6) | **+0.049** [-0.031, +0.148] | [+0.028, +0.070] | 135 |
| daily LFE count, floor de-attenuated | **+0.130** [+0.001, +0.391] | [+0.091, +0.215] | 79 |
| per-family, de-attenuated | **+0.206** | [+0.199, +0.212] | 13,898 |
| **PNSN daily count vs R7 (reference)** | **+0.401** [+0.181, +0.478] | [+0.310, +0.452] | 27 |
| PNSN, same floor removal (control) | +0.412 [+0.166, +0.476] | [+0.351, +0.415] | 145 |
| PNSN, same de-attenuation (control) | +0.528 [+0.174, +0.737] | [+0.428, +0.624] | 141 |

Sub-groups of the daily LFE exponent: north +0.048, central +0.011; 2010-14 +0.077, 2015-19 +0.014,
2020-26 +0.035; pooled with station-event dummies OLS +0.051 / Theil-Sen +0.051 (7,156 day-rows,
r +0.088). Against the 200-draw >= 365-d shifted-LFE null the sign is decisive but the amplitude is
not: null mean +0.006 +- 0.144, z median +0.284, **97/135 positive (sign p 4e-07)**, yet only
16/135 station-events reach null p < 0.05. Lower vs upper half of log10 R7: +0.008 vs +0.107
(upper - lower +0.090, Wilcoxon p 0.008) — no flattening at high rate, the same non-saturating
behaviour S11 found. Paired against the PNSN exponent on the same events:
**b(LFE) - b(PNSN) = -0.290 [-0.445, -0.094], Wilcoxon p 3.4e-19** (event-level -0.353, p 2.8e-07).

### WHY they disagree — the diagnosis (Test 5), and it is a measurement fact, not physics
The densify threshold is **cc >= 0.80 with a 6-s minimum gap** (`fleet_station.sh` ->
`densify_gnw_gpu.py --threshold 0.8 --min-gap-s 6`). At that threshold the matched filter fires
**~108 detections per template per day** (median over stations of the mean daily total,
~32,400 detections/day/station over ~300 templates) in QUIET times. Consequences, all measured here:
- **ETS contrast**: mean daily detections inside event windows / outside = **x1.24 median over
  stations** (range x0.83 SSO2 to x1.81 PR01; WISH 1.47, CARB 1.04, YACT 1.15, WLOO 1.72) against
  **x4.9 for the PNSN count in the SAME 40-km neighbourhood** (x12 on per-event medians).
- **Agreement with PNSN**: Spearman(daily LFE total, daily PNSN count, same neighbourhood) =
  **+0.088 median** (range -0.002 .. +0.179).
- **cc ladder** (b | median detections per window): 0.80 +0.032 | 1.1e6; 0.85 +0.029 | 2.9e5;
  0.90 +0.020 | 3.5e4; 0.95 +0.000 | 1183; 0.98 +0.000 | 149. Raising the threshold thins the
  counts by 4 decades without ever recovering the exponent — at cc >= 0.95 most days are empty and
  the Theil-Sen slope collapses to 0. (Per station the high-cc ETS contrast does improve — WLOO
  1.72 -> 10.1, CARB 1.04 -> 1.87, YACT 1.15 -> 3.74 — but on far too few detections.)
- Restricting to **causality-CERTIFIED families** (the project's own LFE verdict, minus the
  229-family user blocklist) changes nothing: +0.032 vs +0.032 (paired +0.002). Nor does a 60-km
  family radius (+0.032, paired +0.000).
So the `mf_*.csv` files are **a densify product dominated by a nearly constant false-alarm rate,
not an LFE catalogue**. The pipeline extracts LFE signal from them by STACKING
(`build_long_window_3comp.py`), not by counting — and that same fact explains the n_stack
saturation, because the `--min-det 20` gate is met on essentially every day.

### The floor correction, and how far it gets
For `N = F + k R^b` with a constant floor F, `dlogN/dlogR = b (r-1)/r` with
`r = mean(N in window)/mean(N in [start-70, start-10])`, so `b_true = b_obs r/(r-1)`. Applied:
LFE r median 1.25 -> **+0.130 [+0.091, +0.215]**; the identical operation on the PNSN count
(r = 3.78) gives +0.528 against its true +0.401, i.e. the correction over-shoots by ~1.3x, so the
LFE value is really ~0.10-0.16. Stratifying per family by responsiveness gives a clean monotone
ladder (raw | floor-removed): ets_ratio < 1.1 **+0.026 | +0.070** (n 13,379); 1.1-1.5
**+0.070 | +0.116**; 1.5-2.0 **+0.104 | +0.147**; 2.0-3.0 **+0.091 | +0.192**; > 3
**+0.000 | +0.000** (counts too sparse). Only 28.7 % of family-events have ets_ratio > 1.5; those
give raw +0.078, floor-removed +0.117, de-attenuated +0.154. **The exponent rises monotonically
with the responsiveness of the measure and never reaches the PNSN 0.40** — the same
errors-in-variables ladder S11 demonstrated for R1/R3/R7 and S11b for GNSS proxy SNR, here on the
tremor side instead of the slip side.

### Per-LFE amplitude — NOT testable (Test 3)
`mf_*.csv` carries (template, time, cc, station). `cc` is the matched-filter CORRELATION, bounded
in [0.80, 1] and normalised by the running trace energy; it is not an amplitude. **Frank &
Brodsky's per-LFE amplitude ~ rate^0.32 cannot be reproduced from these files** — it needs the
per-detection amplitude, i.e. re-reading the waveforms. Reported as a weak proxy only: the exponent
of the daily MEAN cc on log10 R7 is **-0.000 [-0.001, +0.000]** (mean cc 0.8300) — flat, as
expected for a threshold-defined quantity.

### VERDICT
**S14 does NOT deliver an independent third value for the exponent.** What it delivers:
1. **Sign and sub-linearity are confirmed by a third, fully independent detector**: the daily LFE
   count is a monotone increasing, strongly sub-linear function of the modeled slip rate
   (97/135 station-events positive against a >= 365-d shifted null, p 4e-07; **b < 1 in 100 % of
   fits**; b = 1 excluded by every CI), with no flattening at high rate.
2. **The VALUE is not measurable from these files.** The measured 0.03 is a floor-attenuated lower
   bound; every correction moves it up (0.03 -> 0.05 -> 0.13-0.21) and it is still climbing when the
   corrections run out of signal, so it is consistent with, but cannot confirm, 0.3-0.4.
3. **A new, checkable fact about the project's own data**: the cc >= 0.80 densify catalogues have an
   ETS contrast of only x1.05-1.8 and rho ~ 0.09 with PNSN tremor. Anything that counts these
   detections (rather than stacking them) is measuring the false-alarm rate.
=> For the paper, the third tremor measure must come from a PROPER LFE catalogue (a published
Cascadia LFE catalogue, or a re-detection with a MAD-based per-template threshold + declustering +
a network-coincidence requirement), not from the densify by-product. The Frank & Brodsky comparison
in the memo's "interpretation to develop" cannot be closed with the current files.

### Caveats specific to S14
- **LFE families are station-specific templates.** A "family" is one station's template cluster;
  the same physical source appears under different ids at different stations, so family counts are
  not comparable across stations and everything is per-station with station-event dummies.
- **Detection threshold.** Everything above is at the pipeline's cc >= 0.80 / 6-s-gap setting; the
  ladder to cc >= 0.98 is reported but starves the statistics.
- **30-d smoothing.** `n_stack` and `inc30_mm` are both trailing-30-d quantities, so T1a/T1b are
  low-pass on both axes even before the ceiling; T1c uses a trailing-30-d LFE count against the
  same `inc30`, so it inherits the same bandwidth (and the same floor problem).
- **Events, windows and footprints still come from the slip-patch catalogue** (as in S9/S11/S11b);
  only the tremor measure changed.
- **Coverage is uneven in time**: 24 station-events in 2010-14, 40 in 2015-19, 71 in 2020-26; HEBO
  contributes only 412 days with detections over 2009-2019.
- **The T1a null is subsampled** (2,177 family-events x 50 shifts) — the full null is ~10^8 fits.
- The per-station daily-count matrices are cached OUTSIDE the project, in the session scratchpad
  (`.../scratchpad/s14_cache/<stem>_v2.npz`); nothing outside `tremor_slip_geometry/` was created
  or modified inside the repository.

Files written: `scripts/s14_lfe_scaling.py`, `results/s14_lfe30_exponents.csv` (24 rows =
T1a/T1b/T1c group summaries), `results/s14_lfe_daily_exponents.csv` (145 station-event rows +
102 summary rows x 77 columns), `results/s14_coverage.csv` (12 stations x what was streamed + the
Test-5 diagnosis columns), `figures/s14_lfe_scaling.png` (9 panels: A exponent distributions with
the shifted null, B by segment, C by epoch, D one event as a time series showing the dilution,
E paired LFE vs PNSN, F the cc-threshold ladder, G the n_stack ceiling, H the per-station ETS
contrast vs PNSN, I summary).

## 2026-09-10 orchestrator review — S14
The project's densify detections (data/mf_*.csv, 1.58 billion detections streamed over 12 stations) are NOT an LFE event catalog:
at cc >= 0.80 the matched filter fires ~108 times per template per day in quiet times, ETS/non-ETS ratio only x1.24 (PNSN x4.9
in the same neighbourhood), Spearman with PNSN daily count +0.09; n_stack in the daily dv/v tables is the number of daily STACKS in
30 d (ceiling 30), not a detection count. The LFE-count exponent (0.03, floor-attenuated lower bound) is therefore uninformative
on the value, though qualitatively consistent (b < 1 in 100 % of 135 station-event fits, positive vs null p 4e-7, no flattening
at high rate). The PNSN reference is reproduced digit-for-digit by the same code (validation). A third independent tremor measure
needs a proper LFE catalog (published Cascadia LFE catalog, or re-detection with per-template MAD thresholds + declustering +
network coincidence). FLAG for the main project: the mf_*.csv detections are stacking material, not events — "per-event science
needs no re-detection" (memory 2026-08-31) is wrong for anything that COUNTS detections.

## 2026-09-10 S15 launched (amplitude snapshot check + candidate-count exponent) — see PLAN.md round 5

## 2026-09-10 S15 amplitude snapshot + candidate-count exponent

Two scripts, `OMP_NUM_THREADS=4`, seed 20260910: `scripts/s15_amplitude_check.py` (TASK A, **66 s**)
and `scripts/s15_candidate_scaling.py` (TASK B, **83 s**). Both copy their loaders/fitting verbatim
from `scripts/s11_scaling_robustness.py`; TASK A additionally copies the component labeller verbatim
from `fault_tomography/inversion/sse_slip_patches.py` (that script was NOT modified or re-run).
Nothing outside `tremor_slip_geometry/` was created or modified (TASK A caches `A_new` + `R7` in the
session scratchpad for TASK B).

### TASK A — the "amplitude snapshot mismatch" was a TIME-AXIS ARTEFACT; nothing changes

**The shipped `A_raw` IS the 2026-08-18 solution.** Rebuilding the amplitudes exactly as
`sse_trailing_products.py` does (`component_amplitudes()` -> `groupby(level=0).mean()` ->
daily `reindex` -> `interpolate(limit_area="inside")`) gives `A_new` that equals the stored
`A_raw` to **float32 round-off**: per component Pearson r = **1.000000000** (13/13), rms difference
/ column sd **2.5e-8**, max |difference| **3.6e-6** in amplitude units (max |d|/sd 5e-7), and
**100.0 % of the 93,197 entries are bit-identical after casting `A_new` to float32** — which is the
dtype the npz stores.

**Where S12b's "0.87-0.99" came from.** `sol.component_amplitudes()` has one row per SOLUTION
EPOCH, and those epochs are neither unique nor daily: **7169 epochs but only 5271 distinct dates
(1898 duplicated dates), against a 7169-day daily grid (1898 days missing from the epoch list)**.
The two arrays have the same LENGTH by coincidence — the duplicate count and the missing-day count
happen to be equal, so the running index offset stays within **-1..+1 days** over the whole record —
and an index-by-index comparison therefore compares series shifted by about one day. Reproducing
exactly that comparison here returns **r = 0.870..0.993 across the 13 components**, i.e. the
reported range, digit for digit. So the S12b note ("comes from an EARLIER snapshot and correlates
only 0.87-0.99 ... so it was NOT used") is **WITHDRAWN**: there is no earlier snapshot and no stale
model anywhere in S1-S14. (Side note for S12b itself: it used `Amp` on `sol.dates` AS IF that were
a daily timeline and mapped tremor onto it with `searchsorted`, so its own series carry the same
<= 1-day distortion; that is a sub-sampling-level effect on 30-d trailing quantities and does not
touch its conclusions, but it is the same root cause.)

**(2) The 27-event set is unchanged.** With `A_new`: **118,778 patch-days above 5 mm/30 d, exactly
the same count as `A_raw` and as the shipped `inc30_mm`, with ZERO patch-days disagreeing**
(max |inc30 difference| 1.0e-6 mm); the labeller returns **2196 connected components, the same
number the shipped label array has, with the SAME ids**. All **27/27 confirmed events PASS** the
required test (a > 5 mm/30 d component with start AND end within +-10 d and footprint overlap
>= 70 %): every event matches its own id with **start/end offsets +0/+0 d, footprint overlap
1.0000 and Jaccard 1.0000**, identical triangle counts (23-310).

**(3) The count exponent is unchanged.** b(C, R7) per event with `A_new` versus with `A_raw`
(stored footprints, everything else as S9/S11): median **+0.401307 [+0.1805, +0.4775]**, 90 %
bootstrap CI of the median **[+0.3105, +0.4524]**, **27/27 below 1** — in BOTH cases, and the
recomputed value reproduces the stored `results/s11_exponents.csv` b(C,R7,all) median
**+0.401307** digit for digit. Paired difference **-8e-10 median, max |d| 1.0e-7** (float32
round-off; Wilcoxon p 0.53). Rebuilding the footprints from the new labelling as well changes
nothing (same median, max |d| 1.0e-7).

**PLAIN VERDICT (A).** **YES — the 27-event set and the exponent are unchanged with the current
amplitudes, because the current amplitudes ARE the shipped ones.** The publishability-memo TO DO (d)
is closed with a null result, and the S12b caveat should be struck from RESULTS.md.

### TASK B — the picker-ranked candidates: sub-linear, b ~ 0.68, but an AMPLITUDE-weighted re-read of the PNSN windows

What the measure is (read from `scripts/discover_nllb_pnsn_driven.py` + `lfe_features/score_candidates.py`
+ `scripts/adaptive_cand_threshold.py`): for every PNSN tremor window inside a station's discovery
bbox (+-0.9 deg lat, +-1.35 deg lon) a **300-s** segment of that station's 2-8 Hz record is cut; every
Hilbert-envelope peak above **3 x the window's median envelope**, separated by >= **6 s**, is emitted
as one candidate; `time` is the peak time AT THE STATION, `lat/lon` is the PNSN window's own
epicentre (the source cell, not relocated), `p_lfe` is the picker's P(LFE), and `_cand_filtered`
keeps the top-ranked ~30,000 per station (rank-based, P floor 0.05).
**Loaded: 209 stations, 6,713,290 candidates, 2010-01-10..2026-05-27.** Base product: daily counts
per station and per **0.25-deg source-cell latitude band** — 209 x 44 bands (39.50-50.50 N) x 5991
days, **188,344 occupied (station, band, day) cells**. Candidates per year grow **136k (2010) ->
1,084k (2025)**, stations active 48 (2010) -> 176 (2025); median station record 7.0 yr.

Per event (window [start-10, end+10], footprint from `slip_patches_5mm.npz`, source cell within
60 km of the footprint, station within 150 km of the area-weighted footprint centroid, R7 from
`A_new`): **stations contributing median 26 [11, 70]**, candidates in window median **68,275
[8,496, 360,454]**, PNSN detections in the same neighbourhood median **7,209 [293, 14,456]**,
candidate/PNSN ratio 3.9-57. All 27 events are covered (2010-2025; the only multi-year window is
1317, 2015-2016).

**Exponents** (Theil-Sen of log10(measure + eps) on log10 R7, per event; median [IQR], 90 %
bootstrap CI of the median, pooled = event dummies):

| measure | b median [IQR] | boot90 CI | pooled OLS/TS | b<1 | vs shifted null (pooled / per-station) |
|---|---|---|---|---|---|
| **C_cand** (candidate count) | **+0.677** [+0.324, +0.921] | **[+0.521, +0.770]** | +0.585/+0.676 | 22/27 | 27/27 p 7.5e-09 / 26/27 p 2.1e-07 |
| W_cand (p_lfe-weighted) | +0.801 [+0.483, +1.062] | [+0.703, +0.822] | +0.682/+0.792 | 19/27 | 27/27 p 7.5e-09 / 26/27 p 2.1e-07 |
| C_q75 (top-quartile p_lfe) | +0.808 [+0.527, +0.973] | [+0.571, +0.918] | +0.648/+0.749 | 21/27 | - |
| C_s10 (10 nearest stations) | +0.648 [+0.429, +1.095] | [+0.529, +0.862] | +0.646/+0.737 | 20/27 | - |
| M_cand (multiplicity C_cand/C_pnsn) | **+0.184** [+0.114, +0.323] | [+0.127, +0.248] | +0.267/+0.256 | 27/27 | - |
| C_pnsn (S9/S11 reference) | +0.401 [+0.181, +0.478] | [+0.310, +0.452] | +0.328/+0.397 | 27/27 | - |

Nulls: 200 rigid >= 365-d wrapped shifts of the near-candidate dates with the slip series fixed,
in two flavours — **pooled** (S9/S11's null: null mean +0.001 +- 0.503, z median +1.09, 27/27
positive, sign p 7.5e-09) and **per-station** (each station's dates wrapped inside ITS OWN record
span, so the null cannot borrow another station's active years: null mean +0.003 +- 0.267, 26/27
positive, sign p 2.1e-07, 19/27 events at null p < 0.05). North +0.697 / central +0.542.

**The measure is SEMI-INDEPENDENT and the decomposition says exactly how.** Candidates exist ONLY
inside PNSN tremor windows, and the numbers make that concrete: **Spearman(daily C_cand, daily
C_pnsn) = +0.969** within events (S14's `mf_*.csv` detections gave +0.09), and
**b(C_cand vs C_pnsn) = +1.394 [+1.223, +1.604]** — the candidate count is a super-linear,
essentially monotone transform of the PNSN count in the same neighbourhood. So C_cand is NOT a new
detector; what is new is the station-side **multiplicity**, the number of envelope peaks each
PNSN window produced, and its exponent is **b(M_cand) = +0.184 [+0.127, +0.248], 27/27 events
below 1** (pooled +0.27; b(C_pnsn) + b(M_cand) = 0.401 + 0.184 = 0.585 = the pooled OLS for
C_cand, so the decomposition closes).

**That +0.18 is the useful result: it is a THIRD, amplitude-side confirmation of S11's Test 4.**
S11 found b_E - b_C = **+0.222 [+0.021, +0.387]** (radiated PNSN energy versus count, 2017+, 18
events) and b(E/C vs R7) = +0.231, i.e. "counts under-represent strong tremor". S15B reproduces
that offset from a completely different amplitude proxy, on all 27 events including the pre-2017
ones where PNSN has no energy column: the paired per-event difference is
**b(C_cand) - b(C_pnsn) = +0.225 [+0.160, +0.323], Wilcoxon p 2.0e-06, 25/27 positive**
(Spearman of the two per-event exponents rho +0.83). The candidate exponent therefore lands on the
PNSN **energy** exponent (0.64), not the PNSN **count** exponent (0.40) — which is what an
envelope-peak count above a fixed x3 background threshold should do.

**Controls (all pass).**
- **Count magnitude is not the explanation.** The candidate count is ~10x the PNSN count, so the
  +1 offset compresses low-count PNSN days more. (a) With eps = 0 instead of +1: b(C_cand) +0.677,
  b(C_pnsn) +0.412, paired **+0.204, p 2.5e-06**. (b) Binomially **THINNING** the candidate count
  down to the PNSN mean per event (20 draws): b = **+0.637 [+0.504, +0.688]**, still
  **+0.173 above PNSN (p 6.7e-06, 25/27 positive)**.
- **Network growth is not the explanation**, although it is large (Spearman(n_sta, year) = +0.82,
  p 1.5e-07): b by epoch **+0.612 / +0.700 / +0.687** (Kruskal p 1.00), Spearman(b, n_sta)
  **-0.03 (p 0.87)**, Spearman(b, year) **-0.02 (p 0.93)**, and restricting every event to its
  **10 nearest contributing stations** gives **+0.648** (paired +0.084, p 0.058).
- **The 300-s window ceiling is not reached.** Every window candidate was matched back to its own
  PNSN window (merge_asof on the exact source cell, 300-s tolerance; **100.0 %** matched). A
  (station, window) pair can hold at most **50** candidates (300 s / 6 s); the observed mean
  multiplicity is **3.03**, p99 **11**, max **28**, and the fraction of (station, window) pairs AT
  the ceiling is **0** in every event. Saturation would in any case LOWER b, not raise it.
- **No change of slope across the rate range**: upper minus lower half of log10 R7 is -0.017
  (C_cand, p 1.0), -0.112 (W_cand, p 0.96), -0.118 (C_q75, p 0.75), +0.257 (C_s10, p 0.98),
  +0.099 (C_pnsn, p 0.35) — none significant, 12-14 of 25 positive.
- **Picker rank matters in the expected direction**: keeping only the event's top-quartile p_lfe
  candidates RAISES b to **+0.808** (paired +0.107, p 0.0046), and the p_lfe-weighted count gives
  +0.801. A rank-blind artefact would not order this way.
- **ETS contrast, the S14 diagnostic that condemned `mf_*.csv`**: per-station mean daily candidate
  count inside the 27 windows / outside = **3.99 [0.91, 15.97] (10-90 %), > 1 for 178/206
  stations**, against **3.14** for the PNSN count inside the SAME station bbox on the same days
  (paired candidate - PNSN **+0.56, Wilcoxon p 3.9e-17, higher in 144/206 stations**), margin-wide
  PNSN 2.01, and S14's `mf_*.csv` **1.24**. So these candidates ARE ETS-responsive, unlike the
  cc >= 0.80 densify detections — as they must be, being conditioned on PNSN windows.

**PLAIN VERDICT (B).** **YES — the picker-ranked candidate count reproduces a sub-linear exponent:
b = +0.68, honest range 0.5-0.8 (90 % bootstrap CI of the median [+0.52, +0.77]; per-event IQR
+0.32 to +0.92; 22/27 events below 1; pooled event-dummy fit +0.59; p_lfe-weighted +0.80,
CI [+0.70, +0.82]).** b = 1 is excluded by the bootstrap CI for every candidate variant, and the
exponent is positive against both shifted nulls (27/27 and 26/27 events). **BUT it is not a third
independent tremor measure**: with rho = 0.97 against the daily PNSN count and
b(C_cand vs C_pnsn) = 1.39, it is an amplitude-weighted re-reading of the same PNSN windows. Its
genuine content is the multiplicity exponent **+0.184 [+0.127, +0.248]**, which independently
reproduces S11's energy-minus-count offset (+0.222) on all 27 events. Read for the paper as:
*"a station-side envelope-peak count, which weights strong tremor more heavily than the catalogue
count does, gives b = 0.68 — consistent with the PNSN energy exponent 0.64 and inconsistent with
b = 1 — and confirms that the count under-represents the high-rate end by ~0.2 in exponent."*
It does **NOT** close the memo's TO DO (a) for a genuinely independent third measure (a published
Cascadia LFE catalogue, or a MAD-thresholded declustered re-detection, is still needed), and it
does **NOT** independently confirm the *value* 0.3-0.4 for a pure detection count.

### Caveats specific to S15
- **TASK A is a null result by construction** once the alignment is understood; it cannot detect a
  change of solution DIRECTORY (a different `data/<date>` snapshot), only that `A_raw` matches the
  2026-08-18 one that S1-S14 and S12b both load.
- **The candidate measure is conditioned on PNSN** (PNSN sets WHEN to look). It can never see
  tremor PNSN missed, and its footprint/latitude information is PNSN's, not the station's — the
  station contributes only the peak TIMES and the peak COUNT.
- **`p_lfe` is a rank, not a calibrated probability**: `adaptive_cand_threshold.py` takes whatever
  P gives ~30,000 candidates per station (floor 0.05), and absolute P does not transfer across
  instrument classes. W_cand and C_q75 are therefore rank-weighted, not probability-weighted.
- **`snr` is an envelope ratio to the window median**, so the candidate count is amplitude-sensitive
  by construction; that is the point of the measure but it also means it is NOT a count of events.
- **The per-station null keeps only a median 42 % of the near candidates** (median 25 stations per
  event pass the >= 3-yr span and window-inside-span requirements), so it is stricter but noisier.
- **Events, windows and footprints still come from the slip-patch catalogue**, as in S9/S11/S11b/S14.
- Station coordinates were merged from `data/broadband_fleet_order.csv`,
  `data/station_slab2_depth.csv` and `data/candidate_stations_post2020.csv` (all 209 covered);
  three stations have no ETS/quiet ratio (no quiet or no ETS day inside their record).
- Pandas trap observed: `to_datetime`/parquet give `datetime64[us]` here, so every epoch/int64
  conversion casts to `datetime64[ns]` first.

Files written: `scripts/s15_amplitude_check.py`, `scripts/s15_candidate_scaling.py`,
`results/s15_amplitude_check.csv` (13 rows = components x aligned/misaligned comparison),
`results/s15_events_check.csv` (27 rows), `results/s15_exponent_check.csv` (27 events + summary),
`results/s15_cand_exponents.csv` (27 event rows + 6 measure summaries + 4 paired + 5 growth +
1 magnitude-control row, 119 columns), `results/s15_cand_coverage.csv` (27 event + 937
event-station + 209 station rows), `figures/s15_candidate_scaling.png` (A exponent distributions
with both nulls, B-D three example event tracks, E paired per-event PNSN vs candidates, F coverage
over time, G lower/upper-half saturation check, H per-station ETS/quiet contrast, I the
decomposition boxplot).

## 2026-09-10 orchestrator review — S15
(A) The "earlier amplitude snapshot" caveat (S12b) is WITHDRAWN: the shipped A_raw equals the 2026-08-18 solution to float32 round-off
(r = 1.000000000, 100 % bit-identical after the float32 cast). S12b compared sol.component_amplitudes() rows (7169 epochs with 1898
duplicate dates) index-by-index against the daily grid — a <= 1-day misalignment that reproduces exactly its 0.87-0.99. Nothing in
S1-S14 used a stale model; the 27 events and the exponent (0.4013, digit for digit) are unchanged. S12b's own time-domain numbers
carry the <= 1-day distortion (its verdict does not depend on them at that level).
(B) Picker-ranked candidate counts give b = 0.68 [0.52, 0.77] (22/27 below 1, b = 1 excluded) but are NOT independent of PNSN
(daily Spearman 0.97): they are an envelope-peak count with multiplicity ~ rate^0.18, so they land on the PNSN ENERGY exponent
(0.64) rather than the count exponent (0.40) — an independent reproduction of the energy-count offset (+0.22) on all 27 events, not a
third measure. The truly independent third measure still requires the deep-learning LFE catalog (lfe_catalog_repair E-plan).
## 2026-09-10 S16 launched (F&B aggregation match + paper figures)

## 2026-09-10 S16 (resumed after crash) — F&B aggregation match + the paper figure set

Session context: the pod crashed at ~21:03 with S16 in flight. On disk from the crashed run:
`scripts/s16_fb_aggregation.py` and its two outputs `results/s16_fb_aggregation.csv` (116 rows) and
`results/s16_fb_events.csv` (18 rows), both written 21:02. **Checked first whether the script
intended any output that was missing: it does not** — its OUTPUT block names only
`s16_fb_aggregation.csv` (plus the per-event table written at line 829); there is no log, summary or
figure output in the script, and the CSVs are complete (they include the diagnostic rows written
last, `catalog_EM_relation` and `eps_S11_E_check`, so the run reached the end). Only the stdout
report was lost, and every number in it is recoverable from the CSVs. **The script was therefore NOT
re-run and NOT modified.** All S16 files added in this session are new and prefixed `s16_`.

### PART A — the Frank & Brodsky (2019) aggregation match  (`results/s16_fb_summary.md`, 209 lines)

The like-for-like statistic (F&B's own aggregation: MEDIAN per-detection amplitude per bin) on the
18 amplitude-era events, with S9/S11's windows, R7 regressor, Theil-Sen fitting and 200-draw
>= 365-d shifted-tremor null:

| aggregation / estimator | amplitude exponent | vs F&B 1/3.17 = 0.315 |
|---|---|---|
| forward, 1-d bins (our convention) | **+0.121** [+0.032, +0.151], 18/18 < 1 | 2.6x LOW, CI excludes |
| forward, 3-d bins | +0.166 [+0.100, +0.213] (n 16) | still low |
| forward, 7-d bins | **+0.227 [+0.177, +0.436]** (n 11) | **CI CONTAINS 0.315** |
| cross-event, 1 point per event | **+0.296** (OLS +0.247, r 0.485; no CI) | indistinguishable |
| RMA (direction-symmetric, daily) | +0.271 [+0.046, +0.301]; beta_RMA **2.54** | brackets 3.17 |
| reverse = F&B's own direction (daily) | beta **+1.009** [+0.198, +1.309]; 18/18 events < 2.88 | 3x LOW in beta |

Median per-event |r| of log A1 on log R7 = **0.319 (r^2 0.114)**, so forward and reverse differ by
~1/r^2 ~ 9x; that estimator ambiguity, not a physical disagreement, is what separates 0.12 from
0.99. **VERDICT: the F&B amplitude-slip relation DOES reproduce in Cascadia in sign and
sub-linearity (18/18 events b < 1; positive vs the shifted null, sign p 0.015, and p 6.6e-04 for
the magnitude-based variants), and its VALUE is reproduced once the aggregation granularity and the
estimator are matched (0.23-0.30 vs 0.315); at daily resolution with the forward estimator it is
attenuated to 0.12.** Supporting numbers, all from the CSVs:
- Catalogue constraint measured in the script: `log10 E = 1.781 M + 0.007 log10(dur) + 2.904`,
  residual sd **0.052 dex**, r 0.995 over 482,659 2017+ detections -> `energy` and `magnitude` are
  ONE observable; A2 = 10^(0.5 M) is **0.561x A1 by construction** and must never be quoted as an
  independent value (its 0.059 is 0.561 x 0.106).
- The **S10 "summed-amplitude reading" is retired**: `sqrt(daily summed energy)` = b(E_sum)/2 =
  +0.344 identically, so its apparent match to 0.315 is arithmetic, not a measurement. The LITERAL
  summed amplitude (daily sum of sqrt E) gives **+0.561**, close to b_C + b_A1 = 0.370 + 0.121
  (paired difference between the two summed variants +0.212, p 1.9e-04).
- Whole distribution responds, not just the tail: E_p10 +0.151, p25 +0.211, med +0.243, p75 +0.257,
  p90 +0.257 (`E_p90 - E_med` +0.013 p 0.074; `E_p10 - E_med` -0.058 p 0.004). `E_med` and
  `E_mean` are interchangeable (paired -0.001, p 0.80), so F&B's median statistic == S11's
  energy-per-detection in this catalogue. >= 8-station requirement LOWERS the amplitude exponent
  (-0.031, p 0.007) while it RAISED the count exponent (S11 +0.074) -> not a threshold artefact.
- Against our own exponents on the SAME 18 events (aggregated here from `s16_fb_events.csv`):
  b(A1) - b(C) = **-0.262** [-0.292, -0.155] CI90, Wilcoxon p 1.5e-05, 1/18 positive;
  b(E_med) - b(C) = -0.128 (p 7.6e-05). Per-event rho(b_A1, b_C) = +0.88, rho(b_A1, b_E_sum) =
  +0.92 -> three readings of ONE relation, not independent measurements. **So F&B's 0.32 must NOT
  be compared with our count exponent 0.40** (numerical coincidence between a per-event amplitude
  exponent and a rate exponent).
- Decomposition for the paper: count R^0.37 + per-detection amplitude R^0.12 -> summed amplitude
  R^0.56; count R^0.37 + per-detection energy R^0.23 -> summed energy R^0.69 (18-event medians;
  the 27-event count exponent is 0.40).
- Validation of the crashed run: `C` and `E_sum_epsS11` reproduce `s11_energy_events.csv`
  (`b_C_R7`, `b_E_R7`) per event to **max |diff| 0.0e+00**, and `E_mean` reproduces `b_EperDet_R7`.
- Biggest remaining obstacle (stated in the summary): PNSN `energy`/`magnitude` describe a
  150/300-s tremor window, i.e. a windowed envelope amplitude that mixes source amplitude with the
  number of overlapping bursts; F&B's A is a per-LFE displacement amplitude. A strict comparison
  still needs a per-LFE amplitude catalogue (the S14/S15 LFE-catalogue gap).

### PART B — paper figure set (`scripts/s16_paper_figures.py`, 300 dpi, Okabe-Ito palette)

- `figures/s16_paper_fig1.png` (a-e): tremor-occupancy map + the 27 event centroids; thrust vs
  back-slip Q1/Q3 (4.7/0.9 km vs 18.6/5.3 km, n 27 vs 40); Q3 vs date with the Theil-Sen bound
  (-0.04 km/yr [-0.16, +0.09], MDD 0.2 km/yr); the completeness-growth control; per-year depth
  occupancy (north p10/p50/p90 + slip-weighted depth, trend +0.003 km/yr).
- `figures/s16_paper_fig2.png` (a-f): three example daily tracks; the 27 per-event count exponents
  with the null band and b = 1; the four-aggregation ladder with the F&B line; the model-free S11b
  comparison; segment/epoch constancy; and the F&B aggregation/estimator ladder.
- `figures/s16_paper_fig3.png` (a-f): |L| <= 5 d in 63 % vs 18.5 % null; v_s vs v_t (78 % direction
  agreement, ratio 0.18); the 10-bin depth-vs-phase profile against its null; the event-stacked
  migration; the radius/anchor control with null p; S9's 5-bin version.
- `figures/s16_paper_fig4.png` (a-f): r30 representability map with the central-Oregon cluster;
  r30/misfit/station-distance vs latitude; slip-tremor rank correlation per segment and epoch
  (central 0.40-0.46 vs 0.70-0.77); the GNSS data-vs-model amplitudes (model/data 0.63 central vs
  0.77 north); inter-ETS confirmation by definition (A0 leaky 10/11, A1/A2 strict 0-1 bursts); slip
  per tremor day along strike by epoch.
- `results/s16_table1.csv` (27 x 31): the 27 events with Q1/Q3/potency/log10Q5/z1/z3 (S1), the
  snapshot-check PASS/overlap/jaccard/peak (S15A), b_count + null p (S9), b_energy_sum (S11),
  b_energy_median / b_amplitude_FB / beta_FB_direction (S16) and L/v_t/v_s/dir_agree (S4), rounded
  to paper precision.
- `results/s16_table2.csv` (57 x 23): every exponent with IQR, 90 % bootstrap CI, n_below_1, pooled
  event-dummy fits, shifted-null moments and sign p, plus `source_step`/`source_csv`/`note`.
  Blocks: model-based S11, bandwidth ladder, model-free S11b, station-side S15B (flagged NOT
  independent of PNSN), F&B aggregation S16, F&B bin ladder, F&B direction, stationarity.
- `results/s16_captions.md`: caption drafts naming, panel by panel, exactly what is plotted and the
  source CSV/step, plus the column provenance of both tables.

### The ONE quantity that had to be regenerated (and its validation)

S9/S11/S16 all saved only per-event exponents, never the daily series, and the crashed session's
scratch cache was wiped with the pod, so Fig 2a had nothing to plot.
`scripts/s16_daily_tracks.py` (32 s, CPU, OMP_NUM_THREADS=4) recomputes ONLY the daily tracks with
the loaders/geometry/measures copied verbatim from `s16_fb_aggregation.py` (which copied them from
s11 -> s9); the sole difference is that the near-tremor selection is not pre-filtered to 2017+, so
the same code serves the 27-event count fits and the 18-event amplitude fits (all 18 windows start
after 2017-02-14, so the count series are identical either way). Outputs
`results/s16_daily_tracks.csv` (1549 event-days x 13) and `results/s16_track_validation.csv`
(27 x 27). **Validation: every recomputed per-event exponent reproduces the stored value to
max |diff| 8.3e-17** — b_C vs s9_rate_scaling (27), vs s11_energy_events (18), vs s16_fb_events
(18); b_E_sum(epsS11) vs s11 b_E_R7; b_E_sum(eps 0), b_E_med and b_A1 vs s16_fb_events.

### Decisions taken in this session
1. Did not re-run `s16_fb_aggregation.py` (nothing missing; re-running would only reprint stdout).
2. Regenerated the daily tracks rather than leaving Fig 2a without data, because the memo's
   "27 event panels of log tremor vs log slip rate" cannot be built from saved results at all.
3. Reported A1 = sqrt(median energy) as the primary F&B proxy and kept A2 = 10^(0.5 M) only with
   its 0.561 construction factor stated (as the script itself insists).
4. Did NOT present the picker-ranked candidate exponent as a third independent measure anywhere in
   the figures, summary or tables; in Table 2 it carries the explicit note "NOT independent of PNSN
   (daily rho 0.97)" (orchestrator review of S15).
5. Table 2's segment/epoch rows are medians of SAVED per-event exponents, labelled
   "aggregated here" with no CI, because S9/S11/S16 never bootstrapped the sub-groups.

### Caveats / open items for a human
- **Fig 2a is illustrative**: three events (1863 north 2017, 1601 north 2010, 1361 central 2019)
  chosen to span the exponent range; a 27-panel supplementary figure is possible now that the
  tracks are saved but was not built.
- **No CI for the cross-event fit** (+0.296) — the S16 script does not bootstrap that level, so the
  headline "at F&B's granularity we get 0.30" rests on 18 points with a CI only at the 7-d level
  ([+0.177, +0.436]). If that number goes in the abstract, bootstrap it.
- **A per-LFE amplitude comparison is still not possible** (PNSN energy = windowed envelope, and
  S14 showed `mf_*.csv` cannot give per-detection amplitudes).
- **Table 1 mixes two tremor selections** (`N` from S1's 40-km/event-duration definition; exponents
  from the S9/S11/S16 60-km/padded-window definition) — flagged in the caption.
- Only 4 of the 18 amplitude-era events are central (1361, 1124, 1127, 1434) and two events have
  just 8 fit days (1905, 1434), so the central and per-event amplitude numbers are weak.
- Figures are 7.4 in wide with 9 pt base fonts: legible at double-column, tight at single-column.
