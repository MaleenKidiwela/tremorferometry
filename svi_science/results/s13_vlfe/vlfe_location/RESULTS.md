# vlfe_location — RESULTS (written 2026-10-16 ~10:20 am PT)

Rules: `PLAN.md` (~10:02 am PT, before any number); changes: `NOTES.md` amendment A1 (~10:12 am PT, before any Test-4
number) and one post-hoc diagnostic (Test 5, disclosed as such). Event set: the catalogue's 242 ETS-period tier events
2010–2018 (`../catalogue/tables/prod<Y>_events_v2.csv` 2010–15, `prodr<Y>_events_v2.csv` 2016–18). 188 of them have at
least one located P ≥ 0.75 LFE within ±10 min (SVI families + WA-NW `selected_v2`, at their `families_v1_3.csv`
locations). "Offset" = VLFE catalogue location minus the LFE-weighted centroid of the families that fired in ±10 min.

Raw figures first: `figures/t1_map_offsets.png` (VLFE → LFE-centroid lines coloured by year; offset vectors),
`figures/t1_offset_hist.png` (east/north/distance, real vs scrambled control), `figures/t2_zsurface_examples.png`
(production z surfaces), `figures/t4_gf_variants.png` (distance per Green's function). No error bars.

## Plain answer
- **The user is right that the locations are off, but they are not random.**
  - The VLFEs sit a median **16 km** from the concurrent LFE centroid (IQR 8.5–32 km; mean 23 km, pulled by a tail of
    events placed on the far west of the grid).
  - The offset has a direction: the mean is **11.5 km toward the NW** (az 311°; 8.7 km W, 7.5 km N). Medians: 7.1 km
    W, 4.5 km N.
  - Spread: SD 21 km east–west, 20 km north–south.
  - **Systematic bias: yes.** Bootstrap over episodes, 95 %: east −15.7 to −4.6 km, north +5.0 to +11.1 km; both
    exclude 0.
  - **Tracking: yes.** With the LFE centroid of another event of the same episode (scrambled control), the median
    distance is 26 km (5th percentile of 1000 re-pairings 24 km) against 16 km real (p < 0.001). So each VLFE does sit
    nearer its own LFEs than other events' LFEs. Note: a re-pairing keeps the mean offset by construction, so it tests
    tracking only. The bias is tested by the bootstrap.
- **The grid cannot do better than ~15–20 km.**
  - The Δz ≤ 1 region around the best node has a median of 8–9 nodes, i.e. an equivalent radius of 16–17 km. It spans
    20 km E–W and 20 km N–S, and its farthest node is 22 km from the best one.
  - The LFE centroid's node lies inside the Δz ≤ 1 region for 49 % of events (Δz ≤ 2: 69 %). So the offsets are about
    the size of what the data can resolve.
- **Not simple template copying.**
  - 32 % of events locate within one grid step (≤ 10 km) of the detecting template's own location, against 11 % for a
    random template of the same bank (28 % without the 17 relay self-detections).
  - The best node is nearer the template than the LFE centroid in only 52 % of events. The templates themselves sit
    17 km (median) from these LFE centroids, about the same as the VLFEs.
  - The scan that locates the event uses synthetic grid templates, not the matched-filter template. The excess over
    random is what one expects when similar waveforms carry similar amplitude patterns.
- **A better Green's function does not shrink the offsets** (pre-set bar: median drops ≥ 3 km and ≥ 60 % of events
  closer).
  - 2011–13 events (n = 105): production median 15.1 km; half-space Okada 15.0; Okada + stack terms 14.1; whole-space +
    stack terms 12.9 km (−2.2 km, 33 % of events closer, 20 % farther). No variant passes the bar.
  - The NW bias stays in every variant (mean 6.9–11.1 km W, 3.9–9.7 km N).
  - The 1-D layered GF (Crosson / CVM v1.7) was not run: no layered-GF code or store is on disk offline.
- **What limits the location is the network, and the LFE "truth" itself is broad.**
  - Only ≤ 2 live channels west of the patches in every event (PFB Z, YOUB Z; NLLB and LZB off, horizontals off by tilt
    noise). Everything else is east or south-east.
  - The scan uses 9–18 channels, mostly verticals. The far-west outliers (124.3–124.5 °W, 77–135 km off) are where
    nothing constrains the western side.
  - The LFEs firing in the same ±10 min are themselves spread with an RMS of 11 km about their centroid (18 km when
    ≥ 5 LFEs). So ~10–15 km is about the best agreement one can expect against this reference.

## Test 1 — offsets to concurrent LFEs (`tables/t1_offsets.csv`, `t1_summary_year.csv`, `t1_summary.json`, `t1_bootstrap_by_lfeset.txt`)
| year | n | mean E (km) | mean N (km) | mean vector | SD E | SD N | median distance (IQR) |
|---|---|---|---|---|---|---|---|
| 2010 | 8 | −22.2 | +20.8 | 30.4 km, 313° | 29.7 | 21.3 | 26.4 (15.7–53.4) |
| 2011 | 37 | −11.7 | +7.1 | 13.6 km, 301° | 10.9 | 9.0 | 15.1 (8.4–23.8) |
| 2012 | 51 | −8.3 | +4.8 | 9.6 km, 300° | 16.5 | 16.2 | 14.7 (7.3–21.3) |
| 2013 | 17 | −18.4 | +8.3 | 20.2 km, 294° | 28.1 | 23.4 | 21.5 (7.9–33.1) |
| 2015 | 3 | −0.7 | −18.4 | 18.4 km, 182° | 19.6 | 0.6 | 23.4 |
| 2016 | 7 | −15.9 | +0.2 | 15.9 km, 271° | 7.6 | 16.9 | 17.8 (12.9–28.3) |
| 2017 | 7 | −9.9 | +21.2 | 23.4 km, 335° | 43.1 | 35.8 | 11.7 (7.8–39.2) |
| 2018 | 58 | −2.0 | +8.7 | 8.9 km, 347° | 21.7 | 22.8 | 20.0 (11.6–43.9) |
| **all ±10 min** | **188** | **−8.7** | **+7.5** | **11.5 km, 311°** | **21.0** | **19.7** | **16.2 (8.5–32.2)** |
| all ±2 min | 132 | −8.6 | +5.6 | 10.3 km, 303° | 22.1 | 20.4 | 14.8 (8.2–25.7) |
| SVI families only ±10 min | 180 | −4.7 | +2.6 | 5.4 km, 299° | 18.3 | 17.0 | 14.3 (7.9–23.2) |

- Bias, bootstrap over the 8 episodes (95 %): all ±10 min E [−15.7, −4.6], N [+5.0, +11.1]; ±2 min E [−17.2, −2.8],
  N [+3.1, +11.1]; SVI-only E [−9.8, −0.0], N [+0.6, +4.8] (just excludes 0).
  - The size of the bias depends on the LFE reference. The WA-NW families pull the centroid south-east, so the bias
    halves against SVI families alone. The direction (NW) is the same.
- Scrambled control (within-episode re-pairing, 7 episodes with ≥ 2 events, 187 events): median distance 26.4 km (mean
  of 1000 draws), against 16.0 km real; p < 0.001. Per-component SD: real 19.9 / 18.8 km, scrambled 24.7 / 24.7 km.
- Per episode: 2011-08 15 km, 2012-08 15 km, 2013-09 22 km, 2018-06 20 km (median). The single 2017-02 event is
  135 km off: placed at 124.26 °W against LFEs near 122.86 °W.

## Test 2 — achieved resolution (`tables/t2_zsurface_all.csv`, `t2_summary.json`)
- Code check: the production z surface is reproduced exactly (best node 312/312, z to 1e-3).
- 242 ETS tier events, Δz ≤ 1 region: median 8 nodes (800 km²), equivalent radius 16 km, farthest node 22 km, E–W and
  N–S extents 20 km. Only 1 event has a single node. Δz ≤ 2: 17 nodes, radius 23 km, farthest node 32 km.
- By scan z (events with LFEs): z 5–7 (n = 96) radius 19 km, offset 19 km; z 7–10 (49) 16 km / 15 km; z ≥ 10 (43)
  13 km / 15 km. The surface sharpens with z, but the offset to the LFEs does not shrink below ~15 km.
- LFE centroid's node inside Δz ≤ 1: 49 %; inside Δz ≤ 2: 69 %; median Δz there 1.0.
- Some surfaces have two lobes (e.g. 2018-07-05 11:40 PDT: one lobe at 48.7 °N −123.2, a second at 48.0 °N −124.5).

## Test 3 — template dependence (`tables/t3_template_copy.csv`, `t3_summary.json`)
- Best node within 10 km of the detecting template's node: 32 % (≤ 14.2 km incl. diagonal: 45 %). Either agreeing
  template: 33 %. Scrambled control (a random template of the same bank): 11 % (either: 21 %).
- Without the 17 relay self-detections: 28 % vs 12 %.
- 2011 bank (2010–15): 26 %; 2011 + 2018 relay bank (2016–18): 40 %. The relay templates are 2018 events located by
  the same scan.
- Median distance to the detecting template's node is 20 km; to a random bank template 35 km.
- The best node is nearer the template than the LFE centroid in 52 % of events. Template locations are 17 km (median)
  from the LFE centroids, and the VLFEs 16 km.

## Test 4 — better Green's functions (`tables/t4_reloc.csv`, `t4_summary.csv`, `t4_station_terms.csv`)
Same data windows, weights, live channels, grid and Slab2 depths as production. Statistic for every variant: weighted
VR of a thrust template, free M0 ≥ 0, common lag ±10 s.

| variant | 2011–13 (n = 105) median / mean distance (km) | mean E, N (km) | closer / farther than (0) | vetted 14: median (km) | 2010–18 (n = 188) median (km) |
|---|---|---|---|---|---|
| production catalogue | 15.1 / 19.7 | −11.1, +6.2 | — | 12.5 | 16.2 |
| (0) whole-space, VR | 15.1 / 19.3 | −10.8, +5.7 | — | 12.6 | 15.6 |
| (c0) whole-space + stack terms | **12.9** / 19.8 | −8.6, +5.2 | 33 % / 20 % | **10.4** | 14.5 |
| (b) Okada half-space static × step | 15.0 / 21.1 | −9.5, +7.4 | 33 % / 34 % | 15.7 | 15.8 |
| (b1) Okada, g = 1 | 14.1 / 19.4 | −6.9, +7.5 | 42 % / 27 % | 14.2 | 15.1 |
| (c) Okada + stack terms | 14.1 / 19.9 | −9.4, +3.9 | 36 % / 28 % | 11.7 | 15.0 |
| whole-space static ×2 (diagnostic) | 15.0 / 21.6 | −8.1, +9.7 | 35 % / 37 % | 15.7 | 15.6 |

- (0) reproduces the production node for 71 % of the 2011–13 events (different statistic: VR with a lag search vs
  MAD-normalised z). Its offsets are the same as production.
- No variant passes the pre-set bar. The best, (c0), gains 2.2 km in the median (2.2 km on the vetted 14).
- The NW bias survives every GF.
- Okada vs flat ×2: for static fields they differ by only 10–15 % at these distances (NOTES A1). The static variants
  move individual events (only ~25 % keep the production node) without improving the average.
- Stack terms: the Canadian stations and LRIV are mostly within ×0.7–1.6 of the whole-space prediction (Z 0.95–1.26).
  Several WA stations are far off: WISH N ×14, NLWA ×3–5, D03D Z ×5.3, SP2 opposite sign. These come from only 2–4 WA
  stacks; used terms are clipped to ×0.25–4.
- (a) 1-D layered model: not run. pyrocko is installed, but there are no qseis/qssp binaries or GF stores and no CPS,
  EDGRN or axitra. Building one needs a download.

## Test 5 — post-hoc diagnostics (`tables/t5_diag.csv`, `t5_diag_summary.csv`, `t5_live_channels_by_year.csv`)
- Live scan channels by year: PGC N,Z (Z only from 2015); VGZ E,Z; SNB 3-comp 2012–17; LRIV 3-comp 2010–15, **off
  2016–18**; PFB Z; YOUB Z (to 2014); NLLB Z only 2017–18; LZB never; WA stations (GNW, DOSE, SP2, D03D, WISH, NLWA) Z
  only.
- West of the patches: ≤ 2 live channels in all 188 events. No event has 3 or more western channels, so a
  with/without-west comparison is not possible.
- LFE spread inside the ±10 min window: RMS 11 km about the centroid (18 km when ≥ 5 LFEs). Annual medians range from
  1 to 24 km.

## Scope
- The reference is the location of LFE families active in the same ±10 min. It is not the VLFE's own source, and it is
  itself 11–18 km broad.
- Depth is not free: grid nodes sit on Slab2. So nothing here tests VLFE depth.
- Test 4 is the cheap version: static GFs, a 10-km grid, the production channel set. A dynamic layered GF and a finer
  grid were not tried.
- What would help most is network geometry: horizontals cleaned of tilt, western stations (PFB/YOUB/NLLB 3-comp,
  POLARIS/CAFE eras), and LRIV back in 2016–18. That would be the next test, rather than the GF.
