# S11 ballistic timing, round 2: results (Mon 2026-10-12, 4:02 pm PT)

Rules: `PLAN.md` (written 3:35 pm PT before any round-2 number; Amendments 1–5 time-stamped, each before the number it
affects). Scripts `svi_science/scripts/s11_ball2_{measure,clock,meta,ppass,longterm,locus}.py`, queue
`s11_ball2_queue.sh`. The measurement pass reproduces round 1 exactly (co-located floor 10.18 / 10.02 ms SD, 4.10 / 3.53 ms
robust, as in `../RESULTS.md`). "Scrambled control" = S10 c2 windows (stations shifted to another time of the same day).

## What the data show

**One adaptation works: pooling over years and families gives a new, tight long-term constraint.** Across 2015, the S-time
difference between the PGC-site boreholes and B926 changed by **−0.1 ± 1.3 ms**. On B926's ~19 s S leg that is
**δVs/Vs = 0.00 ± 0.007 %, 2σ bound 0.014 %** (0.024 % taking the least favourable reference). That is ~10× tighter than any
coda-based decade bound from S11 (Q2 < 0.35 %). No long-term change is detected anywhere. The per-episode floor still does
not reach 3 ms, and the locus wander is not detected by the pre-registered rule.

| adaptation | rule | outcome |
|---|---|---|
| 1 clock correction | robust floor ≤ 3 ms after correction (both reference sets) AND far scatter reduced | **partial**: floor 4.8 / 3.8 ms (FAIL); far scatter reduced for 2 of 3 far stations with each reference (PASS). Best S variant (3C, [tS−0.5, tS+4] s): 3.3 / 3.2 ms, still a narrow FAIL |
| 2 better P | P cc ≥ 0.6 for ≥ 3 stations; control < 0.3 | **PASS, marginal**: P3 (3C principal direction) 0.61 / 0.62 / 0.60 at PGC / B010 / B011 on 2-yr stacks; control ≤ 0.27. But σ_P is still 13–20 ms per stack |
| 3 pooled long-term | a 2σ bound < 0.1 % or a detection | **SUCCESS (bounds)**: δVs/Vs on the B926 leg < 0.014–0.024 %; SNB/VGZ legs < 0.04–0.06 %; S−P < 0.075–0.11 % (constant Vp/Vs), < 0.03–0.05 % (Vs only); δVp/Vp on the B926 and VGZ legs < 0.07 %; Vp/Vs on the B926 leg < 0.07 %. No detection |
| 4 locus nuisance | χ² p < 0.01 AND geometry-scramble p < 0.01 | **not detected**: χ² p ≈ 0, but the geometry-scramble control gives p = 0.09 (fits have χ²/dof 2.4). Post hoc, family 87 alone: ~40 m RMS wander (p 0.003) |

## 1. Clock and timing-glitch correction
- **Model.** The 6 co-located pair differences per (family, episode) are solved for 4 instrument lags. The zero is the
  median instrument ("majority clock"). Averaged over the other families (leave-family-out), this is the instrument term.
  Tables: `clock_terms.csv`, `clock_floor_far.csv`, `far_terms.csv`, `frame_drift.csv`. Figures:
  `fig_clock_terms.png`, `fig_floor_before_after.png`.
- **What the correction removes.** The episode-common part of the co-located differences:
  - PGC set: SD 2.6 → 0.0 ms (p 0.0001 → 0.26);
  - B011 set: 2.8 → 1.6 ms (p still 0.0001);
  - the B010 2024-08 step is flagged (−9.3 ± 2.5 ms, 91 % of 35 families agree).
- **Why the floor barely moves.** The robust floor goes 4.7 → 4.8 ms (PGC set) and 4.1 → 3.8 ms (B011 set) on the correctable values.
  - The floor is per-family S-measurement noise (median bootstrap σ 2.4–3.7 ms per value), not clock steps.
  - It is worse for PGC E3 (2017+) and for B009–B010 (robust 4.3–5.5 ms) than for PGC E1 and B009–B011 (2.6–3.4 ms).
- **Round-1 glitches.** PGC 2019-11 (−6.0 ms) and 2021-01 (+3.3) are recovered as PGC instrument terms. They are not flagged at 3σ because the
  family scatter is heavy-tailed (SD 5–11 ms). The B010 2009-03 term is +11 ± 15 ms from 6 families (unresolved).
- **S-estimator variants** (Amendments 4–5; quartet pairs, episode stacks; robust floor after correction, PGC set / B011 set):

  | S window | robust floor (ms) | SD (ms) |
  |---|---|---|
  | round-1 [tS − 0.5, tS + 1.5] s, horizontals | 4.8 / 3.8 | 10.2 / 10.0 |
  | direct S only [tS − 0.3, tS + 0.7] s | 5.5 / 4.5 | — |
  | 3C [tS − 0.5, tS + 2.5] s | 3.6 / 3.1 | 6.5 / 6.1 |
  | 3C [tS − 0.5, tS + 4.0] s | 3.3 / 3.2 | 6.4 / 5.9 (all values) |

  More S energy halves the plain SD and lowers the robust floor. It plateaus just above 3 ms.
- **Far pairs.** Correction lowers the scatter only slightly, e.g. PGC–B926 19.8 → 19.7 ms and PGC–SNB 20.6 → 20.4 ms. The scatter is dominated by the far station's own noise (σ 9–14 ms).
- **Far-station terms** (far station minus PGC site; clock and near-site path cannot be separated): 5 flags, all in few-family
  bursts or low-cc episodes, except VGZ 2020-11 (+7.3 ± 2.0 ms, 18 families). Not interpreted.
- **Metadata** (`tables/meta/`). EarthScope MUSTANG `timing_quality`, `clock_locked`, `suspect_time_tag` and `timing_correction` are empty or zero for these
  stations, so there is no metadata corroboration. B010 had data gaps on days in each of its 2009-03, 2022-03 and 2024-08
  episodes, but gaps occur on 13 % of all B010 days, so this is not specific. CN channel bands (BH vs HH) never mix inside any episode (cut-store index).

## 2. Better P
`P_quality.csv`, `P_splithalf.csv`, `P_verdict.json`, `fig_P_cc.png`. Median P cc per station on episode (EP), 1-year (YR) and 2-year (Y2) stacks:

| estimator | PGC | B009 | B010 | B011 | B926 | VGZ | SNB | scrambled control |
|---|---|---|---|---|---|---|---|---|
| P3 (3C principal direction, headline) | 0.58–0.61 | 0.52–0.57 | 0.58–0.62 | 0.56–0.60 | 0.37–0.39 | 0.41–0.43 | 0.48–0.50 | 0.20–0.27 |
| Pzr (round-1 estimator) | 0.43–0.46 | 0.41–0.43 | 0.44–0.48 | 0.41–0.44 | 0.27–0.28 | 0.29–0.30 | 0.37–0.39 | ≤ 0.19 |

- Pz is close to P3 except at B009, whose P arrives 83° from vertical: its Z carries little P.
- Joint 3C cc (P3j) is worst (0.22–0.39).
- Rule 2 passes on P3 for 3 stations (PGC, B010, B011, Y2 stacks), marginally.
- P precision is still poor: median σ_P 13–20 ms at the PGC site and 20–37 ms at the far stations, against σ_S 4–10 ms.
- Year stacks are not much larger than episode stacks (about one ETS per year), so cc barely rises from EP to YR.

## 3. Pooled long-term constraints (year bins headline; decade contrast = years ≥ 2015 minus ≤ 2014)
`longterm.csv` (all items, YR and Y2, data and control), `longterm_calls.csv`, `longterm_curves.csv`, `fig_longterm_years.png`.

**Method check (co-located, corrected).** B011 vs B009/B010 decade contrast 0.34 ± 0.19 ms, p 0.08. Timing between co-located instruments is stable to ~0.2 ms over a decade once ~50 family series are pooled. The σ of each contrast matches its permutation SD (e.g. B926: 1.29 vs 1.28 ms), so the errors are honest.

**Inter-station S, decade contrast** (borehole-site mean reference, clock-corrected):

| far leg | Δ (ms) | perm p | series | δVs/Vs on far leg, 2σ | uniform on both legs, 2σ | other references |
|---|---|---|---|---|---|---|
| B926 | −0.1 ± 1.3 | 0.91 | 43 | **< 0.014 %** | < 0.07 % | B011 −2.1 ± 1.1 (p 0.10); B009+B010 +0.7 ± 1.4; site mean −0.2 ± 1.0 → conservative < 0.024 % (leg), < 0.11 % (uniform) |
| SNB (E4 only) | +3.9 ± 2.6 | 0.28 | 8 | < 0.056 % | < 0.23 % | — |
| VGZ (E3 only) | −2.0 ± 1.9 | 0.20 | 7 | < 0.040 % | ΔT too small | — |

PGC-referenced decade contrasts are not usable: PGC's epoch changes in 2015 and 2017 zero them by construction. Rows `dS PGC-*` with contrast 0.00 are that artefact.

**S−P, decade contrast** (P3):

| station | Δ (ms) | constant Vp/Vs, 2σ | Vs only, 2σ |
|---|---|---|---|
| B009 | −0.6 ± 2.2 | < 0.075 % | < 0.032 % |
| B010 | −2.3 ± 1.4 | < 0.082 % | < 0.035 % |
| B011 | +3.2 ± 1.2 (p 0.014) | < 0.086 % | < 0.037 % |
| B926 | −2.2 ± 3.2 | < 0.105 % | < 0.045 % |
| VGZ (E3) | +0.3 ± 3.3 | < 0.11 % | — |
| SNB (E4) | −3.1 ± 7.2 | < 0.21 % | — |

- B010 and B011 sit 4 m apart and move in **opposite** directions (Y2: −3.0 vs +3.4 ms, each p 0.01). That is an instrument-level
  P-timing effect, not a path. It sets the practical S−P floor at ~3 ms per decade.

**P between stations.**
- δVp/Vp on the B926 leg: 2σ < 0.067 % (−1.1 ± 3.1 ms on T_P ≈ 10.9 s).
- VGZ leg: < 0.07 % (0.11 % with the permutation SD). SNB leg: < 0.19 %.
- **Vp/Vs on the B926 leg** (from dP and dS): +0.009 ± 0.029 %, 2σ < 0.07 %.

**Year-to-year and trends.**
- No far-pair or S−P family-common year series has permutation p < 0.01.
- The VGZ legs show p 0.01–0.05: a VGZ-site year term of SD ~2.5–3.4 ms, with B009+B010 and with every other reference.
- Co-located sets still carry a year term (SD 1.1–1.7 ms, p ≤ 0.005): residual instrument timing at year scale.
- Linear trends are all within 2σ of zero, e.g. B926 −0.13 ± 0.11 ms/yr.

**Scrambled control.** Few control stacks survive (n matched, ≥ 20): 10–75 values per item. No contrast has p < 0.01. Control scatter is 15–90 ms.

**Detection rule:** nothing qualifies. B011 S−P has p 0.014 > 0.01 and is contradicted by B010.

## 4. Source-locus nuisance (2003–06, families 10, 87, 17; POLARIS + core; Amendment 3)
`locus_lags.csv`, `locus_fits.csv`, `locus_summary.csv`, `fig_locus.png`.
- 18 family-episodes with 7–14 stations at 7–14 sites.
- Horizontal shifts relative to each family's mean:
  - family 87: formal σ 18–48 m;
  - family 10: 54–460 m;
  - family 17: 130–900 m (one burst at cc 0.29 has 3–4 km nonsense shifts).
- Pooled χ² of horizontal Δx is 187 for 36 dof (p ≈ 0). The geometry-scramble control (station geometries permuted) gives median 134 and **p 0.09**, so wander is **not detected** by the rule. The large χ² is mostly underestimated noise: fit χ²/dof 2.4.
- Bound on RMS centroid motion, all three families: < ~1 km. This is dominated by family 17 and is not informative.
- Post hoc, labelled: family 87 alone, with errors inflated by each fit's χ²/dof, gives χ² 32.8 / 14 dof (p 0.003). RMS shift 68 m against 54 m formal → **~40 m of wander**. Family 10: none (p 0.6). Family 17: p 0.26.
- Depth trades completely with origin time (corr −1.00) and is unresolved (σ_Z 300–3,000 m).
- The common velocity term also trades with origin time (corr 0.96–1.0): c = +0.13 % with locus terms vs −0.04 % without; per-fit σ 0.3–2 %. The locus nuisance does absorb velocity information, as Merlin warned.

## What did not work, and what to try next (round 3)
1. **Per-episode floor ≤ 3 ms.** Longer 3C S windows came to within 0.1–0.3 ms. Next:
   - inverse-variance weighting of the floor (the bootstrap σ tracks the scatter; round 1 found ×1.3–2);
   - per-instrument response check of PGC E3 and B009–B010, the worst pairs.
2. **P precision.** σ_P is 13–20 ms even when cc ≥ 0.6. Next:
   - multi-station P alignment (P lags of the 4 co-located instruments solved jointly, as for S);
   - a P window chosen from the split-half cc;
   - the instrument-level P effect at B010/B011 needs a separate look.
3. **Episode-scale path signal.** Still bounded by far-station S noise (σ 9–14 ms). The year/decade pooling is what works. Next:
   - extend it to WA-NW (LRIV + B005/6/7 quartet with B013/B943), as a replication of the 0.014 % B926-leg bound;
   - an ETS-phase composite (inside vs outside ETS) with the same machinery.
4. **Locus wander.** Next:
   - weight fits by their own χ²/dof;
   - drop the velocity term (it trades with origin time);
   - pre-register family 87, the best-constrained family, as a test on POLARIS-era bursts.

## Files
- `tables/lags_station.parquet`: per (family, station, epoch, bin ∈ EP/YR/Y2, label, data/ctrl): S and 4 P estimators, cc, σ, rails, S−P σ.
- `tables/lags_pair.parquet`: 18 common-row pairs × 3 bin types, S and P3.
- `tables/series_geom.csv`: picks, P3 direction, incidence, split-half cc.
- Variants: `*_s3c25`, `*_sh07`, `*_s3c40` (Amendments 4–5).
- `*_test` = family-10 test (Amendment 1).
