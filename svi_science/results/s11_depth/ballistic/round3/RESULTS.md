# S11 ballistic timing, round 3: results (Mon 2026-10-12, 8:20 pm PT)

Rules: `PLAN.md` (written 7:30 pm PT before any round-3 number; Amendments 1–5a time-stamped, each before the number it
affects). Scripts `svi_science/scripts/s11_ball3_{measure,lib,clock,longterm,phase,common,jointp,figures,dh}.py`, queues
`s11_ball3_queue.sh` (headline) and `s11_ball3_queue_dh.sh` (Amendment 5). Tables: `tables/` (both regions),
`svi/tables/`, `wa_nw/tables/`. Figures: `figures/`. "Scrambled control" = c2 windows (each station's data shifted to
another time of the same day). dv/v and δV/V in %. The SVI part of the engine reproduces round 2 exactly (family-10 stacks
identical; B926-leg decade −0.13 ± 1.29 ms as in round 2).

## What the data show

**No arrival-time change is detected anywhere: not over the decade, not through the ETS cycle, not per episode, in
either region. The bounds are now 0.01–0.04 % and are replicated in Washington. Three instrument or catalogue effects
were found and identified, and none of them is a path change.**

| task | pre-registered rule | outcome |
|---|---|---|
| 1 WA-NW replication | co-located controls within 3σ AND a far leg (≥ 20 series) with 2σ ≤ 0.03 % against both references | **year bins FAIL** (WA has too few repeats per year: ≤ 25 series per far station). **Decade-half adaptation (Amendment 5): the leg part PASSES**: B001 leg < 0.018 % (G1, 37 series) and < 0.024 % (G2, 44 series), the same level as SVI's B926 leg (< 0.014 %). The co-located part FAILS because of a B013 instrument step (below) |
| 2 ETS-phase composite | 3σ + label perm + date perm + control + same sign in both regions | **no detection, no candidate**. Front arrival vs between-ETS, S−P: SVI +0.008 ± 0.016 %, WA +0.013 ± 0.018 % (constant Vp/Vs) |
| 3 region-wide change | 2σ bound Vs-only < 0.03 % or constant-Vp/Vs < 0.05 %, or a detection | **PASS (bounds)** in both regions. Joint fit: δVs/Vs < 0.02 %, δVp/Vp < 0.04 %, δ(Vp/Vs) < 0.025 % (2σ, decade). SVI shows a ~1.3 ms catalogue timing shift that S−P removes |
| 4 joint P alignment | σ_P(joint) ≤ 0.7 × best instrument at ≥ 2 of 3 sites | **FAIL** (only G1, 0.69–0.70). But the joint S−P cancels the instrument-specific P effects at the decade scale |

## 1. WA-NW replication of the decade bound
`wa_nw/tables/pairs_decade.csv` (year bins, every pair), `tables/wa_pairs_decade.csv`, `tables/dh_pairs.csv` (decade
halves, both regions), `wa_nw/tables/rule1.json`, `tables/rule1_dh.json`; figures `fig_wa_pairs_decade.png`,
`fig_wa_longterm_years.png`, `fig_dh_decade.png`.

**Year bins (headline, as in round 2).** WA families have ~10 repeats per family per major (SVI ~28). Few family-years reach
20 rows: at most 33 family series per far station have ≥ 3 year stacks, and ≤ 25 have ≥ 2 years on each side of 2015
(SVI 43–88). Best legs: B001 2σ < 0.028 % (G1
reference, 19 series), < 0.031 % (G2, 25 series); GNW < 0.030 %; B943/B013 vs G1 < 0.036 / 0.049 %. Rule 1(b) fails on the
series count.

**Decade halves (Amendment 5: one stack ≤ 2014 and one ≥ 2015 per family; contrast = half the lag difference).**
Per-station 2σ bounds on δVs/Vs along the station's own leg:

| leg | reference G1 (B005/6/7 mean) | reference G2 (B013/B943 mean) | series |
|---|---|---|---|
| B001 | **< 0.018 %** (−0.7 ± 1.2 ms) | < 0.024 % (+3.5 ± 0.7 ms, flip p 0.0001) | 37 / 44 |
| GNW | < 0.030 % | < 0.016 % | 25 / 32 |
| B943 | < 0.025 % | — | 37 |
| B013 | < 0.031 % | — | 37 |
| B005 / B006 / B007 | — | < 0.027 / 0.023 / 0.029 % | 32–35 |
| VGZ | < 0.056 % | < 0.042 % | 13 |
| LRIV (co-located with G1) | < 0.010 % | < 0.045 % | 25 / 26 |
| SVI, for comparison: B926 / SNB / VGZ (borehole mean) | < 0.014 / 0.019 / 0.028 % | | 72 / 38 / 37 |

- The B001 leg passes the replication thresholds with both references. Washington gives the same ~0.015–0.02 % level as SVI.
- **No detection.** G2−B001 +3.5 ± 0.7 ms is not seen against G1 (−0.7 ± 1.2 ms), and the loop G1→B001→G2 does not close
  (1.5σ misfit). It is a G2/B013 instrument term, not a path change.
- **Co-located controls.**
  - G1 (B005, B006, B007, LRIV): all six pairs within 3σ in decade halves (−1.8 to +0.7 ms). G1 is clean.
  - **B013 has an instrument step:** its horizontal-channel S moved −3.9 ± 0.4 ms after 2015 against the detector frame,
    while co-located B943 (40 m away) moved +0.9 ± 0.4 ms. B013's P did not move (+0.6 ± 1.3 ms). So B013 − B943 S changed by
    ~5 ms (decade halves +5.2 ± 0.3 ms; year bins +1.4 ± 0.3 ms) and B013 S−P shortened by 4.5–7.7 ms (p ≤ 0.009).
    B943 S−P is unchanged (+0.8 ± 1.4 ms). It cannot be a path change: two instruments 40 m apart disagree, and only the
    horizontals moved.
  - SVI: B010 against B009/B011 −1.1 / +1.6 ms over the decade (p 0.0001), with S−P at B010 −2.4 and B011 +2.4 ms in
    opposite directions (round-2 finding, now in decade halves).

## 2. ETS-phase composite (both regions)
`tables/phase_contrasts.csv`, `phase_calls.csv`, `phase_curves.csv`, `phase_perm_dists.npz`; figure `fig_phase_composite.png`.
Phases: between-ETS = burst rows; arrival = ≤ 1 day after the family's onset (10 % quantile, the S3 definition); behind the
front = > 1 day. Pooled over all majors, every family and station, 200 episode-date permutations.

| region-wide item | C1 arrival − between | C2 behind − arrival | C3 behind − between |
|---|---|---|---|
| SVI S−P, constant Vp/Vs | +0.008 ± 0.016 % (2σ < 0.039) | −0.047 ± 0.023 (label p 0.010, date p 0.07) | −0.016 ± 0.033 |
| WA S−P, constant Vp/Vs | +0.013 ± 0.018 % (2σ < 0.049) | +0.048 ± 0.017 (label p 0.023, date p 0.19) | +0.068 ± 0.024 (label p 0.022, date p 0.020) |
| SVI S−P, Vs only | +0.003 ± 0.007 % (< 0.017) | −0.020 ± 0.010 | −0.007 ± 0.014 |
| WA S−P, Vs only | +0.005 ± 0.008 % (< 0.021) | +0.021 ± 0.007 | +0.029 ± 0.010 |
| SVI dS far legs (B926/SNB/VGZ) | +0.012 ± 0.007 % (< 0.027) | +0.004 ± 0.006 | +0.008 ± 0.010 |
| WA dS, G2 reference | −0.003 ± 0.007 % (< 0.017) | +0.016 ± 0.005 (label p 0.002, date p 0.07) | +0.016 ± 0.008 |

- **Arrival times do not change systematically through the ETS cycle.** The arrival-vs-between contrast is zero in both
  regions to ±0.02 % (S−P) and ±0.007 % (S legs).
- The behind-vs-arrival S−P contrasts have opposite signs in SVI and WA. None passes the episode-date permutation at 0.01.
- WA G1-referenced dS differs between majors and bursts (C3 +12.8 ± 2.8 ms, label p 0.0003). The date permutation gives the
  same contrast with fake front dates (p 0.30), so it is a bursts-vs-majors difference (which part of a family's patch is
  active), not tied to the front.
- The scrambled control shows nothing (label p 0.06–1.0).
- Finer phases (60-day pre/post-major bursts, first 10 %, 1–3 d, > 3 d) and thirds of the episode window are columns
  (`phase_curves.csv`). All lie within ±0.03 % (S−P, constant Vp/Vs).

## 3. Region-wide (common-mode) change
`tables/common_decade.csv`, `common_year_curves.csv`, `common_joint_ab.csv`, `common_majors*.csv`, `rule3.json`,
`dh_region.csv`, `dh_joint_ab.csv`; figures `fig_common_years.png`, `fig_dh_decade.png`.

**S−P over 20 years, all stations and families pooled** (error = the larger of the formal error and the station-to-station
scatter):

| | year bins: constant Vp/Vs | Vs only | decade halves: constant Vp/Vs | Vs only | Vp only |
|---|---|---|---|---|---|
| SVI (6 stations, 2003–2026) | −0.005 ± 0.021 % (< 0.046) | < 0.020 % | +0.001 ± 0.017 % (< 0.034) | < 0.014 % | < 0.025 % |
| WA (15–22 stations, 2007–2026) | +0.027 ± 0.015 % (p 0.045; B013) | < 0.025 % | +0.009 ± 0.009 % (< 0.027) | < 0.012 % | < 0.020 % |
| WA without B013 | | | −0.008 ± 0.011 % (< 0.029) | < 0.013 % | < 0.022 % |

**Joint uniform δVs/Vs (a) and δVp/Vp (b)** from S−P plus the moveout of dS and dP, decade halves:
- SVI: a = +0.011 ± 0.005 %, b = +0.017 ± 0.011 % (sign-flip p 0.20 / 0.25). 2σ: a < 0.021 %, b < 0.039 %,
  δ(Vp/Vs) < 0.022 %.
- WA, G1 reference, B013 excluded: a = +0.001 ± 0.009 %, b = +0.007 ± 0.017 %. 2σ: a < 0.019 %, b < 0.040 %,
  δ(Vp/Vs) < 0.025 %.
- The G2-referenced WA fit gives a = −0.03 % (p 0.001). That comes from the G2/B013 instrument terms of §1, not from the
  region.

**Catalogue timing drift vs a real change (SVI).**
- In decade halves, the absolute S lag (against the catalogue time) moved +1.10 ± 0.37 ms (flip p 0.002) and the absolute
  P lag +1.52 ± 0.54 ms (p 0.03).
- S−P moved −0.01 ± 0.87 ms.
- S and P moving together with S−P flat is the catalogue (detector-frame) signature: a ~1.3 ms shift of the catalogue
  times across the decade, probably the 2014 switch from the backward to the forward catalogue.
- The pre-registered year-bin call needs 3σ and got 2.7σ (absS +1.45 ± 0.54 ms), so formally no drift is called.
- S−P is immune to it by construction. WA shows no catalogue drift (absS −0.6 ± 0.5, absP +0.8 ± 0.7 ms).

**Per ETS (majors, demeaned over majors):**
- σ per major is 0.036 % (SVI) and 0.046 % (WA), constant Vp/Vs. Majors scatter more than their errors (χ²/dof 2.8 and 2.7;
  scrambled control 0.8 and 1.6).
- **Station-split coherence** (Amendment 3; SVI computed after its χ² was seen): SVI r = 0.48 (p 0.005), WA r = 0.26
  (p 0.07). Restricted to well-measured majors (σ ≤ 0.05 %), SVI −0.23 and WA −0.32.
- **Family-split coherence** (Amendment 4): SVI r = 0.03 (p 0.44), WA r = 0.31 (p 0.06).
- The excess sits in small majors with 1–4 active families. There it is a constant S−P offset at every station (e.g. SVI
  2020-05: −26 ± 7 ms offset, no growth with S−P; 2018-05 +29 ± 9 ms). That pattern points to the source (a depth or locus
  shift of the active patch), not to velocity. **No regional per-ETS S−P change.**

## 4. Joint multi-instrument P alignment
`tables/jointp_quality.csv`, `jointp_decade.csv`, `dh_sitejoint.csv`, `rule4.json`; figure `fig_jointp.png`.
- **Joint σ_P / best single instrument (same rows):**
  - SVI PGC site 0.80–0.97;
  - WA G1 0.69–0.82 (passes on EP and YR bins);
  - WA G2 0.93–1.01.
  - Rule 4 fails (1 of 3 sites). The scrambled control's joint P cc is 0.12–0.23.
- **Why it does not shrink.** The bootstrap σ resamples the LFE rows, and that variability is shared by co-located
  instruments, because the same events (and their source differences) are in every stack.
- **What it does buy.** At the decade scale the joint S−P cancels instrument-specific P/S effects:
  - SVI PGC site: joint +0.07 ± 0.84 ms, while B011 gives +2.5 ± 1.0 and B010 −2.4 ± 1.0 (p ≤ 0.003);
  - WA G2: joint −1.8 ± 1.4 ms, while B013 gives −4.6 ± 1.4 and B943 +1.1 ± 1.5;
  - the robust year-to-year scatter of S−P drops from 20–26 ms to 14 ms (SVI) and from 45–61 ms to 25 ms (WA G1).

## What failed and the next adaptation
1. **WA year bins (rule 1).** Too few repeats per family-year. The decade-half stacks fixed the leg part. Next:
   - the B013 horizontal step needs a date. Q4 bins put it between 2011–14 and 2015–18. Check the B013 channel metadata and
     orientation history (EarthScope) and look at the YR S−P;
   - then use B943 alone as the G2 reference.
2. **Co-located drifts (B013, B010).** Model them as instrument steps with dates, in the same way as the round-2 episode
   glitches, before any joint inversion.
3. **ETS phase.** Bound ~0.02 % per cycle at constant Vp/Vs. Next: pool both regions with a family-level depth/locus offset
   term (the per-ETS excess is source-side), and add a pre-major-burst vs post-major-burst contrast with the bursts stacked
   over several years (BPRE/BPOST are thin per family).
4. **Joint P (rule 4).** More instruments cannot beat the shared source-side σ; more families can. Use the site-joint S−P as
   the instrument-robust column in any decade inversion.
5. **Catalogue drift.** Date the ~1.3 ms SVI shift with the absolute S and P year series (expected 2014-01, the catalogue
   switch). S−P and inter-station differences are unaffected.
