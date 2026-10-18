# S12-LVL, LFE side — results in plain words

Written Tue 2026-10-13 ~10:30 am PT. Rules: `PLAN.md` (§1–5 and §6 written before their numbers; post-hoc items are
labelled). Decision log `NOTES.md`. dv/v in %; "scrambled control" = front dates moved to other episodes (A, task 3),
predicted lags permuted within station (task 2), or time-since-ETS shuffled within series (B). Everything is provisional.

## What the data show (short)
**Update Wed 2026-10-14 ~2:55 pm PT (§7b, §7c):** the 1–3 Hz restack finds no LVL-base reflection (scan p 0.19; it
would need ≥ 0.12–0.18 × the direct S, against ~0.13 for a sharp base). SVI shows no recharge slope in either era
(2003–2013 −0.004 ± 0.017–0.027 %/yr, 2014–2026 −0.005 ± 0.013–0.020 %/yr; target −0.4).
**Update ~10:45 am PT (§7): the WA candidate below did NOT confirm. Over all 50 WA-NW families the slope is
−0.055 ± 0.053 %/yr, and the 39 families outside the original 11 give +0.008 ± 0.059.**
1. **One candidate, in Washington: between ETS, the coda slows down as time since the last ETS grows.** Between-ETS
   LFE days (≥ 60 d after a front), dv/v 2–5 s: **−0.33 %/yr** (bootstrap ± 0.14, 2.3σ; scrambled control p 0.004).
   - Every leave-one-site-out fit stays negative (−0.41 to −0.26). Odd and even families agree (−0.67 / −0.35).
     A calendar trend term does not remove it (−0.37).
   - The binned curve falls from +0.10 % (2–4 months after ETS) to −0.04 % (9–12 months). Most of the fall is in the first
     ~6 months.
   - This is the sign and size the fault valve predicts (target −0.4 %/yr for a +2 % LVL change that decays over the cycle).
   - It fails the pre-set 3σ bootstrap rule (2.3σ), so it is a **candidate**, not a detection.
   - Southern Vancouver Island shows none of it: −0.006 ± 0.008 %/yr, 50× smaller, with 6× more days.
2. **No ETS-locked LVL change in the coda.** The forward model says a Gosselin-sized change (+2 % Vs in a 3.5 km LVL
   below the LFEs) would read **+0.35–0.45 %** in coda dv/v. The contrast between ETS days 0–20 and between-ETS days is:
   - SVI: +0.006 ± 0.011 % (2–5 s; every lapse window within ±0.03 %), which bounds the LVL change to < 0.085 % (2σ),
     about 25× below the target;
   - restricted to Gosselin's footprint (LVL within 15 km of PGC/SNB/VGZ): −0.01 ± 0.10 %, < 0.2 %, 10× below;
   - WA-NW: +0.021 ± 0.028 %.
3. **No LVL-base reflection at 2–8 Hz.** Injection shows the scan would find a reflected S of ≥ 0.06–0.08 × the direct-S
   amplitude. A sharp LVL base predicts ~0.13. So, beneath these families, the LVL base is gradational at 2–8 Hz, or its
   reflection is weaker than ~0.07.
4. **No propagating onset (A) and no step at the front, in either region.** The source scalars (C) change only on the
   busiest days. Episode-size scaling (D) is weakly positive in WA only.

## 1. Forward calibration (task 1) — `tables/fwd_*`, `figures/fig_fwd.png`
- **Kernel.** The S11 Monte-Carlo kernel, with an LVL tally added. The variant with CVM speed in the LVL reproduces
  S11's a_d (0.548 vs 0.547). The footprint run reproduces the nominal coefficients exactly.
- **a_LVL (2–5 s, 3C).** Median 0.18–0.22 over kernels (pairs q10–q90 0.07–0.43).
  - Per window, 1–2 / 2–3 / 3–4 / 4–6 s: 0.14–0.19 / 0.17–0.22 / 0.16–0.22 / 0.17–0.20.
  - The layer sits right where the S-anchored coda dwells (just below the source), and it is slow.
  - Placements: 2 km thick 0.16; source 1 km inside 0.28; top 1 km below the source 0.14. Footprint ≤ 15 km of the three
    stations: median 0.013, > 0.05 for 20 of 101 families.
- **Reflections (flat-layer rays).**
  - LVL-base SrS arrives 1.4–2.3 s after S (median 1.8 s). It spends 2.95 s in the LVL, so a +2 % change → **−59 ms**.
  - Moho SrS (7 km below the top) arrives ~3.0 s after S.
  - PrS arrives within ±1 s of S; SrP sits in the P coda; PrP is 0.05–1 s after P.
- **Source side.** For a source inside the LVL, P/S rises 5 % (Vp fixed) or 2.5 % (δVp = δVs/2). Direct S−P change ≈ 0.

## 2. LVL-base reflection search (task 2) — `tables/refl_*`, `figures/fig_refl_{records,scan}.png`
- **Setup.** 497 (family, station) pairs. Moveout scan over H 1.5–6 km × Vs 2.6–3.8 km/s, against 1000 scrambled draws
  (scan-corrected).
- **Scan results:**

  | phase | best model | peak | scan p | nominal-model p |
  |---|---|---|---|---|
  | SrS LVL base (H envelope) | H 6, Vs 3.6 | 0.088 | 0.62 | 0.49 |
  | SrS Moho | H 2, Vs 3.8 | 0.096 | 0.054 | 0.98 |
  | SrP (P coda, Z) | H 4, Vs 2.6 | 0.108 | 0.68 | 0.43 |
  | PrP (Z) | H 3, Vs 3.2 | 0.099 | 0.59 | 0.12 |
  | SrS, S-wavelet matched filter | H 5, Vs 3.0 | 0.016 | 0.34 | 0.13 |

- **Threshold (`tables/refl_inject.csv`).** An injected S copy at the predicted lag is found (p < 0.01) from
  A = 0.08 × direct S, and at p 0.07–0.11 at A = 0.05. Either polarity gives the same.
- **Station-side phases, not source-side.** SNB shows a later phase at ~1.85 s in every distance bin
  (`fig_refl_records.png`). It does not move out with distance, so it is a station-side reverberation. SNB PrP (p 0.005)
  and B011 SrS (p 0.02) are single-station results among 35 tests.
- **Fixed during the run (post-hoc, labelled).** The first SrP control was invalid: lags were taken after S, so energy
  locked to P was not permuted (it gave p 0.001 at the grid edge). Re-anchored on P and limited to ≥ 0.6 s after P, it is
  p 0.68.
- **Next adaptation.** Restack the existing windows at 1–3 Hz, where a 0.5–1 km gradational boundary still reflects.
  This needs re-reading the cut store, with no new download. Splitting (task 5) would need this phase. Direct S does not
  cross the LVL, so task 5 was not run.

## 3. Superposed epoch on the front: coda dv/v and decay (task 3) — `tables/{svi,wa_nw}/epoch_*`, `figures/fig_epoch_*.png`
- **SVI (464 series, 61 families), ETS days 0–20 minus between-ETS:**

  | window (s) | dv/v | target |
  |---|---|---|
  | 1–2 | +0.020 ± 0.017 | +0.40 |
  | 2–3 | −0.027 ± 0.021 | +0.46 |
  | 3–4 | +0.031 ± 0.018 | +0.46 |
  | 4–6 | +0.004 ± 0.020 | +0.41 |
  | 2–5 | +0.006 ± 0.011 | +0.45 |

  - Implied LVL change: +0.003 ± 0.041 %.
  - Gosselin frame (regional midpoint): 0–15 d +0.000 ± 0.012; −15 to −1 d +0.030 ± 0.016.
  - Pre-front − between: +0.04 ± 0.06.
  - Early (60–240 d) and late (> 240 d) between-ETS differ by 0.04 %, with late faster. That is the opposite sign to a
    decaying post-ETS speed-up.
  - Lapse pattern: flat, with no window near the target.
  - Decay slope 0–20 − between: +0.024 ± 0.025 1/s.
- **WA-NW (11 families with fronts):**
  - 0–20 − between: +0.055 / +0.064 / +0.014 / −0.134 / +0.021 % (same windows, σ 0.03–0.08).
- **Controls.**
  - The time-of-day control fails its 2σ rule at 1–2 s in SVI (−0.22 ± 0.07). No detection depends on it.
  - The front-scrambled control for pre − (0–20) is too narrow to use. Scrambled "pre" classes hold many more days than
    the real ones. It is labelled unusable.

## 4. Source-side Vp/Vs proxy (task 4)
- **ln(P/S) Z, 0–20 − between:** SVI −0.3 ± 4.3 %; WA +3.6 ± 4.2 %. The target is +2.5–5 % for a source inside the
  LVL, and the 2σ bound is ~9 %, so this measure cannot reach the target.
- **ln(P/S) H** is biased by noise on weak horizontal P (|projection|). It is not used.
- **S−P, 0–20 − between:** SVI −4 ± 11 ms (target ≈ 0, as geometry predicts). The ballistic round-3 bound stands.

## A. Propagating onset — `tables/*/novel_A*.csv`, `figures/fig_novel.png` (column 1)
- **Model.** Daily classes from day −10 to day +20. The change switches on with a front moving at 10 km/d (also 5 km/d)
  across an along-strike sensitivity σ of 7 km.

  | region | window | A | target | front-scrambled p |
  |---|---|---|---|---|
  | SVI | 2–5 s | +0.02 ± 0.09 % | +0.45 | 0.74 |
  | SVI | 4–6 s | −0.16 ± 0.08 % | +0.41 | 0.05 (wrong sign) |
  | WA | 2–5 s | −0.12 ± 0.10 % | +0.45 | 0.20 |
  | WA | 4–6 s | −0.27 ± 0.10 % | +0.41 | 0.025 (wrong sign) |

- **Reading.** No speed-up switches on with the front. A is weakly constrained: the pre-front days are few, so the
  2σ bound is 0.4–0.8 of the target.
- **The stronger statement is the ETS − between contrast (§3).** It bounds a change that is already present while the
  family is active.

## B. Recharge between ETS — `novel_B*.csv`, `fig_novel.png` (column 2)

| | SVI | WA-NW |
|---|---|---|
| slope 2–5 s (annual per site) | −0.006 ± 0.008 %/yr | **−0.33 ± 0.14 %/yr** |
| scrambled-control p | 0.47 | **0.004** |
| 4–6 s | −0.012 ± 0.015 | −0.25 ± 0.18 (p 0.03) |
| decay slope | −0.002 ± 0.021 1/s | +0.35 ± 0.32 1/s |
| days / intervals | 4341 / 19 | 738 / 18 |

- **Target.** −0.4 %/yr if a +2 % LVL change decays linearly over a 14-month cycle. The decaying form (τ 45 d) predicts
  +0.45 × e^(−t/45 d). WA's exp-form amplitude is +1.7 ± 0.7, the predicted sign.
- **Checks (post-hoc, labelled; `novel_B_check.csv`).** Leave-one-site-out, calendar trend and scrambled control are
  described in the short summary at the top.
- **Caveats.**
  - The WA geometry uses the SVI kernel scale.
  - The 2–4-month bin has only 52 days.
  - Burst locations may shift through the cycle (source composition), which is untested.
- **Next tests.**
  - The same fit per family, plus a source-locus column (ballistic S−P of the same days).
  - WA-SW, which needs its round-3 store.
  - SVI per zone, because SVI's null may average zones with opposite behaviour.

## C. Source properties across the front — `novel_C.csv`, `fig_novel.png` (column 3)
- **Method.** Per-day scalars, size-adjusted (ln events/day, ln weight), each class minus between-ETS.
- **SVI, behind (2–20 d) − at the front (day 0):**
  - spectral ratio ln(E 5–8 / E 2–4 Hz): +0.04 ± 0.05;
  - S width: +0.013 ± 0.008 s;
  - P/S: +0.001 ± 0.011;
  - decay slope: −0.036 ± 0.032.
- **WA, same contrast:** −0.01 ± 0.03 / +0.003 ± 0.011 s / −0.015 ± 0.013 / −0.037 ± 0.029.
- **Day-0 decay slope (SVI).** It is shallower than between-ETS (+0.075 ± 0.022 1/s, 3.4σ, positive at 4/4 sites).
  - Formally this passes the C rule. But it sits only on the busiest days: day 1 is +0.02 and days 2–3 −0.01.
  - Pre-front days read +0.17 ± 0.08. The same pattern is in the busy-episode tremor background found in S6.
  - Read as a tremor-background effect until the time-of-day control days are run through the same scalar (next test).
- **Pre-front P/S RMS.** +0.10 ± 0.03 in both regions (231 / 138 days). The pre-front stacks are small, and RMS P/S is
  noise-biased upward. Not interpreted.
- **S6 per-detection size, behind − at the front:** +0.031 log10 (CI −0.009 … +0.070), not seen.
- **Expected size (order of magnitude).** A 1–10 MPa Pf drop could change the corner frequency by 10–30 %, which is
  0.1–0.3 in the ln spectral ratio. The 2σ bounds of ~0.1 (SVI) reach the bottom of that range.

## D. Scaling with episode size — `novel_D*.csv`, `fig_novel.png` (column 4)
- **Model.** Per-major contrast of days 0–20 − between vs duration Du (per 10 d).
- **SVI:** +0.019 ± 0.046 (15 majors, 11 with Du).
- **WA:**
  - Du: +0.029 ± 0.023, and +0.042 ± 0.028 with log n̄;
  - tremor count T: +0.030 ± 0.020 / +0.048 ± 0.026;
  - Fs: +0.25 ± 0.12 / +0.33 ± 0.13 (2.1–2.5σ).
- **Reading.** WA is positive in every size measure, matching the round-7 WA hint and the fault-valve sign. None reaches
  3σ. Target if the change scaled with duration: 0.30 %/10 d. The observed slope is ~10 % of that, because the mean
  ETS − between level itself is ~0.

## Overall reading and next adaptations
- **SVI.** The LFE coda bounds an ETS-locked LVL change near the families to < 0.085 % (25× below Gosselin's size),
  including Gosselin's footprint (< 0.2 %).
- **Why that may not contradict Gosselin.**
  - Every LFE day, between-ETS bursts included, is a day of local slip. If each slip event drains the valve, the LFE coda
    always samples the drained state and never the recharged one.
  - Gosselin's change may also sit at RF wavelengths (0.05–0.7 Hz), deeper or broader than the coda of 2–8 Hz LFEs sees.
- **WA.** The coda slows through the inter-ETS interval by about the predicted amount. That is the one place where the
  LFE data look like valve recharge. It is the cheapest strong lead:
  1. per-family fits, with a source-locus control;
  2. WA-SW replication;
  3. the time-of-day control through B;
  4. a WA LVL forward kernel.
- **Other next steps.**
  - 1–3 Hz restack for the reflection.
  - Intraslab earthquakes below the LVL: their up-going S crosses the layer, giving a 25–30 ms target for +2 %. This
    needs waveforms (ask before download). The eq-pick fetch may already give a first look.

## Figures
- `fig_fwd.png`: forward target per window, and predicted later phases.
- `fig_refl_records.png`: record sections per site with predicted SrS curves.
- `fig_refl_scan.png`: model scan, and moveout stacks against the scrambled control.
- `fig_epoch_{svi,wa_nw}.png`, `fig_epoch_lapse_*.png`: superposed epoch and lapse profile, points coloured by weight,
  forward target dotted.
- `fig_novel.png`: A–D for both regions.

## 7. Confirmation of the WA recharge candidate (PLAN §7; written ~10:45 am PT; `tables/wa_nw/r2_*`, `fig_recharge_wa_nw.png`)
**Not confirmed. The slope belongs to the 11 families of §6 and is absent in the other 39 WA-NW families.**
- **All 50 families, every store series (2763 between-ETS day measurements, 18 intervals):** −0.055 ± 0.053 %/yr.
  Scrambled p 0.15. 2σ bound 0.16 %/yr, about 0.4 of the −0.4 %/yr target.
- **C1, families.**
  - The 11 families of §6: −0.33 ± 0.14 (reproduced).
  - The 39 others: **+0.008 ± 0.059**.
  - Per family: 22 of 32 negative (binomial p 0.025). C1 formally passes on this sign count, but the size collapses.
- **C2, source / family mix.**
  - Balanced families: −0.059 ± 0.049.
  - With S−P as a covariate: unchanged (−0.055).
  - S−P drift: −5 ± 6 ms/yr.
  - Family-position drift: ≤ 1.8σ.
- **C3, stations.** Boreholes −0.067 ± 0.058; surface −0.026 ± 0.14.
- **C4, season.** Per-station annual + semi-annual: −0.052 ± 0.052 (−1.0σ). Season-matched scrambled p 0.09.
  **C4 fails.**
- **Lapse 1–2 / 2–3 / 3–4 / 4–6 s:** −0.13 / −0.04 / −0.08 / −0.04 %/yr. Too weak to tell the LVL, shallow and 25–45 km
  patterns apart (χ² 1.5 / 2.6 / 0.7).
- **Shape.**
  - Days ≥ 60 d: the exponential prefers τ 20 d. That is just the earliest bins (ΔAIC −2.1).
  - With 21–60 d days added: τ 120 d, ΔAIC +1.5 (linear preferred).
  - No clear recovery curve.
- **Tie-in: post-ETS level (60–180 d) against the size of the preceding ETS.**
  - Du −0.025 ± 0.042, T −0.044 ± 0.041, Fs −0.26 ± 0.16.
  - Not positive, so it does not follow the round-7 size hint.
- **Replication.**
  - WA-SW is impossible with existing data: 0 between-ETS family-days at P ≥ 0.5–0.75 (counted before measuring).
  - SVI by era (2003–2013 vs 2014–2026): no replication in either era (§7c below).
  - The 1–3 Hz reflection restack: no reflection, and it is less sensitive than 2–8 Hz (§7b below).
- **Are the 11 a spatial patch?** (`scripts/s12_cluster11.py`, `tables/wa_nw/r2_cluster11.json`, `figures/fig_cluster11.png`;
  run ~10:50 am PT at the coordinator's request.) **No, they are scattered.**
  - Compactness: RMS distance to their centroid is 32 km, against a median of 29 km for 1000 random 11-family sets
    (p 0.73).
  - Their slope is extreme among the random sets (p 0.009).
  - Location: all in WA (median 47.80°N, 122.88°W), ~99 km from PGC (the others ~96 km). Not on the PGC-side LVL and not
    at the Strait.
  - Depth: median 42.8 km (37.6–57.5) vs 39.4 km for the others. They lean downdip.
  - The other 13 families deeper than 41 km give −0.03 ± 0.09 %/yr; all families < 41 km give +0.03 ± 0.08.
  - Depth does not carry the slope. Reading: family-set-specific scatter, most likely noise. Not a localised finding.

## 7b. 1–3 Hz LVL-base reflection restack (PLAN §7; run Wed 2026-10-14 2:52–2:53 pm PT; `tables/refl_*_lb*`)
**No reflection at 1–3 Hz either, and this band is less sensitive than 2–8 Hz.**
- **What was run.** SVI near stations (PGC, B009, B010, B011, B926, SNB, VGZ). Per (family, station) up to 100 of the
  highest-P repeats (PLAN said up to 150; the queue used `--nmax 100`), re-read from the existing raw day files, filtered
  1–3 Hz, realigned on S, stacked: 49,681 windows → **380 pairs** (2–8 Hz run: 497). Then the §2 moveout scan, scrambled
  control (lags permuted within station, 1000 draws, scan-corrected) and injection, with the envelope smoothing scaled to the band.
- **Expected (PLAN §7, §1).** SrS from the LVL base ~1.8 s after S (1.4–2.3 s); a sharp base gives ~0.13 × the direct S.
  A ~1 km gradational base still reflects at 1–3 Hz (S wavelength ≈ 1.5 km). A +2 % LVL change would move it by −59 ms.
- **Scan:**

  | phase | best model | peak | scan p | nominal-model p |
  |---|---|---|---|---|
  | SrS LVL base (H envelope) | H 4, Vs 2.6 | 0.095 | 0.19 | 0.90 |
  | SrS Moho | H 2.5, Vs 3.6 | 0.089 | 0.20 | 0.64 |
  | SrP (P coda, Z) | H 3.5, Vs 3.8 | 0.194 | 0.30 | 0.55 |
  | PrP (Z) | H 5, Vs 3.8 | 0.235 | 0.34 | not valid: nominal PrP is 0.6–1.0 s after P, and the wider 1–3 Hz flanks (0.45–0.9 s) reach into the P pulse, so only 8 pairs have a value (10 needed); the printed p 0.001 is an artifact of an empty value |
  | SrS, S-wavelet matched filter | H 6, Vs 3.6 | 0.017 | 0.18 | 0.79 |

- **Threshold (`refl_inject_lb.csv`).** An S copy injected at the nominal lag is found at p 0.02 from A = 0.12 × direct S
  and at p 0.002 from A = 0.18 (either polarity). At 2–8 Hz the same test found it from 0.08. The longer 1–3 Hz S pulse and
  the wider smoothing blur the reflection into the S coda. At the nominal lag the data sit slightly below the scrambled
  level (−0.15 vs q95 −0.07, p 0.90), i.e. no excess energy there.
- **Reading.** A sharp LVL base (~0.13) sits at the 1–3 Hz detection limit, and a gradational base would be weaker still.
  So 1–3 Hz cannot show or exclude a gradational base. Together with 2–8 Hz: a reflection stronger than ~0.08 × direct S
  (2–8 Hz) is excluded beneath these families; weaker ones are not tested. The −59 ms timing change cannot be measured
  without the phase.
- **Single stations (`refl_station_lb.csv`, 70 tests).** PGC SrS best-model p 0.005 and SNB matched filter p 0.005.
  About 0.35 such values are expected by chance among 70. Single-station results, like B011/SNB at 2–8 Hz.
- **Fix (2026-10-14).** The 10-13 scan crashed (rc 1) on 10 of the 380 pairs that have no S−P value in `offsets_sp.csv`.
  For them the P-anchored window was empty (`resid_logenv_P`, `max()` of an empty array). They now get no P-side residual.
  The S-side phases are unchanged, and the 2–8 Hz run (all S−P finite) is unaffected. Rerun only the scan + injection
  (`scripts/s12_queue4.sh`, marker `QUEUE4_DONE`).

## 7c. SVI recharge slope by era (PLAN §7 replication; `tables/svi/r2_verdict.json` era_* keys, check `r2_era_check.json`)
**No replication in SVI, in either era.** The replication rule was a negative slope at ≥ 2σ.
- **What it measures.** Between-ETS LFE days (≥ 60 d after the previous front), coda dv/v 2–5 s, slope against time
  since that ETS, within each family-station series, cc²-weighted, with an annual cycle per site. σ = bootstrap over
  inter-ETS intervals.

  | era | slope (%/yr) | σ (bootstrap; 2 runs) | scrambled control p | days / series / families / intervals |
  |---|---|---|---|---|
  | 2003–2013 (backward catalogue) | −0.004 | 0.017–0.027 | 0.74 | 1,989 / 350 / 74 / 7 |
  | 2014–2026 | −0.005 | 0.013–0.020 | 0.66 | 2,661 / 425 / 70 / 7 |
  | all SVI (§7 headline) | −0.003 | 0.007 | 0.73–0.76 | 4,650 / 603 / 83 / 14 |

- **Size against the target.** The valve target is −0.4 %/yr. The 2σ bounds are about 0.05 %/yr (2003–2013) and
  0.04 %/yr (2014–2026), 8–10× below the target in each era. The two eras agree with each other and with zero.
- **σ range.** There are only 7 intervals per era, so the interval bootstrap moves by ~1.5× between runs (200 draws).
  The table gives both runs. The scrambled-control SDs (0.012 / 0.010 %/yr) are in the same range.
- **B009 veto (2026-10-14; `svi_science/instrument_vetoes.csv`).** Dropping B009 2007-01-08…2009-10-01 removes 28 days.
  2003–2013 moves from −0.0037 to −0.0046 %/yr; the full SVI set moves from −0.0027 to −0.0033. No change in the reading.
- **Scope.** This is the LFE coda volume (the path and ~10 km around the families, weighted to the S-anchored coda),
  on the time scale of one inter-ETS cycle (60 d to ~14 months after a front). It compares days within the same
  family-station series. It does not test a change confined to the post-ETS weeks (< 60 d) or outside the coda's
  sampling volume.
- Post hoc (labelled): the per-era scrambled control and the veto refit (`scripts/s12_era_check.py`) were added 10-14,
  after the era numbers existed.

