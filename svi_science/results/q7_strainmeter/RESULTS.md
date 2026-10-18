# Q7 results — borehole strain on the deep between-episode LFE bursts

Run 2026-10-07, ~5:10–5:45 pm PT. Rules were frozen in `PLAN.md` (4:55 pm PT) before any strain number was computed.
Two changes were made after the first look, and both are declared in "Deviations" below. Scripts:
`svi_science/scripts/q7_download.py`, `q7_process.py`, `q7_strain.py`, `q7_figfix.py` (display only). The log is
`q7.log`; its clock is UTC. No Merlin review yet, so treat this as provisional.

## Plain answer
**No, not at a size the strainmeters can see.** There are 95 tremor-free deep bursts between episodes (SVI 38,
WA 57, 2014–26). Stacking all 20 Olympic, Puget and Vancouver Island borehole strainmeters around those bursts shows
no strain in the thrust-slip pattern. Averaged over the bursts, the stack is +0.05 ± 0.05 of what one Mw 5.2 slip per
burst would give (z = +1.25 against random days, p = 0.11).

The method itself works:
- **It sees ETS.** Every big ETS shows a clear strain step: stack z = +22, and 23 of 30 single episodes are above 2σ.
  The measured size is Mw 6.0–6.7, median 6.4.
- **It sees small episodes.** Even the smallest catalogued episodes (5–20 days) give a stacked Mw 5.7 at z = +14.
- **It recovers injected slip.** An injected Mw 4.8 slip at every burst is recovered in 98 % of stacks.

So the bursts are, on average, **smaller than Mw ≈ 4.7 each** (3σ bound, about 1.3 × 10¹⁶ N·m).

That bound rules out bursts that follow the usual slow-earthquake scaling (Ide et al. 2007) at its central or high end:
- **Central rate (10^12.5 N·m/s, per-burst Mw ≈ 4.6):** injected, it would give stack z ≈ 10. We see z = 1.25.
- **Low end (10^12.0 N·m/s, Mw ≈ 4.3):** only 60 % of injected stacks are recovered. This end is not excluded.

**In plain words.** If these bursts are slow slips, each one is at most about Mw 4.7, and probably smaller than the
scaling relation predicts for its duration. They are not hidden Mw 5–5.5 mini-ETS events. Bursts of that size would
have shown up in the stack, and many of them would have shown up one at a time.

**Pre-registered verdict:**
- (i) Burst stack z > 3: **FAIL** (z = +1.25).
- (ii) Positive control: **PASS** (z = +22, 77 % of single episodes above 2σ).
- (iii) Mw 5.2 injection recovered in > 80 % of stacks: **PASS** (100 %).

So the outcome is **a bound, not a detection.**

## Numbers (`stacks.csv`, `q7_summary.json`)
**Estimator.** For each gauge, the strain step across the burst is fitted as a line plus a step on ±2 days of hourly,
tide- and pressure-corrected data; the burst interval itself is excluded.

**Combination.** The steps from all gauges and stations are combined by GLS onto the expected gauge pattern of thrust
slip on mesh patches within 10 km of the burst. That pattern comes from a cutde half-space model with μ = 40 GPa and
nominal GTSM coupling.

**Scrambled control.** Each burst is moved to 1,000 random days with:
- the same time of day;
- a day of year within ± 30 days;
- no episode within ± 2 days and no burst within ± 3 days;
- similar station coverage.

| set | n | stack (Mw 5.2 = 1) | random days, mean ± σ | z | p | 3σ bound per burst |
|---|---|---|---|---|---|---|
| **primary: SVI + WA tremor-free bursts 2014–26** | 95 | +0.055 | −0.003 ± 0.046 | **+1.25** | 0.11 | **Mw ≤ 4.73** |
| SVI tremor-free | 38 | +0.022 | −0.004 ± 0.058 | +0.44 | | Mw ≤ 4.73 |
| WA tremor-free | 57 | +0.11 | −0.003 ± 0.081 | +1.44 | | Mw ≤ 4.90 |
| WA tremor-free 2009–14 | 34 | +0.048 | −0.011 ± 0.10 | +0.57 | | Mw ≤ 4.91 |
| all tremor-free 2009–26 | 129 | +0.054 | −0.005 ± 0.043 | +1.37 | | Mw ≤ 4.71 |
| big bursts (≥ 3 LFEs) | 34 | +0.073 | −0.002 ± 0.14 | +0.52 | | Mw ≤ 5.00 |
| with PNSN tremor ≤ 100 km ± 6 h | 22 | +0.10 | −0.006 ± 0.13 | +0.84 | | Mw ≤ 5.00 |
| without that tremor | 107 | +0.047 | −0.004 ± 0.045 | +1.12 | | Mw ≤ 4.71 |
| broad set (local tremor allowed), both regions | 156 | +0.12 | −0.003 ± 0.045 | +2.81 | | Mw ≤ 4.81 |
| broad, ≥ 5 LFEs | 39 | +0.41 | −0.03 ± 0.17 | +2.56 | | Mw ≤ 5.18 |
| window ± 1 day (primary) | 97 | | | +0.44 | 0.33 | Mw ≤ 4.68 |
| window ± 5 days (primary) | 83 | | | +2.34 | 0.011 | Mw ≤ 4.92 |

- **Robust re-statement (deviation 2 below).** With fixed empirical weights, the primary z is −0.79 (bound Mw ≤ 4.58).
  The mean of the per-burst z values gives −0.17. None of the three statistics reaches 3 for the primary set.
- **Two descriptive sets come near 3, but neither is a detection:**
  - Broad, ≥ 5 LFEs: GLS z = +2.56, empirical z = +3.15, mean-z z = +2.52.
  - Broad, both regions: GLS z = +2.81, but empirical z = +0.06.
  - About 12 sets × 3 statistics were looked at, so one crossing of 3 is expected by chance.
  - The broad bursts sit inside local tremor bursts. The control days are not matched for time since the last
    episode, so these leads are confounded. If real, they correspond to about Mw 4.9 per burst.
- **Single bursts** (`bursts_strain.csv`):
  - A single burst would need Mw 5.2–5.5 to clear 3σ (median 5.21 SVI, 5.43 WA). The range is 5.0–5.8.
  - None of the biggest bursts shows anything. These are SVI 2015-04-12 and 2018-10-21, and WA 2017-11-10 and
    2018-11-18 (the Q5 GNSS candidate, z here −0.49). All have |z| < 1.1 and 3σ thresholds of Mw 5.4–5.7.
  - Only one tremor-free burst exceeds 3σ: SVI 2024-07-10 17:41 PT (2 LFEs, 21 h, z = +3.0, would be Mw 5.2). It
    starts 3 days before the WA/NW episode of 2024-07-13. The SVI bursts were only screened against SVI episodes, so
    this is most likely the onset of that episode and not a burst signal.

## Controls
**Positive control** (`positive_control_episodes.csv`, `figures/q7_positive_control.png`):
- **Method:** the step from 10 days before to 10 days after each episode, with the same GLS, against same-length
  random windows away from major episodes. The template is thrust slip spread in proportion to the PNSN tremor in the
  episode.
- **P1, big ETS** (in the frozen-major list, ≥ 20 days, 2010–2026): n = 30.
  - Stack z = +22.0 (empirical +26.3). 23 of 30 single episodes are above 2σ.
  - SVI 10 of 11, WA 13 of 19.
  - Strain moment Mw 6.0–6.7, median 6.4. This is close to the published Mw ~6.5–6.8 of the large northern
    Cascadia ETS.
  - **PASS.**
- **P2, smallest catalogued episodes** (minor, 5–20 days): n = 86.
  - Stack z = +13.8, stacked Mw 5.7 per episode.
  - 16 of 86 single episodes are above 2σ.
  - Detected.

**Calibration check (check C).**
- We predicted the gauge steps from the Gualandi (2025) GNSS slip model over each big ETS.
- The GLS amplitude of the observed steps relative to that prediction has a median of 0.92 (WA 0.79, SVI 0.92). The
  nominal calibration is therefore right to within about 20 %.
- Gauge by gauge, the raw steps are noise-dominated over 40–60-day windows: r = 0.03, and the weighted gain is
  0.39–0.63. Some stations agree well: B003 (r 0.47) and B004 (0.45).
- Applying the WA empirical gain moves the bound from Mw 4.58–4.73 to Mw 4.65–4.94.

**Injection** (`injection.csv`, `figures/q7_injection.png`). The expected strain was added as a ramp over the burst
duration to the random-day windows.

| injected Mw per burst | stack recovered (z > 3) | single burst recovered |
|---|---|---|
| 4.6 | 36 % | 1 % |
| 4.8 | 98 % | 2 % |
| 5.0 | 100 % | 6 % |
| 5.2 | **100 %** (median z 22) | 28 % |
| 5.4 | 100 % | 65 % |
| 5.6 | 100 % | 93 % |
| 5.8 | 100 % | 100 % |

- The step remover takes away about 6 % (Mw 5.8) to 26 % (Mw 6.0) of a short injected burst. It is fully
  transparent below Mw 5.8.
- A numerical check, done by adding the ramp to the raw series and re-running the de-stepping and the fit (12 cases at
  Mw 5.6), recovers 1.00–1.03 of the injected moment.

## Figures
- `figures/q7_strain_stack.png`. Left: the stacked, template-projected strain around the 95 bursts (red), in units of
  one Mw 5.2 slip, against 40 random-day stacks (grey 95 % band). Right: stacked moment per burst for each set
  (dots), the random-day distributions (violins), and control + 3σ (bars). The dashed line is Mw 5.2 per burst.
- `figures/q7_positive_control.png`. Left: each big ETS, projected strain in 10¹⁸ N·m; every trace steps up across the
  episode. Right: single-episode z for big and minor episodes.
- `figures/q7_injection.png`: injection–recovery curves (stack and single burst).
- `figures/q7_map_sensitivity.png`: strainmeters and bursts, coloured by the single-burst 3σ threshold in Mw.
- `figures/q7_diag_ets2010.png`: diagnostic plot of raw corrected gauges through the 2010-08 ETS.

## Comparison with earlier work
- **Q5 (GNSS) and Frank et al. (2015, 2016).** Q5 bounded each burst at Mw ≤ 5.9–6.0 with GNSS. The strainmeters
  tighten this by about 1.2–1.3 magnitude units, to Mw ≤ 4.7. That is now below the scaling prediction (Mw 4.6–5.2).
  - Frank's GPS stacks on LFE bursts (Guerrero) found hidden slip only by stacking very many burst days. Here a
    sensitivity ~1.3 units better shows that the deep Cascadia tremor-free bursts carry no slip above Mw ~4.7 each
    on average.
- **Earlier PBO strainmeter studies (e.g. Hawthorne & Rubin 2013; Dragert & Wang 2011).** PBO strainmeters record
  ETS slip. Our positive control reproduces clear ETS strain steps of the right size (Mw ~6.4) with our own raw-gauge
  processing.
- **Wech & Creager (2011).** They describe inter-ETS tremor bursts as small slip events. Our bound does not contradict
  slip in those bursts, but it caps each burst's moment at about 10¹⁶ N·m. That is about 1/1000 of a big ETS, and
  equals roughly 1 mm of slip over a 10-km patch. For the broad bursts (local tremor allowed) there is a hint at the
  Mw ~4.9 level (z 2.5–3.2) that is not significant.

## Deviations and caveats
1. **Data.** EarthScope Level-2 processed strain needs a login (both archive URLs redirect to the login page), so it
   was not used. We used public FDSN raw gauges (RS1–4, 10-min) and RDO pressure, with our own processing:
   linearisation, hourly means, per-year 8-constituent tides plus pressure admittance, spikes, and steps.
   - The download is 269 MB in `data/bsm/` (`SOURCE.txt`, `MANIFEST.csv`).
   - B006 and B011 have unknown gauge azimuths and were dropped. 20 stations were used.
2. **Changes after the first look.** These were made because the pre-registered positive control at first FAILED
   (z = +0.78), not because of the burst result.
   - **(a) Dead-gauge QC.** A flat gauge (B926 g3, B012 g3 from 2015, B943 g1 2010–12) had near-zero variance. GLS
     then gave it enormous weight. A gauge-year is now dropped if its residual RMS is below 0.25 × the station median.
     We also added a 2 nanostrain covariance floor. After this, the positive control passes with the pre-registered
     statistic, and that statistic is the one used for the verdict.
   - **(b) Robust re-statements.** We added two more versions of every stack: fixed empirical inverse-variance weights,
     and the mean of the single-burst z values. These are columns `emp_*` and `mean_single_z*`.
   - The burst conclusion is the same under all three statistics.
   - Bursts were seen in a debug run with 60 draws before the final run. The final run changed nothing in the burst
     rules.
3. **Other changes.**
   - Events with fewer than 100 random days matching the station coverage were dropped: 2 big ETS and a few minor
     episodes, and no primary burst.
   - The positive control starts in 2010-01 because PNSN tremor is needed for its template, so the 2009-05 ETS is not
     included.
   - Check C reports the GLS-scale median as its main number. The pre-registered raw regression gain, 0.48, is given
     too.
4. **Limits.**
   - Tidal calibration is nominal, not tidally calibrated per gauge; check C says it is within about 20 %.
   - Ocean loading is fitted per constituent, not modelled.
   - Rain and groundwater are not corrected. They are in the random-day spread, which makes the WA and B01x stations
     noisy.
   - The template assumes slip at the bursting family's interface point over a ~10-km patch. Slip much farther away,
     for example up-dip in the ETS zone, would be seen less well.
   - Bursts lasting hours that leave no LFEs elsewhere cannot be distinguished from a smaller, shorter event at the
     same moment.
   - The "Mw ≤ 4.7" bound is a stack average. Individual bursts could be larger if most are much smaller, but none of
     the biggest exceeds about Mw 5.4–5.7 (their single 3σ thresholds).

## Files
- `PLAN.md`, `RESULTS.md`, `q7.log`, `run.out`, `q7_summary.json`
- `stacks.csv`, `injection.csv`
- `bursts_strain.csv`: every burst, with m̂ in N·m, the random-day σ, z and single-burst 3σ Mw
- `positive_control_episodes.csv`, `checkC_per_station.csv`, `checkC_obs_pred.csv`
- `controls.npz`, `download*.log`
