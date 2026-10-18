# S11 round 3: kernel-robust deep term, REDPy estimators in Q3, along-strike zones (SVI and WA-NW)

Written Mon 2026-10-12, ~9:30 pm PT. The rules are in `../PLAN.md` §6, written at 7:25 pm PT before any round-3 number, plus
Amendment 1 (~8:00 pm PT, after a code test). The queue ran 7:29–8:59 pm PT. Decision log: `../NOTES.md` (decisions 21–23).
Tables: `{svi,wa_nw}/tables/`. Figures: `figures/fig_r3_{kerr,q3,zones}.png`.

Conventions: dv/v in %. "Scrambled control" = station data shifted to another time of the same day. Cells: SVI 115 families,
n ≥ 20 (3,531 cells); WA-NW n ≥ 20 (2,466 cells). Everything here is provisional.

## In plain words

**1. Making the deep term independent of the kernel: not achieved. The obstacle is now measured, and the data prefer one kernel.**
- **Leave-one-kernel-out.** We averaged five kernel models, made a fake signal with the sixth, and checked recovery. No setup
  passes in every fold, in either region, under either reading of the rules.
  - **What fails is mostly the size of a deep change, not whether it is seen.**
    - A deep step comes back at 0.73–1.32 of its size (SVI) and 0.60–0.91 (WA). It depends on how much deep sensitivity the
      true kernel has compared with the assumed one.
    - So a deep bound is uncertain by about ±30 % from the kernel. It does not appear or disappear.
  - **The second failure is a shallow-to-deep leak of the yearly cycle.** It appears when the true near-surface scattering is
    much stronger than assumed: generated with two3_20 at SVI, 0.43–0.49; generated with two5_40 in WA, 0.36–0.41. The rule
    is < 0.25.
- **The data choose a kernel.** Fitting the early (2–3.5 s) and late (3.5–6 s) coda together, both regions prefer strong
  near-surface scattering:
  - WA-NW: the two3_20 model (scattering length 3 km in the top 3 km) wins by Δχ² ≥ 9 over every other model and by 74 over
    ell10;
  - SVI: two3_20 ≈ two5_40 (Δχ² 0.6), both far better than any uniform-scattering model (Δχ² 5–24).
  - With that kernel, a signal made with the same kernel is recovered cleanly: deep step 0.93 (SVI) and 0.81 (WA); yearly
    leak 0.01 and 0.03.
  - The pre-registered rule demanded robustness to all six kernels, so the lapse-chosen setup does not pass either. Its weak
    point is a true kernel with weak near-surface scattering (ell40/ell80): deep step 0.50–0.64, leak 0.29–0.33.
- **The yearly deep share under the alternative kernels.**

  | region | nominal | 6-kernel average | lapse-chosen |
  |---|---|---|---|
  | SVI | 0.114 ± 0.037 % | 0.071 ± 0.026 % | 0.052 ± 0.015 % |
  | WA-NW | 0.164 ± 0.071 % | 0.160 ± 0.066 % | 0.151 ± 0.043 % |

  - In both regions it is still not called deep. The station pattern rejects the deep form (SVI χ² 31 vs 20 for shallow; WA
    both ~150). The odd/even family halves disagree on the peak date (SVI 260 vs 166 days). The date permutation gives
    p 0.22 (SVI) and 0.15 (WA). WA's lapse ratio (1.24) still matches the shallow prediction.
  - The yearly cycle is a measured value with 2σ bounds: SVI < 0.12 % (6-kernel average), WA < 0.29 %.
  - **No decade-scale deep change under any kernel.** Q2 ranges −0.07 to +0.06 % across all setups (σ 0.03–0.10 %).

**2. REDPy findings in Q3: the all-pairs graph fixes the deep-pulse over-read in both regions. The shallow leak drops only a
little.**
- **Graph estimator (B2).** Episode stacks of ≥ 50 events are compared with each other directly, with no shared reference. A
  deep pulse in one episode now comes back at its true size: SVI 0.98 (0.81–1.00), was 1.14; WA 0.91, was 1.05. Both are inside
  the 0.8–1.1 rule.
- **Weighting each LFE by its waveform match (cc weights, B1).** This lowers the per-major σ by 7–10 %: SVI 0.091 → 0.082 %;
  WA 0.131 → 0.122 %.
  - WA then reaches the precision target: 8 of 32 majors have σ ≤ 0.10 % (was 3).
  - B1 still uses the leave-out reference, so it still over-reads at SVI (1.15).
- **Shallow leak** (a common shallow change during one episode counted as deep):
  - SVI 0.45 → 0.38 (B2); WA 0.60 → 0.47.
  - Still above the 0.35 rule. The drop is mostly the removed over-read; the rest is the model itself. The deep and shallow
    weights differ too little between cells to separate a one-episode shallow change.
- **Adopted (pre-registered rule): B2 in both regions.** It passes deep-pulse recovery, scrambled control, leave-one-site-out,
  and σ ≤ 1.1 × B0. Its composite σ is 0.035 % at SVI and 0.017 % at WA.
  - Caveat: the graph's own scrambled control is measurable for only 1 SVI major and 7 WA majors, because it needs ≥ 50 control
    events per episode.
  - Cost: the graph uses only episodes with ≥ 50 events, so WA keeps 24 of 32 majors and SVI 55 of 115 episodes.
- **Composites stay at zero.** Average major: SVI +0.007 ± 0.035 %, WA −0.014 ± 0.017 %. So an average slow-slip episode
  changes the deep-weighted velocity by less than ~0.07 % (SVI) and 0.03 % (WA), at 2σ.
- **One single episode stands out in the graph: SVI 2019-11, +0.53 ± 0.15 % (3.6σ).**
  - The P-weighted graph gives +0.52 ± 0.17. The old raw path gave +0.19 ± 0.12.
  - It survives leaving out any one instrument (0.43–0.57 %) and either family half (0.42 / 0.52).
  - Without the PGC site it becomes 0.40 ± 0.32 %, i.e. not resolved.
  - It rests on 28 cells, mostly the PGC-site boreholes and B926.
  - It does NOT meet the wording rule: the prior test fails (0.28σ / 0.63σ), the shallow leak is 0.38, and the B2 family-split
    gate fails. So it is a candidate, not a finding.
  - 2019-11 is also the episode with a PGC timing glitch (round-1 ballistic, −6 ms instrument term). That glitch should not
    survive per-event S realignment, but it must be ruled out.
  - The WA graph has no single major beyond 2σ. The 3σ WA values of round 2 (2014-07/08) fall out of the graph: their
    episodes have < 50 events per node.

**3. Along-strike zones: forecast met in 2 SVI zones and 3 WA zones. The data show no difference between zones.**
- **Forecast (design only, no dv/v).**
  - SVI with 3 zones fails the forecast. The south-east zone composite σ is 0.17 %: those families' deep sensitivity lies
    mostly along the paths toward the PGC-site stations, in the next zone. With 2 zones the forecast is 0.100 and 0.068 %, so
    the run went ahead.
  - WA with 3 zones: 0.081 / 0.053 / 0.094 % — run.
- **Result.**

  | region | zone composites | contrast | equal-zones test |
  |---|---|---|---|
  | SVI | SE −0.05 ± 0.10 %, NW +0.05 ± 0.07 % | −0.10 ± 0.13 % | p 0.45 |
  | WA | south −0.07 ± 0.08 %, centre +0.03 ± 0.05 %, north +0.02 ± 0.09 % | ±0.11–0.12 % | p 0.61 |

  - So no part of the fault behaves differently on average during slow slip. Zone-to-zone differences are below ~0.35 % at 2σ.
- **Gates.**
  - Scrambled control and leave-one-site-out pass in both regions.
  - The zone pulses separate cleanly at SVI: recovery 0.86–0.91, cross-zone leak ≤ 0.03. A common shallow pulse barely touches
    the zone contrast (0.05), so zone contrasts are almost immune to the shallow-leak problem.
  - In WA the zone-pulse recovery is only 0.23–0.32. This is a design flaw of the test, not of the zones. The pulse was put
    into every major at once, but ΔD compares each major with its neighbours, and in WA those neighbours are mostly other
    majors. So an every-major pulse cancels itself.
- With the B2 estimator, WA's north zone is +0.16 ± 0.08 % (equal-zones p 0.07). That is not significant.

**Relation to the distant-station agent's "regional shallow per episode" model** (near stations: ±0.18 %/major, leak 0.07).
That model is the only lever so far that removes the one-episode shallow leak. Here the leak stayed at 0.38–0.47 under every
kernel and estimator. The model costs about ×2 in σ compared with B2 (0.085 % per major at SVI). The zone contrasts here are
a second, cheaper route: a common shallow change cancels in the difference between zones (leak 0.05 at SVI).

## Tables

### Task 1: the six injection rules for out-of-sample generating kernels (`kerr_rules.csv`, `kerr_verdict.json`)

| region | in-sample (nominal ell40 → ell40) | leave-one-kernel-out: no fold passes both Q1/Q2 and Q3 (Amendment 1 or strict) | lapse choice (weights) | winner |
|---|---|---|---|---|
| SVI | R1 0.96, R2 0.060 (fail), R3 0.04, R4 1.14 (fail), R5 0.44 (fail), R6 0.02 | worst fold: two3_20 generated (R1 1.32, R3 0.43, R4 1.33) | two3_20 0.55, two5_40 0.41, ell40 0.04 | none |
| WA-NW | R1 0.82, R2 0.024, R3 0.10, R4 1.05, R5 0.60 (fail), R6 0.00 | worst fold: two5_40 generated (R1 0.67, R3 0.36, R4 0.79, R5 1.07) | two3_20 0.99 | none |

- Folds that pass under Amendment 1:
  - SVI: Q1/Q2 only in the ell10 fold; Q3 in the ell20, ell40 and ell80 folds.
  - WA: Q1/Q2 only in the ell20 fold (also strict); Q3 in the ell40 and ell80 folds.
- Q1/Q2 under the nominal and 6-kernel-average setups (`kerr_{NOM,AVG6}_gates.csv`; all setups in `kerr_values.csv`):

| quantity | SVI nominal | SVI average of 6 | WA nominal | WA average of 6 |
|---|---|---|---|---|
| Q1 (%) and peak day | 0.114 ± 0.037, day 229 | 0.071 ± 0.026, day 214 | 0.164 ± 0.071, day 195 | 0.160 ± 0.066, day 202 |
| Q2 (%) | +0.059 ± 0.074 | −0.020 ± 0.053 | −0.003 ± 0.089 | −0.044 ± 0.081 |
| G1 / G2 / G3 | pass / fail / pass | pass / fail / pass | pass / pass / pass | fail / fail / pass |
| G4 deep pattern / G5 lapse | fail / pass (0.55; deep 0.81, shallow 1.15) | fail / pass | fail / fail (1.24 = shallow) | fail / fail |
| date permutation p | 0.33 | 0.22 | 0.08 | 0.15 |
| Q3 composite (%) | −0.000 ± 0.040 | −0.003 ± 0.038 | +0.012 ± 0.021 | +0.011 ± 0.021 |
| SVI v2 column (76 fam, n ≥ 50): Q1 | 0.175 ± 0.052 | 0.104 ± 0.035 | — | — |

### Task 2: Q3 estimators (`r3q3_summary.csv`, `r3q3_pulses.csv`, `r3q3_dD.csv`)

| region | estimator | cells | σ per major (median) | majors ≤ 0.1 % | composite (%) | deep pulse | shallow leak | control / leave-one-site-out | adopt |
|---|---|---|---|---|---|---|---|---|---|
| SVI | B0 raw, P | 3,531 | 0.091 | 12/15 | −0.000 ± 0.040 | 1.14 | 0.45 | pass / pass | — |
| SVI | B0 n ≥ 50 | 1,358 | 0.117 | 0/15 | +0.047 ± 0.049 | 1.15 | 0.40 | pass / pass | — |
| SVI | B1 cc-weighted | 3,774 | 0.082 | 12/15 | −0.025 ± 0.037 | 1.15 | 0.44 | pass / pass | no |
| SVI | **B2 cc-graph** | 1,422 | 0.085 | 11/15 | +0.007 ± 0.035 | **0.98** | 0.38 | pass (1 major) / pass | **yes** |
| SVI | B2p P-graph | 1,369 | 0.094 | 11/15 | +0.042 ± 0.040 | 0.97 | 0.35 | pass (1 major) / pass | yes |
| WA | B0 raw, P | 2,466 | 0.131 | 3/32 | +0.012 ± 0.021 | 1.05 | 0.60 | pass / pass | — |
| WA | B0 n ≥ 50 | 603 | 0.199 | 0/24 | +0.004 ± 0.022 | 1.04 | 0.53 | pass / pass | — |
| WA | B1 cc-weighted | 2,786 | 0.122 | **8/32** | −0.020 ± 0.021 | 1.06 | 0.57 | pass / pass | yes |
| WA | **B2 cc-graph** | 700 | 0.146 | 0/24 | −0.014 ± 0.017 | **0.91** | 0.47 | pass (7 majors) / pass | **yes** |
| WA | B2p P-graph | 630 | 0.171 | 0/24 | −0.010 ± 0.019 | 0.90 | 0.45 | pass / pass | yes |

- The prior test fails for every estimator (0.26–0.97σ).
- Scatter of majors about zero (χ²): SVI B0 1.17, B1 1.99, B2 2.19; WA B0 2.49, B1 3.21, B2 1.25.
- Majors beyond 3σ: SVI B2 1 (2019-11); WA B1 2 (2007-01 −0.70 ± 0.18, 2014-07 +0.53 ± 0.14).

### Task 3: zones (`zones_forecast.csv`, `zones_result.csv`, `zones_gates.csv`, `zones_injection.csv`, `zones_loso.csv`, `zones_peryear.csv`)

| region | zones (along-strike bounds, km) | forecast σ of composite | data composites (%) | equal-zones χ² p | control / leave-one-site-out | zone pulse recovery / cross-leak / shallow leak into contrast |
|---|---|---|---|---|---|---|
| SVI (3 failed: SE zone 0.17) | 2 (−12.7) | 0.100, 0.068 | −0.047 ± 0.100, +0.052 ± 0.068 | 0.45 | pass / pass | 0.86–0.91 / ≤ 0.03 / 0.05 |
| WA-NW | 3 (−22.0, +16.9) | 0.081, 0.053, 0.094 | −0.074 ± 0.081, +0.026 ± 0.053, +0.022 ± 0.094 | 0.61 | pass / pass (contrast signs flip at 5 sites; contrasts ≈ 0) | **0.23–0.32 (fail; design flaw)** / ≤ 0.03 / ≤ 0.15 |

## What failed, why, and the next adaptation
1. **Kernel robustness (task 1).**
   - Obstacle: the deep-step size scales with the assumed deep sensitivity (±30 % across kernels). Strong near-surface
     scattering, if not modelled, moves a shallow yearly cycle into the deep term.
   - The data strongly prefer strong near-surface scattering (two3_20).
   - Next: (a) a selection-recovery test — synthetic data sets made with each kernel, to check that the joint-lapse χ² picks
     the generating kernel. If it does, the lapse-chosen in-sample test applies (deep step 0.93 / 0.81, leak 0.01 / 0.03, both
     pass). (b) A finer family of near-surface variants around two3_20 (top scattering length 1–5 km, layer 1–5 km). (c) Report
     every deep bound inflated by the worst kernel recovery, rather than demanding robustness to physically disfavoured kernels.
2. **Shallow leak in Q3 (task 2).**
   - Obstacle: a structural property of ΔD under a smooth shallow model. It is unchanged by kernels and estimators (0.35–0.60).
   - Next: combine the B2 graph with the distant-station agent's regional-shallow-per-episode term (leak 0.07). Also use zone
     contrasts, where a common shallow change cancels.
3. **SVI 2019-11 candidate (+0.53 ± 0.15 %, B2).**
   - Next:
     - per-instrument graph values for 2019-11 at B009/B010/B011/B926;
     - repeat with the round-2 clock-corrected S picks;
     - check whether the change sits in the deep-weighted families (high a_d) and not the low-a_d ones;
     - check WA's closest majors (2019-08, B0 +0.21 ± 0.13).
4. **WA zone-pulse test (task 3).**
   - Obstacle: the every-major pulse cancels inside ΔD because WA's neighbouring episodes are mostly majors.
   - Next: one-major-at-a-time zone pulses (as the Q3 rule does), or a zone ΔD referenced to inter-ETS bursts only. The 3-zone
     SVI split needs either more south-east stations or a forecast with the distant stations.

## Files
- Scripts `svi_science/scripts/`:
  - `s11_r3_common.py`, `s11_r3_fit.py`, `s11_r3_inject.py`;
  - `s11_r3_kernel.py`, `s11_r3_measure.py`, `s11_r3_kerr.py`, `s11_r3_q3.py`, `s11_r3_zones.py`, `s11_r3_figures.py`;
  - queue `s11_r3_queue.sh`, kernel chain `s11_r3_kernel_chain.sh`.
- Kernel checks:
  - the ell40 re-run reproduces the stored kernels exactly;
  - along-strike bins sum to a_d within 3e-15;
  - new variants: |Σa − 1| ≤ 2e-14, and a_d < 0 for 0 % of pairs.
- Measurement: 767 SVI and 860 WA causal series. Store-rebuilt values reproduce the measurement within 1e-8 %.
- Code-test outputs (not results): `svi/tables/test/`, `*_test.*`.
