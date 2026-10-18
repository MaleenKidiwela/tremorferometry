# S11 pilot (SVI): separating deep from shallow dv/v with S-pinned residence kernels

Written Mon 2026-10-12, 1:42 pm PT. Rules: `PLAN.md` (written before any number). Decision log: `NOTES.md`. Design:
`MERLIN_DESIGN.md`. dv/v in %. "Scrambled control" = station data shifted to another time of the same day. All results
are provisional pilot results (76 families × 7 stations, 2003–2026). Steps 1–6 done. Step 7: WA-NW done (section below); POLARIS not started.

**Instrument vetoes (Wed 2026-10-14, ~2:55 pm PT).** B013 2012-07-13 … 2013-07-24 and B009 2007-01-08 … 2009-10-01 are
now dropped (`svi_science/instrument_vetoes.csv`; reasons in `../s14_cc/RESULTS.md`). The tables below are pre-veto.
Refits: SVI Q3 composite −0.002 ± 0.040 %, Q1 0.177, Q2 +0.109; WA-NW composite +0.008 ± 0.020 %, Q1 0.162, Q2 −0.005.
No reading changes.

## Step 7: WA-NW, the same pipeline on an independent region (Mon 2026-10-12, ~4:05 pm PT). Read this first.
Rules: `PLAN.md` §5, written before any WA number. Decisions: `NOTES.md` 16–19. Outputs: `wa_nw/tables`, `wa_nw/logs`.
Figures: `figures/fig_wa_{compare,q3_dD,annual,sites,kerr}.png`. Raw path. Cells: n ≥ 20, cc ≥ 0.35, 1/σ² with control
inflation. "Scrambled control" = station data shifted to another time of the same day. dv/v in %.

**In plain words**
- **The annual cycle replicates in Washington.**
  - The kernel-free cycle of the deep-LFE coda is 0.090 ± 0.007 %, fastest around day 268 (late September). SVI's is
    0.079 ± 0.009 % (day 226).
  - When the smooth model splits it into deep (25–45 km) and shallow (0–5 km), the deep share is 0.164 ± 0.071 %
    (2.3σ), peaking at day 195. SVI gives 0.114 ± 0.037 % (day 229) on the same n ≥ 20 rule. The two agree within errors,
    and their peaks are 34 days apart.
- **Where the cycle sits is still not shown to be deep, and WA leans shallow.**
  - **Lapse test.** WA's cycle grows slightly later in the coda (3.5–6 / 2–3.5 s = 1.24 ± 0.16). That is exactly the
    shallow prediction (1.23); the deep prediction is 0.79. SVI's 0.52 sat nearer the deep prediction.
  - **Station pattern.** Neither the deep nor the shallow kernel explains the per-station amplitudes (both χ² ≈ 153 on
    24 dof).
  - **Kernel dependence.** The deep share moves from 0.08 % to 0.20 % across kernel models.
  - So Q1 is a measured value with a 2σ upper bound of 0.31 %, not a deep seasonal change.
  - It rests on cells with 20–49 repeats. With n ≥ 50 (603 cells) it is 0.02 ± 0.11 %.
- **No decade-scale deep change.** Q2 (2015–25 minus 2007–14) = −0.003 ± 0.089 %, so |Q2| < 0.18 % (2σ). That is half
  of SVI's bound (0.35 %). An injected +0.3 % deep step comes back at 0.83, or at 0.66–0.71 if the true kernel differs,
  so the bound after calibration is about 0.27 %.
- **Slow-slip episodes: no deep-weighted change on average, the tightest bound so far.**
  - Over the 32 WA majors, ΔD averages +0.012 ± 0.021 %, so an average episode changes the deep-weighted velocity by less
    than 0.05 % (2σ).
  - The 13 majors co-timed with SVI majors give −0.026 ± 0.039 %. SVI's 15 majors gave −0.000 ± 0.041 %.
  - Single majors have σ 0.09–0.19 % (median 0.13 %). Only 3 have σ ≤ 0.1 %, so the per-episode precision target
    (≥ 5 majors) is missed.
- **Single episodes scatter more than their errors in WA, but that scatter cannot be called deep.**
  - Across the 32 majors, χ² = 2.5, and 9 lie beyond 2σ. The largest are 2014-08 +0.60 ± 0.19 % (3.1σ), 2014-07
    +0.41 ± 0.14 % (2.9σ) and 2007-01 −0.51 ± 0.19 % (−2.6σ).
  - The scrambled control scatters as its errors predict (χ² 1.13), so the excess is in the data, not the noise model.
  - Two tests stop a deep reading. A common shallow change in one episode leaks into ΔD at 0.60 × (rule ≤ 0.35). Doubling
    or halving the prior moves single-episode ΔD by up to 0.84σ (rule 0.1σ).
  - No WA episode passes the wording rule. The excess scatter is most likely episode-locked shallow or site change, or
    model-conditional.
- **Confirming tests that pass.**
  - **Scrambled control:** 29/32 majors within 2σ; composite +0.004 ± 0.048 %.
  - **Leave-one-site-out:** passes for every one of the ~50 sites. The heaviest site, B013+B943, holds 37 % of the weight;
    without it the composite is +0.034 ± 0.026 %.
  - **Odd/even family halves:** agree on 97 % of majors; composites −0.011 ± 0.052 and +0.021 ± 0.022 %.
  - **Deep pulse:** an injected deep pulse comes back at 1.06 × (rule 0.8–1.1; SVI 1.14).
- **No single station hijacks the result.** This was checked because a bursty station (W060) once hijacked a WA result.
  - No site exceeds 50 % of the weight (B013+B943 37 %, LRIV group 20 %, B001 19 %).
  - DOSE has an outlying annual amplitude (0.46 %, robust z 4.1). It carries only 39 cells (2.9 % of the weight) and is
    kept as a flagged column.
- **The two regions' year-to-year deep-weighted values do not track each other** (correlation 0.08 over 20 calendar
  years).
- **A wrong kernel matters more in WA. One adaptation passes; it is post hoc.**
  - The pre-registered coherent kernel-error injections fail in WA. Data generated with the two-layer kernel (two5_40)
    give deep-step recovery 0.66 and a shallow-to-deep leak of 0.41. With ell80 the leak is 0.27. At SVI the same tests
    passed.
  - Re-inverting the same waveform injections (`s11_wa_kerr.py`, `fig_wa_kerr.png`):
    - **M2, kernel-averaged coefficients (mean over the 4 kernel models), passes all six cases**: leak 0.13–0.18,
      recovery 0.73–0.85. Under M2, Q1 = 0.152 ± 0.081 % (1.9σ, day 194) and Q2 = −0.068 ± 0.097 %.
    - **M3, kernel-free shallow term,** removes the leak (≤ 0.10) but under-recovers a deep step (0.69 with two5_40).
    - **M1, one shallow coefficient per site,** halves the leak (0.24) but has the same recovery problem.
  - **Caveat:** M2 is tested partly with kernels that are inside its own average.
  - **Next adaptation to try:**
    - (i) Leave-one-kernel-out: generate with a variant outside the average (e.g. a 3 km/20 km two-layer, ell 10 km)
      and invert with the average of the rest.
    - (ii) Choose the kernel from the data, using the lapse dependence (joint-lapse χ² per variant) and the per-station
      amplitudes.
    - (iii) Run M2 on SVI for a like-for-like Q1.
- **WA-SW could not be measured.** Of 311 SW repeats at P ≥ 0.75 on tremor days, only 51 lie on active family-days,
  which gives no measurable series.

| quantity | SVI headline | SVI, n ≥ 20 (115 fam) | WA-NW, n ≥ 20 |
|---|---|---|---|
| Q1 deep annual (smooth) | 0.174 ± 0.053 %, day 225 | 0.114 ± 0.037 %, day 229 | 0.164 ± 0.071 %, day 195 (M2: 0.152 ± 0.081) |
| kernel-free pooled annual | 0.079 ± 0.009 %, day 226 | 0.063 ± 0.006 %, day 231 | 0.090 ± 0.007 %, day 268 |
| Q2 decade contrast | +0.12 ± 0.11 % (2015–25 − 2004–14) | +0.06 ± 0.07 % | −0.003 ± 0.089 % (2015–25 − 2007–14) |
| Q3 composite, all majors | −0.000 ± 0.041 % (15) | — | +0.012 ± 0.021 % (32) |
| Q3 composite, co-timed majors | (the same 15) | — | −0.026 ± 0.039 % (13) |
| σ_ΔD per major (median) | 0.091 % (12/15 ≤ 0.1) | — | 0.131 % (3/32 ≤ 0.1) |
| lapse ratio 3.5–6 / 2–3.5 s | 0.52 ± 0.13 (deep 0.80, shallow 1.20) | — | 1.24 ± 0.16 (deep 0.79, shallow 1.23) |

| WA-NW gate | rule | result |
|---|---|---|
| kernels | \|Σa − 1\| ≤ 0.02; a_d < 0 ≤ 10 % | **pass** (1.7e-14; 0 %); 684 pairs; a_d 0.49, a_s 0.29; ell spread of a_s ×2.1 |
| Q1 G1 leave-one-site-out / G2 family split / G3 scrambled control | as §3 | **pass** / **pass** / **pass** (0.22 ± 0.14 %); date permutation p 0.08 |
| Q1 G4 station pattern / G5 lapse | as §3 | **fail** (neither fits) / **fail** (lapse = shallow) |
| A deep step (ell40 / two5_40 / ell80 generated) | Q2 rec 0.7–1.3, Q1 leak < 0.03 % | 0.83 / **0.66** / 0.71; leak 0.022–0.026 % |
| B shallow annual (same three) | Q1 leak < 0.25 | 0.10 / **0.41** / **0.27** |
| prior ×2/×0.5 | < 0.1σ | Q1 **0.32σ/0.40σ**; Q2 ≤ 0.10σ; Q3 majors **0.62σ/0.84σ** (composite 0.03σ/0.15σ) |
| Q3 acceptance | σ ≤ 0.1 % for ≥ 5 majors | **fail** (3) |
| Q3 scrambled control | ≥ 90 % within 2σ; composite within 2σ | **pass** (29/32; 0.004 ± 0.048) |
| Q3 deep pulse / shallow pulse | 0.8–1.1 / ≤ 0.35 | **pass** 1.06 (0.98–1.15) / **fail** 0.60 (all-majors composite 0.22) |
| Q3 leave-one-site-out / family split | ≥ 90 % within 2σ | **pass** (all sites) / **pass** (97 %) |
| site dominance | no site > 50 % of Σ1/σ² | **pass** (max 37 %, B013+B943) |

**Scale of the run.** 108 families (28 Slab2 depths), 25,400 repeats on 2,336 active family-days, 90 stations in 76
sites. 1,006,851 windows (601,647 data + 405,204 scrambled control; 11.6 GB) were read from 16,041 raw station-days.
2,466 headline cells in 288 series and 47 episodes. Per-episode D (episode model): 2σ bound 0.52 %.

---
## Revision 2 (Q3, Mon 2026-10-12 ~2:25 pm PT) — read this first
Rules: `PLAN.md` §4, written before these numbers. Source: `MERLIN_FORECAST.md`. Raw path. Script `s11_q3.py`.
Tables `tables/q3_*.csv`. Figure `figures/fig_q3_dD.png`.

**Correction to v1 / Revision 1.** The per-episode deep σ of ~0.4 % came mostly from the common level of D, which
the data cannot fix. Splitting one episode's change into deep and shallow (attribution) is not identifiable at any
cell count (σ > 1 %, 14 % for the 15-major composite), because the series offsets absorb a_d. Merlin's numbers are
reproduced with our code: 1.76 % per major, 14.05 % composite, at prior 1000.

**In plain words**
- **What we can measure:** ΔD(e), how much faster or slower the deep-weighted part of the coda was in one slow-slip
  episode compared with its neighbours within ±3 years. The shallow part is held to a smooth seasonal-plus-slow curve
  per site. With all 115 SVI families and cells of ≥ 20 repeats, the precision target is reached: σ ≈ 0.09 % per
  major episode (0.07–0.12 %). 12 of the 15 majors have σ ≤ 0.1 %. The 15-major average is known to ± 0.04 %.
- **What the data show: no episode stands out.** Every major's ΔD is within 2σ of zero. The largest are
  2004-07 +0.18 %, 2008-05 −0.19 % and 2019-11 +0.19 %, each about 2σ. Across the 15 majors the scatter is what the
  errors predict (χ² 1.17). The 15-major average is −0.000 ± 0.041 %. So a slow-slip episode does not change the
  deep-weighted velocity by more than about 0.15–0.23 % (single episode, 2σ), or 0.08 % on average.
- **Confirming tests that pass.** The scrambled control gives zero (all 13 measurable majors within 2σ). Leaving out
  any site, PGC included, gives the same answers (15/15 within 2σ). So ΔD is not a station or instrument artefact.
- **Limits (tests that fail).**
  - A common shallow change in one episode leaks into ΔD at 0.45 × (rule ≤ 0.35).
  - An injected deep pulse comes back at 1.14 × (rule 0.8–1.1). The leave-episode-out reference makes the raw path
    over-read slightly.
  - Doubling or halving the prior moves ΔD by up to 0.5σ (rule 0.1σ).

  ΔD is therefore a deep-weighted measure: deep change plus about half of any shallow change common to all sites in
  that episode. It is calibrated to ~15 % and its values are model-conditional below 0.5σ.
- **The 39 missing families added little.** 37 have repeats, mostly 2003–2006: +84 cells. The jump in precision came
  from using cells with ≥ 20 repeats and from restoring 32 B010 cells. Those cells had been silently given zero
  weight by a control-inflation overflow, now fixed.

| configuration | cells | σ_ΔD majors (median, range) | majors ≤ 0.1 % | 15-major composite |
|---|---|---|---|---|
| (a) 76 fam, n ≥ 50 | 1,304 | 0.114 (0.098–0.200) | 1 | +0.072 ± 0.048 % |
| (b) 76 fam, n ≥ 20, B010 fixed | 3,447 | 0.091 (0.073–0.116) | 12 | +0.017 ± 0.041 % |
| (c) 115 fam, n ≥ 20 (headline) | 3,531 | 0.091 (0.073–0.116) | 12 | −0.000 ± 0.041 % |

| gate (headline c) | rule | result |
|---|---|---|
| acceptance | σ ≤ 0.1 % for ≥ 5 majors | **pass** (12) |
| scrambled control | ≥ 90 % within 2σ, composite within 2σ | **pass** (13/13; 0.07 ± 0.15 %) |
| leave-one-site-out (PGC, B926, VGZ, SNB) | ≥ 90 % within 2σ | **pass** (15/15 each; without PGC composite +0.03 ± 0.07 %) |
| deep pulse +0.3 % in one major | median recovery 0.8–1.1 | **fail** 1.14 (1.10–1.55) |
| common shallow pulse 0.1 % in one major | median leak ≤ 0.35 | **fail** 0.45 (0.28–0.62); all-15 composite 0.57 |
| prior ×2 / ×0.5 | < 0.1σ | **fail** (≤ 0.27σ / 0.49σ) |

---
## Revision 1 (v2, Mon 2026-10-12 ~2:10 pm PT)
Rules: `PLAN.md` §3, written before these numbers. Source: `MERLIN_REVIEW_1.md`. No new measurement. Raw path.
v2 cells: n ≥ 50, cc ≥ 0.35, 1/σ² weights → 1,304 cells, 250 series, 50 episodes. Script `s11_v2.py`. Tables
`tables/v2_*.csv`. Figure `figures/fig_v2_annual.png`.

**In plain words**
- **The dv/v from deep LFE coda has a clear annual cycle.** The cycle is fastest in mid-August. When the smooth model
  splits it between a deep layer (25–45 km) and shallow layers (0–5 km) at each site, the deep share comes out at
  **Q1 = 0.17 ± 0.05 %**. The check tests do not confirm that this part is really deep.
  - **Two family halves.** Odd and even families give the same size of deep term but peak 100 days apart (G2 fails).
  - **Shuffling episode dates.** A deep term this large appears 15 % of the time (p = 0.15).
  - **Station pattern.** Under a deep cause, the stations that see the deep layer most (SNB, B926) should show the
    biggest cycle. They show the smallest. PGC and VGZ see the deep layer least and show the largest. This fits a
    shallow cause better (G4: deep pattern p 0.001, shallow p 0.07).
  - **Kernel model.** The deep term changes from 0.08 to 0.17 % with the choice of kernel model.
  - **Leak from a deep step.** A deep step in 2015 leaks 0.06 % into the deep annual term, because the slow-slip
    episodes recur at a regular time of year.

  So Q1 is reported as a measurement with an upper bound (≤ 0.28 %, 2σ), not as a deep seasonal change.
- **One test points the other way.** The cycle is weaker later in the coda (3.5–6 s / 2–3.5 s = 0.52 ± 0.13). A
  deep source predicts 0.80 and a shallow one 1.20, so this sits closer to deep. But it is steeper than either
  prediction, so the kernels do not fully explain it.
- **No decade-scale change in the deep layer.** The 2015–2025 minus 2004–2014 contrast is **Q2 = +0.12 ± 0.11 %**, so
  any deep decade change is smaller than about 0.35 % (2σ). A +0.3 % deep step inserted into the waveforms comes back
  at 0.81–1.06 of its size, including when the data are generated with another kernel model, so this bound is
  calibrated.
- **Single-episode deep values are still not resolved** (2σ bound 0.82 %). Revision 2 corrects the reason: that σ was mostly the unidentifiable common level, and per-episode attribution is unidentified at any cell count.
- **SNB tension, investigated.** SNB records the deep layer most strongly but has almost no cycle in 2003–2006
  (0.016 ± 0.023 %; 137 of its 198 cells). Its own cycle peaks two months earlier (day 168). Leaving SNB out raises
  Q1 to 0.20 ± 0.08 %. SNB behaves like a site with its own weak shallow cycle, not like a deep signal.

| quantity / gate | rule | result |
|---|---|---|
| Q1 deep annual (smooth) | — | 0.174 ± 0.053 %, peak day 225; slow+annual 0.164; joint-lapse 0.147; ell20/ell80/two5_40/seed2 0.148/0.075/0.091/0.110; cc ≥ 0.8 0.05 ± 0.09 |
| Q2 decade contrast | — | +0.12 ± 0.11 % (joint-lapse −0.05 ± 0.09) |
| prior ×2 / ×0.5 | \|Δ\| < 0.1σ | Q1 +0.08σ / −0.20σ (flag); Q2 ≤ 0.08σ (pass) |
| G1 leave-one-site-out | within 2σ, peak ±60 d | pass (0.15–0.24 %) |
| G2 family split | within 2σ, peak ±60 d | **fail** (peaks 273 vs 172 d) |
| G3 scrambled control | Q1 < 2σ | pass (0.39 ± 0.35 %, 106 cells); date permutation p 0.15 |
| G4 station pattern | ∝a_d p > 0.05 and Δχ² > 4 | **fail** (∝a_d χ² 22.8, p 0.001; ∝a_s 11.5, p 0.07) |
| G5 lapse ratio | ≥ 1σ closer to deep | pass (0.52 ± 0.13; z −2.2 deep, −5.3 shallow) |
| A deep step (ell40 / two5_40 / ell80 generated) | Q2 rec 0.7–1.3, Q1 leak < 0.03 % | Q2 0.89 / 1.06 / 0.81; Q1 leak 0.06–0.07 % → **fail** on the leak clause |
| B shallow annual (same three) | Q1 leak < 0.25 | 0.06 / 0.04 / 0.18 → pass |

**Next (suggestion).** Per-station shallow amplitudes need a lever that does not depend on the deep coefficient. The
lapse result (steeper than either kernel) should be checked against the straight-ray approximation (refraction in
the top 3 km), before Q1 is worth replicating with WA or POLARIS.

---
## Pilot v1 (first run, superseded where it conflicts with Revision 1)
### What the data show (v1)
- **The kernels work and are well behaved.** In the 2–5 s window after S, a dv/v change in the 25–45 km layer counts
  with weight a_d ≈ 0.55 (10–90 %: 0.43–0.83). A change in the top 5 km counts with weight a_s ≈ 0.33. The 5–25 km
  layer counts with ≈ 0.08. The unit test passes exactly: a uniform change gives the right answer (max |Σa − 1| = 2e-14).
  No pair has a negative deep weight. The deep weight changes ×4–6 between shallow and deep families at each site, and
  that spread is what can separate deep from shallow. (`figures/fig_kernels.png`)
- **The raw episode stacks are the calibrated data path.** A step inserted into the waveforms comes back at
  1.14 × its size on the raw path. On the user's SVD-Wiener day-stack path it comes back at 0.26 × (median over ~500
  series). By the rule fixed in the plan, the raw path is the headline. The SVD path is reported beside it.
- **A slow deep change can be recovered.** A +0.3 % deep step from 2015, inserted into the waveforms and run through
  the whole chain, comes back as 0.91 × in D. At most 0.21 of it leaks into any site's shallow term (0.98 and 0.15 for
  the post-hoc column without the cc cut). Single-site shallow steps leak ≤ 0.19 into D. A ×/÷1.5 kernel error leaks
  0.28. A common shallow annual cycle leaks 0.30 into D, just above the 0.25 bar. On the SVD path the same deep step
  comes back at only 0.32 ×.
- **No long-term deep change is seen.** The difference in D after vs before 2015 is +0.04 ± 0.12 % on the raw headline,
  −0.01 ± 0.09 % on raw without the cc cut, and −0.05 ± 0.07 % on SVD (≈ −0.15 ± 0.2 % after correcting for its 0.32
  recovery).
- **The deep value for each single episode is not resolved.** Its σ is 0.4 % on the raw path and 0.25 % on SVD. It
  moves by 34–66 % when the prior is doubled or halved (the rule is < 10 %), so by the plan it is reported as
  unresolved. The reason is the expected trade-off between the deep term and the mean of the site shallow terms
  (correlation −0.9 in the weak eigen-directions). Because of this, the split-network, per-site-sign and lapse gates
  (d, e, f) fail. They test the per-episode D, and that D is mostly noise and prior here.
- **The annual cycle.** Without kernels, the cells carry a 0.07–0.10 % annual cycle. The co-located PGC-site instruments
  (surface vault and three boreholes at 193–218 m) see the same amplitude within errors for 3C, Z and H (χ² p 0.23–0.88).
  After the inversion, the deep annual amplitude is 0.07 ± 0.09 % on the raw headline (gate g passes). It is
  0.13 ± 0.06 % on the post-hoc column without the cc cut, and part of that is leakage (injection B). The SVD path gives
  0.02 ± 0.05 %. Shallow S at PGC is 0.15 ± 0.14 %. The pilot does not yet decide which layer carries the cycle.

## Gates
| gate | rule | raw (headline) | SVD |
|---|---|---|---|
| unit test (step 3) | \|Σa−1\| ≤ 0.02, a_d<0 ≤ 10 % | pass (2e-14, 0 %) | same kernels |
| step-4 precision | borehole σ(n≥100) ≤ 0.25 % | pass (0.13–0.15 %) | — |
| (a) co-located | report | equal within errors (3C p 0.31) | — |
| (b) A deep step | recovery 0.7–1.3, S leak < 0.3 | pass (0.91, 0.21) | fail (0.32) |
| (b) B shallow annual | D leak < 0.25 | fail (0.30) | pass (0.05) |
| (b) C one-site step | D leak < 0.2 | pass (≤ 0.19) | pass |
| (b) D kernel error | D leak < 0.3 | pass (0.28) | pass |
| (c) scrambled control | ≥ 90 % within 2σ | not measurable (control cells railed) | pass (15/15) |
| (d) split network / families | slope 0.7–1.3, r > 0.5 | fail (0.15, r 0.07 / 0.19, r 0.29) | fail |
| (e) site sign / cell ρ | ≥ 70 % / \|ρ\| < 0.2 | fail 56 % / pass (≤ 0.09) | fail 49 % / pass |
| (f) lapse | slope 0.7–1.3 | fail (−0.20) | fail (−0.22) |
| (g) seasonal D | amp < 2σ | pass (0.07 ± 0.09) | pass (0.02 ± 0.05) |
| D4 prior test | Δ < 10 % | fail → per-episode D unresolved | fail |

## Main limitations
- The pre-registered raw cell rule (n ≥ 50 and cc3c ≥ 0.8) keeps only 21 % of the n ≥ 50 cells. That leaves 257 cells
  over 38 episodes, about 7 per episode, where the design expected 100–250. Dropping the cc cut (post hoc) gives
  1,345 cells, but the per-episode D is still unresolved (σ 0.41 %).
- The kernels use straight legs: no refraction in the low-velocity top layers. a_s changes ×1.9 across the
  mean-free-path variants.
- Only four sites are used, and three of the PGC-site instruments sit at the same place.

## Files
- Tables:
  - `kernel_coef.parquet`, `kernel_unit_test.csv`, `kernel_summary.csv`, `kernel_mc.npz`;
  - `episode_raw.parquet` (raw cells, all columns), `episode_svd*.parquet`, `series.csv`;
  - `recovery.csv`, `headline.json`;
  - `inv_<path>_<run>.csv` and `cells_<path>_<run>.csv`, `inv_summary*.csv`;
  - `inject_summary*.csv`, `inject_deltas*.csv`;
  - `gates.csv`, `gate_a.csv`, `sensor_depths.csv`;
  - `store/` (day stacks).
- Figures: `fig_kernels`, `fig_D_S_series`, `fig_cells_by_site`, `fig_seasonal`, `fig_injection`, `fig_gates_scatter`.
- Scripts: `svi_science/scripts/s11_{common,gate_a,kernel_mc,measure,paths,svdpath,invert,inject,gates,figures}.py`,
  `s11_queue.sh`.
