# S10 follow-up — tidal-phase dv/v (A) and sub-daily dv/v inside big episodes (B): results

Rules: `PLAN.md` (written ~10:20 am PT Sun 2026-10-11, before any number; Addendum 1 before its post-hoc number).
Scripts `svi_science/scripts/s10_tide_{quantities,common,prep,phase,phase_summary,subday,figures}.py`, queue
`s10_tide_queue.sh`. Tables `tables/`, figures `figures/`. dv/v = −eps. Provisional: three families, one pass.

## A. Tidal-phase dv/v

**Plain words.**
- dv/v does not change with the tidal cycle at a level we can measure.
- Repeats taken at different phases of the M2 up-dip shear stress (UDSS) give the same coda velocity to within a few
  hundredths of a percent, and that holds for every family and design. The pooled curve for all three families
  (`figures/s10_tide_phase_all_24.png`) is flat: amplitude 0.001 % (D1), 0.023 % (D2), 0.027 % (D3). The
  random-shift control gives the same size or larger (99th percentile 0.023 / 0.043 / 0.045 %).
- Binning by UDSS value (kPa) instead of phase gives the same answer. The pooled slope is −0.005 %/kPa (D1), 0.000
  (D2) and −0.001 (D3).
- **No combination passes the pre-registered rule.** That covers 3 families plus pooled, 3 designs, tidal and fake
  cycles, and the scrambled control.
- The tides do show up in *which* LFEs are caught. Repeats near the peak of the rising stress rate (phase 270–315°)
  are ~10 % bigger than those near phase 180° (within-day log amplitude ±0.05). That is a population effect, not a
  velocity effect. It is consistent with S5's finding that the LFE rate peaks near the stress-rate peak.
- **What this rules out, after scaling by the injection recovery** (99 % upper bounds on the M2 tidal dv/v amplitude,
  pooled families):

  | design | bound | per kPa |
  |---|---|---|
  | D3 | 0.08 % | 0.08 %/kPa |
  | D2 | 0.10 % | 0.11 %/kPa |
  | D1 | 0.07 % | 0.07 %/kPa |

  From the UDSS value bins the bound is 0.04–0.07 %/kPa. Single families are about 1.5–2× looser (0.09–0.21 %).
- **Comparison, from memory, not checked against the papers:** shallow (< 1 km) tidal dv/v is reported at
  ~10⁻⁴–10⁻³ (0.01–0.1 %) for tidal strains ~10⁻⁸, i.e. roughly 0.01 %/kPa.
  - Our bound sits at the top of that range in amplitude and ~4–7× above it per kPa.
  - So a deep source zone as stress-sensitive as shallow rock is not excluded. Anything stronger than the most
    sensitive shallow sites is.

**Why the injection scaling matters.**
- We injected eps = a·cos θ_UDSS into every family-10 window and re-ran everything.
- The fraction of a true tidal pattern that comes back:

  | design | 2–4 s, a = 0.1 % | 2–4 s, a = 0.03 % | 3–5 s, a = 0.1 % |
  |---|---|---|---|
  | D1 within-day | 0.30 | 0.25 | 0.15 |
  | D2 day-matched | 0.62 | 0.55 | 0.54 |
  | D3 pooled | 0.84 | 0.79 | 0.78 |

- D1 loses most of the signal. LFEs of one day bunch in tidal phase, so the day's own all-bin reference already holds
  part of the tidal pattern, and the differential removes it along with the day-to-day changes. D2 has the same
  problem to a lesser degree.
- So D1's small raw amplitudes are not more sensitive than D3's. After scaling, the three designs give similar bounds.
  The scaling uses family 10's pooled recovery for all families. It is mildly amplitude-dependent (D1: 0.30 vs 0.25),
  so the D1 bound would be 0.08 % with the 0.03 % recovery.
- Recovery of the 90° (quadrature) component is ≤ 0.05: the fit does not invent a phase.

**Numbers (2–4 s, M2 UDSS phase, control (i) = 1,000 independent random phase shifts per repeat;
`tables/phase_rules.csv`, `phase_bins.csv`).** A = fitted cosine amplitude of eps (= |dv/v|), unscaled.

| family | D1 A (99th ctrl), p | D2 A (99th), p | D3 A (99th), p | sites passing (D1/D2/D3) |
|---|---|---|---|---|
| 10 | 0.012 % (0.036), 0.61 | 0.010 % (0.070), 0.94 | 0.016 % (0.075), 0.81 | 1 / 0 / 1 |
| 87 | 0.011 % (0.037), 0.69 | 0.030 % (0.078), 0.50 | 0.044 % (0.081), 0.28 | 0 / 0 / 2 |
| 17 | 0.009 % (0.041), 0.80 | 0.063 % (0.083), 0.07 | 0.065 % (0.086), 0.08 | 2 / 0 / 4 |
| pooled | 0.001 % (0.023), 0.995 | 0.023 % (0.043), 0.32 | 0.027 % (0.045), 0.24 | 1 / 0 / 3 |

- The pooled per-bin values (2–4 s, D3) lie within −0.04 … +0.05 %. The control SD of one pooled bin is 0.014 %
  (D1) and 0.028 % (D2/D3), which is the per-bin uncertainty. The plots carry no error bars.
- Family 17 in D2/D3 comes closest (p 0.07–0.08: slower around phase 135–180°, faster around 270°). The 3–5 s window
  does not repeat it (p 0.32–0.92). It also has the strongest population effect of the three (see (iv)).
- The one small p in the whole table, family 10 D1 at 3–5 s (A 0.037 %, p 0.004), has 2–4 s p 0.61. The fake 12.0-h
  cycle produces a value just as small (family 10 D1 3–5 s, p 0.008). With ~100 such numbers in the table, a p near
  0.005 is expected by chance.
- Same-day common-offset control (column): same verdicts; the smallest is family 17 D2 p 0.046.
- Leave-bin-out day reference (D1 column, observed only): pooled amplitude 0.007 % (2–4 s), family values
  0.03–0.05 %. All are within the D1 control range.
- Normal stress: M2 σn phase gives pooled 0.002–0.022 %, p 0.35–0.94. The σn value bins give 0.001–0.005 %/kPa,
  p 0.17–0.64. σn and UDSS rates: nothing passes; the smallest p is σn rate D1 pooled at 0.04 (3–5 s same sign, but p > 0.01
  and only 2 sites).

**Controls.**
- (i) Random phase shift: all observed amplitudes sit inside the control distribution (table above;
  `figures/s10_tide_controls.png`).
- (ii) Fake 12.0 h, 13.5 h and 10.5 h cycles: same size as the tide (pooled 0.002–0.017 %), and none passes. The real
  tide is indistinguishable from a fake cycle, which is consistent with nothing being there.
- (iii) Scrambled control (c2 noise windows, same path): no pass. Its D3 amplitudes are large (0.14–0.25 %) because
  those stacks are noise against the data reference, as expected.
- (iv) Population check, per family, against the same 1,000 draws (`tables/phase_population.csv`):
  - Mean P does not vary with tidal phase (p 0.24 and 0.42 for families 10 and 87). Family 17 sits right at the bar
    (p 0.010; ±0.01 in P).
  - Within-day S amplitude varies with M2 UDSS phase in all three families (p 0.001), by ±0.03–0.06 log10 units,
    biggest at 270–315°. S6 dM shows the same in families 10 and 17 (p 0.001; 87 p 0.015).
  - Waveform similarity to the grand stack (realignment cc) varies by ±0.01 (p ≤ 0.007).
  - Realignment lag, a sub-location proxy, does not vary (p 0.20–0.74). So the bins hold the same sources, but the
    LFEs in them are slightly bigger or smaller depending on tidal phase. The equal dv/v across bins means this
    amplitude/similarity difference did not create a velocity artefact.
- No passes, so the bound in the plain-words section applies.

Figures:
- `figures/s10_tide_phase_all_24.png`, the key figure: pooled families, the three designs, per-station bins plus the
  pooled points and cos fit.
- `figures/s10_tide_phase_fam{10,87,17}_{24,35}.png`: per family, phase row plus UDSS-value row.
- `figures/s10_tide_controls.png`: observed amplitudes against the control spread for tide, σn, the fakes and the
  scrambled control.
- All dv/v points are coloured by CC, marker = station, no error bars.

## B. Sub-daily dv/v inside the biggest episodes

**Plain words.** Within an episode, 6-hour and 12-hour stacks scatter by about ±0.4 % from bin to bin, at every kind of
station. There is no consistent change from the start to the end of an episode: the pooled first-third vs
last-third difference is −0.03 % (6 h) and −0.005 % (12 h), well inside what shuffling the bins in time gives, and
stations split about evenly in sign. A change larger than about 0.2 % (99 %) between the first and last third of
an episode would have been seen when the three families are pooled. The bins with ≥ 20 repeats come almost only
from the two largest episodes (2004-07, 2005-09) and from POLARIS-era stations, and almost none fall before the front,
so the before/after-front comparison cannot be made as pre-registered.

Numbers:
- Bins with ≥ 20 repeats: 518 at 6 h (median 29 repeats, median cc 0.64) and 483 at 12 h (34; 0.70) over all
  stations and the three families; 18–20 % (6 h) and 15–16 % (12 h) give no correlation peak within ±1 % and are
  dropped. By episode (6 h): 2004-07 262, 2005-09 162, 2013-09 34, 2003-02 18, 2011-08 14, 2024-08 12, 2008-05 9,
  2019-08 7, 2023-04 0. By family: 10 → 281, 17 → 207, 87 → 30.
- Scatter of a 6-h bin about its series' episode mean (series with ≥ 2 bins): 0.40 % (boreholes 0.33 %, CN 0.39 %,
  POLARIS 0.41 %). This is the per-bin noise level; the per-bin uncertainty is not drawn on the plots.
- Front: S3 tremor front at the family for 2019-08, 2023-04, 2024-08; elsewhere the family's onset (10th-percentile
  time of its high-P repeats). `tables/subday_fronts.csv`.
- **First vs last third** (`tables/subday_tests.csv`; Δ = last − first, dv/v, CC²-weighted; 2,000 time shuffles):

  | bins | families | units | Δ dv/v 2–4 s | shuffle SD | p | Δ dv/v 3–5 s | 99 % bound |
  |---|---|---|---|---|---|---|---|
  | 6 h | 10 | 28 | +0.054 % | 0.079 | 0.50 | +0.059 % | 0.26 % |
  | 6 h | 17 | 29 | −0.173 % | 0.105 | 0.10 | +0.062 % | 0.44 % |
  | 6 h | pooled | 57 | −0.030 % | 0.063 | 0.65 | +0.061 % | 0.19 % |
  | 12 h | 10 | 47 | +0.026 % | 0.069 | 0.73 | −0.008 % | 0.20 % |
  | 12 h | 17 | 25 | −0.070 % | 0.141 | 0.63 | +0.018 % | 0.43 % |
  | 12 h | pooled | 72 | −0.005 % | 0.063 | 0.97 | 0.000 % | 0.17 % |

  Family 87 has no series with ≥ 3 bins. No row passes (p < 0.01, same sign at 3–5 s, ≥ 2 sites): the 2–4 s and 3–5 s
  windows disagree in sign in most rows, and stations split about evenly (pooled 6 h: 10 of 19 stations share the
  pooled sign; per-station values −0.67 … +0.62 %).
- **Before vs after the front:** 0 units have ≥ 2 bins on each side (Addendum 1). Post-hoc with ≥ 1 bin each side
  (labelled post-hoc, not a pass rule): pooled 6 h −0.038 % (38 units, shuffle SD 0.080, p 0.65), 3–5 s +0.037 %;
  family 10 −0.063 % (p 0.58); family 17 −0.019 % (p 0.87).
- Scrambled control: only 7 bins have ≥ 3 scrambled-control windows; their stacks give |eps| ~0.7 % (noise against
  the reference), as expected.
- Figures: `figures/s10_subday_fam{10,87,17}_{24,35}.png` (one panel per episode, dv/v minus the series' episode
  mean vs time in PT, colour = CC, marker = station, front dashed).
