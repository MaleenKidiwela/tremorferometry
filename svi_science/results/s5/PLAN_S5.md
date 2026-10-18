# S5 — tidal sensitivity of LFE families, 2003–2026 (pre-registered)

Written 2026-10-07 2:05 pm PT, BEFORE any modulation number was computed. Follows PLAN.md Revision 1 point 6 and
MERLIN_BIG_PICTURE_2026-10-07.md (ranking item 2, prescription step 5). Any later change is appended at the end with
time and reason; nothing above that line is edited after computing.

## 0. What Royer et al. (2015) found (the positive-control target)
Royer, Thomas & Bostock (2015), JGR 120, 384–405, doi:10.1002/2014JB011430. Full text not accessible here; the
abstract was read from Crossref (verbatim in RESULTS). From the abstract:
- LFE families (SVI + Washington, ETS 2003–2013) are sensitive to tidal **up-dip shear stress (UDSS)**.
- Sensitivity to UDSS "rises dramatically from near zero on the first day of strong activity to a maximum ~4 days later".
- The peak LFE rate moves from a correlation with peak tidal **shear-stress rate** (early) to peak **shear stress** (late).
- RTRs start a few days after the main front; 72 % occur at large positive UDSS.
Secondary source (Zhou et al. 2026, arXiv 2602.06703, Table 4 and text): Royer 2015 tidal shear ~2 kPa,
aσ ≈ 7.6 kPa (Beeler 2013 model); phase moves from about −90° early to ~0° late.
From memory only (flagged, not used as a pass rule): modulation amplitude of order tens of % (consistent with
exp(2/7.6) − 1 ≈ 30 %).

## 1. Step 1 — tidal stress 2003-01 → 2013-12 at Bostock's families
- All 130 Bostock hypocentres (`e5r_families_2012.csv`, bkey = his family = our fid for the 115 production families).
  Event→family key is the exact (lat, lon, depth) match already made by P3 (`results/p3_bostock/bostock_events_mapped.parquet`
  holds the 115; the 15 non-production families are keyed the same way here, no new mapping logic).
- Same model as `results/s5_setup` (body tide pyTMD IERS2010 plane stress; ocean load EOT20 + inland tide-gauge field,
  Boussinesq half-space, cells to 700 km; Slab2 plane; slip direction 236°). Ocean cell fields are re-used from the
  setup's saved `ocean_cells_EOT20_gauge_{inner,outer}.npz`; geometry for the 15 extra families by the setup's own code.
- Hourly UTC 2003-01-01 → 2014-01-31 (one month of overlap with the setup).
- Acceptance (before any LFE is touched): (a) for the 115 shared families the recomputed ocean constituent stresses
  equal `ocean_constituent_stress.csv` (EOT20/gauge/cat, inner+outer) to < 0.1 %; (b) hourly τ, σn in 2014-01
  equal the setup's `tidal_stress_svi_hourly.parquet` (max |diff| < 1 Pa). Fail → fix before step 2.

## 2. Observables (both eras, identical code)
For family f: hourly total τ_f (236° slip direction, + promotes thrust), σn_f (tension +), τ_updip_f (Royer's UDSS, column).
A harmonic fit per family over each era's hourly span gives the constituent complex amplitudes; the phase of
constituent c at time t is θ_c(t) = V_c(t) + u_c(t) − G_c,f, so that τ_c = A f cos θ_c (θ = 0 at peak thrust-promoting
shear of that constituent; θ = −90° at peak shear-stress RATE).
- **Headline phase statistic (M2):** for a set of events, R = mean exp(iθ_M2). Directional projection C = Re R
  (= mean cos θ_M2). Modulation amplitude a = 2|R − R_surr| (rate ∝ 1 + a cos(θ − θ0)), preferred phase θ0 = arg(R − R_surr),
  where R_surr is the surrogate mean.
- **O1** separately (lunar diurnal; phase shifts 25°/day relative to the solar day, so it is testable with whole-day
  shifts). **K1** and **S2** are locked to the solar day (K1 drifts 1°/day, S2 0°/day): whole-day surrogates cannot
  randomise them, so they are reported descriptively only, flagged as confounded with time-of-day detection.
- **Stress-level statistic (total τ):** Δτ = mean τ(events) − mean τ(surrogate events), in Pa, and the implied
  exponential sensitivity aσ_hat = Var(τ at surrogate times)/Δτ (rate ∝ exp(τ/aσ) ⇒ Δτ ≈ Var/aσ), comparable to Royer's
  ~7.6 kPa. Columns: σn, Coulomb μ′ = 0.4, τ_updip, body-only τ, ocean-only τ.
- **Rate-vs-stress curve (figure):** events per τ decile / expected per decile (expected from the surrogate events).

## 3. Expected distribution and significance (bursts respected)
- **Active windows:** a family's active days = UTC days with ≥ 1 of its events in the set analysed (Bostock: inside
  his 10 ETS windows; ours: any day).
- **Surrogate = whole-day shifts:** for each family independently, its active days are randomly permuted and every
  event moves with its whole day to the assigned day, keeping its time of day. Within-day bursts move together; time-of-day
  detection effects are identical in data and surrogate. Families with < 3 active days are kept in the observed statistic
  but get no shift (contribute identically to both → bias toward no effect, conservative). 1,000 surrogates (seed 5).
- One-sided p (expected sign) = (1 + #{surr ≥ obs}) / 1001.
- **Random times** (column): uniform continuous random times on the same active days; their M2 phase must be flat
  (Rayleigh p > 0.01).
- **Scrambled control** (ours only; Bostock has none): series c1 at P ≥ 0.75 counted once (tiny, ~65), P ≥ 0.5, and all c1
  candidates (detector-level check). Pass = no M2 modulation at p < 0.01 two-sided at P ≥ 0.5 and all-P.

## 4. Step 2 — POSITIVE CONTROL on Bostock 2003–2013 (130 families; 115 as a column)
Within-episode time: per family per Bostock ETS window, onset t0 = 10th-percentile time of the family's events in that
window (column: first event). dt = (t − t0) in days. Bins: < 1, 1–2, 2–3, 3–4, 4–6, 6–10, ≥ 10 d (dt < 0 goes in "< 1").
**PASS requires both:**
- P1 (sign): pooled M2 projection C > surrogate with p < 0.01 (one-sided, positive = more LFEs at peak thrust-promoting
  shear, i.e. |θ0| < 90°), AND Δτ > 0 with p < 0.01.
- P2 (rise): C(late) − C(first day) > 0 with p < 0.01, late = 2 ≤ dt < 6 d, first day = dt < 1 d. Test = permute the
  late/first-day labels between whole family-days within each family-episode (1,000 permutations; only family-episodes
  with both kinds of day enter), so the overall modulation is kept and only its timing is broken.
Reported, not pass rules: θ0 per dt bin (Royer: ~−90° early → ~0° late), aσ_hat vs ~7.6 kPa, a in % vs "tens of %",
σn and τ_updip, per-episode C, per-depth-third C.
**FAIL → stop**, report, do not run step 3.

## 5. Step 3 — only if step 2 passes: ours 2014–2026 (115 families)
Set: `results/s0/analysis_table.parquet`, series data, P_same ≥ 0.75, keep_xfam8 (counted once). Columns: P ≥ 0.5, P ≥ 0.85, raw.
- Overall modulation (same statistics, same surrogates).
- **Within-episode rise:** front arrival per family from PNSN tremor for the 10 migrating SVI episodes
  (`results/s3/episodes_kinematics.csv` migrates = True): arrival at the family's 5-km along-strike bin = 10th-percentile
  time of tremor in that bin (≥ 20 epicentres), exactly as `scripts/s3b_deep_early.py`. Column: family's own onset
  (10th pct, as Bostock) for the same 10 episodes and for all major episodes. Same bins, same permutation test.
- **Eras:** Bostock 2003–2013 vs ours 2014–2026 — SHAPE only (detectors differ): compare θ0 and the normalised
  rise curve C(dt)/C(all); amplitude ratio reported but not interpreted as change. Within ours: 2014–2019 vs 2020–2026
  (same detector) — difference in C, a and aσ_hat; two-sided p from the pooled surrogate distribution of the difference
  (difference statistic recomputed on each surrogate).
- **Depth:** families split in thirds of down-dip coordinate (`families.csv` downdip_km, as S3b); deep vs shallow C,
  both eras; also by the 12 location groups (descriptive).
- **Power table:** for each comparison, minimum detectable difference MDD = 2.8 × SE(difference) (80 % power,
  α = 0.05 two-sided), SE from the surrogate SD of each side combined in quadrature; expressed in C and in a (%).
- Change "detected" only at p < 0.01; otherwise stated as "a change larger than MDD would have been seen".
- Flags: episodes 2017-03 and 2025-10 (detection problems) as a drop-column.

## 6. Outputs
`results/s5/` tables + `RESULTS_S5.md`; figures `figures/s5_*.png` (modulation vs tidal stress; sensitivity vs days
since front; per-era comparison); scripts `scripts/s5_*.py`. CPU only, < 30 GB.

---
## Appendix A (added 2026-10-07 2:52 pm PT, AFTER step 2 and the first step-3 numbers were seen) — post-hoc, labelled as such
1. The random-continuous-time M2 phase check (section 3) gave Rayleigh p = 0.0035 for Bostock (R = 0.0046, i.e. a 0.9 %
   artefact, because 24 h is not a whole number of M2 cycles and events sit on few days). This fails the letter of the
   "p > 0.01" check but is 40× smaller than the signal; the surrogates (the comparison basis) have R = 0.002. Kept, reported.
2. Raw era and depth contrasts mix the tidal response with the share of events that fall on the first day of activity
   (no response) vs days 2–6 (strong response). Added `scripts/s5_step3_posthoc.py`: the same contrasts STANDARDISED over
   days since the family's own onset (common weights), same surrogates, same p and MDD definitions. These are post-hoc;
   any p < 0.01 found here is "suggestive", not a pre-registered detection.
3. (Appendix B, 2:53 pm PT) The surrogate SD only measures tidal-sampling noise. Episodes may differ from each other by
   more than that (Bostock per-episode a ranged 28–57 %). Added `scripts/s5_step3_episodes.py`: per-episode excess C,
   heterogeneity chi², between-episode SD, and era contrasts tested by permuting EPISODE labels between halves (exact,
   all splits) with weights 1/(sd² + between²). If episodes are heterogeneous, this episode-level test is the one that
   counts for "has sensitivity changed between eras"; the surrogate-only p is then too optimistic.
