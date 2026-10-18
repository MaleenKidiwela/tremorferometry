# Q7 PLAN — borehole strain stack on deep between-episode LFE bursts (pre-registered)

Written 2026-10-07 4:55 pm PT, before any strain number was computed. Merlin big-picture step 6
(`svi_science/MERLIN_BIG_PICTURE_2026-10-07.md`). Question in plain words: are the deep between-episode LFE bursts small
slow-slip events? Borehole strainmeters see slip ~1.5 magnitude units smaller than GNSS (Q5 GNSS bound Mw ≤ 5.9–6.0).

## Data (deviation from the brief, stated up front)
- Level-2 processed strain at `data.earthscope.org/archive/strain/bsm/processed/` (and the GAGE archive) **redirects to
  a login page**. Per the user rule (no accounts, no credentials) it was not used.
- Instead: the **public** EarthScope FDSN dataselect service, PB network, channels RS1–RS4 (the four raw gauge readings,
  10-min samples, counts) and RDO (barometric pressure, 30-min). ~1 MB per station-year. Script `q7_download.py`,
  manifest `data/bsm/MANIFEST.csv`, files `data/bsm/raw/<STA>/<STA>.<YEAR>.mseed`. 22 PB strainmeters 46.3–49.3 N,
  2009-01-01 to 2026-10-06.
- So we do our own Level-2-like processing (below). Gauge orientations come from the FDSN channel azimuths; stations whose
  four azimuths are all 0 (unknown) are dropped.

## Processing (per station, per gauge)
1. Counts → linearised gauge strain: e = (0.2 mm / 87 mm) · [r/(1−r)], r = counts/1e8 (GTSM gap/diameter; UNAVCO
   formula). Hourly means (≥ 4 of 6 samples).
2. Offsets: hourly first differences > 20 × MAD (of that station-year's hourly differences) and > 5 nanostrain are
   treated as instrument steps and removed (difference set to the local median). A 1–25 h slow event changes strain
   by < 1 nanostrain per hour at Mw 5.5 and > 30 km, so this cannot remove the signal sought (checked by injection).
3. Tides and air pressure: per station-year, least squares on the high-passed hourly series (minus a 3-day running
   median) of 8 tidal constituents (M2 S2 N2 K2 K1 O1 P1 Q1) + pressure admittance (high-passed RDO). The fitted tides
   and pressure term are subtracted from the full series. Years with < 60 days of data are dropped.
4. Gauge → model: nominal isotropic GTSM response, e_i = ½[A·εa + B·(−γ1 cos2φ_i + γ2 sin2φ_i)], εa = εEE+εNN,
   γ1 = εEE−εNN, γ2 = 2εEN, φ_i = gauge azimuth, A = 1.5, B = 3 (nominal coupling; real tidal calibrations differ by up to
   ~2× and show cross-coupling — the detection test does not depend on this, only the Mw scale does; see check C below).

## Step estimator (one event, one gauge)
- Hourly corrected series y(t). Event [t0, t1] (burst first/last LFE, UTC). Half-window W = 2 days (primary).
- Fit y = a + b·(t − tm) + s·H(t > t1) on hours in [t0 − W, t0] ∪ [t1, t1 + W]; the burst interval itself is excluded.
  Each side needs ≥ 70 % of its hours. s = strain change across the burst.
- Secondary (descriptive): W = 1 day and W = 5 days.

## Expected pattern and combination (matched filter)
- Template g: gauge steps per 1 N·m of pure thrust slip on the Slab2-based fault mesh (`cascadia-sse-strain`
  `data/cascadia_mesh.txt`, read-only) on all patches within 10 km of the burst's interface point (as Q5), uniform
  slip, μ = 40 GPa, cutde half-space strain at 200 m depth below each station, ν = 0.25. Thrust sign pinned by the
  same check as `slip_build.sign_check` (+dip-slip on this mesh = normal sense).
- Per station-year, the 4 × 4 covariance C of gauge steps is estimated from the random-day windows of that year
  (robust: steps clipped at 5 × MAD). Per event: m̂ = Σ_s gᵀC⁻¹s / Σ_s gᵀC⁻¹g (in N·m), σ_m = (Σ_s gᵀC⁻¹g)^−½;
  stations within 200 km with valid data.
- Stack over bursts: inverse-variance weighted mean of m̂ (N·m per burst) → Mw.

## Event sets
- **BE bursts** (no PNSN tremor within 30 km that day) from `results/q5_bursts/bursts.csv`: SVI primary (40, 2014–26),
  WA primary (62, 2015–26), WA E1/E2 secondary (41, 2009–14). **Primary test set = SVI + WA primary.** Subsets
  (descriptive): big bursts (≥ 3 LFEs), with / without PNSN tremor ≤ 100 km ± 6 h.
- **Broad bursts** (between episodes, local tremor allowed): `s7/between_broad_events.parquet`,
  `a_wa_between/between_broad_events.parquet`, chained by gaps < 24 h as Q5 (secondary test).
- Exclusion: event windows [t0 − W, t1 + W] must not overlap any catalogued episode ± 2 days
  (`episodes_svi_final.csv`, `episodes_nw_final.csv`; SVI bursts vs SVI episodes, WA vs NW).

## Scrambled control (random days)
- For each burst, 1,000 draws: the same template and duration moved to a random start time on another day with the
  same time of day, day-of-year within ± 30 days, any year of that set's span, window fully outside episodes ± 2 d and
  outside every burst ± 3 d, and the same stations having valid data (data-gap matched: m̂ is only computed with the
  stations that have data on the draw; a draw is redrawn if σ_m is > 1.5 × the burst's own σ_m).
- Stack z = (stack − mean of control stacks) / σ of control stacks.

## Positive control
- P1 **big ETS**: episodes with in_frozen_major and span ≥ 20 days, 2009–2026 (SVI and NW lists). Step from
  [start − 10 d, start] to [end, end + 10 d] (shared linear trend). Template = slip ∝ PNSN tremor epicentres in the
  episode per mesh patch (unit total moment). Control = same-length windows at random times outside major episodes.
  **Pass:** stacked m̂ > 3σ of its control AND ≥ 50 % of single episodes > 2σ.
- P2 **smallest catalogued episodes**: major = False (minor), same estimator; descriptive (detected or not).
- Check C (calibration, descriptive): for P1, also predict the gauge steps from the Gualandi (2025) GNSS slip solution
  over the same window (`lfe_extension/results/slip_compare/slip_cache.npz`, L·ΔA) and regress observed on predicted per
  station; slope ≈ 1 means the nominal calibration is right; the Mw bound is also given with this empirical gain.

## Injection test
- Into the actual burst set's scrambled-control windows, add the template × M0 as a ramp over the burst duration
  (step if duration 0) to the corrected hourly series before the step fit, for Mw 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8,
  6.0 (each burst gets the same Mw). Recovery = stack z > 3 against the un-injected control distribution; 200 injected
  stacks per Mw. Also per-burst recovery (single m̂ > 3σ_m) for the curve.

## Pass rule (frozen)
- "Bursts are small slow slips (Mw ≈ X)" only if ALL: (i) primary BE stack z > 3 (positive, thrust sense);
  (ii) P1 passes; (iii) Mw 5.2 injection recovered in > 80 % of stacks.
- Otherwise: report the bound = Mw of (max(stack − control mean, 0) + 3σ_control) per burst, with nominal and
  empirical (check C) calibration, and the smallest Mw recovered in ≥ 80 % of injected stacks.
- Mw = (2/3)(log10 M0 − 9.1).

## Limits stated now
- Raw-gauge processing is ours, not EarthScope's Level 2; tidal calibration is nominal; ocean loading is absorbed by
  the per-year harmonic fit (constituent-wise) but not modelled; rain/groundwater are not corrected (they land in the
  scrambled-control spread).
- Times: event files are UTC; text gives Pacific time.
