# Temporal resolution & stability of dv/v — methods program (2026-06-09)

**The user's question:** improve temporal resolution AND stability of the dv/v series without
median-over-time. Current production = 30-calendar-day causal trailing waveform stack (1-day steps,
≥5 stacks) + full-record SVD-Wiener + stretch 2–4 s. Sometimes works, sometimes doesn't; unclear
whether the denoising makes earthquake-produced sudden drops look non-sudden.

---
## 1. What is actually at stake (why "ambient-noise people see drops with 30-day averages" ≠ all clear)

1. **Shape/timing, not detection.** A causal trailing boxcar turns a true step into a ~30-day ramp
   *starting at the event*. Ambient-noise coseismic drops are large (0.05–0.5%, shallow kernel +
   strong shaking) so the ramp still reads as a drop. Our deep-fault signals are ~50× smaller — a
   small step smeared into a ramp can sink below noise where a sharp step would be detectable.
2. **Backward leakage.** Centered moving averages / acausal filters make the drop start BEFORE the
   event — a known source of spurious "precursors" in the literature. Our boxcar is causal (safe);
   the **SVD-Wiener basis spans the whole record (acausal)** and can in principle move energy across
   the event. This is the specific worry to quantify.
3. **Waveform-stack vs series-average.** Ambient-noise studies average the dv/v *series*; we stack
   *waveforms* then measure. During the 30 post-event days the stack MIXES unstretched and stretched
   waveforms — the stretch fit on a mixture is not the mean dv/v: it can lock onto one population,
   sit in between with depressed cc, or jitter. Empirical question → benchmark.

## 2. ⚠ DISCOVERY (benchmark calibration, 2026-06-09): absolute dv/v scale is ~0.44×

`stretch_dvv` stretches about **sample 0** of the trace (t = −3 s), but physical stretch is about
the **LFE origin** (t = 0). In the 2–4 s window (= 5–7 s from sample 0) a true stretch a reads as
≈ a·(lapse_from_origin/lapse_from_sample0) ≈ 0.44·a (measured: injected −0.150% reads +0.0657%;
sign is the project's convention). **Uniform scale factor: all temporal structure, cross-station
comparisons, tomography patterns unaffected; ABSOLUTE numbers are ~2.3× larger in true dv/v** (ETS
bound 0.019% → ~0.044% true; site terms, checkerboard amplitudes likewise). Fix for future: trim
trace to t ≥ 0 before stretching (origin-anchored lapse axis). DO NOT change production silently —
it rescales every product; do it as one coordinated re-measure + note.

## 3. Method candidates (ranked by expected value)

| # | Method | Core idea | Resolution | Causal? | Step-preserving? | Status |
|---|---|---|---|---|---|---|
| M-A | **Pairwise doublet matrix** | stretch between PAIRS of daily stacks at staggered lags; invert pair graph (WLS/Huber) for v(t). Reference-free, no smoothing kernel; stability from graph redundancy | ~1 d | yes (per-day estimates) | yes (each day keeps identity) | Opus agent prototyping (`scripts/dvv_pairwise.py`) |
| M-B | **Robust causal Kalman with jump acceptance** | random-walk state on per-day stretch; gated innovation confirmed by next sample = accepted step (heavy-tailed process noise) | ~1–3 d | YES (filter) | yes by construction | in benchmark (`daily_kal`) — pre-noise 0.045% < prod30 0.071% on B928 |
| M-C | **TV / ℓ1 trend filtering** | total-variation denoise of per-day series — the canonical step-preserving smoother (L2 smoothers always smear) | ~1–5 d | no (but symmetric, no SVD-style basis leakage) | yes (exactly) | in benchmark (`daily_tv`); watch outlier-latching |
| M-D | **Known-event step fitting** | we KNOW earthquake times: fit dvv(t)=smooth + Σ aₖ·H(t−tₖ) (+ optional log-healing) jointly; turns "is the drop sudden" into a parameter estimate | exact at known tₖ | n/a | exact | cheap to add after M-A/M-B |
| M-E | **Equal-information adaptive windows** | accumulate detections to target SNR instead of fixed 30 d → days-resolution during ETS bursts, honest widening in droughts | variable (reported per point) | yes | partial | design only |
| M-F | **Causal SVD-Wiener** | trailing-only basis (or drop SVD for timing work, keep for background fields) | as prod | yes | n/a | trivial variant if SVD proves leaky |

**Composition that likely wins:** M-A (pairwise graph) as the measurement layer → M-B/M-C as the
state layer → M-D for known events. Cross-family per-day median (spatial) stays allowed.

## 4. Benchmark design (the referee): differential step injection
`scripts/bench_step_recovery.py` — inject a known stretch step into REAL daily stacks (post-t0
traces resampled u(t·(1+a)) about the LFE origin), run every estimator on **injected AND control
arms, subtract** → natural variability cancels exactly, leaving each estimator's pure step response.
Metrics: pre-step noise (raw series), differential amplitude fraction, rise-to-80% days, pre-event
leakage fraction. First (non-differential) run on B928 taught: the injected step (0.066% as-read) is
*comparable to natural monthly variability* → differential design is mandatory. Nonlinear filters
(TV/Kalman/SVD) are applied per-arm so their data-dependence is honestly tested.

## 5. Open questions
- Does SVD-Wiener actually leak the step backward / shrink it? (differential run answers this)
- Mixture behavior of the rolling stack across t0: average, bimodal lock-on, or cc collapse?
- Pairwise method: lag-set, weighting, and drift control (long-lag pairs anchor the datum).
- Kalman/TV robustness to multi-day glitches (2024 spike latching seen in v1) → Huber measurement.
- Per-point uncertainty: pairwise graph gives formal errors; Kalman gives P(t) — report both.
- After a winner emerges: re-run ETS superposed-epoch + B201-type cases with the sharp estimator —
  transients shorter than 30 d (the caveat in ETS_ANALYSIS §5) become testable for the first time.
