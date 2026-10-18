---
name: dvv-temporal-estimators
description: "Benchmark-settled estimator strategy — boxcar for background fields, pairwise doublets for events; TV fakes precursors; per-window absolute-scale fix (--origin-anchor, _calT)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 18c98477-b5c6-4af0-9788-f240c207a136
---

Settled 2026-06-10 by differential step-injection benchmark on real B928 stacks
(`scripts/bench_step_recovery.py`, `scripts/dvv_pairwise.py`; notes/2026-06-10_Notes.md; design: inject
known stretch step into daily stacks, run estimator on injected AND control arms, subtract).

**TWO-ESTIMATOR STRATEGY:**
- **30-d causal boxcar + SVD-Wiener (production)** = BACKGROUND/tomography estimator. Best stability
  (0.071% on B928). Turns any true step into a ~24-d ramp STARTING at the event; SVD-Wiener does NOT
  leak backward (≤4%) and shrinks amplitude only ~5% — the denoiser is innocent; the boxcar is the smearer.
- **Pairwise doublet graph (`scripts/dvv_pairwise.py`)** = EVENT/transient estimator. Staggered-lag (≤14 d)
  day-pair stretches inverted by robust WLS; no temporal waveform averaging. Earthquake-size step: rise50
  1 d / rise80 2 d, amp 0.96, leak nil; sub-noise step: rise80 12 d. Noise 0.12% (~1.7× boxcar). Its
  vectorized spline stretch ≡ stretch_dvv (3e-17) and ~6× faster.
- A real sudden change must be sharp in pairwise AND a 24-d ramp starting the SAME date in the boxcar.
- **TV/ℓ1 smoothing DISQUALIFIED for event timing: leaks 27–34% of a step BEFORE the event** (precursor-faker).
- Evaluations of step recovery MUST be differential (inj−ctrl); raw metrics are contaminated by natural
  variability (burned an Opus agent's report and my own v1).

**ABSOLUTE-SCALE FIX:** stretch about sample 0 (t[0]=−3 s) under-reads dv/v per window — 1–3 s 0.30×,
2–4 s 0.44×, 3–5 s 0.56× of true (lapse-ratio + S-dilution in 1–3 s). So the multi-window inversion mixed
inconsistent scales (needs redo) and all absolute claims need rescale (ETS bound 0.019%→~0.044% true).
Fix: `dvv_roll30cal.py --origin-anchor` (trim t≥0) → 0.97–1.02× all windows; validated B928 (corr 0.997
to old shape, ratio 2.29, cc unchanged). True-scale generation = `data/daily_dvv_<STA>_<win>_calT.csv`
(re-measure launched 2026-06-10; `_cal` files kept). See [[coda-window-2to4-correction]], [[success-criteria]].
