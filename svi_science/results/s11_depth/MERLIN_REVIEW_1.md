# S11 — Merlin review 1 of the SVI pilot (2026-10-12) — summary saved by the main session
(The full file text was returned only in Merlin's transcript; this is the hand-back summary, verbatim in substance.)

Verdict: pipeline sound; per-episode deep D is structurally unresolvable at SVI (σ_D(e) ≈ 0.4 %, set by the within-site a_d spread ≈ 0.13 × √cells, not by the cc cut, site count, WA or POLARIS). Headline quantities change to:
- Q1: deep annual amplitude (σ ≈ 0.04–0.05 %);
- Q2: decade contrast of D (σ ≈ 0.13 %), under a pre-registered smooth S/D parameterisation.

Per-episode D is reported as a bound. Do not scale now.

1. Cell rule v2: n ≥ 50 AND cc3c ≥ 0.35 (99th percentile of scrambled-control cc; noise-derived), 1/σ² weights with control inflation. Justification: bootstrap σ tracks the series-demeaned scatter in every cc bin (rms/σ 1.1–1.3, χ²_red 1.2–1.8 flat); corr(σ, cc) = −0.83. cc ≥ 0.8 stays as a column. About 1,315 cells. σ_D 0.44 → 0.41; Q2 σ 0.17 → 0.13.
2. Parameterisation: S per site = 2 annual harmonics + linear + 4-yr hinges (10 dof per site, no per-year steps; the per-year steps are why "harm" failed). D = same slow basis + annual. Rejected: a regional S per episode (no gain, 0.45 %); H/Z (no lever, kernels identical); a hydrology proxy (a forcing prior). Joint-lapse (a_d/a_s 2.01 → 1.46) is kept as a column.
3. Gate (a) carries no depth information: the kernel predicts vault = 200-m borehole for every layer (paired Δa_s ≤ 0.009, Δa_d ≤ 0.004). It only rules out instrument or vault artefacts. Depth is assigned by an annual-only inversion using the distance lever (corr(a_d, hypocentral distance) +0.6–0.9), gated by:
   - per-site replicates;
   - the station-amplitude pattern (deep ∝ a_d of the site);
   - the lapse ratio (deep 0.83, shallow 1.14).
   Current evidence for a deep annual term: 0.114 ± 0.048 % (2.4σ), leave-one-site-out 0.08–0.15, with an SNB tension (largest a_d but smallest amplitude, 0.025 ± 0.020). Fragile; report honestly.
4. Scaling does not break the degeneracy. WA NW has its own D (no effect on SVI); POLARIS adds sites with their own S (√cells only). Scale later only to replicate Q1/Q2, after a design-only resolution forecast.
5. Minimal next run (< 1 h, no new measurement):
   - amend PLAN first;
   - re-run s11_invert.py with models episode / smooth / slow+annual / joint-lapse on the v2 cells;
   - re-run injections A+B under the new model, plus coherent kernel-error runs (generate with two5_40/ell80 a_s, invert with ell40; design-level leak 0.02/0.11, recovery 0.87–0.96);
   - gates for Q1: leave-one-site-out within 2σ and peak ±60 d; family split; scrambled control; station-amplitude χ² against ∝a_d and ∝a_s; lapse ratio.
Flaw labels: gate (a) has no depth content (MAJOR); the lever is distance-based and model-conditional (MAJOR); per-episode D is noise-limited by geometry (MAJOR); the cc cut was a mistake but not the limiting factor (MINOR); annual–slow D correlation via ETS timing (MINOR); SNB tension (SPECULATIVE).
