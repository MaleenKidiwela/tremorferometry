# S11 — design-only resolution forecast for per-episode deep change (Merlin, 2026-10-12)
Model: Claude Fable 5.1. Design-matrix analysis on `tables/v2_cells.csv` (1,304 cells, 50 episodes) and `kernel_coef.parquet`; no new measurement. "Scrambled control" = station data shifted to another time of the same day.

## Finding that reframes the problem
The per-episode σ_D of 0.19–0.21 % (big episodes; 0.41 % median over all 50) is dominated by the unidentifiable common level of D (1 %/√50 = 0.14 % from the prior), not by the deep/shallow trade-off. The failed prior test is that level moving.

With a free shallow term per site per episode, any per-episode deep contrast is data-unidentified: σ 1.7 % per episode, 14 % for the 15-major composite, at prior 1000. The reason is that a_d is constant within each series, so the free series offsets absorb it. No lever changes this. The same-station family-pair double difference is this model written pairwise (a_s varies only ±0.05 within a site), so it is equivalent.

## What is resolvable: ΔD(e) under a smooth shallow term
Temporal model: S per site = 2 harmonics + linear + hinges, as v2 "smooth"; D free per episode. ΔD(e) = D(e) − mean D of the other episodes within ±3 yr.
- σ: 0.10–0.13 % per big episode; 0.08–0.11 % against the all-episode mean; D(2012-08) − D(2011-08) 0.08 %. Prior-insensitive (< 10 % change at prior 1000). Composite σ 0.05 %.
- Leak:
  - 0.1 % common shallow pulse in one episode → 0.029 (0.29); 15-major composite 0.41;
  - PGC-site-only 0.3 % pulse → 0.22; other sites ≤ 0.03;
  - annual 0.1 % → ≤ 0.003.
- Deep pulse recovery 0.89–0.97.

## Levers
| rank | lever | σ_ΔD per big episode | leak | cost |
|---|---|---|---|---|
| 1 | temporal model (enabling) | 0.10–0.13 % | 0.29 (0.41 composite) | code, ~1 h |
| 2 | n ≥ 20 cells, 1/σ² (already measured) + restore 32 B010 cells (infl overflow) | ×0.75–1.0 | same | minutes |
| 3 | 39 missing SVI families (no new depth range) | ×0.85–0.9 | same | windows ~1 h, 273 kernels ~6 min |
| 4 | distant stations 100–150 km | ×0.7 (speculative) | + mid-crust −0.4, > 45 km +0.2 | new picks/windows, new kernels, least-trusted geometry |
| 5 | P ≥ 0.3 repeats | ≤ 5 % | — | low-P source needed |
| — | double difference | = unidentified model | — | — |
| best | 1 + 2 + 3 | 0.07–0.10 %; composite 0.03–0.04 % | 0.3 / 0.4 | ~3 h |

## Recommendation
Yes, ≤ 0.1 % per big episode is reachable for ΔD(e), the deep-weighted, episode-locked contrast, with a declared ≤ 0.3 shallow leak. It is not reachable for a per-episode deep/shallow attribution, which is impossible at any cell count.

Make Q3 = ΔD(e) and its 15-major composite the target. Gates:
- prior ×2/×0.5 moves ΔD by < 0.1σ;
- scrambled-control ΔD within 2σ;
- deep-pulse recovery 0.8–1.1;
- leave-PGC-out within 2σ. PGC holds 64 % of cells, and an ETS-locked estimator systematic maps 100 % into D here, so the control gate is the guard.

Lever 4 only if steps 1–3 leave σ > 0.1 %, and then with a_m and a_b free.

## Prescription (executor steps)
1. Correct NOTES/RESULTS: the per-episode σ was dominated by the common level, and per-episode attribution is unidentified (σ > 1 %). Reproduce the prior-1000 numbers with `s11_invert.design/demean`; if they cannot be reproduced, stop.
2. Pre-register Q3 = ΔD(e) (temporal model) + the 15-major composite, with these pass rules:
   - prior ×2/×0.5 moves ΔD by < 0.1σ;
   - scrambled-control ΔD within 2σ in ≥ 90 % of majors;
   - injected one-episode deep pulse recovery 0.8–1.1;
   - injected common shallow pulse leak ≤ 0.35;
   - leave-one-site-out ΔD within 2σ.
   Failure of the control or leave-PGC-out gate means ΔD is a station systematic.
3. Re-run with nmin 20 and the 32 B010 cells restored; then add the 39 families (measure first, then kernels). Acceptance: σ_ΔD ≤ 0.1 % for ≥ 5 majors.
4. Lever 4 only if step 3 leaves σ_ΔD > 0.1 %.
