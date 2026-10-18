# SVI LFEs vs modeled slow slip — Merlin design, 2026-10-07 (~7 am PT). Model output; verify before acting.

Verdict: viable only as a regional / along-strike study. Family-scale slip comparison is below geodetic resolution.
- Resolution: 14 GNSS stations in the SVI box (ALBH, PGC5, SC02, SC04, NEAH, PTRF, PTAA, NANO, P435/P436/P439, CLRS, COUP, +1); 36 Slab2 patches under SVI (29–56 km, ~16.7 km), smoothing 21 km, 13 margin-wide components → family cloud (48.12–48.93 N, 28.8–46.6 km) spans ~5×4 patches, expect 2–4 resolved slip parameters. The LFE residual coherence (0.73 at < 3 km, gone by 20–40 km) is below that scale. Depth (down-dip) gradient < 2 patches → prior-driven; "LFEs per mm vs depth" forbidden unless resolution test passes.
- Slip reproduction unvalidated (per-component misfit 0.52–0.94; fig08/fig09 peaks look like daily noise spikes; IC 26 −90 excursion 2024.8–2025.0 coincides with 2025 episodes).
- Low-rank leakage: SVI slip = Σ of 13 margin-wide patterns × shared time functions; a Washington-only ETS can leak into SVI.
- Detrended convention: inter-ETS slip negative by construction; |slip| biased; use signed rake-projected slip.
- Circularity via tremor: none (component selection uses Hurst/frequency criteria only). 31-day smoothing only in the strain netCDF.

Prescription (freeze PLAN.md first):
1. Tier A, model-free GNSS offsets per episode at the 14 stations (5-day medians, ±2 d, own ETS direction, station bootstrap; null = random windows ≥ 10 d from episodes). Accept: ≥ 15 of 23 major episodes > 3σ; fail < 8 → only the detection curve (step 7).
2. Tier B gate: resolution trace over 36 SVI patches; N/S and up/down-dip spike tests; modelled vs observed offsets at ALBH/PGC5/SC04 (slope 0.7–1.3) and per-episode Mw vs Michel 2019 / Gualandi 2025; IC 26 share. Along-strike use if trace ≥ 2 and N/S spike ≥ 50 % recovered with < 30 % cross-talk.
3. Headline A: log N = a + b log D across major episodes; ratios N/D, X/D, N/X trend tests. New if b ≠ 1 by > 2σ, or N/D trends while N/X does not (p < 0.01). LFE moment vs geodetic moment as a column.
4. Headline B: north/south family groups (≥ 25 km apart): LFE share vs slip share per episode (Tier B and Tier A subsets); controls shuffle and 180° rotation; pass ρ > 0.5, p < 0.01 in both tiers.
5. Per-family GLM local vs regional slip only after 4 (expected: unresolvable at 21 km).
6. Timing (superposed epoch) as validation; claim only |lag| ≥ 3 d vs random-day control.
7. Inter-ETS: signed Tier-A rate on 262 quiet-day LFE days vs 1,000 season-matched random sets (> 99th pct); detection curve for 56 non-major episodes.
8. Cut: depth dependence, migration vs slip-front speed, per-family moment comparison.
Run first: steps 1 and 2. Would change mind: trace ≥ 6 with up/down-dip separated (per-family + depth come back), or Tier-A S/N < 3 for most episodes (only step 7 survives).
