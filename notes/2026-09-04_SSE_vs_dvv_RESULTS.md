# Slow-slip model vs LFE-coda dv/v — consolidated results (2026-09-02 → 09-04)

Goal (user): find a pattern in how the slow-slip model's variables relate to δβ/β (and to coda dv/v generally).
Model: Gualandi (2025) vbICA slow slip via github MaleenKidiwela/cascadia-sse-strain, re-evaluated in OUR measurement
frame (30-d trailing mean on the 13 amplitude functions, gauge 2019-24) at our 872 interface cells and on the model grid;
22 strain/stress fields + slip. Products: `fault_tomography/inversion/res_catalog_sse/`, map `figures/sse_model_map.html`.
Working ledger (every number, in order): `notes/2026-09-03_Notes.md`. Merlin rulings: `notes/2026-09-03_merlin_sse_plan_ruling.md`
+ two audits recorded in the ledger. Scripts: `fault_tomography/inversion/sse_*.py`.

## Bottom line
1. **No detectable response of the interface δβ/β to modeled slow slip.** Pre-registered footprint test over the 21
   tremor-confirmed north/central episodes: +0.03 ± 0.07 % (rigid-shift and seasonal-phase nulls; ring, segment split,
   controls, leave-one-out all clean). 95 % bound on |footprint response| ≈ 0.17 % (production), 0.22 % (drift-gated
   product), in estimator units. Family-level 2-4 s product: bound 0.065 %.
2. **Stress/strain hypotheses are untestable with point evaluation.** Stress 10 km above the interface is a near-
   cancellation of vbICA components; footprint-median shear change is negative in only 10-14 % of episodes at 0.5/10 km.
   Slip-rate is the only usable predictor. Near-station dv/v vs modeled strain at the station (169 stations): null;
   sensitivities |d(dv/v)/dε| ≳ 1e4 excluded.
3. **What IS in the data: an ETS-locked, station-specific apparent dv/v shift.** Nearly every station's LFE stacks
   change during modeled slip (median |shift| 0.18 % in 2-4 s; 93/167 stations > 3 SE in a direct ON-vs-PRE stack test).
   It is not a medium velocity change: sign set by the station, not the source cell (stations sharing cells disagree
   0/10, opposite to the kernel); sign flips between the 2-4 and 3-5 s windows at 30 stations; the pre-arrival noise
   window carries an equal or larger shift; amplitude scales with the pre-arrival noise-level change (ρ +0.34); it does
   not replicate across time halves (18 % same-sign, below chance), even for the 26 stations that were lapse-consistent.
   Excluded mechanisms: SVD-Wiener shrinkage (identical without the filter), detection-count change, tremor intensity,
   modeled site strain, regression to the mean, depth/year/distance. Remaining: an ETS-locked change of the stacked
   LFE waveform (source population / noise field / equal-weight stacking of ETS-day waveforms) → a SYSTEMATIC of the
   S-pinned stretch estimator, to be treated as such in any ETS-scale claim.
4. **Leads that arose and were killed (with the killing test):**
   - 49.2 N deep patch, +0.34 % during slip in every product, FWER-safe (z 4.9): station-side (same-cell sign
     inversion B011 vs B928; lapse flip; mirror window). Merlin audit F1-F3.
   - Regional sign pattern (north +, south −, +0.11 % p 0.01): reverses in the 3-5 s window; same class.
   - 3-5 s footprint DECREASE (−0.05 %, mean p 0.03/0.005): carried by the July-2018 northern episode whose driver
     stations (B010, PGC, SNB, BFSB, CBB) flip sign vs 2-4 s of the same stacks with 7-17 % rail fraction and 0.1-0.2
     coherence loss in the ON window; median/sign-count not significant; station split-half 43 % agreement; south
     opposite sign. Merlin audit: not defensible. Decisive block-gated lapse test (15,274 family-episode pairs, 113
     stations): Theil-Sen slope of the 3-5 s on the 2-4 s contrast 0.45 [0.33, 0.56] (a velocity change needs ≥ 0.6),
     sign agreement 62 % (needs ≥ 70 %), gated pooled mean −0.03 %, Null-A p 0.18 (needs ≤ −0.04 % at p < 0.01) → all three fail →
     waveform-change class, BOUND. Closed.
5. **Model-side findings (independent of dv/v):** the model's episodes at 42.5-45 N coincide with tremor at chance
   (29 %), are winter-locked and sign-symmetric → not slow slip (hydrological/seasonal leakage through vbICA); north/
   central episodes are 76-80 % tremor-supported. Earthquake offsets (Dec-2024 M7.0 etc.) sit inside the amplitude
   functions and must be gated.

## Round 4 (2026-09-04, after the user asked to keep investigating without Merlin; large events; cumulative results)
6. **Tremor-defined episodes, margin-wide (211 episodes incl. the south):** interface +0.004 % (105/211 decrease,
   p 0.94); family 2-4 s −0.014 % (p 0.46); family 3-5 s −0.036 % (p 0.17); near-station −0.007 % (p 0.67). With five
   times more episodes the model-catalog 3-5 s "decrease" vanishes (N+C −0.010 %, p 0.65). Tighter bounds: interface
   < ~0.13 %, family 2-4 s < ~0.04 %, family 3-5 s < ~0.05 %.
7. **Large events:** per-episode significance against each episode's own footprint null: interface 2/32 beyond 2σ
   (1.6 expected), no relation between response and episode size (Spearman −0.03). Top-quartile and top-5 northern
   composites: interface −0.08 % (p 0.39) and −0.22 % (p 0.09); the family 3-5 s top-5 (−0.16 %, p 0.003) still rests
   on the artifact-flagged July-2018 episode. Largest events do not respond either.
8. **Cumulative slip / stress state, rise and fall (14 regions × 5 products):** lagged level correlations with
   cumulative slip, Coulomb, shear, normal and dilatation show ~7-12 of 59 cells at p < 0.05 per variable but with
   contradictory signs between the 2-4 and 3-5 s windows in the same region → no coherent relation. Phase test
   (dv/v drift during slip release vs loading): regional 48.2N family products −0.03 %/30 d (p 0.005/0.008), but at
   station level it is a weak 60 %-negative lean that does not replicate in 3-5 s across station halves; fleet-wide
   station phase test null (9/154 beyond 2σ, 7.7 expected). Daily spatial correlation between the slip-rate map and
   the dv/v-change map: null (p 0.2-0.35).
9. **Kernel-weighted model fields (each pair's own coda sensitivity applied to the model strain/stress, 4,940 pairs):**
   cumulative slip, dilatation, Coulomb, shear, normal and all six strain components correlate with pair dv/v at chance
   level (5.2-5.9 % of pairs beyond their own null, mean r ≈ 0); pooled regression of 30-d changes R² = 0.00003.
   Slip-in-window shows 19.8 % of pairs beyond null with zero mean and 46 % positive: the station-specific two-signed
   shift again. Null.

**End state (2026-09-04):** every lens exhausted; no relation between modeled slow-slip variables and δβ/β found;
the only robust pattern is the station-specific ETS-locked estimator systematic.

## Round 5 (2026-09-04, user: "I know there's a pattern" / hint: slip patches across time)
10. **Spatial dipole across dip.** The change maps of the three largest northern episodes show up-dip positive /
    down-dip negative structure, which footprint means cancel. Tested formally (per-episode spatial correlation of the
    change map with 33 model change maps at 3 depths, the slip map and a geometry template; both catalogs; shifted-
    episode null; family-wise): nothing survives (best FWER p 0.10; interface vs east-west strain change r +0.02-0.04,
    same sign in both catalogs). Across-dip profile aggregated over 21 episodes: shape present but inside the null.
11. **Patch-cycle folding** (each slip patch followed through its own inter-slip cycle): interface inside the null;
    family 2-4 s drifts down through the cycle (slope p 0.000, only 4 distinct patches), not seen in 3-5 s.
12. **Moving slip patch, inside vs outside at the same time** (718 northern days): interface −0.13 % at +40 d after
    the slip map (z −2.8; lag-family-wise p 0.10; negative in both time halves and 11/15 years; −0.08 in the gated
    product; ~−0.02 in family products). The most hypothesis-consistent residual of the whole search; suggestive, not
    significant.
13. **Year-scale before/after per patch:** interface 15/21 episodes lower in the following year (−0.07 %, p 0.39).
14. Station map vs model spatial variables; cross-segment transfer (48 cells, mixed signs); long lags to ±2 yr; scatter
    change during slip; long-term trend map vs total slip map; nonlinear gradient boosting (max CV R² 0.006); other
    coda observables (energy ratio, noise level): all null.
15. **Re-inversion with the S-pinned differential kernel** (`splits/Gdiff24`; λ 0.215, checkerboard recovery 0.93 vs
    0.60 for the production operator, aliasing probe 1.00): pre-registered footprint test +0.011 ± 0.007 % → BOUND
    |response| < 0.026 %; across-dip profile flat; the northern inside-vs-outside drop of item 12 is not reproduced
    (north +0.02 %, central +0.03 % at lag 0); whole-record patch levels: 48.2N vs cumulative slip r +0.21 (p 0.01),
    49.2N vs slip-in-window r +0.10 (p 0.01, the station-side effect), others null.

**End state (2026-09-04, 16:45):** the sharper corrected-kernel product tightens the bound six-fold and removes the last
suggestive residual. No pattern linking the modeled slow-slip variables to δβ/β was found in any framing tried.

## Standing products for future work
- Gated episode catalog (32; N 11, C 10, S 7): `res_catalog_sse/episodes.csv`, footprints + local onsets.
- Model in our frame: `components_at_cells_trailing30.npz`, `megathrust_strain_trailing30.nc`, `strain_gf_offsets.npz`.
- 3-5 s dv/v products for 184 stations: `data/daily_dvv_<TAG>_Z_3to5.csv` (map toggle in `figures/borehole_dvv_map.html`).
- Fleet direct stack test table: `res_catalog_sse/stack_lapse_fleet_classified.csv` (per station: e24, e35, ratio,
  noise ratio, time halves) — the reference table for the estimator systematic.
- Round-4 scripts: `sse_tremor_episodes.py`, `sse_large_events.py`, `sse_cumulative_regional.py`, `sse_kernel_weighted.py`;
  tables `tremor_episodes.csv`, `large_events.csv`, `cumulative_regional.csv`, `station_phase_test.csv`, `kernel_weighted_pairs.csv`.
- Test scripts: `sse_catalog.py`, `sse_primary_test.py` (DVV_PRODUCT env), `sse_pattern_search.py`, `sse_explore2.py`,
  `sse_explore3.py`, `sse_cycle_coherence.py`, `sse_verify_492.py`, `sse_stack_lapse_fleet.py`, `sse_lapse_gated_test.py`.

## Wording to use (Merlin-approved)
"Across the interface, gated, and family products in both coda windows, the slow-slip footprint contrast is consistent
with zero; the one look that gave a mean-significant result (3-5 s family product, −0.05 %, Null-A p 0.03 on the mean
only) is carried by the July-2018 northern episode, whose 3-5 s decrease at B010, PGC and SNB is opposite in sign to the
2-4 s measurement of the same stacks and coincides with a 0-15 % rail fraction and a 0.1-0.2 loss of coda coherence,
i.e. an estimator response to a change in the stacked waveform. We report no detection and 95 % bounds on |footprint
response| of ~0.065 % (2-4 s) and ~0.09 % (3-5 s), in estimator units, uncalibrated; window-dependent apparent shifts of
this kind are a systematic of the method."
