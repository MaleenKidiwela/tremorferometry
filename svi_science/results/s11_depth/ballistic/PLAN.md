# S11 ballistic test — do P and S arrival times move from episode to episode, and is it the path or the source? (rules fixed first)

Written Mon 2026-10-12, 2:45 pm PT, before any lag was computed. User-approved cheapest test of MERLIN_IDEAS.md "Idea 3"
(ranked test 1): extend `s10_ps_timing.py` (family 10, PGC + B011, ±500 ms search → P cycle-skips) to every SVI family
and core station with a ±100 ms search. "Scrambled control" = the detector's c1 rows: each station's data shifted to
another time of the same day, read at the detector frame time + offset (the S10 `__c2` windows), i.e. noise through the
same path.

## 1. Data (reused, nothing re-read)
- Windows: `s10_dvv/multi/win/<NET.STA>.npz` (data rows, P ≥ 0.5, 2–8 Hz, 40 Hz, despiked, RMS-normalised, components
  Z/N(1)/E(2), τ = −14…+10 s about the catalogue-predicted S = detection time + fixed (family, station) offset) and
  `<NET.STA>__c2.npz` (scrambled control). **Windows are NOT realigned** (the S10 dv/v realign would erase the timing).
- Stations: PB.B009, PB.B010, PB.B011, PB.B926, CN.PGC, CN.VGZ, CN.SNB. NLLB has no windows → not used.
- Families: the 76 S10-multi families. Episode of a row = `multi/tables/active_days.csv` (fam, date → episode unit; 52
  majors + 123 bursts). Rows on non-active days have no episode and are used only in the grand stacks.
- Sensor epochs: `s10_common.EPOCHS`; an epoch is never crossed (each series has its own zero).

## 2. Lags (per family, station, epoch, episode)
- Episode stack: P-weighted mean of the episode's rows, **n ≥ 20 rows**, 3C. Reference: the station-epoch
  **leave-episode-out** P-weighted grand stack (all rows of the series minus the episode's).
- **S:** S time t_S picked on the full grand stack (`s10_measure.pick_S`, 3C envelope, ±2 s of predicted S). Window
  [t_S − 0.5, t_S + 1.5] s on N and E jointly (one normalised cc summed over both components).
- **P:** t_P = maximum of the Z envelope of the grand stack within t_S − SP ± 1.0 s (SP = `multi/tables/offsets_sp.csv`).
  Window = full P wavetrain [t_P − 0.5, min(t_P + 2.5, t_S − 1.0)] s on **Z plus the polarised horizontal** R (the
  horizontal direction of largest grand-stack P energy in that window, from the 2×2 covariance), joint cc. Z-only lag kept
  as a column.
- Search ±100 ms (±4 samples), parabolic sub-sample refinement. Convention: positive lag = arrival later than in the
  reference. A peak on the search edge = **railed = not a measurement** (missing, counted).
- Errors: 200 bootstrap draws of the episode's rows (P-weighted, fixed seed); S and P lag from the same draws, so σ(S−P)
  is direct.
- Self-check (not a gate): the episode stack Fourier-shifted by +5 ms must read +5 ms ± 0.5 ms in the median.

## 3. Observables
- (a) **S−P** per (family, station, episode) = lag_S − lag_P (origin-time free).
- (b) **Inter-station S difference** dS = lag_S(X) − lag_S(PGC) per (family, episode), computed on the **rows present at
  both stations** (same rows → the catalogue-time drift and the 20 Hz detector-frame jitter cancel exactly), each station
  against its own leave-episode-out reference, bootstrap with common draws. Pairs: PGC–B926, PGC–SNB, PGC–VGZ (far),
  PGC–B009, PGC–B010, PGC–B011 (co-located floor), plus B011–B009 and B011–B010 (extra floor). dP the same way (column).
- (c) **Floor** = scatter of dS for PGC–B009/B010/B011.
- Every series ((family, station, epoch) for S−P; (family, pair, epoch pair) for dS) is demeaned over its episodes;
  series with < 3 episodes are dropped. Scatter = √(Σr² / Σ(K−1)) (dof-corrected) and a robust 1.4826 × MAD version.

## 4. Key test (headline set: all episode units with n ≥ 20; majors-only reported as a second column)
Per far pair (and per station for S−P): one-way decomposition of the demeaned values by episode (episodes with ≥ 3
families). MSB, MSW → **episode-common variance** σ²_ep = (MSB − MSW)/n₀; **within-episode variance** MSW, split into
**noise** (mean bootstrap σ²) and **family-specific** σ²_fs = max(0, MSW − noise). Statistic F = MSB/MSW;
**permutation**: 10,000 draws, each series' values shuffled among its own episodes (seed 20261012),
p = (1 + #{F* ≥ F}) / 10,001.

## 5. Pass rules (pre-registered)
- **Floor ≤ 2 ms** (dof-corrected scatter of the co-located dS) → timing precision confirmed. > 2 ms → not confirmed,
  and every other number is read against the measured floor.
- Far-pair scatter "well above the floor" = > 2 × floor.
- **PATH**: σ²_ep > σ²_fs AND p < 0.01 → real path/velocity signal worth a ballistic inversion.
- **WANDER**: scatter > 2 × floor and σ²_fs ≥ σ²_ep (or p ≥ 0.01) → source-locus wander, needs a nuisance model.
- **MIXED**: p < 0.01 but σ²_ep ≤ σ²_fs → a common part exists but wander dominates.
- **FLAT**: scatter ≤ 2 × floor → nothing to measure.
The verdict is per pair (B926, SNB, VGZ) and per station for S−P; the overall call = the majority of the three far pairs.

## 6. Conversion to δβ/β (any common signal)
Hypocentral S time T_S = SP × 1.75/0.75 per (family, station) (Vp/Vs 1.75). dS common term a_e: uniform δβ/β on both
legs → δβ/β = −a_e / (T_S,X − T_S,PGC); change on the X leg only → −a_e / T_S,X (both reported). S−P common term:
constant Vp/Vs → δβ/β = −δ(S−P)/SP; Vs only → −δ(S−P)/(2.33 SP). The floor is converted the same way (sensitivity).

## 7. Scrambled control
The c2 rows through the identical path (per-station episode stacks, n matched to data ≤ available, ≥ 20, equal weights,
against the data leave-episode-out reference; pair dS on c1 rows present at both stations). Reported: railed fraction,
cc distribution, scatter and the same decomposition/permutation on its non-railed values. If the control shows an
episode-common part with p < 0.01, the data PATH call is withdrawn as a path artefact.

## 8. Figures (user style: points coloured by cc, no error bars, line graphs)
- `fig_dS_series.png`: per pair, demeaned dS per episode vs date; thin grey line per family, points coloured by cc,
  thick line = episode-common mean.
- `fig_SmP_series.png`: the same for S−P per station.
- `fig_variance_decomp.png`: bars (episode-common, family-specific, noise) per pair/station, data and control.

## 9. Conduct
CPU only, ≤ 8 workers, RAM < 40 GB (another queue, WA-NW S11, is running); watchdog kept alive; kill by PID only;
nothing deleted. Scripts `svi_science/scripts/s11_ball_*.py`; outputs here (tables/, figures/, logs/).

## Amendment 1 (Mon 2026-10-12, 2:50 pm PT — after a test run on family 10 only, BEFORE the 76-family run; no key-test number seen)
What the family-10 test showed (tables/*_test.*):
1. **P pick on the search edge.** The catalogue SP (offsets_sp.csv) is 10–20 % shorter than the observed S−P (PGC 5.2 vs
   6.2 s), so the ±1 s P search hit its edge. Fix: t_P = maximum of the Z envelope of the station's **all-epoch** pooled
   grand stack within [t_S − 1.3·SP, t_S − 1.0·SP] s (one t_P per family–station, the same window in every epoch; a first try with [t_S − 1.35·SP − 0.5, t_S − 0.85·SP] put VGZ's P on the early edge of the search, 4.8 s instead of the 3.6 s seen in its envelope; observed/catalogue S−P is 1.08–1.19 at PGC, SNB, B926, VGZ).
2. **PGC timing in E3.** PGC–B009/B010/B011 dS jump together by ~10 ms in some 2017–2025 episodes (e.g. 2021-10) while
   B011–B009 stays within ~2 ms: a PGC-side timing effect that is episode-common by construction and is NOT path (the
   sites are co-located). With PGC as reference it would masquerade as a far-pair "path" signal. Added before the run:
   - far pairs also referenced to **B011** (B011–B926, B011–SNB, B011–VGZ), and the floor pairs get the same
     decomposition/permutation;
   - **PATH additionally requires** (i) the far-pair episode-common variance to exceed the co-located episode-common
     variance of the same reference station, and (ii) the call to hold with both references (PGC and B011).
   - The floor is reported per reference: PGC–boreholes and B011–B009/B010; the pre-registered floor rule (≤ 2 ms)
     is applied to each.

## Amendment 2 (Mon 2026-10-12, 2:58 pm PT — while the 76-family run is in flight, no output read)
Two definitions the plan left open, fixed before analysis:
- **S−P floor** (for the per-station S−P verdicts) = dof-corrected scatter of d(S−P) = dS − dP over the co-located pairs
  (PGC–B009/B010/B011, B011–B009/B010), demeaned per series like everything else.
- **SP in the δβ/β conversion** = the observed grand-stack S−P (t_S − t_P picks, `tables/series_geom.csv`), not the
  catalogue SP; T_S = 1.75/0.75 × SP.
