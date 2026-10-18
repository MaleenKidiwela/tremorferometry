# S11 lever 4 — distant stations (80–150 km): rules fixed before any number

Written Mon 2026-10-12, 8:20–8:45 pm PT, before any far-station window, pick, SNR, kernel, dv/v or inversion number.
Source: `../MERLIN_FORECAST.md` lever 4 ("rank 4": far pairs > 80 km carry a_d 0.91 ± 0.26, a_s 0.23, a_m −0.34,
a_b 0.20; σ ×0.7 speculative; a_m and a_b must then be free). User (via main session): persevere; coarse results are
fine if new. Design by the executor (no Merlin call: Fable-quota rule); every choice below is fixed now.
Scripts `svi_science/scripts/s11_dist_*.py`; outputs here (`distant/<svi|wa_nw>/{tables,logs}`, `distant/figures`).
No file of round3/, tomo/ or of the earlier S11 runs is written. ε = −dv/v; plots show dv/v = −ε.
"Scrambled control" = the detector's c1 rows (each station's data shifted to another time of the same day), read at
frame time + the S offset, through the identical path.

## 0. Status
- 2026-10-12 8:45 pm PT — PLAN written. Nothing computed.
- 8:31 / 8:38 pm PT — step 1 inventory done for both regions (no STOP): SVI 104 usable pairs at 19 stations, 1,459 expected
  usable cells; WA-NW 136 usable pairs at 29 stations, 2,135 expected cells. Tables `<region>/tables/far_pairs.csv`,
  `far_stations_summary.csv`, `inventory_summary.json`.
- ~8:40 pm PT — code tests on one SVI far station (CN.LZB, tag _test; not results): far measurement, kernels (unit sum 4e-15),
  invert (N0 reproduces the stored SVI Q3 composite to 1e-17 and the stored pulse responses 1.14 / 0.45), gates run.
  Amendment 1 written. Queue launched (§9).
- 9:10 pm PT — QUEUE_DONE (both regions). ~9:35 pm PT RESULTS.md written (incl. post-hoc B2 column, `s11_dist_b2.py`).

## 1. Why far pairs, and what would count as the lever working
Near pairs (≤ 80 km) have a_d ≈ 0.5–0.55, a_s ≈ 0.3. Far pairs have a different mix (forecast: a_d ×1.6, a_s ×0.7,
a_m of opposite sign). Three expected effects:
1. the common-shallow → deep leak of ΔD (SVI 0.45, WA 0.60) and of the annual term should fall, because the shallow
   and deep columns of the design matrix are less collinear;
2. σ_ΔD per major should shrink (forecast ×0.7);
3. a kernel-free deep/shallow discriminator for the annual cycle: the annual amplitude at far pairs relative to near
   pairs (predicted ≈ a_s ratio ≈ 0.7 if shallow, ≈ a_d ratio ≈ 1.6 if deep), and the lapse ratio at far pairs.

**"The lever works" (per region, pre-registered):** in the headline configuration F1 (§6) versus N1 (near cells, same
model): median common-shallow-pulse leak ≤ 0.8 × N1's, AND composite σ_ΔD ≤ N1's. **"Deep annual ambiguity reduced":**
B-injection Q1 leak and the coherent-kernel-error B leaks each ≤ 0.8 × N1's, AND σ(Q1) ≤ N1's. Both reported as factors.

## 2. Pairs and stations
- **Regions.** SVI: the 115 families (76 S10-multi + 39 added; Bostock hypocentres, `data/catalogue/families_v1_3.csv`);
  repeats, c1 rows, active days and episode units exactly as the near SVI run (S10-multi tables + `tables/ext`).
  WA-NW: the 108 families of `wa_nw/tables/inv` (Slab2-replaced depths as there); rows, active days, units from there.
- **Far pair** = (family, station) with 80 < hypocentral distance ≤ 150 km (epicentral distance ⊕ family depth), and
  not already in the region's near cell set (SVI: the 7 core stations; WA-NW: the pairs of `wa_nw/tables/inv/offsets_sp.csv`).
- **Candidate stations:** permanent networks only (temporary arrays XG, XU, ZK, X4, ZZ, YW, YG, XD, 7D, Z2, X9 excluded:
  ≤ 3 years, too few episodes per series); a raw day file (any root of `wa_prod_common.prod_roots()` +
  `lfe_extension/traces_backward`) on ≥ 150 of the region's active dates (SVI 702, WA 440); ≥ 1 far pair. No downloads.
- Station coordinates: WA production station table, `backward_scope/inventory_stations.csv`, else EarthScope FDSN.
  Sensor depths: `tables/sensor_depths.csv` ∪ `wa_nw/tables/sensor_depths.csv`, else FDSN (channel level), else 0 (flagged).

## 3. Step 1 — inventory: expected S, picks, SNR, expected cells
- **Expected S.** 1-D flat-layer travel times (S and P, first arrival = min of the direct ray and the head waves along every
  layer top below the source) in the source region's CVM v1.7 profile (SVI `south_VI`, WA `royer_NW_WA`). Origin offset per
  family = median over its near stations of (template S offset − TT_S). Predicted far offset = origin offset + TT_S(far).
  S−P(far) = TT_S − TT_P (first arrivals) — used only for the pre-P window.
- **Sample.** Per calendar year the 4 active dates with the most data rows of the region (all years with active dates). On
  those dates every data row of the region's families and their c1 rows (as stored, ≤ n data per family-day) are read at
  every candidate station at the predicted S, span −30…+12 s, with the unchanged S10 processing (2–8 Hz Butterworth-4
  zero-phase, 40 Hz, despike 8 × MAD, RMS-normalised, 3C).
- **Per far pair** (P-weighted 3C sample stack; c1 equal-weight stack):
  - S pick = maximum of the summed 3C smoothed envelope² within ±6 s of the predicted S; pick shift = pick − predicted.
  - pre-P window = [pick − S−P − 9 s, pick − S−P − 1 s] ∩ span.
  - **pick-SNR** = peak envelope / median envelope of the pre-P window (Z+N+E summed envelope²; square root taken).
  - **causality** = 3C RMS 2–5 s past the pick / 3C RMS of the pre-P window (the series rule of the near runs).
  - control: the same pick-SNR on the c1 stack (its own maximum in the same ±6 s band).
  - **cc20** (pairs with ≥ 40 sample data windows): median over 20 seeded draws of the joint 3C cc, 2–5 s past the pick, of a
    random 20-window stack against the P-weighted stack of the other sample windows.
- **Usable pair:** causality > 1.5 AND pick-SNR > the 99th percentile of the control pick-SNR over all far pairs of the
  region AND |pick shift| < 5.5 s (not at the search edge).
- **Expected usable cells:** over usable pairs, the potential cells (family, station, epoch, episode) with ≥ 20 data rows on
  dates with a file, whose predicted cc(n) ≥ 0.35, cc(n) = √ρ(n), ρ(n) = 1 / (1 + (1/ρ20 − 1) · 20/n), ρ20 = max(cc20, 0.01).
  Pairs without cc20 are not counted.
- **STOP rule (per region):** expected usable cells < 100 → that region stops after the inventory (`DIST_STOP` marker in
  its folder) and is reported with the measured obstacle (pick-SNR / causality vs distance and station). The other region
  continues.
- Full-read S offset = predicted offset + the pair's pick shift.

## 4. Step 2 — windows and measurement (identical recipe)
- Rows: every data and c1 row of the usable pairs' families on the region's active days, at usable pairs only; read at the
  corrected S on a wider grid **−28…+10 s** (the near grid −14…+10 s holds no pre-P noise at 100–150 km, S−P 12–18 s);
  identical processing (the RMS normalisation then spans the wider window — logged, harmless for equal-weight stacking).
- Epochs: `s10_common.EPOCHS` where the station has an entry, else the WA node/epoch-change table (`s11_wa_inventory.epochs_table`
  logic: node start dates ∪ `wa_pre/w2/epoch_changes.csv`). Series = (family, station, epoch); never referenced across epochs.
  'dead' station-days of `wa_pre/w2/r2_bad_for_dvv.csv` not read.
- Measurement = `s11_measure.run_series` logic unchanged (grand stack, S re-pick ±2 s, per-event 3C realignment ±0.5 s,
  P-weighted day and episode stacks, leave-episode-out reference, joint 3C S-anchored stretch 2–5 (headline), 3–6, 2–3.5,
  3.5–6 s, Z/H, 2-parameter fit, 200-draw bootstrap σ, c2 control stacks), except the pre-P (causality) window, which is
  re-expressed relative to the grid start: [tS − 26.5 s, tS − S−P − 1 s] (near grid: [tS − 12.5, tS − S−P − 1]).
- Cells: the region's headline rule — n ≥ 20, cc3c(2–5) ≥ 0.35, causality > 1.5, finite ε and σ, series ≥ 3 cells, 1/σ²
  with the series' scrambled-control χ² inflation (same code; station median fallback).
- Sites: a far station joins an existing site if < 3 km from one of its stations (single linkage); else its own site.

## 5. Kernels (new pairs)
`s11_kernel_mc.simulate` unchanged physics, source region's CVM, sensor depths as §2, 100,000 particles, variants ell20,
ell40 (primary), ell80, two5_40, ell40_s2; lapse weights from the far store grand stacks on their own τ grid; windows as
near. **Unit test |Σa − 1| ≤ 0.02 for every pair/variant/window/component; stop if a_d < 0 for > 10 % of far pairs
(ell40, 2–5 s, 3C).** Reported: far medians of a_s, a_m, a_d, a_b vs the forecast (0.23, −0.34, 0.91, 0.20), neff.

## 6. Step 3 — inversion with far pairs
Model: d = o_series + a_s S_site(t) + a_d D + **a_m M(t) + a_b B(t)**; offsets free (eliminated exactly), prior σ = 1 % on
all coefficients (conditioning only), weights 1/σ² with control inflation, σ scaled by max(1, χ²_red).
- **smooth** (Q1, Q2): S per site = 2 annual harmonics + slow basis (linear + hinges 2007/2011/2015/2019/2023; 10 dof);
  D = slow + 1 harmonic (8 dof); M and B = slow + 1 harmonic each (8 dof each), one per region.
- **temporal** (Q3): D free per episode; S, M, B as above. Q3(e) = D(e) − mean D over episodes within ±3 yr; composite =
  equal-weight mean over the region's majors (SVI the 15 of PLAN §4; WA the 32 of `wa_nw/tables/maj_list.csv`).
- Configurations (per region): **N0** near cells, M and B absent — must reproduce the stored near results (SVI §4 (c):
  composite −0.000 ± 0.041 %; WA: +0.012 ± 0.021 %; Q1/Q2 of the n ≥ 20 smooth fits) to ≤ 0.002 % — else stop and debug;
  **N1** near cells, M and B free; **F1 HEADLINE** near + far cells, M and B free; F0 near + far, M and B absent (column);
  FO far cells only (column). Kernel setups: NOM = ell40 (headline), AVG4 = mean of ell20/40/80/two5_40 (column).
- Quantities: Q1 (deep annual amplitude, peak day), Q2 (decade contrast: SVI 2015–25 − 2004–14, WA 2015–25 − 2007–14),
  Q3 per major and composite, σ_ΔD per major; M and B annual amplitudes reported.

## 7. Injections (waveform level, full raw chain, near AND far stores; injected − natural)
- Ps: common shallow pulse 0.1 % × a_s in one major (each major in turn) → median |leak| (rule ≤ 0.35); Ps(all majors) composite.
- Pd: deep pulse 0.3 % × a_d in one major → median recovery (rule 0.8–1.1).
- A: deep step 0.3 % × a_d from 2015-01-01 → Q2 recovery (0.7–1.3) and |ΔQ1| (< 0.03 %).
- B: shallow annual 0.1 % × a_s (peak day 233) → Q1 leak (< 0.25) and composite leak.
- Kernel error: A, B and Ps(all) generated with the two5_40 and with the ell80 coefficients, inverted with NOM (and AVG4).
- Mid-crust: annual 0.1 % × a_m (peak day 233) → Q1 leak (column; expected ≈ 0 with M free, reported for N0/F0).
- Each configuration N0, N1, F1 (F0 column) inverts the same injected cell values (cells common to all configurations).

## 8. Step 4 — gates on F1 (pass rules as PLAN §3–§5)
- prior ×2 / ×0.5: Q1, Q2, Q3 majors and composite move < 0.1σ.
- scrambled control: Q1_ctrl < 2σ; Q3 ≥ 90 % of majors with ≥ 5 control cells within 2σ, composite within 2σ.
- leave-one-site-out (every site incl. each far site): Q1 within 2σ (peak ±60 d where ≥ 1σ); Q3 ≥ 90 % of majors within 2σ.
  Leave-all-far-out = N1 (reported).
- family split (odd/even): Q1 G2 rule; Q3 composites within 2√(σ₁² + σ₂²) and ≥ 90 % of common majors within that bound.
- **Lapse test (far pairs):** kernel-free pooled in-phase annual amplitude ratio 3.5–6 / 2–3.5 s on the far cells (and on
  the near cells, as before), at the near pooled peak day; predictions = ratio of the mean a_d (deep) and of the mean a_s
  (shallow) over those cells. "Deep-like" iff the observation is ≥ 1σ closer to the deep prediction than to the shallow one.
- **Far/near amplitude test (new, kernel-free):** pooled in-phase annual amplitude (2–5 s) of far cells / of near cells, at
  the near peak day; predictions = (mean a_d far / mean a_d near) for deep, (mean a_s far / mean a_s near) for shallow;
  "deep supported" iff ≥ 1σ closer to the deep prediction; z to each reported.
- site dominance: no site > 50 % of Σ1/σ² (flag + leave-out guard).
- Wording: Q1 is a "deep annual change" only if Q1/σ ≥ 3, control/LOSO/family-split pass and the B leak (incl. coherent)
  < 0.25. A major's ΔD is a "deep-weighted episode change" only if |ΔD| ≥ 3σ and all Q3 gates pass. Otherwise values with
  2σ bounds.

## 9. Staging and resources
nohup queue `svi_science/scripts/s11_dist_queue.sh` (both regions; markers in `distant/<region>/`: INVENTORY_DONE,
[DIST_STOP], WINDOWS_DONE, MEASURE_DONE, KERNEL_DONE|KERNEL_STOP, INVERT_DONE, GATES_DONE; `distant/QUEUE_DONE`).
≤ 10 CPU workers, 1 BLAS thread each, RAM < 60 GB, no GPU. Restart line in `scripts/resume_after_reboot.sh`; row in
`notes/RESUME_2026-10-12_dvv_depth.md` §2. Plots: dv/v points coloured by CC or weight, no error bars, line graphs.

## Amendment 1 (Mon 2026-10-12 ~8:45 pm PT; after the inventory counts and the one-station code test, before any far dv/v,
## kernel or inversion number of the full run)
1. **Expected-cell formula.** §3 sets ρ20 = cc20. cc20 is itself the cc of a 20-window stack against a high-count reference,
   so the consistent form is ρ20 ≈ cc20² (cc(n) = √ρ(n)); the written form over-counts. Recomputed for SVI: 887 expected cells
   instead of 1,459 (WA logged in NOTES). Both ≥ 100, so the stop decision is unchanged; the written rule stays the record.
2. **Column model T+Sc** (temporal + an episode-locked REGIONAL shallow term Sc(e) × a_s, free per episode, prior 1 %), fitted
   for N1 and F1. Reason: the forecast showed per-SITE per-episode shallow terms make ΔD unidentified, but a single regional
   Sc(e) is separable from D(e) wherever the cells span different a_s/a_d ratios — which is exactly what far pairs add
   (a_s/a_d ≈ 0.2–0.3 vs ≈ 0.6 near). Reported: σ_ΔD per major and composite, deep-pulse recovery, one-major shallow-pulse
   leak and Ps(all) composite leak, prior ×2/×0.5, for N1 vs F1. No pass rule beyond §7's; it is the direct test of whether
   far pairs make a per-episode deep/shallow attribution measurable (the code test with near cells alone gave σ ≈ 0.18 %
   per major and a 0.07 shallow leak — so the near network already identifies it weakly; the question is how much far
   pairs sharpen it).
