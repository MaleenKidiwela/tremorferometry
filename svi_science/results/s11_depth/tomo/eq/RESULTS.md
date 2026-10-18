# Earthquake κ and the differential Wadati — RESULTS

> **Update ~11:30 am PT (read first):** the synthetic test, the SVI calibration and the volume check below the picker-test section supersede the WA reading in §3 and the picker-test section. Near-interface earthquakes (the same start depth as the LFEs) give WA κ 1.74–1.78, and the picker is unbiased on synthetics and on SVI. So the WA LFE κ excess of +0.17 is a measurement effect specific to the LFE waveforms (not yet identified), not evidence for a fluid-rich high-κ layer.
>
> **Update 2026-10-14 ~3:05 pm PT:** the network split and the S-pick scatter correction (last section) rule out a network/instrument effect and pick scatter. Temporary-array stations alone give NW κ 1.94–1.95; the scatter-corrected all-station κ is unchanged (NW 1.95, SW 1.89). The WA excess is still unexplained.

2026-10-13. The fetch ran 9:31–10:26 am PT and the analysis 10:26–10:51 am PT; written up ~10:55 am PT.
Rules: `../PLAN.md` §6, written before any pick was fetched. Scripts: `svi_science/scripts/s11_eq_{select,fetch,parse,wadati,diff,figures}.py`,
with the runners `s11_eq_fetch.sh` and `s11_eq_analysis.sh`.
Data: `data/eq_picks/picks.parquet` holds 105,646 picks (75,089 P, 30,557 S) from 2,273 events at 849 stations; `data/eq_picks/download_log.md` has the counts.
All κ below are the clipped tP-on-tS pooled fit with one κ per zone and per-event intercepts, the same estimator as the LFE headline values.
The headline rows are stations ≤ 100 km from the event. ±σ is the larger of the event and station bootstraps.

## What the data show

1. **Earthquake κ is low and very uniform everywhere: 1.71–1.79.** Every crustal value passes every gate.
   The scrambled-control z-score is 30–72, leave-one-station-out shifts are ≤ 0.002, and station terms are ≤ 0.011.

   | region | crustal ≤ 25 km | intraslab | events at 25–45 km |
   |---|---|---|---|
   | SVI | **1.768 ± 0.005** (185 events, passes) | 1.746 ± 0.006 (65, passes) | 1.787 ± 0.007 (46, passes) |
   | WA NW | **1.747 ± 0.004** (588, passes) | 1.736 ± 0.005 (221; station-term shift 0.020, just over the rule) | 1.735 ± 0.005 (224; station terms 0.042, fails) |
   | WA SW | **1.778 ± 0.003** (336, passes) | 1.724 ± 0.006 (109; station terms 0.036, fails) | 1.711 ± 0.006 (70; station terms 0.051, fails) |

   - **Crustal κ per 30-km zone:**
     - SVI: 1.73–1.78, 4 zones, all pass.
     - WA NW: 1.73–1.77, 8 zones, 5 pass.
     - WA SW: 1.76–1.79, 6 zones, 3 pass.
   - **Intraslab κ per zone:** 1.70–1.80.
   - **Robustness:** the inlier profile, the OLS fit, ≤ 60 km, ≤ 150 km and weight > 0 picks all move the values by ≤ 0.02. One exception: at ≤ 60 km the intraslab and 25–45 km values rise to 1.75–1.81.
   - Figure: `figures/fig_eq_kappa_zones.png`.
2. **The differential (LFE κ minus crustal κ, ray-weighted in the CVM v1.7 profiles) gives κ of the layer from 25 km down to the LFE sources.**
   - The ray tracing matters. In the Wadati slope this layer carries a weight of 0.45 (SVI) to 0.53 (WA). The vertical-ray time fraction is only 0.28–0.33.
   - Whole-region values:
     - **SVI: 1.78 (16–84 % range 1.76–1.80).** Both inputs pass.
     - **WA SW: 2.02 (1.99–2.04).** Both inputs pass. With the frozen-P LFE picks it is 1.97.
     - **WA NW: 2.12 (2.07–2.17).** Both inputs pass. With the eye-checked families it is 2.05; with frozen-P picks, 1.87.
   - Per zone (full ray solution after the zone-label fix, `tables/diff_layer.csv`; 16–84 % in brackets), only two zones have both inputs passing:
     - SVI 180–210 km: **1.81 (1.79–1.83)**;
     - WA NW 60–90 km: **2.12 (2.07–2.19)**.
     - Every other zone has an LFE input that fails a gate. Their values, read with that caveat:
       - SVI 150–180 km 1.76;
       - WA NW 90–120 km 1.79;
       - WA SW −30–0 km 1.66;
       - WA SW 0–30 km 1.99 (whose crustal input also fails a gate).
     - The linear form (`tables/diff_layer_zones_linear.csv`) agrees within 0.01.
3. **The cross-checks support SVI and contradict WA.**
   - **SVI agrees with itself.** The ray model predicts 1.77 for the 25–45 km events; they give 1.787. It predicts 1.75–1.79 for intraslab events; they give 1.746. So the SVI layer κ ≈ 1.78 stands.
   - **WA does not.** If the 25 km → source layer had κ 2.0–2.1, earthquakes in that same layer should show κ 1.80–1.81, and intraslab events 1.80–1.86. Instead they give 1.71–1.74, lower than the crust.
   - The cleanest case is WA NW zone 60–90 km, where every input passes all gates:
     - LFE κ 1.95;
     - crustal earthquakes 1.768;
     - earthquakes at 25–45 km 1.766;
     - intraslab 1.751.
   - Earthquakes sampling the LFE depth range therefore see no high-κ layer in WA. The ~0.15–0.2 excess in WA LFE κ is more likely in the LFE picks than in the rock.
     - The LFE κ already moves by 0.15 with the P picker (V3 vs frozen). The frozen-P values (1.815–1.88) give a layer κ of 1.87–1.97, closer to the earthquakes.
4. **Compared with receiver functions (Audet/Bostock κ 2.35–2.8): no body-wave path here sees such high κ.**
   - **SVI:** the layer κ 1.78 sits just above the CVM value for the same layer (1.70). It would allow a κ_RF layer only 0.7–1.2 km thick inside 25–35 km.
   - **WA:** the face-value LFE differential (2.0–2.1) would need a 3.4–8.8 km κ_RF layer. The earthquakes do not support it.
   - The likely reason is geometry, not a contradiction. LFE rays start at the plate interface and go up, so they barely cross a low-velocity layer at or below the sources. The earthquakes we used are above it (25–45 km) or cross it only steeply (intraslab events).
5. **Along-strike contrast.**
   - **Crust:** WA SW 1.778 > SVI 1.768 > WA NW 1.747. The SW − NW difference is 0.031 ± 0.005; all three crustal fits pass their gates.
   - **Deep paths:** SVI's 25–45 km events give 1.787, above WA's 1.71–1.74, i.e. 0.05–0.08 higher. The WA fits fail the station-term rule, so this is a candidate, not a result.
   - The deep-path contrast runs opposite to the LFE ordering (WA 1.9 > SVI 1.77). That is another sign that the WA LFE κ is biased high.

## Picker test (§6 Amendment 1b; run ~10:58–11:10 am PT): does the LFE picker push κ up?

**What the data show.** On earthquake records the LFE picker gives κ only about **+0.06 higher** than the catalogue picks of the same events. That is a real but small bias, about one third of the WA LFE excess. It does not explain the whole gap.

**Method.** The frozen LFE picker and its V3 variant were run on 1–8 Hz, 40 Hz records from the same cache the LFE stacks came from. Window: origin − 4 s to + 40 s, at the LFE stations.
- Events: 89 deep earthquakes (25–45 km plus intraslab) on days that are in the cache (64 in WA NW, 25 in WA SW, mostly 2006–2011 ETS days).
- Only stations with rP, rS ≥ 1.2 are used, and a station is kept only if it is ≤ 100 km from the event.
- Most LFE stations are temporary-array stations with no catalogue picks, so the comparison is made on the same events but different stations: the catalogue κ comes from those events' own picks at stations ≤ 100 km.
- Both fits use the same clipped tP-on-tS estimator. Files: `tables/pickertest_rows.csv`, `pickertest_same_events.csv`, `pickertest_kappa.csv`, `pickertest_bias.csv`.

| region | records | picker | κ picker | κ catalogue, same events | bias |
|---|---|---|---|---|---|
| WA NW | raw | V3 P | 1.798 ± 0.027 (29 ev) | 1.736 ± 0.008 | **+0.06 ± 0.03** |
| WA NW | raw | frozen P | 1.795 ± 0.022 (23 ev; passes all gates) | 1.732 ± 0.008 | **+0.06 ± 0.02** |
| WA SW | raw | V3 / frozen | 1.815 / 1.823 ± 0.045 (8 / 7 ev) | 1.750 ± 0.022 | +0.07 ± 0.05 |
| WA NW | noise-matched | V3 / frozen | 1.688 ± 0.054 / 1.760 ± 0.060 (18 / 11 ev) | 1.735 / 1.734 | −0.05 ± 0.06 / +0.03 ± 0.06 |
| WA SW | noise-matched | V3 | 1.91 ± 0.14 (6 ev) | 1.735 ± 0.021 | +0.17 ± 0.14 |

**Pick offsets (picker minus catalogue, rows where both exist).**
- P is late by a median of 0.26–0.55 s, and 32–58 % of P picks are more than 0.5 s off.
- S is late by 0.13–0.31 s, and ≤ 8 % of S picks are more than 0.5 s off.
- The extremum picker marks the peak, not the onset, and its P is the noisy half.

**Reading against the pre-written rule.**
- On clean records the bias is +0.06–0.07 in both regions (2–3σ in WA NW). That is the "partly" case: the picker can explain about +0.06 of the WA LFE excess.
- The noise-matched headline settles nothing: NW −0.05 ± 0.06, SW +0.17 ± 0.14 on 6 events. Adding noise leaves too few events and too much scatter.
- **Corrected layer κ.** Removing +0.06 from the LFE κ:
  - WA NW: LFE 1.89, layer κ ≈ 2.02;
  - WA SW: LFE 1.84, layer κ ≈ 1.91.
  These are still above the 1.71–1.74 that earthquakes in the same layer show. So about two thirds of the WA LFE excess remains unexplained by picking on earthquake records.

**Obstacles.**
- Earthquake records are impulsive; LFE stacks are emergent and low-frequency, so a picker bias could be larger on LFEs than this test can show.
- The station sets differ: the picker runs on the LFE (mostly temporary-array) stations, the catalogue κ on PNSN stations.

**Next adaptation.**
- Put a synthetic LFE-like source (an emergent 2–4 Hz pulse) on the real station noise, with known κ and the LFE geometry, and measure the picker's κ there.
- Alternatively, repeat this test with SVI POLARIS-era earthquakes. There the LFE κ (1.773) already matches the earthquakes, so the picker bias should come out near 0. That would calibrate the method.


## Synthetic test, SVI calibration and volume check (§6 Amendment 2; run ~11:16–11:30 am PT)

**What the data show.**
- The LFE picker does not manufacture a high κ. On synthetic emergent pulses at the real WA geometry and noise levels it recovers the true κ within −0.02 to +0.002.
- On SVI earthquakes the picker matches the catalogue exactly (bias 0.00).
- Earthquakes at the LFE source level show WA κ 1.74–1.78. That is 0.17 below the WA LFE κ.
- So the WA LFE excess is not explained by picker mechanics, and it is not seen by earthquakes starting at the same place. It is specific to the LFE waveforms.

**1. Synthetic emergent LFE-like records** (`tables/synth_kappa.csv`, `synth_rows.csv`; script `s11_eq_synth.py`).
- Geometry: 3,902 real WA (family, station) rows at 25–50 km depth. P times come from the CVM v1.7 ray model, with tS = κ·tP.
- Pulses: 3 Hz carrier with a t²e^(−t/0.25 s) envelope. S is 4× P, with S leaking onto Z at 0.3× (and a no-leak column).
- Noise: real cached station noise, matched to each row's LFE-stack SNR. 3 noise draws plus a noise-free run.
- Pickers and fit: frozen and V3 pickers, rP, rS ≥ 1.2, the same clipped fit.

| true κ | recovered, WA NW (V3 / frozen; LFE noise, S leak 0.3) | recovered, WA SW |
|---|---|---|
| 1.75 | 1.750 / 1.748 | 1.744 / 1.746 |
| 1.90 | 1.892 / 1.891 | 1.895 / 1.896 |
| 2.10 | 2.095 / 2.094 | 2.093 / 2.087 |

- Bias is −0.021 to +0.002 in every case (σ 0.004–0.007; 900–1,000 rows in NW, 550–660 in SW).
- The picker marks the wavelet peak 0.5 s after onset. That offset is common to P and S, so it cancels in the slope.
- About 75 % of noisy P picks are more than 0.5 s off the onset, yet κ is still recovered: the clipping and the rP, rS ≥ 1.2 cut handle it.
- **Rule outcome:** |bias| ≤ 0.05 at all three κ. The picker cannot turn 1.75 into 1.9.

**2. SVI picker test (calibration)** (`tables/pickertest_svi_rows.csv`, `pickertest_svi_same_events.csv`).
- Data: 64 deep SVI earthquakes from 2003–2013, records from the 20-Hz backward cache.
- Raw records: picker 1.773 ± 0.021 (V3) and 1.780 ± 0.016 (frozen, passes all gates). The catalogue gives 1.774 / 1.777 for the same events, so the **bias is 0.00 ± 0.02**.
- Noise-matched, frozen picker: +0.01 ± 0.03. The V3 fit collapsed onto a mixture (κ 1.14, rms 2.5 s) and gives no value.
- The WA raw-record +0.06 is therefore not a general picker property. It may come from the WA station set (mostly temporary-array stations) or from WA waveforms.

**3. Which rock do the "earthquakes in the layer" sample?** (`tables/volume_split_kappa.csv`)
- Most of them are overriding-plate events, not source-level events. In WA NW 170 of 247 are a median 22 km above Slab2 at ~27 km depth; the rest sit at the source level (table below).

| subset | SVI | WA NW | WA SW |
|---|---|---|---|
| overriding plate, 17–22 km above Slab2, ~26 km deep | 1.806 ± 0.005 (station terms fail) | 1.725 ± 0.005 (station terms fail) | 1.710 ± 0.007 (station terms fail) |
| **near interface (within 5 km of Slab2), ~36–41 km = LFE source level** | **1.755 ± 0.011 (passes)** | **1.775 ± 0.008 (passes; 68 events)** | **1.735 ± 0.015 (passes; 11 events)** |
| intraslab 0–15 km below Slab2 | 1.743 ± 0.006 (passes) | 1.736 ± 0.005 (passes) | 1.725 ± 0.006 (station terms fail) |

- The fair comparison for the LFEs is the near-interface events. Their rays start where the LFE rays start.
  - WA NW: LFE 1.95 ± 0.03 vs near-interface earthquakes 1.775 ± 0.008, a difference of **+0.17 ± 0.03**.
  - WA SW: LFE 1.904 ± 0.015 vs 1.735 ± 0.015, a difference of **+0.17 ± 0.02**.
  - SVI: LFE 1.773 vs 1.755, a difference of **+0.02 ± 0.015**.

**Verdict on the WA NW layer κ ≈ 2.0: not supported as a property of the rock. It is most likely a measurement effect specific to the LFE waveforms, and is still unexplained.**
- The fault-valve, fluid-rich reading is not supported, for four reasons:
  1. Earthquakes starting at the same level and travelling the same upward paths give κ 1.775, which implies a 25 km → source layer κ of about 1.80.
  2. A high-κ zone confined near the sources would barely change the Wadati slope. All rays share the near-source path, so the slope's weight comes from the 25 km → source layer as a whole.
  3. The picker cannot make the excess. On synthetics the bias is ≤ 0.02, and on SVI earthquakes it is 0.00.
  4. SVI LFEs agree with their earthquakes to within 0.02.
- **What remains.**
  - A waveform effect not in the synthetics. The leading candidate is LFE S envelopes that broaden with distance: at far stations the extremum picker lands later on a long, scattered S wavetrain, which inflates tS − tP with distance.
  - Or a WA station-set effect, consistent with the +0.06 seen on WA earthquake records but not on SVI ones.
  - A thin, very-high-κ zone right at the WA sources (the receiver-function layer) is still allowed, because it is nearly invisible to the slope; it just cannot be the explanation.
- **Next (cheap, existing data).**
  - Measure the S envelope width vs distance on the WA LFE stacks and on the WA earthquake records, and feed the measured broadening into the synthetics.
  - Re-fit the WA LFE κ on stations ≤ 40–50 km only, where broadening is small. If it drops toward 1.78, the excess is a measurement effect.

## The obstacle and the next adaptation (WA) — written before the picker test above

The WA differential cannot be trusted because the LFE κ and the earthquake κ do not come from like-for-like measurements: different pickers, stacked vs catalogue picks, and different station sets.

Next step (cheap; needs no new download):
- Run the LFE picker (V3, the frozen picker) on waveform stacks of the 25–45 km earthquakes located within 10 km of LFE families (WA NW zone 60–90 km has 43).
- Use the same stations and ≤ 60 km rows as the LFE fit, and compare with the catalogue-pick κ of the same events.
- If the picker returns ~1.9 where the catalogue picks give ~1.77, the WA LFE κ is picker bias and the corrected layer κ ≈ crust κ.
- If it returns ~1.77, the excess is real and sits in the LFE source region itself, not on the earthquake paths.

## Caveats
- **Per-zone differential, FIXED ~10:57 am PT:** `s11_eq_diff.py` had matched zones by label, and the labels differed ("+30" vs "30"), so only WA SW −30 had a per-zone row in v1. The label is now formatted as "+30". `tables/diff_layer.csv` and `figures/fig_eq_layer_kappa.png` (zone squares now drawn) were rerun; the whole-region values are unchanged.
- **Coverage:** NRCan-only events (419, mostly SVI) have no public picks. SVI therefore rests on PNSN-located events: 185 crustal and 65 intraslab, which is still enough.
- **Station terms:** the WA 25–45 km and intraslab fits fail the station-term rule (shifts 0.020–0.051). Their κ stays 1.71–1.75 under every variant.
- **Test files:** `tables/test_partial_*` and `figures/test_partial_*` are a code smoke test on 66 events, not results.

## Files
- `tables/eq_kappa_zones.csv`: every cell × estimator (clip_inv, clip_ols, profile_0.25, profile_0.5) × variant (d100, d60, d150, d100_w), with gates.
- `tables/eq_kappa_rows.csv`: 27.9k event–station P+S rows.
- `tables/diff_layer.csv`: whole-region ray-weighted layer κ for every LFE × earthquake estimator, plus the simple-layer column and the κ_RF thickness.
- `tables/diff_layer_zones_linear.csv`: per-zone layer κ (linear form).
- `tables/diff_crosscheck.csv`: predicted vs observed κ for intraslab and 25–45 km events.
- `tables/diff_model.json`: ray-model geometry and layer weights.
- `selection.csv`, `selection_counts.csv`: the event selection.
- Figures: `figures/fig_eq_kappa_zones.png` and `figures/fig_eq_layer_kappa.png`.

## Envelope broadening with distance and the distance-limited κ (§6 Amendments 3 + 4; run 10-13 ~11:32–11:42 am and 10-14 12:34–12:50 pm PT)

**What the data show.**
- LFE S envelopes do get wider with distance, but just as much in SVI, where the LFE κ matches the earthquakes. The broadening is not specific to WA.
- Feeding the measured broadening into the synthetics gives a κ bias of only **+0.002 to +0.03**, far below the +0.17 WA excess.
- Dropping the far stations does not lower the WA κ; it raises it (≤ 50 km: NW 2.02, SW 1.98; ≤ 40 km: 2.0–2.26).
- So the WA LFE excess is **not** a far-station broadening / peak-picking effect. It is still unexplained.

**1. Width vs epicentral distance** (`tables/swidth2_slopes.csv`; script `s11_eq_swidth2.py`).
- Method: the same 1–8 Hz records and picks as the κ fits. Width = full width at half maximum (FWHM) of the Hilbert envelope at the picker's S (horizontals) or P (vertical) peak. Slope by OLS, ± bootstrap over families/events. Scrambled control: distances shuffled among the rows of each family/event (1,000), p = how often the shuffled slope is as far out.
- LFE stacks, slopes in ms per km (median width ≤ 30 km → ≥ 60 km):

  | | S width | P width | S − P width | S rise time (onset → peak) |
  |---|---|---|---|---|
  | WA NW (1,319 S rows, 140 fam.) | 1.2 ± 0.5 (0.73 → 0.85 s; scrambled p 0.30) | 0.0 ± 0.6 | 2.4 ± 0.9 | 1.0 ± 0.3 |
  | WA SW (942, 81) | 3.9 ± 0.5 (0.50 → 0.78 s; p < 0.001) | 1.2 ± 0.3 | 2.0 ± 0.7 | 2.3 ± 0.3 |
  | SVI (1,773, 122) | 3.7 ± 0.4 (0.48 → 0.73 s; p < 0.001) | 1.7 ± 0.3 | 2.1 ± 0.5 | 2.0 ± 0.2 |

  - 2–8 Hz gives the same picture (S 1.3 / 3.4 / 3.0; S − P 2.8 / 2.1 / 2.1 ms/km).
  - The quantity that could bias κ is the S − P difference. It is 2.0–2.4 ms/km in all three regions, equal within 1σ.
  - Amendment 3 rule ("a finding" only if WA > SVI at 2σ): **not met.** WA SW − SVI S slope = 0.2 ± 0.6; WA NW is below SVI.
  - The picker sits on the envelope peak: pick − peak median 0.00 s, slope ≤ 0.3 ms/km.
- Earthquake records at the same stations (461 WA rows, 484 SVI rows, deep events):
  - S widths are narrower: 0.33–0.48 s (LFE stacks 0.48–0.85 s).
  - S slopes are flat to weak: NW −0.2 ± 0.6, SW 1.6 ± 1.1, SVI 1.0 ± 0.4 ms/km. S − P slopes are ≤ 0 (−2.4 ± 1.2, −0.9 ± 1.1, −0.5 ± 1.1).
  - So LFE stacks broaden more than earthquakes. That fits stacking jitter and emergent sources, and it happens in SVI too.
  - Picker − catalogue S delay vs distance: NW 8 ± 8 ms/km (78 rows), SVI 4.5 ± 1.3 (112 rows); SW has 9 rows only. The SVI picker κ is still unbiased (0.00 ± 0.02, Amendment 2).

**2. Synthetics with the measured broadening** (`tables/synthw_kappa.csv`, `synthw_closure.csv`; script `s11_eq_synth_width.py`).
- Same rows, ray geometry, real noise, SNR matching and fit as Amendment 2b. The pulse time constant is set from the measured FWHM (P and S separately; FWHM ≈ 3.5 τ).
- Bias = recovered − true κ (true 1.75; mean of 3 noise draws; V3 / frozen P picker):

  | variant | WA NW | WA SW |
  |---|---|---|
  | W0 constant width (repeat of Amendment 2b) | +0.000 / +0.002 | −0.009 / −0.012 |
  | **W1 measured trend (headline)** | **+0.002 / +0.007** | **+0.014 / +0.015** |
  | W2 each row's own measured width | +0.020 / +0.024 | +0.019 / +0.022 |
  | W3 3× the trend slope | +0.016 / +0.019 | +0.027 / +0.030 |
  | W2s scrambled control (widths shuffled across rows) | −0.001 / −0.001 | +0.004 / +0.007 |

  - Bootstrap σ 0.003–0.008. True κ 1.90 gives W1 bias +0.004 to +0.011. Noise-free W1: +0.007 (NW), +0.012 (SW).
  - Closure: noise flattens the measured width trend in the synthetics. W1 synthetics show S slopes 0.1 (NW) / 1.5 (SW) ms/km, and W3 shows 2.1 / 7.1. The observed 1.2 / 3.9 lie between them, so the realistic bias is between the W1 and W3 values: **+0.002 to +0.03**.
  - Rule (|W1 bias| ≤ 0.05 → broadening does not explain the excess): **met.** Broadening cannot turn 1.75 into 1.9.

**3. LFE κ from near stations only** (`tables/swidth2_dist_kappa.csv`, `swidth_kappa.csv`; epicentral distance; families with ≥ 3 rows, ≥ 4 in brackets).

| | all stations | ≤ 50 km | ≤ 40 km | scrambled control at ≤ 50 km (median, 16–84 %) |
|---|---|---|---|---|
| WA NW | 1.950 ± 0.029 (passes) | **2.022 ± 0.048** (passes) [1.957 ± 0.044, passes] | **2.234 ± 0.095** (fails σ, LOSO) [2.258] | 1.957 (1.87–1.98); real at 96th pct [53rd] |
| WA SW | 1.904 ± 0.016 (passes) | **1.977 ± 0.064** (station terms −0.06, fails) [1.972] | **2.047 ± 0.119** (fails) [2.013] | 1.895 (1.88–1.91); real at 100th pct |
| SVI | 1.773 ± 0.010 (passes) | 1.785 ± 0.016 (passes) [1.783] | 1.809 ± 0.026 (passes) [1.800] | 1.772 (1.76–1.79); real at 82nd pct |

- The scrambled control draws random station subsets of the same size per family (distances shuffled within family; 300 draws). Random subsets give the all-station κ. Near-station subsets give the same or **higher** κ, never lower.
- Hypocentral ≤ 60 km: NW 1.98 ± 0.05, SW 1.98 ± 0.045 (station terms fail), SVI 1.79.
- Amendment 3 rule (drops to within 0.05 of the near-interface earthquakes, NW 1.775 / SW 1.735): **not met in either WA region.** The ≤ 50 km values stay ≥ 1.96.
- The rise at short range occurs in SVI too (+0.01–0.04). A likely, untested reason: near stations have a short tS lever arm, and S-pick scatter then inflates the tP-on-tS κ. WA S picks scatter more than SVI ones (rms 0.15–0.20 s vs 0.08 s), so the rise is larger there. This works against, not for, the broadening hypothesis.
- Half-max onset picks instead of peaks (Amendment 3b): NW 1.83 (all; station terms 0.029, fails) and 1.85–1.87 (≤ 50 km, passes); SW 1.87–1.88; SVI 1.74–1.75. NW falls between the two Amendment 3 thresholds ("partly"); SW stays 0.14 above the earthquakes.
  - The NW drop with onset picks (0.12) is larger than the rise-time trend can produce (S − P rise-time slope ≈ 0.8 ms/km → ≤ 0.03 by the synthetics). So it is not the distance trend; which rows survive the clipping changes too.

**Verdict.** The far-station broadening is real, but it is the same in SVI and gives only +0.002 to +0.03 of bias. Near-station κ does not drop. **The WA excess is not a broadening or peak-picking measurement effect.**
- What remains open:
  - a WA station-set effect: the +0.06 seen on WA earthquake records (Amendment 1b), absent in SVI;
  - near-station S-pick scatter;
  - a real property of the rock near the WA sources.
- Next (cheap, existing data):
  - Refit WA LFE κ on PB boreholes plus permanent stations vs temporary-array stations separately. This tests the station set.
  - Apply the E3 S-noise correction (Amendment 7) to the ≤ 50 km subsets.

Figure: `figures/fig_swidth_kappa.png` (width vs distance for LFE stacks and earthquakes; κ by distance limit with the scrambled control; synthetic bias per variant). Rows: `tables/swidth2_lfe_rows.csv`, `swidth2_eq_rows.csv`, `synthw_rows.csv`; τ lookup `synthw_tau_lookup.csv`, trend `synthw_trend.csv`. Decision log: `NOTES.md`.

## Network split, station residuals, B013 and the S-pick scatter correction (§6 Amendment 5; run 10-14 ~2:48–2:55 pm PT)

**What the data show.**
- **The excess is in the temporary-array data on its own.** Temporary stations alone (CAFE XU + TA + FACES YW) give WA NW
  κ 1.94–1.95, and XU alone gives 1.95–1.96. In WA SW the temporary stations are essentially the whole data set (1.90–1.91).
  The permanent stations cannot be fitted on their own with the clipped fit. Their start-free profile κ is 2.04, not lower.
  So no network carries the excess, and neither does a timing, orientation or response problem in one network.
- **S-pick scatter does not explain it.** The scatter-corrected estimator (instrumental variable, IV) gives the same
  all-station κ as before: NW 1.949 ± 0.029, SW 1.891 ± 0.019, SVI 1.769 ± 0.011. Synthetics show that scatter of the LFE size
  raises the near-station κ by only ≤ 0.02. The real data rise by 0.14–0.28.
- **The near-station rise in WA NW persists after correction.** IV at ≤ 40 km gives 2.15 ± 0.10, which is +0.20 ± 0.06 above
  all stations. That fit still fails the σ and leave-one-station-out gates, and it rests on a few stations. In SW the rise
  shrinks to +0.075 ± 0.043 (1.7σ). Most of it comes from one station, S020.
- **B013 is irrelevant.** All WA stack data date from 2007-01 to 2011-12, before both B013 instrument periods. A 5 ms
  B013 step moves κ by 0.0002.

Script `svi_science/scripts/s11_eq_netsplit.py` (3 workers), figures `s11_eq_netsplit_fig.py`. Rows = the V3 LFE picks of
Amendments 3–4 (`swidth2_lfe_rows.csv`, rP, rS ≥ 1.2). Fit = the same clipped pooled tP-on-tS (E1) and gates as before.
Sanity check passed: all-station NW 1.950, SW 1.904. Classes: permanent = UW, US, PB; temporary = XU, TA, YW.
No XG (Array of Arrays) station is in the 25-station NW/SW pools, so XG could not be tested.

**1. κ by station class** (`tables/netsplit_class_kappa.csv`; families with ≥ 4 class rows, ≥ 3 in brackets).
Profile = the start-free inlier-profile κ (Amendment 7), shown as a column.

| | stations | clipped κ (E1) | profile κ | gates |
|---|---|---|---|---|
| NW all | 21 | 1.950 ± 0.028 [1.936] | 2.000 | pass |
| NW temporary (XU, TA, YW) | 12 | 1.941 ± 0.033 [**1.950 ± 0.024**] | 1.905 | fail at ≥ 4 (station terms −0.020); [pass] |
| NW XU only | 10 | 1.951 ± 0.038 [1.956 ± 0.029] | 2.02–2.03 | [pass] |
| NW permanent (UW + PB) | 9 | 1.21 ± 0.10: **collapsed fit** (rms 2.2 s, unphysical) | **2.04** | fail |
| NW PB only / UW only | 4 / 5 | 1.29–1.43 / 1.09 (5–20 families, collapsed) | 1.66–2.04 / 1.42 | fail |
| NW temporary ≤ 50 km | 11 | [1.987 ± 0.074] | 1.95 | fail (σ, LOSO) |
| SW temporary | 20 | 1.907 ± 0.014 [1.897] | 1.895 | pass |
| SW permanent (NLWA, WISH, SQM, B017–B019) | 6, median 58 km | no family has ≥ 3 rows | — | — |

- **Scrambled control:** class labels were shuffled across stations (300 draws). Half-networks drawn at random also often
  collapse: their median κ is 1.39 / 1.47. The real NW Δκ (permanent − temporary) = −0.73 sits at the 6th percentile,
  inside the control's 2.5–97.5 % range. Shuffling within distance terciles gives the same result. So the permanent-class
  collapse is what a random 9-station subset of this network does. It is not a property of the permanent stations.
- **Common families** (≥ 3 rows in both classes, 23–46 families):
  - The temporary rows give 1.75 ± 0.05–0.06 (profile 1.75–1.85).
  - That fit fails leave-one-station-out (0.054–0.066) and station terms (0.040–0.056). It is a lead, not a result.
  - The permanent rows of the same families collapse (profile 1.99).
- **Earthquake records, same picker, per class** (`tables/netsplit_eq_class_kappa.csv`):
  - NW temporary stations: 1.83 ± 0.05 (V3 picker) and 1.78 ± 0.04 (frozen picker).
  - NW permanent stations: not measurable (σ 1.6).
  - Catalogue picks at the permanent stations of the same events: 1.815 ± 0.033.
  - Earthquakes show no class difference.
- **Rule outcome (Amendment 5a):** neither class is within 0.05 of the earthquake value. The temporary class alone gives
  ≥ 1.88, and the permanent class gives profile κ 2.04 (its clipped fit failed). **The excess is not a network effect.**

**2. Per-station residuals of tS − tP vs the SVI-calibrated model** (κ 1.773 × CVM v1.7 ray tP from the catalogue location;
family-demeaned; `tables/netsplit_station_residuals.csv`, `netsplit_residual_class.csv`; figure `figures/fig_netsplit_residual_map.png`).
- **Residuals grow with distance in both classes alike.** On the NW kept rows, the 10-km medians of the PB stations and the
  temporary stations agree within ~0.1 s out to 55 km. Beyond ~45 km the temporary medians rise to +0.1 to +0.3 s, and in
  SW to +0.5 s at 60–70 km. This is what κ above 1.773 looks like. It is not a step between networks.
- **Stations with the largest positive median residual** (more S − P than the model):
  - NW: SQM +0.84 (UW, 8 rows), PL11 +0.54 (XU), GNW +0.43 (UW), B001 +0.36 (PB), BS11 +0.33 (XU), C04A +0.33 (XU).
  - SW: a western cluster of XU stations, S015 +1.02, N015 +0.94, S020 +0.75, N030 +0.64, plus E03A +0.70 (TA) and S055 +0.42.
  - The positive stations are mixed between permanent and temporary in NW. In SW almost all stations are temporary anyway.
- **Moveout of each phase against the model** (within-family slope on model tP; both slopes are lowered equally by location
  error, so only their ratio is robust; column only):

  | | P slope | S slope |
  |---|---|---|
  | NW | 0.80 | 1.55 |
  | SW | 0.95 | 1.80 |
  | SVI | 0.96 | 1.71 |

  - Against SVI, NW P moveout is smaller, while SW S moveout is larger.
  - This is a hint about which phase drives the excess in each region, not a result.

**3. B013** (`tables/netsplit_b013.csv`).
- **The data predate both B013 instrument periods.** The stack data cover 2007-01-13 to 2011-12-14 (NW) and to 2008-07-17 (SW).
  The ballistic round-3 S step (after 2015) and the S14 B013 CC period (2012-08 to 2013-05) come later.
- **Shifting B013 tS** by −5 ms moves κ by +0.0001 (all stations) and +0.0002 (≤ 50 km). Shifting it by −50 ms gives −0.027 and
  +0.002; this changes which rows the clipping keeps.
- **Co-located picks on the same families:**
  - B013 − B943: ΔtP 0.00 s (robust SD 0.04), ΔtS 0.00 s (0.07), 54 families.
  - B013 − W030: ΔtP +0.10 s (0.04), ΔtS 0.00 (0.19). The W030 P is 0.1 s earlier. This is a constant site offset, which the
    family intercepts and station terms absorb.
- **Dropping B013 collapses the clipped fit.** Without B013, or without the whole B013/B943/W030 site, the fit collapses
  (rms 1.7–2.1 s, κ 1.30–1.41).
  - The cause: this one site supplies 40 % of the kept NW rows (152 of 382) and anchors the clipping.
  - The start-free profile does not move without B013 (2.000). Without the whole site it is 1.905 (all stations) and 1.94
    (≤ 50 km, down from 2.105).
- **Rule outcome (Amendment 5c):** the pre-registered rule is triggered only by a collapsed fit. The step itself changes
  κ by ≤ 0.0002. The three co-located instruments agree, so the ~0.1 pull of the site on the profile κ comes from that site's
  path or geometry, not from an instrument.

**4. S-pick scatter correction for near stations** (`tables/netsplit_scatter_kappa.csv`; IV = instrumental variable with the
model tP as instrument. Pick errors are independent of geometry, so IV needs no σ_S).

| | E1 (before) | **IV (headline)** | E3 / SIMEX / Deming | IV rise vs all stations | E1 rise | random same-size subsets (IV median; real pct) |
|---|---|---|---|---|---|---|
| NW all | 1.950 | **1.949 ± 0.029** | 1.946 / 1.945 / 1.945 | — | — | — |
| NW ≤ 50 km | 2.022 | **2.015 ± 0.049** | 2.006 / 2.008 / 2.007 | +0.066 ± 0.034 (≥ 4 rows: +0.009) | +0.071 | 1.948; 94th |
| NW ≤ 40 km | 2.234 | **2.151 ± 0.101** (LOSO 0.062, B001) | 2.205 / 2.208 / 2.212 | **+0.202 ± 0.064** | +0.283 | 1.847; 100th |
| SW all | 1.904 | **1.891 ± 0.019** | 1.903 / 1.904 / 1.894 | — | — | — |
| SW ≤ 50 km | 1.977 | **1.948 ± 0.058** | 1.973 / 1.973 / 1.940 | +0.058 ± 0.035 | +0.073 | 1.885; 100th |
| SW ≤ 40 km | 2.047 | **1.965 ± 0.109** | 2.040 / 2.042 / 1.989 | +0.075 ± 0.043 | +0.143 | 1.888; 100th |
| SVI all | 1.773 | **1.769 ± 0.011** | 1.770 / 1.768 / 1.767 | — | — | — |
| SVI ≤ 40 km | 1.809 | **1.793 ± 0.028** | 1.797 / 1.797 / 1.789 | +0.024 ± 0.016 | +0.036 | 1.770; 88th |

- **E3, SIMEX and Deming hardly move.** The measured S-pick scatter (σ_S 0.13 s at the NW co-located site) is small next to the
  1.4–2.6 s within-family tS spread at ≤ 40–50 km, so R_S = 0.99. Doubling σ_S gives NW ≤ 40 km 2.12–2.15.
- **Synthetic closure** (`tables/netsplit_synth_kappa.csv`): the synthetics use the real WA geometry, LFE noise and pickers,
  and the instrument carries 5-km location errors.
  - With κ_true 1.75, E1 rises at ≤ 40 km by only +0.018 (NW) and −0.009 (SW). The real E1 rise is +0.283 and +0.143.
  - IV recovers the true κ within 0.014 at every distance cut and every κ_true (1.75, 1.90, 2.10).
- **Rule outcomes (Amendment 5d):**
  - The NW ≤ 40 km rise **persists** (+0.20 > 2σ and > 0.05). At ≤ 50 km it counts as **removed** (1.9σ).
  - SW: **removed** at both cuts (≤ 1.7σ).
  - Synthetic closure: scatter **does not explain** the E1 rise (< 10 % of it).
  - The corrected all-station κ (NW 1.949, SW 1.891) is still 0.16–0.17 above the near-interface earthquakes. **The
    correction does not remove the excess.**
- **Post hoc, which stations make the near rise** (`tables/netsplit_posthoc_drop.csv`; checked after seeing the leave-one-station-out stations):
  - SW without S020: IV 1.83 (≤ 50 km) and 1.84 (≤ 40 km), i.e. no rise.
  - NW without B001: the rise grows (IV ≤ 40 km 2.24).
  - NW without the B013 site: the clipped fit collapses, and the profile gives 2.29 at ≤ 40 km.
  - So the NW near-station rise is not one station. It is carried by a 64-family, 278-row fit that still fails the gates.

**Verdict.**
- **What the WA excess is not:**
  - a network effect: temporary arrays alone give 1.94–1.95, and the permanent profile κ is 2.04;
  - a B013 instrument effect: the data predate both B013 periods;
  - S-pick scatter: the IV-corrected κ is unchanged, and synthetics show ≤ 0.02 of scatter bias.
- **What remains:** the near-station rise in NW (+0.20 ± 0.06 at ≤ 40 km after correction), which is real but fails the
  σ and leave-one-station-out gates. Taken with the earlier tests (picker, synthetics, broadening), every measurement check
  on the WA LFE picks has now come back negative.
- **What is left as an explanation:** a property of the LFE sources or the rock near them that earthquakes starting at the
  same level do not sample, or an LFE-waveform effect not in any synthetic so far.
- **Lead (column, not a result):** against the model, the NW excess comes with a smaller P moveout, and the SW excess with a
  larger S moveout.
- **Cheapest next tests:**
  1. A per-phase moveout check with earthquake travel times (no model), to see whether WA LFE P or S times are the odd ones.
  2. The common-family temporary κ of 1.75, refit with the profile estimator over more families, as a pre-registered test.

Figures: `figures/fig_netsplit_kappa.png` (κ by class with the profile κ; scrambled control; near-station E1 vs IV/E3/SIMEX;
synthetic bias) and `figures/fig_netsplit_residual_map.png` (per-station residual maps and residual vs distance by class).
