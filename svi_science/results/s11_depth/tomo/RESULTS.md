# Static Vp/Vs pilot (LFE P and S picks) — RESULTS

2026-10-12, ~7:30–8:00 pm PT. Rules: `PLAN.md` (written before the first κ; Amendments 1–7 time-stamped). Spec:
`../MERLIN_TOMO.md` §4. Scripts `svi_science/scripts/s11_tomo_{picks,wadati,cluster,pooled,repick,svi_picks,profile,figures}.py`.

## Amendment 8 (added 2026-10-12 ~8:20 pm PT): SVI 2003–06 POLARIS-era stacks, ≤ 40 km rows, WA boreholes

**What the data show**
- **Era stacks are the cleanest κ data so far, but not because of the ≤ 40 km stations.**
  - Build: 39,522 era-detector detections (P ≥ 0.75, 2003-02 to 2006-06), 125 families, 26 stations
    (12 POLARIS, LZB/PGC/VGZ/SNB/NLLB, 9 PB boreholes from 2005), 2,696 family-station stacks in 6.6 min.
  - Clipped tP-on-tS, all stations: **κ = 1.773 ± 0.010**. 79 families, 555 rows, rms **0.096 s** (WA 0.17–0.21 s).
    Leave-one-station-out 0.007, station terms +0.012, 25.6 σ from the scrambled control. Passes every
    Amendment 3 rule.
  - Same, stations ≤ 60 km: 1.782 ± 0.021 (68 families, rms 0.079 s), passes.
  - Start-free inlier profile: 1.825 ± 0.042 for all stations (118 families), 1.845 ± 0.057 for ≤ 60 km.
    E1/E2/E3 give 1.80–1.81 and station terms ≤ 0.011. It fails only leave-one-station-out (0.09, limit 0.05).
    The profile is broad (≥ 90 % of peak over ~1.7–2.0, `figures/fig_tomo_profiles.png`).
- **Zones (clipped fit, all stations):**
  - along-strike 180–210 km: 1.782 ± 0.009, passes;
  - along-strike 150–180 km: 1.776 ± 0.014, fails station terms (+0.052);
  - depth 30–36 km: 1.767 ± 0.013, passes;
  - depth > 42 km: 1.628 ± 0.037 (9 families), passes;
  - depth 36–42 km and along-strike 120–150 km: the fit converges on a mixture (κ ~1.1–1.3, rms ~2 s), so no value.
  - The profile estimator gives 1.66–1.93 for the same zones, with leave-one-station-out 0.02–0.2.
    **So no along-strike or depth contrast survives both estimators.**
- **Per family:** still fails G1. Four families pass G1 (κ 1.74–1.79), against ≥ 50 required.
  - The scrambled control yields 0 good families, so G2b passes.
  - χ²/ν 1.38 (p 0.25): no between-family spread beyond the errors.
- **The requested ≤ 40 km (hypocentral) rows cannot carry κ.**
  - Only 13 stations and 22 families have ≥ 3 such rows; none has 6.
  - Pooled σκ 0.3–4, profile σ 0.31, leave-one-station-out 0.9; the profile peak is 46 against a scrambled-control p99 of 42.
  - Reason: within 40 km of a 30–40 km deep source, tP spans < 1 s. The precision comes from stations at
    40–60 km once their P is picked well, which the POLARIS-era stacks allow.
- **WA boreholes within 35 km:**
  - NW has only B001, B007 and B013 + B943, i.e. 3 sites. SW has B017/B018/B019 with 0 families of ≥ 3 rows.
  - NW result: κ 2.0–2.1 ± 0.16–0.2, leave-one-station-out 0.7, station terms degenerate, scrambled-control z 0.8. **Fails.**
- **Reading:** the SVI path-averaged κ (≈ 1.77–1.83) is below WA SW (1.91) and WA NW (1.8–2.0). The gap is
  ~0.1 at about 2–4 σ of the clipped fits. The two regions use different stacks, pickers' station sets and
  geometry, so this is **not** a claimed regional contrast without a like-for-like test.

Files: `tables/era_stacks/` (+ `era_stack_index.csv`), `picks_E{all,60,40}V3.csv`, `kappa_E*`, `gates_E*`,
`pooled_E60V3_inv.csv`, `pooled_EallV3_inv.csv`, `pooled_E40bh_inv.csv`, `profile_E_*`, `profile_E40bh_*`,
`picks_{BV3,CV3,B,C}bh.csv`. Scripts `s11_tomo_era_stacks.py`, `s11_tomo_era_picks.py`, `s11_tomo_subset.py`.
Earthquake-pick fetch: still NOT run (user decides).

## What the data show (first pass, before Amendment 8)

1. **Per-family κ: the pilot fails as specified.** No WA or SVI family reaches σκ ≤ 0.03 and rms ≤ 0.25 s (G1 needs
   ≥ 50). The spec set has 104 families with ≥ 6 P picks at rP ≥ 1.2, but only 17 with ≥ 6 stations where P **and**
   S both pass. Median fit rms is 1.4–1.7 s and median σκ 0.26–0.61. The fail signature (rms > 0.3 s) holds in
   every variant: clipping, ≤ 60 km, 3/5-km cluster stacks aligned on S, the new P picker, tP-on-tS. G2–G4 cannot
   be evaluated: they need ≥ 3 good families.
2. **Why.** P and S are clean only at stations within ~35 km. There the tP span for a 35–45 km deep source is only
   0.5–2 s. Co-located pick scatter is 0.04 s (P) and 0.09 s (S) robust, but 6–50 % of co-located pairs differ by
   > 0.5 s. This puts a floor of σκ ≳ 0.08 on any single family. Beyond ~40 km the frozen picker's P often
   lands on S energy arriving on the vertical (tS − tP ≈ 1.5–2.5 s at every distance, `figures/diag_sp_vs_distance.png`).
   That drives the OLS κ below 1, which is physically impossible. Stacking families within 3 km does not help:
   most clusters are singletons (NW 79/106), and the 22-member 5-km stack still shows P only within ~10 km.
3. **Zone κ does exist in WA, at about 1.9.** Fixes: a polarisation-weighted P picker (V3) plus pooling all
   families of a zone with one κ and per-family intercepts. Results (`tables/profile_*_zones.csv`, `pooled_*_inv.csv`):
   - **WA SW (81 families):** start-free inlier profile κ = **1.91 ± 0.024** (72 families, 648 rows, 43 % inliers).
     The peak is 2.6× the scrambled control's 99th percentile, leave-one-station-out shift 0.035, and the
     E1/E2/E3 estimators give 1.906–1.909. The clipped tP-on-tS fit gives 1.904 ± 0.015 (45 families, rms 0.21 s).
     With the frozen P picker it gives 1.88 ± 0.024 (38 families). **Only failed rule: station terms move κ by
     +0.034 (rule ≤ 0.02).**
   - **WA NW (141 families):** the profile is broad, flat from 1.8 to 2.1 (`figures/fig_tomo_profiles.png`).
     Peak κ = 2.00 ± 0.056 (V3 picks, LOSO 0.125 → fails) and 1.815 ± 0.047 (frozen P picks). The clipped
     tP-on-tS fit gives 1.95 ± 0.03 (V3, 72 families, rms 0.17 s, passes Amendment 3) and 1.85 (frozen P).
     The 33 eye-checked real run-1 families give 1.89–1.91 (clipped fit passes; the profile is broad, σ 0.14).
     NW κ is about 1.8–2.0; the picker choice alone moves it by ~0.15.
   - **SVI** (2012 e5r stacks, 116 families × 15 stations): no coherent κ. The profile is barely above the
     scrambled control (1.2–1.35×), peaks wander 1.6–1.9, and LOSO shifts reach 0.9. The S10 pilot grand stacks
     (3 families × 21 stations incl. POLARIS) give per-family κ 0.69/0.77/0.89. P is picked well at the
     PGC/VGZ/SILB distances (S−P 4.7–6.5 s) and wrongly at most far POLARIS stations (S−P 1.5–3.9 or 10–12 s).
     These stacks also mix detector eras, so the shared time frame is not guaranteed.
4. **No along-strike or depth structure is resolved.** Zone values differ by up to ~0.3, but no zone passes all
   gates, and the differences change with the picker (`figures/fig_tomo_zones.png`). For per-family κ vs depth,
   along-strike position and firing rate, Spearman p ≥ 0.25 everywhere except SW κ vs along-strike
   (ρ −0.41, p 0.0005, n 58). That one is not trusted: per-family σκ is 0.03–1, and the SW zone values run the
   other way (1.73 → 1.88 → 1.89). ETS segmentation cannot be tested: all WA families lie in one Brudzinski &
   Allen (2007) segment, and NW vs SW (1.8–2.0 vs 1.91) are equal within the picker spread.
5. **Physics, read plainly:** κ ≈ 1.9 is a path average from the LFE source (~35–40 km) up through the forearc
   crust. It is higher than the 1.73 the GSC locator assumes, and lower than the 2.35–2.8 receiver-function values
   for the thin low-velocity layer. A thin very-high-κ layer under a ~1.8 crust would give this kind of average.
   It does not isolate the layer.

## `sp` provenance (step 3)
`sp` in the template npz is **modelled, not measured**. It is set at `lfe_catalog_repair/scripts/e5r_templates.py:172`:
`SP = (ttS − ttP)[Crosson 1-D at the catalogue location] + (S pick − nominal) − 0.5 s`.
`bw_gate_g3_templates.py` builds the era/POLARIS templates with the same function. No P was ever picked for it.
`s10_dvv/pilot/tables/offsets_sp.csv:sp` inherits this. Using it for κ would just return the Crosson model's Vp/Vs.

## WA SW picks (step 2)
`SW_mute/hypo/picks.csv` and `SW_lag15/hypo/picks.csv` exist; they come from the older all-station runs. The
catalogue SW families (81; `stack_file` = SW_mute/SW_lag15 `hypo_28/stacks_loc.npz`) were re-picked here with
the frozen picker (`tables/picks_C.csv`) and with V3 (`picks_CV3.csv`). Sanity check: the same code on
`NW/hypo_28/stacks_loc.npz` reproduces the on-disk NW picks.csv exactly (max |Δ| 0).

## Gate table (headline rows)
| set / estimator | families | κ | σκ | LOSO | Δκ station terms | control | verdict |
|---|---|---|---|---|---|---|---|
| A spec, per family OLS | 17 measured | — | median 0.61 | — | — | — | G1 FAIL (0/50) |
| B NW, per family (any variant) | 22–95 measured | — | ≥ 0.33 median | — | — | — | G1 FAIL (≤ 1) |
| C SW, per family (any variant) | 26–63 measured | — | ≥ 0.23 median | — | — | — | G1 FAIL (≤ 2) |
| SW V3, inlier profile | 72 | 1.91 | 0.024 | 0.035 | +0.034 | 276 vs p99 107 | fails station-term rule only |
| SW V3, clipped tP-on-tS | 45 | 1.904 | 0.015 | 0.013 | +0.010 | z 13.6 | pass (Amdt 3) |
| SW frozen P, clipped tP-on-tS | 38 | 1.88 | 0.024 | 0.027 | +0.006 | z 12.5 | pass (Amdt 3) |
| NW V3, inlier profile | 130 | 2.00 | 0.056 | 0.125 | −0.014 | 361 vs 200 | FAIL (σ, LOSO) |
| NW frozen P, inlier profile | 98 | 1.815 | 0.047 | 0.07 | +0.052 | 323 vs 176 | FAIL |
| NW V3, clipped tP-on-tS | 72 | 1.95 | 0.029 | 0.019 | −0.006 | z 13.6 | pass (Amdt 3) |
| NW eye-real 33, clipped tP-on-tS | 19 | 1.909 | 0.039 | 0.012 | +0.001 | z 11.1 | pass (Amdt 3) |
| SVI 2012 V3, inlier profile | 79 | 1.645 | 0.15 | 0.945 | — | 122 vs 90 | FAIL |
| SVI pilot grand stacks, per family | 3 | 0.69–0.89 | 0.09–0.5 | 0.06–1.0 | — | — | FAIL |

Caveat on the clipped fits: they converge on a coherent subset (32–50 % of rows). In some zones of the same data
they converge on a mixture instead (rms ~2 s, κ 1.1–1.3). That dependence on the starting point is why the
start-free profile (Amendment 7) is the stricter reference.

## Step 4 — earthquake-pick fetch (PLANNED, NOT RUN; step 1 failed as specified)
The fetch is worth it only for the zone-level route, and SW WA is the one candidate. Volume, from
`g2_local/tables/eq_local.parquet` (non-duplicate rows, ≤ 50 km epicentral from a family):
- crustal (≤ 25 km): WA 6,019 at M ≥ 1 (1,805 at M ≥ 1.5); SVI 1,341 (261).
- 25–45 km: WA 1,449 (646); SVI 375 (101).
- intraslab (> Slab2 + 5 km): WA 315 (157); SVI 136 (40).
- union: 7,969 at M ≥ 1 (ComCat 6,805, NRCan 1,164); 2,595 at M ≥ 1.5; 508 at M ≥ 2. Covers 8–10 30-km zones per region.

Requests and time:
- ComCat: one event-detail call plus one `phase-data` QuakeML per event, so ~13.6k requests at M ≥ 1.
- NRCan: `fdsnws/event/1/query?eventid=…&includearrivals=true`, ~1.2k requests.
- At ~1 request/s that is ~4 h for M ≥ 1, or ~1.3 h for M ≥ 1.5. Payload ~30–80 kB per QuakeML, < 1 GB in total.

Needs the user's OK. Suggested start: M ≥ 1.5 in the SW and NW WA zones only, ~1,900 events, ~1 h.

## Files
- `tables/picks_{A,B,C,AV3,BV3,CV3,S1V3,S2V3}.csv`: picks (frozen picker; V3 = polarisation-weighted P).
- `kappa_*.csv` and `gates_*.json`: per-family κ, with gates for every variant.
- `pooled_*.csv`: Amendment 3 zone fits.
- `profile_*_{zones.csv,curves.npz,inliers.csv}`: Amendment 7.
- `colocated_pick_errors.csv`, `clusters_*`, `clstacks_*`: supporting tables.
- Figures: `fig_tomo_family.png` (per-family κ vs σκ, colour = rms), `fig_tomo_profiles.png`,
  `fig_tomo_zones.png` (colour = bootstrap spread), `fig_tomo_sp_distance.png`, `diag_*`.
