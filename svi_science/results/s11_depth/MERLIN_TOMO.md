# MERLIN_TOMO — LFE-based Vp/Vs imaging: feasibility and pilot (2026-10-13, Fable 5.1; saved by main session)

**Verdict.** Viable in a coarse form; 3-D tomoDD is premature. The decision that matters: start with per-family Wadati Vp/Vs (κ). It needs no location, no velocity model and no new software. Let its scrambled-control separation decide whether DD tomography is worth building.

## 1. Picks on disk (verified)
- **WA NW:** `lfe_array_discovery/results/ab_royer/NW/hypo_28/picks.csv`, 169 groups × 25 stations, 4225 tP / 3549 tS from Hilbert-rotated stacks (`scripts/a4f_hypo.py`).
  - Quality is weak: unambiguity ratio rP median 1.10.
  - At rP ≥ 1.2: 146 families have ≥ 4 P picks (104 have ≥ 6). At ≥ 1.5, only 26.
  - S: 61 families have ≥ 4 picks at 1.2.
  - The SW picks.csv was not found under `ab_royer/`; locate it before counting on it.
- **SVI:**
  - `lfe_catalog_repair/results/e5r_families_2012.csv` has 130 families with S offsets only (`offset_s`).
  - `s10_dvv/pilot/tables/offsets_sp.csv` has an `sp` for 3 families × 21 stations, inherited from the template npz `sp` dict (`bw_gate_common.py:118-130`). Provenance (measured or modelled) is unverified.
  - SVI P must be measured on stacked vertical templates, boreholes B009/B010/B011/B926 first.
- **Earthquakes:** `svi_science/results/g2_local/tables/eq_local.parquet`, 41,207 events M ≥ 1 (45.9–50.5N). By depth: 30.7k ≤ 15 km; 6.6k at 15–25; 2.5k at 25–35; 763 at 35–45; 662 > 45 km.
  - No phase picks on disk. Fetch per-event QuakeML (USGS phase-data product; NRCan FDSN `includearrivals`).
- **Software:** NLLoc/Vel2Grid/Grid2Time (`~/.local/bin`), `~/HypoDD_3d_build/src/hypoDD_3d`, pykonal 0.4.1, pyocto, obspy.
  - No tomoDD (not pip-installable; request the source, or write a pykonal + LSQR Vp/Vs inversion).

## 2. Resolution
- **LFE rays alone:** an upward fan. Each family yields one path-averaged κ, with zero vertical resolution.
- **Vertical resolution comes from:**
  - (a) families at 29–46 km;
  - (b) crustal earthquakes (paths above 25 km);
  - (c) about 1,400 intraslab events whose rays cross the LFE zone.
- **Differential Wadati** (LFE κ minus co-located crustal-earthquake κ) gives the κ of the 25–45 km layer per ~30-km along-strike zone. Coarse, but a body-wave number independent of receiver functions.
- **POLARIS 2003–06** (13 stations, `polaris_coverage.csv`) and the WA arrays raise the station count, not the crossing-ray geometry.

## 3. Novelty
- Audet/Bostock receiver-function κ (2.35–2.8) is a thickness–κ trade-off averaged along the profile.
- Merrill et al. 2020 (WA forearc DD tomography) is the closest prior art. Check whether it used LFEs as sources.
- An along-strike κ profile of the LFE source layer from the LFE families themselves, compared with family firing rate and ETS segmentation, is new at that resolution.

## 4. Pilot (WA NW, < 1 day)
1. **Per family κ:** tS−tP vs tP, picks with rP, rS ≥ 1.2, ≥ 6 stations; slope + 1 = κ; bootstrap over stations. Acceptance: ≥ 50 families with κ σ ≤ 0.03 and residual rms ≤ 0.25 s.
2. **Controls:** scrambled control (permute station labels within each family) — the real κ spread across families must exceed the control's at 95 %. Leave-one-station-out: |Δκ| ≤ 0.05.
3. **Mapping:** κ vs family depth and along-strike position. Test against station-set (site-term) confounding by regressing κ on station subsets.

**Fail signature (infeasible):** rms > 0.3 s or leave-one-station-out shifts > 0.05. Then, before abandoning, try stacking vertical stacks of families within 3 km, aligned on S, to sharpen P.

**If it passes:** fetch earthquake picks and do the differential Wadati; only then hypoDD_3d/tomography.

**What would change my mind:** if WA rP ≥ 1.2 P picks give κ scatter indistinguishable from the scrambled control, the LFE P data cannot carry κ and the project must rely on earthquake tomography alone.
