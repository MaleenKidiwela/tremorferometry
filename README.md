# tremorferometry

**LFE coda-wave interferometry for 4-D shear-velocity change (δβ/β) on the Cascadia plate interface.**

Low-Frequency Earthquake (LFE) families are used as *repeating sources* for coda-wave interferometry.
Because LFEs sit on the megathrust / transition zone, the dv/v sensitivity kernel is biased toward the
slipping patch — the opposite of ambient-noise dv/v, which is dominated by the shallow crust. The novel
axis is **depth**.

---

## Network

| tier | stations | certified families |
|---|---|---|
| PB boreholes | 30 | 2,095 |
| broadband fleet | 124 INCLUDE / 55 FLAG (178 processed) | 10,608 (with anchors) |
| anchors (PGC, SHB, CLRS) | 3 | — |
| **total** | **187** | **12,703** |

Span 39.7–50.4°N, 2009–2026, everything processed at 40 Hz. A family is *causality-certified* when
`RMS(coda 2–4 s) / RMS(pre-arrival mirror −2..0 s) > 1.5`. Frozen station-inclusion gate:
**≥20 certified families AND ≥15% survival** (certified ÷ densified).

Validated two independent ways: borehole↔co-located broadband (B011↔PGC, per-family p=0.0005) and the
no-borehole discovery path (CLRS, +120 s shift-null) — the latter is what licenses the fleet.

---

## The finalized pipeline

Per station, driven by `scripts/fleet_station.sh`:

| # | stage | script | output |
|---|---|---|---|
| 1 | download | `scripts/download_broadband.py` / `download_station.py` | day files (resumable) |
| 2 | band pick | `scripts/pick_band.py` | one vertical high-gain band (2nd SEED char `H`) |
| 3 | detect | `scripts/discover_nllb_pnsn_driven.py` | PNSN-tremor-driven candidates |
| 4 | score | `lfe_features/score_candidates.py --target-fs 40` | P(LFE) per candidate |
| 5 | threshold | `scripts/adaptive_cand_threshold.py` | rank-based top-30k |
| 6 | cluster | `scripts/discover_gpu.py` (GPU, cc≥0.80, `--min-years 1`) | matched-filter families |
| 7 | select | `scripts/select_families_coverage.py 'snr'` | ≤300 densify budget, by SNR |
| 8 | densify | `scripts/densify_gnw_gpu.py` (forward only, 2–8 Hz) | detections |
| 9 | daily stacks | `scripts/build_long_window_3comp.py` (≥20 det/day) | `long_window_daily_<STA>_Z.npz` |
| 10 | **dv/v** | `scripts/dvv_roll30cal.py` | `daily_dvv_<STA>_Z_2to4*.csv` |
| 11 | certify | `scripts/finalize_causality.py` | `<stem>_causality_cert.csv` |

**dv/v measurement (stage 10).** 30-calendar-day *trailing* rolling stack → peak-normalize → SVD-Wiener
filter → stretch the **2–4 s** coda against a per-family **all-time reference** (`ref = Rf.mean(0)`),
origin-anchored. The 2–4 s window is load-bearing: 1–4 s dilutes the stretch ~50× toward zero because the
pinned direct-S dominates the window energy. Never take a median dv/v.

`--cert-csv` restricts the measurement to causality-certified families. Family measurement is fully
independent (own stacks, own SVD, own reference), so this is bit-identical for the kept families and
skips ~77% of wasted work.

### The joint 4-D inversion

```
assemble_res_catalog.py   →  pairs.csv, cells.csv, pair_months.parquet   (the finalized dv/v tensor)
build_G_captured.py       →  G.npz          capture-weighted single-scatter operator
build_era_table.py        →  era_table.csv  instrument/response eras
invert_dvv_4d.py          →  inversion_4d.npz
plot_4d_interface_maps.py →  DEEP    (interface δβ/β, from MODEL)
plot_4d_surface_maps.py   →  SHALLOW (near-receiver field, from SITES)
```

One system solves the interface model **and** per-station site terms jointly:

    m = argmin ‖W(d − [Gc S]m)‖² + λ_f²‖Lm‖² + λ_s²‖m_site‖²

`inversion_4d.npz` carries both halves — `MODEL` (n_windows × n_cells, the deep interface) and `SITES`
(n_windows × n_stations, the shallow near-receiver field) — plus the gate statistics (`idx`, `null_idx`,
`idx_pctile`, `VR`, `closure_*`, `mda_model`, `bound_vals`).

Grids: `res_catalog_g20` (0.2°, 293 cells) and `res_catalog_g40` (0.4°, 101 cells), both 189 stations,
2009-01 → 2026-07. Cells are represented by their **LFE family centroid**, not the geometric grid centre.

---

## Layout

```
scripts/                        the 11-stage pipeline + drivers + infra watchdogs
fault_tomography/inversion/     the 6 scripts of the joint 4-D inversion
  res_catalog_g20/  g40/        assembled tensors, operators, results, figures
lfe_features/                   score_candidates.py (the LFE picker)
src/tremorferometry/            the library (stretch_dvv, matched filter, io, qc)
notes/                          methodology, resume notes, the adversarial audit ledger
archive/                        223 superseded / exploratory / one-off scripts (nothing deleted)
```

`archive/MOVES.json` records every move. Start with `notes/METHODOLOGY_END_OF_RUN_2026-07-21.md`
(methods + coverage) and `notes/POST_DVV_ANALYSIS_2026-07-21.md` (the audit trail).

---

## What is and isn't in git

Versioned: all code, the **assembled inversion tensors and results** (`res_catalog*/`), the 191
causality-cert files, and the small shared inputs (slab geometry, station rosters). ~40 MB — enough to
reproduce **both inversions from a clone**.

Not versioned (regenerable, far past GitHub's limits):

| product | size |
|---|---|
| daily stacks `long_window_daily_*.npz` | 283 GB |
| per-family daily dv/v CSVs | 7.8 GB (17 files >100 MB) |
| `figures/borehole_dvv_map.html` | 176 MB |
| `catalogs/pnsn_tremor_cascadia_full.csv` | 55 MB — **required by `invert_dvv_4d.py`**, re-fetch from the PNSN API |

---

## Results

**Resolution (signed off).** Network *geometry* resolves the deep interface down to the 44 km cell scale
noise-free — there is no geometric obstruction; the limit is SNR. At measurement noise the finest resolved
scale is ~70 km (recent), sharpened from ~180 km as the fleet grew. Deep-index precision improved
0.64% → 0.12% and above-noise deep modes ~3 → ~24. The 3-D **volume** checkerboard reproduces ~0.85 at
250 km; only the older 2-D *interface* 0.77 is retired (a thin interface captures a median ~5% of the
volumetric coda sensitivity, so unit-sum kernels inflated it ~12×).

**Inversion (preliminary — not signed off).** The deep megathrust is **velocity-stable within resolution**
— a null-gated, spatially-resolved *bound*, not a detection. A planted ~1% coherent deep patch is
recovered at +0.89% (correct sign), so this is a genuine bound rather than a coverage failure. The ETS
composite is sign-correct and leave-one-station-out robust but sub-threshold (p≈0.055 against a
pre-registered p<0.01).

Both results come from an adversarial review loop (fix → re-run → re-audit) that has so far caught a sign
flip, a rigged null, a λ transplant, a baseline-drift bug, and an instrument-era contamination path.

### Known issues

- **Search-range clipping.** dv/v was measured with `eps_max=0.02`, so any true |dv/v| > 2% returned
  *exactly* 2%. Fleet-wide this affects 1.8% of inversion input, but it is concentrated (NEMA 37%, B036
  31%, B035 29%). Ten stations have been re-measured at a uniform `eps_max=0.05` (`scripts/rerun_unif05.sh`,
  `*_unif05.csv`); **the inversions on disk still use the clipped input.** Whether the freed 2.4–3.0%
  values are real or cycle-skips is unresolved — the discriminator is lapse-proportionality.
- **Mirror v2.** A fleet-wide mirror-corrected pass would cut the noise floor ~65%. Since geometry already
  resolves 44 km, that buys resolution directly. It is the single biggest lever and is not yet run.
