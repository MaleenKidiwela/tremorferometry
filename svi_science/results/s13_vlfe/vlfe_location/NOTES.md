# vlfe_location — NOTES (work log, PT)

## 2026-10-16
- ~10:00 am PT — folder created; PLAN.md written (~10:02) before any offset number.
- ~10:04 — t1_lfe_centroids.py: final event set = 312 rows (242 ETS-period tier events = the catalogue's "242"; 52 between
  raw; 18 'other_episode'); 188 have ≥ 1 located LFE in ±10 min (the catalogue's 189 "with LFEs" minus one whose only LFE
  is an unlocated family).
  - Note on the scrambled control: a within-episode re-pairing is a permutation, so it keeps the MEAN offset vector exactly
    (by construction). It therefore tests tracking (does each VLFE sit nearer its own LFE centroid than other events'?),
    not the bias. The bias test is the bootstrap over episodes, as pre-registered.
- ~10:06 — t2_zsurface.py per year: the production z surface is reproduced exactly (best node 312/312, z to 1e-3).
  Diagnostic seen in passing: in 2012 the production run has NLLB/PFB/YOUB/LZB off and the WA stations Z-only, so the
  western side of the grid has no live station (tabulated in t4 / RESULTS).
- ~10:08 — t23_summary.py: figure bug (8 surfaces into 8 axes, overwriting the width panel) fixed; relay self-detections
  (17 events in 2018 that built RC0–RC7) flagged in Test 3.
- ~10:12 — Amendment A1 (before any Test-4 number):
  - Stack terms use only pooled S13 stacks (svi + wa, per episode scalef > refit > scale, as scan_calib.pick) with
    thrust z ≥ 4 and ≥ 4 usable channels; per-stack M0 refit with the variant's own GF (lag ±8 s); per station-component
    term = inverse-variance mean of (template·stack)/(M0·template²) over stacks.
  - Term use: SNR ≥ 2 and positive → used, clipped to [0.25, 4]; SNR ≤ −2 (opposite sign) → channel dropped;
    otherwise 1. Channels without stack coverage keep 1.
  - Columns added: (b1) Okada with g = 1; 'wsst' = whole-space STATIC ×2 × step (separates the dynamic-vs-static effect
    from the free-surface effect); all 2010–2018 ETS events as an extra subset.
  - Okada check: unit-moment thrust at 35 km, half-space (same Vp, Vs, ρ as s13_green) vs whole-space static ×2 at 0, 20 km
    (both azimuths) and 40 km N: same signs everywhere, amplitudes 10–15 % lower. So for STATIC fields the flat ×2 is
    already close to the half-space; any large change in Test 4 must come from the dynamic vs static difference.
  - (a) layered GF: no offline code — pyrocko 2026.6.2 present but no qseis/qssp backend binaries and no GF store
    (gf_store_superdirs = []), no CPS/EDGRN/axitra. Not run (would need a download/build: ask the user).
- ~10:13 — t4_relocate.py (2 threads, 92 s): first run failed on a two-part stack name (parse fix); then done. No GF
  variant passes the pre-set bar (best: whole-space + stack terms, −2.2 km median). t4_fig.py → figures/t4_gf_variants.png.
- ~10:15 — Test 5 (post-hoc, disclosed): t5_diag.py — LFE spread inside the window (RMS 11 km) and live western channels
  (≤ 2 in every event, so no with/without-west split possible); live-channel table by year.
- ~10:20 — RESULTS.md written. Catalogue folder untouched (only read). Peak memory < 2 GB; ≤ 2 threads throughout.
