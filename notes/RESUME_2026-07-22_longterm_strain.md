# RESUME 2026-07-22 — Long-term interface STRAIN analysis (OPEN, user thinking on the design)

## POLICY (user, 2026-07-25): NO DEPTH CUT — use ALL illuminated interface cells
Drop the `depth>30` "deep" mask everywhere. The grid represents the WHOLE subduction interface; include
shallow cells (e.g. N. California LFEs ~20–25 km) even though not deep. "Deep index/map" → "INTERFACE index/map".
Family-centroid depth is still used for the KERNEL source depth (unchanged) — the change is only that we no
longer SELECT cells by depth for the index/maps/inversion summaries.


**Where we stopped:** user wants to fully think through the "slow-field" idea before we build it. Resume by
asking if they've settled on the design below, then decide whether to run it.

## The target (reframed this session)
NOT ETS-locked episodic change. The goal is **long-term (multi-year) strain on the deep (>30 km) interface** =
slow evolution of interseismic locking. The ETS composite work (ledger R6) was the wrong hypothesis for this.

## What's recoverable (Merlin, vetted)
- Signal is NOT destroyed upstream: `dvv_roll30cal.py:53` uses a fixed all-time reference → per-pair slow
  variation survives on disk.
- My feared "per-era demean breaks trends" **does not exist** in `invert_coarse_cm.py` (era only drives the
  ±35 d straddle drop + a no-op grouping key). The 2019–24 baseline demean removes a constant → can't change a trend.
- **The real signal-killer = common-mode PRE-subtraction** (`invert_coarse_cm.py:55`): subtracting each
  station's family-mean before the inversion spatially high-passes the slow field at ~50–100 km. FIX = do NOT
  pre-subtract; estimate per-(station-era,month) **site terms JOINTLY** with the deep model.
- **Identifiability:** only the spatially-VARYING (differential) slow field is recoverable = differential
  locking evolution, modulo one margin-wide time function (one connected station–cell graph component). A
  margin-UNIFORM slow change is unidentifiable mirror-free (degenerate with climate/site drift). Bonus lever:
  deep-vs-shallow trend contrast within a station's families (depth-differential).

## KEY design refinement from the user (the reason we paused)
A single 15-yr **slope assumes MONOTONIC strain**. If strain builds → relaxes → rebuilds, the net slope ≈ 0 and
we'd miss it. So DO NOT fit a straight line. Fit a **smooth low-frequency CURVE per cell** (spline ~2–3 yr knots,
or low-pass) that can wander up/down/up while rejecting the fast (30-day-stack) noise. Map the slow FIELD through
time = a **de-noised, slow-only version of the monthly movie**. The slope is just the 1-param special case.
- Smoothness rule: keep variation slower than the noise decorrelation time (lag-1 ρ≈0.64 monthly), smooth away faster.

## The honest referee (same whether signal is a trend or a wiggle)
**Split the network into two independent station-halves → two slow-field maps → do they AGREE on the same
wandering at the same place & time?** Null = permute pair→cell labels WITHIN station (≥500). Year-scramble null
is INVALID for trends (whitens persistence → false positives) — do not reuse it here. Block-bootstrap only as a
sanity cross-check.

## Biggest confound to gate FIRST (cheap, decisive)
**Family-specific detection-rate drift** — LFE catalog grows over 15 yr; the event population feeding each
family's stack drifts; it's family-specific so it SURVIVES common-mode removal and masquerades as differential
deep signal. Gate: do per-pair dv/v slopes track their `n_stack`/`cc_max` slopes? If yes → it's bookkeeping drift,
stop. Physical gate: real velocity change is LAPSE-PROPORTIONAL across the coda (use the dt-vs-lapse machinery).
Direct control: the **mirror** (does the pre-arrival noise window trend too?) — this is mirror-as-DIAGNOSTIC,
distinct from mirror CLEANING/subtraction the supervisor paused, so it's available.

## Honest prior (Merlin)
~15–25% for a spatially-resolved secular DETECTION; **~70%+ for a defensible spatially-resolved BOUND**
(differential deep secular change < ~0.1% / 15 yr at 100–200 km). The BOUND is itself a novel publishable
constraint on locking heterogeneity. 15-yr lever = this dataset's unique asset. Do it, expecting the bound.
Merlin: do it mirror-corrected (confound gate); supervisor paused mirror cleaning → use mirror-as-control only,
or resolve the mirror-method decision.

## Build plan (Merlin's 7 steps, adapted to the smooth-curve refinement)
1. Per-pair SLOW-CURVE table (smooth spline/low-pass, NOT single slope; + era-step offsets; NO common-mode
   pre-subtraction, NO baseline demean; AR(1)/Newey-West errors on the curve coefficients).
2. **Detection-rate-drift gate** (n_stack/cc_max) — run BEFORE any map; kills or clears the idea in one pass.
3. Joint slow-field inversion: cell slow-fields + station site-term nuisances (min-norm gauge, Laplacian, oracle λ — do NOT transplant λ).
4. Significance = split-network coherence + within-station label permutation (≥500).
5. Replicate on long-record subnetwork (310 pairs ≥12 yr, 48 stns, 27 multi-station deep cells) — guards network-growth bias.
6. If null → publish the per-cell 95% BOUND on differential deep secular change (state gauge + mirror-free caveats).
7. Mirror-corrected v2 rerun (boreholes) = the direct control on confound #1.

Full Merlin transcript logic is in this session; ledger rows R5 (family-centroid) + R6 (ETS composite sub-threshold)
are in `notes/POST_DVV_ANALYSIS_2026-07-21.md`. Scripts to reuse: `invert_coarse_cm.py` (data prep),
`ets_composite_confirm.py` (per-pair machinery + null pattern), `pair_months.parquet` (the data).
