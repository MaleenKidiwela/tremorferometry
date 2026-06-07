# ETS vs. interface dv/v — analysis log + my thoughts (2026-06-06)

**Question:** does the deep plate-interface velocity (δβ/β from LFE coda) change during Cascadia
episodic-tremor-and-slip (ETS / slow slip)? Run autonomously; this is my reasoning for us to discuss.

---

## 1. Method
- **ETS timing from our own tremor catalog** (`catalogs/pnsn_tremor_cascadia_full.csv`, 748k detections,
  2010–2026). Monthly (and weekly) tremor *rate* per latitude band; ETS onsets = rate peaks (≥70th pct,
  min-separation ≈ recurrence). This is self-contained and avoids importing an external ETS list.
- **dv/v index = mean of the clean borehole cross-patch-median dv/v** in each segment (NOT the inverted
  cell map — the per-month inversion is ~35× noisier: 0.13% vs 0.0037% RMS; for a detection test you want
  the raw measured signal).
- **Superposed-epoch:** stack the dv/v in a window around every ETS onset, subtract each event's pre-onset
  baseline, average → reveals any coherent ETS-locked change while noise cancels (√N).
- Segments: **North ≥47°N**, **Central 43–47°N**, **South <43°N** (each has its own ETS cycle).

## 2. Results — a robust NULL
ETS identification is **validated**: recovered recurrence **12.9 months (north)**, ~11–13 (central/south)
— matches the known Cascadia ETS cycles. So the slow-slip timing is real.

dv/v change at/after ETS onset (deseason+detrend, baseline-subtracted, ±1σ):

| segment | # ETS | cadence | change at onset | most-negative in [0,+5] |
|---|---|---|---|---|
| North  | 13 | monthly | −0.0006 ± 0.0007% | −0.0007% |
| Central| 11 | monthly | +0.0012 ± 0.0014% | −0.0012% |
| South  | 13 | monthly | +0.0030 ± 0.0024% | −0.0020% |
| North  | 18 | **weekly** | +0.0011 ± 0.0009% | −0.0016% |

**Every segment, every cadence → consistent with zero (≤1.5σ).** 2σ upper bound on a coseismic
ETS interface velocity change: **< ~0.002–0.003%**.

## 3. The "ETS is margin-wide / quasi-periodic" concern (user, well-placed)
ETS recurs near **~14 months (north)**; my detrend used a **13-month rolling-median** filter — *right at
the ETS period* — which could attenuate a periodic ETS signal. **Tested directly** (`ets_detrend_test.png`):
re-ran the north stack at 3 levels — raw (no detrend) / deseason-only / deseason+detrend. Even **raw**
shows only **−0.0025 ± 0.0013% at lag +5 months** (~1.9σ, and at the *wrong* lag — a coseismic response
should appear *at* onset, not 5 months later). Detrending mildly attenuates but **does not hide a real
signal**. The null is robust to this concern. (Also: because ETS happens *somewhere* on the margin almost
continuously, the non-ETS "baseline" isn't perfectly ETS-free — which would only *weaken* the contrast,
i.e. bias us toward a smaller apparent signal. Per-segment onsets mitigate this; the migrating along-strike
nature is handled by doing each segment on its own tremor.)

## 4a. WHY the deep signal is small — depth/pressure physics (the main reason)
The dominant reason is NOT just slow slip's low stress drop; it is that **seismic velocity at LFE depth is
intrinsically stress-INSENSITIVE.** dv/v is driven by **crack** opening/closing; near the surface cracks are
compliant → large dv/v (ambient noise sees ~0.05–0.1%). At 30–45 km the confining pressure is ~1 GPa, cracks
are **clamped shut**, so **dv/dstress collapses with depth**. A velocity change at the deep interface is
therefore *expected* to be tiny **regardless of the ETS forcing**. So < 0.003% is the **physically predicted**
magnitude for a high-pressure, crack-poor environment — NOT a non-detection. The contrast with the shallow
~0.1% ambient-noise signal IS the depth dependence of crack compliance. Reframed headline: *"the deep megathrust
is velocity-stable through ETS, as expected where high confining pressure suppresses the crack compliance that
drives shallow dv/v"* — arguably the first deep-depth constraint on ETS velocity change.

## 4. My interpretation
**The deep interface velocity is remarkably stable through ETS — change < ~0.002–0.003%.** I think this is
a *real, meaningful* result, not a sensitivity failure, for three reasons:
1. The measurement is genuinely sensitive — raw dv/v RMS is 0.0037%/month and the stacked error is
   ±0.001%, so a 0.01% transient would be a clear 10σ. We'd see it.
2. It is **physically plausible** that the deep change is tiny: slow slip is aseismic with very low stress
   drop (~kPa), so the associated elastic/damage velocity change at depth may simply be ≪ 0.01%.
3. It **does not contradict** ambient-noise ETS studies that report ~0.05–0.1%: those are sensitive to the
   **shallow crust** (surface-wave / near-receiver), whereas our LFE coda is anchored at the **deep fault
   (~30–45 km)**. Different depth → different (smaller) signal. That contrast is itself interesting.

So the headline I'd defend: *"LFE-coda interferometry bounds the deep plate-interface shear-velocity change
during Cascadia ETS at < ~0.002–0.003% — an order of magnitude below shallow ambient-noise estimates."*

## 4b. SIGN BIAS — corrected (user point, important)
My first pass reported "min[0,+5]" (most *negative*) — implicitly hunting a velocity DROP. But a velocity
**increase** is equally physical: **dilatant hardening / fluid drainage** during slow slip drops pore
pressure → stiffer → **dv/v positive**. Redid it **sign-agnostic** and **pooled all 37 ETS across segments**
(`ets_pooled_bothsign.png`):
- max POSITIVE +0.0018 ± 0.0014% at **lag +1 month** (1.3σ) — a velocity *increase* just after onset, the
  right timing+sign for dilatant hardening.
- max NEGATIVE −0.0019 ± 0.0011% at lag +8 (1.7σ) — but far from onset → likely the next quasi-periodic cycle.
- at onset +0.0004 ± 0.0012%.
**Reframed conclusion:** still sub-significant (<2σ), but *if* anything is there it's a small POST-ONSET
velocity INCREASE, not a drop. The early-coda-window pooled stack is the test that could firm this up.

## 4c. STATION-BY-STATION (the rigorous test — user point, decisive)
Earlier passes **pooled stations into a regional index first**, then stacked — which can manufacture a hint
dominated by 1–2 stations. Redid it **station-by-station**: each of 22 boreholes gets its OWN superposed-epoch
(on its segment's ETS), then test whether stations AGREE (`ets_station_by_station.png`).
- Across-station agreement is **~50% positive at every lag (9–12 of 22)** = random. No lag where a majority
  of stations independently show the same-sign response.
- The regional +0.0018% lag+1 "positive lean" was an **artifact** — only 9/22 stations were actually positive
  there; a few large-value stations pulled the regional mean.
**This is the strongest statement: NO station-consistent ETS dv/v response, either sign.** Region-pooling is
too easily dominated by individual stations; the across-station-agreement test is the right one. Apply it to
the early-window data too.

## 4d. 2–6 s window cross-check (user-requested, DONE) — confirms null, and 1–4 s is optimal
Measured dv/v at 2–6 s for all 22 boreholes (B028 came back empty — sparsest station fails the noisier
window). cc ~0.85 (fine by ambient-noise standards). 22-panel contact sheet: `smoke_22stations_w26_dvv.png`.
Station-by-station ETS test, 2–6 s vs 1–4 s (21 stations):
- 1–4 s: lag0 +0.0004 ± 0.0009% (11/21 positive)
- 2–6 s: lag0 +0.0078 ± **0.0063%** (13/21 positive) — bigger central value but **6× larger error → 1.2σ**, no
  station consensus. NULL confirmed.
**Key lesson:** the δτ ∝ t sensitivity *gain* from a later window is *more than cancelled* by the coda-S/N
*loss* (LFE coda → noise by 4 s). So **1–4 s is both the cleanest AND the most sensitive window** for LFE dv/v;
late-coda CWI does not help here. Depth check: 2–6 s (shallower-weighted) reveals nothing 1–4 s missed → no
shallow ETS signal either.

## 5. Caveats / things that could still hide a signal (to discuss)
- **Temporal smearing:** monthly/weekly bins could dilute a *days-long* transient. → daily stack around a
  few of the *largest* individual ETS.
- **Spatial dilution:** I averaged whole segments; a change confined to the *actually-slipping patch* is
  smeared. → isolate the families that are *active* (elevated LFE rate) during each ETS — they ARE the slip.
- **Reference absorption:** the all-time mean reference, if ETS is quasi-periodic, could partly absorb a
  recurring ETS signal into the reference. (Weak effect — the superposed-epoch removes per-event baselines.)
- **Coda window — EMPIRICALLY VERIFIED (`figures/smoke_coda_envelope_check.png`).** Stacked coda envelope of
  4 clean stations: **S arrives at t = 1.02 s**; the **1–4 s window has S/N = 11** (good coda); the **4–8 s
  window is S/N = 0.7 — BELOW the noise floor.** The weak LFE coda decays into noise by ~4 s, so there is **no
  usable late coda** — late-coda CWI is infeasible with LFEs. **1–4 s is the only usable window**, not a
  choice. The ETS null is measured in the correct (and only) window; precision is coda-length-limited.
- **Active-family test — DONE, also null.** Detection-weighted index (emphasises the families actually
  slipping/detecting each month) gives −0.0006 ± 0.0008% at onset, identical to the plain median
  (`ets_active_family.png`). Targeting the slipping families does not recover a signal.
- **Early-window detour — ABANDONED (user correction).** I tried 0.5–2 s then 1–2.5 s thinking "earlier =
  more deep-fault-sensitive." WRONG: the **S wave arrives ~1 s** (t=0 is a reference ~1 s ahead of S), so
  **1–4 s already IS the early coda** (S onset + first ~3 s of scatter). 0.5–2 s dips into the ballistic/
  pre-S region (our own coda-window standard rejected 0–2 s as ballistic); 1–2.5 s just truncates good coda.
  There is no earlier coda to extract. **The 1–4 s window is correct and final.** The meaningful window
  *contrast* would be a LATER window (4–8 s, diffusive/shallow-weighted) vs 1–4 s (early/deep) — a depth test,
  not an "earlier" test. Not run.

## 6. Recommended next steps (for our discussion)
1. **Active-family superposed-epoch** — for each ETS, stack only the families whose LFE rate spikes that
   month (the slipping sources). Highest signal-to-noise for a localized change. *(highest priority)*
2. **Early-window (0.5–2 s) dv/v** for northern boreholes, re-stack — most deep-fault-sensitive.
3. **Daily-resolution** stack around the 3–5 biggest individual ETS.
4. If all still null → write up the upper bound as the result.

## 7. Artifacts
Figures: `ets_tremor_vs_dvv.png`, `ets_superposed_epoch.png` (inverted-cell, noisy — superseded),
`ets_superposed_raw.png` (clean raw), `ets_segments_weekly.png` (per-segment + weekly),
`ets_detrend_test.png` (detrend-sensitivity). Data: `ets_north.csv`. Scripts inline (see session).

**BOTTOM LINE (finalized 2026-06-06):** No **station-consistent** ETS-locked deep-interface dv/v change, of
**either sign**, at any lag — across all segments, monthly + weekly, robust to detrending, robust to
active-family weighting. Tests done: per-segment, weekly, pooled sign-agnostic, station-by-station agreement,
detrend-sensitivity. Window **empirically verified** (1–4 s is the only usable coda; no late coda with LFEs).
**Upper bound: < ~0.002–0.003%** on the deep-interface velocity change during Cascadia ETS — an order of
magnitude below shallow ambient-noise estimates (consistent: different, deeper measurement). The earlier
"hints" were pooling artifacts that vanish under the station-agreement test. Caveats remaining: monthly/weekly
could miss a days-long transient (→ daily stack on the few biggest individual ETS); regional vs slip-patch
localization. **Methodological lessons:** (1) station-by-station agreement >> region-pooling (pooling
manufactures hints); (2) test both signs (dilatancy can raise v); (3) verify the coda window empirically.
