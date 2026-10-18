---
name: nisqually-drop-artifact-check
description: "The old GNW Nisqually-2001 dv/v \"precursor\" drop does NOT survive standard reprocessing — likely a processing artifact; investigation DEFERRED"
metadata: 
  node_type: memory
  type: project
  originSessionId: 400c26f9-17ff-46fc-914c-9c6a7671906d
---

**DEFERRED (2026-06-06, user: "forget about era 1, come back once everything else is done"). The project's founding observation — a dv/v DROP at GNW around the 2001 Nisqually M6.8 — appears to be a PROCESSING ARTIFACT of the old non-standard pipeline. Must revisit + finish the response-removal test before any precursor claim.**

**Context:** GNW was rebuilt 2026-06-06 with the standard pipeline (no response removal / `--no-deconv`, coverage+top-10% selection floor 6 = 81 families, despike 8, per-era ref bounds 2010-12-01 / 2019-05-07). The OLD GNW (deleted) used response removal + deseasonalization + a merged hand-curated ~96-family set, and its Nisqually-drop figures were the response-removed/deseason per-family ones (`c963_nisqually_drop`, `c334_resp_deseason`, `nisqually_dropmap/dropdist`).

**Evidence so far that the drop is an ARTIFACT** (new 81-family standard pipeline, era 1):
1. Spring-2001 change only **−0.009%** (within ±0.048% scatter) — essentially flat.
2. **Deseasonalization RULED OUT** — raw and deseasonalized both give −0.009% (`figures/smoke_GNW_nisqually_check.png`).
3. **No near-source signature** — families <40 km from the Nisqually epicentre (47.15,−122.73) change −0.011%, families >60 km change −0.015% (far slightly LARGER). A real Nisqually-induced change would be strongest near the source; it isn't (`figures/smoke_GNW_nisqually_distsplit.png`).

**Remaining suspect = RESPONSE REMOVAL** (only GNW had it; family count 96→81 is too small a change to matter). **THE definitive test (not yet run):** rebuild the NEW 81 families' stacks WITH deconvolution (raw waveforms KEPT at `data/waveforms/UW.GNW` 58G; re-fetch the GNW response XML — it was deleted; stack via `build_long_window_resp.py` WITHOUT `--no-deconv` + `--inv`), re-measure era-1 dv/v, and check if the Nisqually drop REAPPEARS. Same families, only deconvolution flips → isolates the cause. The old 96 cannot be faithfully recreated (deleted + non-standard merge of gnw/gnwcircle/gnwgap/cg runs).

**GNW handling now:** canonical dv/v = full 1995-2026 (3 eras, kept). Tomography uses **era 2+3 only** (`data/daily_dvv_GNW_tomo_e23.csv`, 2010-12+). See [[fault-tomography-next-phase]], [[nisqually-precursor-project]].
