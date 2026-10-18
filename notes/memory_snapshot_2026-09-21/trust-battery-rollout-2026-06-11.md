---
name: trust-battery-rollout-2026-06-11
description: "10-station family-trust-battery rollout result — 341 GOLD certified LFE families, two-gate scoreboard, the inversion backbone"
metadata: 
  node_type: memory
  type: project
  originSessionId: 18c98477-b5c6-4af0-9788-f240c207a136
---

**UPDATE: extended to 29/35 stations, 937 GOLD total (2026-06-11 PM, scripts/rollout_queue2.sh).** Added the remaining PB boreholes (B004/B005/B011/B013/B020/B023/B026/B028/B030/B032/B033/B036/B039/B045/B049/B204/B927/B932/B941). The 6 ANCHOR stations (PGC/GNW/NLLB/HDW/COLT/COR) are DEFERRED per user (not yet run; networks CN.PGC/CN.NLLB/UW.GNW/UW.HDW/UW.COLT/IU.COR; NLLB npz symlinked _100km). Test-station waveforms KEPT on disk: B927 (58G, family-rich), B941 (52G, ambient-coherent), B935 (41G) — for the pending ambient-ACF/shallow-monitor + deep-sensitivity tests; their janitor keep-alarms are benign. Original 10-station detail below.

10-station family-trust-battery rollout COMPLETE (2026-06-11), via scripts/rollout_queue.sh → run_trust_battery.sh → family_trust_tier1.py (v2). Each station: reversed-template (time-reversed = guaranteed-fake) calibration densified on 2 sample years, then real families scored by coda coherence (+2..+4s, outside the 2-s selection window) in σ above an nwin-matched random-time null, plus day/night source-independence.

**Two gates (use BOTH; standardized this run):**
- GOLD = coda_sigma > the station's own fake-MAX (strict tail rule the user emphasized) → the INVERSION set.
- TRUSTED = coda_sigma > fake-95th + daynight_cc>0.3 → broader; the extra families over GOLD are a watch list (sit in the fake tail).

**Result: 341 GOLD / 433 TRUSTED / 114 FAIL across 781 scored families (10 stations).**

Per-station GOLD | TRUSTED | FAIL | real-median-σ | fake-max:
- B022 67|87|4  16.9σ  14.3 (elevated floor)
- B040 62|65|3  17.5σ  9.3
- B001 60|66|7  13.7σ  9.6
- B018 40|44|3  18.0σ  11.0
- B928 33|38|0  8.1σ  6.3 (clean, 0 fail)
- B201 23|29|6  14.2σ  15.4 (elevated floor)
- B014 19|35|8  9.3σ  15.6 (elevated floor)
- B003 16|17|40 0.5σ  6.6 (marginal: real med = fake med)
- B017 15|32|2  3.1σ  5.8 (clean separation but modest)
- B935 6|20|41  1.3σ  6.8 (marginal; positive-control, waveforms KEPT)

Tiers: RICH (B022/B040/B001/B018/B201, real med 13-18σ) = SW-WA/Olympics core, does the heavy lifting; SOLID (B928/B014); MARGINAL (B003/B017/B935, real med ≈ fake med, 40+ fail, only tail usable).

CAVEATS: (1) Elevated-floor stations B201/B022/B014 (fake max 14-16σ) have a coherent persistent site signal — even reversed fakes get high σ + day/night coherence; per-station reversed calibration auto-raised their bar (this is WHY we calibrate per-station, not a fixed σ). Sanity-check their dv/v vs the resonance before inversion; trust their GOLD tail, not the marginal TRUSTED band. (2) Every station over-detects ~45-50/day/family (drowned-list inflation) — GOLD is about whether detections STACK into real coda, independent of counts.

Verdict CSVs: data/family_trust_tier1_<STA>.csv + merged data/family_trust_master.csv. See [[coda-window-2to4-correction]] for the measurement band, [[fault-tomography-next-phase]] for what GOLD feeds. Next: Phase B = GOLD-only dv/v + re-derive ETS null + inversion on certified families; plus B935 detection-level cleaning experiment (waveforms retained, data/.keep_wf_B935).
