# Family verdicts — "which LFE families are good" (snapshot 2026-06-12)

Git-tracked backup of the family quality/trust outputs (the canonical data lives in the
gitignored `data/`; this is a small, versioned second copy on GitHub so the verdicts survive a
home-volume loss). 29 of 35 densified stations certified; 6 anchor stations (PGC/GNW/NLLB/HDW/COLT/COR)
deferred. ~937 GOLD families total.

## Files
- `family_trust_tier1_<STA>.csv` (29) — per-station Tier-1 trust verdicts (the core output).
- `family_trust_master.csv` — all stations merged.
- `family_quality_flags.csv` — older per-family quality flags (gappy/glitch/sensor-era).
- `family_fingerprint.csv` — physical-identity labels from detection timing (CULTURAL/BLAST/NATURAL-like/…).
- `exp_t2a_<STA>_fixed.csv` (B018/B927/B941) — T2a shallow-share **site filter** (block-bootstrap CIs).
- `daily_dvv_<STA>_2to4_CERTIFIED.csv` (B003/B018) — certified per-family 2–4 s dv/v (real families only).

## Tier-1 verdict columns (`family_trust_tier1_*.csv`)
`kind` F=real family / R=reversed-template FAKE (calibration). `coda_sigma` = coda(+2..+4 s) coherence
in σ above a random-time null (primary discriminator). `daynight_cc` = day vs night coda correlation
(source-independence). `verdict`:
- **GOLD** = `coda_sigma` > the station's reversed-fake **MAX** (strict tail; the inversion set).
- **TRUSTED** = > fake-95th + daynight_cc>0.3 (broader; watch list).
- **FAIL** = < fake-median (behaves like a reversed-template fake = noise).

## Critical caveats (do not skip)
1. **GOLD/TRUSTED = real coherent LFE, NOT "deep".** T1 certifies a genuine repeating source vs noise;
   it does not certify depth. Even high-σ GOLD families can be shallow site-carriers.
2. **T2a (`exp_t2a_*_fixed.csv`) is a SITE filter, not a depth certificate.** `cls`=SHALLOW-COUPLED
   (candidate site-carrier), DECOUPLED (deep-or-noise), ANOMALOUS. T1 σ is ORTHOGONAL to shallow-coupling
   → can't σ-weight out site-carriers; the deep inversion input must be T2a-filtered.
3. **Deep velocity CHANGE = an upper bound (~0.02%)**, not a positive detection. Deep SENSITIVITY
   (geometry) is separately confirmable via lapse-dilation (pending `_calS` per-window scale regen).
4. The ambient-ACF shallow monitor is usable at only ~1/4 stations (station-idiosyncratic) → the shallow
   map / T2a depth-filter can only be applied cleanly at good-ACF stations.

## Regenerate
`scripts/run_trust_battery.sh <NET> <STA>` per station (via `rollout_queue*.sh`); T2a via
`scripts/exp_t2a_fixed.py <STA>` (needs that station's `autocorr_dvv_<STA>.csv`). See
`notes/MASTER_PLAN_2026-06-11.md` §7 and `notes/FAMILY_TRUST_TEST.md`.
