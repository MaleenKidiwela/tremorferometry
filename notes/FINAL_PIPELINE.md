# ★★ FINALIZED PIPELINE — Cascadia borehole LFE-coda dv/v (settled 2026-07-10)

**This is the authoritative, settled pipeline. No reverse densify anywhere — one forward densify per
station; everything is computed from the forward stacks. The mirror window (the part of the stack BEFORE
the LFE arrival) does double duty: it certifies reliable families AND cleans the noise.**

Supersedes all earlier reverse-densify / trust-battery discussion. If a future session sees notes about
reverse densify, sampled reverse, fwd-vs-rev certification, or the trust battery — those are RETIRED.

---

## PER STATION

| # | Step | What | Tool | Resource |
|---|------|------|------|----------|
| 1 | **Download** | 3-comp waveforms 2007–2026 | `download_borehole_3comp.py` | net |
| 2 | **Detect candidates** | detect near station + picker score; keep an **adaptive ~15k** (P for 15k cands, floor 0.2, cap 0.7 — so weak stations aren't starved) | `discover_nllb_pnsn_driven.py --candidates-only`, `score_candidates.py`, `adaptive_cand_threshold.py` | CPU |
| 3 | **Cluster** | group candidates into families | `discover_gpu.py --fs 100 --cc-threshold 0.80 --min-family-members 3 --min-years 1` | GPU ~1 min |
| 4 | **Select** | coverage-balanced ~300 families | `select_families_coverage.py <disc> <picker_csv> 300` | CPU |
| 5 | **FORWARD densify** | matched-filter search of the whole record for every occurrence. **ONE densify. NO reverse.** | `densify_gnw_gpu.py --fs 40 --top-n 0 --despike-mad 8 --threshold 0.8` (templates = 100→40 resampled) | GPU ~1–2 h |
| 6 | **Daily stacks** | per family, 30-cal-day rolling stack aligned on the arrival; each day has a **coda window (after)** and a **mirror window (before)** | `build_long_window_3comp.py --fs 40 --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 --despike-mad 8` | CPU |
| 7 | **Certify = CAUSALITY** | on the grand stack, `RMS(coda 2–4 s) / RMS(mirror −2..0 s) > 1.5` → real causal LFE. Keep those (~1/3 of families). | `finalize_causality.py` | CPU |
| 8 | **dv/v** | coda stretch per day (2–4 s window, origin-anchored at t_s=1.0), median over reliable families (cc_max≥0.6, ≥3 families/day, 15-day rolling median) | `dvv_roll30cal.py --window 2 4 --origin-anchor` | CPU |
| 8b | **Battery gate (after dv/v)** — user-mandated 2nd gate (reinstated 2026-07-11) | Tier-1 stack-vs-random battery on the **causality-certified subset** (NOT the full cap-off set, where it false-negatives). Keep families that are causality-reliable AND not battery-FAIL; recompute Z dv/v over the doubly-gated set. Forward-only (no `--rev-mf` → script's f50=1.0/f95=6.0 fallback). **NEEDS TRACES → run before trace deletion.** | `battery_gate.sh` → `family_trust_tier1.py` + `battery_gate_dvv.py` | CPU |
| 9 | **Noise-clean = MIRROR** | measure the dv/v of the **mirror window** (front of the stacks = noise) and subtract: `residual = fwd_coda − β·mirror` (β = regression, on deseasoned series). Free, from the same stacks. **No reverse.** | `build_mirror_npz.py` → `dvv_roll30cal.py` on the mirror npz → subtract | CPU |
| 10 | **Finalize + map** | corrected dv/v per station; refresh the map | `build_borehole_dvv_map.py` | CPU |

## AFTER ALL STATIONS

| 11 | **Inversion** | combine every station's corrected dv/v into the 4-D δβ/β map on the plate interface, with an explicit **shallow/station nuisance term** to separate shallow noise from the deep megathrust signal | `fault_tomography/inversion/` | — |

---

## WHY THIS IS THE PIPELINE (the settled logic)

- **The mirror window is the key trick.** Each daily stack is aligned on the LFE arrival. AFTER = coda
  (the earthquake). BEFORE ("mirror") = pure background noise. From the SAME forward stacks:
  - **Certify** (step 7): coda ≫ mirror ⇒ real causal LFE (a real quake's energy comes after it arrives;
    ringing/noise is symmetric). Verified equivalent to the old fwd-vs-reverse certification at 99% (corr 1.00).
  - **Clean** (step 9): the mirror's own dv/v = the noise-field's apparent stretch ⇒ subtract it.
- **Reverse densify is NOT needed.** It was the old certification method AND a candidate noise reference,
  but the mirror does both jobs from the forward stacks. Cost: one densify per station, not two.
- **Validated:** the raw dv/v is ~65% noise-field artifact (reversed-control test). The MIRROR-corrected
  residual recovers a real velocity signal locked to ETS episode onsets: **+0.043% co-onset, p=0.019** at
  B926+B011 (reverse-corrected gave +0.039%, p=0.003 — reverse is marginally cleaner but not worth 2× GPU).
  Figures: `figures/mirror_corrected_dvv_ets.png`.

## KEY PARAMETERS (do not drift)
- Densify/stacks/dv/v at **fs=40 Hz**; discovery+picker stay 100 Hz. Cap OFF (`--top-n 0`).
- Coda window **2–4 s**, `--origin-anchor` (t_s=+1.0 s). Mirror window = the pre-arrival slice (~−3..−1 s).
- Causality threshold **1.5**. Adaptive candidate target **~15k** (**floor 0.2**, cap 0.7 — lowered from 0.5 on
  2026-07-13; 0.5 starved low-cal stations to 0 families). ~**300** families/station (hard cap = densify BUDGET).
- Flag any station finalizing with **<20 certified families** (starved or genuinely weak) — don't silently pass.

## ★ FLEET INCLUSION RULE (frozen 2026-07-13, BEFORE any inversion — Merlin) — rule-before-result
A station enters the **inversion + any fleet-pooled statistic** ONLY if: **z_certified ≥ 20** AND **causality
survival ≥ ~15%** of densified families (ratio p50 ≥ 1.0). Flagged stations stay on disk + on the map
labeled **FLAGGED**, but are excluded from analysis and never presented as a result.
- **★ KEY LESSON (B045):** lowering the candidate floor recovers **candidates, not certified families**. June
  trust-battery **GOLD count does NOT predict causality-certified count**. B045 (GOLD 52) recovered to 15k
  candidates → 172 families → 19 LFE-label → **4 causality-certified** (median causality 0.79 = noise-dominated).
  Mendocino = UNRESOLVED. The strong stations are the FIRST ones done; the tail is inherently weak. **Do not
  force a station count** by padding with weak/low-cal stations — finalize at the honest count.
- **Early-exit** (`borehole_gpu_phase_fwd.sh`): <60 LFE-label families → skip densify (can't reach ~20 cert).
- **HONEST LEDGER (2026-07-13):** 25 stations ≥20 certified (min above line B022=24); 3 FLAGGED
  (B045=4, B009=10, B027=13 — clean gap 13→24); B030 = last pending slot. B932 killed (P_for_15k=0.070,
  stale cand marked). Junk refused (B935/B204/B028: duplicates/no new interface coverage).

## STATUS (2026-07-11, after server crash + restart)
**Finalized (causality Z dv/v):** B926(110), B011(119), B001(98), B004(64), B927(55), B928(106). ⚠ these 6
were finalized BEFORE the battery gate was reinstated → their battery gate is NOT yet applied; 5 of them
(B926/B011/B001/B927/B928) also have their **traces deleted** (rolling buffer), so backfilling the battery
gate on those needs re-download. B004 traces still on disk. **DECISION (user 2026-07-11): come back to these
6 AFTER the current wave finishes** — B004 backfill is free (traces present); the other 5 need re-download.

**Orchestration (rolling buffer, downloads HELD):** `scripts/rolling_wave.sh` — serial; per station
prep→gpu_phase→cpu_tail→**battery_gate**→**DELETE TRACES**. Queue = complete on-disk stations
`B009 B003 B010 B012 B024 B027 B035`. OOM guard = `scripts/oom_guard.sh` (persistent; SIGSTOP/CONT densify
if anon>150). Old drivers `process_wave.sh`/`borehole_cpu_tail.sh` still exist but the tail no longer
overlaps (rolling buffer is serial so traces free immediately).

**In flight at restart:** B003 stacks, B009 stacks (both finishing their pre-crash tail), B010 prep
(scoring). rolling_wave waits for each via `wait_station_idle`.

**GAP AUDIT (user 2026-07-11: "don't skip gaps"):** `scripts/gap_audit.py` flags INTERIOR years with <200
dv/v days (a truncated year fakes a dv/v discontinuity) and reports whether traces are on disk (download
truncation = re-downloadable) vs deleted. Wired into the map daemon → runs on every finalize, logs FLAG
lines to rollout_status. **Current backlog:**
- B011 2017 (59 days, download truncation) — FIXING NOW (`fix_2017_b011.sh`: re-download 2017 → re-densify
  reliable → merge into npz → recompute dv/v → delete 2017 traces → rebuild map).
- B926 2017 (59 days, download truncation) — deferred with the 6 (same fix, needs re-download).
- B012 2025 (132 days, traces deleted) — likely download not caught up to 2025; re-download 2017-style to confirm.
- B004 2020/2023/2024 — NOT a download gap: 2023/2024 have FULL traces (365/364) but sparse dv/v =
  genuine low LFE activity (B004 is the tremor-poor P≥0.55 re-run). Do NOT "fix" by re-download.

**Map auto-refresh:** `scripts/map_refresh_daemon.sh` (nohup, decoupled) watches `data/*_3comp_summary.json`
and rebuilds `figures/borehole_dvv_map.html` (via `build_borehole_dvv_map.py`) within ~2 min of any station
finalizing. RESTART IT after a crash (alongside oom_guard + rolling_wave). Battery step removed from the
rolling tail (`battery_gate.sh`/`battery_gate_dvv.py` still exist for the deferred backfill, just not called).

**CONTINUE-THROUGH-29 (user 2026-07-12): downloads UN-held; driving to 29 completed boreholes.**
Driver = `scripts/continue_29.sh` (nohup): throttled download feeder (≤4 buffered) + serial rolling
processor (prep→gpu_phase→cpu_tail→delete traces), battery-free, stops at 29. RESUMABLE after a crash —
just re-launch it (uses `data/.dl_done_<STA>` markers + guarded steps + skip-if-`<sta>_3comp_summary.json`).
Queue = B031 B943 then prepped boreholes ranked by GOLD (B036 B033 B039 B022 B040 B045 B032 B013 B018 B005
B932 B201 B030 B023 …). 13 done at launch → target 29. Keep oom_guard + map_refresh_daemon alive alongside.

**SCOPE (user 2026-07-11): BOREHOLES ONLY — skip broadbands/land for now.** Wave is 100% PB boreholes;
CN.PGC (46 G, only broadband on disk) left untouched, not reprocessed. Broadband/land stations revisited
only after all boreholes done. Battery gate also skipped (causality-only product).

**Downloads HELD (user 2026-07-11):** "downloads way ahead; finish+delete traces before downloading new."
Complete on disk: B003 B004 B009 B010 B012 B024 B027 B035. Partial (deferred): **B031 (~85%), B943 (~57%)** —
resume their download only after the backlog clears and traces free. Then B012/B024/B027/B035 already local.

**Battery gate (step 8b) is being VALIDATED on B009/B003 first** — if it guts the certified set like the old
full-set run did (notes/2026-07-10), revisit the fake-calibration (f95=6.0 fallback) before fleet trust.

Inversion is the final step after all stations.
