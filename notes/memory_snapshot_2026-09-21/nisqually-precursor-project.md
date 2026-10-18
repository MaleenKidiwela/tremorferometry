---
name: nisqually-precursor-project
description: Current work — dual-reference coda dv/v at GNW for the Nisqually 2001 EQ precursor question
metadata: 
  node_type: memory
  type: project
  originSessionId: be83f309-c83d-4336-b8a2-99a689098c27
---

**Authoritative resume doc lives in the repo: `notes/NISQUALLY_GNW_RESUME.md` — read it first.** It has the full, current plan. This file is just the index-level summary.

Active line (2026-05-29): **GNW coda 1–3 s dv/v, NLLB-style** (`figures/smoke_dvv_nllb58_coda_1to3.png` is the template), to look for a **dv/v precursor to the 2001 Nisqually M6.8 EQ** at station GNW (UW network, 1995–2026, 10,662 day-files). 53 high-SNR seed templates (SNR≥15) already built. Map done.

**The whole blocker has been the matched-filter DENSIFY of the 53 seeds.** This session diagnosed and fixed WHY it was slow/crashing — see [[pod-cpu-cap]] and the resume doc. Summary of the locked design the user approved:
- **workers=30** (pod is capped at **32 CPUs** — `os.cpu_count()` lies, reports 176; the old `--workers 48` thrashed the quota → ~13 h instead of ~1 h). Thread env pinned to 1.
- **resample_poly** not obspy FFT-resample (25× faster; GNW needs resampling at 50/100 Hz, NLLB was mostly 40 Hz native so it was fast).
- **reuse day FFT** across templates (scan_day_multi recomputed it 53×/day).
- **forward-chronological, YEARLY chunks** `mf_gnw_YYYY.csv` + resume (skip done years).
- **top-N=100 cap per (family,day)** AFTER per-day QC; keep the cc>1 fix + QC (drop median-count>2000 or median-cc>0.96 days). cc threshold stays **0.8** (not 0.85+: stacker min-det=20 is per-family; cc≥0.85 drops ~half the family-days, worse pre-2001).
- NO streaming, NO bounded-submission (the cap removes the firehose; ~10–30M rows total, few GB).

**✅ DONE 2026-05-29 ~08:27** (ran autonomously while user slept): `figures/smoke_dvv_GNW_coda_1to3.png` + `data/daily_dvv_GNW_coda_1to3.csv` (371,330 measurements, 53 patches, mean cc 0.986, 1995–2026). Densify all 32 years on the GPU in ~17 min; stack → 382,646 daily stacks; dv/v via new parallel `scripts/dvv_coda_parallel.py` (~2.5 min). **NISQUALLY VERDICT (after peeling away artifacts):** a real, **REGION-WIDE ~0.04–0.05% coseismic-shaped velocity drop** at 2001 (min at the EQ, recovery; survives deseasonalizing — GNW seasonal cycle only ~0.006%), in clean resample_poly+per-era data. Initially looked NW-clustered/"wrong place," but that was UNDER-SAMPLING the south (only 6 near-epicenter families at SNR≥15); re-running the 22 SNR≥12 families within 40 km of the epicenter (`run_gnw_circle.sh`) showed they ALSO dip → the drop is region-wide, **consistent with a coseismic origin**. NOT yet established: significance vs the noisy early-era background (needs the null/false-positive test). Lesson: GNW LFE families are 5:1 N:S; use SNR≥12 to sample the Nisqually source region. Getting there required removing 4 artifacts that each faked a bigger signal: instrument steps (per-era ref), **obspy-FFT-resample drift (use scipy `resample_poly`, NOT obspy — the original `build_long_window_daily_all51.py` is contaminated; `build_long_window_resp.py` is clean)**, all-time-ref inflation (use per-era), and season. Also: response removal alone doesn't flatten cross-era steps (per-era referencing does). Full chain + methods + figures in the resume doc. Remaining test if revisited: null/false-positive rate + dual-ref with error bars on resample_poly stacks.

**User directive 2026-05-29: run autonomously — "don't stop until you finish densifying every year, then run the dv/v setup then plot."** "Every year" = full record (all years 1995–2026), which also resolved the earlier 2010–2026-vs-1999 span question in favor of the full span.

This crash = **pod EVICTION at 77 GB** (not OOM; full pod restart). Recurring historical crashes were OOMs (firehose), now killed by the cap. memory.max≈187 GB, cpu.max=32.

Downstream (scripts generic, just `--station GNW --network UW`): densify → `build_long_window_daily_all51.py` → `dvv_coda_51.py --window 1.0 3.0`. `redpy.py` NOT in pipeline.

Sessions crash often — see [[crash-recovery-context]]; env can be wiped — see [[pod-env-ephemeral]].
