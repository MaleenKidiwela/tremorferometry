---
name: codex-lfe-parallel-effort
description: "User handed the LFE-catalog problem to ChatGPT/Codex (folder codex-LFE/, 2026-09-11..15) as an independent re-implementation + audit of my lfe_catalog_repair work; state, verdicts and the corrections it made to my tests"
metadata: 
  node_type: memory
  type: project
  originSessionId: d6a923e6-f9f1-4ee1-8d37-955b6b2e0f35
  modified: 2026-09-16T05:20:28.333Z
---

**2026-09-11..15: the user ran OpenAI Codex in `codex-LFE/` (untracked, ~21 GB: 16 GB cache, 3.7 GB results) as an independent, clean-room
attempt at the LFE catalog ([[lfe-catalog-goal-2026-09-10]]) and as an audit of my `lfe_catalog_repair/` E2-E6 work.** Read
`codex-LFE/docs/CLAUDE_TEST_AUDIT.md`, `docs/ACTIVE_WORK.md`, `docs/CATALOG_BUILD_PLAN.md`, `results/ACTIVE_STATUS.json` first. 146 tests pass.

**Corrections to MY work it found (accept these):** (1) E2 picker scorer `match_frac` was many-to-one (recall inflated up to 10.9 pp; 6.11%->4.19%
one-to-one); (2) E2 P rows were scored against S references; (3) E2 clock offset fitted on the evaluation period; (4) E2 borehole streams
(B011 az 11/101, B926 az 299/29) went into LFEDetect via 1->N/2->E aliasing WITHOUT rotation; (5) E4 divided by 34 days for a 35-day window;
(6) E6 GPU results (site-merge, one-per-site, noVGZSNB) completed 2026-09-11 13:24-14:02 UTC but were NEVER summarised in RESULTS.md
(only e6_checkpoint.txt + e6_* files). It also retracted nothing of E4/E5 that I had not already retracted.

**Its own verdicts (all "unvalidated", deliberately conservative):** roster = 183 stations from figures/borehole_dvv_map.html (NOT 198);
NLLB/LRIV/LZB/PFB/YOUB are OFF-roster (so my 15-station 2012 E5 catalog is not a roster result). Matched filter on 3-5 roster stations,
9 MAD, training-only 8-s templates: 3-9% recovery of Bostock references, controls <=5% mostly; "north candidate catalog v1" = 61 candidates /
17 days (2010). Frozen 2012 47-family bank FAILED its 5% gate (6.67%). Neural (jms/cascadia LFEDetect) + PyOcto on roster: ~0 reference
recovery; two fine-tuning attempts failed natural checks; 2026 central/northern panels do NOT separate from shift controls (one "seed" =
M6.0 Japan P arrival). Southern Ducellier templates on roster stations work on selected hours (7/8, 9/10 matches) — best regional lead.
Waveform-availability probes: Jan 2011 only 1-7 full-coverage stations/region; Aug 2026 6-43. Storage: monthly download/process/delete
runner built (not integrated); user wants 2011-01..2026-08, fallback 2019 start.

**Why:** the user said ChatGPT "gave a try to fix your mistakes" — treat codex-LFE as a peer review of my LFE work; its folder-boundary
rule (never write into lfe_catalog_repair) and provenance discipline (hashes, frozen thresholds, one-to-one scoring) are the bar to meet.
**How to apply:** before any further LFE work, reconcile the two efforts (my 15-station E5 result vs its roster-only result differ mainly
by station count/roster, not by contradiction); write up E6; adopt one-to-one scoring + rotation + roster restriction; do not re-run tests
listed in `results/prior_picker_test_audit_v2/reuse_ledger.json`.
