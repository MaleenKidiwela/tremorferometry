# RESUME — LFE array discovery (Cascadia) — written 2026-09-21 03:55 UTC
Read this first after any reset or account change. Everything below is on disk under `/home/jovyan/tremorferometry/`.

## 1. One-paragraph state
We reproduced Bostock et al. (2012/2015) LFE-family discovery on his own 2005 POLARIS stations using his PUBLISHED recipe
(network autocorrelation → grouping → stack → iterative matched filter), blind: **58/96 open ("kept") families and 22/30 sealed
("hidden") families rediscovered** (ceilings 89/96 and 29/30 from catalogue-built templates). Every earlier failure was OUR
machinery (added gates, a toned station, midnight-only quiet hours, whole-second catalogue times, a 0.5-s grid that cancels stacks,
a merge that fused neighbours). The same recipe on PERMANENT stations discovers 3/77 (2012) and 5/60 (2010) while templates from
known events run 33/77 and 50/60 → arrays FIND families, permanent stations only RUN them. The reproduced recipe is now running,
unchanged, on the Washington 2010 Array-of-Arrays/CAFE data (N-Olympic) with location-based verification. Seven conservative
candidate NEW northern families exist (one spatially distinct), unconfirmed.

## 2. What is running right now (check first)
**UPDATE 07:45Z:** the WA 2010 run FINISHED (`results/WA_2010_RESULTS.md` + addendum): raw screen 50 real vs 17 null candidates, but 46/50 rail on a grid edge; interior 4 vs 3; strict (≥5 days, co-timing excess ≥0.05) **2 real (+1 marginal) vs 0 null**: grp 1656 (48.04 N, −123.92, 30 km), grp 1193 (48.16 N, −123.52, 30 km), grp 1152. Agent M launched 07:50Z for the sealed-2011 re-fire + wider-grid relocation (markers `results/M_REFIRE2011_DONE|FAILED`, log `logs/m_refire2011_queue.log`). If M's marker is absent and no `m_refire2011_queue.sh` process exists, relaunch it the same way as L.

- `lfe_array_discovery/scripts/l_wa2010_queue.sh` — pid 1222985 (setsid-detached; flock `locks/l_wa2010.lock`); log
  `lfe_array_discovery/logs/l_wa2010_queue.log`; markers `lfe_array_discovery/results/L_WA2010_DONE|FAILED`; writes
  `results/WA_2010_RESULTS.md`. Started 2026-09-21 ~02:40Z, expected ≈6 h. Stage at 03:50Z: growth merged (204,403 detections);
  remaining: locator → screen → whole-chain site-shift null → report. If the marker is absent and no `l_wa2010_queue.sh`
  process exists, relaunch: `cd lfe_array_discovery && nohup setsid -f bash scripts/l_wa2010_queue.sh >> logs/l_wa2010_queue.log 2>&1 < /dev/null`
  (per-stage guards make it resume).
- Nothing else. All downloads complete (`results/DL_ALL_DONE`). No agents alive (agents from the old session cannot be resumed from
  another account — spawn fresh ones with the file pointers below).

**UPDATE 09:15Z:** the sealed-2011 re-fire FINISHED (`results/WA_2011_REFIRE_RESULTS.md`): **0/50 real and 0/17 null confirmed**; the real templates fire 3× more on non-tremor than tremor days and relocate unstably → cultural/site repeaters, not LFE families. N-Olympic 2010 yields no Washington family. Agent N launched ~09:20Z for the last full-array WA bed: Puget/S-Olympic 2008 (CAFE 44 sites + PB), with a per-station tremor-hour audit (the 2008 network-energy hours were 86 % daytime) and 2007-01 as the confirmation re-fire bed; markers `results/N_PUGET2008_DONE|FAILED`, log `logs/n_puget2008_queue.log`. **Launched 09:50Z, pid 1271745, ≈2–3 h**; its hour audit reconstructed the May-2008 ETS from waveforms alone (115 tremor hours on 2008-05-05..26, 0.9 % daytime; 0/456 quiet hours fire); tremor/quiet pair ratio 8.2× (2010: 4.2×). If its marker is absent and no `n_puget2008_queue.sh` process exists, relaunch the same way (resumable).

**UPDATE 10:20Z — GOAL RESTATED BY THE USER:** reproduce the PUBLISHED Washington families (Royer & Bostock 2014 down to the Oregon border on CAFE/PB; Chestler & Creager 2017 ≈40 Olympic families on the Array of Arrays; 69 derived coordinates in `codex-LFE/results/olympic_reference_discovery_v1/`). The 2010 null is therefore a reproduction FAILURE to debug (suspect: 64/96 stations in 7 tight clusters treated independently → cluster noise dominates; Chestler used the arrays as arrays). Agent O launched 10:20Z: array-stacked 2010 re-run + Chestler check + 2011 re-fire (markers `results/O_OLY2010_DONE|FAILED`, log `logs/o_oly2010_queue.log`). Agent N (Puget 2008 + 2007 confirmation) still running. PLAN.md §17.

**UPDATE 10:30Z:** Puget 2008 FINISHED (`results/WA_2008_RESULTS.md`): 0 real / 0 null candidates; pair-stage tremor/quiet ratio 0.89; grown templates fire equally on quiet and tremor days → noise repeaters. WA surface stations are noise-repeater dominated. Agent P launched ~10:35Z: per-station noise-repeater census (quiet-hour pair detections per station) → quiet-subset 2008 bed (PB boreholes + quietest CAFE) → recipe → screen → 2007 confirmation (markers `results/P_PUGET2008Q_DONE|FAILED`, log `logs/p_puget2008q_queue.log`). **Launched 11:05Z, pid 1298829, ≈1.5–2 h.** Smoke: CAFE surface stations ~14 dB louder than PB boreholes on quiet days; a provisional 16-station quiet subset gives a pair tremor/quiet ratio of 10.8× (vs 0.89× on all 51; north 11×). Agent O launched 10:55Z: pid 1297832, ≈7–8 h, markers `results/O_OLY2010_DONE|FAILED`, output `results/WA_2010_ARRAY_RESULTS.md` (Chestler's own chain from her UW thesis runs as Arm B; PLAN.md §18).

**UPDATE 11:20Z:** Puget 2008 quiet-subset FINISHED (`results/WA_2008Q_RESULTS.md`): pair ratio 0.86 (boreholes-only 0.45) — station quality is NOT the lever; 0 candidates. Agent Q (screen validation on the 58 Bostock families → `results/SCREEN_VALIDATION_NORTH.md`) and agent O (2010 arrays-as-arrays + Chestler's own chain, pid 1297832) still running. Open: what exactly Royer & Bostock 2014 found in WA (fetch Royer's UBC PhD thesis, open access at UBC cIRcle).

**UPDATE 12:00Z — MAJOR CORRECTION:** `results/SCREEN_VALIDATION_NORTH.md`: the WA screen as coded keeps only 2 of the 58 blind-rediscovered Bostock families; interior-grid, re-stack stability, 8-km Slab2, shuffled control and co-timing clauses are UNJUSTIFIED. Justified tier A = ≥5-station location + |Δz| ≤ 19 km + ≥3 days + quiet FA ≤ 0.25/template-day (& ≥3× below tremor). Re-screened: WA 2008 → 17 real vs 2 null (tier A), 4 vs 0 (tiers B/C); WA 2010 → 234 vs 59 (B), 127 vs 27 (C). The 2007/2011 confirmation rules (interior + 10 km) were equally unjustified (a true family's re-stack moves a median 35 km). Agent R launched ~12:05Z: calibrate confirmation on the north (half-split re-fire, true vs noise), re-confirm 2008→2007 and 2010→2011 with justified tiers, Chestler check → `results/WA_CANDIDATES_v2.md` (markers `results/R_CONFIRM_DONE|FAILED`). PLAN.md §19.

## 3. The recipe (binding; do not add clauses — see §6)
From Bostock 2015 §2 / Plourde 2015 §2 / Brown 2008, parameter table in `lfe_array_discovery/results/REPRO_BOSTOCK.md` §0:
6-s windows, 1–8 Hz (frozen 20-Hz cache is equivalent to raw 100 Hz — measured), per-component cc summed over ALL admitted 3-C
stations (no site merge, no per-station z, no clip), threshold 8×MAD of the block's network-cc distribution (≈5σ), 4-s decluster,
Bungum–Husebye realignment before every stack (REQUIRED: any 0.5-s/1-s time grid gives |sinc| = 0 at 2/4/6/8 Hz), transitive
grouping (4-s node merge, ≥3 members), stack → template cut at stack onset → 8×MAD network matched filter over tremor days →
iterate (pass 0 is best; iterating loses families — §7), NO waveform-similarity merge (families 0.1–1.6 km apart co-detect; resolve
duplicates by location), scoring vs truth = one-to-one ±2 s, ≥30 % of the grown family's detections (Bostock's criterion).
Refuted by measurement: 6×MAD first pass (0 families), survival gates, quiet-day FA-only thresholds, coincidence clauses, site merge.
Code: `lfe_array_discovery/scripts/k_*.py` (k_auto, k_group, k_grow, k_locate, k_common, k_s1_loud for sub-second catalogue times);
WA bed binding `l_common.py`/`l_run.py`. Older `w_*`, `x_*`, `z_*`, `y_*`, `q_*` scripts are the superseded gated engine — do not reuse.

## 4. Key results files
- `lfe_array_discovery/results/NORTH_DISCOVERY_SUMMARY.md` — the north result for the user (plain language + tables).
- `lfe_array_discovery/results/REPRO_BOSTOCK.md` — the full ledger §0–§11 (parameters, every defect, every measurement).
- `lfe_array_discovery/results/repro/candidate_new_families.csv` — 18 candidates; conservative 7 (grp 534, 1215, 912, 889, 247, 745,
  822); only grp 1215 (47.90 N, −123.66, 32 km) is >20 km from any catalogued family. NOTHING CLAIMED.
- `lfe_array_discovery/PLAN.md` §13–§14 (mandate change + WA design), `NOTES.md` (dated ledger of every launch/result).
- `lfe_picker/results/P4_0_RESULTS.md`, `P4_0B_RESULTS.md` (+ corrections) — the closed permanent-network Puget attempt.
- `lfe_extension/results/trackD/RESULTS.md` — the 130-family extension to 2026 on the permanent roster (v1.1), the one earlier product.
- Merlin reviews (verbatim) in `lfe_array_discovery/results/merlin_*.md` and `lfe_picker/results/merlin_*.md` — historical; the
  gated framework they prescribe is superseded by the user's mandate (§6).

## 5. Data on disk (lfe_array_discovery/)
`traces/` + `cache20hz/` (20 Hz 1–8 Hz ZNE npz per station-day): WA N-Olympic 2010-07-20..09-10 (T2A–T2D; XG 7 sites @50 Hz =
chain=rational, XU, PB, UW, CN) + 20 quiet days (T2B); Puget/S-Olympic 2008-04-15..06-10 (T4; ZK is vertical-only in the data);
2007-01..02 (T5); intraslab EQ windows (`traces_eq/`, `results/intraslab_eq_2006_2012.csv`); the 2007–2026 extension record for 17
Olympic permanents (T7). `traces_sealed/` = 2011-07-20..09-10, hashed (`results/sealed_2011_manifest.sha256`), NEVER read — the
re-fire test for WA families. North beds: `traces/bed2005_polaris/` + cache (2005-08-29..09-21 + 20 quiet days, 14 stations; SHVB has
a 4.76-Hz tone — exclude), `lfe_catalog_repair/cache20hz` (2010/2012 permanents). Station tables:
`results/channels_2006_2012.csv`, `station_flags.csv`, `sites_array_era_final.csv`, `site_probe_2010.csv` (admitted stations),
`wa_station_hour_veto_2010.csv`/`_2008.csv`, `wa_quiet_hours_2010.csv`/`_2008.csv`, `P_STEP2_SITES.md`, `WA_V2_INPUTS.md`.
Bostock catalogue: `lfe_catalog_repair/results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv` (2004–2005 times are
whole-second). Hidden families: `lfe_picker/results/holdouts.json` (30; sealed — one blind check done, do not iterate on them).

## 6. Rules the user set (binding)
- Bostock is reproducible on his stations; a failed rediscovery = OUR bug. Never say "not reproducible" / "method doesn't work".
- Use the published recipe literally; do not invent clauses, bars or gates; change one thing at a time, justified by a measurement
  against known families. Merlin only for a final blind check, never mid-development.
- Agents: launch long jobs as `nohup setsid -f bash queue.sh >> log 2>&1 < /dev/null` from a FOREGROUND call, exit immediately,
  no idle waiting/polling/interim turns; per-stage output guards; flock; signal traps writing FAILED markers; `env VAR=… python`
  per stage (an exported K_YEAR once leaked into a blind score). Never `pkill -f` (it killed the orchestrator's relaunches).
  Spawn fresh agents once a transcript passes ~200k tokens instead of resuming. Quota is the scarce resource (Opus and Fable).
- No author contact/emails. Don't delete results without explicit OK. Document as you go (PLAN/NOTES/results + memory).

## 7. Next steps, in order
1. When `L_WA2010_DONE` appears: read `results/WA_2010_RESULTS.md` — candidate WA families (located ≥5 stations, shuffled control
   pass, |Δz| ≤ 8 km Slab2, ≥3 days, tremor co-timing > base rate) vs the site-shift null count; Chestler-derived coordinates
   (`codex-LFE/results/olympic_reference_discovery_v1/derived_coordinate_rows_not_a_catalog.csv`) as a descriptive comparison only.
2. Confirm candidates by RE-FIRE, never by re-tuning: WA → run the unchanged templates on the SEALED 2011 window (open it once;
   ≥3 days at 8×MAD with a reproducing location); north candidates → download 2004-07 POLARIS (~900 station-days; POLARIS ran
   2004-04→2005-10; 59,659 Bostock events) and re-fire the 2005 templates.
3. Puget/S-Olympic 2008 (CAFE 44 sites + PB; ZK is vertical-only): audit the hour selection first — the 2008 array-energy candidate
   hours are 86 % daytime (cultural noise); use a tremor-timing source or a quiet-day contrast per station before running.
4. Extension: register confirmed WA families on the permanents (stack at detection times; templates from each permanent's own
   stack), then run 2007→2026 with the lfe_extension frozen chain (FA-only tiers per network epoch, completeness by station-set
   ablation) — the machinery that produced the 130-family v1.1 catalogue.
5. Parked: front estimator v2 (`lfe_extension/results/trackD/RESULTS.md` §3b-2; de-duplicate co-detected events across families
   ≤3 km first); far-south existence test (Shelly/Ducellier families on B045/B935/RVIT/PETY; ~55 GB download, needs user OK).

## 8. How to spawn the next agent (template)
Opus, fresh context. Give it: this file, `lfe_array_discovery/PLAN.md` §13–§14, `results/REPRO_BOSTOCK.md` §0 (recipe),
`results/NORTH_DISCOVERY_SUMMARY.md`, the k_/l_ scripts, the marker/log paths in §2, and the rules in §6. Ask for a ≤10-line report.

## 9. Lessons / memory (in case the auto-memory does not follow the account)
The auto-memory lives at `/home/jovyan/.claude/projects/-home-jovyan-tremorferometry/memory/` (keyed by project path on this
machine). A full copy as of 2026-09-21 is in `notes/memory_snapshot_2026-09-21/` — read `MEMORY.md` there first, then the files it
indexes. The lessons that matter most for this line: `reproduce-known-answer-first.md`, `agents-must-not-idle-wait.md`
(+ its detach addendum), `minimize-fable-delegate-opus.md` (+ fresh-agent addendum), `scoring-truth-colocated-families.md`,
`dont-over-conclude.md`, `check-geometry-before-density-claims.md`, `lfe-extension-study-2026-09-16.md` (the dated status ledger).
If the new session's memory directory is empty, copy the snapshot back into it.

## 10. FINAL BOOKKEEPING (2026-09-21 11:20 UTC — account at 98 % of the weekly limit; session may stop abruptly)

### 10.1 Running now (no agent attached; all survive the session)
- `scripts/o_oly2010_queue.sh` pid 1297832 — Olympic 2010 arrays-as-arrays run + Chestler's own chain (Arm B) + sealed-2011 re-fire. Started 10:55Z, ≈7–8 h. Log `lfe_array_discovery/logs/o_oly2010_queue.log`; markers `results/O_OLY2010_DONE|FAILED`; output `results/WA_2010_ARRAY_RESULTS.md`. NOTE: its screen is agent N's UNJUSTIFIED one — re-screen its templates with the tier-A/B/C rules of `results/SCREEN_VALIDATION_NORTH.md` before reading any count.
- Agent R was STOPPED by the orchestrator at 11:27Z (user: let jobs run without agents; a new agent reads the results) while still building `scripts/r_confirm_queue.sh` — its queue was NOT launched. A new agent must redo its brief (below) (calibrate the confirmation rule on the north half-split; re-confirm 2008→2007 and 2010→2011 under tiers A/B/C; Chestler check). If `results/R_CONFIRM_DONE|FAILED` and `scripts/r_confirm_queue.sh` are both absent, the agent was cut off before launching — a new agent should redo its brief (PLAN.md §19 + the agent-R paragraph in NOTES.md) from `results/SCREEN_VALIDATION_NORTH.md`. If the script exists but no process and no marker: relaunch it detached (resumable).

### 10.2 What is established (numbers traceable to files)
- North: Bostock's published recipe reproduces his families blind — kept 58/96 (ceiling 89), hidden 22/30 (ceiling 29), isolated 6/10; permanent-station beds discover 3/77 (2012) and 5/60 (2010) vs ceilings 33 and 50 → arrays find families, permanents only run them. 7 conservative candidate NEW northern families (grp 534, 1215, 912, 889, 247, 745, 822; 1215 spatially distinct) — unconfirmed; confirmation = re-fire in 2004-07 POLARIS (download ~900 station-days). `results/NORTH_DISCOVERY_SUMMARY.md`, `results/REPRO_BOSTOCK.md`.
- Washington: under the JUSTIFIED screen (tier A) Puget 2008 has 17 real vs 2 null candidates (4 vs 0 tiers B/C) and Olympic 2010 has 127–234 vs 27–59 (thresholds to re-derive per bed; quiet-hour FA never run on the 2010 templates). Earlier "0 candidates" verdicts (WA_2010_RESULTS, WA_2011_REFIRE, WA_2008, WA_2008Q) were produced by unjustified clauses and are SUPERSEDED as verdicts (their pair-stage and census measurements stand). Pair-stage tremor/quiet ratio: north 11×, Olympic 2010 4.2×, Puget 2008 0.86–0.89× (boreholes 0.45×) — the WA pair stage is noise-repeater rich; the recipe still grows candidates that pass the justified screen at 4–8× the null rate.
- Royer & Bostock 2014's WA family count/region/stations could not be obtained (all full-text sources blocked); Chestler & Creager 2017's method is in `results/WA_2010_ARRAY_RESULTS.md` §0; 69 derived Chestler coordinates are the only WA location check.

### 10.3 Next actions, in order (for whoever resumes)
1. Read `results/O_OLY2010_DONE` / `WA_2010_ARRAY_RESULTS.md`; re-screen with tiers A/B/C; compare member vs beam vs siteavg pair ratios; read the Arm-B (Chestler chain) result and the Chestler-coordinate check.
2. Read `results/R_CONFIRM_DONE` / `WA_CANDIDATES_v2.md` (or redo agent R's brief): calibrated confirmation rule; confirmed WA families real vs null for 2008 (via 2007) and 2010 (via 2011).
3. Any confirmed WA family → register on the permanents (stack at detection times; template from each permanent's own stack) → extension 2007→2026 with the lfe_extension frozen chain (data already on disk: T7).
4. North candidates → download 2004-07 POLARIS and re-fire (the only north step needing a download).
5. Parked: front estimator v2; far-south existence test (Shelly/Ducellier; ~55 GB, needs user OK).
6. Puget 2007 (CAFE full) as a further discovery bed only if 2008 confirms families.

### 10.4 Rules reaffirmed today (see §6 and `notes/memory_snapshot_2026-09-21/`)
Literal published recipe; every added criterion must pass ≥85 % of the 58 rediscovered Bostock families and reject noise, else dropped (`justify-criteria-on-known-families.md`); a failed rediscovery is our bug, never "not reproducible"; agents launch detached queues and exit; never `pkill -f`; fresh agents over resuming >200k-token transcripts; document as you go.

## 11. STATE AS OF 2026-09-23 05:00 UTC (supersedes §10 where they differ)
**Running (detached queues, no agent attached):** `lfe_south/scripts/t_south_queue.sh` pid 1412422 (Ducellier reproduction: download → score event-by-event; markers T_SOUTH_REPRO_A_DONE | T_SOUTH_STAGE1_DONE | T_SOUTH_REPRO_DONE | T_SOUTH_DONE | T_SOUTH_FAILED); `lfe_array_discovery/scripts/u_locate_queue.sh` pid 1428881 (rebuilt locator → re-locate WA real/null/ArmB → Chestler comparison; markers U_LOCATE_DONE|FAILED; report results/WA_LOCATION_REBUILD.md). `s_final_queue.sh` FINISHED (S_FINAL_DONE) — its report `results/WA_FINAL_ASSESSMENT.md` has a WRONG verdict, see PLAN §21.
**Headline findings of the night:** see `lfe_array_discovery/PLAN.md` §21 and NOTES.md 2026-09-23 entries. Washington's 2010 Olympic array bed contains repeating sources the null does not (280 templates with ≥20 dedupe partners vs 0); the locator was the failure (27.6 km → 4.0 km rebuilt); all WA location columns must be recomputed.
**Binding rule added today:** rank, do not discard (`memory/rank-dont-discard.md`) — never filter candidates; every criterion is a column; the null gets the identical table; evidence is the real-vs-null separation curve. The ≥85 %-of-known-families test is a VETO on proposed cuts, never a licence.
**Next three actions:** (1) on U_LOCATE_DONE, recompute every WA table with the 4-km locator; (2) audit the direction flag of every column in sf_separation.py and rewrite the verdict; (3) redo the Chestler-coordinate comparison with bootstrap radii — the test of whether those 280 repeating templates are Chestler's 43 families, new families, or something else. Then the southern reproduction score decides whether the chain is trustworthy at all.

## 12. HANDOFF BEFORE A CLAUDE CODE UPDATE — written 2026-09-25 ~04:00 UTC
Read §12 first; it supersedes §10 and §11 where they differ. Written because the user is updating Claude Code
and the session/agents may not survive. Everything needed to continue is on disk.

### 12.1 THE ONE THING THAT CHANGED EVERYTHING: the papers were on disk all along
Eight full-text papers sit in `notes/*.pdf` (Royer & Bostock 2014 + its supplement, Bostock 2012,
Bostock 2015, Chestler & Creager 2017, Thomas & Bostock 2015, Plourde 2015, Shelly 2025, Sweet 2014),
plus Ducellier's complete manuscript source in `catalogs/ducellier_catalog/paper/`. For two days the
orchestrator asserted they were unavailable and reconstructed parameters from secondary sources.
**Two documents now hold the extraction and MUST be read before touching any reproduction:**
- `notes/PAPER_METHODS_TABLE.md` — 12 parameter sections, verbatim quotes with section cites, a column
  for what OUR code does with file:line, and a ranked mismatch list (M1…M12).
- `notes/CATALOGUE_SUMMARY.md` — reader-facing, one section per catalogue north→south, what they found,
  whether they LOCATED families, what they released, whether we use their method, status.
- `notes/MARGIN_CATALOG_PLAN_2026-09-24.md` — the phase plan (phases 1–4, decision points D1–D4) and the
  rule that a published POSITIVE result is itself testable without the authors' catalogue.

### 12.2 WHO LOCATED WHAT (the orchestrator was wrong about this and must not be again)
Royer & Bostock 2014: **122 NW + 54 SW + 47 NVI** templates, **ALL located** (Hyp2000 → hypoDD, Table S1;
retained 77 % NW, 72 % SW, 91 % NVI, 92 % SVI); depths **31–46 km NW, 32–49 km SW, 33–39 km NVI**.
Bostock 2012: **140**, all located (hypoinverse), 27–37 km, deliberately no dedupe.
Chestler & Creager 2017: **57 → 43**, 32,387 LFEs, all located (3-D Preston grid search + Waldhauser–Ellsworth DD).
Thomas & Bostock 2015: **18**, >15,000 detections, only **13 located**, just 4 with free depths (31–34 km).
Plourde 2015: 59 templates → **34 families**, all located, 28–47 km. Bostock 2015: 130, did NOT locate.
Ducellier 2022: Plourde's 66 → 60 kept, did NOT locate (adopted Plourde's).

### 12.3 WHY WASHINGTON FAILED — four faults, all ours, now specified
1. **52 stations vs their 7 NAMED.** NW: B001, SQM, BS11, W020, W040, GNW, PL11. SW: N050, S030, W070,
   N060, PL11, S040, N070. **SQM, GNW, W020 (NW) and N070 (SW) are NOT ON DISK — download first.**
2. **6-s windows at 20 sps vs their 15-s at 40 sps** (`k_auto.py:42`, `k_common.py:28`); our docstring cited
   the 6 s to Bostock 2012, a different region.
3. **Connected-component grouping (`k_group.py:51-65,96-99`) vs their hierarchical clustering** on a
   super-trace similarity matrix (Rowe et al. 2002). Components chain; clustering does not.
4. **Their second stage was never run at all** — network cross-correlation with subsets of 7–10 high-SNR
   3-C stations, 8×MAD over a 24-HOUR period, restack, iterate to convergence. That is what turns 224
   rough clusters into 122 final families.
Their episodes: NW 2007-01-17..30, 2008-05-01..23, 2010-08-07..24, 2011-08-04..23; SW 2007-01-14..02-01,
2008-05-01..25. Stage yields to score against: NW 4,915 detections → 224 clusters → 122 families;
SW 3,306 → 70 → 54.
**Therefore the two Puget zeros (`WA_2008_RESULTS.md`, `WA_2008Q_RESULTS.md`) measured a method nobody
published, and the `y_wa_finish` run was killed mid-flight for the same reason. Washington has never been
attempted with the method that produced its 176 families.**

### 12.4 A BUG CLASS THAT HAS BITTEN TWICE AND IS NOT YET AUDITED IN WASHINGTON
The southern Shelly run produced **12,848,734 detections against his 61,441**. Three faults, all ours:
(a) `lfe_south/scripts/t11_mf.py` `ncc()` divides by `sd=sqrt(max(M*S2-S1*S1,1e-30))` while
`t3_cache.py:31-32` **zero-fills gaps with no validity mask** → constant runs give sd≈1e-15 and
"normalised" cc of 2.4e15 (one family/day: 128,922 detections vs his 71); (b) cc series summed at a
**common lag** although per-station offsets span 1.49–5.85 s within a family; (c) the recall denominator
used all 61,441 events against 385 days holding 5,021, so the cap was 0.082 and "0.026" was 0.32 of it.
Fixed in COPIES (`t11_mf_aa.py` etc.): same day went **128,922 → 21 detections, 11 of them his**.
Ducellier's path is immune (validity mask, gapless blocks, one absolute cut time) — hence 0.836 vs 0.026.
**OWED: `t11_mf.py` still carries both bugs, so `SOUTH_RESULTS.md` §5 extension counts are inflated.**
**OWED: the same audit on the WASHINGTON chain** (`k_common.py`, `k_auto.py`, `k_grow.py`, `o_auto.py`,
`o_beam.py`, `uu_*.py`). If present, the 2010 Olympic "61 of 69", the 280 dedupe-rich templates and the
2008 zeros are ALL suspect. Agent AB was asked to do this read-only audit.

### 12.5 CORRECTIONS TO OUR OWN CLAIMS
- The Olympic "61 of 69 coordinates matched vs 3 for the null" is coverage against
  `codex-LFE/results/olympic_reference_discovery_v1/derived_coordinate_rows_not_a_catalog.csv` — it is
  **NOT 61 of Chestler's 43 families** and must be restated.
- Our 4-s declustering comes from Bostock **2015**'s magnitude-catalogue cull, not from any detector's
  published rule.
- The rebuilt locator is real and validated: **4.0 km** median epicentral error against Bostock's
  published coordinates vs **27.6 km** for the old one (chance 39.8 km), rail 0.75→0.09, 11→35 of 58
  within 10 km. Adopted `L9_xcdiff` + `L12_comb`. Double-difference measured WORSE (7.40 vs 6.44 km) and
  was NOT adopted. `results/WA_LOCATION_REBUILD.md`.
- Bostock's published family depths sit a median **9.7 km below** this project's Slab2 field, so the old
  "|Δz| ≤ 8 km of Slab2" clause was wrong on its face.

### 12.6 WHAT IS RUNNING AND WHAT SURVIVES AN UPDATE
Both are `setsid`-detached shell queues and should survive a CLI update but NOT a pod restart. Both are
per-stage guarded, so relaunching resumes rather than redoing.
- `lfe_south/scripts/t_south_queue.sh` pid **1509294** — Ducellier matched filter (66 families, 34
  stations, 460 templates). Markers `T_SOUTH_DONE|FAILED`; heartbeat `T_SOUTH_HEARTBEAT`.
- `lfe_south/scripts/aa_shelly_queue.sh` pid **1800120** — the fixed Shelly reproduction, ablation rows
  A/B/C/D. Markers `AA_SHELLY_DONE|FAILED`; report `results/aa_shelly/AA_SHELLY_REPORT.md`.
- Relaunch either with `cd <workspace> && nohup setsid -f bash scripts/<queue>.sh >> logs/<queue>.log 2>&1 < /dev/null`.
- A standing watch writes to `notes/queue_watch.log` when a queue dies or finishes; re-arm it after an update.
- **Agent AB is writing the Washington code but was told NOT to run it** — see `results/AB_ROYER_PLAN.md`
  (its launch command is the last line of that file). Nothing Washington executes until the user says so.

### 12.7 THE NEXT ACTIONS, IN ORDER, AFTER THE UPDATE
1. Read `results/AB_ROYER_PLAN.md` and its bug audit; if the bug is present in the Washington chain, fix
   it before any reproduction and re-run what it invalidates.
2. Download the 4 missing Royer stations; run their method on NW and SW, scoring at every stage checkpoint.
3. Then NVI (47 families, same paper, never attempted).
4. Re-run the southern extension with the fixed matched filter.
5. Oregon (Thomas & Bostock) last — and read `notes/PAPER_METHODS_TABLE.md`'s Thomas row first.
