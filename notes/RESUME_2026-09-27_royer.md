# RESUME — Royer & Bostock (2014) reproduction (written 2026-09-27 ~03:45Z, before an account switch)

READ THIS FIRST. Supersedes RESUME_2026-09-26_royer.md. Full log: notes/2026-09-25_Notes.md (sections A–W).
Agent log: lfe_array_discovery/results/ab_royer/AB_ROYER_LOG.md. Work dir: lfe_array_discovery/.

## STATE 09-29 ~16:30Z (latest — read this block first; full log notes/2026-09-25_Notes.md 09-29 entries)
- CATALOGUE v1.1: data/catalogue/families_v1_1.csv — 610 families, 586 distinct, 508 located (A 371/B 83/C 54), 461 distinct with ≥1 permanent carrier.
  NW 141 (union of 5 runs, NW_union/), SW 81 (SW_mute+SW_lag15, SW_union/), NVI 19 (NVI_union/), Olympic Arm A 84 (either band) + Arm B 47 new
  (Chestler chain, scrambled control 0/50), Oregon 26 + 4 candidates, SVI 115, Ducellier 66, Shelly 27. Union script: lfe_array_discovery/scripts/nw_union.py <runs> <TAG>.
  v1.2 DONE: data/catalogue/families_v1_2.csv — 551 located, 496 distinct with ≥1 permanent carrier; Arm B 43/47 located.
- FILTER SLIP: Olympic site-count scores were 1–4 Hz (butter on 20 sps with fs-40 design). True 2–8 Hz: Arm A 62 real / control 1; only 31 pass both bands.
- SHELLY (lfe_south/): reproduced with his Text S1 recipe (notes/shelly_sources/) and all 17 of his stations (v3b): 9×MAD any-family recall 0.68, precision 0.78;
  same-family 0.49/0.55; results/SH_V3B_LOG.md. 27 families = his (none new).
- MENDOCINO DISCOVERY: STOPPED 09-30 by user (6 real, none of Shelly's 27 rediscovered, 5 "new" look like EQ repeaters). Use Shelly's 27 (v3b templates).
  markers results/MENDO_{DONE,FAILED,HEARTBEAT}, ETA ~2026-10-01 06Z. Location headline = MEN model (Yoon & Shelly 2024 Table S3; data/velocity_models/yoon_shelly_2024/,
  NCSN men.crt in data/velocity_models/ncsn_hypoinverse/) at Vp/Vs 1.78 (scripts/ab_hypo_mendo.py).
- ab_hypo.py "ZTR 35" → "ZTR 35 F" fixed 09-29 (was ignored by Hyp2000 1.40; tested on SW: no material change, locations stand).
- NEXT (GPU, after MENDO): extend Shelly's 27 families 2024-11 → 2026, or the WA trial extension episode — user to choose. Z5 batch-2 waits for user's Drive link.

## STATE 09-28 ~22:10Z (latest — read this block first)
- NW FINAL (user 09-29) = NW_lag15: 65 real, 64 located (1 unlocated: 2338, approx 47.86 N −123.00). NW_mute (66) superseded.
- SW 54 real (Royer 54); USGS model test: tighter errors, 4.5 km shallower than Royer (notes §AC).
- NVI 15 real (Royer 47).
- Oregon (lfe_oregon/): 8–26 families (paper 18), 5 located; per-event (Savard) and eye-pick routes failed (only 3–4 recording stations);
  looser station admission (chain_admit95/80) made it WORSE (5 and 3 families) → original rule kept. Better data = 3F nodes 2025 (embargoed).
- Olympic 2010: 53 real families (scrambled control 0), 46 located with the USGS NW model (median 34 km) → results/ab_royer/OLY2010_cvm/real_USGS.csv; 33 have a permanent carrier (B001, LRIV, GNW, VGZ).
- Extension plan (tree section 'Extension plan'): dedupe across catalogues → permanent carriers → permanent templates → trial episode 2011/12 → scan to 2026. Washington first.
- Velocity model: data/velocity_models/usgs_cascadia_v1.7/ (headline for new work; GSC comparison). Pod memory limit 200 GB; mem_watchdog running.
- Decision tree: notes/decision_tree/royer_path.html → https://claude.ai/artifact/WUf2L64jjdf2iQpU9cX6av (new account).

## User rules (also in memory)
- Follow Royer's method and GET FAMILIES; no side tests/extras unless asked. New locatable families count as success.
- Order: NW (done for now) → SW (RUNNING) → NVI → extension to 2026 (Phase 3).
- Say "scrambled control", never "null". Don't gate detections on ETS/tremor windows (labels only).
- Agents: launch-and-exit, TaskStop on hand-back. No git commits/pushes. Use the GPU. Be concise.

## RESULT SO FAR — NW (burst-muted run)
- results/ab_royer/NW_mute/: 299 grown families → **66 real LFE families** (review/class_final.csv), 65 located
  (hypo/loc_GSC_1.73_R1.5.csv), depth median 35.4 km, 41 inside 31–46 km, 36 with erh<5 km.
  vs Royer truth: 45/122 within 10 km (p=0.014), 27/122 within 5 km (p=0.066). Figure NW_mute/NW_mute_catalogue.png;
  locations real_R15.csv. Earlier runs gave 33 / 29 real (noise from W060/W020 bursts took over growth).
- Real-vs-noise classifier = number of independent sites with a short impulsive arrival in the 2–8 Hz stack
  (scratchpad famscore.py logic; B013/W030/B943 = one site); ≥3 sites = real. Matched eye labels 60/62 real, 149/149 noise.
  Saved as scripts/ab_famscore.py (score(npz) -> DataFrame grp,n_imp_sites; logic:
  per station, Hilbert envelope summed over 3 comps; impulsive if peak/median > 5 and half-max duration < 3 s; count distinct sites).

## SW RESULT (09-27 05:08Z)
- results/ab_royer/SW_mute/: 693 families → **54 real (Royer 54), all located** (real_R15.csv), depth median 35.9 km, 43 erh<5 km;
  vs Royer SW: 14/54 within 5 km (p=0.035), 28/54 within 10 km (p=0.028); depths agree to 0.6 km. Map SW_mute/WA_catalogue_NW_SW.png.
  Classifier with co-located site groups (stations <3 km apart) — see notes §Y.

## NVI RESULT (09-27 07:05Z)
- results/ab_royer/NVI_mute/: burst muting fixed the −52 % first-detection gap (4,235 vs 4,267). 40 families → **15 real (Royer 47), all located**
  (real_R15.csv), depth median 35.3 km; 10/47 within 5 km (p=0.033). TOTAL Royer reproduction: 66 NW + 54 SW + 15 NVI = 135.

## RUNNING NOW (superseded: SW and NVI finished)
- **SW muted run**: scripts/sw_mute_queue.sh, pid in results/SW_MUTE_PID (2762174), started 03:40Z, log logs/sw_mute_queue.log,
  markers results/SW_MUTE_{DONE,FAILED,HEARTBEAT}. Steps: SW mute masks → autocorr → cluster → ab_grow_fast (14 it) → ab_hypo →
  ab_restack28 → review sheets in results/ab_royer/SW_mute/hypo_28/review/. Expected ~3–4 h.

## WHEN SW_MUTE_DONE APPEARS
1. Score families with the site-count classifier (≥3 sites = real); eye-check the borderline (2–4 sites) and a random sample from sheets.
2. Take real families' locations from SW_mute/hypo/loc_GSC_1.73_R1.5.csv; score with
   `python scripts/ab_truth_score.py <csv with lon,lat,depth_km> --subarrays SW --nnull 1000` (Royer SW = 54 templates, 32–49 km).
3. Make a map + depth histogram like NW_mute/NW_mute_catalogue.png. Log in notes/2026-09-25_Notes.md.
4. If it FAILS: read the log; stages are guarded, relaunch with
   `cd lfe_array_discovery && nohup setsid -f bash scripts/sw_mute_queue.sh > logs/sw_mute_queue.log 2>&1 < /dev/null`
   (check nothing holds locks/sw_mute.lock first; kill only exact pids, never pkill -f).

## VELOCITY MODEL (09-28)
USGS Cascadia CVM v1.7 regional 1-D models: data/velocity_models/usgs_cascadia_v1.7/profiles/ (headline from now on; GSC = comparison; NVI stays GSC, outside grid).
Oregon chain running: lfe_oregon/results/chain_paper (paper reading, 189 candidates).

## MERLIN PLAN (09-28, user-approved; notes §AB)
1 no-merge growth test on 176 seeds → 2 full NW no-merge + final dedupe → 3 full-lag clustering → 4 SNR cull; then SW, NVI. Test output results/ab_royer/NW_nomerge/.

## NEXT LEVERS (older list)
- Most families are still noise (NW: 233/299), mainly low-amplitude W020 packets muting can't catch. Idea: use noise families as
  detectors and drop starting pairs that match them before growth.
- Tighten locations: eye-pick P/S on the 66 NW families (skill seismologist-pick).
- NVI (POLARIS PO, 2006-09 + 2007-06) with the same muted chain; then extension to 2026 (78 velocity stations in the WA box in 2025).

## Key code
ab_auto.py (autocorr), ab_cluster.py, ab_grow_fast.py (GPU, verified identical to ab_grow.py), ab_hypo.py (Hyp2000+hypoDD, /home/jovyan/tools/),
ab_restack28.py (2–8 Hz restack), ab_review_sheets.py, ab_mute.py (burst mute K=15 + line notch), ab_truth_score.py (Royer truth in results/ab_royer/truth/).
Main queue ab_royer_queue.sh is PAUSED; don't relaunch it (its end has an extension block that must not run unless asked).
- Decision tree: notes/lfe_decision_tree.json → python notes/decision_tree/build.py → republish https://claude.ai/artifact/Hw5t5SnuEYpX8UPgQQmFq1 (update SW node when it finishes).
