# lfe_extension_wa — WA LFE families on PERMANENT stations (trial extension)

Trial of step 4 of the extension plan for the Washington families (Royer NW/SW unions, Olympic 2010
Arm A + B): permanent-station templates from each family's own array-era detections, a validation
episode inside the array era scored against those detections (the known answer), and a trial episode
after the arrays (2012-09 ETS). Design and every deviation from the SVI recipe: `PLAN.md`. Log: `NOTES.md`.

## Run
```bash
nohup bash scripts/wa_trial_queue.sh > logs/queue.out 2>&1 &
cat results/WA_TRIAL_HEARTBEAT        # "<utc> <stage>", refreshed every 60 s
ls results/WA_TRIAL_{DONE,FAILED}     # end markers
kill $(cat results/WA_TRIAL_PID)      # stop (exact pid); resumable: finished stages have results/stage_<name>.done
```

## Scripts (in queue order)
| script | stage |
|---|---|
| `wa_common.py` | constants (frozen chain, windows, thresholds), helpers (clear-arrival rule, one-to-one matching, halves) |
| `wa_setup.py` | families (329 representatives), 42 permanent stations + 3-km sites, station-days, truth, scan days → `results/setup/` |
| `wa_fetch.py` | download missing station-days (reuses lfe_array_discovery / lfe_extension / lfe_catalog_repair traces) → `traces/` |
| `wa_dl.py` + `wa_dl_queue.sh` | PRODUCTION downloads, tiers A+B of `results/WA_PROD_PLAN.md` §4 (2012-01-01 → 2026-09-30, approved 2026-10-02) → NEW trace root `traces_dl/` (appended to `wa_common.TRACE_ROOTS`; only days missing from every other root; days held by other projects' roots are symlinked); manifest `results/wa_dl/manifest.csv`, plan `results/wa_dl/plan.csv`, per-station `results/wa_dl/completeness_by_station.csv`; log `logs/wa_dl.log` |
| `wa_dl2.py` + `wa_dl2_queue.sh` | PRODUCTION downloads, second queue (2026-10-02): tiers C (SW carriers) + D (UW 3-C ≥ 2017, ≤ 10 km) of `WA_PROD_PLAN.md` §4, the useful TEMPORARY-network stations (tier T, `results/wa_temp_networks/selected_for_download.csv`) and short UW deployments (D2), in that order, auto-stop at 2.0 TB → NEW trace root `traces_dl2/` (same `<NET.STA>/<YYYY>/<DOY>.mseed` layout; NOT yet in `wa_common.TRACE_ROOTS`: append it after the GO/NO-GO run; it holds only days missing from every other root; never a day of the A+B plan); manifest `results/wa_dl2/manifest.csv`, per-station `results/wa_dl2/completeness_by_station.csv`; markers `results/wa_dl2/WA_DL2_{PID,HEARTBEAT,ALIVE,DONE,FAILED}`; log `logs/wa_dl2.log` |
| `wa_temp_inv.py` | Temporary-network inventory 2012-01 → 2026-09 within 75 km of a WA family (FDSN EarthScope/NCEDC/NRCan; one 10-min sample per station-year) → `results/wa_temp_networks/{INVENTORY.md,inventory.csv,by_network.csv,channels.csv,samples.csv,selected_for_download.csv}` |
| `wa_cache.py` | 20 Hz 1–8 Hz midnight-aligned cache with validity mask and node id → `cache20hz/` |
| `wa_templates.py` | permanent templates per (family, station node) from build-half detections; clear-arrival QC; SNR weights → `results/templates/` |
| `wa_mf.py` | GPU network matched filter + site-rolled scrambled control → `results/mf_<tag>/` |
| `wa_hf.py` | HF ratio per detection at the family's HF station, native rate → `results/mf_<tag>/hf.parquet` |
| `wa_score.py` | flags (tiers, HF veto, cross-family), HF cut calibration, validation / trial tables → `results/score_<tag>/` |
| `wa_newnodes.py` | trial station epochs that need re-templating + pass-1 events → `results/trial/` |
| `wa_report.py` | `results/WA_TRIAL_RESULTS.md`, `figures/wa_trial_summary.png` |

Smoke-test outputs (6 days of 2010-08, not part of the result): `results/mf_smoke`, `results/score_smoke*`,
`results/templates/*_smoke*`.

## Step 4b — added stations (NW + Olympic)
```bash
nohup bash scripts/step4b_queue.sh > logs/step4b_queue.out 2>&1 &
cat results/STEP4B_HEARTBEAT; ls results/STEP4B_{DONE,FAILED}; kill $(cat results/STEP4B_PID)
```
Scripts: `step4b_setup.py` (candidates, selection, stations, station-days, truth → `results/step4b/`),
`step4b_events.py` (bridge templating events + nodes), `step4b_report.py` (`results/STEP4B_RESULTS.md`),
`step4b_queue.sh`; reuses `wa_fetch/cache/templates/mf/hf/score` with `WA_SETUP=results/step4b/setup`,
`WA_EXTRA_NODES=results/step4b/setup/extra_nodes.csv`. Stage guards `results/step4b/stage_<name>.done`.
Smoke outputs (5 days; not the result): `results/mf_step4b_smoke_*`, `score_step4b_smoke_*`, `results/step4b/smoke/`,
`results/templates/*step4b_smoke*`. Design: `PLAN.md` §Step 4b.

## Stack activity (family activity by stacking; runs after / beside step 4b)
```bash
nohup setsid -f bash scripts/stack_activity_queue.sh > logs/stack_activity_queue.out 2>&1 < /dev/null
cat results/STACK_ACT_HEARTBEAT; ls results/STACK_ACT_{DONE,FAILED}; kill $(cat results/STACK_ACT_PID)
```
Scripts: `stack_activity.py` (engine: leave-one-site-out top-K candidates, 2–8 Hz daily stacks for data / 3 scrambled
controls / 6 random-times draws), `stack_activity_score.py` (windows W = 1/5/10/30 d, z, A, p_ctl, truth),
`stack_activity_select.py` (continuous-series families + stations + setup), `stack_activity_report.py`
(`results/stack_activity/STACK_ACTIVITY_RESULTS.md`, `figures/stack_activity_*.png`). Guards
`results/stack_activity/stage_<name>.done`. Smoke outputs (not results): `results/stack_activity/run_smoke*`, `score_smoke*`.
Design: `PLAN.md` §Stack activity.

## Merlin tests (T1 all templates, T2 withheld-station template recovery; 2026-09-30)
```bash
nohup setsid -f bash scripts/merlin_tests_queue.sh > logs/merlin_tests_queue.out 2>&1 < /dev/null
cat results/MERLIN_T_HEARTBEAT; ls results/MERLIN_T_{DONE,FAILED}; kill $(cat results/MERLIN_T_PID)
```
Scripts: `merlin_t1_setup.py` (merged setup incl. A2011 truth → `results/merlin_tests/t1_setup/`), `wa_mf.py --all-templates
--drop-dead 0.05`, `wa_hf.py --all-templates`, `merlin_t1_score.py` (own / ≤3 km / any levels, separation), `merlin_t2.py`
(`--stage setup|round|score`), `merlin_tests_report.py` → `results/MERLIN_TESTS_RESULTS.md`. MF outputs `results/mf_mt1_*`;
everything else `results/merlin_tests/`. Guards `results/merlin_tests/stage_<name>.done`; chain failures `T1_FAILED` / `T2_FAILED`.
Smoke outputs (not results): `results/mf_mt1smoke_all`, `results/merlin_tests/t1/score_mt1smoke_all`, `results/merlin_tests/t2_smoke`.
Design: `PLAN.md` §Merlin's two deciding tests.

## STEP 2 — group detection + member attribution (Merlin, 2026-09-30)
```bash
nohup setsid -f bash scripts/step2_groups_queue.sh > logs/step2_queue.out 2>&1 < /dev/null
cat results/STEP2_HEARTBEAT; ls results/STEP2_{DONE,FAILED}; kill $(cat results/STEP2_PID)
```
Scripts: `step2_groups.py` (groups, alignment, group setup, group templates), `wa_mf.py`/`wa_hf.py` (tag step2grp, step2grp_c34),
`step2_score.py` (group vs T1 summed over members), `step1_features.py`/`step1_fit.py` (STEP1_* env → group P(real)),
`step2_attrib.py` (cc at the nearest sites, softmax attribution, expected counts), `step2_report.py`. Outputs `results/step2_groups/`,
report `results/STEP2_GROUPS_RESULTS.md`. Smoke outputs (not the result): `results/mf_step2grp_smoke`, `results/step2_groups/smoke/`.

## PNET as per-repeat evidence (Merlin step 6, 2026-09-30; queue `scripts/pnet_queue.sh`)
```bash
nohup setsid -f bash scripts/pnet_queue.sh > logs/pnet_queue.out 2>&1 < /dev/null
cat results/PNET_HEARTBEAT; ls results/PNET_{DONE,FAILED}; kill $(cat results/PNET_PID)
```
Scripts: `pnet_common.py` (ZNE rotation, label rule, held-out stations), `pnet_cf.py` (PNET CF per station-day), `pnet_feat.py`
(PNET features per step-1 candidate), `pnet_auc.py` (stage A/B feature AUCs), `pnet_labels.py` + `pnet_finetune.py` (stage B),
`pnet_eval.py` (per-station AUC / quiet-day FP on held-out stations and days), `step1_fit.py` with `STEP1_EXTRA_FEATS` (stage C),
`pnet_t2.py` + `merlin_t2.py --stage score` with `MT2_OD` (stage D), `pnet_report.py` → `results/PNET_RESULTS.md`.
Outputs `results/pnet/` (cf_<tag>/, feat_<tag>/, train/, models/, step1_<tag>/, t2_<tag>/). Guards `results/pnet/stage_<name>.done`;
failed / blocked stages listed in `results/pnet/failed_stages.txt` (copied to `results/PNET_FAILED`). Smoke outputs (not results):
`results/pnet/{cf_smoke,train_smoke,t2_smoke,step1_smoke}`, `results/pnet/models/smoke.pt`. Design: `PLAN.md` §PNET.
