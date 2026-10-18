#!/bin/bash
# Run once after a pod reboot (written 2026-10-12, before the 10-13 reboot). Idempotent: safe to run twice.
# 1) restarts the memory watchdog, 2) checks the Python environment, 3) restarts the S11 WA queue if it had not
# finished (stages with DONE markers are skipped), 4) prints the state of every S11 job and where to read next.
# Agent-driven jobs (ballistic test, REDPy pilot) cannot be restarted from bash: the next Claude session re-launches
# them from their PLAN.md (see notes/RESUME_2026-10-12_dvv_depth.md §2).
R=/home/jovyan/tremorferometry; cd $R
S=$R/svi_science/results/s11_depth
TS() { TZ=America/Los_Angeles date '+%Y-%m-%d %-I:%M %P PT'; }
echo "== resume_after_reboot $(TS)"

# 1) watchdog
if pgrep -f "scripts/mem_watchdog.sh" >/dev/null; then echo "watchdog: already running"
else nohup bash scripts/mem_watchdog.sh >/dev/null 2>&1 & echo "watchdog: started pid $!"; fi

# 2) environment (pod rebuilds can wipe /opt/conda; user packages live in ~/.local)
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')$LD_LIBRARY_PATH
python3 -c "import numpy, scipy, pandas, pyarrow, obspy, matplotlib, cartopy; print('python imports: ok')" \
  || { echo "python imports: FAILED -> pip install --user <missing> (see memory pod-env-ephemeral). STOPPING: nothing restarted."; exit 1; }
python3 -c "import obspy.clients.fdsn, decorator" || { echo "obspy import FAILED (2026-10-14: 'decorator' was wiped) -> pip install --user decorator. STOPPING."; exit 1; }

# 3) S11 WA-NW queue (+ WA-SW chain)
if grep -q QUEUE_DONE $S/wa_nw/logs/queue.log 2>/dev/null; then echo "WA-NW queue: finished (QUEUE_DONE)"
elif pgrep -f s11_wa_queue.sh >/dev/null; then echo "WA-NW queue: running"
else
  echo "WA-NW queue: NOT finished -> restarting (DONE stages skipped)"
  (cd svi_science/scripts && S11_REGION=WA_NW nohup bash s11_wa_queue.sh >> $S/wa_nw/logs/queue.log 2>&1 &)
fi
if grep -q SW_CHAIN_DONE $S/wa_sw/logs/queue.log 2>/dev/null; then echo "WA-SW column: finished"
elif pgrep -f s11_wa_sw_chain.sh >/dev/null; then echo "WA-SW column: running/waiting"
else
  echo "WA-SW column: not finished -> restarting chain (waits for WA-NW)"
  (cd svi_science/scripts && nohup bash s11_wa_sw_chain.sh >> $S/wa_sw/logs/queue.log 2>&1 &)
fi

# 4) agent-driven pilots: report only
for d in ballistic redpy; do
  n=$(ls $S/$d/*DONE 2>/dev/null | wc -l); r=$([ -f $S/$d/RESULTS.md ] && echo yes || echo no)
  echo "S11 $d: DONE markers $n, RESULTS.md $r -> if no RESULTS.md, re-launch an agent from $S/$d/PLAN.md"
done
# 5) S11 round 3 queue (PLAN.md §6; added 2026-10-12 7:35 pm PT): kernels -> measure -> kerr -> r3q3 -> zones, SVI + WA-NW
R3=$S/round3
if [ -f $R3/QUEUE_DONE ]; then echo "S11 round 3: finished (QUEUE_DONE)"
elif pgrep -f s11_r3_queue.sh >/dev/null; then echo "S11 round 3: running (log $R3/queue.log)"
else
  echo "S11 round 3: NOT finished -> restarting (DONE stages skipped; log $R3/queue.log)"
  (cd svi_science/scripts && nohup bash s11_r3_queue.sh > /dev/null 2>&1 &)
fi
# 6) S11 ballistic round 3 queue (ballistic/round3/PLAN.md §5; added 2026-10-12 7:42 pm PT): measure SVI -> measure WA-NW
#    -> analysis (s11_ball3_analysis.sh), each stage with a DONE marker; idempotent
B3=$S/ballistic/round3
if [ -f $B3/QUEUE_DONE ]; then echo "S11 ballistic round 3: finished (QUEUE_DONE)"
elif pgrep -f s11_ball3_queue.sh >/dev/null; then echo "S11 ballistic round 3: running (log $B3/logs/queue.log)"
else
  echo "S11 ballistic round 3: NOT finished -> restarting (DONE stages skipped; logs $B3/logs/)"
  (cd svi_science/scripts && nohup bash s11_ball3_queue.sh > /dev/null 2>&1 &)
fi
if [ -f $B3/QUEUE_DH_DONE ]; then echo "S11 ballistic round 3 DH adaptation: finished"
elif pgrep -f s11_ball3_queue_dh.sh >/dev/null; then echo "S11 ballistic round 3 DH adaptation: running"
else echo "S11 ballistic round 3 DH adaptation: NOT finished -> restarting"; (cd svi_science/scripts && nohup bash s11_ball3_queue_dh.sh > /dev/null 2>&1 &); fi
# 7) S11 tomo Vp/Vs pilot, SVI 2003-06 era stacks (tomo/PLAN.md Amendment 8; added 2026-10-12 8:06 pm PT);
#    per-station npz files are skipped when present; afterwards an agent runs s11_tomo_era_picks + wadati/profile
TM=$S/tomo
if [ -f $TM/ERA_STACKS_DONE ]; then echo "S11 tomo era stacks: finished (ERA_STACKS_DONE)"
elif pgrep -f s11_tomo_era_stacks.py >/dev/null; then echo "S11 tomo era stacks: running (log $TM/logs/era_stacks.log)"
else
  echo "S11 tomo era stacks: NOT finished -> restarting (done stations skipped)"
  (cd svi_science/scripts && OMP_NUM_THREADS=1 nohup python s11_tomo_era_stacks.py --workers 6 >> $TM/logs/era_stacks.log 2>&1 &)
fi
echo "READ NEXT: notes/RESUME_2026-10-12_dvv_depth.md, then $S/NOTES.md (decision log)"
# 8) S11 distant-station lever (s11_depth/distant/PLAN.md; added 2026-10-12 8:42 pm PT): per region inventory -> windows ->
#    measure -> kernels -> invert -> gates, then figures; DONE markers in distant/<svi|wa_nw>/; idempotent restart
DS=$S/distant
if [ -f $DS/QUEUE_DONE ]; then echo "S11 distant lever: finished (QUEUE_DONE)"
elif pgrep -f s11_dist_queue.sh >/dev/null; then echo "S11 distant lever: running (log $DS/queue.log)"
else
  echo "S11 distant lever: NOT finished -> restarting (DONE stages skipped; log $DS/queue.log)"
  (cd svi_science/scripts && nohup bash s11_dist_queue.sh > /dev/null 2>&1 &)
fi
# 9) S11 round 4 (PLAN.md §8 + Amendment 1; added 2026-10-12 ~9:50 pm PT): per region measure (near/ext/far) -> far kernels ->
#    lapse-selection recovery -> combined model + injection suite; then the 2019-11 dive (SVI) and figures; DONE markers in
#    s11_depth/round4/; idempotent restart (stages with markers are skipped)
R4D=$S/round4
if [ -f $R4D/QUEUE_DONE ]; then echo "S11 round 4: finished (QUEUE_DONE)"
elif pgrep -f s11_r4_queue.sh >/dev/null; then echo "S11 round 4: running (log $R4D/queue.log)"
else
  echo "S11 round 4: NOT finished -> restarting (DONE stages skipped; log $R4D/queue.log)"
  (cd svi_science/scripts && nohup bash s11_r4_queue.sh > /dev/null 2>&1 &)
fi
# 10) S11 round 5 (PLAN.md §9 + Amendment 1; added 2026-10-12 ~10:45 pm PT): task-1 forms per region -> selection -> smooth
#     Sc + Pcal + Q3 band -> Sc adoption -> Q1 (WA yearly cycle) per region -> 2019-11 per instrument (SVI) -> figures;
#     DONE markers in s11_depth/round5/; idempotent restart (stages with markers are skipped)
R5D=$S/round5
if [ -f $R5D/QUEUE_DONE ]; then echo "S11 round 5: finished (QUEUE_DONE)"
elif pgrep -f s11_r5_queue.sh >/dev/null; then echo "S11 round 5: running (log $R5D/queue.log)"
else
  echo "S11 round 5: NOT finished -> restarting (DONE stages skipped; log $R5D/queue.log)"
  (cd svi_science/scripts && nohup bash s11_r5_queue.sh > /dev/null 2>&1 &)
fi
# 11) S11 round 6 (PLAN.md §10; added 2026-10-12 ~11:09 pm PT): task-1 forms (E4y/E4l/E4h/E4y+Scs2 + E3/E4 columns) and the
#     T1 cause analysis per region -> adoption -> task-2 ΔD-vs-episode-size regressions per region -> figures;
#     DONE markers in s11_depth/round6/; idempotent restart (stages with markers are skipped)
R6D=$S/round6
if [ -f $R6D/QUEUE_DONE ]; then echo "S11 round 6: finished (QUEUE_DONE)"
elif pgrep -f s11_r6_queue.sh >/dev/null; then echo "S11 round 6: running (log $R6D/queue.log)"
else
  echo "S11 round 6: NOT finished -> restarting (DONE stages skipped; log $R6D/queue.log)"
  (cd svi_science/scripts && nohup bash s11_r6_queue.sh > /dev/null 2>&1 &)
fi
# 11b) S11 round 6 §10 Amendment 1 (added ~11:25 pm PT): forms E4lu/E4yu -> adoption -> task-2 column; markers in round6/
if [ -f $R6D/AMEND_QUEUE_DONE ]; then echo "S11 round 6 amendment: finished (AMEND_QUEUE_DONE)"
elif pgrep -f s11_r6a_queue.sh >/dev/null; then echo "S11 round 6 amendment: running (log $R6D/queue.log)"
else
  echo "S11 round 6 amendment: NOT finished -> restarting (log $R6D/queue.log)"
  (cd svi_science/scripts && nohup bash s11_r6a_queue.sh > /dev/null 2>&1 &)
fi
# 13) S11 tomo/eq earthquake-pick fetch for the differential Wadati (tomo/PLAN.md §6; added 2026-10-13 ~9:45 am PT):
#     ComCat phase-data QuakeML -> data/eq_picks/raw/comcat/ (cached events skipped, .fail retried) -> parse -> picks.parquet;
#     markers data/eq_picks/{FETCH_DONE,PARSE_DONE}, PID data/eq_picks/fetch.pid, log tomo/eq/logs/fetch.log; idempotent
EQ=$R/data/eq_picks
if [ -f $EQ/PARSE_DONE ]; then echo "S11 eq-pick fetch: finished (PARSE_DONE) -> analysis s11_eq_wadati.py / s11_eq_diff.py"
elif pgrep -f s11_eq_fetch.sh >/dev/null; then echo "S11 eq-pick fetch: running (log $S/tomo/eq/logs/fetch.log)"
else
  echo "S11 eq-pick fetch: NOT finished -> restarting (cached events skipped)"
  (cd svi_science/scripts && nohup bash s11_eq_fetch.sh >> $S/tomo/eq/logs/fetch.log 2>&1 &)
fi
#     analysis chain (waits for PARSE_DONE): κ per zone -> differential -> figures; markers tomo/eq/{KAPPA,DIFF,ANALYSIS}_DONE
if [ -f $S/tomo/eq/ANALYSIS_DONE ]; then echo "S11 eq analysis: finished (ANALYSIS_DONE) -> write tomo/eq/RESULTS.md"
elif pgrep -f s11_eq_analysis.sh >/dev/null; then echo "S11 eq analysis: waiting/running (log $S/tomo/eq/logs/analysis.log)"
else echo "S11 eq analysis: restarting chain"; (cd svi_science/scripts && nohup bash s11_eq_analysis.sh >> $S/tomo/eq/logs/analysis.log 2>&1 &); fi
# 12) S11 round 7 (PLAN.md §11 + Amendment 1; added 2026-10-13 ~9:45 am PT): per region front timing (day stacks vs days since
#     the front) -> zones with E4 (per major, slipped vs not) -> joint ballistic S−P + coda -> WA hint halves / SVI
#     replication; then the hint verdict and figures; DONE markers in s11_depth/round7/; idempotent restart
R7D=$S/round7
if [ -f $R7D/QUEUE_DONE ]; then echo "S11 round 7: finished (QUEUE_DONE)"
elif pgrep -f s11_r7_queue.sh >/dev/null; then echo "S11 round 7: running (log $R7D/queue.log)"
else
  echo "S11 round 7: NOT finished -> restarting (DONE stages skipped; log $R7D/queue.log)"
  (cd svi_science/scripts && nohup bash s11_r7_queue.sh > /dev/null 2>&1 &)
fi
# 14) S12 rf_gosselin (receiver-function reproduction of Gosselin et al. 2020; svi_science/results/s12_lvl/rf_gosselin/PLAN.md;
#     added 2026-10-13 ~10:05 am PT): teleseismic P-window download (CN from NRCAN, UW/US from EarthScope, <=3 parallel),
#     cached per event in data/rf_teleseis/raw/<NET>.<STA>/, per-chunk markers _chunks/*.done, station markers
#     raw/<NET>.<STA>.DONE, all done = raw/DOWNLOAD_DONE; log rf_gosselin/logs/download_*.log; idempotent restart
RFW=$R/svi_science/results/s12_lvl/rf_gosselin
if [ -f $R/data/rf_teleseis/raw/DOWNLOAD_DONE ]; then echo "S12 RF download: finished (DOWNLOAD_DONE)"
elif pgrep -f rf01_download.py >/dev/null; then echo "S12 RF download: running (logs $RFW/logs/download_*.log)"
else
  echo "S12 RF download: NOT finished -> restarting (done chunks skipped)"
  (cd $RFW && nohup python -P scripts/rf01_download.py --workers 3 > logs/download_main.out 2>&1 & echo $! > logs/download.pid)
fi
#     14b) S12 rf_gosselin processing queue (waits for the download markers): RFs -> repro (1993-2019) -> 0.7 Hz -> 2019-26
#     out-of-sample -> full -> recovery curves -> extension SVI/WA -> slip-front timing; markers rf_gosselin/*_DONE,
#     all done = rf_gosselin/QUEUE_DONE; log rf_gosselin/logs/queue.log; idempotent restart
if [ -f $RFW/QUEUE_DONE ]; then echo "S12 RF queue: finished (QUEUE_DONE) -> write RESULTS.md + target.json"
elif pgrep -f rf_queue.sh >/dev/null; then echo "S12 RF queue: running (log $RFW/logs/queue.log)"
else echo "S12 RF queue: NOT finished -> restarting"; (cd $RFW && nohup bash scripts/rf_queue.sh > /dev/null 2>&1 & echo $! > logs/queue.pid); fi
# 15) S12 LVL, LFE side (svi_science/results/s12_lvl/lfe/PLAN.md; added 2026-10-13 ~10:10 am PT): LVL-base reflection scan
#     -> superposed-epoch coda dv/v / decay / P-S per lapse window vs a between-ETS baseline, SVI then WA-NW -> pooling;
#     DONE markers in s12_lvl/lfe/, log s12_lvl/lfe/logs/queue.log; idempotent restart (stages with markers are skipped)
L12=$R/svi_science/results/s12_lvl/lfe
if [ -f $L12/QUEUE_DONE ]; then echo "S12 LVL LFE side: finished (QUEUE_DONE) -> figures: python $L12/scripts/s12_figures.py"
elif pgrep -f s12_queue.sh >/dev/null; then echo "S12 LVL LFE side: running (log $L12/logs/queue.log)"
else
  echo "S12 LVL LFE side: NOT finished -> restarting (DONE stages skipped; log $L12/logs/queue.log)"
  (cd $L12/scripts && nohup bash s12_queue.sh > /dev/null 2>&1 &)
fi
if [ -f $L12/QUEUE2_DONE ]; then echo "S12 LVL LFE side A-D: finished (QUEUE2_DONE)"
elif pgrep -f s12_queue2.sh >/dev/null; then echo "S12 LVL LFE side A-D: running"
elif [ -f $L12/QUEUE_DONE ]; then echo "S12 LVL A-D: restarting"; (cd $L12/scripts && nohup bash s12_queue2.sh > /dev/null 2>&1 &)
fi
#     15b) PLAN §7 (added ~10:40 am PT): WA recharge confirmation, SVI by era, 1-3 Hz reflection restack, figures
if [ -f $L12/QUEUE3_DONE ]; then echo "S12 LVL recharge/low-band: finished (QUEUE3_DONE)"
elif pgrep -f s12_queue3.sh >/dev/null; then echo "S12 LVL recharge/low-band: running (log $L12/logs/queue.log)"
elif [ -f $L12/QUEUE2_DONE ]; then echo "S12 LVL recharge/low-band: restarting"; (cd $L12/scripts && nohup bash s12_queue3.sh > /dev/null 2>&1 &)
fi
# 16) S13 VLFE-band (20-50 s) stacks on LFE times (svi_science/results/s13_vlfe/PLAN.md; added 2026-10-13 ~10:45 am PT):
#     pilot queue = 2011-08 + 2012-08 pooled/family stacks with 500 scrambled-control draws -> figures + pilot_verdict.json;
#     markers s13_vlfe/STACK_pilot_*_DONE, PILOT_QUEUE_DONE; log s13_vlfe/logs/stack_pilot.log; idempotent restart.
#     Scale-up queue (if launched) = scripts/s13_scale_queue.sh, markers SCALE_*_DONE, all done = SCALE_QUEUE_DONE.
S13=$R/svi_science/results/s13_vlfe
if [ -f $S13/PILOT_QUEUE_DONE ]; then echo "S13 VLFE pilot: finished (PILOT_QUEUE_DONE) -> tables/pilot_verdict.json"
elif pgrep -f s13_pilot_queue.sh >/dev/null; then echo "S13 VLFE pilot: running (log $S13/logs/stack_pilot.log)"
else echo "S13 VLFE pilot: NOT finished -> restarting"; (cd $S13 && nohup bash scripts/s13_pilot_queue.sh > logs/pilot_queue.out 2>&1 & echo $! > logs/pilot_queue.pid); fi
if [ -f $S13/scripts/s13_scale_queue.sh ]; then
  if [ -f $S13/SCALE_QUEUE_DONE ]; then echo "S13 VLFE scale-up: finished (SCALE_QUEUE_DONE)"
  elif pgrep -f s13_scale_queue.sh >/dev/null; then echo "S13 VLFE scale-up: running (log $S13/logs/scale_queue.log)"
  else echo "S13 VLFE scale-up: NOT finished -> restarting"; (cd $S13 && nohup bash scripts/s13_scale_queue.sh > logs/scale_queue.out 2>&1 & echo $! > logs/scale_queue.pid); fi
fi
# 17) S12c ACF under PGC (svi_science/results/s12_lvl/acf_pgc/PLAN.md; added 2026-10-13 ~11:14 am PT): hourly whitened
#     auto/cross-correlations from raw mseed on disk (no download). Pilot 2018-2021: scripts/acf_pilot_queue.sh (PILOT_PRE_DONE)
#     + scripts/acf_pilot_ana.sh (PILOT_ANA_DONE). Scale-up 2003-2026: scripts/acf_full_queue.sh (FULL_PRE_DONE, FULL_VIS_DONE,
#     EPOCH_full_*_DONE, all done = FULL_QUEUE_DONE; log logs/full_queue.log). Finished station-years are skipped on restart.
A12=$R/svi_science/results/s12_lvl/acf_pgc
#     17b) amendment B echo amplitude/similarity stage: scripts/acf_amp_queue.sh (waits FULL_QUEUE_DONE; AMP_full_*_DONE,
#     all done = AMP_QUEUE_DONE, appends an auto section to RESULTS.md)
for Q in acf_pilot_queue:PILOT_PRE_DONE acf_pilot_ana:PILOT_ANA_DONE acf_full_queue:FULL_QUEUE_DONE acf_amp_queue:AMP_QUEUE_DONE; do
  N=${Q%%:*}; M=${Q##*:}
  if [ -f $A12/$M ]; then echo "S12c ACF $N: finished ($M)"
  elif pgrep -f $N.sh >/dev/null; then echo "S12c ACF $N: running"
  else echo "S12c ACF $N: NOT finished -> restarting"; (cd $A12 && nohup bash scripts/$N.sh > logs/$N.restart.log 2>&1 &); fi
done
#     16b) S13 gap fill (Amendment 3, user OK 10-13): fetch -> wait SCALE_QUEUE_DONE -> prep -> scalef_* restacks -> summary;
#     markers s13_vlfe/data/raw/FETCH_fill_DONE, PREP_FILL_DONE, STACK_scalef_*_DONE, FILL_QUEUE_DONE; log logs/fill_queue.log
if [ -f $S13/FILL_QUEUE_DONE ]; then echo "S13 gap fill: finished (FILL_QUEUE_DONE) -> tables/scale_verdict.json"
elif pgrep -f s13_fill_queue.sh >/dev/null; then echo "S13 gap fill: running (log $S13/logs/fill_queue.log)"
else echo "S13 gap fill: NOT finished -> restarting"; (cd $S13 && nohup bash scripts/s13_fill_queue.sh > logs/fill_queue.out 2>&1 & echo $! > logs/fill_queue.pid); fi
