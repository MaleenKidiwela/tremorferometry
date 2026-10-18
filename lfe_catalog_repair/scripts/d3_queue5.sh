#!/bin/bash
# D3 GPU queue (part 5) — final ordering after the stage-3 sensitivity budget.
# 1. STAGE 5 network matched filter on the 2010 ground-truth window (the decisive experiment):
#    templates from the project's REFINED Lin picks, 9 stations / 21 channels, L = 2/4/6 s,
#    k = 4..10. Predicted from the measured budget (s/n = 0.078, sigma_1ch(6 s) = 0.109):
#    sigma_net = 0.109/sqrt(21) = 0.0238 -> a median LFE reaches k ~ 4.9 MAD, so k = 4-5 should
#    detect and k >= 8 should not.
# 2. YACT (cheap, completes the single-station leg with D1's 0 %-responsive station).
# 3. B943, restricted to the 1,927 families in the 160 source cells it SHARES with WISH -- the only
#    ones the surface/borehole coincidence test uses (1.4 h instead of 2.5 h).
cd /home/jovyan/tremorferometry
R=lfe_catalog_repair
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
export CUDA_PATH=/opt/conda/targets/x86_64-linux
PY=/home/jovyan/envs/tremorferometry/bin/python
wait_gpu(){ while pgrep -f "d3_redetect.py|d3_stage5_network.py" >/dev/null; do sleep 45; done; }

wait_gpu
$PY $R/scripts/d3_stage5_network.py \
    --stations 'CN:PGC,PB:B011,CN:VGZ,CN:SNB,PB:B926,UW:LRIV,UW:OCP,CN:GOBB,UW:STW' \
    --start 2010-08-10 --end 2010-09-13 --lengths 2,4,6 --workers 5 \
    > $R/logs/d3_stage5.log 2>&1
echo "[d3 queue5 STAGE5 DONE]"

wait_gpu
bash $R/scripts/d3_run_redetect.sh --net UW --sta YACT --templates-z data/yact_disc.npz \
    --templates-3c $R/results/d3_templates3c_YACT.npz --summary data/yact_disc.summary.csv \
    --comps Z,N,E --start 2017-01-01 --end 2019-01-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_YACT > $R/logs/d3_redet_yact.log 2>&1
echo "[d3 queue5 YACT DONE]"

wait_gpu
bash $R/scripts/d3_run_redetect.sh --net PB --sta B943 \
    --templates-z data/b943_disc_p70_2010_2026_m3_40hz.npz \
    --templates-3c $R/results/d3_templates3c_B943.npz \
    --summary data/b943_disc_p70_2010_2026_m3.summary.csv \
    --family-list $R/results/d3_b943_wishshared_families.txt \
    --comps Z,1,2 --start 2017-01-01 --end 2019-01-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_B943 > $R/logs/d3_redet_b943.log 2>&1
echo "[d3 queue5 ALL DONE]"
