#!/bin/bash
# Stage 5 rerun in --single-cell mode (the 0.1-deg cells gave only ~35-pick fit stacks, too shallow
# to clear any gate), then B943.
cd /home/jovyan/tremorferometry
R=lfe_catalog_repair
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
export CUDA_PATH=/opt/conda/targets/x86_64-linux
PY=/home/jovyan/envs/tremorferometry/bin/python
wait_gpu(){ while pgrep -f "d3_redetect.py|d3_stage5_network.py" >/dev/null; do sleep 45; done; }
wait_gpu
$PY $R/scripts/d3_stage5_network.py --single-cell --min-picks 30 --snr-ratio 1.3 \
    --stations 'CN:PGC,PB:B011,CN:VGZ,CN:SNB,PB:B926,UW:LRIV,UW:OCP,CN:GOBB,UW:STW' \
    --start 2010-08-10 --end 2010-09-13 --lengths 2,4,6 --workers 5 \
    > $R/logs/d3_stage5b.log 2>&1
echo "[queue6 STAGE5B DONE]"
wait_gpu
bash $R/scripts/d3_run_redetect.sh --net PB --sta B943 \
    --templates-z data/b943_disc_p70_2010_2026_m3_40hz.npz \
    --templates-3c $R/results/d3_templates3c_B943.npz \
    --summary data/b943_disc_p70_2010_2026_m3.summary.csv \
    --family-list $R/results/d3_b943_wishshared_families.txt \
    --comps Z,1,2 --start 2017-01-01 --end 2019-01-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_B943 > $R/logs/d3_redet_b943.log 2>&1
echo "[queue6 ALL DONE]"
