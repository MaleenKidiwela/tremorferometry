#!/bin/bash
cd /home/jovyan/tremorferometry
R=lfe_catalog_repair
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
export CUDA_PATH=/opt/conda/targets/x86_64-linux
PY=/home/jovyan/envs/tremorferometry/bin/python
wait_gpu(){ while pgrep -f "d3_redetect.py|d3_stage5_network.py" >/dev/null; do sleep 45; done; }
# Stage 5c: the network matched filter itself, summing ALL 9 stations / 21 channels.
$PY $R/scripts/d3_stage5_network.py --single-cell --min-picks 30 --use-all --min-net-sta 3 \
    --stations 'CN:PGC,PB:B011,CN:VGZ,CN:SNB,PB:B926,UW:LRIV,UW:OCP,CN:GOBB,UW:STW' \
    --start 2010-08-10 --end 2010-09-13 --lengths 2,4,6 --workers 5 \
    > $R/logs/d3_stage5c.log 2>&1
echo "[queue7 STAGE5C DONE]"
# PGC-only control: the same machinery on 1 station / 3 channels, to isolate what the network buys.
wait_gpu
$PY $R/scripts/d3_stage5_network.py --single-cell --min-picks 30 --use-all --min-net-sta 1 \
    --stations 'CN:PGC' --start 2010-08-10 --end 2010-09-13 --lengths 6 --workers 4 \
    --outdir $R/results/pgconly > $R/logs/d3_stage5c_pgconly.log 2>&1
echo "[queue7 PGCONLY DONE]"
wait_gpu
bash $R/scripts/d3_run_redetect.sh --net PB --sta B943 \
    --templates-z data/b943_disc_p70_2010_2026_m3_40hz.npz \
    --templates-3c $R/results/d3_templates3c_B943.npz \
    --summary data/b943_disc_p70_2010_2026_m3.summary.csv \
    --family-list $R/results/d3_b943_wishshared_families.txt \
    --comps Z,1,2 --start 2017-01-01 --end 2019-01-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_B943 > $R/logs/d3_redet_b943.log 2>&1
echo "[queue7 ALL DONE]"
