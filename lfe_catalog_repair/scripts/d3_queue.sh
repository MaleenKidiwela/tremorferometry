#!/bin/bash
# D3 stage 2 queue: one station at a time on the single GPU, after its 3-C templates exist.
cd /home/jovyan/tremorferometry
R=lfe_catalog_repair
wait_gpu(){ while pgrep -f "d3_redetect.py" >/dev/null; do sleep 60; done; }
wait_file(){ while [ ! -f "$1" ]; do sleep 60; done; while pgrep -f "d3_build_3c_templates.py --net $2" >/dev/null; do sleep 60; done; }
run(){ bash $R/scripts/d3_run_redetect.sh "$@"; }

# ---- WISH (redone: the tb=256 pass OOM'd on the L40S at day ~250; tb=64 is the same speed) ----
wait_gpu
bash lfe_catalog_repair/scripts/d3_run_redetect.sh --net UW --sta WISH \
    --templates-z data/wish_disc.npz --templates-3c lfe_catalog_repair/results/d3_templates3c_WISH.npz \
    --summary data/wish_disc.summary.csv --comps Z,N,E --start 2017-01-01 --end 2019-01-01 \
    --tb 64 --workers 6 --out-prefix lfe_catalog_repair/results/d3_raw_WISH \
    > lfe_catalog_repair/logs/d3_redet_wish.log 2>&1
echo "[d3 queue WISH DONE]"

# YACT (templates ready)
wait_gpu
run --net UW --sta YACT --templates-z data/yact_disc.npz \
    --templates-3c $R/results/d3_templates3c_YACT.npz --summary data/yact_disc.summary.csv \
    --comps Z,N,E --start 2017-01-01 --end 2019-01-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_YACT > $R/logs/d3_redet_yact.log 2>&1

# B943 borehole (components Z,1,2)
wait_file $R/results/d3_templates3c_B943.npz PB
wait_gpu
run --net PB --sta B943 --templates-z data/b943_disc_p70_2010_2026_m3_40hz.npz \
    --templates-3c $R/results/d3_templates3c_B943.npz \
    --summary data/b943_disc_p70_2010_2026_m3.summary.csv \
    --comps Z,1,2 --start 2017-01-01 --end 2019-01-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_B943 > $R/logs/d3_redet_b943.log 2>&1

# PGC 2010 ground-truth window (Lin/Bostock overlap)
wait_file $R/results/d3_templates3c_PGC.npz CN
wait_gpu
run --net CN --sta PGC --templates-z data/pgc_pnsn_families_100km.npz \
    --templates-3c $R/results/d3_templates3c_PGC.npz \
    --summary data/pgc_pnsn_families_100km.summary.csv \
    --comps Z,N,E --start 2010-05-01 --end 2010-11-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_PGC2010 > $R/logs/d3_redet_pgc2010.log 2>&1

# PGC 2017 BH era (stage-4 reference detections + tremor test)
wait_gpu
run --net CN --sta PGC --templates-z data/pgc_pnsn_families_100km.npz \
    --templates-3c $R/results/d3_templates3c_PGC.npz \
    --summary data/pgc_pnsn_families_100km.summary.csv \
    --comps Z,N,E --start 2017-01-01 --end 2017-08-10 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_PGC2017 > $R/logs/d3_redet_pgc2017.log 2>&1

# PGC 2015-16 ETS (bonus era)
wait_gpu
run --net CN --sta PGC --templates-z data/pgc_pnsn_families_100km.npz \
    --templates-3c $R/results/d3_templates3c_PGC.npz \
    --summary data/pgc_pnsn_families_100km.summary.csv \
    --comps Z,N,E --start 2015-11-01 --end 2016-04-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_PGC1516 > $R/logs/d3_redet_pgc1516.log 2>&1
echo "[d3 queue ALL DONE]"
