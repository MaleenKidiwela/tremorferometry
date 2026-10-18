#!/bin/bash
# D3 stage 2 queue (part 2) — runs after the in-flight WISH pass, one station at a time on the GPU.
# SCOPE NOTE: PGC is restricted to its 2,129 families (of 5,640) whose source cell lies within 5 km
# of a Lin/Bostock family (`results/d3_pgc_lin_families.txt`, 287 distinct cells). Both reasons PGC
# is in this pilot are confined to that set: the Lin ground-truth test, and the stage-4 locations
# ("the region with ground truth"). It also cuts the two PGC passes from ~2.6 h to ~1 h of GPU.
# The 2015-16 bonus era is dropped for the same budget reason.
cd /home/jovyan/tremorferometry
R=lfe_catalog_repair
wait_gpu(){ while pgrep -f "d3_redetect.py" >/dev/null; do sleep 45; done; }
wait_file(){ while [ ! -f "$1" ]; do sleep 60; done; }
run(){ bash $R/scripts/d3_run_redetect.sh "$@"; }

wait_gpu
run --net UW --sta YACT --templates-z data/yact_disc.npz \
    --templates-3c $R/results/d3_templates3c_YACT.npz --summary data/yact_disc.summary.csv \
    --comps Z,N,E --start 2017-01-01 --end 2019-01-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_YACT > $R/logs/d3_redet_yact.log 2>&1
echo "[d3 queue2 YACT DONE]"

wait_file $R/results/d3_templates3c_B943.npz
wait_gpu
run --net PB --sta B943 --templates-z data/b943_disc_p70_2010_2026_m3_40hz.npz \
    --templates-3c $R/results/d3_templates3c_B943.npz \
    --summary data/b943_disc_p70_2010_2026_m3.summary.csv \
    --comps Z,1,2 --start 2017-01-01 --end 2019-01-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_B943 > $R/logs/d3_redet_b943.log 2>&1
echo "[d3 queue2 B943 DONE]"

wait_file $R/results/d3_templates3c_PGC.npz
wait_gpu
run --net CN --sta PGC --templates-z data/pgc_pnsn_families_100km.npz \
    --templates-3c $R/results/d3_templates3c_PGC.npz \
    --summary data/pgc_pnsn_families_100km.summary.csv \
    --family-list $R/results/d3_pgc_lin_families.txt \
    --comps Z,N,E --start 2010-05-01 --end 2010-11-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_PGC2010 > $R/logs/d3_redet_pgc2010.log 2>&1
echo "[d3 queue2 PGC2010 DONE]"

wait_gpu
run --net CN --sta PGC --templates-z data/pgc_pnsn_families_100km.npz \
    --templates-3c $R/results/d3_templates3c_PGC.npz \
    --summary data/pgc_pnsn_families_100km.summary.csv \
    --family-list $R/results/d3_pgc_lin_families.txt \
    --comps Z,N,E --start 2017-01-01 --end 2017-08-10 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_PGC2017 > $R/logs/d3_redet_pgc2017.log 2>&1
echo "[d3 queue2 ALL DONE]"
