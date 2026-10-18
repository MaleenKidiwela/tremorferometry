#!/bin/bash
# D3 stage 2 queue (part 3) — REORDERED after the WISH result.
# WISH came back with ~0 ETS response and ~0 tremor correlation under every rule, and the stage-4
# control showed the STORED times stack to nothing. The decisive question is therefore no longer
# "does the new rule beat the old one at WISH" but "does the new rule find the LFEs we KNOW are
# there" -> PGC over the Lin/Bostock ground-truth window is run FIRST.
cd /home/jovyan/tremorferometry
R=lfe_catalog_repair
wait_gpu(){ while pgrep -f "d3_redetect.py" >/dev/null; do sleep 45; done; }
run(){ bash $R/scripts/d3_run_redetect.sh "$@"; }

wait_gpu
run --net CN --sta PGC --templates-z data/pgc_pnsn_families_100km.npz \
    --templates-3c $R/results/d3_templates3c_PGC.npz \
    --summary data/pgc_pnsn_families_100km.summary.csv \
    --family-list $R/results/d3_pgc_lin_families.txt \
    --comps Z,N,E --start 2010-05-01 --end 2010-11-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_PGC2010 > $R/logs/d3_redet_pgc2010.log 2>&1
echo "[d3 queue3 PGC2010 DONE]"

wait_gpu
run --net UW --sta YACT --templates-z data/yact_disc.npz \
    --templates-3c $R/results/d3_templates3c_YACT.npz --summary data/yact_disc.summary.csv \
    --comps Z,N,E --start 2017-01-01 --end 2019-01-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_YACT > $R/logs/d3_redet_yact.log 2>&1
echo "[d3 queue3 YACT DONE]"

wait_gpu
run --net CN --sta PGC --templates-z data/pgc_pnsn_families_100km.npz \
    --templates-3c $R/results/d3_templates3c_PGC.npz \
    --summary data/pgc_pnsn_families_100km.summary.csv \
    --family-list $R/results/d3_pgc_lin_families.txt \
    --comps Z,N,E --start 2017-01-01 --end 2017-08-10 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_PGC2017 > $R/logs/d3_redet_pgc2017.log 2>&1
echo "[d3 queue3 PGC2017 DONE]"

wait_gpu
run --net PB --sta B943 --templates-z data/b943_disc_p70_2010_2026_m3_40hz.npz \
    --templates-3c $R/results/d3_templates3c_B943.npz \
    --summary data/b943_disc_p70_2010_2026_m3.summary.csv \
    --comps Z,1,2 --start 2017-01-01 --end 2019-01-01 --tb 64 --workers 6 \
    --out-prefix $R/results/d3_raw_B943 > $R/logs/d3_redet_b943.log 2>&1
echo "[d3 queue3 ALL DONE]"
