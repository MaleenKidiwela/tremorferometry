#!/bin/bash
# D3 stage 3 — validation, network coincidence, Lin ground truth, figure.
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
R=lfe_catalog_repair
set -x
for t in WISH YACT B943 PGC2010 PGC2017 PGC1516; do
  [ -f $R/results/d3_raw_${t}_detections.parquet ] || continue
  $PY $R/scripts/d3_validate.py --sta $t --prefix $R/results/d3_raw_${t} 2>&1 | tail -n 20
done
$PY $R/scripts/d3_network.py --a WISH --b B943 2>&1 | tail -n 12
$PY $R/scripts/d3_lin.py --det $R/results/d3_detections_PGC2010.parquet 2>&1 | tail -n 20
$PY $R/scripts/d3_figures.py --stations WISH,YACT,B943,PGC2010,PGC2017 \
   --coinc $R/results/d3_coincidence_WISH_B943.csv 2>&1 | tail -n 3
