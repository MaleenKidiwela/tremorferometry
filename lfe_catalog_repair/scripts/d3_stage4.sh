#!/bin/bash
# D3 stage 4 — network-consistent stacks -> family locations -> figure.
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
R=lfe_catalog_repair
# PGC itself is included: it anchors the frame (its predicted-S offset is 0 by
# construction, so its stack also checks that S really lands where the model puts it).
STNS='CN:PGC,CN:SYMB,PB:B011,CN:VGZ,UW:LRIV,CN:CLRS,CN:PTRF,CN:SNB,UW:OCP'
set -x
$PY $R/scripts/d3_stage4_stacks.py --pgc-det $R/results/d3_detections_PGC2017.parquet \
    --rule "${RULE:-new9}" --min-det "${MINDET:-200}" --stations "$STNS" \
    --start 2017-02-20 --end 2017-04-21 --workers 6 --min-snr 3.0 \
    --max-fams "${MAXFAMS:-0}" --out $R/results/d3_stage4_stacks 2>&1 | tail -n 25
$PY $R/scripts/d3_stage4_locate.py 2>&1 | tail -n 25
$PY $R/scripts/d3_stage4_fig.py 2>&1 | tail -n 3
