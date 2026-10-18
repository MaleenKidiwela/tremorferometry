#!/bin/bash
# P2 — the whole evaluation of one year in a single pass (so the 2011 referee budget is spent once).
#   bash scripts/p2_run_year.sh <model.pt> <tag> <year>
# Produces: CF per station-day, the 15-station back-projection, its moveout-shuffle null, and the
# 12-station (withheld B926/NLLB/LRIV) run used for the independent calibration check.
set -e
MODEL=$1; TAG=$2; YEAR=$3
PY=/opt/conda/bin/python
S12=PGC,VGZ,SNB,B001,B003,B004,B005,B006,B007,B009,B010,B011
cd /home/jovyan/tremorferometry/lfe_picker
$PY scripts/p2_cf.py --model "$MODEL" --year "$YEAR" --out "results/p2/cf_${TAG}_${YEAR}"
$PY scripts/p2_backproj.py --cf "results/p2/cf_${TAG}_${YEAR}" --year "$YEAR" --tag "${TAG}_${YEAR}" \
    --norm mad --site-delay results/p2/site_delay.json --maps-days "$4"
$PY scripts/p2_backproj.py --cf "results/p2/cf_${TAG}_${YEAR}" --year "$YEAR" --tag "${TAG}_${YEAR}_shuf" \
    --norm mad --site-delay results/p2/site_delay.json --shuffle --seed 11
$PY scripts/p2_backproj.py --cf "results/p2/cf_${TAG}_${YEAR}" --year "$YEAR" --tag "${TAG}_${YEAR}_s12" \
    --norm mad --site-delay results/p2/site_delay.json --stations "$S12"
echo "[p2_run_year] done $TAG $YEAR"
