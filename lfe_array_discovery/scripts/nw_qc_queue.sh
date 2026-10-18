#!/usr/bin/env bash
# NW growth rerun with station QC (2026-09-26, user-approved): W060 and FACB removed from the growth/location
# pool (126/136 eye-classified noise families had W060 among their growth stations; FACB 61/136). Same
# cluster stacks, same ab_grow (duplicate merge, gap guard), then ab_hypo. Markers results/NW_QC_{DONE,FAILED,HEARTBEAT}.
set -u -o pipefail
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
exec 9>locks/nw_qc.lock; flock -n 9 || exit 0
R=results; OUT=results/ab_royer/NW/grow_qc; HY=results/ab_royer/NW/hypo_qc
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
POOL=STW,B007,OPC,B04A,LRIV,B001,BLN,BS11,W020,B943,B013,W030,W040,DOSE,C04A,HDW,GNW,SQM,PL11,SMW,N050,N060,N080
( while true; do date -u +%FT%TZ > $R/NW_QC_HEARTBEAT; sleep 120; done ) 9>&- &
HB=$!
trap 'kill $HB 2>/dev/null' EXIT
fail(){ echo "stage=$1 $(date -u +%FT%TZ)" > $R/NW_QC_FAILED; exit 1; }
rm -f $R/NW_QC_FAILED
python scripts/ab_grow_fast.py --sub NW --pool "$POOL" --out "$OUT" --max-iters 14 --tol 0.01 --prep-workers 6 --restack-workers 8 || fail grow
LAST=$(ls $OUT/det_it*.parquet | sed 's/.*det_it//; s/\.parquet//' | sort -n | tail -1)
python scripts/ab_hypo.py --sub NW --pool "$POOL" --det "$OUT/det_it${LAST}.parquet" --out "$HY" --boot 60 || true
[ -f "$HY/loc_headline.csv" ] || fail hypo
date -u +%FT%TZ > $R/NW_QC_DONE
