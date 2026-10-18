#!/bin/bash
# A4 location-seeded positive control queue.  Usage: a4_queue.sh <pilot|full>
# Detached launch:  nohup setsid -f bash scripts/a4_queue.sh full > logs/a4_full.log 2>&1
# Markers (full): results/A4_SEED_{DONE,FAILED,HEARTBEAT}; pilot: results/a4_seed/pilot/PILOT_{DONE,FAILED}
# Jobs: (NW,SW,NVI) x (real,nullT,nullR), at most $PAR (=4) at once, 3 threads each (<=12 CPUs).
set -u
MODE=${1:-full}
ST=/home/jovyan/tremorferometry/lfe_array_discovery
cd $ST || exit 1
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=3 MKL_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 PYTHONUNBUFFERED=1
PAR=${PAR:-4}
if [ "$MODE" = pilot ]; then
  OUT=results/a4_seed/pilot; LOCS=results/a4_seed/locs_pilot.csv
  DONE=$OUT/PILOT_DONE; FAILED=$OUT/PILOT_FAILED; HB=$OUT/PILOT_HEARTBEAT
else
  OUT=results/a4_seed/full; LOCS=results/a4_seed/locs_all.csv
  DONE=results/A4_SEED_DONE; FAILED=results/A4_SEED_FAILED; HB=results/A4_SEED_HEARTBEAT
fi
mkdir -p $OUT locks logs
exec 9>locks/a4_seed_$MODE.lock
flock -n 9 || { echo "another a4 $MODE queue holds the lock"; exit 1; }
echo $$ > $OUT/QUEUE_PID
rm -f $DONE $FAILED
OK=0
trap 'if [ $OK -ne 1 ]; then date -u +%FT%TZ" queue exited abnormally (see logs/a4_${MODE}_*.log)" > $FAILED; fi; kill $HBPID 2>/dev/null' EXIT
( while true; do date -u +%FT%TZ > $HB; sleep 60; done ) & HBPID=$!

FAILS=0
run() {  # sub cond
  python scripts/a4_seed.py --sub $1 --locs $LOCS --cond $2 --out $OUT > logs/a4_${MODE}_$1_$2.log 2>&1
  echo "$1 $2 exit $?" >> $OUT/job_status.txt
}
: > $OUT/job_status.txt
# biggest (NW) first so the long jobs overlap
JOBS="NW:real NW:nullT NW:nullR SW:real SW:nullT SW:nullR NVI:real NVI:nullT NVI:nullR"
for j in $JOBS; do
  while [ $(jobs -rp | wc -l) -ge $((PAR + 1)) ]; do sleep 20; done   # +1 = the heartbeat loop
  run ${j%%:*} ${j##*:} &
done
while [ $(jobs -rp | wc -l) -gt 1 ]; do sleep 20; done
python scripts/a4_report.py $OUT --title "A4 location-seeded positive control ($MODE)" > logs/a4_${MODE}_report.log 2>&1
REP=$?
if grep -qv "exit 0" $OUT/job_status.txt || [ $REP -ne 0 ]; then
  date -u +%FT%TZ" job failures: $(grep -v 'exit 0' $OUT/job_status.txt | tr '\n' ';') report exit $REP" > $FAILED
  OK=1; exit 1
fi
[ "$MODE" = full ] && cp $OUT/REPORT.md results/a4_seed/REPORT.md
date -u +%FT%TZ > $DONE
OK=1
