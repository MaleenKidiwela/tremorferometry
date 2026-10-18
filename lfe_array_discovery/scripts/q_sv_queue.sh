#!/usr/bin/env bash
# Q_SV: the two missing measurements for the screen-validation table on the 2005 POLARIS bed.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/q_screenval
exec 9>locks/q_sv.lock
flock -n 9 || { echo "already running"; exit 0; }
rm -f results/Q_SV_DONE results/Q_SV_FAILED
STAGE=start
trap '[ -f results/Q_SV_DONE ] || echo "stage=$STAGE" > results/Q_SV_FAILED' EXIT
E2005="K_YEAR=2005 K_DAY0=2005-08-29 K_DAY1=2005-09-21"
echo "[qsv] start $(date -u +%FT%TZ) pid=$$"

STAGE=quiet
env $E2005 python scripts/q_sv_data.py --stage quiet > logs/q_sv_quiet.log 2>&1 || exit 1
tail -3 logs/q_sv_quiet.log

STAGE=restack_h2
env $E2005 python scripts/q_sv_data.py --stage restack --half 2 > logs/q_sv_h2.log 2>&1 || exit 1
tail -3 logs/q_sv_h2.log

STAGE=restack_h1
env $E2005 python scripts/q_sv_data.py --stage restack --half 1 > logs/q_sv_h1.log 2>&1 || exit 1
tail -3 logs/q_sv_h1.log

echo "OK $(date -u +%FT%TZ)" > results/Q_SV_DONE
echo "[qsv] DONE $(date -u +%FT%TZ)"
