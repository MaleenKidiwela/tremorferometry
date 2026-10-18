#!/usr/bin/env bash
# Detached queue: corrected times for all 100 kept families, then the LITERAL network
# autocorrelation over the whole 24-day ETS bed, then the pair-stage score.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/repro results/markers
exec 9>locks/k_repro.lock
flock -n 9 || { echo "already running"; exit 0; }
rm -f results/K_REPRO1_DONE results/K_REPRO1_FAILED
trap 'echo "[queue] TRAP exit $?"; [ -f results/K_REPRO1_DONE ] || echo "stage=$STAGE" > results/K_REPRO1_FAILED' EXIT
STAGE=start
echo "[queue] start $(date -u +%FT%TZ) pid=$$"

STAGE=s1_all100
if [ ! -s results/repro/s1all_times.csv ]; then
  python scripts/k_s1_loud.py --nfam 100 --out results/repro/s1all || exit 1
fi
echo "[queue] s1all done $(date -u +%FT%TZ)"

STAGE=s2a_all100
python scripts/k_s2a_stacks.py --times results/repro/s1all_times.csv \
  --stas GLBC,JRBC,LCBC,MGCB,SHVB,TWBB,LZB,PGC,SNB,VGZ,B001,B003,B004,B005,B006,B007,SILB,SSIB,TSJB,TWKB \
  --out results/repro/s2a_all100 || exit 1
echo "[queue] s2a done $(date -u +%FT%TZ)"

STAGE=autocorr_full
DATES=$(python - <<'PY'
import pandas as pd
print(",".join(d.strftime("%Y-%m-%d") for d in pd.date_range("2005-08-29","2005-09-21")))
PY
)
python scripts/k_auto.py --dates "$DATES" --hours all --out results/repro/auto_full || exit 1
echo "[queue] autocorr done $(date -u +%FT%TZ)"

STAGE=score
python scripts/k_score_pairs.py --det results/repro/auto_full \
  --times results/repro/s1all_times.csv > logs/k_score_full.log 2>&1 || exit 1
tail -40 logs/k_score_full.log

STAGE=auto7
python scripts/k_auto.py --dates 2005-09-13,2005-09-14 --hours all --nsta 7 \
  --out results/repro/auto7 || true
python scripts/k_auto.py --dates 2005-09-13,2005-09-14 --hours all \
  --out results/repro/auto11_2d || true
python scripts/k_score_pairs.py --det results/repro/auto7 \
  --times results/repro/s1all_times.csv > logs/k_score_auto7.log 2>&1 || true
python scripts/k_score_pairs.py --det results/repro/auto11_2d \
  --times results/repro/s1all_times.csv > logs/k_score_auto11_2d.log 2>&1 || true

echo "OK $(date -u +%FT%TZ)" > results/K_REPRO1_DONE
echo "[queue] DONE $(date -u +%FT%TZ)"
