#!/usr/bin/env bash
# SH_V2 (2026-09-28): Shelly (2025) Mendocino reproduction, detector v2 -- see results/SH_V2_LOG.md.
# prep (mute/lines/dead floors) -> it0 templates at his published times (own build days) ->
# for k=0..3: MF on held-out TEST days + score (source & fam) + report; MF on BUILD days -> restack it(k+1).
# Outputs only under results/sh_v2/.  Markers results/SH_V2_{DONE,FAILED,HEARTBEAT}.
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south
cd "$WS" || exit 1
mkdir -p results/sh_v2 logs locks
exec 9>locks/sh_v2.lock
flock -n 9 || { echo "[v2] another instance holds the lock; exiting"; exit 0; }
STAGE="init"; FINISHED=0; HB=""
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ)" > results/SH_V2_HEARTBEAT; }
hb_start(){ ( exec 9>&-; while sleep 300; do [ -f results/SH_V2_HEARTBEAT ] && touch results/SH_V2_HEARTBEAT; done ) & HB=$!; }
cleanup(){ [ -n "$HB" ] && kill "$HB" 2>/dev/null; }
fail(){ echo "[v2] FAILED stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > results/SH_V2_FAILED; cleanup; exit "$1"; }
onexit(){ rc=$?; cleanup; [ "$FINISHED" = 1 ] && return 0; [ -f results/SH_V2_FAILED ] && return 0
          echo "[v2] EXITED unexpectedly stage=$STAGE rc=$rc $(date -u +%FT%TZ)"
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > results/SH_V2_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[v2] === $STAGE $(date -u +%FT%TZ)"; "$@" || fail $?; }
soft(){ STAGE="$1"; shift; beat; echo "[v2] === $STAGE (non-fatal) $(date -u +%FT%TZ)"; "$@" || echo "[v2] WARN $STAGE rc=$? -- continuing"; }
rm -f results/SH_V2_FAILED results/SH_V2_DONE
hb_start
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
V=results/sh_v2
TEST=results/aa_shelly/stage2_days_test.json
BUILD=results/aa_shelly/stage2_days_build.json
T0=2018-01-01; T1=2024-12-01
SW=6,7,8,9,10,11,12

famfilter(){ python - "$1" "$2" <<'EOF'
import sys, json, pandas as pd
src, out = sys.argv[1], sys.argv[2]
D = pd.read_parquet(src); D["fam"] = D.fam.astype(str)
dd = json.load(open("results/aa_shelly/stage2_days_test.json"))
FD = {k.split("|")[1]: {pd.Timestamp(x) for x in v} for k, v in dd.items() if k.startswith("Shelly2025|")}
day = D.t.dt.floor("D")
D[[f in FD and d in FD[f] for f, d in zip(D.fam, day)]].to_parquet(out)
EOF
}

run prep python scripts/sh_v2_prep.py --workers 12 --t0 $T0 --t1 $T1
[ -f $V/tpl_it0.csv ] || run tpl0 python scripts/sh_v2_restack.py --base results/aa_shelly/templates6_build \
    --times shelly --days-file $BUILD --out $V/tpl_it0.npz
for k in 0 1 2 3; do
  [ -f $V/mf_it${k}_test.parquet ] || run mf${k}_test python scripts/t11_mf_v2.py --templates $V/tpl_it$k \
      --days-file $TEST --truth shelly --t0 $T0 --t1 $T1 --out $V/mf_it${k}_test.parquet
  [ -f $V/mf_it${k}_test_score.csv ] || run score${k} python scripts/t8_score_aa.py --ours $V/mf_it${k}_test.parquet \
      --truth shelly --t0 $T0 --t1 $T1 --days-file $TEST --days-mode source --sweep $SW --out $V/mf_it${k}_test_score.csv
  [ -f $V/mf_it${k}_test_fam.parquet ] || soft famfilt${k} famfilter $V/mf_it${k}_test.parquet $V/mf_it${k}_test_fam.parquet
  [ -f $V/mf_it${k}_test_score_fam.csv ] || soft score${k}_fam python scripts/t8_score_aa.py --ours $V/mf_it${k}_test_fam.parquet \
      --truth shelly --t0 $T0 --t1 $T1 --days-file $TEST --days-mode fam --sweep $SW --out $V/mf_it${k}_test_score_fam.csv
  soft report$k python scripts/sh_v2_report.py
  if [ $k -lt 3 ]; then
    n=$((k+1))
    [ -f $V/mf_it${k}_build.parquet ] || run mf${k}_build python scripts/t11_mf_v2.py --templates $V/tpl_it$k \
        --days-file $BUILD --t0 $T0 --t1 $T1 --out $V/mf_it${k}_build.parquet
    [ -f $V/tpl_it$n.csv ] || run tpl$n python scripts/sh_v2_restack.py --base $V/tpl_it$k \
        --det $V/mf_it${k}_build.parquet --days-file $BUILD --nstack 200 --thr 9 --out $V/tpl_it$n.npz
  fi
done
run report python scripts/sh_v2_report.py
date -u +%FT%TZ > results/SH_V2_DONE
STAGE="done"; FINISHED=1; beat
echo "[v2] DONE $(date -u +%FT%TZ)"
