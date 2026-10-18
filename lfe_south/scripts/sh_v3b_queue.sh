#!/usr/bin/env bash
# SH_V3B (2026-09-29): v3 (scripts/sh_v3_queue.sh) rerun EXACTLY, with the full 17-station Text S1 set:
# the 13 v3 stations + BK.GHOP, NC.KMPB, PB.B932, BK.HUNT (see results/SH_V3B_LOG.md).
# data (guarded): plan -> t2_fetch (NCEDC for BK/NC, IRIS for PB; unchanged) -> t3_cache (unchanged, 20 sps 1-8 Hz)
# -> prep (sh_v3b_data.py: v2 prep entries kept, new stations added) -> the v3 chain with SH_V3_DIR / SH_PREP_DIR
# pointing at results/sh_v3b.  Markers results/SH_V3B_{DONE,FAILED,HEARTBEAT}; pid results/SH_V3B_PID.
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south
cd "$WS" || exit 1
mkdir -p results/sh_v3b logs locks
exec 9>locks/sh_v3b.lock
flock -n 9 || { echo "[v3b] another instance holds the lock; exiting"; exit 0; }
echo $$ > results/SH_V3B_PID
STAGE="init"; FINISHED=0; HB=""
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ)" > results/SH_V3B_HEARTBEAT; }
hb_start(){ ( exec 9>&-; while sleep 300; do [ -f results/SH_V3B_HEARTBEAT ] && touch results/SH_V3B_HEARTBEAT; done ) & HB=$!; }
cleanup(){ [ -n "$HB" ] && kill "$HB" 2>/dev/null; }
fail(){ echo "[v3b] FAILED stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > results/SH_V3B_FAILED; cleanup; exit "$1"; }
onexit(){ rc=$?; cleanup; [ "$FINISHED" = 1 ] && return 0; [ -f results/SH_V3B_FAILED ] && return 0
          echo "[v3b] EXITED unexpectedly stage=$STAGE rc=$rc $(date -u +%FT%TZ)"
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > results/SH_V3B_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[v3b] === $STAGE $(date -u +%FT%TZ)"; "$@" || fail $?; }
soft(){ STAGE="$1"; shift; beat; echo "[v3b] === $STAGE (non-fatal) $(date -u +%FT%TZ)"; "$@" || echo "[v3b] WARN $STAGE rc=$? -- continuing"; }
rm -f results/SH_V3B_FAILED results/SH_V3B_DONE
hb_start
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
V=results/sh_v3b
export SH_V3_DIR=$WS/$V SH_PREP_DIR=$WS/$V/prep SH_V3_LOG=$WS/results/SH_V3B_LOG.md SH_V3_BASE=$WS/results/sh_v3 SH_V3_NAME=SH_V3B
TEST=results/aa_shelly/stage2_days_test.json
BUILD=results/aa_shelly/stage2_days_build.json
T0=2018-01-01; T1=2024-12-01

datasum(){ python - <<'PY'
import pandas as pd, os
P = pd.read_csv("results/sh_v3b/download_plan.csv")
M = pd.read_csv("results/sh_v3b/download_manifest.csv").drop_duplicates(["net","sta","day"], keep="last")
C = pd.read_parquet("results/sh_v3b/cache_v3b.parquet")
S = pd.read_csv("results/sh_v3b/prep/../prep_summary_new.csv") if os.path.exists("results/sh_v3b/prep_summary_new.csv") else None
L = ["", f"### data obtained (auto, {pd.Timestamp.now('UTC'):%Y-%m-%d %H:%MZ})", "",
     "| station | planned days | miniSEED days | no-data days | fetch errors | cached (built/cached) | prepped | station-dead days | band |", "|---|---|---|---|---|---|---|---|---|"]
for (n, s), g in P.groupby(["net", "sta"]):
    ns = f"{n}.{s}"; m = M[(M.net == n) & (M.sta == s)]; c = C[C.netsta == ns]
    sp = S[S.netsta == ns] if S is not None else pd.DataFrame()
    L.append(f"| {ns} | {len(g)} | {(m.status == 'data').sum()} | {(m.status == 'nodata').sum()} | {(m.status == 'error').sum()} | "
             f"{c.status.isin(['built', 'cached']).sum()} | {len(sp)} | {int(sp.station_dead.sum()) if len(sp) else 0} | "
             f"{','.join(sorted(c.band.dropna().unique())) if 'band' in c else ''} |")
s = open("results/SH_V3B_LOG.md").read()
B0 = "<!-- SH_V3_RESULTS_BEGIN -->"
s = (s[:s.index(B0)].rstrip() + "\n" + "\n".join(L) + "\n\n" + s[s.index(B0):]) if B0 in s else s.rstrip() + "\n" + "\n".join(L) + "\n"
open("results/SH_V3B_LOG.md", "w").write(s); print("\n".join(L))
PY
}

# ---- data (guarded) ------------------------------------------------------------------------------
[ -f $V/download_plan.csv ] || run plan python scripts/sh_v3b_data.py plan
if [ ! -f $V/FETCH_DONE ]; then
  run fetch python scripts/t2_fetch.py $V/download_plan.csv --workers 8 --manifest $V/download_manifest.csv
  run fetch_retry python scripts/t2_fetch.py $V/download_plan.csv --workers 4 --manifest $V/download_manifest.csv
  date -u +%FT%TZ > $V/FETCH_DONE
fi
if [ ! -f $V/CACHE_DONE ]; then
  run cache python scripts/t3_cache.py $V/download_plan.csv --workers 24 --out $V/cache_v3b.parquet
  date -u +%FT%TZ > $V/CACHE_DONE
fi
if [ ! -f $V/PREP_DONE ]; then
  run prep python scripts/sh_v3b_data.py prep --workers 12
  date -u +%FT%TZ > $V/PREP_DONE
fi
soft datasum datasum

# ---- v3 chain, unchanged --------------------------------------------------------------------------
[ -f $V/tpl_it0.csv ] || run tpl0 python scripts/sh_v3_restack.py --mode his --days-file $BUILD --workers 4 --out $V/tpl_it0.npz
for k in 0 1; do
  [ -f $V/mf_it${k}_test.parquet ] || run mf${k}_test python scripts/t11_mf_v2.py --templates $V/tpl_it$k \
      --days-file $TEST --truth shelly --t0 $T0 --t1 $T1 --out $V/mf_it${k}_test.parquet
  [ -f $V/rescore_it$k.csv ] || run rescore$k python scripts/sh_v3_rescore.py --ours $V/mf_it${k}_test.parquet --tag it$k --out $V/rescore_it$k.csv
  soft report$k python scripts/sh_v3_report.py
  if [ $k -eq 0 ]; then
    [ -f $V/mf_it0_build.parquet ] || run mf0_build python scripts/t11_mf_v2.py --templates $V/tpl_it0 \
        --days-file $BUILD --t0 $T0 --t1 $T1 --out $V/mf_it0_build.parquet
    [ -f $V/tpl_it1.csv ] || run tpl1 python scripts/sh_v3_restack.py --mode det --det $V/mf_it0_build.parquet \
        --picks $V/tpl_it0.csv --days-file $BUILD --thr 8 --workers 4 --out $V/tpl_it1.npz
  fi
done
run mix python scripts/sh_v3_mix.py
[ -f $V/rescore_mix.csv ] || run rescore_mix python scripts/sh_v3_rescore.py --ours $V/mf_mix_test.parquet --tag mix --out $V/rescore_mix.csv
run report python scripts/sh_v3_report.py
date -u +%FT%TZ > results/SH_V3B_DONE
STAGE="done"; FINISHED=1; beat
echo "[v3b] DONE $(date -u +%FT%TZ)"
