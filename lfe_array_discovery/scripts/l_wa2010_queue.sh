#!/usr/bin/env bash
# L_WA2010 -- Bostock's published recipe (agent K's literal engine, unedited) on the Washington
# N-Olympic 2010 array bed, as ONE detached queue.
#
#   nohup setsid -f bash scripts/l_wa2010_queue.sh >> logs/l_wa2010_queue.log 2>&1 < /dev/null
#
# Every stage pins its OWN environment with `env VAR=... python` and nothing is exported across
# stages -- the K_YEAR leak that voided the first blind run (ledger sec.9) cannot happen here.
# Every stage has an output guard, so the queue resumes from where it stopped.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/l_wa2010 results/l_wa2010_null

exec 9>locks/l_wa2010.lock
flock -n 9 || { echo "[L] already running"; exit 0; }
rm -f results/L_WA2010_DONE results/L_WA2010_FAILED
STAGE=start
fail() { [ -f results/L_WA2010_DONE ] || \
  printf 'stage=%s signal=%s at=%s pid=%s\n' "$STAGE" "${1:-exit}" "$(date -u +%FT%TZ)" "$$" \
  > results/L_WA2010_FAILED; }
trap 'rc=$?; echo "[L] TRAP EXIT rc=$rc stage=$STAGE"; fail "exit$rc"' EXIT
trap 'echo "[L] TRAP TERM stage=$STAGE"; fail TERM; exit 143' TERM
trap 'echo "[L] TRAP INT  stage=$STAGE"; fail INT;  exit 130' INT
trap 'echo "[L] TRAP HUP  stage=$STAGE"; fail HUP;  exit 129' HUP

# CUDA-12 runtime libs (memory: cupy/torch die with libcufft.so.11 missing after a pod restart)
CU=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')
export LD_LIBRARY_PATH="${CU}${LD_LIBRARY_PATH:-}"

R=results/l_wa2010
N=results/l_wa2010_null
GROW_MAX=1500          # the north grew 1,514 groups; the GPU is shared with k_repro5
CHUNK=750              # k_grow allocates one day-length accumulator row per template
echo "[L] start $(date -u +%FT%TZ) pid=$$ ppid=$PPID"

# ---------------------------------------------------------------- bed + hours ------------------
STAGE=bed
if [ ! -s "$R/bed.json" ]; then
  env L_OUT="$R" python - <<'PY' > logs/l_bed.log 2>&1 || exit 1
import sys, json, os; sys.path.insert(0, "scripts")
import l_common as L
os.makedirs("results/l_wa2010", exist_ok=True); os.makedirs("results/l_wa2010_null", exist_ok=True)
b = L.bed_summary()
H = L.candidate_hours()
b["n_candidate_hours"] = int(len(H)); b["n_tremor_days"] = int(H.date.nunique())
b["n_quiet_hours"] = int(len(L.quiet_hours()))
corr, flagged = L.station_corrections()
b["n_intraslab_corr"] = int(len(corr) - len(flagged)); b["n_corr_flagged"] = int(len(flagged))
json.dump(b, open("results/l_wa2010/bed.json", "w"), indent=1)
json.dump(L.null_offsets(L.stations()), open("results/l_wa2010_null/null_offsets.json", "w"),
          indent=1)
print(json.dumps({k: v for k, v in b.items() if k != "stations"}, indent=1))
open("results/l_wa2010/tremor_dates.txt", "w").write(",".join(sorted(set(H.date))))
PY
fi
tail -30 logs/l_bed.log
[ -s "$R/tremor_dates.txt" ] || { echo "[L] no tremor_dates.txt"; exit 1; }
TREMOR_DATES=$(cat "$R/tremor_dates.txt")
echo "[L] tremor days: $(echo "$TREMOR_DATES" | tr ',' '\n' | wc -l)"

# ---------------------------------------------------------------- one chain --------------------
# $1 = output dir, $2 = 1 for the site-shift null, $3 = tag
chain() {
  local OUT="$1" NUL="$2" TAG="$3"
  local E="L_OUT=$OUT"
  if [ "$NUL" = "1" ]; then E="$E L_NULL_SHIFT=1 L_NULL_SEED=11"; fi
  mkdir -p "$OUT"

  STAGE="${TAG}_auto_tremor"
  if [ ! -s "$OUT/auto/_params.json" ]; then
    env $E python scripts/l_auto.py --hours tremor --out "$OUT/auto" \
      > "logs/l_${TAG}_auto.log" 2>&1 || exit 1
  fi
  tail -3 "logs/l_${TAG}_auto.log"

  if [ "$NUL" != "1" ]; then
    STAGE="${TAG}_auto_quiet"
    if [ ! -s "$OUT/auto_quiet/_params.json" ]; then
      env $E python scripts/l_auto.py --hours quiet --out "$OUT/auto_quiet" \
        > "logs/l_${TAG}_auto_quiet.log" 2>&1 || exit 1
    fi
    tail -3 "logs/l_${TAG}_auto_quiet.log"
  fi

  STAGE="${TAG}_pairstats"
  if [ ! -s "$OUT/pair_pairstats.json" ]; then
    if [ "$NUL" = "1" ]; then
      env $E python scripts/l_pairstats.py --tremor "$OUT/auto" --out "$OUT/pair" \
        > "logs/l_${TAG}_pairstats.log" 2>&1 || exit 1
    else
      env $E python scripts/l_pairstats.py --tremor "$OUT/auto" --quiet "$OUT/auto_quiet" \
        --out "$OUT/pair" > "logs/l_${TAG}_pairstats.log" 2>&1 || exit 1
    fi
  fi
  tail -40 "logs/l_${TAG}_pairstats.log"

  STAGE="${TAG}_group"
  if [ ! -s "$OUT/gb_summary.csv" ]; then
    env $E python scripts/l_run.py k_group --detfile "$OUT/pair_all.parquet" \
      --budget 100000000 --min-members 3 --times /nonexistent --out "$OUT/gb" \
      > "logs/l_${TAG}_group.log" 2>&1 || exit 1
  fi
  grep -E "nodes,|clean detections" "logs/l_${TAG}_group.log" || true

  STAGE="${TAG}_split"
  if [ ! -s "$OUT/chunks/_split.json" ]; then
    env $E python scripts/l_chunk.py split --nodes "$OUT/gb_nodes.csv" \
      --summary "$OUT/gb_summary.csv" --out "$OUT/chunks" --size "$CHUNK" \
      --max-groups "$GROW_MAX" > "logs/l_${TAG}_split.log" 2>&1 || exit 1
  fi
  head -1 "logs/l_${TAG}_split.log"
  local NC
  NC=$(python -c "import json;print(json.load(open('$OUT/chunks/_split.json'))['n_chunks'])")
  case "$NC" in (''|*[!0-9]*) echo "[L] bad chunk count '$NC'"; exit 1;; esac

  local i CH DIRS=""
  for ((i=0; i<NC; i++)); do
    CH=$(printf "%02d" "$i")
    STAGE="${TAG}_grow_c${CH}"
    if [ ! -s "$OUT/grow_c${CH}/det_it0.parquet" ]; then
      env $E python scripts/l_run.py k_grow --nodes "$OUT/chunks/chunk${CH}_nodes.csv" \
        --summary "$OUT/chunks/chunk${CH}_summary.csv" --dates "$TREMOR_DATES" \
        --iters 0 --out "$OUT/grow_c${CH}" > "logs/l_${TAG}_grow_c${CH}.log" 2>&1 || exit 1
    fi
    tail -2 "logs/l_${TAG}_grow_c${CH}.log"
    DIRS="${DIRS}${OUT}/grow_c${CH},"
  done

  STAGE="${TAG}_merge"
  if [ ! -s "$OUT/grow/det_it0.parquet" ]; then
    env $E python scripts/l_chunk.py merge --dirs "$DIRS" --out "$OUT/grow" \
      > "logs/l_${TAG}_merge.log" 2>&1 || exit 1
  fi
  cat "logs/l_${TAG}_merge.log"

  STAGE="${TAG}_locate"
  if [ ! -s "$OUT/loc_it0.csv" ]; then
    env $E python scripts/l_locate.py --tpl "$OUT/grow/tpl_it0.npz" \
      --det "$OUT/grow/det_it0.parquet" --out "$OUT/loc_it0.csv" \
      > "logs/l_${TAG}_locate.log" 2>&1 || exit 1
  fi
  tail -3 "logs/l_${TAG}_locate.log"

  STAGE="${TAG}_verify"
  if [ ! -s "$OUT/verify_it0.csv" ]; then
    env $E python scripts/l_verify.py --loc "$OUT/loc_it0.csv" \
      --det "$OUT/grow/det_it0.parquet" --out "$OUT/verify_it0.csv" --tag "$TAG" \
      > "logs/l_${TAG}_verify.log" 2>&1 || exit 1
  fi
  tail -30 "logs/l_${TAG}_verify.log"
}

echo "[L] ===================== REAL CHAIN ====================="
chain "$R" 0 real
echo "[L] ===================== SITE-SHIFT NULL ================"
chain "$N" 1 null

STAGE=report
env L_OUT="$R" python scripts/l_report.py --real "$R" --null "$N" \
  --out results/WA_2010_RESULTS.md > logs/l_report.log 2>&1 || exit 1
cat logs/l_report.log

echo "OK $(date -u +%FT%TZ)" > results/L_WA2010_DONE
echo "[L] DONE $(date -u +%FT%TZ)"
