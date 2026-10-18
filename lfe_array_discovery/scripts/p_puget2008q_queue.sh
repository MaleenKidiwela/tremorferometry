#!/usr/bin/env bash
# P_PUGET2008Q -- station QUALITY as the lever on the Puget / S-Olympic 2008 bed.
#
# Agent N's 51-station CAFE+PB run (results/WA_2008_RESULTS.md) returned a pair-stage
# tremor/quiet ratio of 0.89 and grown templates that fire equally on quiet and tremor days: the
# pair stage is dominated by repeating cultural/site noise, not LFE pairs.  Royer & Bostock
# (2014) built COMPOSITE ARRAYS out of the stations that worked.  This queue
#   1. measures, per station, who manufactures those pairs (p_census.py / p_subset.py),
#   2. re-runs THE UNCHANGED RECIPE on the quiet subset the census defines,
#   3. runs the same five-clause screen real vs a whole-chain site-shift null,
#   4. re-fires on 2007-01,
# reusing agent N's audited tremor hours, quiet days, health veto, screen and confirmation rule.
# The ONLY change to any n_/k_ file is the optional `N_STATION_SUBSET` whitelist in
# n_common.admitted(): unset, every agent-N number is reproduced byte-for-byte.
#
#   nohup setsid -f bash scripts/p_puget2008q_queue.sh >> logs/p_puget2008q_queue.log 2>&1 < /dev/null
#
# Every stage pins its OWN environment with `env VAR=... python`; nothing is exported across
# stages.  Every stage has an output guard, so the queue resumes from where it stopped.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/p_2008

exec 9>locks/p_puget2008q.lock
flock -n 9 || { echo "[P] already running"; exit 0; }
rm -f results/P_PUGET2008Q_DONE results/P_PUGET2008Q_FAILED
STAGE=start
fail() { [ -f results/P_PUGET2008Q_DONE ] || \
  printf 'stage=%s signal=%s at=%s pid=%s\n' "$STAGE" "${1:-exit}" "$(date -u +%FT%TZ)" "$$" \
  > results/P_PUGET2008Q_FAILED; }
trap 'rc=$?; echo "[P] TRAP EXIT rc=$rc stage=$STAGE"; fail "exit$rc"' EXIT
trap 'echo "[P] TRAP TERM stage=$STAGE"; fail TERM; exit 143' TERM
trap 'echo "[P] TRAP INT  stage=$STAGE"; fail INT;  exit 130' INT
trap 'echo "[P] TRAP HUP  stage=$STAGE"; fail HUP;  exit 129' HUP

# CUDA-12 runtime libs (memory: torch dies with libcufft.so.11 missing after a pod restart)
CU=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')
export LD_LIBRARY_PATH="${CU}${LD_LIBRARY_PATH:-}"

N=results/n_puget2008          # agent N's outputs: hours, quiet days, 51-station pair stage
P=results/p_2008               # this run
R=$P/2008
NUL=$P/2008_null
BHD=$P/bh
C7=$P/2007
SUB=$P/subset_stations.txt
BHL=$P/boreholes_stations.txt
GROW_MAX=1500
CHUNK=750
MAX_SCREEN_PASS=200
mkdir -p "$R" "$NUL" "$C7" "$BHD"
echo "[P] start $(date -u +%FT%TZ) pid=$$ ppid=$PPID"

for f in "$N/2008/auto/_stats.csv" "$N/2008/auto_quiet/_stats.csv" "$N/tremor_dates_2008.txt" \
         "$N/quiet_dates_2008.txt" "$N/hour_audit_2007.csv"; do
  [ -s "$f" ] || { echo "[P] missing agent-N input $f"; exit 1; }
done
TREMOR_DATES=$(cat "$N/tremor_dates_2008.txt")
QUIET_DATES=$(cat "$N/quiet_dates_2008.txt")

# ================================================================ 1. the census ================
for H in tremor quiet; do
  STAGE="census_$H"
  A="$N/2008/auto"; [ "$H" = "quiet" ] && A="$N/2008/auto_quiet"
  if [ ! -s "$P/census_${H}.parquet" ]; then
    env N_BED=2008 python scripts/p_census.py --hours "$H" --auto "$A" \
      --out "$P/census_${H}.parquet" > "logs/p_census_${H}.log" 2>&1 || exit 1
  fi
  tail -2 "logs/p_census_${H}.log"
done

STAGE=subset
if [ ! -s "$P/subset_meta.json" ]; then
  env N_BED=2008 python scripts/p_subset.py --tremor "$P/census_tremor.parquet" \
    --quiet "$P/census_quiet.parquet" --outdir "$P" > logs/p_subset.log 2>&1 || exit 1
fi
cat logs/p_subset.log
[ -s "$SUB" ] || { echo "[P] no subset station list"; exit 1; }
echo "[P] subset: $(cat "$SUB")"

# ================================================================ 2. bed summaries =============
STAGE=bed
if [ ! -s "$P/bed_2008q.json" ] || [ ! -s "$P/bed_2007q.json" ]; then
  for B in 2008 2007; do
    env N_BED="$B" N_STATION_SUBSET="$SUB" P_DIR="$P" python - <<'PY' >> logs/p_bed.log 2>&1 || exit 1
import sys, json, os
sys.path.insert(0, "scripts")
import n_common as N
P = os.environ["P_DIR"]
b = N.bed_summary()
corr, flagged = N.station_corrections()
b["n_intraslab_corr"] = int(len(corr) - len(flagged))
b["n_corr_flagged"] = int(len(flagged))
H = N.candidate_hours()
b["n_tremor_hours"] = int(len(H)); b["n_tremor_days"] = int(H.date.nunique())
json.dump(b, open(f"{P}/bed_{N.BED}q.json", "w"), indent=1)
if N.BED == "2008":
    os.makedirs(f"{P}/2008_null", exist_ok=True)
    json.dump(N.null_offsets(N.stations()),
              open(f"{P}/2008_null/null_offsets.json", "w"), indent=1)
print(json.dumps({k: v for k, v in b.items() if k != "stations"}, indent=1))
PY
  done
fi
tail -40 logs/p_bed.log

# ================================================================ 3. one discovery chain =======
# $1 = output dir, $2 = 1 for the site-shift null, $3 = tag, $4 = station list
chain() {
  local OUT="$1" NU="$2" TAG="$3" STL="$4"
  local E="N_BED=2008 N_OUT=$OUT N_STATION_SUBSET=$STL"
  if [ "$NU" = "1" ]; then E="$E N_NULL_SHIFT=1 N_NULL_SEED=11"; fi
  mkdir -p "$OUT"

  STAGE="${TAG}_auto_tremor"
  if [ ! -s "$OUT/auto/_params.json" ]; then
    env $E python scripts/n_auto.py --hours tremor --out "$OUT/auto" \
      > "logs/p_${TAG}_auto.log" 2>&1 || exit 1
  fi
  tail -3 "logs/p_${TAG}_auto.log"

  if [ "$NU" != "1" ]; then
    STAGE="${TAG}_auto_quiet"
    if [ ! -s "$OUT/auto_quiet/_params.json" ]; then
      env $E python scripts/n_auto.py --hours quiet --out "$OUT/auto_quiet" \
        > "logs/p_${TAG}_auto_quiet.log" 2>&1 || exit 1
    fi
    tail -3 "logs/p_${TAG}_auto_quiet.log"
  fi

  STAGE="${TAG}_pairstats"
  if [ ! -s "$OUT/pair_pairstats.json" ]; then
    if [ "$NU" = "1" ]; then
      env $E python scripts/n_pairstats.py --tremor "$OUT/auto" --out "$OUT/pair" \
        > "logs/p_${TAG}_pairstats.log" 2>&1 || exit 1
    else
      env $E python scripts/n_pairstats.py --tremor "$OUT/auto" --quiet "$OUT/auto_quiet" \
        --out "$OUT/pair" > "logs/p_${TAG}_pairstats.log" 2>&1 || exit 1
    fi
  fi
  tail -40 "logs/p_${TAG}_pairstats.log"
  [ "${5:-}" = "paironly" ] && return 0

  STAGE="${TAG}_group"
  if [ ! -s "$OUT/gb_summary.csv" ]; then
    env $E python scripts/n_run.py k_group --detfile "$OUT/pair_all.parquet" \
      --budget 100000000 --min-members 3 --times /nonexistent --out "$OUT/gb" \
      > "logs/p_${TAG}_group.log" 2>&1 || exit 1
  fi
  grep -E "nodes,|clean detections" "logs/p_${TAG}_group.log" || true

  STAGE="${TAG}_split"
  if [ ! -s "$OUT/chunks/_split.json" ]; then
    env $E python scripts/n_chunk.py split --nodes "$OUT/gb_nodes.csv" \
      --summary "$OUT/gb_summary.csv" --out "$OUT/chunks" --size "$CHUNK" \
      --max-groups "$GROW_MAX" > "logs/p_${TAG}_split.log" 2>&1 || exit 1
  fi
  head -1 "logs/p_${TAG}_split.log"
  local NC
  NC=$(python -c "import json;print(json.load(open('$OUT/chunks/_split.json'))['n_chunks'])")
  case "$NC" in (''|*[!0-9]*) echo "[P] bad chunk count '$NC'"; exit 1;; esac
  if [ "$NC" -eq 0 ]; then echo "[P] $TAG: no groups to grow"; return 0; fi

  local i CH DIRS=""
  for ((i=0; i<NC; i++)); do
    CH=$(printf "%02d" "$i")
    STAGE="${TAG}_grow_c${CH}"
    if [ ! -s "$OUT/grow_c${CH}/det_it0.parquet" ]; then
      env $E python scripts/n_run.py k_grow --nodes "$OUT/chunks/chunk${CH}_nodes.csv" \
        --summary "$OUT/chunks/chunk${CH}_summary.csv" --dates "$TREMOR_DATES" \
        --iters 0 --out "$OUT/grow_c${CH}" > "logs/p_${TAG}_grow_c${CH}.log" 2>&1 || exit 1
    fi
    tail -2 "logs/p_${TAG}_grow_c${CH}.log"
    DIRS="${DIRS}${OUT}/grow_c${CH},"
  done

  STAGE="${TAG}_merge"
  if [ ! -s "$OUT/grow/det_it0.parquet" ]; then
    env $E python scripts/n_chunk.py merge --dirs "$DIRS" --out "$OUT/grow" \
      > "logs/p_${TAG}_merge.log" 2>&1 || exit 1
  fi
  cat "logs/p_${TAG}_merge.log"

  STAGE="${TAG}_quiet_fa"
  if [ ! -s "$OUT/det_fa.parquet" ]; then
    env $E python scripts/n_refire.py --tpl "$OUT/grow/tpl_it0.npz" --dates "$QUIET_DATES" \
      --out "$OUT" --tag fa > "logs/p_${TAG}_quietfa.log" 2>&1 || exit 1
  fi
  tail -12 "logs/p_${TAG}_quietfa.log"

  STAGE="${TAG}_locate"
  if [ ! -s "$OUT/loc_it0.csv" ]; then
    env $E python scripts/n_locate.py --tpl "$OUT/grow/tpl_it0.npz" \
      --det "$OUT/grow/det_it0.parquet" --out "$OUT/loc_it0.csv" \
      > "logs/p_${TAG}_locate.log" 2>&1 || exit 1
  fi
  tail -3 "logs/p_${TAG}_locate.log"

  STAGE="${TAG}_screen"
  if [ ! -s "$OUT/screen.csv" ]; then
    env $E python scripts/n_screen.py --stage screen --loc "$OUT/loc_it0.csv" \
      --det "$OUT/grow/det_it0.parquet" --det-quiet "$OUT/det_fa.parquet" \
      --out "$OUT/screen.csv" --tag "$TAG" > "logs/p_${TAG}_screen.log" 2>&1 || exit 1
  fi
  tail -60 "logs/p_${TAG}_screen.log"
}

echo "[P] ===================== QUIET-SUBSET REAL CHAIN ========"
chain "$R" 0 real "$SUB"
echo "[P] ===================== SITE-SHIFT NULL ================"
chain "$NUL" 1 null "$SUB"

# --- boreholes-only bed: pair stage only.  Its job is to answer ONE question -- whether ANY
# --- Washington station set has an LFE-dominated pair stage -- so it stops at the ratio.
echo "[P] ===================== BOREHOLES-ONLY PAIR STAGE ======"
chain "$BHD" 0 bh "$BHL" paironly

# ================================================================ 4. clause (e): stability =====
stability() {
  local OUT="$1" NU="$2" TAG="$3" STL="$4"
  local E="N_BED=2008 N_OUT=$OUT N_STATION_SUBSET=$STL"
  if [ "$NU" = "1" ]; then E="$E N_NULL_SHIFT=1 N_NULL_SEED=11"; fi
  [ -s "$OUT/screen.csv" ] || { echo "[P] $TAG: no screen table"; return 0; }

  STAGE="${TAG}_passlist"
  if [ ! -s "$OUT/screen_pass.csv" ]; then
    env $E python - "$OUT" "$MAX_SCREEN_PASS" <<'PY' || exit 1
import sys, pandas as pd
out, cap = sys.argv[1], int(sys.argv[2])
V = pd.read_csv(f"{out}/screen.csv")
P = V[V.pass_abcd.astype(bool)].sort_values("n_det", ascending=False).head(cap)
P[["grp"]].to_csv(f"{out}/screen_pass.csv", index=False)
print(f"{len(P)} templates pass clauses (a)-(d) of {len(V)}", flush=True)
PY
  fi
  local NP
  NP=$(python -c "import pandas as pd;print(len(pd.read_csv('$OUT/screen_pass.csv')))")
  echo "[P] $TAG: $NP templates carried into the stability re-stack"

  if [ "$NP" -gt 0 ]; then
    for H in 2 1 0; do
      STAGE="${TAG}_reloc_h${H}"
      local NAME="reloc_h${H}"; [ "$H" = "0" ] && NAME="reloc_all"
      if [ ! -s "$OUT/${NAME}.csv" ]; then
        env $E python scripts/n_relocate.py --det "$OUT/grow/det_it0.parquet" \
          --grps "$OUT/screen_pass.csv" --half "$H" --out "$OUT/${NAME}" \
          > "logs/p_${TAG}_${NAME}.log" 2>&1 || exit 1
      fi
      tail -2 "logs/p_${TAG}_${NAME}.log"
    done
  fi

  STAGE="${TAG}_stability"
  if [ ! -s "$OUT/screen_final.csv" ]; then
    env $E python scripts/n_screen.py --stage stability --screen "$OUT/screen.csv" \
      --reloc2 "$OUT/reloc_h2.csv" --reloc1 "$OUT/reloc_h1.csv" \
      --reloc-all "$OUT/reloc_all.csv" --out "$OUT/screen_final.csv" --tag "$TAG" \
      > "logs/p_${TAG}_stability.log" 2>&1 || exit 1
  fi
  tail -12 "logs/p_${TAG}_stability.log"
}

echo "[P] ===================== STABILITY (clause e) ==========="
stability "$R" 0 real "$SUB"
stability "$NUL" 1 null "$SUB"

# ================================================================ 5. 2007-01 confirmation ======
STAGE=refire_dates
if [ ! -s "$C7/refire_meta.json" ]; then
  env N_BED=2007 N_STATION_SUBSET="$SUB" python - "$C7" "$N" <<'PY' > logs/p_refire_dates.log 2>&1 || exit 1
import sys, os, json, pandas as pd
sys.path.insert(0, "scripts")
import n_common as N
c7, nroot = sys.argv[1], sys.argv[2]
A = pd.read_csv(f"{nroot}/hour_audit_2007.csv")
W = A[~A.quiet_day]
sel = sorted(set(W[W.tremor_hour].date))
note = ""
if len(sel) < 3:
    sel2 = sorted(set(W[W.clause2].date))
    note = (f"the 2007 audit selects only {len(sel)} tremor days at clause 2 AND 3; "
            f"the clause-2-only day set ({len(sel2)} days) is used for the re-fire and RECORDED")
    sel = sel2
if len(sel) < 3:
    sel = N.window_days()
    note = (note + "; " if note else "") + \
        f"fewer than 3 days pass clause 2 either, so ALL {len(sel)} window days are used"
q = N.quiet_days()
open(f"{c7}/refire_dates.txt", "w").write(",".join(sel))
open(f"{c7}/quiet_dates.txt", "w").write(",".join(q))
json.dump(dict(n_tremor_days=len(sel), tremor_days=sel, n_quiet_days=len(q), quiet_days=q,
               note=note), open(f"{c7}/refire_meta.json", "w"), indent=1)
print(json.dumps(dict(n_tremor_days=len(sel), n_quiet_days=len(q), note=note), indent=1))
PY
fi
cat logs/p_refire_dates.log
D7T=$(cat "$C7/refire_dates.txt")
D7Q=$(cat "$C7/quiet_dates.txt")

confirm() {
  local OUT="$1" TAG="$2"
  local NP=0
  [ -s "$OUT/screen_pass.csv" ] && \
    NP=$(python -c "import pandas as pd;print(len(pd.read_csv('$OUT/screen_pass.csv')))")
  if [ "$NP" -eq 0 ]; then
    echo "[P] $TAG: nothing survived the 2008 screen -- no 2007 re-fire"
    return 0
  fi
  for W in tremor fa; do
    local DD="$D7T"; [ "$W" = "fa" ] && DD="$D7Q"
    STAGE="refire_${TAG}_${W}"
    if [ ! -s "$C7/det_${TAG}_${W}.parquet" ]; then
      env N_BED=2007 N_STATION_SUBSET="$SUB" python scripts/n_refire.py \
        --tpl "$OUT/grow/tpl_it0.npz" --grps "$OUT/screen_pass.csv" --dates "$DD" \
        --out "$C7" --tag "${TAG}_${W}" > "logs/p_refire_${TAG}_${W}.log" 2>&1 || exit 1
    fi
    tail -12 "logs/p_refire_${TAG}_${W}.log"
  done
  STAGE="reloc2007_${TAG}"
  if [ ! -s "$C7/loc2007_${TAG}.csv" ]; then
    env N_BED=2007 N_STATION_SUBSET="$SUB" python scripts/n_relocate.py \
      --det "$C7/det_${TAG}_tremor.parquet" --grps "$OUT/screen_pass.csv" \
      --out "$C7/loc2007_${TAG}" > "logs/p_reloc2007_${TAG}.log" 2>&1 || exit 1
  fi
  tail -2 "logs/p_reloc2007_${TAG}.log"
  STAGE="confirm_${TAG}"
  if [ ! -s "$C7/confirm_${TAG}.csv" ]; then
    env N_BED=2007 N_STATION_SUBSET="$SUB" python scripts/n_confirm.py \
      --screen "$OUT/screen_final.csv" --grps "$OUT/screen_pass.csv" \
      --det-tremor "$C7/det_${TAG}_tremor.parquet" --det-fa "$C7/det_${TAG}_fa.parquet" \
      --loc2007 "$C7/loc2007_${TAG}.csv" --tremor-dates "$C7/refire_dates.txt" \
      --quiet-dates "$C7/quiet_dates.txt" --out "$C7/confirm_${TAG}.csv" --tag "$TAG" \
      > "logs/p_confirm_${TAG}.log" 2>&1 || exit 1
  fi
  tail -20 "logs/p_confirm_${TAG}.log"
}

echo "[P] ===================== 2007-01 CONFIRMATION ==========="
confirm "$R" real
confirm "$NUL" null

STAGE=refire_meta
env N_BED=2007 python - "$C7" "$R" "$NUL" <<'PY' >> logs/p_refire_dates.log 2>&1 || true
import sys, json, os, pandas as pd
c7, r, n = sys.argv[1], sys.argv[2], sys.argv[3]
m = json.load(open(f"{c7}/refire_meta.json"))
for k, d in (("n_real", r), ("n_null", n)):
    p = f"{d}/screen_pass.csv"
    m[k] = int(len(pd.read_csv(p))) if os.path.exists(p) else 0
json.dump(m, open(f"{c7}/refire_meta.json", "w"), indent=1)
PY

# ================================================================ 6. the report ================
STAGE=report
env P_OUTROOT="$P" N_OUTROOT="$N" python scripts/p_report.py \
  --out results/WA_2008Q_RESULTS.md > logs/p_report.log 2>&1 || exit 1
cat logs/p_report.log

STAGE=notes
env P_OUTROOT="$P" python scripts/p_notes.py "$P" "$R" "$NUL" "$C7" "$BHD" \
  > logs/p_notes.log 2>&1 || true
tail -5 logs/p_notes.log

echo "OK $(date -u +%FT%TZ)" > results/P_PUGET2008Q_DONE
echo "[P] DONE $(date -u +%FT%TZ)"
