#!/usr/bin/env bash
# Y_WA_FINISH -- finish Washington.  Run the Puget / S-Olympic 2008 bed and the 2007 (CAFE at its
# fullest) bed AS ARRAYS, in the mode the data chooses, with the frozen recipe, a whole-chain
# site-shift null, the REBUILT locator, and the identical table for real and null.
#
#   nohup setsid -f bash scripts/y_wa_finish_queue.sh >> logs/y_wa_finish_queue.log 2>&1 < /dev/null
#
# Every stage pins its OWN environment with `env VAR=... python`; nothing is exported across
# stages.  Every stage has an output guard, so the queue resumes from where it stopped.
# Markers: results/Y_WA_DONE | results/Y_WA_FAILED; interim results/Y_WA_MODE_DONE after the
# decisive early number (stage 2); heartbeat results/Y_WA_HEARTBEAT rewritten each stage.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/y_wa_arrays beamcache_y

exec 9>locks/y_wa_finish.lock
flock -n 9 || { echo "[Y] already running"; exit 0; }
rm -f results/Y_WA_DONE results/Y_WA_FAILED
STAGE=start
beat() { STAGE="$1"; printf 'stage=%s at=%s pid=%s\n' "$STAGE" "$(date -u +%FT%TZ)" "$$" \
  > results/Y_WA_HEARTBEAT; echo "[Y] == $STAGE $(date -u +%FT%TZ)"; }
fail() { [ -f results/Y_WA_DONE ] || \
  printf 'stage=%s signal=%s at=%s pid=%s\n' "$STAGE" "${1:-exit}" "$(date -u +%FT%TZ)" "$$" \
  > results/Y_WA_FAILED; }
trap 'rc=$?; echo "[Y] TRAP EXIT rc=$rc stage=$STAGE"; fail "exit$rc"' EXIT
trap 'echo "[Y] TRAP TERM stage=$STAGE"; fail TERM; exit 143' TERM
trap 'echo "[Y] TRAP INT  stage=$STAGE"; fail INT;  exit 130' INT
trap 'echo "[Y] TRAP HUP  stage=$STAGE"; fail HUP;  exit 129' HUP

CU=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')
export LD_LIBRARY_PATH="${CU}${LD_LIBRARY_PATH:-}"

YR=results/y_wa_arrays
U=results/u_locate
BEDS="2008 2007"
GROW_MAX=1500          # a cap, not a target
CHUNK=750
MAX_SCREEN_PASS=200
MODE_HOURS=12
VAR=L9_xcdiff          # agent U's adopted locator
VAR2=L12_comb          # its companion
echo "[Y] start $(date -u +%FT%TZ) pid=$$ ppid=$PPID"

# =============================================================== 0. bed + sites ================
for B in $BEDS; do
  beat "bed_$B"
  if [ ! -s "$YR/bed_${B}.json" ] || [ ! -s "$YR/sites_${B}.csv" ]; then
    env Y_BED="$B" python - <<'PY' >> logs/y_bed.log 2>&1 || exit 1
import sys, json, os
sys.path.insert(0, "scripts")
import y_common as Y
b = Y.bed_summary()
json.dump(b, open(f"{Y.OUTROOT}/bed_{Y.BED}.json", "w"), indent=1, default=str)
os.makedirs(f"{Y.OUTROOT}/{Y.BED}_null", exist_ok=True)
json.dump(Y.null_offsets(Y.stations()),
          open(f"{Y.OUTROOT}/{Y.BED}_null/null_offsets.json", "w"), indent=1)
open(f"{Y.OUTROOT}/tremor_dates_{Y.BED}.txt", "w").write(",".join(Y.tremor_dates()))
open(f"{Y.OUTROOT}/quiet_dates_{Y.BED}.txt", "w").write(",".join(Y.quiet_days()))
print(json.dumps({k: v for k, v in b.items()
                  if k not in ("arrays", "quiet_days", "singles")}, indent=1, default=str))
PY
  fi
done
tail -80 logs/y_bed.log
TD08=$(cat "$YR/tremor_dates_2008.txt"); QD08=$(cat "$YR/quiet_dates_2008.txt")
TD07=$(cat "$YR/tremor_dates_2007.txt"); QD07=$(cat "$YR/quiet_dates_2007.txt")
echo "[Y] 2008 tremor days $(echo "$TD08"|tr ',' '\n'|wc -l) quiet days $(echo "$QD08"|tr ',' '\n'|wc -l);\
 2007 tremor days $(echo "$TD07"|tr ',' '\n'|wc -l) quiet days $(echo "$QD07"|tr ',' '\n'|wc -l)"

# =============================================================== 1. array traces ===============
# O_STEER=none: the zero-slowness delay-and-sum beam.  Agent O measured that the per-window
# slowness search raises the semblance by the SAME factor on hour-matched QUIET hours as on
# tremor hours (WA_2010_ARRAY_RESULTS sec.0c) -- it is the bias of taking a maximum over the
# grid, not signal -- so the beam is the equal-weight mean of the rms-normalised members.  No
# statics are applied: agent O's leave-one-event-out test showed they do not cross-validate.
for B in $BEDS; do
  beat "beam_$B"
  env Y_BED="$B" O_STEER=none python scripts/y_run.py o_beam --stage build \
    >> "logs/y_beam_${B}.log" 2>&1 || exit 1
  tail -3 "logs/y_beam_${B}.log"
done

# the slowness control, re-measured on THESE beds -- a diagnostic, never a decision (non-fatal)
for B in $BEDS; do
  beat "ctrl_slowness_$B"
  if [ ! -s "$YR/slowness_control_${B}.csv" ]; then
    env Y_BED="$B" python scripts/y_run.py o_controls --stage slowness --hours 4 \
      > "logs/y_ctrl_slowness_${B}.log" 2>&1 \
      && { [ -s "$YR/slowness_control.csv" ] && \
           mv "$YR/slowness_control.csv" "$YR/slowness_control_${B}.csv"; \
           [ -s "$YR/smoke_beam.json" ] && mv "$YR/smoke_beam.json" "$YR/smoke_beam_${B}.json"; } \
      || echo "[Y] slowness control failed on $B (diagnostic only, continuing)"
  fi
  tail -6 "logs/y_ctrl_slowness_${B}.log" 2>/dev/null || true
done

# =============================================================== 2. MODE BY MEASUREMENT ========
for B in $BEDS; do
  beat "mode_keep_$B"
  if [ ! -s "$YR/${B}/mode_smoke.json" ]; then
    mkdir -p "$YR/${B}"
    env Y_BED="$B" python scripts/y_mode.py --stage keep --out "$YR/${B}" \
      --hours "$MODE_HOURS" > "logs/y_mode_keep_${B}.log" 2>&1 || exit 1
    cat "logs/y_mode_keep_${B}.log"
    KEEP=$(cat "$YR/${B}/mode_keep.txt")
    for M in beam siteavg member; do
      for W in tremor quiet; do
        beat "mode_${B}_${M}_${W}"
        F="$YR/${B}/mode_hours_${W}.txt"
        KP="--keep $KEEP"; [ "$M" = "beam" ] && KP=""
        for HH in $(tr ';' ' ' < "$F"); do
          env Y_BED="$B" python scripts/y_run.py o_auto --mode "$M" --hours tremor \
            --only "$HH" $KP --out "$YR/${B}/mode/${M}_${W}" \
            >> "logs/y_mode_${B}_${M}.log" 2>&1 || exit 1
        done
      done
      tail -2 "logs/y_mode_${B}_${M}.log"
    done
    beat "mode_pick_$B"
    env Y_BED="$B" python scripts/y_mode.py --stage pick --out "$YR/${B}" \
      --hours "$MODE_HOURS" > "logs/y_mode_pick_${B}.log" 2>&1 || exit 1
  fi
  cat "logs/y_mode_pick_${B}.log"
done
MODE08=$(cat "$YR/2008/mode_chosen.txt")
MODE07=$(cat "$YR/2007/mode_chosen.txt")
echo "[Y] site treatment: 2008 -> $MODE08 ; 2007 -> $MODE07"
beat mode_done
env PYTHONUNBUFFERED=1 python scripts/y_report.py > logs/y_report.log 2>&1 || true
{ echo "mode_2008=$MODE08"; echo "mode_2007=$MODE07"; date -u +%FT%TZ; } > results/Y_WA_MODE_DONE
cat results/Y_WA_MODE_DONE

# =============================================================== 3. DISCOVERY ==================
chain() {
  local B="$1" NU="$2" MODE="$3" TD="$4" QD="$5"
  local TAG="$B" OUT="$YR/$B"
  if [ "$NU" = "1" ]; then TAG="${B}_null"; OUT="$YR/${B}_null"; fi
  local E="Y_BED=$B O_OUT=$OUT"
  if [ "$NU" = "1" ]; then E="$E O_NULL_SHIFT=1 O_NULL_SEED=11"; fi
  mkdir -p "$OUT"

  beat "${TAG}_auto_tremor"
  if [ ! -s "$OUT/auto/_params.json" ]; then
    env $E python scripts/y_run.py o_auto --mode "$MODE" --hours tremor --out "$OUT/auto" \
      > "logs/y_${TAG}_auto.log" 2>&1 || exit 1
  fi
  tail -3 "logs/y_${TAG}_auto.log"

  if [ "$NU" != "1" ]; then
    beat "${TAG}_auto_quiet"
    if [ ! -s "$OUT/auto_quiet/_params.json" ]; then
      env $E python scripts/y_run.py o_auto --mode "$MODE" --hours quiet --out "$OUT/auto_quiet" \
        > "logs/y_${TAG}_auto_quiet.log" 2>&1 || exit 1
    fi
    tail -3 "logs/y_${TAG}_auto_quiet.log"
  fi

  beat "${TAG}_pairstats"
  if [ ! -s "$OUT/pair_pairstats.json" ]; then
    if [ "$NU" = "1" ]; then
      env $E python scripts/y_run.py n_pairstats --tremor "$OUT/auto" --out "$OUT/pair" \
        > "logs/y_${TAG}_pairstats.log" 2>&1 || exit 1
    else
      env $E python scripts/y_run.py n_pairstats --tremor "$OUT/auto" \
        --quiet "$OUT/auto_quiet" --out "$OUT/pair" \
        > "logs/y_${TAG}_pairstats.log" 2>&1 || exit 1
    fi
  fi
  tail -45 "logs/y_${TAG}_pairstats.log"

  beat "${TAG}_group"
  if [ ! -s "$OUT/gb_summary.csv" ]; then
    env $E python scripts/y_run.py k_group --detfile "$OUT/pair_all.parquet" \
      --budget 100000000 --min-members 3 --times /nonexistent --out "$OUT/gb" \
      > "logs/y_${TAG}_group.log" 2>&1 || exit 1
  fi
  grep -E "nodes,|components" "logs/y_${TAG}_group.log" | tail -3 || true

  beat "${TAG}_split"
  if [ ! -s "$OUT/chunks/_split.json" ]; then
    env $E python scripts/n_chunk.py split --nodes "$OUT/gb_nodes.csv" \
      --summary "$OUT/gb_summary.csv" --out "$OUT/chunks" --size "$CHUNK" \
      --max-groups "$GROW_MAX" > "logs/y_${TAG}_split.log" 2>&1 || exit 1
  fi
  head -1 "logs/y_${TAG}_split.log"
  local NC
  NC=$(python -c "import json;print(json.load(open('$OUT/chunks/_split.json'))['n_chunks'])")
  case "$NC" in (''|*[!0-9]*) echo "[Y] bad chunk count '$NC'"; exit 1;; esac
  if [ "$NC" -eq 0 ]; then echo "[Y] $TAG: no groups to grow"; return 0; fi

  local i CH DIRS=""
  for ((i=0; i<NC; i++)); do
    CH=$(printf "%02d" "$i")
    beat "${TAG}_grow_c${CH}"
    if [ ! -s "$OUT/grow_c${CH}/det_it0.parquet" ]; then
      env $E python scripts/y_run.py k_grow --nodes "$OUT/chunks/chunk${CH}_nodes.csv" \
        --summary "$OUT/chunks/chunk${CH}_summary.csv" --dates "$TD" \
        --iters 0 --out "$OUT/grow_c${CH}" > "logs/y_${TAG}_grow_c${CH}.log" 2>&1 || exit 1
    fi
    tail -2 "logs/y_${TAG}_grow_c${CH}.log"
    DIRS="${DIRS}${OUT}/grow_c${CH},"
  done

  beat "${TAG}_merge"
  if [ ! -s "$OUT/grow/det_it0.parquet" ]; then
    env $E python scripts/n_chunk.py merge --dirs "$DIRS" --out "$OUT/grow" \
      > "logs/y_${TAG}_merge.log" 2>&1 || exit 1
  fi
  cat "logs/y_${TAG}_merge.log"

  beat "${TAG}_quiet_fa"
  if [ ! -s "$OUT/det_fa.parquet" ]; then
    env $E python scripts/y_run.py n_refire --tpl "$OUT/grow/tpl_it0.npz" --dates "$QD" \
      --out "$OUT" --tag fa > "logs/y_${TAG}_quietfa.log" 2>&1 || exit 1
  fi
  tail -10 "logs/y_${TAG}_quietfa.log"

  beat "${TAG}_screen_locate"
  if [ ! -s "$OUT/loc_it0.csv" ]; then
    env $E python scripts/y_run.py n_locate --tpl "$OUT/grow/tpl_it0.npz" \
      --det "$OUT/grow/det_it0.parquet" --out "$OUT/loc_it0.csv" \
      > "logs/y_${TAG}_locate.log" 2>&1 || exit 1
  fi
  tail -3 "logs/y_${TAG}_locate.log"

  beat "${TAG}_screen"
  if [ ! -s "$OUT/screen.csv" ]; then
    env $E python scripts/y_run.py n_screen --stage screen --loc "$OUT/loc_it0.csv" \
      --det "$OUT/grow/det_it0.parquet" --det-quiet "$OUT/det_fa.parquet" \
      --out "$OUT/screen.csv" --tag "$TAG" > "logs/y_${TAG}_screen.log" 2>&1 || exit 1
  fi
  tail -40 "logs/y_${TAG}_screen.log"

  # clause (e), stability -- the same LEDGER column the 2010 bed carries; non-fatal
  beat "${TAG}_passlist"
  if [ ! -s "$OUT/screen_pass.csv" ]; then
    env $E python - "$OUT" "$MAX_SCREEN_PASS" <<'PY' || return 0
import sys, pandas as pd
out, cap = sys.argv[1], int(sys.argv[2])
V = pd.read_csv(f"{out}/screen.csv")
P = V[V.pass_abcd.astype(bool)].sort_values("n_det", ascending=False).head(cap)
P[["grp"]].to_csv(f"{out}/screen_pass.csv", index=False)
print(f"{len(P)} templates pass clauses (a)-(d) of {len(V)}", flush=True)
PY
  fi
  local NP
  NP=$(python -c "import pandas as pd;print(len(pd.read_csv('$OUT/screen_pass.csv')))" 2>/dev/null) || NP=0
  echo "[Y] $TAG: $NP templates carried into the stability re-stack"
  if [ "${NP:-0}" -gt 0 ]; then
    for H in 2 1 0; do
      local NAME="reloc_h${H}"; [ "$H" = "0" ] && NAME="reloc_all"
      beat "${TAG}_${NAME}"
      if [ ! -s "$OUT/${NAME}.csv" ]; then
        env $E python scripts/y_run.py n_relocate --det "$OUT/grow/det_it0.parquet" \
          --grps "$OUT/screen_pass.csv" --half "$H" --out "$OUT/${NAME}" \
          > "logs/y_${TAG}_${NAME}.log" 2>&1 || break
      fi
      tail -2 "logs/y_${TAG}_${NAME}.log"
    done
    beat "${TAG}_stability"
    if [ ! -s "$OUT/screen_final.csv" ] && [ -s "$OUT/reloc_all.csv" ]; then
      env $E python scripts/y_run.py n_screen --stage stability --screen "$OUT/screen.csv" \
        --reloc2 "$OUT/reloc_h2.csv" --reloc1 "$OUT/reloc_h1.csv" \
        --reloc-all "$OUT/reloc_all.csv" --out "$OUT/screen_final.csv" --tag "$TAG" \
        > "logs/y_${TAG}_stability.log" 2>&1 || true
      tail -10 "logs/y_${TAG}_stability.log" 2>/dev/null || true
    fi
  fi
  return 0
}

echo "[Y] ===================== 2008 REAL ======================"
chain 2008 0 "$MODE08" "$TD08" "$QD08" || exit 1
echo "[Y] ===================== 2008 SITE-SHIFT NULL ==========="
chain 2008 1 "$MODE08" "$TD08" "$QD08" || exit 1
env PYTHONUNBUFFERED=1 python scripts/y_report.py >> logs/y_report.log 2>&1 || true
echo "[Y] ===================== 2007 REAL ======================"
chain 2007 0 "$MODE07" "$TD07" "$QD07" || exit 1
echo "[Y] ===================== 2007 SITE-SHIFT NULL ==========="
chain 2007 1 "$MODE07" "$TD07" "$QD07" || exit 1
env PYTHONUNBUFFERED=1 python scripts/y_report.py >> logs/y_report.log 2>&1 || true

# =============================================================== 4. THE REBUILT LOCATOR ========
# uu_stacks -> uu_picks -> uu_locate (L9_xcdiff adopted + L12_comb companion), with the 40-draw
# station bootstrap radius and the 20-draw shuffled-station control as COLUMNS.  NOT k_locate.
ulocate() {
  local B="$1" NU="$2"
  local TAG="$B" OUT="$YR/$B"
  if [ "$NU" = "1" ]; then TAG="${B}_null"; OUT="$YR/${B}_null"; fi
  [ -s "$OUT/grow/det_it0.parquet" ] || { echo "[Y] $TAG: nothing grown, no relocation"; return 0; }
  local E="Y_BED=$B O_OUT=$OUT"
  if [ "$NU" = "1" ]; then E="$E O_NULL_SHIFT=1 O_NULL_SEED=11"; fi

  beat "uu_${TAG}_stacks"
  if [ ! -s "$YR/uu_${TAG}_stacks.npz" ]; then
    env $E PYTHONUNBUFFERED=1 python scripts/y_run.py uu_stacks --bed wa \
      --det "$OUT/grow/det_it0.parquet" --tpl "$OUT/grow/tpl_it0.npz" \
      --out "$YR/uu_${TAG}_stacks.npz" > "logs/y_uu_${TAG}_stacks.log" 2>&1 \
      || { echo "[Y] $TAG stacks failed"; return 1; }
  fi
  tail -2 "logs/y_uu_${TAG}_stacks.log"

  beat "uu_${TAG}_picks"
  if [ ! -s "$YR/uu_${TAG}/picks.parquet" ]; then
    env $E PYTHONUNBUFFERED=1 python scripts/y_run.py uu_picks \
      --stacks "$YR/uu_${TAG}_stacks.npz" --out "$YR/uu_${TAG}/picks.parquet" \
      > "logs/y_uu_${TAG}_picks.log" 2>&1 || { echo "[Y] $TAG picks failed"; return 1; }
  fi
  tail -2 "logs/y_uu_${TAG}_picks.log"

  beat "uu_${TAG}_locate"
  if [ ! -s "$YR/uu_${TAG}/loc_${VAR2}.csv" ]; then
    env $E PYTHONUNBUFFERED=1 python scripts/y_run.py uu_locate --bed wa \
      --picks "$YR/uu_${TAG}/picks.parquet" --stacks "$YR/uu_${TAG}_stacks.npz" \
      --outdir "$YR/uu_${TAG}" --only "$VAR,$VAR2" --boot 40 --nshuf 20 \
      > "logs/y_uu_${TAG}_locate.log" 2>&1 || { echo "[Y] $TAG locate failed"; return 1; }
  fi
  tail -4 "logs/y_uu_${TAG}_locate.log"
  return 0
}

for B in $BEDS; do
  ulocate "$B" 0 || exit 1
  ulocate "$B" 1 || exit 1
done

# =============================================================== 5. SCORE ======================
for B in $BEDS; do
  for V in "$VAR" "$VAR2"; do
    beat "score_${B}_${V}"
    [ -s "$YR/uu_${B}/loc_${V}.csv" ] || { echo "[Y] no $B $V solution to score"; continue; }
    [ -s "$YR/uu_${B}_null/loc_${V}.csv" ] || { echo "[Y] no $B null $V solution"; continue; }
    if [ ! -s "$YR/score_${B}_${V}/y_score_${B}_${V}.json" ]; then
      env Y_BED="$B" PYTHONUNBUFFERED=1 python scripts/y_score.py \
        --real "$YR/uu_${B}/loc_${V}.csv" --null "$YR/uu_${B}_null/loc_${V}.csv" \
        --det-real "$YR/${B}/grow/det_it0.parquet" \
        --det-null "$YR/${B}_null/grow/det_it0.parquet" \
        --screen-real "$YR/${B}/screen.csv" --screen-null "$YR/${B}_null/screen.csv" \
        --out "$YR/score_${B}_${V}" --tag "${B}_${V}" \
        > "logs/y_score_${B}_${V}.log" 2>&1 || { echo "[Y] score $B $V failed"; continue; }
    fi
    cat "logs/y_score_${B}_${V}.log"
  done
done

# the agent-U Chestler table too, identical code, so the Olympic comparison is byte-identical to
# the way the 2010 bed was scored (sf_* ranking columns are 2010-only and are passed as empty)
for B in $BEDS; do
  beat "uu_wa_${B}"
  [ -s "$YR/uu_${B}/loc_${VAR}.csv" ] || continue
  [ -s "$YR/uu_${B}_null/loc_${VAR}.csv" ] || continue
  if [ ! -s "$YR/chest_${B}/wa_chestler_summary.csv" ]; then
    env Y_BED="$B" PYTHONUNBUFFERED=1 python scripts/y_run.py uu_wa \
      --real "$YR/uu_${B}/loc_${VAR}.csv" --null "$YR/uu_${B}_null/loc_${VAR}.csv" \
      --det-real "$YR/${B}/grow/det_it0.parquet" \
      --det-null "$YR/${B}_null/grow/det_it0.parquet" \
      --sf-real "" --sf-null "" --out "$YR/chest_${B}" \
      > "logs/y_uu_wa_${B}.log" 2>&1 || echo "[Y] uu_wa $B failed (non-fatal)"
  fi
  tail -20 "logs/y_uu_wa_${B}.log" 2>/dev/null || true
done

# =============================================================== 6. report =====================
beat report
env PYTHONUNBUFFERED=1 python scripts/y_report.py > logs/y_report.log 2>&1 \
  || { echo "[Y] report failed"; exit 1; }
tail -3 logs/y_report.log

date -u +"%Y-%m-%dT%H:%M:%SZ" > results/Y_WA_DONE
{ echo "mode_2008=$MODE08"; echo "mode_2007=$MODE07"; } >> results/Y_WA_DONE
beat done
echo "[Y] DONE $(date -u +%FT%TZ)"
