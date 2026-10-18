#!/usr/bin/env bash
# O_OLY2010 -- reproduce the PUBLISHED Olympic LFE families (Chestler & Creager 2017 / Sweet
# et al. 2014, 2019) by running the 2010 N-Olympic Array-of-Arrays bed the way those papers ran
# it: the 1-km arrays AS ARRAYS.  PLAN.md sec.17.
#
#   nohup setsid -f bash scripts/o_oly2010_queue.sh >> logs/o_oly2010_queue.log 2>&1 < /dev/null
#
# Every stage pins its OWN environment with `env VAR=... python`; nothing is exported across
# stages (the K_YEAR leak that voided the first blind run, REPRO_BOSTOCK sec.9, cannot happen).
# Every stage has an output guard, so the queue resumes from where it stopped.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/o_oly2010

exec 9>locks/o_oly2010.lock
flock -n 9 || { echo "[O] already running"; exit 0; }
rm -f results/O_OLY2010_DONE results/O_OLY2010_FAILED
STAGE=start
fail() { [ -f results/O_OLY2010_DONE ] || \
  printf 'stage=%s signal=%s at=%s pid=%s\n' "$STAGE" "${1:-exit}" "$(date -u +%FT%TZ)" "$$" \
  > results/O_OLY2010_FAILED; }
trap 'rc=$?; echo "[O] TRAP EXIT rc=$rc stage=$STAGE"; fail "exit$rc"' EXIT
trap 'echo "[O] TRAP TERM stage=$STAGE"; fail TERM; exit 143' TERM
trap 'echo "[O] TRAP INT  stage=$STAGE"; fail INT;  exit 130' INT
trap 'echo "[O] TRAP HUP  stage=$STAGE"; fail HUP;  exit 129' HUP

CU=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')
export LD_LIBRARY_PATH="${CU}${LD_LIBRARY_PATH:-}"

O=results/o_oly2010
R=$O/2010
NUL=$O/2010_null
C11=$O/2011
ARMB=$O/chestler_arm
GROW_MAX=1500
CHUNK=750
MAX_SCREEN_PASS=200
mkdir -p "$R" "$NUL" "$C11" "$ARMB"
echo "[O] start $(date -u +%FT%TZ) pid=$$ ppid=$PPID"

# =============================================================== 0. bed + sites ================
STAGE=bed
if [ ! -s "$O/bed_2010.json" ] || [ ! -s "$O/bed_2011.json" ]; then
  for B in 2010 2011; do
    env O_BED="$B" python - <<'PY' >> logs/o_bed.log 2>&1 || exit 1
import sys, json, os
sys.path.insert(0, "scripts")
import o_common as O
b = O.bed_summary()
json.dump(b, open(f"{O.OUTROOT}/bed_{O.BED}.json", "w"), indent=1)
if O.BED == "2010":
    os.makedirs(f"{O.OUTROOT}/2010_null", exist_ok=True)
    json.dump(O.null_offsets(O.stations()),
              open(f"{O.OUTROOT}/2010_null/null_offsets.json", "w"), indent=1)
    open(f"{O.OUTROOT}/tremor_dates_2010.txt", "w").write(",".join(O.tremor_dates()))
    open(f"{O.OUTROOT}/quiet_dates_2010.txt", "w").write(",".join(O.quiet_days()))
else:
    open(f"{O.OUTROOT}/tremor_dates_2011.txt", "w").write(",".join(O.tremor_dates()))
    open(f"{O.OUTROOT}/quiet_dates_2011.txt", "w").write(",".join(O.quiet_days()))
print(json.dumps({k: v for k, v in b.items()
                  if k not in ("arrays", "quiet_days", "singles")}, indent=1))
PY
  done
fi
tail -60 logs/o_bed.log
TD10=$(cat "$O/tremor_dates_2010.txt"); QD10=$(cat "$O/quiet_dates_2010.txt")
TD11=$(cat "$O/tremor_dates_2011.txt"); QD11=$(cat "$O/quiet_dates_2011.txt")
echo "[O] 2010 tremor days $(echo "$TD10"|tr ',' '\n'|wc -l), quiet days $(echo "$QD10"|tr ',' '\n'|wc -l);\
 2011 tremor days $(echo "$TD11"|tr ',' '\n'|wc -l), quiet days $(echo "$QD11"|tr ',' '\n'|wc -l)"

# =============================================================== 1. the two beam controls ======
STAGE=ctrl_bands
if [ ! -s "$O/array_coherence_bands.csv" ]; then
  env O_BED=2010 python scripts/o_controls.py --stage bands --events 25 \
    > logs/o_ctrl_bands.log 2>&1 || exit 1
fi
tail -12 logs/o_ctrl_bands.log

STAGE=ctrl_statics
if [ ! -s "$O/statics_eq.json" ]; then
  env O_BED=2010 python scripts/o_beam.py --stage eqstatics --out "$O/statics_eq.json" \
    > logs/o_ctrl_statics.log 2>&1 || exit 1
fi
tail -12 logs/o_ctrl_statics.log

STAGE=ctrl_slowness
if [ ! -s "$O/smoke_beam.json" ]; then
  env O_BED=2010 python scripts/o_controls.py --stage slowness --hours 4 \
    > logs/o_ctrl_slowness.log 2>&1 || exit 1
fi
tail -12 logs/o_ctrl_slowness.log

# =============================================================== 2. array traces ===============
# O_STEER=none: the zero-slowness delay-and-sum beam.  The slowness search is retained in
# o_beam.py and is what O_STEER=grid runs; the control above is why it is not used.
for B in 2010 2011; do
  STAGE="beam_$B"
  env O_BED="$B" O_STEER=none python scripts/o_beam.py --stage build \
    >> "logs/o_beam_${B}.log" 2>&1 || exit 1
  tail -2 "logs/o_beam_${B}.log"
done

# =============================================================== 3. site-treatment smoke =======
STAGE=mode_smoke
if [ ! -s "$O/mode_smoke.json" ]; then
  env O_BED=2010 python - "$O" <<'PY' > logs/o_mode_keep.log 2>&1 || exit 1
import sys, os, json, pandas as pd
sys.path.insert(0, "scripts")
import o_common as O
# 12 tremor hours and 12 quiet hours, and the MEMBER set present on both populations, so the
# three treatments are compared on the same ground (the June quiet days carry 78 of the 96).
T = O.candidate_hours(); Q = O.quiet_hours()
T = T.iloc[[int(i) for i in pd.Series(range(len(T))).sample(12, random_state=0)]]
Q = Q.iloc[[int(i) for i in pd.Series(range(len(Q))).sample(12, random_state=0)]]
idx = O.member_cache_index()
def have(dates):
    s = None
    for d in dates:
        ts = pd.Timestamp(d)
        k = {m for (m, y, doy) in idx if (y, doy) == (ts.year, ts.dayofyear)}
        s = k if s is None else (s | k)
    return s
keep = have(sorted(set(T.date))) & have(sorted(set(Q.date)))
out = sys.argv[1]
open(f"{out}/mode_keep.txt", "w").write(",".join(sorted(keep)))
open(f"{out}/mode_hours_tremor.txt", "w").write(
    ";".join(f"{r.date}:{int(r.hour)}" for r in T.itertuples()))
open(f"{out}/mode_hours_quiet.txt", "w").write(
    ";".join(f"{r.date}:{int(r.hour)}" for r in Q.itertuples()))
print(f"{len(keep)} member stations common to both populations")
PY
  cat logs/o_mode_keep.log
  KEEP=$(cat "$O/mode_keep.txt")
  for M in beam siteavg member; do
    for W in tremor quiet; do
      STAGE="mode_${M}_${W}"
      F="$O/mode_hours_${W}.txt"
      KP="--keep $KEEP"; [ "$M" = "beam" ] && KP=""
      for HH in $(tr ';' ' ' < "$F"); do
        env O_BED=2010 python scripts/o_auto.py --mode "$M" --hours tremor --only "$HH" \
          $KP --out "$O/mode/${M}_${W}" >> "logs/o_mode_${M}.log" 2>&1 || exit 1
      done
    done
    tail -2 "logs/o_mode_${M}.log"
  done
  STAGE=mode_pick
  env O_BED=2010 python - "$O" <<'PY' > logs/o_mode_pick.log 2>&1 || exit 1
import sys, os, glob, json, pandas as pd, numpy as np
out = sys.argv[1]
rows = []
for m in ("member", "beam", "siteavg"):
    r = {}
    for w in ("tremor", "quiet"):
        d = f"{out}/mode/{m}_{w}"
        S = pd.read_csv(f"{d}/_stats.csv") if os.path.exists(f"{d}/_stats.csv") else None
        if S is None:
            continue
        ok = S[S.nmad_in_sd > 3]
        r[f"n_blocks_{w}"] = int(len(ok))
        r[f"det_per_hour_{w}"] = float(ok.n_det.sum() / max(len(ok), 1))
        r[f"traces_{w}"] = float(ok.n_sta.median())
        r[f"thr_sd_{w}"] = float(ok.nmad_in_sd.median())
    if "det_per_hour_quiet" in r and r["det_per_hour_quiet"] > 0:
        r["tremor_quiet_ratio"] = r["det_per_hour_tremor"] / r["det_per_hour_quiet"]
    else:
        r["tremor_quiet_ratio"] = float("inf") if r.get("det_per_hour_tremor", 0) > 0 else 0.0
    rows.append(dict(mode=m, **r))
T = pd.DataFrame(rows)
best = max(rows, key=lambda x: (x["tremor_quiet_ratio"]
                                if np.isfinite(x["tremor_quiet_ratio"]) else -1))
why = ("highest tremor/quiet pair-count ratio at 8 x MAD on 12 + 12 hours restricted to the "
       "member set present in both populations; the 2010 run's whole-bed value was 4.2x and "
       "the north's truth-based value 11x")
json.dump(dict(table=rows, chosen=best["mode"], why=why),
          open(f"{out}/mode_smoke.json", "w"), indent=1, default=float)
open(f"{out}/mode_chosen.txt", "w").write(best["mode"])
print(T.to_string()); print("CHOSEN:", best["mode"])
PY
fi
cat logs/o_mode_pick.log
MODE=$(cat "$O/mode_chosen.txt")
echo "[O] site treatment for the full chain: $MODE"

# =============================================================== 4. one discovery chain ========
chain() {
  local OUT="$1" NU="$2" TAG="$3"
  local E="O_BED=2010 O_OUT=$OUT"
  if [ "$NU" = "1" ]; then E="$E O_NULL_SHIFT=1 O_NULL_SEED=11"; fi
  mkdir -p "$OUT"

  STAGE="${TAG}_auto_tremor"
  if [ ! -s "$OUT/auto/_params.json" ]; then
    env $E python scripts/o_auto.py --mode "$MODE" --hours tremor --out "$OUT/auto" \
      > "logs/o_${TAG}_auto.log" 2>&1 || exit 1
  fi
  tail -3 "logs/o_${TAG}_auto.log"

  if [ "$NU" != "1" ]; then
    STAGE="${TAG}_auto_quiet"
    if [ ! -s "$OUT/auto_quiet/_params.json" ]; then
      env $E python scripts/o_auto.py --mode "$MODE" --hours quiet --out "$OUT/auto_quiet" \
        > "logs/o_${TAG}_auto_quiet.log" 2>&1 || exit 1
    fi
    tail -3 "logs/o_${TAG}_auto_quiet.log"
  fi

  STAGE="${TAG}_pairstats"
  if [ ! -s "$OUT/pair_pairstats.json" ]; then
    if [ "$NU" = "1" ]; then
      env $E python scripts/o_run.py n_pairstats --tremor "$OUT/auto" --out "$OUT/pair" \
        > "logs/o_${TAG}_pairstats.log" 2>&1 || exit 1
    else
      env $E python scripts/o_run.py n_pairstats --tremor "$OUT/auto" \
        --quiet "$OUT/auto_quiet" --out "$OUT/pair" \
        > "logs/o_${TAG}_pairstats.log" 2>&1 || exit 1
    fi
  fi
  tail -45 "logs/o_${TAG}_pairstats.log"

  STAGE="${TAG}_group"
  if [ ! -s "$OUT/gb_summary.csv" ]; then
    env $E python scripts/o_run.py k_group --detfile "$OUT/pair_all.parquet" \
      --budget 100000000 --min-members 3 --times /nonexistent --out "$OUT/gb" \
      > "logs/o_${TAG}_group.log" 2>&1 || exit 1
  fi
  grep -E "nodes,|components" "logs/o_${TAG}_group.log" | tail -3 || true

  STAGE="${TAG}_split"
  if [ ! -s "$OUT/chunks/_split.json" ]; then
    env $E python scripts/n_chunk.py split --nodes "$OUT/gb_nodes.csv" \
      --summary "$OUT/gb_summary.csv" --out "$OUT/chunks" --size "$CHUNK" \
      --max-groups "$GROW_MAX" > "logs/o_${TAG}_split.log" 2>&1 || exit 1
  fi
  head -1 "logs/o_${TAG}_split.log"
  local NC
  NC=$(python -c "import json;print(json.load(open('$OUT/chunks/_split.json'))['n_chunks'])")
  case "$NC" in (''|*[!0-9]*) echo "[O] bad chunk count '$NC'"; exit 1;; esac
  if [ "$NC" -eq 0 ]; then echo "[O] $TAG: no groups to grow"; return 0; fi

  local i CH DIRS=""
  for ((i=0; i<NC; i++)); do
    CH=$(printf "%02d" "$i")
    STAGE="${TAG}_grow_c${CH}"
    if [ ! -s "$OUT/grow_c${CH}/det_it0.parquet" ]; then
      env $E python scripts/o_run.py k_grow --nodes "$OUT/chunks/chunk${CH}_nodes.csv" \
        --summary "$OUT/chunks/chunk${CH}_summary.csv" --dates "$TD10" \
        --iters 0 --out "$OUT/grow_c${CH}" > "logs/o_${TAG}_grow_c${CH}.log" 2>&1 || exit 1
    fi
    tail -2 "logs/o_${TAG}_grow_c${CH}.log"
    DIRS="${DIRS}${OUT}/grow_c${CH},"
  done

  STAGE="${TAG}_merge"
  if [ ! -s "$OUT/grow/det_it0.parquet" ]; then
    env $E python scripts/n_chunk.py merge --dirs "$DIRS" --out "$OUT/grow" \
      > "logs/o_${TAG}_merge.log" 2>&1 || exit 1
  fi
  cat "logs/o_${TAG}_merge.log"

  STAGE="${TAG}_quiet_fa"
  if [ ! -s "$OUT/det_fa.parquet" ]; then
    env $E python scripts/o_run.py n_refire --tpl "$OUT/grow/tpl_it0.npz" --dates "$QD10" \
      --out "$OUT" --tag fa > "logs/o_${TAG}_quietfa.log" 2>&1 || exit 1
  fi
  tail -10 "logs/o_${TAG}_quietfa.log"

  STAGE="${TAG}_locate"
  if [ ! -s "$OUT/loc_it0.csv" ]; then
    env $E python scripts/o_run.py n_locate --tpl "$OUT/grow/tpl_it0.npz" \
      --det "$OUT/grow/det_it0.parquet" --out "$OUT/loc_it0.csv" \
      > "logs/o_${TAG}_locate.log" 2>&1 || exit 1
  fi
  tail -3 "logs/o_${TAG}_locate.log"

  STAGE="${TAG}_screen"
  if [ ! -s "$OUT/screen.csv" ]; then
    env $E python scripts/o_run.py n_screen --stage screen --loc "$OUT/loc_it0.csv" \
      --det "$OUT/grow/det_it0.parquet" --det-quiet "$OUT/det_fa.parquet" \
      --out "$OUT/screen.csv" --tag "$TAG" > "logs/o_${TAG}_screen.log" 2>&1 || exit 1
  fi
  tail -50 "logs/o_${TAG}_screen.log"
}

echo "[O] ===================== REAL CHAIN ====================="
chain "$R" 0 real
echo "[O] ===================== SITE-SHIFT NULL ================"
chain "$NUL" 1 null

# =============================================================== 5. clause (e): stability ======
stability() {
  local OUT="$1" NU="$2" TAG="$3"
  local E="O_BED=2010 O_OUT=$OUT"
  if [ "$NU" = "1" ]; then E="$E O_NULL_SHIFT=1 O_NULL_SEED=11"; fi
  [ -s "$OUT/screen.csv" ] || { echo "[O] $TAG: no screen table"; return 0; }
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
  echo "[O] $TAG: $NP templates carried into the stability re-stack"
  if [ "$NP" -gt 0 ]; then
    for H in 2 1 0; do
      local NAME="reloc_h${H}"; [ "$H" = "0" ] && NAME="reloc_all"
      STAGE="${TAG}_${NAME}"
      if [ ! -s "$OUT/${NAME}.csv" ]; then
        env $E python scripts/o_run.py n_relocate --det "$OUT/grow/det_it0.parquet" \
          --grps "$OUT/screen_pass.csv" --half "$H" --out "$OUT/${NAME}" \
          > "logs/o_${TAG}_${NAME}.log" 2>&1 || exit 1
      fi
      tail -2 "logs/o_${TAG}_${NAME}.log"
    done
  fi
  STAGE="${TAG}_stability"
  if [ ! -s "$OUT/screen_final.csv" ]; then
    env $E python scripts/o_run.py n_screen --stage stability --screen "$OUT/screen.csv" \
      --reloc2 "$OUT/reloc_h2.csv" --reloc1 "$OUT/reloc_h1.csv" \
      --reloc-all "$OUT/reloc_all.csv" --out "$OUT/screen_final.csv" --tag "$TAG" \
      > "logs/o_${TAG}_stability.log" 2>&1 || exit 1
  fi
  tail -10 "logs/o_${TAG}_stability.log"
}
echo "[O] ===================== STABILITY (clause e) ==========="
stability "$R" 0 real
stability "$NUL" 1 null

# =============================================================== 6. the Chestler check =========
STAGE=chestler
if [ ! -s "$O/chestler_check.json" ]; then
  env O_BED=2010 python scripts/o_chest_check.py --screen "$R/screen_final.csv" \
    --loc "$R/loc_it0.csv" --loc-null "$NUL/loc_it0.csv" --out "$O/chestler_check.json" \
    > logs/o_chestler.log 2>&1 || exit 1
fi
cat logs/o_chestler.log

# =============================================================== 7. 2011 confirmation ==========
confirm() {
  local OUT="$1" TAG="$2"
  local NP=0
  [ -s "$OUT/screen_pass.csv" ] && \
    NP=$(python -c "import pandas as pd;print(len(pd.read_csv('$OUT/screen_pass.csv')))")
  if [ "$NP" -eq 0 ]; then echo "[O] $TAG: nothing survived the screen -- no 2011 re-fire";
    return 0; fi
  for W in tremor fa; do
    local DD="$TD11"; [ "$W" = "fa" ] && DD="$QD11"
    STAGE="refire_${TAG}_${W}"
    if [ ! -s "$C11/det_${TAG}_${W}.parquet" ]; then
      env O_BED=2011 python scripts/o_run.py n_refire --tpl "$OUT/grow/tpl_it0.npz" \
        --grps "$OUT/screen_pass.csv" --dates "$DD" --out "$C11" --tag "${TAG}_${W}" \
        > "logs/o_refire_${TAG}_${W}.log" 2>&1 || exit 1
    fi
    tail -10 "logs/o_refire_${TAG}_${W}.log"
  done
  STAGE="reloc2011_${TAG}"
  if [ ! -s "$C11/loc2011_${TAG}.csv" ]; then
    env O_BED=2011 python scripts/o_run.py n_relocate --det "$C11/det_${TAG}_tremor.parquet" \
      --grps "$OUT/screen_pass.csv" --out "$C11/loc2011_${TAG}" \
      > "logs/o_reloc2011_${TAG}.log" 2>&1 || exit 1
  fi
  tail -2 "logs/o_reloc2011_${TAG}.log"
  STAGE="confirm_${TAG}"
  if [ ! -s "$C11/confirm_${TAG}.csv" ]; then
    env O_BED=2011 python scripts/o_run.py n_confirm --screen "$OUT/screen_final.csv" \
      --grps "$OUT/screen_pass.csv" --det-tremor "$C11/det_${TAG}_tremor.parquet" \
      --det-fa "$C11/det_${TAG}_fa.parquet" --loc2007 "$C11/loc2011_${TAG}.csv" \
      --tremor-dates "$O/tremor_dates_2011.txt" --quiet-dates "$O/quiet_dates_2011.txt" \
      --out "$C11/confirm_${TAG}.csv" --tag "$TAG" \
      > "logs/o_confirm_${TAG}.log" 2>&1 || exit 1
  fi
  tail -18 "logs/o_confirm_${TAG}.log"
}
echo "[O] ===================== 2011 CONFIRMATION =============="
confirm "$R" real
confirm "$NUL" null

# =============================================================== 8. report (arm A) =============
STAGE=report_a
env O_BED=2010 python scripts/o_report.py --out results/WA_2010_ARRAY_RESULTS.md \
  > logs/o_report.log 2>&1 || exit 1
cat logs/o_report.log
echo "OK-ARM-A $(date -u +%FT%TZ)" > results/O_OLY2010_ARMA_DONE

# =============================================================== 9. ARM B (non-fatal) ==========
# Chestler's OWN chain (beam 18xMAD peaks -> cross-array boxcar stack -> 15 s 2-8 Hz seeds ->
# 2nd-root network stack at 45 x MAD -> 150-window templates, 3 iterations -> R_com dedupe).
# Failures here must not cost the arm-A report, so every stage is guarded with `|| ARMB_OK=0`.
ARMB_OK=1
armb() { [ "$ARMB_OK" = "1" ] || return 0; STAGE="armB_$1"; shift; "$@" || { ARMB_OK=0;
  echo "[O] ARM B stage failed, continuing"; }; }

echo "[O] ===================== ARM B: Chestler's chain ========"
if [ ! -s "$ARMB/seeds.csv" ]; then
  armb seeds env O_BED=2010 O_OUT="$ARMB" python scripts/o_chestler.py --stage seeds \
    --out "$ARMB" > logs/o_armb_seeds.log 2>&1
fi
tail -4 logs/o_armb_seeds.log 2>/dev/null || true
if [ "$ARMB_OK" = "1" ] && [ ! -s "$ARMB/det_it3.parquet" ]; then
  armb grow env O_BED=2010 O_OUT="$ARMB" python scripts/o_chestler.py --stage grow \
    --out "$ARMB" --max-seeds 1200 > logs/o_armb_grow.log 2>&1
fi
tail -6 logs/o_armb_grow.log 2>/dev/null || true
if [ "$ARMB_OK" = "1" ] && [ ! -s "$ARMB/keep.csv" ]; then
  armb dedupe env O_BED=2010 O_OUT="$ARMB" python scripts/o_chestler.py --stage dedupe \
    --out "$ARMB" > logs/o_armb_dedupe.log 2>&1
fi
tail -4 logs/o_armb_dedupe.log 2>/dev/null || true
if [ "$ARMB_OK" = "1" ] && [ ! -s "$ARMB/det_fa.parquet" ]; then
  armb fa env O_BED=2010 O_OUT="$ARMB" python scripts/o_chestler.py --stage fa \
    --out "$ARMB" > logs/o_armb_fa.log 2>&1
fi
if [ "$ARMB_OK" = "1" ] && [ ! -s "$ARMB/loc.csv" ]; then
  armb locate env O_BED=2010 O_OUT="$ARMB" python scripts/o_run.py n_locate \
    --tpl "$ARMB/tpl_final.npz" --det "$ARMB/det_final.parquet" --out "$ARMB/loc.csv" \
    > logs/o_armb_locate.log 2>&1
fi
tail -3 logs/o_armb_locate.log 2>/dev/null || true
if [ "$ARMB_OK" = "1" ] && [ ! -s "$ARMB/screen.csv" ]; then
  armb screen env O_BED=2010 O_OUT="$ARMB" python scripts/o_run.py n_screen --stage screen \
    --loc "$ARMB/loc.csv" --det "$ARMB/det_final.parquet" --det-quiet "$ARMB/det_fa.parquet" \
    --out "$ARMB/screen.csv" --tag armB > logs/o_armb_screen.log 2>&1
fi
tail -40 logs/o_armb_screen.log 2>/dev/null || true
if [ "$ARMB_OK" = "1" ] && [ ! -s "$ARMB/chestler_check.json" ]; then
  armb chest env O_BED=2010 python scripts/o_chest_check.py --screen "$ARMB/screen.csv" \
    --loc "$ARMB/loc.csv" --loc-null "$NUL/loc_it0.csv" \
    --out "$ARMB/chestler_check.json" > logs/o_armb_chest.log 2>&1
fi
tail -30 logs/o_armb_chest.log 2>/dev/null || true
STAGE=armB_summary
env O_BED=2010 python - "$ARMB" <<'PY' >> logs/o_armb_sum.log 2>&1 || true
import sys, os, json, pandas as pd
d = sys.argv[1]
s = {}
for k, p in (("seeds", "seeds.csv"), ("keep", "keep.csv"), ("screen", "screen.csv")):
    if os.path.exists(f"{d}/{p}"):
        s[f"n_{k}"] = int(len(pd.read_csv(f"{d}/{p}")))
for it in range(4):
    p = f"{d}/det_it{it}.parquet"
    if os.path.exists(p):
        D = pd.read_parquet(p)
        s[f"iter{it}"] = dict(n_det=int(len(D)), n_templates=int(D.grp.nunique()) if len(D) else 0)
for k, p in (("dedupe", "dedupe.json"), ("screen_summary", "screen_summary.json"),
             ("chestler", "chestler_check.json")):
    if os.path.exists(f"{d}/{p}"):
        s[k] = json.load(open(f"{d}/{p}"))
json.dump(s, open(f"{d}/summary.json", "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in s.items() if k != "screen_summary"}, indent=1)[:3000])
PY
tail -40 logs/o_armb_sum.log

# =============================================================== 10. report + notes ============
STAGE=report
env O_BED=2010 python scripts/o_report.py --out results/WA_2010_ARRAY_RESULTS.md \
  > logs/o_report.log 2>&1 || exit 1
cat logs/o_report.log

STAGE=notes
env O_BED=2010 python - <<'PY' >> logs/o_notes.log 2>&1 || true
import sys, json, os, datetime
sys.path.insert(0, "scripts")
import o_common as O
R = O.OUTROOT
def jl(p, d=None):
    try: return json.load(open(p))
    except Exception: return d
bed = jl(f"{R}/bed_2010.json", {}); mode = jl(f"{R}/mode_smoke.json", {})
ps = jl(f"{R}/2010/pair_pairstats.json", {})
sr = jl(f"{R}/2010/screen_final_summary.json", {}); sn = jl(f"{R}/2010_null/screen_final_summary.json", {})
ch = jl(f"{R}/chestler_check.json", {}); cr = jl(f"{R}/2011/confirm_real_summary.json", {})
cn = jl(f"{R}/2011/confirm_null_summary.json", {}); ab = jl(f"{R}/chestler_arm/summary.json", {})
txt = f"""
## {datetime.date.today()} -- agent O: the 2010 N-Olympic bed run with the arrays AS ARRAYS
Goal (PLAN.md sec.17): reproduce Chestler & Creager 2017 / Sweet et al.  Sites defined by geometry:
{bed.get('n_member_stations')} stations -> {bed.get('n_sites')} sites ({bed.get('n_array_sites')} arrays
covering {bed.get('n_members_in_arrays')} stations + {bed.get('n_single_sites')} singles).
Beam = zero-slowness delay-and-sum; the Ghosh slowness grid was measured to give the SAME semblance
gain on quiet hours as on tremor hours (max-over-grid bias), and intraslab-earthquake S waves show the
1-km arrays are only ~25 % coherent at 1-2 Hz and incoherent above 4 Hz.
Site treatment chosen by measurement: {mode.get('chosen')} ({[r for r in mode.get('table', [])]}).
Pair stage: tremor/quiet ratio {ps.get('tremor_quiet_ratio')} (2010 member-wise run 4.2x, north 11x).
Screen: real candidates {sr.get('n_candidate')} vs null {sn.get('n_candidate')}.
Chestler check: {json.dumps(ch.get('A', {}), default=float)} / {json.dumps(ch.get('B', {}), default=float)}.
2011 confirmation: real {cr.get('n_CONFIRMED')} / null {cn.get('n_CONFIRMED')}.
Arm B (Chestler's own chain): {json.dumps({k: v for k, v in ab.items() if not isinstance(v, dict)}, default=float)}.
Full ledger results/WA_2010_ARRAY_RESULTS.md; workspace results/o_oly2010/.
"""
open(f"{O.STUDY}/NOTES.md", "a").write(txt)
print(txt)
PY
tail -25 logs/o_notes.log

echo "OK $(date -u +%FT%TZ)" > results/O_OLY2010_DONE
echo "[O] DONE $(date -u +%FT%TZ)"
