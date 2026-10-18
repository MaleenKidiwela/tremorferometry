#!/usr/bin/env bash
# STAGE 3 — forward extension to 2026-05, month by month, with per-month guards so the queue
# resumes exactly where it stopped.  Shelly's roster first (the southernmost located families),
# then Ducellier's.  Each month: fetch -> cache -> matched filter.  Tiers and the per-family
# per-year catalogue are re-derived after each source finishes its months.
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south; cd "$WS" || exit 1
mkdir -p results/ext logs
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
run_source(){
  SRC="$1"; M0="$2"; M1="$3"
  for M in $(python - "$M0" "$M1" <<'PY'
import sys,pandas as pd
print("\n".join(str(x)[:7] for x in pd.date_range(sys.argv[1],sys.argv[2],freq="MS")))
PY
); do
    OUTP="results/ext/${SRC}_${M}.parquet"
    [ -f "$OUTP" ] && continue
    PLAN="results/ext/plan_${SRC}_${M}.csv"
    python scripts/t16_plan.py --source "$SRC" --month "${M}-01" --out "$PLAN" || continue
    python scripts/t2_fetch.py "$PLAN" --workers 8 --manifest results/ext_download_manifest.csv || continue
    python scripts/t3_cache.py "$PLAN" --workers 24 --out "results/ext/cache_${SRC}_${M}.parquet" || continue
    NEXT=$(python - "$M" <<'PY'
import sys,pandas as pd; print(str((pd.Timestamp(sys.argv[1]+"-01")+pd.offsets.MonthBegin(1)).date()))
PY
)
    python scripts/t11_mf.py --source "$SRC" --t0 "${M}-01" --t1 "$NEXT" --out "$OUTP" || continue
    echo "[extend] $SRC $M done $(date -u +%FT%TZ)"
  done
  python scripts/t12_tiers.py --source "$SRC" --det "results/ext/${SRC}_*.parquet" || true
  python scripts/t13_report.py || true
}
run_source Shelly2025 2018-01 2026-05
run_source Ducellier2022 2008-01 2026-05
