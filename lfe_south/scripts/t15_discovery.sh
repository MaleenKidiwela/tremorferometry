#!/usr/bin/env bash
# STAGE R2 — blind discovery on the SOUTHERN permanent bed with the reproduced Bostock recipe
# (ported unchanged as ks_*.py; the only edit is the 1-C/2-C loader fidelity fix).
# Bed = the Ducellier 2008-04 window already on disk; truth = his 66 families.
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south; cd "$WS" || exit 1
OUT=results/discovery_south; mkdir -p "$OUT"
STAS=$(python - <<'PY'
import pandas as pd
A=pd.read_csv('results/duc_bedA_stations.csv')
print(",".join(sorted(A.sta.unique())))
PY
)
DATES=$(python - <<'PY'
import pandas as pd
print(",".join(str(d.date()) for d in pd.date_range('2008-04-15','2008-05-05')))
PY
)
echo "[discovery] stations: $STAS"
if [ ! -f "$OUT/auto.parquet" ] && [ ! -d "$OUT/auto" ]; then
  python scripts/ks_auto.py --dates "$DATES" --stas "$STAS" --out "$OUT/auto" || exit 1
fi
[ -f "$OUT/groups.parquet" ] || python scripts/ks_group.py --auto "$OUT/auto" --out "$OUT/groups.parquet" || exit 1
[ -f "$OUT/grown.parquet" ] || python scripts/ks_grow.py --groups "$OUT/groups.parquet" --out "$OUT/grown.parquet" || exit 1
python scripts/t15_score_discovery.py --grown "$OUT/grown.parquet" --out results/south_discovery_summary.csv
