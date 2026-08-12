#!/bin/bash
# Battery gate AFTER causality + dv/v (user-mandated second gate, 2026-07-11). Runs the Tier-1 trust
# battery (stack-vs-random coda-sigma) forward-only (NO --rev-mf -> script's f50=1.0/f95=6.0 fallback,
# since reverse densify is retired), then recomputes the Z dv/v over families that are BOTH causality-
# reliable AND not battery-FAIL. NEEDS TRACES (harvests raw windows) -> MUST run before trace deletion.
# Usage: scripts/battery_gate.sh <STA>
cd /home/jovyan/tremorferometry
STA=$(echo "$1" | tr '[:lower:]' '[:upper:]'); sta=$(echo "$STA" | tr '[:upper:]' '[:lower:]')
PY=/home/jovyan/envs/tremorferometry/bin/python
STATUS=logs/rollout_status.log; log(){ echo "$(date +%H:%M) $*" | tee -a "$STATUS"; }
TAG=${STA}p90f40

[ -f data/${sta}_causality_cert.csv ] || { log "$STA battery SKIP: no causality_cert (finalize first)"; exit 0; }

# 1. Tier-1 battery on all families (forward-only) -> data/family_trust_tier1_<STA>.csv
if [ ! -f data/family_trust_tier1_${STA}.csv ]; then
  log "$STA battery start (Tier-1 stack-vs-random, forward-only)"
  $PY scripts/family_trust_tier1.py --net PB --sta "$STA" --mf "data/mf_${sta}p90f40_*.csv" \
     --workers 12 > logs/battery_${sta}.log 2>&1 || { log "$STA ABORT: battery failed"; exit 20; }
fi
# 2. Doubly-gated dv/v (causality ∩ battery-not-FAIL); never destroys the causality dv/v
$PY scripts/battery_gate_dvv.py "$STA" "$TAG" >> logs/battery_${sta}.log 2>&1 \
   || { log "$STA ABORT: battery_gate_dvv failed"; exit 21; }
gs=$($PY -c "import json;d=json.load(open('data/${sta}_battery_gate_summary.json'));print('gated %d/%d std %s%% (TRUSTED %d UNDET %d FAIL %d not-scored %d)'%(d['gated_n'],d['causality_cert'],d['gated_std_pct'],d['battery_trusted'],d['battery_undet'],d['battery_fail'],d['cert_not_scored']))" 2>/dev/null)
log "$STA BATTERY GATE done: $gs"
