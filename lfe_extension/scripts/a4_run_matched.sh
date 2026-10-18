#!/bin/bash
# Track A-2 STEP 4 -- MATCHED-N control leg: re-stack the 2012 (previous-epoch) reference events with
# each family capped at the number of events the NEW leg used.  cc between two stacks rises with N on
# either side, and the 2018/2021 SSE windows yield fewer detections per family than the 2012 ETS, so
# the raw *_old control is optimistic; this leg removes that asymmetry.
cd /home/jovyan/tremorferometry/lfe_extension
S=/opt/conda/bin/python
m () {  # m <key> <sta> <netsta>
  local key=$1 sta=$2 ns=$3
  [ -f results/e5r_detections_a4_${key}_old.parquet ] || { echo "[wait] $key old detections"; return 0; }
  [ -f results/a4_epoch_identity_${key}_new.csv ] || { echo "[wait] $key new leg"; return 0; }
  $S scripts/a4_epoch_stack.py --sta $sta --netsta $ns --det results/e5r_detections_a4_${key}_old.parquet \
     --year 2012 --epoch 2012match --boundary matched --leg old \
     --match-n results/a4_epoch_identity_${key}_new.csv \
     --out results/a4_epoch_identity_${key}_oldm.csv > logs/a4_stack_${key}_oldm.log 2>&1
  echo "  $key/oldm exit=$?  $(grep 'families stacked' logs/a4_stack_${key}_oldm.log)"
}
m PGC PGC CN.PGC; m VGZ VGZ CN.VGZ; m NLLB NLLB CN.NLLB; m SNB SNB CN.SNB
m LRIV LRIV UW.LRIV; m VGZ17 VGZ CN.VGZ; m B926 B926 PB.B926
