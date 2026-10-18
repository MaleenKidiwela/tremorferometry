#!/bin/bash
# Track A-2 STEP 4 -- stack each withheld station at the predicted moveout of its reference events.
cd /home/jovyan/tremorferometry/lfe_extension
S=/opt/conda/bin/python
stack () {  # stack <key> <sta> <netsta> <leg> <year> <epoch> <boundary>
  local key=$1 sta=$2 ns=$3 leg=$4 yr=$5 ep=$6 bnd=$7
  local det=results/e5r_detections_a4_${key}_${leg}.parquet
  [ -f "$det" ] || { echo "[wait] $det missing"; return 1; }
  [ -f "results/a4_epoch_identity_${key}_${leg}.csv" ] && { echo "[skip] $key $leg"; return 0; }
  $S scripts/a4_epoch_stack.py --sta $sta --netsta $ns --det $det --year $yr --epoch $ep \
     --boundary $bnd --leg $leg --out results/a4_epoch_identity_${key}_${leg}.csv \
     > logs/a4_stack_${key}_${leg}.log 2>&1
  echo "  $key/$leg exit=$?"
}
case "${1:-}" in
 PGC)  stack PGC  PGC  CN.PGC  new 2018 2018_HH 2017-08-10; stack PGC  PGC  CN.PGC  old 2012 2012_BH 2017-08-10 ;;
 VGZ)  stack VGZ  VGZ  CN.VGZ  new 2018 2018_HH 2017-03-03; stack VGZ  VGZ  CN.VGZ  old 2012 2012_HH 2017-03-03 ;;
 NLLB) stack NLLB NLLB CN.NLLB new 2018 2018_HH 2017-02-23; stack NLLB NLLB CN.NLLB old 2012 2012_BH 2017-02-23 ;;
 SNB)  stack SNB  SNB  CN.SNB  new 2018 2018_HH 2017-06-08; stack SNB  SNB  CN.SNB  old 2012 2012_HH 2017-06-08 ;;
 LRIV) stack LRIV LRIV UW.LRIV new 2021 2021_HH 2019-05-30; stack LRIV LRIV UW.LRIV old 2012 2012_BH 2019-05-30 ;;
 VGZ17) stack VGZ17 VGZ CN.VGZ new 2017 2017_HH 2017-03-03; stack VGZ17 VGZ CN.VGZ old 2012 2012_HH 2017-03-03 ;;
 B926) stack B926 B926 PB.B926 new 2018 2018_EH 2017-08-10; stack B926 B926 PB.B926 old 2012 2012_EH 2017-08-10 ;;
 *) echo "usage: a4_run_stacks.sh PGC|VGZ|NLLB|SNB|LRIV|VGZ17" ;;
esac
