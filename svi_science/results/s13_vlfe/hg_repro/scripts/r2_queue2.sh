#!/bin/bash
# Amendment R2/R2a stages 2-5 (Table-S1 targets) (chained after R2_VAR_DONE): calibrate + choose on 1 Jul-6 Sep 2011 -> [V7 combination run if
# the rule asks] -> held-out test + detection match vs Table S1 -> matched filter at the chosen setting (real + scrambled controls + 100 noise templates;
# skipped when the chosen setting is V0/V8, whose R1 outputs come from the identical code path) -> stage-2 analysis -> Table-S1 own-detection tests (hg_r2_s1.py).
# Markers: R2_CALIB_DONE R2_HELDOUT_DONE R2_MF_DONE R2_DONE.  <= 8 workers, OMP 1.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while [ ! -f R2_VAR_DONE ]; do
  if ! pgrep -f "bash scripts/r2_queue1.sh" > /dev/null; then echo "queue1 ended without R2_VAR_DONE"; exit 1; fi
  sleep 60
done
python scripts/hg_r2_calib.py > logs/r2_calib.log 2>&1 || exit 1
V7=$(python -c "import json;print(json.load(open('tables/r2_choice.json'))['V7_needed'])")
if [ "$V7" = "True" ]; then
  read SPEC FILE <<< "$(python scripts/hg_r2_calib.py --v7spec)"
  echo "V7 spec $SPEC file $FILE"
  for K in 1 2 3 4 5; do python scripts/hg_r2.py real $K "$SPEC" & done; wait
  python scripts/hg_r2_calib.py "$FILE" > logs/r2_calib_v7.log 2>&1 || exit 1
fi
touch R2_CALIB_DONE
python scripts/hg_r2_heldout.py > logs/r2_heldout.log 2>&1 || exit 1
python scripts/hg_r2_match.py > logs/r2_match.log 2>&1 || exit 1
touch R2_HELDOUT_DONE
SPEC=$(python -c "import json;print(json.load(open('tables/r2_choice.json'))['chosen_spec'].split('+')[0])")
CH=$(python -c "import json;print(json.load(open('tables/r2_choice.json'))['chosen_variant'])")
if [ "$CH" != "V0" ] && [ "$CH" != "V8" ]; then
  NORM=$(echo $SPEC | cut -d, -f3)
  if [ "$NORM" = "net" ]; then P=3; else P=5; fi
  ls data/r2/real_*_T?.npz 2>/dev/null
  for K in 1 2 3 4 5; do python scripts/hg_r2.py real $K "$SPEC" & 
    if [ $K -eq $P ]; then wait; fi; done; wait
  for c in 0 13 26 39 52 65 78 91; do e=$((c+13)); [ $e -gt 100 ] && e=100; python scripts/hg_r2.py noise $c $e "$SPEC" > logs/r2_noise_$c.log 2>&1 & done; wait
fi
touch R2_MF_DONE
python scripts/hg_r2_stage2.py > logs/r2_stage2.log 2>&1 || exit 1
python scripts/hg_r2_s1.py > logs/r2_s1.log 2>&1 || exit 1
touch R2_DONE
echo "$(TZ=America/Los_Angeles date) queue2 done"
