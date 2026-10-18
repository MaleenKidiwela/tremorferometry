#!/bin/bash
# S11 tomo/eq analysis chain (tomo/PLAN.md §6.4-6.5): waits for data/eq_picks/PARSE_DONE, then earthquake κ per
# zone -> differential layer κ -> figures. Markers in tomo/eq/: KAPPA_DONE, DIFF_DONE, ANALYSIS_DONE. Idempotent.
R=/home/jovyan/tremorferometry; E=$R/data/eq_picks; SC=$R/svi_science/scripts; Q=$R/svi_science/results/s11_depth/tomo/eq
TS() { TZ=America/Los_Angeles date '+%Y-%m-%d %-I:%M %P PT'; }
echo $$ > $Q/analysis.pid
until [ -f $E/PARSE_DONE ]; do sleep 60; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd $SC
[ -f $Q/KAPPA_DONE ] || { echo "== kappa $(TS)"; python s11_eq_wadati.py --workers 8 > $Q/logs/kappa.log 2>&1 && touch $Q/KAPPA_DONE; }
[ -f $Q/DIFF_DONE ] || { echo "== diff $(TS)"; python s11_eq_diff.py > $Q/logs/diff.log 2>&1 && touch $Q/DIFF_DONE; }
echo "== figures $(TS)"; python s11_eq_figures.py > $Q/logs/figures.log 2>&1 && touch $Q/ANALYSIS_DONE
echo "== end $(TS)"; rm -f $Q/analysis.pid
