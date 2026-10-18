#!/bin/bash
# S10 tide / sub-daily follow-up queue (rules: svi_science/results/s10_dvv/tide_subday/PLAN.md). DONE markers in that folder.
cd /home/jovyan/tremorferometry/svi_science/scripts
OUT=/home/jovyan/tremorferometry/svi_science/results/s10_dvv/tide_subday
[ -f $OUT/QUANT_DONE ] || python s10_tide_quantities.py > $OUT/logs/quant.out 2>&1
[ -f $OUT/PREP_DONE ] || python s10_tide_prep.py --workers 6 > $OUT/logs/prep.out 2>&1
[ -f $OUT/PHASE_DONE ] || python s10_tide_phase.py --workers 6 --fams 10,87,17 > $OUT/logs/phase.out 2>&1
[ -f $OUT/SUBDAY_DONE ] || { [ -f s10_tide_subday.py ] && python s10_tide_subday.py > $OUT/logs/subday.out 2>&1; }
[ -f $OUT/SUMMARY_DONE ] || { [ -f s10_tide_phase_summary.py ] && python s10_tide_phase_summary.py > $OUT/logs/summary.out 2>&1; }
[ -f $OUT/FIGURES_DONE ] || { [ -f s10_tide_figures.py ] && python s10_tide_figures.py > $OUT/logs/figures.out 2>&1; }
date > $OUT/QUEUE_DONE
