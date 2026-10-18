#!/bin/bash
# S11 pilot queue (PLAN.md D8). Stages skip if their DONE marker exists. Stops on KERNEL_STOP.
cd /home/jovyan/tremorferometry/svi_science/scripts
O=/home/jovyan/tremorferometry/svi_science/results/s11_depth
L=$O/logs/queue.log
say(){ echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') $*" >> $L; }
say "queue start (pid $$)"
while [ ! -f $O/KERNEL_DONE ] && [ ! -f $O/KERNEL_STOP ]; do sleep 20; done
if [ -f $O/KERNEL_STOP ]; then say "KERNEL_STOP: unit test or a_d sign rule failed -> queue halted"; exit 1; fi
run(){ m=$1; shift; if [ -f $O/$m ]; then say "skip $m"; return 0; fi; say "start $*"; "$@" >> $O/logs/queue_stage.out 2>&1; rc=$?; say "end $* rc=$rc"; [ -f $O/$m ] || { say "missing $m -> halt"; exit 1; }; }
run MEASURE_DONE python s11_measure.py --workers 12
python s11_gate_a.py >> $O/logs/gate_a.out 2>&1; say "gate_a part 2 rc=$?"
run SVD_DONE python s11_svdpath.py --workers 12
run INVERT_DONE python s11_invert.py
run INJECT_DONE python s11_inject.py --workers 12
run GATES_DONE python s11_gates.py
run FIGURES_DONE python s11_figures.py
say "QUEUE_DONE"; touch $O/QUEUE_DONE
