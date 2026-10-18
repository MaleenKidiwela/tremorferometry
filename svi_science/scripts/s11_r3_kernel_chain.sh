#!/bin/bash
# S11 round 3 kernel stage (PLAN §6.1/§6.3): SVI main, SVI ext, WA-NW. Skips parts whose marker exists.
cd /home/jovyan/tremorferometry/svi_science/scripts
L=/home/jovyan/tremorferometry/svi_science/results/s11_depth/round3
say(){ echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') $*" >> $L/queue.log; }
[ -f $L/svi/KERNEL_DONE ] || { say "start kernel SVI main"; S11_REGION=SVI python s11_r3_kernel.py --workers 12 >> $L/svi/logs/kernel_stage.out 2>&1; say "end kernel SVI main"; }
[ -f $L/svi/KERNEL_ext_DONE ] || { say "start kernel SVI ext"; S11_REGION=SVI python s11_r3_kernel.py --workers 12 --ext >> $L/svi/logs/kernel_stage.out 2>&1; say "end kernel SVI ext"; }
[ -f $L/wa_nw/KERNEL_DONE ] || { say "start kernel WA_NW"; S11_REGION=WA_NW python s11_r3_kernel.py --workers 12 >> $L/wa_nw/logs/kernel_stage.out 2>&1; say "end kernel WA_NW"; }
