#!/bin/bash
# Track A-2: remaining step-4 legs (the first a4_finish.sh self-matched its own pgrep guard and hung;
# no run was affected, only the wait).  NLLB_old is re-run here because it died of a CUDA OOM when six
# matched-filter passes shared the H100.
cd /home/jovyan/tremorferometry/lfe_extension
bash scripts/a4_run_selection.sh g1
for k in NLLB SNB VGZ17 B926; do bash scripts/a4_run_stacks.sh $k; done
bash scripts/a4_run_matched.sh
echo "A4 FINISH2 DONE $(date -u +%FT%TZ)"
