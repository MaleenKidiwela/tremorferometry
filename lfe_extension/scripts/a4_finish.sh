#!/bin/bash
# Track A-2: finish the remaining step-4 legs after the first GPU wave (NLLB_old died of a CUDA OOM
# when six matched-filter passes shared the H100).  Runs them ONE AT A TIME.
cd /home/jovyan/tremorferometry/lfe_extension
while pgrep -f "e4_template_test.py --lengths 12" > /dev/null; do sleep 30; done
bash scripts/a4_run_selection.sh g1        # re-runs only the missing a4_NLLB_old
bash scripts/a4_run_stacks.sh NLLB
bash scripts/a4_run_stacks.sh SNB
bash scripts/a4_run_stacks.sh VGZ17
bash scripts/a4_run_stacks.sh B926
bash scripts/a4_run_matched.sh
echo "A4 FINISH DONE $(date -u +%FT%TZ)"
