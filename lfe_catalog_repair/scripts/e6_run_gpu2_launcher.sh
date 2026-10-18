#!/bin/bash
# waits for GPU chain 1 to finish, then runs chain 2 with variant B (adopted 2026-09-11 13:35 UTC: higher own recall than A at the
# smallest T meeting the quiet + post-SSE FA bars, 8.5 MAD: B 0.145 vs A 0.108; the ETS site-rolled <= 5 % bar is unmet by both)
cd /home/jovyan/tremorferometry/lfe_catalog_repair
until grep -q "CHAIN 1 DONE" logs/e6_gpu_chain.log; do sleep 15; done
bash scripts/e6_run_gpu2.sh B
