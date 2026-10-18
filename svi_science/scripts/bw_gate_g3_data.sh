#!/usr/bin/env bash
# G3 data: POLARIS + LZB 2004-06-04..06-17 (quiet) and 2004-07-01..07-28 (ETS); CN4 quiet days; repair pass; 20 Hz caches.
set -u -o pipefail
S=/home/jovyan/tremorferometry/svi_science/scripts
O=/home/jovyan/tremorferometry/svi_science/results/backward_gates
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=2
rm -f $O/G3_DATA_DONE $O/G3_DATA_FAILED
fail(){ echo "FAILED $1 $(date)"; date > $O/G3_DATA_FAILED; exit 1; }
ERA=SHVB,MGCB,TWBB,SHDB,KHVB,TWKB,TSJB,KLNB,SILB,SSIB,TWGB,LZB
python $S/bw_gate_download.py --stations $ERA,PGC,VGZ,NLLB,SNB --start 2004-06-04 --end 2004-06-17 --tag g3_quiet || fail dl1
python $S/bw_gate_download.py --stations $ERA --start 2004-07-01 --end 2004-07-28 --tag g3_ets || fail dl2
python $S/bw_gate_download.py --stations $ERA,PGC,VGZ,NLLB,SNB --start 2004-06-04 --end 2004-06-17 --repair --tag repair_g3q || fail rp1
python $S/bw_gate_download.py --stations $ERA --start 2004-07-01 --end 2004-07-28 --repair --tag repair_g3e || fail rp2
NS=C8.SHVB,C8.MGCB,C8.TWBB,C8.SHDB,C8.KHVB,PO.TWKB,PO.TSJB,PO.KLNB,PO.SILB,PO.SSIB,PO.TWGB,CN.LZB,CN.PGC,CN.VGZ,CN.NLLB,CN.SNB
python $S/bw_gate_cache.py --stations $NS --start 2004-06-04 --end 2004-07-28 --workers 10 || fail cache
date > $O/G3_DATA_DONE
