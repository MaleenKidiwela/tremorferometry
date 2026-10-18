#!/usr/bin/env bash
# G0 pilot (+ the G2 windows): download -> 20 Hz cache -> forward detector with scrambled control on the pilot days.
# Pilot: 2008-05-01..06-04 (whole 2008 ETS) + 2008-08 (quiet), all template stations alive then (B3 = forward S14).
# PB 2008-05 ETS days are already cached on disk (lfe_picker/cache20hz, forward chain) and are not re-downloaded.
# G2 windows: CN4 2004-07-01..28, CN4 + Olympic PB 2005-09-02..26.  Markers OUT/G0_{DONE,FAILED}.  Nothing deleted.
set -u -o pipefail
S=/home/jovyan/tremorferometry/svi_science/scripts
O=/home/jovyan/tremorferometry/svi_science/results/backward_gates
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
rm -f $O/G0_DONE $O/G0_FAILED
fail(){ echo "FAILED at $1 $(date)"; date > $O/G0_FAILED; exit 1; }
CN=PGC,VGZ,NLLB,SNB; PB=B001,B003,B004,B005,B006,B007,B009,B010,B011,B926; OLY=B001,B003,B004,B005,B006,B007
echo "== download $(date)"
python $S/bw_gate_download.py --stations $CN --start 2008-05-01 --end 2008-06-04 --tag pilot_ets || fail dl1
python $S/bw_gate_download.py --stations $CN,$PB --start 2008-08-01 --end 2008-08-31 --tag pilot_quiet || fail dl2
date > $O/G0_DOWNLOAD_PILOT_DONE
python $S/bw_gate_download.py --stations $CN --start 2004-07-01 --end 2004-07-28 --tag g2_2004 || fail dl3
python $S/bw_gate_download.py --stations $CN,$OLY --start 2005-09-02 --end 2005-09-26 --tag g2_2005 || fail dl4
for spec in "2008-05-01 2008-06-04" "2008-08-01 2008-08-31" "2004-07-01 2004-07-28" "2005-09-02 2005-09-26"; do
  set -- $spec; python $S/bw_gate_download.py --stations $CN --start $1 --end $2 --repair --tag repair_$1 || fail repair
done
date > $O/G0_DOWNLOAD_DONE
echo "== cache $(date)"
NS=CN.PGC,CN.VGZ,CN.NLLB,CN.SNB,PB.B001,PB.B003,PB.B004,PB.B005,PB.B006,PB.B007,PB.B009,PB.B010,PB.B011,PB.B926
python $S/bw_gate_cache.py --stations $NS --start 2008-05-01 --end 2008-08-31 --workers 10 || fail cache08
python $S/bw_gate_cache.py --stations $NS --start 2004-07-01 --end 2005-09-26 --workers 10 || fail cache0405
date > $O/G0_CACHE_DONE
echo "== detector $(date)"
python - <<EOF
import json, pandas as pd
d = [str(x.date()) for x in pd.date_range("2008-05-01", "2008-06-04")] + [str(x.date()) for x in pd.date_range("2008-08-01", "2008-08-31")]
json.dump(d, open("$O/tables/pilot_days.json", "w"))
EOF
while [ ! -f $O/G1_RUN_DONE ] && [ ! -f $O/G1_RUN_FAILED ]; do sleep 60; done     # <= 3 GPU passes at once
export SVI_REF_DIR=$O/mf
for i in 0 1; do
  python $S/bw_gate_mf.py --backward -- --run g0_pilot_B3 --kind prod --sset B3 --clean 0 --days-file $O/tables/pilot_days.json \
     --cache-dirs /home/jovyan/tremorferometry/lfe_extension/cache20hz_backward,/home/jovyan/tremorferometry/lfe_picker/cache20hz \
     --shard $i --nshard 2 > $O/logs/g0_pilot_B3_s$i.log 2>&1 &
done
wait
python -c "import sys; sys.path.insert(0,'/home/jovyan/tremorferometry/lfe_extension/scripts'); import svi_ref_common as R; R.guard('mf_g0_pilot_B3','mf_g0_pilot_B3_s*')" \
  > $O/logs/g0_guard.log 2>&1; echo "guard rc=$?" >> $O/logs/g0_guard.log
date > $O/G0_DONE
echo "== done $(date)"
