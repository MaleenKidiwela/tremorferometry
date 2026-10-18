#!/bin/bash
# Track D, D1(d) -- the WITHHELD-STATION network on production months, so the winter and summer
# quiet-day detections can be confirmed (or refuted) by stations that did NOT contribute to them.
# Exactly Track A step 3b's configuration (frozen mergeB, 12-s 2012 templates, 2012 weights/moveout,
# CN.VGZ + CN.SNB + PB.B926 EXCLUDED from the network sum), run month by month on Track C's 20 Hz cache.
# Usage: bash scripts/d1d_run_withheld.sh 2015-02 2014-02 ...
cd /home/jovyan/tremorferometry/lfe_extension
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=3
I=/home/jovyan/tremorferometry/lfe_catalog_repair/results
O=/home/jovyan/tremorferometry/lfe_extension/results
ST=PGC,NLLB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011
NS="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B005=PB.B005,B006=PB.B006,B007=PB.B007,B009=PB.B009,B010=PB.B010,B011=PB.B011,B926=PB.B926"
SITES="PGC:PGC,B009,B010,B011;LRIV:LRIV,B005,B006,B007"
for YM in "$@"; do
  Y=${YM%%-*}
  read D0 D1 <<< $(/opt/conda/bin/python -c "
import pandas as pd,sys
ym='$YM'; p=pd.Period(ym,freq='M')
print(p.start_time.dayofyear, p.end_time.dayofyear)")
  TAG=d1d_wh_${YM}
  echo "[$(date -u +%FT%TZ)] $TAG  year $Y doy $D0..$D1"
  python scripts/e4_template_test.py --year $Y --lengths 12 --n-random 200 --doy-from $D0 --doy-to $D1 \
    --stations $ST --netsta "$NS" --no-seisbench --families-csv $I/e5r_families_2012.csv \
    --stations-csv $I/e5r_stations.csv --gt-parquet $I/e5r_gt_bostock2012.parquet \
    --offsets-csv $I/e5r_family_offsets_2012.csv --extra-tpl $I/e5r_templates_2012.npz --skip-z --half all \
    --moveout-in $I/e4_moveout_e5r_mv_build.csv --budget --budget-thr 8,8.25,8.5,8.75,9 --catalog-out \
    --weights $I/e5r_weights_snr_2012.csv --sites "$SITES" --site-merge B --out-dir $O --tag $TAG \
    --cache-dirs /home/jovyan/tremorferometry/lfe_extension/cache20hz \
    > logs/$TAG.log 2>&1
  echo "$TAG exit=$? $(date -u +%FT%TZ)" >> logs/d1d_withheld.log
done
