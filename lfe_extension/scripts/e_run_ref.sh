#!/bin/bash
# Track E — the per-site pass on the two REFEREE years (labels live here).
# 2012 = training year (Track A frozen run a1_mergeB_L12, 15 stations / 9 sites)
# 2010 = held-out validation year (frozen transfer a2_xfer2010, 9 stations / 6 sites)
# Both also generate the SITE-ROLLED chance control (negative class / chance-coincidence reference).
cd /home/jovyan/tremorferometry/lfe_extension
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages OMP_NUM_THREADS=4
R=/home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz
ST12=PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
NS12="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B005=PB.B005,B006=PB.B006,B007=PB.B007,B009=PB.B009,B010=PB.B010,B011=PB.B011,B926=PB.B926"
ST10=PGC,VGZ,SNB,NLLB,B926,LRIV,B009,B010,B011
NS10="PGC=CN.PGC,VGZ=CN.VGZ,SNB=CN.SNB,NLLB=CN.NLLB,B926=PB.B926,LRIV=UW.LRIV,B009=PB.B009,B010=PB.B010,B011=PB.B011"
python scripts/e_persite.py --mode ref --year 2012 --doy-from 228 --doy-to 279 \
  --det results/a1_hf_ratio_a1_mergeB_L12.parquet --stations $ST12 --netsta "$NS12" \
  --cache-dirs $R --control --tag ref2012 > logs/e_ref2012.log 2>&1
echo "2012 exit=$? $(date -u +%FT%TZ)"
python scripts/e_persite.py --mode ref --year 2010 --doy-from 222 --doy-to 256 \
  --det results/a1_hf_ratio_a2_xfer2010.parquet --stations $ST10 --netsta "$NS10" \
  --cache-dirs $R --control --tag ref2010 > logs/e_ref2010.log 2>&1
echo "2010 exit=$? $(date -u +%FT%TZ)"
