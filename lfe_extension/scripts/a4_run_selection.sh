#!/bin/bash
# Track A-2 STEP 4 -- WITHHELD-STATION event selection for the epoch-identity test (Merlin step 4).
#
# For every (station, epoch boundary) the reference events are selected by the network with THAT
# STATION EXCLUDED from the sum, at the frozen T_core = 8.75 MAD, frozen config (mergeB, 12-s
# templates, 2012 weights/moveout/templates, two-site coincidence gate), on the tremor days of the
# SSE window.  Each boundary gets TWO passes with the SAME withheld station:
#   *_new : the first cached SSE of the NEW epoch  (2018 doy 171-196 / 2021 doy 276-302 / 2017 doy 77-93)
#   *_old : the 2012 ETS tremor days (doy 247-270), which lie in the PREVIOUS epoch for all five
#           boundaries -- this is the CONTROL that measures how much cc / pick shift the
#           "select-with-the-station-withheld, stack at detection times" procedure costs by itself.
# Nothing here is scored against Bostock ground truth (there is none after 2013); the GT parquet is
# passed only because e4_template_test.py needs one, and its match_* columns are ignored.
set -u
cd /home/jovyan/tremorferometry/lfe_extension
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=4
I=/home/jovyan/tremorferometry/lfe_catalog_repair/results
O=/home/jovyan/tremorferometry/lfe_extension/results
C=/home/jovyan/tremorferometry/lfe_extension/cache20hz
NS="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B005=PB.B005,B006=PB.B006,B007=PB.B007,B009=PB.B009,B010=PB.B010,B011=PB.B011,B926=PB.B926"
SITES="PGC:PGC,B009,B010,B011;LRIV:LRIV,B005,B006,B007"
BASE="--lengths 12 --n-random 50 --netsta $NS --no-seisbench --families-csv $I/e5r_families_2012.csv \
  --stations-csv $I/e5r_stations.csv --gt-parquet $I/e5r_gt_bostock2012.parquet \
  --offsets-csv $I/e5r_family_offsets_2012.csv --extra-tpl $I/e5r_templates_2012.npz --skip-z --half all \
  --moveout-in $I/e4_moveout_e5r_mv_build.csv --weights $I/e5r_weights_snr_2012.csv --sites $SITES \
  --site-merge B --budget --budget-thr 8.75 --catalog-out --amp-ref-csv $I/e5r_amp_ref_2012.csv \
  --coinc-gate 2,2.5 --cache-dirs $C --out-dir $O"

run () {   # run <tag> <year> <doyfrom> <doyto> <stations>
  local tag=$1 yr=$2 d0=$3 d1=$4 st=$5
  if [ -f "$O/e5r_detections_$tag.parquet" ]; then echo "[skip] $tag exists"; return; fi
  echo "[$(date -u +%FT%TZ)] $tag  year $yr doy $d0-$d1  stations $st"
  /opt/conda/bin/python scripts/e4_template_test.py $BASE --year $yr --doy-from $d0 --doy-to $d1 \
      --stations "$st" --tag "$tag" > logs/$tag.log 2>&1
  echo "exit=$? $(date -u +%FT%TZ)" >> logs/$tag.log
}

ALL=PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
A21=PGC,VGZ,NLLB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926   # SNB dead 2021-03-28
A17=PGC,VGZ,SNB,LRIV,B011                                                # only these are on disk for 2017

case "${1:-all}" in
 g1)
  run a4_PGC_new  2018 171 196 VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
  run a4_PGC_old  2012 247 270 VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
  run a4_NLLB_new 2018 171 196 PGC,VGZ,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
  run a4_NLLB_old 2012 247 270 PGC,VGZ,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
  ;;
 g2)
  run a4_VGZ_new  2018 171 196 PGC,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
  run a4_VGZ_old  2012 247 270 PGC,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
  run a4_SNB_new  2018 171 196 PGC,VGZ,NLLB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
  run a4_SNB_old  2012 247 270 PGC,VGZ,NLLB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
  ;;
 g3)
  run a4_LRIV_new 2021 276 302 PGC,VGZ,NLLB,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
  run a4_LRIV_old 2012 247 270 PGC,VGZ,NLLB,SNB,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
  run a4_VGZ17_new 2017 77 93  PGC,SNB,LRIV,B011
  run a4_VGZ17_old 2012 247 270 PGC,SNB,LRIV,B011
  ;;
 g4)
  # NEGATIVE CONTROL: B926 has NO epoch boundary (PB EH 100 Hz throughout), so the same withheld-station
  # procedure across the same 2012 -> 2018 gap measures what the YEAR costs, with the instrument fixed.
  run a4_B926_new 2018 171 196 PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011
  run a4_B926_old 2012 247 270 PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011
  ;;
 *) echo "usage: a4_run_selection.sh g1|g2|g3|g4" ;;
esac
echo "GROUP ${1:-} DONE $(date -u +%FT%TZ)"
