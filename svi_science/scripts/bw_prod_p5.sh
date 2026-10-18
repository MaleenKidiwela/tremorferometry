#!/usr/bin/env bash
# P5 matched filter, every day 2003-01..2013-12, data + scrambled control (forward svi_ref_mf.py via bw_gate_mf.py):
#   frozen network (sset S15: whatever of the 15 frozen-template stations is alive that day) on every day;
#   era-augmented detector (S15 + POLARIS/LZB era templates, P3) on 2003-02-01..2006-06-30 -> both variants on the whole
#   era = the overlap that measures every era-boundary step; SNB era template override before 2006-04-20 only if P3 adopted it.
# Starts each year as soon as P2_<Y>_DONE exists (era years also wait for P3_DONE).  <= 3 GPU passes at once.
# Markers P5F_<Y>_DONE / P5E_<Y>_DONE (or _FAILED, > 3 % station-day errors by the forward guard).  Outputs mf/feat/<run>/.
set -u -o pipefail
source /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_env.sh
echo $$ > $O/P5_PID
export SVI_REF_DIR=$O/mf
CD=$CB,$FWD
guard(){ python -c "import sys; sys.path.insert(0,'$LFE/scripts'); import svi_ref_common as R; R.guard('mf_$1','mf_$1_s*')" >> $O/logs/p5_guard.log 2>&1; }
mf(){ # run sset daysfile extra-env
  local run=$1 sset=$2 days=$3; shift 3
  [ "$(python -c "import json; print(len(json.load(open('$days'))))")" = 0 ] && return 0
  for i in 0 1 2; do
    env "$@" python $S/bw_gate_mf.py --backward -- --run $run --kind prod --sset $sset --clean 0 --days-file $days \
      --cache-dirs $CD --shard $i --nshard 3 > $O/logs/${run}_s$i.log 2>&1 &
  done
  wait; guard $run
}
for Y in 2010 2011 2012 2013 2003 2004 2005 2006 2007 2008 2009; do
  while [ ! -f $O/P2_${Y}_DONE ] && [ ! -f $O/P2_${Y}_FAILED ]; do sleep 300; done
  python $S/bw_prod_p5_days.py $Y $O >> $O/logs/p5.log
  if [ ! -f $O/P5F_${Y}_DONE ]; then
    mf p5f_$Y S15 $O/tables/p5_days_${Y}_frozen.json BW_NOP=1 && date -u +%FT%TZ > $O/P5F_${Y}_DONE || date -u +%FT%TZ > $O/P5F_${Y}_FAILED
  fi
  if [ $Y -le 2006 ] && [ ! -f $O/P5E_${Y}_DONE ]; then
    while [ ! -f $O/P3_DONE ] && [ ! -f $O/P3_FAILED ]; do sleep 300; done
    [ -f $O/P3_DONE ] || { date -u +%FT%TZ > $O/P5E_${Y}_FAILED; continue; }
    OV=$(python -c "import json; print(json.load(open('$O/templates/p3_config.json'))['snb_override'])")
    ok=1
    mf p5e_$Y S15ERA $O/tables/p5_days_${Y}_era_pre.json BW_SITELIST_EXT=1 BW_ERA_TPL=$O/templates/g3_era_templates.npz BW_TPL_OVERRIDE=$OV || ok=0
    mf p5e2_$Y S15ERA $O/tables/p5_days_${Y}_era_post.json BW_SITELIST_EXT=1 BW_ERA_TPL=$O/templates/g3_era_templates.npz || ok=0
    [ $ok = 1 ] && date -u +%FT%TZ > $O/P5E_${Y}_DONE || date -u +%FT%TZ > $O/P5E_${Y}_FAILED
  fi
done
date -u +%FT%TZ > $O/P5_DONE
