#!/usr/bin/env bash
# P5 v2 (replaces bw_prod_p5.sh from 2026-10-09).  Same detector runs; era runs now use the era template set WITH KELB
# (templates/kelb/).  Order: (1) p5e_2003_kelb = the 2003 era run again with KELB (all 2003 era days; it REPLACES p5e_2003,
# which is kept on disk); (2) each ETS window as soon as its caches exist (frozen p5f_<Y> + era p5e_<Y> on the window days);
# (3) the rest of each year when P2_<Y>_DONE exists -- the same run names, so window days already done are skipped.
# P3 did not adopt the SNB era template, so the era days are one run per year (no 2006-04-20 split).  <= 3 GPU passes.
set -u -o pipefail
source /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_env.sh
echo $$ > $O/P5_PID
export SVI_REF_DIR=$O/mf
CD=$CB,$FWD
waitm(){ while [ ! -f $O/$1 ] && [ ! -f ${O}/${1/DONE/FAILED} ]; do sleep 300; done; }
guard(){ python -c "import sys; sys.path.insert(0,'$LFE/scripts'); import svi_ref_common as R; R.guard('mf_$1','mf_$1_s*')" >> $O/logs/p5_guard.log 2>&1; }
mf(){ # run sset daysfile extra-env...
  local run=$1 sset=$2 days=$3; shift 3
  [ "$(python -c "import json; print(len(json.load(open('$days'))))")" = 0 ] && return 0
  for i in 0 1 2; do
    env "$@" python $S/bw_gate_mf.py --backward -- --run $run --kind prod --sset $sset --clean 0 --days-file $days \
      --cache-dirs $CD --shard $i --nshard 3 >> $O/logs/${run}_s$i.log 2>&1 &
  done
  wait; guard $run
}
days(){ # out start end mode(all|era)
  python - "$@" <<'PY'
import sys, json, pandas as pd
o, a, b, m = sys.argv[1:5]
d = [str(x.date()) for x in pd.date_range(a, b)]
if m == "era": d = [x for x in d if "2003-02-01" <= x <= "2006-06-30"]
json.dump(d, open(o, "w")); print(o, len(d))
PY
}
ERA="BW_SITELIST_EXT=1 BW_ERA_TPL=$ERA_TPL"
# (1) 2003 era again, with KELB
waitm P2_KELB_DONE
if [ ! -f $O/P5E_2003_KELB_DONE ]; then
  days $O/tables/p5_days_2003_era_kelb.json 2003-01-01 2003-12-31 era >> $O/logs/p5.log
  mf p5e_2003_kelb S15ERA $O/tables/p5_days_2003_era_kelb.json $ERA && date -u +%FT%TZ > $O/P5E_2003_KELB_DONE || date -u +%FT%TZ > $O/P5E_2003_KELB_FAILED
fi
# (2) ETS windows
for w in $WINDOWS; do
  IFS=: read n a b <<< "$w"; [ -f $O/P5W_${n}_DONE ] && continue
  waitm P2_${n}_DONE; ok=1
  for Y in $(seq ${a:0:4} ${b:0:4}); do
    s=$a; [ "$Y" != "${a:0:4}" ] && s=$Y-01-01; e=$b; [ "$Y" != "${b:0:4}" ] && e=$Y-12-31
    days $O/tables/p5_days_${n}_${Y}_frozen.json $s $e all >> $O/logs/p5.log
    days $O/tables/p5_days_${n}_${Y}_era.json $s $e era >> $O/logs/p5.log
    mf p5f_$Y S15 $O/tables/p5_days_${n}_${Y}_frozen.json BW_NOP=1 || ok=0
    mf p5e_$Y S15ERA $O/tables/p5_days_${n}_${Y}_era.json $ERA || ok=0
  done
  [ $ok = 1 ] && date -u +%FT%TZ > $O/P5W_${n}_DONE || date -u +%FT%TZ > $O/P5W_${n}_FAILED
done
# (3) full years
for Y in 2004 2005 2006 2007 2008 2009; do
  waitm P2_${Y}_DONE
  days $O/tables/p5_days_${Y}_frozen.json $Y-01-01 $Y-12-31 all >> $O/logs/p5.log
  days $O/tables/p5_days_${Y}_era.json $Y-01-01 $Y-12-31 era >> $O/logs/p5.log
  [ -f $O/P5F_${Y}_DONE ] || { mf p5f_$Y S15 $O/tables/p5_days_${Y}_frozen.json BW_NOP=1 && date -u +%FT%TZ > $O/P5F_${Y}_DONE || date -u +%FT%TZ > $O/P5F_${Y}_FAILED; }
  if [ $Y -le 2006 ] && [ ! -f $O/P5E_${Y}_DONE ]; then
    mf p5e_$Y S15ERA $O/tables/p5_days_${Y}_era.json $ERA && date -u +%FT%TZ > $O/P5E_${Y}_DONE || date -u +%FT%TZ > $O/P5E_${Y}_FAILED
  fi
done
date -u +%FT%TZ > $O/P5_DONE
