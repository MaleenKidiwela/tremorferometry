#!/usr/bin/env bash
# P4 envelope-CC quiet flag queue (CPU only, 8 workers): each year's window table as soon as its 20 Hz caches exist
# (P2_<Y>_DONE), then the flag file is rewritten with the frozen (C, k); the one-time 2012-13 validation runs inside
# bw_p4env_flag.py at the first run where P2_2012/2013 and P5F_2012/2013 all exist.  Rules: p4/envflag/ENVFLAG.md.
# Markers in p4/envflag/: P4ENV_<Y>_DONE | P4ENV_<Y>_FAILED, P4ENV_VAL_DONE, P4ENV_ALL_DONE.  Nothing deleted.
set -u -o pipefail
source /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
E=$O/p4/envflag
echo $$ > $E/Q_ENV_PID
L=$E/logs/queue.log
ts(){ TZ=America/Los_Angeles date "+%a %F %H:%M %Z"; }
w2(){ while [ ! -f $O/P2_$1_DONE ] && [ ! -f $O/P2_$1_FAILED ]; do sleep 600; done; }
for Y in 2012 2013 2003 2004 2005 2006 2007 2008 2009; do
  [ -f $E/P4ENV_${Y}_DONE ] && continue
  w2 $Y
  [ -f $O/P2_${Y}_FAILED ] && echo "[envq] P2 $Y FAILED marker present; computing on whatever is cached $(ts)" >> $L
  echo "[envq] detect $Y start $(ts)" >> $L
  if python $S/bw_p4env_detect.py --years $Y --workers 8 >> $E/logs/detect_$Y.log 2>&1 \
     && python $S/bw_p4env_flag.py >> $L 2>&1; then
    ts > $E/P4ENV_${Y}_DONE
  else
    ts > $E/P4ENV_${Y}_FAILED; echo "[envq] $Y FAILED $(ts)" >> $L
  fi
  ls $E/tables/envflag_validation_CN4.json >/dev/null 2>&1 && [ ! -f $E/P4ENV_VAL_DONE ] && ts > $E/P4ENV_VAL_DONE
done
# validation, if the P5F markers came after the 2012/2013 envelopes
while [ ! -f $E/P4ENV_VAL_DONE ]; do
  if [ -f $O/P5F_2012_DONE ] && [ -f $O/P5F_2013_DONE ]; then
    python $S/bw_p4env_flag.py >> $L 2>&1
    ls $E/tables/envflag_validation_CN4.json >/dev/null 2>&1 && ts > $E/P4ENV_VAL_DONE || { echo "[envq] validation not written $(ts)" >> $L; break; }
  elif [ -f $O/P5F_2012_FAILED ] || [ -f $O/P5F_2013_FAILED ]; then
    echo "[envq] P5F 2012/2013 FAILED -> no sigmaP validation $(ts)" >> $L; break
  else
    sleep 900
  fi
done
ts > $E/P4ENV_ALL_DONE
echo "[envq] all done $(ts)" >> $L
