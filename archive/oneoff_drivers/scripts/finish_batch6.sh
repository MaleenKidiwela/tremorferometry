#!/bin/bash
# Finish the 6-station PB batch the blocked agent left mid-flight (main-session Python works).
# Resumable: every step skips if its output exists. One GPU job at a time (script is sequential).
cd /home/jovyan/tremorferometry
export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export CUDA_PATH=/opt/conda/targets/x86_64-linux
PY=/home/jovyan/envs/tremorferometry/bin/python
L=logs/finish_batch6.log; : > $L
say(){ echo "=== $(date +%H:%M:%S) $* ===" | tee -a $L; }

# coords for selection (only B933 needs it; others pre-selected)
declare -A LAT=( [B933]=40.0571 ) LON=( [B933]=-123.9728 ) BBOX=( [B933]="39.1561 40.9581 -125.1490 -122.7966" )

mindet_for(){  # echo min-det from median det/fam-day of an mf_all csv
  $PY - "$1" <<'EOF'
import pandas as pd,sys
df=pd.read_csv(sys.argv[1],usecols=['template','time']); df['d']=pd.to_datetime(df.time).dt.date
m=df.groupby(['template','d']).size().median()
print(20 if m>=15 else 8 if m>=6 else 5)
EOF
}

concat(){ s=$1; [ -f data/mf_${s}_all.csv ] && return
  $PY - $s <<'EOF'
import pandas as pd,glob,sys; s=sys.argv[1]
fs=sorted(glob.glob(f'data/mf_{s}_[12]*.csv')); d=pd.concat([pd.read_csv(f) for f in fs],ignore_index=True)
d.to_csv(f'data/mf_{s}_all.csv',index=False); print('concat',s,len(d))
EOF
}

finish_cpu(){  # concat(if needed) -> stack -> all-time dvv -> plots -> 3 rolling windows
  S=$1; s=$(echo $S|tr A-Z a-z); U=$S
  say "$S CPU finish"
  concat $s
  if [ ! -f data/long_window_daily_${U}.npz ]; then
    MD=$(mindet_for data/mf_${s}_all.csv); say "$S min-det $MD"
    $PY scripts/build_long_window_resp.py --mf-csv data/mf_${s}_all.csv --network PB --station $U \
        --no-deconv --min-det $MD --despike-mad 8 --workers 16 --out data/long_window_daily_${U}.npz >> $L 2>&1
  fi
  if [ ! -f data/daily_dvv_${U}_coda_1to4.csv ]; then
    $PY scripts/dvv_coda_parallel.py --npz data/long_window_daily_${U}.npz --window 1.0 4.0 --station $U \
        --workers 24 --out-csv data/daily_dvv_${U}_coda_1to4.csv --out-fig figures/smoke_dvv_${U}_coda_1to4.png >> $L 2>&1
  fi
  $PY scripts/plot_dvv_metadata.py --dvv-csv data/daily_dvv_${U}_coda_1to4.csv --network PB --station $U \
      --out figures/smoke_dvv_${U}_metadata.png >> $L 2>&1 || true
  for W in "1.0 3.0 1to3" "2.0 4.0 2to4" "3.0 5.0 3to5"; do set -- $W
    [ -f data/daily_dvv_${U}_${3}_cal.csv ] || $PY scripts/dvv_roll30cal.py --station $U \
        --npz data/long_window_daily_${U}.npz --window $1 $2 --out data/daily_dvv_${U}_${3}_cal.csv --workers 16 >> $L 2>&1
  done
  say "$S DONE"
}

gpu_densify(){  # densify selected families (GPU) then CPU finish
  S=$1; s=$(echo $S|tr A-Z a-z)
  if [ $(ls data/mf_${s}_[12]*.csv 2>/dev/null|wc -l) -eq 0 ]; then
    say "$S densify (GPU)"
    $PY scripts/densify_launcher.py --templates-npz data/${s}_pnsn_families_100km.npz \
        --summary-csv data/${s}_coverage_selection.summary.csv --min-snr 0 --network PB --station $S \
        --out-prefix mf_${s}_ --workers 20 --top-n 100 --max-raw-det 3000000 --despike-mad 8 >> $L 2>&1
  fi
  finish_cpu $S
}

# ---- B933 full discovery (no families yet) ----
b933_full(){
  S=B933; s=b933
  if [ ! -f data/${s}_pnsn_candidates_100km.parquet ]; then
    say "B933 stage-1"
    $PY scripts/discover_nllb_pnsn_driven.py --station $S --pnsn catalogs/pnsn_tremor_cascadia_full.csv \
        --bbox ${BBOX[$S]} --candidates-out data/${s}_pnsn_candidates_100km.parquet --candidates-only --workers 16 >> $L 2>&1
  fi
  if [ ! -f data/${s}_pnsn_families_100km.npz ]; then
    say "B933 stage-2 (GPU)"
    $PY scripts/discover_gpu.py --station $S --candidates data/${s}_pnsn_candidates_100km.parquet \
        --out data/${s}_pnsn_families_100km.npz --max-bin-candidates 2000 --workers 24 >> $L 2>&1
  fi
  if [ ! -f data/${s}_coverage_selection.summary.csv ]; then
    say "B933 floor-tune + select"
    FLOOR=$($PY - <<'EOF'
import pandas as pd
df=pd.read_csv('data/b933_pnsn_families_100km.summary.csv')
# pick floor so eligible pool ~300-1000 (top ~1%); borehole SNR compressed
for f in [4,5,6,7,8,9,10,12]:
    n=(df.snr>=f).sum()
    if n<=1000: print(f); break
else: print(12)
EOF
)
    say "B933 floor $FLOOR"
    $PY scripts/select_coverage_families.py --summary data/${s}_pnsn_families_100km.summary.csv \
        --station-lat ${LAT[$S]} --station-lon ${LON[$S]} --min-snr $FLOOR --az-sectors 12 \
        --dist-rings "0,20,40,65,100" --k 2 --out data/${s}_coverage_selection.summary.csv \
        --out-fig figures/smoke_${s}_coverage_selection.png >> $L 2>&1
    $PY scripts/plot_family_map.py --station $S --station-lat ${LAT[$S]} --station-lon ${LON[$S]} \
        --summary data/${s}_pnsn_families_100km.summary.csv --min-snr $FLOOR \
        --out figures/smoke_${s}_family_map.png >> $L 2>&1 || true
  fi
  gpu_densify $S
}

say "BATCH FINISH START (GPU free, py OK)"
finish_cpu B932          # stacked already -> just dvv+plots+rolling
finish_cpu B045          # densified -> concat+stack+dvv+plots+rolling
gpu_densify B005
gpu_densify B003
gpu_densify B049
b933_full
# refresh path map with all new stations
$PY scripts/plot_pb_path_map.py >> $L 2>&1 || true
say "ALL SIX DONE"
