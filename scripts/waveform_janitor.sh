#!/bin/bash
# Backstop for the storage rule: delete a NEW station's raw waveforms as soon as its daily stacks
# (.npz) exist AND no pipeline job is using that station. Never touches B018 (user wants it kept) or
# the partial-download sample stations. Emits one stdout line per cleanup (Monitor turns it into an event).
DONE_KEEP="B018"
PARTIALS="UW.ERW UW.LRIV CN.LZB CN.SNB"
WF=/home/jovyan/tremorferometry/data/waveforms
DATA=/home/jovyan/tremorferometry/data
while true; do
  for d in "$WF"/PB.B* "$WF"/UW.* "$WF"/CN.* "$WF"/CC.* "$WF"/NC.* "$WF"/UO.* "$WF"/BK.*; do
    [ -d "$d" ] || continue
    base=$(basename "$d")                 # e.g. PB.B017
    sta=${base#*.}                        # e.g. B017
    [ "$sta" = "$DONE_KEEP" ] && continue
    case " $PARTIALS " in *" $base "*) continue;; esac
    npz="$DATA/long_window_daily_$sta.npz"
    [ -f "$npz" ] || continue             # stacks not built yet -> waveforms still needed
    # require the stacks file to be settled (>30 min old) and no active job referencing this station
    if [ -n "$(find "$npz" -mmin +30 2>/dev/null)" ] && \
       ! pgrep -f "densify.*$sta|build_long.*$sta|discover.*$sta|download_station.*$sta|dvv_.*$sta" >/dev/null 2>&1; then
      sz=$(du -sh "$d" 2>/dev/null | cut -f1)
      rm -rf "$d" && echo "janitor: deleted $base ($sz) — stacks $npz present, station idle"
    fi
  done
  sleep 300
done
