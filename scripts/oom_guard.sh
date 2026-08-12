#!/bin/bash
# Persistent OOM guard (pod cgroup cap ~187 GB; `free` lies/shows 1.5 TB). Logs anon continuously and,
# on danger, SIGSTOPs the densify (the big job — per-year atomic checkpoints make it safe to pause) until
# memory drains, then SIGCONTs. Never touches stacks/dvv/finalize/battery/downloads. Replaces the old
# one-shot oom_watch.sh (which `break`ed after the first alert). Usage: nohup bash scripts/oom_guard.sh &
CAP=187; WARN=125; STOP=150; RESUME=110
paused=""
while true; do
  anon=$(awk '/^anon /{print int($2/1073741824)}' /sys/fs/cgroup/memory.stat)
  ts=$(date +%H:%M:%S)
  if [ "$anon" -gt "$STOP" ] && [ -z "$paused" ]; then
    dp=$(ps -eo pid,cmd | awk '/[d]ensify_gnw_gpu/{print $1; exit}')
    if [ -n "$dp" ]; then kill -STOP "$dp" 2>/dev/null; paused="$dp"
      echo "$ts !!! anon=${anon}GB > ${STOP} — SIGSTOP densify pid $dp (will resume < ${RESUME})"
    else echo "$ts !!! anon=${anon}GB > ${STOP} but no densify to pause — MAIN SESSION INTERVENE"; fi
  elif [ -n "$paused" ] && [ "$anon" -lt "$RESUME" ]; then
    kill -CONT "$paused" 2>/dev/null; echo "$ts anon=${anon}GB < ${RESUME} — SIGCONT densify pid $paused"; paused=""
  elif [ "$anon" -gt "$WARN" ]; then echo "$ts anon=${anon}GB/${CAP} WARN"
  else echo "$ts anon=${anon}GB/${CAP}"; fi
  sleep 30
done
