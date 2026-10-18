#!/bin/bash
# Pod memory watchdog: anon (real, OOM-relevant) from cgroup. WARN >=160G, at >=175G kill largest python; polls every 5 s (2026-10-06: two pod deaths jumped past 185G inside one 60-s poll) (pod cgroup limit 200 GB, 2026-09-28).
LOG=/home/jovyan/tremorferometry/logs/mem_watchdog.log
echo "$(date -Is) watchdog started" >> $LOG
while true; do
  A=$(awk '/^anon /{print $2}' /sys/fs/cgroup/memory.stat)
  G=$((A/1073741824))
  if [ "$G" -ge 175 ]; then
    P=$(ps -eo pid,rss,args --sort=-rss | awk 'NR>1 && /python/ {print $1; exit}')
    echo "$(date -Is) anon=${G}G CRITICAL -> kill -9 $P ($(ps -o args= -p $P 2>/dev/null | cut -c1-140))" >> $LOG
    kill -9 $P 2>/dev/null
  elif [ "$G" -ge 160 ] && [ $(( $(date +%s) - ${LW:-0} )) -ge 60 ]; then
    LW=$(date +%s); echo "$(date -Is) anon=${G}G WARN" >> $LOG
  fi
  sleep 5
done
