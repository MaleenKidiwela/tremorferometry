#!/bin/bash
# Pod memory watchdog: anon (real, OOM-relevant) from cgroup. WARN >=90G, at >=118G kill largest python.
LOG=/home/jovyan/tremorferometry/logs/mem_watchdog.log
echo "$(date -Is) watchdog started" >> $LOG
while true; do
  A=$(awk '/^anon /{print $2}' /sys/fs/cgroup/memory.stat)
  G=$((A/1073741824))
  if [ "$G" -ge 118 ]; then
    P=$(ps -eo pid,rss,args --sort=-rss | awk 'NR>1 && /python/ {print $1; exit}')
    echo "$(date -Is) anon=${G}G CRITICAL -> kill -9 $P ($(ps -o args= -p $P 2>/dev/null | cut -c1-140))" >> $LOG
    kill -9 $P 2>/dev/null
  elif [ "$G" -ge 90 ]; then
    echo "$(date -Is) anon=${G}G WARN" >> $LOG
  fi
  sleep 60
done
