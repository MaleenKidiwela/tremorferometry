#!/bin/bash
while true; do
  anon=$(awk '/^anon /{print int($2/1073741824)}' /sys/fs/cgroup/memory.stat)
  echo "$(date +%H:%M:%S) anon=${anon}GB/187"
  if [ "$anon" -gt 155 ]; then echo "!!! OOM DANGER anon=${anon}GB — main session should intervene"; break; fi
  sleep 30
done
