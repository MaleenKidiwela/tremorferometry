#!/bin/bash
# Stage-4 network download + PGC 2017 BH era.
# WINDOW NOTE: 9 stations x 2 years would be ~146 GB and 9 x 120 d ~29 GB would put the study
# folder at ~89 GB, over the ~80 GB cap. The window is therefore the 60 days spanning the
# 2017-02/03/04 northern ETS (PNSN tremor within 60 km of PGC: Feb 197, Mar 3,511, Apr 326),
# which is what the network-consistent template stacks need. 9 stations kept (azimuthal
# coverage matters more for LOCATION than record length).
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
DL="$PY lfe_catalog_repair/scripts/d3_download.py --outroot lfe_catalog_repair/traces --workers 3"
while pgrep -f "d3_download.py --outroot" | grep -qv $$ && pgrep -f "sta PGC" >/dev/null; do sleep 30; done
$DL --net CN --sta PGC --chan 'BH?' --start 2017-01-01 --end 2017-08-10
S=2017-02-20; E=2017-04-21
$DL --net CN --sta SYMB --chan 'HH?' --start $S --end $E
$DL --net PB --sta B011 --chan 'EH?' --start $S --end $E
$DL --net CN --sta VGZ  --chan 'HH?' --start $S --end $E
$DL --net UW --sta LRIV --chan 'BH?' --start $S --end $E
$DL --net CN --sta CLRS --chan 'HH?' --start $S --end $E
$DL --net CN --sta PTRF --chan 'HH?' --start $S --end $E
$DL --net CN --sta SNB  --chan 'HH?' --start $S --end $E
$DL --net CN --sta GOBB --chan 'HH?' --start $S --end $E
$DL --net UW --sta OCP  --chan 'EH?' --start $S --end $E
echo "[stage4-dl ALL DONE]"
du -sh lfe_catalog_repair/traces
