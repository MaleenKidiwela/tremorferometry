#!/bin/bash
# Stage-5 Lin leg: the 2010 network around PGC over the DENSEST 34 days of the Lin ETS-2010
# catalogue (2010-08-10..2010-09-13 holds 80,341 of its 106,644 events). 5 three-component
# neighbours + PGC = 18 channels, so the network-sum cc noise should fall to ~0.17/sqrt(18) ~ 0.04.
# Budget: ~3.8 GB, paid for by deleting the unused PGC 2015-11..2016-04 window (2.7 GB).
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
DL="$PY lfe_catalog_repair/scripts/d3_download.py --outroot lfe_catalog_repair/traces --workers 3"
S=2010-08-10; E=2010-09-13
$DL --net CN --sta VGZ  --chan 'HH?' --start $S --end $E
$DL --net PB --sta B011 --chan 'EH?' --start $S --end $E
$DL --net PB --sta B926 --chan 'EH?' --start $S --end $E
$DL --net CN --sta SNB  --chan 'BH?' --start $S --end $E
$DL --net UW --sta LRIV --chan 'BH?' --start $S --end $E
echo "[stage5-dl2010 ALL DONE]"
du -sh lfe_catalog_repair/traces

# Z-only neighbours, added for azimuthal coverage in the 2010 location work (~1 GB).
S=2010-08-10; E=2010-09-13
$DL --net UW --sta OCP  --chan 'EH?' --start $S --end $E
$DL --net CN --sta GOBB --chan 'EH?' --start $S --end $E
$DL --net UW --sta STW  --chan 'EH?' --start $S --end $E
echo "[stage5-dl2010 ZONLY DONE]"
