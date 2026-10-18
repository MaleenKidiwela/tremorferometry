#!/bin/bash
# E2 step 3 — five extra 3-C Vancouver Island stations for the 2010 association/location test.
# Window = the 34 densest days of the Lin ETS-2010 catalogue, matching the neighbours already on
# disk (2010-08-10..2010-09-13).  Stations chosen to (a) be 3-component and (b) reproduce Lin's own
# network: LZB, PFB, YOUB, NLLB were 4 of his 14 training stations; C8.SOKB fills the SW azimuth.
# Budget: ~4.3 GB (HH 100 Hz ~32 MB/day, BH 40 Hz ~14 MB/day).
cd /home/jovyan/tremorferometry
DL="python lfe_catalog_repair/scripts/d3_download.py --outroot lfe_catalog_repair/traces --workers 3"
S=2010-08-10; E=2010-09-13
$DL --net CN --sta LZB  --chan 'BH?' --start $S --end $E
$DL --net CN --sta PFB  --chan 'HH?' --start $S --end $E
$DL --net CN --sta YOUB --chan 'HH?' --start $S --end $E
$DL --net CN --sta NLLB --chan 'BH?' --start $S --end $E
$DL --net C8 --sta SOKB --chan 'HH?' --start $S --end $E
echo "[e2-dl-vi2010 DONE]"
du -sh lfe_catalog_repair/traces
