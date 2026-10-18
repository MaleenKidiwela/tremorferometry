#!/bin/bash
# VLFE scan (PLAN_SCAN.md §1): ComCat event METADATA only (no waveforms). Global M>=4.5 2003-2026 per year;
# regional M>=2.5 within 300 km and M>=4.0 within 1500 km of the SVI centroid (48.37 N, -123.61), 2003-2026.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/data/eq
U="https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&orderby=time-asc"
for y in $(seq 2003 2026); do
  f=comcat_global_M45_$y.csv
  [ -s $f ] || curl -sS --retry 3 -o $f "$U&starttime=$y-01-01&endtime=$((y+1))-01-01&minmagnitude=4.5"
done
[ -s comcat_reg300_M25_2003_2026.csv ] || curl -sS --retry 3 -o comcat_reg300_M25_2003_2026.csv "$U&starttime=2003-01-01&endtime=2027-01-01&minmagnitude=2.5&latitude=48.37&longitude=-123.61&maxradiuskm=300"
[ -s comcat_reg1500_M40_2003_2026.csv ] || curl -sS --retry 3 -o comcat_reg1500_M40_2003_2026.csv "$U&starttime=2003-01-01&endtime=2027-01-01&minmagnitude=4.0&latitude=48.37&longitude=-123.61&maxradiuskm=1500"
wc -l *.csv | tail -3
