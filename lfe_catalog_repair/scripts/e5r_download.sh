#!/bin/bash
# E5 step 5 STEP A -- THE ONE APPROVED DOWNLOAD (user "Go", 2026-09-11 ~06:05 UTC).
# 15 stations, 2012-08-15 .. 2012-10-05 inclusive (52 days), all three components at the native rate.
# Channel band per station from FDSN station metadata queried 2026-09-11 (e5r ledger): CN PGC/NLLB + UW LRIV = BH? (40 Hz),
# CN VGZ/SNB = HH? (100 Hz; SNB switched BH->HH 2011-03-17, VGZ HH since 2009-07-31), PB B* = EH? (100 Hz).
# Layout traces/<NET>.<STA>/2012/<JJJ>.mseed via scripts/d3_download.py (resumable, .nd no-data markers).
# Concurrency: 3 stations at a time x 2 day-workers each = <= 6 concurrent FDSN requests.  NOTHING ELSE IS FETCHED.
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages
S0=2012-08-15; S1=2012-10-06   # end exclusive -> through 10-05
run() {  # NET STA CHAN
  python scripts/d3_download.py --net $1 --sta $2 --start $S0 --end $S1 --chan "$3" --outroot /home/jovyan/tremorferometry/lfe_catalog_repair/traces --workers 2 \
      > logs/e5_dl_$2.log 2>&1
  echo "[$(date -u +%FT%TZ)] $1.$2 done: $(tail -1 logs/e5_dl_$2.log)"
}
export -f run; export S0 S1
printf "%s\n" "CN PGC BH?" "PB B011 EH?" "PB B009 EH?" "CN VGZ HH?" "CN SNB HH?" "PB B010 EH?" "CN NLLB BH?" "UW LRIV BH?" \
              "PB B926 EH?" "PB B001 EH?" "PB B003 EH?" "PB B004 EH?" "PB B005 EH?" "PB B006 EH?" "PB B007 EH?" \
  | xargs -P 3 -L 1 bash -c 'run $0 $1 $2'
echo "[$(date -u +%FT%TZ)] ALL DONE"
