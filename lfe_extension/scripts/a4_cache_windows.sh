#!/bin/bash
# Track A-2 step 4: 20 Hz / 1-8 Hz caches for the epoch-test windows (frozen chain, e4_cache20 recipe).
#   2018 doy 161-206 (SSE 2018-06-10..07-25)  -> PGC / NLLB / SNB / VGZ boundary tests
#   2021 doy 274-309 (SSE 2021-10-01..11-05)  -> LRIV boundary test
#   2017 doy 70-100  (SSE 2017-03-18..04-03)  -> VGZ secondary test (traces live in lfe_catalog_repair)
cd /home/jovyan/tremorferometry/lfe_extension
O=/home/jovyan/tremorferometry/lfe_extension/cache20hz
EXT=/home/jovyan/tremorferometry/lfe_extension
REP=/home/jovyan/tremorferometry/lfe_catalog_repair
S18=CN.PGC,CN.VGZ,CN.NLLB,CN.SNB,UW.LRIV,PB.B001,PB.B003,PB.B004,PB.B005,PB.B006,PB.B007,PB.B009,PB.B010,PB.B011,PB.B926
S21=CN.PGC,CN.VGZ,CN.NLLB,UW.LRIV,PB.B001,PB.B003,PB.B004,PB.B005,PB.B006,PB.B007,PB.B009,PB.B010,PB.B011,PB.B926
S17=CN.PGC,CN.VGZ,CN.SNB,UW.LRIV,PB.B011
python scripts/a4_cache20.py --stations $S18 --year 2018 --doy-from 161 --doy-to 206 --outdir $O --traces-root $EXT --workers 10 --skip-existing
python scripts/a4_cache20.py --stations $S21 --year 2021 --doy-from 274 --doy-to 309 --outdir $O --traces-root $EXT --workers 10 --skip-existing
python scripts/a4_cache20.py --stations $S17 --year 2017 --doy-from  70 --doy-to 100 --outdir $O --traces-root $REP --workers 10 --skip-existing
echo "CACHE DONE $(date -u +%FT%TZ)"
