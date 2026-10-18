#!/bin/bash
# S11 tomo/eq earthquake-pick fetch (tomo/PLAN.md §6). Idempotent: cached events are skipped, failures retried.
# Up to 4 passes; then parse (python -I) -> data/eq_picks/picks.parquet; markers in data/eq_picks/.
R=/home/jovyan/tremorferometry; E=$R/data/eq_picks; SC=$R/svi_science/scripts
mkdir -p $E/raw/comcat
echo $$ > $E/fetch.pid
TS() { TZ=America/Los_Angeles date '+%Y-%m-%d %-I:%M %P PT'; }
if [ ! -f $E/FETCH_DONE ]; then
  for pass in 1 2 3 4; do
    echo "== pass $pass $(TS)"
    python -I $SC/s11_eq_fetch.py --workers 2 && { touch $E/FETCH_DONE; break; }
    echo "pass $pass left failures; waiting 5 min"; sleep 300
  done
  [ -f $E/FETCH_DONE ] || { ls $E/raw/comcat/*.fail 2>/dev/null | wc -l > $E/FETCH_PARTIAL; echo "fetch PARTIAL $(TS)"; }
fi
echo "== parse $(TS)"
python -I $SC/s11_eq_parse.py && touch $E/PARSE_DONE
echo "== end $(TS)"
rm -f $E/fetch.pid
