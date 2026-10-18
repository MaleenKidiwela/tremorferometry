#!/usr/bin/env python
"""bw_prod_p5_days -- P5 day lists for one year: frozen network every day of the year (2003-01-01..2013-12-31);
era-augmented detector on 2003-02-01..2006-06-30, split at 2006-04-20 (SNB CMG40T epoch end, where an adopted SNB
era-template override stops).  Writes OUT/tables/p5_days_<Y>_{frozen,era_pre,era_post}.json (empty lists allowed)."""
import sys, json
import pandas as pd
Y = int(sys.argv[1]); O = sys.argv[2]
d = [str(x.date()) for x in pd.date_range(f"{Y}-01-01", f"{Y}-12-31")]
era = [x for x in d if "2003-02-01" <= x <= "2006-06-30"]
for k, v in dict(frozen=d, era_pre=[x for x in era if x < "2006-04-20"], era_post=[x for x in era if x >= "2006-04-20"]).items():
    json.dump(v, open(f"{O}/tables/p5_days_{Y}_{k}.json", "w"))
print(Y, len(d), len(era))
