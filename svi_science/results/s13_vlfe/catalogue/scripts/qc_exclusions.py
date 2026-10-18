#!/usr/bin/env python
"""Amendment W2b: station-component-months with median 10-min RMS > 3x the station-component's median of monthly medians
over 2011-2018 -> tables/qc_exclusions.csv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import cat_common as K
D = pd.concat([pd.read_csv(f'{K.TB}/diag_noise_{y}.csv') for y in range(2011, 2019)], ignore_index=True)
ref = D.groupby(['sta', 'comp']).med_rms.median().rename('ref')
if len(sys.argv) > 1:                                     # Amendment X: append the exclusions of a forward year
    y = int(sys.argv[1])
    F = pd.read_csv(f'{K.TB}/diag_noise_{y}.csv')
    own = F.groupby(['sta', 'comp']).med_rms.median().rename('own')
    F = F.join(ref, on=['sta', 'comp']).join(own, on=['sta', 'comp'])
    F['ref'] = F.ref.fillna(F.own)                        # stations without 2011-2018 data: their own year's median
    Xn = F[F.med_rms > 3 * F.ref][['year', 'sta', 'comp', 'month', 'med_rms', 'ref']]
    old = pd.read_csv(f'{K.TB}/qc_exclusions.csv'); old = old[old.year != y]
    pd.concat([old, Xn], ignore_index=True).to_csv(f'{K.TB}/qc_exclusions.csv', index=False)
    print(y, Xn.groupby('sta').size().to_dict()); sys.exit(0)
D = D.join(ref, on=['sta', 'comp'])
X = D[D.med_rms > 3 * D.ref][['year', 'sta', 'comp', 'month', 'med_rms', 'ref']]
X.to_csv(f'{K.TB}/qc_exclusions.csv', index=False)
print(X.groupby(['year', 'sta']).size().unstack(fill_value=0).to_string())
