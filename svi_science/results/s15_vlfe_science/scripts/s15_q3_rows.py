"""Q3 step 0: filtered S13 row files (pooled unit only) for 2012-08, 2018-06, 2024-08: 'ex' = LFEs farther than 60 s from
every catalogued VLFE; 'in' = the complement (within 60 s). Same rows otherwise."""
import numpy as np, pandas as pd
import s15_common as C
V = C.vlfes(); ev = np.sort(V.epoch.values.astype(float))
out = []
for e in ('2012-08', '2018-06', '2024-08'):
    R = pd.read_parquet(f'{C.RES}/s13_vlfe/data/rows/svi_{e}.parquet')
    P = R[R.unit == 'pooled'].copy()
    t = P['to'].values.astype(float)
    j = np.searchsorted(ev, t)
    d = np.minimum(np.abs(t - ev[np.clip(j - 1, 0, len(ev) - 1)]), np.abs(ev[np.clip(j, 0, len(ev) - 1)] - t))
    near = d <= 60
    P[~near].to_parquet(f'{C.S15}/data/rows/s15ex_svi_{e}.parquet', index=False)
    P[near].to_parquet(f'{C.S15}/data/rows/s15in_svi_{e}.parquet', index=False)
    out.append(dict(episode=e, n_pooled=len(P), n_ex=int((~near).sum()), n_in=int(near.sum()),
                    n_vlfe=int((V.episode == e).sum())))
pd.DataFrame(out).to_csv(f'{C.OUT_T}/q3_restack_rows.csv', index=False)
print(pd.DataFrame(out))
