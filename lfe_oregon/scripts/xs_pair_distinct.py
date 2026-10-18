"""Are co-detecting Oregon families duplicates or neighbours? (2026-09-28)
For every pair of the 53 families (>=100 detections) sharing >10 % of detections within 5 s (net of a shifted control),
compare their final templates on the stations both have:
  - per-station lag of the max cross-correlation (horizontal components, 2-8 Hz), and its cc;
  - lag spread = max - min of the per-station lags. A duplicate (same source) needs only one common lag on every
    station (spread ~0); a neighbouring source moves the S arrivals differently at different stations.
Output results/chain_paper/scan46/pair_distinct.csv + a summary.
"""
import sys, numpy as np, pandas as pd
from scipy.signal import butter, sosfiltfilt, correlate
sys.path.insert(0, 'scripts'); import xs_scan46 as X
R = 'results/chain_paper'
it = X.last_iter(f'{R}/grow'); T = X.load_tpl(f'{R}/grow/tpl_it{it}.npz')
FS = 40.0; sos = butter(4, [2 / 20, 8 / 20], 'band', output='sos'); MAXLAG = int(3 * FS)
pairs = pd.read_csv(f'{R}/scan46/dup_pairs_shiftctrl.csv'); pairs = pairs[pairs.excess > 0.10]
def lagcc(x, y):
    x = sosfiltfilt(sos, x); y = sosfiltfilt(sos, y)
    c = correlate(x, y, 'full') / (np.linalg.norm(x) * np.linalg.norm(y) + 1e-12)
    mid = len(y) - 1; w = c[mid - MAXLAG: mid + MAXLAG + 1]; k = np.argmax(w)
    return (k - MAXLAG) / FS, w[k]
rows = []
for a, b in pairs[['a', 'b']].values:
    common = sorted(set(T[a]) & set(T[b]))
    lags, ccs = [], []
    for s in common:
        la, ca = [], []
        for comp in (0, 1):                         # N, E
            l, c = lagcc(T[a][s][0][comp], T[b][s][0][comp]); la.append(l); ca.append(c)
        k = int(np.argmax(ca)); lags.append(la[k]); ccs.append(ca[k])
    lags, ccs = np.array(lags), np.array(ccs)
    good = ccs >= 0.5
    spread = float(np.ptp(lags[good])) if good.sum() >= 2 else np.nan
    rows.append(dict(a=a, b=b, n_common=len(common), n_good=int(good.sum()), cc_med=float(np.median(ccs)) if len(ccs) else np.nan,
                     lag_spread_s=spread))
d = pd.DataFrame(rows).merge(pairs, on=['a', 'b'])
d['verdict'] = np.select([(d.n_good >= 2) & (d.lag_spread_s <= 0.15),
                          (d.n_good >= 2) & (d.lag_spread_s > 0.30)], ['duplicate', 'neighbour'], 'unclear')
d.to_csv(f'{R}/scan46/pair_distinct.csv', index=False)
print(d.verdict.value_counts().to_dict()); print(d[['cc_med', 'lag_spread_s', 'n_good']].describe().round(2).to_string())
import networkx as nx
fam = sorted(set(pairs.a) | set(pairs.b) | set(pd.read_csv(f'{R}/scan46/dup_pairs_shiftctrl.csv').a))
g = nx.Graph(); g.add_nodes_from(set(pd.read_csv(f'{R}/scan46/dup_pairs_shiftctrl.csv')[['a','b']].values.ravel()))
g.add_edges_from(d[d.verdict == 'duplicate'][['a', 'b']].values.tolist())
print('families after merging only true duplicates (connected groups):', nx.number_connected_components(g), 'of', g.number_of_nodes())
