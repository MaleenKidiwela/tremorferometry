"""Does any earlier PGC single-station repeating-LFE family (fleet run, data/long_window_daily_PGCbhz_Z.npz) match
SVI family 10's grand stack at PGC? Grand stack per single-station family = n_det-weighted mean of its daily stacks (Z, 2-8 Hz, 40 Hz).
Normalised cross-correlation over all lags vs family 10's P>=0.5 grand stack at CN.PGC epoch E1 (BHZ, same era). User request 2026-10-11."""
import numpy as np, pandas as pd
R = "/home/jovyan/tremorferometry"; PIL = f"{R}/svi_science/results/s10_dvv/pilot"; FAM = 10
d = np.load(f"{R}/data/long_window_daily_PGCbhz_Z.npz", allow_pickle=True)
S = d["stacks"]; pt = d["patches"]; nd = d["n_det"].astype(np.float64); t1 = d["t"]
up, inv = np.unique(pt, return_inverse=True)
G = np.zeros((len(up), S.shape[1])); W = np.zeros(len(up))
ok = np.isfinite(S).all(1)
np.add.at(G, inv[ok], S[ok] * nd[ok, None]); np.add.at(W, inv[ok], nd[ok])
G /= np.maximum(W, 1)[:, None]
# family 10 grand stack at PGC, Z, epoch E1 (BHZ era)
rows = pd.read_parquet(f"{PIL}/tables/rows.parquet"); rows = rows[(rows.fam == FAM) & (rows.series == "data") & (rows.P >= 0.5)]
rs = pd.read_parquet(f"{PIL}/tables/row_station.parquet"); rs = rs[(rs.netsta == "CN.PGC") & (rs.epoch == "E1")]
z = np.load(f"{PIL}/win/CN.PGC.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}; tau = z["tau"]
g = rows[rows.row_id.isin(set(rs.row_id)) & rows.row_id.isin(pos)]
Wf = z["W"][[pos[r] for r in g.row_id], 0, :].astype(np.float64); Wf[~np.isfinite(Wf)] = 0
f10 = (Wf * g.P.to_numpy()[:, None]).sum(0) / g.P.sum()
def nxc(a, b):  # a long (family 10), b short (template): best normalised cc over all full-overlap lags
    b = (b - b.mean()) / (np.linalg.norm(b - b.mean()) + 1e-12); n = len(b); best, lag = 0.0, 0
    for i in range(len(a) - n + 1):
        s = a[i:i + n] - a[i:i + n].mean(); c = float(s @ b / (np.linalg.norm(s) + 1e-12))
        if abs(c) > abs(best): best, lag = c, i
    return best, lag
res = []
for k, name in enumerate(up):
    if W[k] < 50: continue
    c, i = nxc(f10, G[k]); lat, lon = [float(x) for x in name.split("__")[0].split("_")]
    dist = np.hypot((lat - 48.22) * 111.2, (lon + 123.43) * 111.2 * np.cos(np.radians(48.2)))
    res.append(dict(family=name, n_det=int(W[k]), cc=round(c, 3), lag_s=round(tau[i] - t1[0], 2), patch_lat=lat, patch_lon=lon, dist_km=round(dist, 1)))
T = pd.DataFrame(res).sort_values("cc", key=abs, ascending=False)
T.to_csv(f"{PIL}/tables/fam10_vs_pgc_single_station.csv", index=False)
print(f"family 10 grand stack: {len(g)} repeats (PGC E1, Z). single-station families compared: {len(T)}")
print(T.head(12).to_string(index=False))
print("\ncc distribution (|cc|):", np.percentile(np.abs(T.cc), [50, 90, 99]).round(2))
np.savez(f"{PIL}/tables/fam10_vs_pgc_single_station_stacks.npz", f10=f10, tau=tau, G=G, names=up, t1=t1, W=W)
