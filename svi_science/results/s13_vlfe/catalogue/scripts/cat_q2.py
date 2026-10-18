#!/usr/bin/env python
"""Amendment Q2: triage of the 2014 between-ETS tier events (after the Q0 re-summary).
Class (pre-stated): real-looking iff vr_thrust >= 0.12 AND far max |Z| (5 far LHZ, [t0, t0 + 180 s]) < 0.5 x near max
(scan near stations, same window) AND (ratio_GNW <= 0.551 if the best node is north of 48.1 N); else 'likely noise'.
-> tables/q2_triage.csv, figures/q2_wave_*.png"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

NEAR = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV']
FAR = ['IU.COR', 'US.NEW', 'US.HAWA', 'US.MSO', 'US.WVOR']
SD = f'{K.S13}/data/scan_2014'
Y0 = pd.Timestamp('2014-01-01').value / 1e9


def arr(path):
    return np.load(path, mmap_mode='r') if os.path.isfile(path) else None


def main():
    D = pd.read_parquet(f'{K.DATA}/p14/p14_rows_scored.parquet')
    ev = D[D.tier & (D.source == 'real') & (D.period == 'between')].sort_values('epoch')
    # one row per event (+-30 s): keep the highest MF z
    ev = ev.assign(grp=(ev.i.diff().fillna(1e9) > 30).cumsum()).sort_values('zmad').groupby('grp').tail(1).sort_values('epoch')
    rows = []
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scan_refs as SR
    cc = SR.all_coords()
    for r in ev.itertuples():
        ts = int(round(r.epoch - Y0))
        a, n = ts - 300, 800
        tr = []
        for s in NEAR:
            X = arr(f'{SD}/bp_2014g/{s}.npy')
            if X is not None:
                tr.append(('near', s, np.asarray(X[2, a:a + n], float) * 1e9))
        for s in ['UW.GNW'] + FAR:
            X = arr(f'{SD}/bp2_ref/{s}.npy')
            if X is not None:
                tr.append(('far' if s in FAR else 'GNW', s, np.asarray(X[2, a:a + n], float) * 1e9))
        w0, w1 = 300, 480
        near_max = max([np.nanmax(np.abs(x[w0:w1])) for g, s, x in tr if g == 'near' and np.isfinite(x[w0:w1]).any()] or [np.nan])
        far_max = max([np.nanmax(np.abs(x[w0:w1])) for g, s, x in tr if g == 'far' and np.isfinite(x[w0:w1]).any()] or [np.nan])
        north = r.lat_best > 48.1
        c1 = r.vr_thrust >= 0.12
        c2 = far_max < 0.5 * near_max
        c3 = (r.ratio_GNW <= 0.551) if north else True
        cls = 'real-looking' if (c1 and c2 and c3) else 'likely noise'
        tpt = pd.Timestamp(r.epoch, unit='s').tz_localize('UTC').tz_convert('America/Los_Angeles').strftime('%Y-%m-%d %H:%M:%S %Z')
        rows.append(dict(time_pt=tpt, template=r.template, zmad=r.zmad, agree=r.agree, z_best=r.z_best, lat_best=r.lat_best,
                         lon_best=r.lon_best, vr_thrust=r.vr_thrust, ratio_far=r.ratio_far, n_far=r.n_far, ratio_GNW=r.ratio_GNW,
                         near_max_nm=near_max, far_max_nm=far_max, far_over_near=far_max / near_max,
                         vr_ok=c1, far_ok=c2, gnw_ok=c3, lfe_status=r.lfe_status, n_lfe_10m=r.n_lfe_10m,
                         head2_maxmag=r.head2_maxmag, cls=cls))
        fig, ax = plt.subplots(figsize=(10, 0.5 * len(tr) + 1.5))
        t = np.arange(n) - 300
        for k, (g, s, x) in enumerate(tr):
            if not np.isfinite(x).any():
                ax.text(-300, -k, f'{s} no data', fontsize=7); continue
            pk = np.nanmax(np.abs(x))
            ax.plot(t, x / (2.2 * pk) - k, color={'near': 'k', 'GNW': 'tab:blue', 'far': 'tab:red'}[g], lw=0.7)
            try:
                dkm = K.dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], r.lat_best, r.lon_best)
            except KeyError:
                dkm = np.nan
            ax.text(402, -k, f'{s} {dkm:.0f} km peak {pk:.0f} nm', fontsize=7)
        ax.axvspan(0, 180, color='0.9', zorder=0); ax.set_xlim(-300, 600); ax.set_yticks([])
        ax.set_xlabel('s after t0'); ax.set_title(f'{tpt}  {r.template}+[{r.agree}]  MF {r.zmad:.1f}  scan z {r.z_best:.1f}  '
                                                  f'vr_thrust {r.vr_thrust:.2f}  far/near {far_max / near_max:.2f}  -> {cls}\n'
                                                  'Z, 0.02-0.05 Hz displacement, each trace scaled to its own peak', fontsize=8)
        fig.tight_layout(); fig.savefig(f'{K.FIG}/q2_wave_{tpt[:19].replace(" ", "_").replace(":", "")}.png', dpi=120); plt.close(fig)
    T = pd.DataFrame(rows); T.to_csv(f'{K.TB}/q2_triage.csv', index=False)
    pd.set_option('display.width', 250); print(T.round(2).to_string())


if __name__ == '__main__':
    main()
