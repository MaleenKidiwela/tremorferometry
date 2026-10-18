#!/usr/bin/env python
"""Flagged production episodes (Amendment U item 3): raw waveform figures (scan near stations + regional + 5 far LHZ, Z,
0.02-0.05 Hz displacement) and battery columns for the tier events of an episode. usage: cat_flag.py YEAR EPISODE
-> tables/flag_<YEAR>_<EPISODE>.csv, figures/flag_<YEAR>_<EPISODE>_*.png"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

NEAR = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV', 'CN.LZB']
REG = ['UW.GNW', 'UW.DOSE', 'UW.SP2', 'TA.D03D', 'TA.B05D', 'UW.WISH', 'US.NLWA']
FAR = ['IU.COR', 'US.NEW', 'US.HAWA', 'US.MSO', 'US.WVOR']


def main(y, ep, nm='prod'):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scan_refs as SR
    import hg_r3 as X
    cc = SR.all_coords()
    SD = f'{K.S13}/data/scan' if y == 2011 else f'{K.S13}/data/scan_{y}'
    E = pd.read_csv(f'{K.TB}/{nm}{y}_events_v2.csv')
    E = (E[E.time_pt.str.startswith(ep)] if ep[:4].isdigit() and len(ep) > 10 else E[E.episode == ep]).sort_values('epoch')
    D = pd.read_parquet(f'{K.DATA}/{nm}{y}/rows_scored_v2.parquet')
    LF = X.lfe_lists()
    y0 = pd.Timestamp(f'{y}-01-01').value / 1e9
    out = []
    for k, e in enumerate(E.itertuples()):
        r = D[(D.source == 'real') & (D.i == e.i)].sort_values('zmad').iloc[-1]
        ts = int(round(e.epoch - y0)); a, n = ts - 300, 800
        tr = []
        for s in NEAR + REG + FAR:
            f = f'{SD}/bp_{y}g/{s}.npy' if s in NEAR else f'{SD}/bp2_ref/{s}.npy'
            if os.path.isfile(f):
                x = np.asarray(np.load(f, mmap_mode='r')[2, a:a + n], float) * 1e9
                if np.isfinite(x).sum() > 400:
                    tr.append(('far' if s in FAR else 'near', s, x))
        amp = {s: np.nanmax(np.abs(x[300:480])) for g, s, x in tr}
        nm = max([amp[s] for g, s, x in tr if g == 'near'] or [np.nan]); fm = max([amp[s] for g, s, x in tr if g == 'far'] or [np.nan])
        n_lfe_svi_1h = int(X.cnt(LF['svi075'], [e.epoch], 3600)[0]); n_pnsn = int(X.cnt(LF['pnsn100_T2'], [e.epoch], 3600)[0])
        near2 = sorted([s for g, s, x in tr if g == 'near' and s in cc.index], key=lambda s: K.dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], r.lat_best, r.lon_best))[:2]
        o = dict(time_pt=e.time_pt, template=e.template, agree=e.agree, zmad=e.zmad, z_best=e.z_best, Mw_scan=e.Mw_scan,
                 ratio_far=r.ratio_far, n_far=r.n_far, ratio_GNW=r.get('ratio_GNW', np.nan), vr_thrust=r.vr_thrust, z_loso=r.z_loso,
                 pass_a=r.pass_a, pass_b=r.pass_b, pass_c=r.pass_c, pass_e2=r.pass_e2, pass_f=r.pass_f, pass_S5=r.pass_S5,
                 lat=r.lat_best, lon=r.lon_best, nearest_family=r.get('nearest_family', ''), lfe_status=e.lfe_status,
                 svi_lfe075_pm1h=n_lfe_svi_1h, pnsn_tremor_pm1h=n_pnsn, near_max_nm=nm, far_max_nm=fm, far_over_near=fm / nm,
                 head2_maxmag=r.head2_maxmag, in_head2=r.in_head2, near2=','.join(near2),
                 far_over_near2=fm / max([amp[s] for s in near2] or [np.nan]))
        out.append(o)
        fig, ax = plt.subplots(figsize=(10, 0.42 * len(tr) + 1.6)); t = np.arange(n) - 300
        for j, (g, s, x) in enumerate(tr):
            pk = np.nanmax(np.abs(x)); ax.plot(t, x / (2.2 * pk) - j, color='k' if g == 'near' else 'tab:red', lw=0.7)
            try:
                dk = K.dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], r.lat_best, r.lon_best)
            except KeyError:
                dk = np.nan
            ax.text(402, -j, f'{s} {dk:.0f} km peak {pk:.0f} nm', fontsize=7)
        ax.axvspan(0, 180, color='0.9', zorder=0); ax.set_yticks([]); ax.set_xlim(-300, 600); ax.set_xlabel('s after t0')
        ax.set_title(f'{e.time_pt}  {e.template}+[{e.agree}]  MF {e.zmad:.1f}  scan z {e.z_best:.1f}  vr_thrust {r.vr_thrust:.2f}  '
                     f'far/near {fm / nm:.2f}  LFE ±10 min: {e.lfe_status}\nZ, 0.02–0.05 Hz displacement, each trace scaled to its own peak',
                     fontsize=8)
        fig.tight_layout(); fig.savefig(f'{K.FIG}/flag_{y}_{ep.replace(" ", "_").replace(":", "")}_{k:02d}.png', dpi=110); plt.close(fig)
    T = pd.DataFrame(out); T.to_csv(f'{K.TB}/flag_{y}_{ep.replace(" ", "_").replace(":", "")}.csv', index=False)
    pd.set_option('display.width', 250); print(T.round(2).to_string())


if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else 'prod')
