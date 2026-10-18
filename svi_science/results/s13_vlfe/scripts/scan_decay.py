#!/usr/bin/env python
"""VLFE scan pilot 2 (PLAN_SCAN.md S7.2): far-field decay of the LFE-locked VLF signal with the S13 stacker.
  scan_decay.py rows     -> data/rows_scan2/svi_{2011,2012}-08_pooled.parquet (unit == 'pooled' rows only)
  scan_decay.py summary  -> tables/scan2_decay.csv, figures/scan2_decay_stacks.png, tables/scan2_decay_verdict.json
Stacks (run by scan_pilot2_queue.sh): scan2dec_svi_2011-08 (S13_DISP v3), scan2dec_svi_2012-08 (S13_DISP v2).
Per station, Z: A = sqrt(max(RMS^2 - sigma^2, 0)) over the stacker's fit window; R_obs = A / A_PGC; A <= 2 sigma -> not
detected, R_up = RMS / A_PGC. R_pred = RMS(template Z) / RMS(PGC template Z). a_s = <stk_Z, tpl_Z>/<tpl_Z, tpl_Z>."""
import os, sys, json, pickle
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
from s13_stack import FITM, TREL

EPIS = {'2011-08': 'scan2dec_svi_2011-08', '2012-08': 'scan2dec_svi_2012-08'}
RD = f'{C.DATA}/rows_scan2'


def rows():
    os.makedirs(RD, exist_ok=True)
    for e in EPIS:
        R = pd.read_parquet(f'{C.DATA}/rows/svi_{e}.parquet')
        R = R[R.unit == 'pooled']
        R.to_parquet(f'{RD}/svi_{e}_pooled.parquet')
        print(e, len(R))


def summary():
    import scan_refs as RF
    cc = RF.all_coords()
    out = []
    V = {}
    figdat = {}
    for e, tag in EPIS.items():
        f = f'{C.DATA}/stacks/{tag}.pkl'
        if not os.path.isfile(f):
            print('missing', f); continue
        P = pickle.load(open(f, 'rb'))['units']['pooled']
        U = pd.read_csv(f'{C.TB}/{tag}_units.csv').set_index('unit').loc['pooled']
        src = (U.src_lat, U.src_lon)
        A = {}
        for s, stk in P['stk'].items():
            z = np.asarray(stk[2])[FITM]
            tpl = np.roll(P['tpl'][s], P['lag'], axis=1)[2][FITM]
            sig = float(P['sig'][s][2])
            if not np.isfinite(z).all() or not np.isfinite(sig):
                continue
            rms = float(np.sqrt(np.mean(z ** 2)))
            A[s] = dict(rms=rms, sig=sig, amp=float(np.sqrt(max(rms ** 2 - sig ** 2, 0))), det=bool(np.sqrt(max(rms ** 2 - sig ** 2, 0)) > 2 * sig),
                        tpl_rms=float(np.sqrt(np.mean(tpl ** 2))), a_s=float((z * tpl).sum() / (tpl * tpl).sum()) * 1e-12)
            figdat[(e, s)] = (np.asarray(stk[2]), np.roll(P['tpl'][s], P['lag'], axis=1)[2] * U.thrust_m0)
        if 'CN.PGC' not in A:
            continue
        ap, tp = A['CN.PGC']['amp'], A['CN.PGC']['tpl_rms']
        for s, d in A.items():
            la, lo = cc.loc[s, 'lat'], cc.loc[s, 'lon']
            kx = 111.19 * np.cos(np.radians(src[0]))
            dist = float(np.hypot((lo - src[1]) * kx, (la - src[0]) * 111.19))
            out.append(dict(episode=e, sta=s, dist_from_centroid_km=dist, rms_nm=d['rms'] * 1e9, sigma_nm=d['sig'] * 1e9,
                            amp_nm=d['amp'] * 1e9, detected=d['det'], R_obs=d['amp'] / ap if d['det'] else np.nan,
                            R_up=np.nan if d['det'] else d['rms'] / ap, R_pred=d['tpl_rms'] / tp,
                            station_m0_1e12=d['a_s'], pooled_m0_1e12=U.thrust_m0 / 1e12))
    T = pd.DataFrame(out).sort_values(['episode', 'dist_from_centroid_km'])
    T.to_csv(f'{C.TB}/scan2_decay.csv', index=False)
    print(T.round(3).to_string())
    g = T[T.sta == 'UW.GNW'].set_index('episode')
    rr = {e: float(r.R_obs if np.isfinite(r.R_obs) else r.R_up) for e, r in g.iterrows()}
    V['gnw_ratio'] = rr
    V['gnw_detected'] = {e: bool(r.detected) for e, r in g.iterrows()}
    V['gnw_R_pred'] = {e: float(r.R_pred) for e, r in g.iterrows()}
    V['decay_pass'] = bool(len(rr) == 2 and all(v <= 0.5 for v in rr.values()))
    cfac = [float(r.R_obs / r.R_pred) for _, r in g.iterrows() if np.isfinite(r.R_obs) and r.R_pred > 0]
    V['veto_pred_factor_c'] = 1.0 if V['decay_pass'] else (max(cfac) if cfac else 1.0)
    V['gnw_cannot_separate'] = bool(any(v > 0.8 for v in rr.values()))
    json.dump(V, open(f'{C.TB}/scan2_decay_verdict.json', 'w'), indent=1)
    print(json.dumps(V, indent=1))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    eps = sorted({k[0] for k in figdat})
    fig, ax = plt.subplots(1, len(eps), figsize=(6 * len(eps), 9), squeeze=False)
    for j, e in enumerate(eps):
        tt = T[T.episode == e].sort_values('dist_from_centroid_km')
        for i, r in enumerate(tt.itertuples()):
            d, m = figdat[(e, r.sta)]
            sc = 1.0 / max(np.nanmax(np.abs(figdat[(e, 'CN.PGC')][0])), 1e-15)
            ax[0, j].plot(TREL, d * sc * 0.8 - i, color='k', lw=0.8)
            ax[0, j].plot(TREL, m * sc * 0.8 - i, color='tab:red', lw=0.8, alpha=0.8)
            ax[0, j].text(TREL[0], -i + 0.25, f'{r.sta} {r.dist_from_centroid_km:.0f} km  R={r.R_obs if r.detected else r.R_up:.2f}'
                          + ('' if r.detected else ' (upper)'), fontsize=7)
        ax[0, j].set_title(f'SVI {e} pooled LFE-time Z stacks (black) and thrust synthetic x pooled M0 (red);\n'
                           'common scale = PGC peak; sorted by distance from the pooled centroid', fontsize=8)
        ax[0, j].set_xlabel('s after LFE origin'); ax[0, j].set_yticks([])
    fig.tight_layout(); fig.savefig(f'{C.FIG}/scan2_decay_stacks.png', dpi=130); plt.close(fig)


if __name__ == '__main__':
    {'rows': rows, 'summary': summary}[sys.argv[1]]()
