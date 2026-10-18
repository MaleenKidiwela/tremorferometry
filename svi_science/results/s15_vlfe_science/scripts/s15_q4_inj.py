#!/usr/bin/env python
"""S15 Q4 (PLAN Amendment 5): Green's-function thrust sources at 7 scan-grid nodes + C0/C1-node GF references, injected
through a production run's own detector, bank and battery (cat_pilot.inj_row unchanged). Run with the run's environment
(see s15_q4_queue.sh).  usage: s15_q4_inj.py sources|inj YEAR"""
import os, sys
os.environ['OMP_NUM_THREADS'] = '1'
CATS = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue/scripts'
S13S = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/scripts'
sys.path[:0] = [CATS, S13S, '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro/scripts']
import numpy as np, pandas as pd

OUT = '/home/jovyan/tremorferometry/svi_science/results/s15_vlfe_science/data/q4'
os.makedirs(OUT, exist_ok=True)
NODES = {'G0': 77, 'G1': 94, 'N1': 112, 'N2': 109, 'N3': 126, 'N4': 143, 'N5': 176, 'W1': 14, 'W2': 80}
MWREF = 3.6
NPER = 150
NW = int(os.environ.get('CAT_NW', '4'))


def build_sources(year):
    import cat_pilot as P
    import cat_inj as I
    import scan_common as S
    import s13_green as G
    import scan_refs as SR
    import hg_cont
    from scipy.signal import sosfiltfilt
    Gr = np.load(os.environ['SCAN_GRID']); gst = list(Gr['stations'])
    cc = SR.all_coords()
    m0 = 10 ** (1.5 * MWREF + 9.1)
    M = G.mt_ned(*S.MECHS['thrust'], m0)
    t = np.arange(-400, 601).astype(float); o = 400
    z, meta, chk = {}, [], []
    mf_sta = sorted({ch.rsplit('.', 1)[0] for ch in P.CH})
    for src, node in NODES.items():
        lat, lon, dep = float(Gr['lat'][node]), float(Gr['lon'][node]), float(Gr['dep'][node])
        nmf = 0
        for s in mf_sta:
            if s not in cc.index:
                continue
            raw = S.raw_disp(M, (lat, lon, dep), (float(cc.loc[s, 'lat']), float(cc.loc[s, 'lon'])), t)      # (3, n) N, E, Z
            for ci, c in enumerate('NEZ'):
                ch = f'{s}.{c}'
                if ch not in P.CH:
                    continue
                x = raw[ci] - raw[ci][:50].mean()
                xp = np.r_[np.zeros(600), x, np.full(600, x[-1])]
                f = sosfiltfilt(hg_cont.SOS, xp)[600:-600]
                v = np.zeros_like(f); v[1:-1] = 0.5 * (f[2:] - f[:-2])
                z[f'{src}|mf|{ch}'] = v[o + I.W0:o + I.W1].copy(); nmf += 1
        for j, s in enumerate(gst):
            z[f'{src}|sc|{s}'] = m0 * Gr['thrust'][node, j]
        meta.append(dict(source=src, members='', n=0, Mw_ref=MWREF, Mw_members='', node=node, lat=lat, lon=lon, dep=dep,
                         n_mf_ch=nmf, scan_stations=','.join(gst)))
    # check: GF reference vs the empirical C0/C1 sources of this run's source file
    Se, Me = I.load_sources()
    for g, e in (('G0', 'C0'), ('G1', 'C1')):
        for kind in ('mf', 'sc'):
            for ch, w in Se[e][kind].items():
                k = f'{g}|{kind}|{ch}'
                if k not in z:
                    continue
                a = np.asarray(w, float).ravel(); b = np.asarray(z[k], float).ravel()
                sc_e = 10 ** (1.5 * (MWREF - Me.loc[e, 'Mw_ref']))
                cc_ = float(np.dot(a - a.mean(), b - b.mean()) / (np.linalg.norm(a - a.mean()) * np.linalg.norm(b - b.mean()) + 1e-30))
                chk.append(dict(year=year, gf=g, empirical=e, kind=kind, ch=ch, cc=cc_,
                                amp_ratio_gf_over_emp=float(np.linalg.norm(b) / (np.linalg.norm(a) * sc_e + 1e-30))))
    np.savez(f'{OUT}/sources_{year}.npz', **z)
    pd.DataFrame(meta).to_csv(f'{OUT}/sources_{year}.csv', index=False)
    C = pd.DataFrame(chk); C.to_csv(f'{OUT}/check_{year}.csv', index=False)
    print(pd.DataFrame(meta)[['source', 'node', 'lat', 'lon', 'dep', 'n_mf_ch']].to_string())
    print(C.groupby(['gf', 'kind']).agg(n=('cc', 'size'), cc_med=('cc', 'median'), amp_med=('amp_ratio_gf_over_emp', 'median')).round(3))


_SRC = {}


def my_init():
    import cat_pilot as P
    import cat_inj as I
    I.bat_init(early=False)
    year = int(os.environ['SCAN_YEAR'])
    zz = np.load(f'{OUT}/sources_{year}.npz')
    S = {}
    for k in zz.files:
        src, kind, ch = k.split('|')
        S.setdefault(src, {'mf': {}, 'sc': {}})[kind][ch] = zz[k]
    SM = pd.read_csv(f'{OUT}/sources_{year}.csv').set_index('source')
    full = P.bank_local([])
    P._P.update(I=I, S=S, SM=SM, get=P.getter(), bank={src: full for src in SM.index})


def run_inj(year):
    import cat_pilot as P
    from multiprocessing import Pool
    L = pd.read_parquet(f'{P.PD}/inj_list.parquet')
    base = L[(L.source == 'C0') & (L.Mw == 3.6)].sort_values('inj_id')
    picks = pd.concat([base[base.period == p].head(NPER) for p in ('ETS', 'between')])
    rows, k = [], 0
    for src in NODES:
        for mw in (3.6, 3.8):
            for r in picks.itertuples():
                rows.append(dict(inj_id=k, source=src, Mw=mw, i=int(r.i), epoch=float(r.epoch), period=r.period, base_inj_id=int(r.inj_id)))
                k += 1
    Lst = pd.DataFrame(rows); Lst.to_parquet(f'{OUT}/inj_list_{year}.parquet')
    recs = Lst.to_dict('records')
    nch = NW * 6; out = []
    with Pool(NW, initializer=my_init) as p:
        for r in p.imap_unordered(P.inj_chunk, [recs[j::nch] for j in range(nch)]):
            out += r
    pd.DataFrame(out).to_parquet(f'{OUT}/inj_{year}.parquet')
    print('inj done', year, len(out))


if __name__ == '__main__':
    mode, year = sys.argv[1], int(sys.argv[2])
    assert os.environ.get('SCAN_YEAR') == str(year)
    if mode == 'sources':
        build_sources(year)
    elif mode == 'inj':
        run_inj(year)
