#!/usr/bin/env python
"""Step 4 (PLAN.md row 10 + Amendment 1 A6): stack-and-invert per template, as H&G did, plus the selection-bias controls.
Stacks: linear mean of the band-passed displacement over detections, per channel of the template's primary config,
window [t-10, t+100) for the inversion and [t-205, t+295) for the Fig.-2-style figure.
Sets per template: 'template' (the event itself), 'all' (1 Jul 2011-1 Dec 2012), 'ETS' (2011 + 2012 ETS), 'inter'
(6 Sep 2011-30 Aug 2012); controls: noise templates j = 0-4 ('noiseJ_single' = the window itself, 'noiseJ_stack' = its
detections above the T-config threshold, self excluded, 1 Jul 2011-1 Dec 2012), 'randomR' (R = 0-4: random live times,
as many as the real 'all' detections).
Outputs tables/stackinv.csv, data/mf/stacks.npz, figures/fig2_repro.png, figures/stackinv_map.png.
usage: hg_stackinv.py gf     (GF caches for tau in {2,5,10,20} on the coarse grid; ~5 min each)
       hg_stackinv.py run"""
import os, sys
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import hg_common as H
import hg_cmt as M
import hg_analyze as A

TAUS = [2, 5, 10, 20]


def cont_all():
    return {s: np.load(f'{H.MF}/cont_{s}.npy', mmap_mode='r') for s in H.STAS}


def stack(X, chans, idx, a=-10, b=100):
    out = np.full((len(chans), b - a), np.nan)
    nused = np.zeros(len(chans), int)
    for q, ch in enumerate(chans):
        s, c = ch.rsplit('.', 1)
        Y = X[s][H.COMPS.index(c)]
        acc = np.zeros(b - a); n = 0
        for i in idx:
            w = np.asarray(Y[i + a:i + b], np.float64)
            if w.size == b - a and np.all(np.isfinite(w)):
                acc += w; n += 1
        if n:
            out[q] = acc / n
        nused[q] = n
    return out, nused


def do_inv(data, chans, G, nodes, tau):
    ok = np.all(np.isfinite(data), axis=1)
    ch = [c for c, k in zip(chans, ok) if k]
    d = data[ok]
    r = M.invert(d, ch, nodes, G)
    fn = M.fine_nodes(r[0])
    Gf = M.gf_cache(fn, tau)
    rf = M.invert(d, ch, fn, Gf)
    b = rf[0]
    out = {'lat': b['lat'], 'lon': b['lon'], 'dep': b['dep'], 'shift': b['shift'], 'VR': b['vr'],
           'coarse_VR': r[0]['vr'], 'n_chan_used': len(ch), 'm': b['m']}
    out.update(M.describe(b['m']))
    return out


def km(a, b, c, d):
    return 6371 * np.hypot(np.radians(c - a), np.radians(d - b) * np.cos(np.radians(a)))


def main_gf():
    nodes = M.coarse_nodes()
    from multiprocessing import Pool
    with Pool(4) as p:
        p.starmap(M.gf_cache, [(nodes, t, 'coarse') for t in TAUS])
    print('gf done', len(nodes))


def choose_tau():
    f = f'{H.TB}/stackinv_tau.json'
    import json
    if os.path.isfile(f):
        return json.load(open(f))['tau']
    nodes = M.coarse_nodes(); X = cont_all()
    R = {K: dict(np.load(f'{H.MF}/real_T{K}.npz', allow_pickle=True)) for K in range(1, 6)}
    tv = {}
    for t in TAUS:
        G = M.gf_cache(nodes, t, 'coarse')
        v = []
        for K in range(1, 6):
            chans = list(R[K]['prim']); i0 = int(R[K]['i0'])
            d, _ = stack(X, chans, [i0])
            v.append(M.invert(d, chans, nodes, G)[0]['vr'])
        tv[t] = float(np.mean(v))
    tau = max(tv, key=tv.get)
    json.dump({'tau': tau, 'mean_VR_templates': tv}, open(f, 'w'))
    print('tau mean VR', tv, '-> tau', tau, flush=True)
    return tau


def run_template(K, thk, tau):
    nodes = M.coarse_nodes(); X = cont_all(); TT = H.templates()
    R, NZ = A.load()
    TH = pd.read_csv(f'{H.TB}/thresholds.csv')
    CAT = pd.read_csv(f'{H.TB}/catalogue_prim_{thk}.csv')
    G = M.gf_cache(nodes, tau, 'coarse')
    rng = np.random.default_rng(H.SEED + 7 + K)
    a18, b18 = H.per_idx('all18')
    rows, ST = [], {}
    r = TT[TT.template == K].iloc[0]
    chans = list(R[K]['prim']); i0 = int(R[K]['i0'])
    c = CAT[CAT.template == K]
    sets = {'template': [i0],
            'all': c[(c.i >= a18) & (c.i < b18)].i.values,
            'ETS': c[((c.i >= H.per_idx('ETS2011')[0]) & (c.i < H.per_idx('ETS2011')[1])) |
                     ((c.i >= H.per_idx('ETS2012')[0]) & (c.i < H.per_idx('ETS2012')[1]))].i.values,
            'inter': c[(c.i >= H.per_idx('interETS')[0]) & (c.i < H.per_idx('interETS')[1])].i.values}
    thr = float(TH[(TH.template == K) & (TH.variant == 'prim')][thk].iloc[0])
    for j in range(5):
        z = NZ[j]
        pi, pv, _ = A.noise_peaks(z, f'T{K}_prim')
        sel = (pv > thr) & (np.abs(pi - int(z['i0'])) > H.L) & (pi >= a18) & (pi < b18)
        sets[f'noise{j}_single'] = [int(z['i0'])]
        sets[f'noise{j}_stack'] = pi[sel]
    nall = len(sets['all'])
    for q in range(5):
        sets[f'random{q}'] = np.sort(rng.integers(a18, b18 - 300, size=max(nall, 1)))
    tsol = None
    for name, idx in sets.items():
        if len(idx) == 0:
            rows.append({'template': K, 'set': name, 'n_det': 0}); continue
        d, nused = stack(X, chans, idx)
        dl, _ = stack(X, chans, idx[:2000] if len(idx) > 2000 else idx, -205, 295)
        ST[f'T{K}_{name}'] = dl
        o = do_inv(d, chans, G, nodes, tau)
        o.update({'template': K, 'set': name, 'n_det': len(idx), 'median_traces_per_chan': int(np.median(nused)),
                  'dist_to_table1_km': km(r.lat, r.lon, o['lat'], o['lon']), 'ddep_table1_km': o['dep'] - r.depth_km,
                  'tau': tau, 'thr_kind': thk})
        if name == 'template':
            tsol = o
        if tsol is not None:
            o['dist_to_template_inv_km'] = km(tsol['lat'], tsol['lon'], o['lat'], o['lon'])
            o['kagan_to_template_deg'] = M.kagan(o['m'], tsol['m'])
            o['cos_m_to_template'] = float(np.dot(o['m'], tsol['m']) / np.linalg.norm(o['m']) / np.linalg.norm(tsol['m']))
        rows.append(o)
        print(K, name, len(idx), round(o['VR'], 3), round(o['lat'], 3), round(o['lon'], 3), o['dep'],
              round(o.get('dist_to_template_inv_km', 0), 1), flush=True)
    return rows, ST


def main_run(thk='thrA_3rd_incl_self'):
    tau = choose_tau()
    from multiprocessing import Pool
    with Pool(4) as p:
        out = p.starmap(run_template, [(K, thk, tau) for K in range(1, 6)])
    rows = [r for o in out for r in o[0]]
    ST = {k: v for o in out for k, v in o[1].items()}
    R = {K: dict(np.load(f'{H.MF}/real_T{K}.npz', allow_pickle=True)) for K in range(1, 6)}
    D = pd.DataFrame(rows)
    D['m'] = D['m'].apply(lambda m: ' '.join(f'{x:.4g}' for x in m) if isinstance(m, np.ndarray) else '')
    D.to_csv(f'{H.TB}/stackinv_{thk}.csv', index=False)
    np.savez_compressed(f'{H.MF}/stacks_{thk}.npz', **ST)
    figs(D, ST, R, thk)


def figs(D, ST, R, thk):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    X = cont_all()
    CAT = pd.read_csv(f'{H.TB}/catalogue_prim_{thk}.csv')
    a18, b18 = H.per_idx('all18')
    panels = [(1, 'UW.SP2.E'), (2, 'TA.B05D.N'), (4, 'TA.A04D.Z'), (1, 'TA.A04D.N'), (2, 'TA.B05D.E'), (3, 'TA.D03D.N')]
    fig, ax = plt.subplots(4, 3, figsize=(13, 9), gridspec_kw={'height_ratios': [1, 2, 1, 2]})
    for p, (K, ch) in enumerate(panels):
        r0, c0 = 2 * (p // 3), p % 3
        c = CAT[(CAT.template == K) & (CAT.i >= a18) & (CAT.i < b18)].sort_values('S', ascending=False)
        s, comp = ch.rsplit('.', 1)
        if s not in X or ch not in list(R[K]['prim']):
            ax[r0, c0].set_title(f'Template {K} - {ch}: not in config', fontsize=9); continue
        Y = X[s][H.COMPS.index(comp)]
        W = [np.asarray(Y[i - 205:i + 295], float) for i in c.i.values]
        W = np.array([w / np.max(np.abs(w)) for w in W if w.size == 500 and np.all(np.isfinite(w))])
        t = np.arange(-205, 295) + 205
        if len(W):
            st = W.mean(0)
            ax[r0, c0].plot(t, st / np.max(np.abs(st)), 'b', lw=1)
            ax[r0 + 1, c0].imshow(W, aspect='auto', cmap='gray', vmin=-1, vmax=1, extent=[0, 500, len(W), 0])
        ax[r0, c0].set_title(f'Template {K} - {s.split(".")[1]} - BH{comp} (n={len(W)})', fontsize=9)
        ax[r0 + 1, c0].set_xlabel('Time (s)'); ax[r0 + 1, c0].set_ylabel('Trace #')
    fig.suptitle('H&G Fig. 2 reproduction (2011 templates; each detection normalised, white = positive; stack on top)')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/fig2_repro_{thk}.png', dpi=120); plt.close(fig)
    # map of solutions
    TT = H.templates()
    fig, ax = plt.subplots(figsize=(7, 7))
    C = M.coords()
    for s, (la, lo) in C.items():
        ax.plot(lo, la, 'D', color='orange'); ax.text(lo, la, s.split('.')[1], fontsize=7)
    cols = dict(zip(range(1, 6), ['C0', 'C1', 'C2', 'C3', 'C4']))
    for _, r in TT.iterrows():
        ax.plot(r.lon, r.lat, '*', ms=14, color=cols[r.template], mec='k')
    for _, d in D.dropna(subset=['VR']).iterrows():
        mk = {'template': 'o', 'all': 's', 'ETS': '^', 'inter': 'v'}.get(d.set, 'x' if 'noise' in d.set else '.')
        ax.scatter(d.lon, d.lat, marker=mk, s=40 if mk != '.' else 15, color=cols[d.template],
                   alpha=1 if mk in 'os^v' else 0.5)
    ax.set_title('stars: Table 1; o template inv; s all; ^ ETS; v inter; x noise-template stacks; . random', fontsize=8)
    ax.set_xlabel('lon'); ax.set_ylabel('lat')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/stackinv_map_{thk}.png', dpi=120); plt.close(fig)


if __name__ == '__main__':
    if sys.argv[1] == 'gf':
        main_gf()
    else:
        import json
        best = json.load(open(f'{H.TB}/best_threshold.json'))['best_thr_kind_K1']
        for thk in (sys.argv[2:] or list(dict.fromkeys(['thrA_3rd_incl_self', best]))):
            print('=====', thk, flush=True)
            main_run(thk)
