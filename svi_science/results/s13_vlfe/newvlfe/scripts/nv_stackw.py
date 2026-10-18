#!/usr/bin/env python
"""Point 5 for the 2018 western group (coordinator 2026-10-16 ~11:05 am PT).
Part A (compare): stack the 6 LFE-coincident western candidates (W6) and all 17 western events (W17) with the catalogue
relay recipe (cat_relay.stack_frame logic: 90-s velocity windows 0.02-0.05 Hz from origin - 9 s, unit-normalised per
channel, +-5 s alignment to the stack of the others, 2 passes; channel weights = mean template-RMS / pre-origin RMS),
on the 2018 production station set (catalogue/data/prod2018/cont, read only). Then:
  - channel-sum CC (mean over common channels) at the best lag (-20..20 s) against every bank template (T1-T5, C0-C5,
    RC0-RC7) + reference pairs (relay LOO vs primary = same source; RC pairs = different clusters);
  - station amplitude pattern (scan-band 20-50 s Z displacement, max |Z| in [t0, t0 + 180 s], per member, median),
    W events vs RC4 members (R12, R13) vs RC0 members;
  - stacked scan location: per node sum_e num_e / sqrt(sum_e den_e) on the 2018 v2 scan (variant 0), for W6, W17, RC4, RC0.
-> data/stackw/bank.npz, tables/nv_stackw_cc.csv, nv_stackw_amp.csv, nv_stackw_loc.csv, figures/nv_stackw_section.png,
   nv_stackw_amp.png, nv_stackw_loc.png. Nothing is written under catalogue/."""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
os.environ.setdefault('PROD_YEAR', '2018'); os.environ.setdefault('SCAN_YEAR', '2018')
S13 = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe'
NV = f'{S13}/newvlfe'; CAT = f'{S13}/catalogue'
for p in (f'{S13}/scripts', f'{S13}/hg_repro/scripts', f'{CAT}/scripts', os.path.dirname(os.path.abspath(__file__))):
    sys.path.insert(0, p)
import numpy as np, pandas as pd
import cat_pilot as P
import cat_common as K
OUT = f'{NV}/data/stackw'; os.makedirs(OUT, exist_ok=True)


def stack_frame(starts):
    W = [P.cut_frame(int(s)) for s in starts]
    out = {}
    for ch in P.CH:
        ws = [w[ch] for w in W if w[ch] is not None]
        if not ws or len(ws) < 0.5 * len(W):
            out[ch] = None; continue
        s = np.mean(ws, 0); s = s - s.mean(); n = np.linalg.norm(s); out[ch] = s / n if n > 0 else None
    return out


def sumcc(a, b):
    v = [float(a[c] @ b[c]) for c in P.CH if a.get(c) is not None and b.get(c) is not None]
    return (sum(v), len(v))


def aligned(starts):
    starts = list(map(int, starts)); sh = np.zeros(len(starts), int)
    if len(starts) >= 2:
        for _ in range(2):
            for k, s0 in enumerate(starts):
                ref = stack_frame([starts[j] + sh[j] for j in range(len(starts)) if j != k])
                sh[k] = max(range(-5, 6), key=lambda d: sumcc(P.cut_frame(s0 + d), ref)[0])
    st = [s + d for s, d in zip(starts, sh)]
    return stack_frame(st), st, sh


def weights(starts):
    rr = [P.rms_ratio_frame(int(s)) for s in starts]
    return {ch: float(np.nanmean([r[ch] for r in rr])) if np.isfinite([r[ch] for r in rr]).any() else np.nan for ch in P.CH}


def as_dict(arr):
    return {ch: (arr[k] if np.isfinite(arr[k]).all() else None) for k, ch in enumerate(P.CH)}


def lag_cc(starts, tpl, lags=range(-20, 21)):
    best = (-9, 0, 0)
    for d in lags:
        s = stack_frame([x + d for x in starts]); v, n = sumcc(s, tpl)
        if n and v / n > best[0]:
            best = (v / n, d, n)
    return best


def main():
    Wd = pd.read_csv(f'{NV}/tables/nv_w2b_west.csv')
    Wd = Wd[(Wd.run == 'after W2b') & (Wd.kind == 'real')].sort_values('epoch')
    sets = {'W17': Wd, 'W6': Wd[Wd.patch_p <= 0.05]}
    B = np.load(f'{CAT}/data/prodr2018/bank.npz')           # T1-T5, C0-C5, RC0-RC7 (+ noise) on the 2018 channels
    assert list(B['chans']) == P.CH
    RB = np.load(f'{CAT}/data/prod2018/relay/bank.npz')     # relay primaries + LOO stacks
    bank = {n: as_dict(B[n]) for n in [f'T{k}' for k in range(1, 6)] + [f'C{k}' for k in range(6)] + [f'RC{k}' for k in range(8)] if n in B.files}
    REL = pd.read_csv(f'{CAT}/tables/prod2018_relay_events.csv')
    stacks, members = {}, {}
    for nm, g in sets.items():
        starts = (g.epoch.values - P.EP0 - P.PRE).round().astype(int)
        st, al, sh = aligned(starts)
        stacks[nm] = st; members[nm] = al
        print(nm, len(al), 'shifts', sh.tolist(), 'live ch', sum(v is not None for v in st.values()))
    # save the candidate templates (newvlfe only)
    arr = {nm: np.array([st[ch] if st[ch] is not None else np.full(P.L, np.nan) for ch in P.CH]) for nm, st in stacks.items()}
    wts = {f'w__{nm}': np.array([weights(members[nm])[ch] for ch in P.CH]) for nm in stacks}
    np.savez(f'{OUT}/bank.npz', chans=np.array(P.CH), **arr, **wts,
             **{f'starts__{nm}': np.array(members[nm]) for nm in members})
    # CC against the bank at the best lag (stack rebuilt from the aligned members shifted by the lag)
    rows = []
    for nm in stacks:
        for t, tpl in bank.items():
            c, d, n = lag_cc(members[nm], tpl)
            rows.append(dict(stack=nm, template=t, cc_mean=c, lag_s=d, n_ch=n))
    # reference pairs: LOO stack vs its primary (same source) and RC primaries vs each other (different clusters)
    rel = [n for n in RB.files if n.startswith('RC') and '_loo' not in n]
    for a in rel:
        for b in rel:
            if a < b:
                A, Bb = as_dict(RB[a]), as_dict(RB[b])
                best = (-9, 0, 0)
                for d in range(-20, 21):
                    v, n = sumcc({c: (np.roll(A[c], d) if A[c] is not None else None) for c in P.CH}, Bb)
                    if n and v / n > best[0]:
                        best = (v / n, d, n)
                rows.append(dict(stack=a, template=b, cc_mean=best[0], lag_s=best[1], n_ch=best[2], ref='RC pair (roll)'))
    for n in RB.files:
        if '_loo' in n and not n.startswith('w__'):
            prim = n.split('_loo')[0]; v, k = sumcc(as_dict(RB[n]), as_dict(RB[prim]))
            rows.append(dict(stack=n, template=prim, cc_mean=v / max(k, 1), lag_s=0, n_ch=k, ref='LOO vs primary'))
    # single members of RC4 / RC0 vs their own primary's LOO (how well one event matches its own source stack)
    for c in (0, 4):
        g = REL[REL.cluster == c]
        for r in g.itertuples():
            loo = f'RC{c}_loo{r.rid}'
            if loo in RB.files:
                cc_, d, n = lag_cc([int(r.start0 + r.shift_s)], as_dict(RB[loo]))
                rows.append(dict(stack=f'{r.rid} (single)', template=loo, cc_mean=cc_, lag_s=d, n_ch=n, ref='member vs own LOO'))
    CC = pd.DataFrame(rows); CC.to_csv(f'{NV}/tables/nv_stackw_cc.csv', index=False)
    pd.set_option('display.width', 200)
    print(CC[CC['stack'].isin(['W6', 'W17'])].sort_values(['stack', 'cc_mean'], ascending=[True, False]).groupby('stack').head(8).round(3).to_string())
    print(CC[CC.get('ref').notna()].groupby('ref').cc_mean.describe().round(3))
    # amplitude pattern (scan band)
    SD = f'{S13}/data/scan_2018'
    stas = [s for s in json.load(open(f'{SD}/runs/real_v2/meta.json'))['stations'] if s not in ('UW.LRIV', 'CN.LZB')]
    y0 = (pd.Timestamp('2018-01-01') - pd.Timestamp('1970-01-01')).total_seconds()
    X = {s: np.load(f'{SD}/bp_2018v2/{s}.npy', mmap_mode='r') for s in stas}
    groups = {'W6': sets['W6'].epoch.values, 'W17': sets['W17'].epoch.values,
              'RC4': (REL[REL.cluster == 4].epoch.values), 'RC0': REL[REL.cluster == 0].epoch.values,
              'RC2': REL[REL.cluster == 2].epoch.values}
    arows = []
    for gname, eps in groups.items():
        for ep in eps:
            t = int(round(ep - y0))
            for s in stas:
                x = np.asarray(X[s][2, t:t + 180]) * 1e9; pre = np.asarray(X[s][2, t - 420:t - 240]) * 1e9
                if np.isfinite(x).sum() > 150 and np.isfinite(pre).sum() > 150:
                    arows.append(dict(group=gname, epoch=ep, sta=s, amp=float(np.nanmax(np.abs(x))), pre=float(np.nanmax(np.abs(pre)))))
    A = pd.DataFrame(arows); A.to_csv(f'{NV}/tables/nv_stackw_amp.csv', index=False)
    Am = A.groupby(['group', 'sta'])[['amp', 'pre']].median().reset_index()
    Am['snr'] = Am.amp / Am.pre
    piv = Am.pivot(index='sta', columns='group', values='amp')
    print(Am.pivot(index='sta', columns='group', values='snr').round(2).to_string())
    lp = np.log(piv); lp = lp - lp.mean()
    print('log-amplitude pattern correlation\n', lp.corr().round(2).to_string())
    # stacked scan location
    import nv_battery as NB
    NB.winit('real'); SC = NB.SC; W = SC.W
    nn = len(W['nodes']) - 1
    def numden(tsec):
        d, u = divmod(int(tsec), 86400); g = W['g'][d]
        num = np.zeros(nn); den = np.zeros(nn)
        for i, s in enumerate(W['stas']):
            x = SC.cleaned_window(i, s, int(d) * 86400 + u, tsec, NB.S.FIT0, NB.S.FIT1 + 1)
            sig = W['DS'][s]['sig'][d] * 1e9
            for c in range(3):
                if W['won'][d, i, c] and np.isfinite(x[c]).all() and np.isfinite(sig[c]) and sig[c] > 0:
                    T = W['T']['thrust'][:nn, i, c] * g[i, c]; w = 1 / sig[c] ** 2
                    num += w * (T @ x[c]); den += w * (T * T).sum(1)
        return num, den
    N = W['nodes'].iloc[:nn]
    lrows, Z = [], {}
    for gname, eps in groups.items():
        nu = np.zeros(nn); de = np.zeros(nn)
        for ep in eps:
            t = int(round(ep - y0))
            # best scan time within +-10 s of the listed origin (relay events: origin-aligned MF index)
            cand = []
            for k in range(-10, 11, 2):
                a, b = numden(t + k)
                cand.append((np.nanmax(np.where(b > 0, a / np.sqrt(np.maximum(b, 1e-30)), -np.inf)), a, b))
            _, a, b = max(cand, key=lambda q: q[0]); nu += a; de += b
        z = np.where(de > 0, nu / np.sqrt(de), np.nan); Z[gname] = z
        k = int(np.nanargmax(z))
        lrows.append(dict(group=gname, n=len(eps), lat=N.lat.values[k], lon=N.lon.values[k], dep=N.dep.values[k], zrel_max=float(z[k])))
    Lc = pd.DataFrame(lrows); Lc.to_csv(f'{NV}/tables/nv_stackw_loc.csv', index=False); print(Lc.round(3).to_string())
    np.savez(f'{OUT}/loc_maps.npz', lat=N.lat.values, lon=N.lon.values, **Z)
    figs(stacks, bank, piv, Z, N, Lc)


def figs(stacks, bank, piv, Z, N, Lc):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 9), sharey=True)
    for a, nm in zip(ax, ('W6', 'W17', 'RC4')):
        st = stacks.get(nm, bank.get(nm))
        for k, ch in enumerate(P.CH):
            if st[ch] is not None:
                a.plot(np.arange(P.L) - P.PRE, k + 2.5 * st[ch], lw=0.8, color='tab:red' if nm != 'RC4' else 'k')
        a.set_title(nm + (' (relay template)' if nm == 'RC4' else ' stack')); a.set_xlabel('s after origin')
    ax[0].set_yticks(range(len(P.CH))); ax[0].set_yticklabels(P.CH, fontsize=7)
    fig.suptitle('Raw stacks: unit-normalised 0.02–0.05 Hz velocity per channel (catalogue V4 recipe), 2018 station set')
    fig.tight_layout(); fig.savefig(f'{NV}/figures/nv_stackw_section.png', dpi=110); plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 5))
    for g, mk in (('W6', 'o-'), ('W17', 's--'), ('RC4', '^-'), ('RC0', 'v:'), ('RC2', 'd:')):
        if g in piv:
            ax.semilogy(range(len(piv)), piv[g], mk, label=g)
    ax.set_xticks(range(len(piv))); ax.set_xticklabels([s.split('.')[1] for s in piv.index], rotation=45)
    ax.set_ylabel('median max |Z| (nm), 20–50 s, [t0, t0+180 s]'); ax.legend(); ax.set_title('Station amplitude pattern')
    fig.tight_layout(); fig.savefig(f'{NV}/figures/nv_stackw_amp.png', dpi=110); plt.close(fig)
    fig, ax = plt.subplots(1, len(Z), figsize=(4 * len(Z), 4.2), sharey=True)
    for a, (g, z) in zip(ax, Z.items()):
        sc = a.scatter(N.lon, N.lat, c=z / np.nanmax(z), s=30, cmap='viridis', vmin=0, vmax=1)
        r = Lc[Lc.group == g].iloc[0]; a.plot(r.lon, r.lat, 'r*', ms=14); a.set_title(f'{g}: max {r.lat:.2f} N {r.lon:.2f}')
    fig.colorbar(sc, ax=ax, label='stacked scan statistic / max'); fig.savefig(f'{NV}/figures/nv_stackw_loc.png', dpi=110); plt.close(fig)


if __name__ == '__main__':
    main()
