#!/usr/bin/env python
"""Amendment R2 (PLAN.md): reverse-engineered reproduction of H&G (2019). Summed-CC matched filter with the pre-listed
processing variants, plus the controls and noise templates at the chosen setting.

A setting ("spec") = kind,off,norm,chset,hw
  kind  disp | vel      displacement (R1 cont arrays) or its time derivative (data/r2/contv_*.npy)
  off   0 | -30         template window [origin+off, origin+off+90]; detection index = window start - off (origin time)
  norm  sum | mean | net  plain sum of per-channel CC (V0) / sum / n_live (V1) / one normalised CC over all live
                        channels concatenated (V2: sum_c |w_c| cc_c den_c / sqrt(sum_live |w_c|^2 * sum_live den_c^2))
  chset all | snr15 | near5   channels live at the template / template SNR >= 1.5 / 5 stations nearest the epicentre
  hw    90 | 30         declustering half-width (s): a peak is the max of S over +-hw
Modes
  vel STA[,STA]            build data/r2/contv_STA.npy (float32 (3,N), NaN dead) and denv_STA.npy
  var K                    stage 1: template K under V0..V6, V5b (window origin-9 s), V6b (decluster 20 s) (R2a) (one CC pass per kind/off) -> data/r2/var_T{K}.npz
  real K SPEC              stage 2: template K under SPEC + 8 station-day-shift + 8 polarity-scrambled draws (R1 draws)
                           -> data/r2/real_{tag}_T{K}.npz
  noise J0 J1 SPEC         stage 2: R1 noise templates j in [J0, J1) under each template's channel set
                           -> data/r2/noise_{tag}_JJJ.npz
Peaks: all positive declustered peaks above median(S) are stored (index = origin-aligned sample, value, n_live)."""
import os, sys
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.signal import oaconvolve
from scipy.ndimage import maximum_filter1d
import hg_common as H
import hg_cont
import hg_mf

R2 = f'{H.HG}/data/r2'
os.makedirs(R2, exist_ok=True)
COORD = {'CN.NLLB': (49.225, -124.087), 'CN.VGZ': (48.413, -123.325), 'UW.LRIV': (48.0575, -123.504),
         'UW.DOSE': (47.7172, -122.9718), 'UW.GNW': (47.5641, -122.8250), 'UW.SP2': (47.5563, -122.2492),
         'TA.A04D': (48.7201, -122.7063), 'TA.B05D': (48.2641, -122.0960), 'TA.D03D': (47.5347, -123.0894)}
_C = {}


def cont(sta, kind):
    k = (sta, kind)
    if k not in _C:
        if kind == 'disp':
            _C[k] = (np.load(f'{H.MF}/cont_{sta}.npy', mmap_mode='r'), np.load(f'{H.MF}/den_{sta}.npy', mmap_mode='r'))
        else:
            _C[k] = (np.load(f'{R2}/contv_{sta}.npy', mmap_mode='r'), np.load(f'{R2}/denv_{sta}.npy', mmap_mode='r'))
    return _C[k]


def build_vel(sta):
    out = f'{R2}/contv_{sta}.npy'
    if os.path.isfile(out) and os.path.isfile(f'{R2}/denv_{sta}.npy'):
        return f'{sta} exists'
    X = np.load(f'{H.MF}/cont_{sta}.npy', mmap_mode='r')
    V = np.full(X.shape, np.nan, np.float32)
    for c in range(3):
        x = np.asarray(X[c], np.float64)
        v = np.full(x.shape, np.nan)
        v[1:-1] = 0.5 * (x[2:] - x[:-2])          # central difference (m/s at 1 sps); NaN propagates at gaps
        V[c] = v.astype(np.float32)
    np.save(out + '.part.npy', V); os.replace(out + '.part.npy', out)
    D = hg_cont.dens(V)
    np.save(f'{R2}/denv_{sta}.part.npy', D); os.replace(f'{R2}/denv_{sta}.part.npy', f'{R2}/denv_{sta}.npy')
    return f'{sta} vel live {(D > 0).mean(axis=1).round(3)}'


def tag(spec):
    return spec.replace(',', '_').replace('-', 'm')


def parse(spec):
    kind, off, norm, chset, hw = spec.split(',')
    return {'kind': kind, 'off': int(off), 'norm': norm, 'chset': chset, 'hw': int(hw), 'spec': spec}


def cut(i0, kind):
    """27-channel window at sample i0: ch -> (unit w, raw norm, snr) or None."""
    out = {}
    for s in H.STAS:
        X, D = cont(s, kind)
        for c in range(3):
            ch = f'{s}.{H.COMPS[c]}'
            w = np.array(X[c, i0:i0 + H.L], np.float64)
            if not np.all(np.isfinite(w)) or D[c, i0] <= 0:
                out[ch] = None
                continue
            pre = np.array(X[c, i0 - 1800:i0], np.float64); pre = pre[np.isfinite(pre)]
            snr = np.sqrt(np.mean(w ** 2)) / np.sqrt(np.mean(pre ** 2)) if pre.size > 600 else np.nan
            w = w - w.mean(); n = np.linalg.norm(w)
            out[ch] = (w / n, n, snr)
    return out


def chans_for(K, T, chset):
    """channel list of template K's configuration (T = its cut at the template origin + off)."""
    live = [ch for ch in H.CHANS if T[ch] is not None]
    if chset == 'all':
        return live
    if chset == 'snr15':
        return [ch for ch in live if np.isfinite(T[ch][2]) and T[ch][2] >= 1.5]
    if chset == 'near5':
        r = H.templates().set_index('template').loc[K]
        d = {s: np.hypot(COORD[s][0] - r.lat, (COORD[s][1] - r.lon) * np.cos(np.radians(r.lat))) for s in H.STAS}
        near = sorted(d, key=d.get)[:5]
        return [ch for ch in live if ch.rsplit('.', 1)[0] in near]
    raise ValueError(chset)


def ccser(ch, w, kind):
    s, comp = ch.rsplit('.', 1)
    c = H.COMPS.index(comp)
    X, D = cont(s, kind)
    x = np.nan_to_num(np.asarray(X[c], np.float32))
    num = oaconvolve(x, w[::-1].astype(np.float32), mode='valid')
    del x
    d = np.asarray(D[c])
    live = d > 0
    cc = np.zeros(H.M, np.float32)
    np.divide(num, d, out=cc, where=live)
    np.clip(cc, -1, 1, out=cc)
    return cc, live, d


class Acc:
    """accumulates one summed statistic (any norm) with an optional per-station shift and per-channel sign."""
    def __init__(self, norm):
        self.norm = norm
        self.S = np.zeros(H.M, np.float32)
        self.NL = np.zeros(H.M, np.int16)
        if norm == 'net':
            self.D2 = np.zeros(H.M, np.float64); self.W2 = np.zeros(H.M, np.float32)

    def add(self, cc, live, d, wn, sign=1.0, shift=0):
        r = (lambda a: np.roll(a, shift)) if shift else (lambda a: a)
        if self.norm == 'net':
            self.S += r(cc * (d * wn)) * sign
            self.D2 += r(np.where(live, d.astype(np.float64) ** 2, 0.0))
            self.W2 += r(live.astype(np.float32) * np.float32(wn * wn))
        else:
            self.S += r(cc) * sign if sign != 1.0 else r(cc)
        self.NL += r(live).astype(np.int16)

    def final(self):
        if self.norm == 'sum':
            return self.S
        if self.norm == 'mean':
            out = np.zeros(H.M, np.float32)
            np.divide(self.S, self.NL, out=out, where=self.NL > 0)
            return out
        den = np.sqrt(self.D2 * self.W2)
        out = np.zeros(H.M, np.float32)
        np.divide(self.S, den, out=out, where=den > 0)
        return out


def peaks(S, NL, hw, off):
    """positive declustered peaks above the median of S over live samples; index shifted to origin time."""
    v = S[NL > 0]
    med = float(np.median(v)); mad = float(np.median(np.abs(v - med)))
    mf = maximum_filter1d(S, 2 * hw + 1, mode='constant', cval=-np.inf)
    i = np.where((S == mf) & (S > med) & (NL > 0))[0]
    if i.size > 1:
        i = i[np.concatenate([[True], np.diff(i) > hw])]
    return (i - off).astype(np.int64), S[i].astype(np.float32), NL[i].astype(np.int16), med, mad


def store(out, tag, S, NL, hw, off):
    pi, pv, pn, med, mad = peaks(S, NL, hw, off)
    out[f'{tag}_pi'] = pi; out[f'{tag}_pv'] = pv; out[f'{tag}_pnl'] = pn
    out[f'{tag}_med'] = med; out[f'{tag}_mad'] = mad


def run_var(K):
    """stage 1: V0-V6 for template K."""
    f = f'{R2}/var_T{K}.npz'
    if os.path.isfile(f):
        return f'T{K} exists'
    row = H.templates().set_index('template').loc[K]
    io = H.tidx(row.time)
    out = {'K': K, 'i0': io}
    # pass 1: displacement, off 0 -> V0 (sum, hw90), V1 (mean), V2 (net), V3a (snr15), V3b (near5), V6 (sum, hw30)
    T = cut(io, 'disp')
    ch_all, ch_snr, ch_near = chans_for(K, T, 'all'), chans_for(K, T, 'snr15'), chans_for(K, T, 'near5')
    acc = {'V0': Acc('sum'), 'V2': Acc('net'), 'V3a': Acc('sum'), 'V3b': Acc('sum')}
    for ch in ch_all:
        cc, live, d = ccser(ch, T[ch][0], 'disp')
        acc['V0'].add(cc, live, d, T[ch][1]); acc['V2'].add(cc, live, d, T[ch][1])
        if ch in ch_snr:
            acc['V3a'].add(cc, live, d, T[ch][1])
        if ch in ch_near:
            acc['V3b'].add(cc, live, d, T[ch][1])
        del cc, live, d
    acc['V1'] = acc['V0']
    for tag, hw, a in [('V0', 90, 'V0'), ('V6', 30, 'V0'), ('V6b', 20, 'V0'), ('V2', 90, 'V2'), ('V3a', 90, 'V3a'),
                       ('V3b', 90, 'V3b')]:
        store(out, tag, acc[a].final(), acc[a].NL, hw, 0)
    a = Acc('mean'); a.S, a.NL = acc['V0'].S, acc['V0'].NL
    store(out, 'V1', a.final(), a.NL, 90, 0)
    out['chans_V0'] = np.array(ch_all); out['chans_V3a'] = np.array(ch_snr); out['chans_V3b'] = np.array(ch_near)
    del acc, a
    # pass 2: velocity (V4); pass 3: window [origin-30, origin+60] (V5)
    for tag, kind, off in [('V4', 'vel', 0), ('V5', 'disp', -30), ('V5b', 'disp', -9)]:
        T = cut(io + off, kind)
        chs = chans_for(K, T, 'all')
        a = Acc('sum')
        for ch in chs:
            cc, live, d = ccser(ch, T[ch][0], kind)
            a.add(cc, live, d, T[ch][1]); del cc, live, d
        store(out, tag, a.final(), a.NL, 90, off)
        out[f'chans_{tag}'] = np.array(chs)
        del a
    np.savez_compressed(f + '.part.npz', **out); os.replace(f + '.part.npz', f)
    return f'T{K} var done: n_ch V0 {len(ch_all)} V3a {len(ch_snr)} V3b {len(ch_near)}'


def run_real(K, spec):
    """stage 2: template K under spec, plus A1 (station-day-shift) and A2 (polarity-scrambled) draws (R1 draws)."""
    P = parse(spec)
    f = f'{R2}/real_{tag(spec)}_T{K}.npz'
    if os.path.isfile(f):
        return f'T{K} exists'
    row = H.templates().set_index('template').loc[K]
    io = H.tidx(row.time)
    T = cut(io + P['off'], P['kind'])
    chs = chans_for(K, T, P['chset'])
    offs, signs = hg_mf.draws()
    names = ['real'] + [f'A1_{r}' for r in range(hg_mf.NDRAW)] + [f'A2_{r}' for r in range(hg_mf.NDRAW)]
    acc = {n: Acc(P['norm']) for n in names}
    for ch in chs:
        ci = H.CHANS.index(ch); si = ci // 3
        cc, live, d = ccser(ch, T[ch][0], P['kind'])
        acc['real'].add(cc, live, d, T[ch][1])
        for r in range(hg_mf.NDRAW):
            acc[f'A1_{r}'].add(cc, live, d, T[ch][1], shift=int(offs[r, si]) * 86400)
            acc[f'A2_{r}'].add(cc, live, d, T[ch][1], sign=float(signs[r, ci]))
        del cc, live, d
    out = {'K': K, 'i0': io, 'spec': spec, 'chans': np.array(chs)}
    for n in names:
        store(out, n, acc[n].final(), acc[n].NL, P['hw'], P['off'])
        acc[n] = None
    np.savez_compressed(f + '.part.npz', **out); os.replace(f + '.part.npz', f)
    return f'T{K} real+controls done ({len(chs)} ch)'


def run_noise(J0, J1, spec):
    P = parse(spec)
    NT = pd.read_csv(f'{H.MF}/noise_times.csv')
    TT = H.templates()
    cfg = {}
    for K in TT.template:
        io = H.tidx(TT.set_index('template').loc[K].time)
        cfg[K] = chans_for(K, cut(io + P['off'], P['kind']), P['chset'])
    uniq = {}
    for K, v in cfg.items():
        uniq.setdefault(tuple(v), []).append(K)
    for j in range(J0, J1):
        f = f'{R2}/noise_{tag(spec)}_{j:03d}.npz'
        if os.path.isfile(f):
            continue
        i0 = int(NT.i[j])
        T = cut(i0, P['kind'])
        cl = list(uniq.items())
        acc = [Acc(P['norm']) for _ in cl]
        need = sorted({ch for v, _ in cl for ch in v})
        for ch in need:
            if T[ch] is None:
                continue
            cc, live, d = ccser(ch, T[ch][0], P['kind'])
            for q, (v, _) in enumerate(cl):
                if ch in v:
                    acc[q].add(cc, live, d, T[ch][1])
            del cc, live, d
        out = {'j': j, 'i0': i0, 'ncfg': len(cl)}
        for q, (v, Ks) in enumerate(cl):
            store(out, f'c{q}', acc[q].final(), acc[q].NL, P['hw'], 0)
            out[f'c{q}_templates'] = np.array(Ks); out[f'c{q}_n'] = len(v)
            acc[q] = None
        np.savez_compressed(f + '.part.npz', **out); os.replace(f + '.part.npz', f)
        print(f'{pd.Timestamp.now(tz="America/Los_Angeles"):%H:%M:%S} PT noise {j} done', flush=True)
    return f'noise {J0}-{J1} done'


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'vel':
        from multiprocessing import Pool
        with Pool(min(5, len(sys.argv[2].split(',')))) as p:
            r = '; '.join(p.map(build_vel, sys.argv[2].split(',')))
    elif mode == 'var':
        r = run_var(int(sys.argv[2]))
    elif mode == 'real':
        r = run_real(int(sys.argv[2]), sys.argv[3])
    elif mode == 'noise':
        r = run_noise(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    line = f'{pd.Timestamp.now(tz="America/Los_Angeles"):%Y-%m-%d %H:%M:%S} PT r2 {" ".join(sys.argv[1:])}: {r}'
    print(line, flush=True)
    open(f'{H.LOGD}/r2.log', 'a').write(line + '\n')
