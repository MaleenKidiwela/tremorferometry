#!/usr/bin/env python
"""H&G (2019) summed-CC matched filter, literal (PLAN.md rows 4-8) + our controls A1 (station-day-shift) and A2
(polarity-scrambled), computed from the same per-channel CC series.

Per channel c with template w (90 s from the origin, demeaned, unit norm):  cc_c(t) = sum_i w_i x(t+i) / den_c(t),
den_c(t) = ||x[t:t+90] - mean||  (hg_cont.py); cc_c = 0 and not live where den_c = 0.  Summed CC S(t) = sum over live
channels of the configuration; t = window start = origin time of a matched event.  Peaks: S(t) = max of S over
t +- 90 s (declustering, PLAN row 6); negative peaks likewise on -S.

Modes
  real  K      template K (1-5): sums 'real' (primary config = channels live at the template), 'snr' (V-snr: template
               RMS >= 3x RMS of the 30 min before the origin), A1_r (r = 0..7: each station's CC series rolled by its own
               whole-day offset), A2_r (r = 0..7: random +-1 per channel).  Writes data/mf/real_TK.npz and per-channel CC
               at the 40,000 highest real positive peaks.
  noise J0 J1  noise templates j in [J0, J1): the 27-channel 90-s window at noise time t_j (data/mf/noise_times.csv,
               made by 'noisesel') as a template, summed for every distinct configuration (primary and V-snr of T1-T5).
               Writes data/mf/noise_JJJ.npz (positive peaks above median + 3 MAD, self peak flagged, stats).
  noisesel     draws the 100 noise times (PLAN row 7) after the real runs exist.
usage: hg_mf.py real K | noise J0 J1 | noisesel"""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.signal import oaconvolve
from scipy.ndimage import maximum_filter1d
import hg_common as H

NDRAW = 8
SNR_MIN = 3.0


def cont(sta):
    return np.load(f'{H.MF}/cont_{sta}.npy', mmap_mode='r')


def den(sta):
    return np.load(f'{H.MF}/den_{sta}.npy', mmap_mode='r')


def cut_template(i0):
    """27-channel template at sample i0: dict ch -> (w unit-norm or None, snr)."""
    out = {}
    for si, s in enumerate(H.STAS):
        X, D = cont(s), den(s)
        for c in range(3):
            ch = f'{s}.{H.COMPS[c]}'
            w = np.array(X[c, i0:i0 + H.L], np.float64)
            if not np.all(np.isfinite(w)) or D[c, i0] <= 0:
                out[ch] = (None, np.nan)
                continue
            pre = np.array(X[c, i0 - 1800:i0], np.float64)
            pre = pre[np.isfinite(pre)]
            snr = np.sqrt(np.mean(w ** 2)) / np.sqrt(np.mean(pre ** 2)) if pre.size > 600 else np.nan
            w = w - w.mean()
            out[ch] = (w / np.linalg.norm(w), snr)
    return out


def channel_cc(ch, w):
    s, comp = ch.rsplit('.', 1)
    c = H.COMPS.index(comp)
    x = np.nan_to_num(np.asarray(cont(s)[c], np.float32))
    num = oaconvolve(x, w[::-1].astype(np.float32), mode='valid')
    del x
    d = np.asarray(den(s)[c])
    live = d > 0
    cc = np.zeros(H.M, np.float32)
    np.divide(num, d, out=cc, where=live)
    np.clip(cc, -1, 1, out=cc)
    return cc, live


def peaks(S, floor, sign=1):
    Z = S if sign > 0 else -S
    mf = maximum_filter1d(Z, 2 * H.L + 1, mode='constant', cval=-np.inf)
    i = np.where((Z == mf) & (Z > floor))[0]
    if i.size > 1:                                   # plateaus: keep the first of equal neighbours within 90 s
        keep = np.concatenate([[True], np.diff(i) > H.L])
        i = i[keep]
    return i.astype(np.int64), S[i].astype(np.float32)


def stats(S, NL):
    v = S[NL > 0]
    med = float(np.median(v)); mad = float(np.median(np.abs(v - med)))
    return med, mad


def draws():
    rng = np.random.default_rng(H.SEED)
    offs = []
    for r in range(NDRAW):
        while True:
            o = rng.choice(np.arange(30, H.NDAY - 30), size=len(H.STAS), replace=False)
            if np.min(np.diff(np.sort(o))) >= 2:
                break
        offs.append(o.astype(int))
    signs = []
    for r in range(NDRAW):
        while True:
            sg = rng.choice([-1.0, 1.0], size=len(H.CHANS))
            if 0 < (sg > 0).sum() < len(sg):
                break
        signs.append(sg)
    return np.array(offs), np.array(signs)


def config(T):
    prim = [ch for ch in H.CHANS if T[ch][0] is not None]
    snr = [ch for ch in prim if np.isfinite(T[ch][1]) and T[ch][1] >= SNR_MIN]
    return prim, snr


def run_real(K):
    TT = H.templates()
    row = TT[TT.template == K].iloc[0]
    i0 = H.tidx(row.time)
    T = cut_template(i0)
    prim, snr = config(T)
    offs, signs = draws()
    names = ['real', 'snr'] + [f'A1_{r}' for r in range(NDRAW)] + [f'A2_{r}' for r in range(NDRAW)]
    S = {n: np.zeros(H.M, np.float32) for n in names}
    NLn = ['real', 'snr'] + [f'A1_{r}' for r in range(NDRAW)]
    NL = {n: np.zeros(H.M, np.int8) for n in NLn}
    for ch in prim:
        ci = H.CHANS.index(ch); si = ci // 3
        cc, live = channel_cc(ch, T[ch][0])
        S['real'] += cc; NL['real'] += live
        if ch in snr:
            S['snr'] += cc; NL['snr'] += live
        for r in range(NDRAW):
            sh = int(offs[r, si]) * 86400
            S[f'A1_{r}'] += np.roll(cc, sh); NL[f'A1_{r}'] += np.roll(live, sh)
            if signs[r, ci] > 0:
                S[f'A2_{r}'] += cc
            else:
                S[f'A2_{r}'] -= cc
        del cc, live
    out = {'K': K, 'i0': i0, 'prim': np.array(prim), 'snr_chans': np.array(snr),
           'tsnr': np.array([T[ch][1] for ch in H.CHANS]), 'offs': offs, 'signs': signs}
    for n in names:
        nl = NL[n if n in NL else 'real']
        med, mad = stats(S[n], nl)
        ip, vp = peaks(S[n], 0.0, +1)
        ineg, vneg = peaks(S[n], 0.0, -1)
        out[f'{n}_pi'] = ip; out[f'{n}_pv'] = vp; out[f'{n}_pnl'] = nl[ip]
        out[f'{n}_ni'] = ineg; out[f'{n}_nv'] = vneg; out[f'{n}_nnl'] = nl[ineg]
        out[f'{n}_med'] = med; out[f'{n}_mad'] = mad
        out[f'{n}_livefrac'] = float((nl > 0).mean())
    # per-channel CC at the 40,000 highest real positive peaks (pass 2)
    top = out['real_pi'][np.argsort(out['real_pv'])[::-1][:40000]]
    top = np.sort(top)
    C = np.full((len(H.CHANS), top.size), np.nan, np.float32)
    for ch in prim:
        cc, live = channel_cc(ch, T[ch][0])
        v = cc[top]; v[~live[top]] = np.nan
        C[H.CHANS.index(ch)] = v
    out['top_i'] = top; out['top_cc'] = C
    del S, NL
    np.savez_compressed(f'{H.MF}/real_T{K}.part.npz', **out)
    os.replace(f'{H.MF}/real_T{K}.part.npz', f'{H.MF}/real_T{K}.npz')
    return f'T{K} prim {len(prim)} snr {len(snr)} real med {out["real_med"]:.3f} mad {out["real_mad"]:.3f}'


def noise_select():
    """100 noise times in the 2011 ETS: all 27 channels live; >= 180 s from the five template origins and from every
    positive real-run peak (primary) above median + 6 MAD of that template's real S."""
    rng = np.random.default_rng(H.SEED + 1)
    a, b = H.per_idx('ETS2011')
    TT = H.templates()
    excl = [H.tidx(t) for t in TT.time]
    for K in TT.template:
        z = np.load(f'{H.MF}/real_T{K}.npz')
        thr = float(z['real_med'] + 6 * z['real_mad'])
        ii = z['real_pi'][z['real_pv'] > thr]
        excl += list(ii[(ii > a - 1000) & (ii < b + 1000)])
    excl = np.array(sorted(excl))
    D = {s: den(s) for s in H.STAS}
    out = []
    tries = 0
    while len(out) < 100 and tries < 200000:
        tries += 1
        i = int(rng.integers(a, b - H.L))
        if excl.size and np.min(np.abs(excl - i)) < 180:
            continue
        if any(np.min(np.abs(np.array(out) - i)) < 180 for _ in [0] if out):
            continue
        if all(D[s][c, i] > 0 for s in H.STAS for c in range(3)):
            out.append(i)
    df = pd.DataFrame({'j': range(len(out)), 'i': out, 'time_utc': H.itime(out)})
    df.to_csv(f'{H.MF}/noise_times.csv', index=False)
    return f'noise times {len(df)} (tries {tries}, excluded peaks {excl.size})'


def all_configs():
    cfg = {}
    for K in H.templates().template:
        z = np.load(f'{H.MF}/real_T{K}.npz')
        cfg[f'T{K}_prim'] = tuple(z['prim']); cfg[f'T{K}_snr'] = tuple(z['snr_chans'])
    uniq = {}
    for k, v in cfg.items():
        uniq.setdefault(v, []).append(k)
    return uniq


def run_noise(J0, J1):
    NT = pd.read_csv(f'{H.MF}/noise_times.csv')
    uniq = all_configs()
    cfgs = list(uniq.items())
    for j in range(J0, J1):
        f = f'{H.MF}/noise_{j:03d}.npz'
        if os.path.isfile(f):
            continue
        i0 = int(NT.i[j])
        T = cut_template(i0)
        need = sorted({ch for v, _ in cfgs for ch in v})
        S = [np.zeros(H.M, np.float32) for _ in cfgs]
        NL = [np.zeros(H.M, np.int8) for _ in cfgs]
        for ch in need:
            if T[ch][0] is None:
                continue
            cc, live = channel_cc(ch, T[ch][0])
            for q, (v, _) in enumerate(cfgs):
                if ch in v:
                    S[q] += cc; NL[q] += live
            del cc, live
        out = {'j': j, 'i0': i0}
        for q, (v, names) in enumerate(cfgs):
            med, mad = stats(S[q], NL[q])
            ip, vp = peaks(S[q], med + 3 * mad, +1)
            tag = f'c{q}'
            out[f'{tag}_names'] = np.array(names); out[f'{tag}_n'] = len(v)
            out[f'{tag}_pi'] = ip; out[f'{tag}_pv'] = vp; out[f'{tag}_pnl'] = NL[q][ip]
            out[f'{tag}_med'] = med; out[f'{tag}_mad'] = mad
        out['ncfg'] = len(cfgs)
        del S, NL
        np.savez_compressed(f + '.part.npz', **out)
        os.replace(f + '.part.npz', f)
        print(f'{pd.Timestamp.now(tz="America/Los_Angeles"):%H:%M:%S} PT noise {j} done', flush=True)
    return f'noise {J0}-{J1} done'


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'real':
        r = run_real(int(sys.argv[2]))
    elif mode == 'noisesel':
        r = noise_select()
    elif mode == 'noise':
        r = run_noise(int(sys.argv[2]), int(sys.argv[3]))
    line = f'{pd.Timestamp.now(tz="America/Los_Angeles"):%Y-%m-%d %H:%M:%S} PT mf {" ".join(sys.argv[1:])}: {r}'
    print(line, flush=True)
    open(f'{H.LOGD}/mf.log', 'a').write(line + '\n')
