#!/usr/bin/env python
"""Stage 3 of the literal recipe: stack each group, cut templates, run the network matched filter
at 8 x MAD over the ETS days, and iterate the stack.

Literal elements (see results/REPRO_BOSTOCK.md sec.0):
  * the group's member windows are already mutually aligned (the autocorrelation lag IS the
    relative time), so the first stack is a plain linear stack of a long window at every station;
  * fine alignment follows Bungum & Husebye (1971) as Bostock 2015 sec.2 describes: each member is
    re-shifted by correlation against the composite stack, iterated;
  * the template is cut about the maximum of the stack envelope AT EACH STATION (Plourde et al.
    2015's station-dependent taper), 12 s long starting 2 s before the peak;
  * detection = the network mean of the per-station normalised 3-C cc exceeding
    median + 8 x MAD of that day's own function;  strongest detection kept per 4 s;
  * the stack is rebuilt from the detections and the whole thing repeated (iterations).
No survival gate, no scramble ratio, no quiet-day calibration, no burst guard.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k_common as K                                                        # noqa: E402

LONG_PRE = 12.0        # long stack window starts this far before the node time
LONG_LEN = 44.0        # ... and is this long (covers the ~10 s network moveout at both ends)
TPL_LEN = 12.0         # Plourde et al. 2015: a short window of 12 s
TPL_PRE = 2.0          # cut 2 s before the stack's envelope peak at that station
ALIGN_LAG = 1.0        # Bungum-Husebye re-shift half-width (s)
ALIGN_ITERS = 2
NMAD = 8.0
DECL_S = 4.0
MAX_STACK = 200
MIN_PEAK_MAD = 3.0     # a station is admitted if its stack shows an arrival at all (data QC)


def envelope(x, fs, smooth_s=0.4):
    from scipy.signal import hilbert
    e = np.sqrt((np.abs(hilbert(x, axis=-1)) ** 2).sum(0))
    k = max(int(round(smooth_s * fs)), 1)
    return np.convolve(e, np.ones(k) / k, mode="same")


def build_templates(dc, stas, groups, fs=K.FS, log=True, align_rounds=2, pad_s=1.0):
    """groups = {g: node times}.  Returns {g: {sta: (tpl(3,n), off_s, peak_mad)}}.

    The autocorrelation row grid is 0.5 s, so a node time carries +-0.25 s of jitter, and a plain
    stack of 1-8 Hz waveforms with that jitter has coherent gain |sinc(f*0.5 s)| -- zero at 2, 4,
    6, 8 Hz.  Bostock 2015 sec.2 (after Bungum & Husebye 1971) fixes exactly this: every member is
    re-shifted by correlation against the composite stack of the previous iteration, repeatedly.
    Round 0 is the plain stack; each further round realigns every member on the NETWORK sum of its
    per-station cc against the current stack, then restacks."""
    nL = int(round(LONG_LEN * fs))
    nP = int(round(pad_s * fs))
    byd = {}
    for g, t in groups.items():
        for x in t:
            byd.setdefault(pd.Timestamp(float(x), unit="s").strftime("%Y-%m-%d"),
                           []).append((g, float(x)))
    ST = {}
    shifts = {}
    for r in range(align_rounds + 1):
        acc = {}
        cnt = {}
        for date in sorted(byd):
            dc.day(date)
            for g, x in byd[date]:
                seg = {}
                for s in stas:
                    c = dc.cut(s, x - LONG_PRE - pad_s, LONG_LEN + 2 * pad_s, date)
                    if c is None or (c ** 2).sum() <= 0:
                        continue
                    seg[s] = K.norm3(c)
                if not seg:
                    continue
                k = 0
                if r > 0 and g in ST:
                    ccs = np.zeros(2 * nP + 1)
                    nn = 0
                    for s, c in seg.items():
                        ref = ST[g].get(s)
                        if ref is None:
                            continue
                        rn = np.sqrt((ref ** 2).sum())
                        if rn <= 0:
                            continue
                        V = np.lib.stride_tricks.sliding_window_view(c, nL, axis=1)
                        V = V[:, :2 * nP + 1, :]                     # (3, 2nP+1, nL)
                        num = np.einsum("clt,ct->l", V, ref)
                        den = np.sqrt(np.einsum("clt,clt->l", V, V))
                        ccs += np.where(den > 0, num / (rn * np.maximum(den, 1e-30)), 0.0)
                        nn += 1
                    if nn:
                        k = int(np.argmax(ccs)) - nP
                for s, c in seg.items():
                    i0 = nP + k
                    if i0 < 0 or i0 + nL > c.shape[1]:
                        continue
                    A = acc.setdefault((g, s), np.zeros((3, nL)))
                    A += c[:, i0:i0 + nL]
                    cnt[(g, s)] = cnt.get((g, s), 0) + 1
                shifts[(g, x)] = k / fs
        ST = {}
        for (g, s), A in acc.items():
            if cnt[(g, s)] >= 2:
                ST.setdefault(g, {})[s] = A / cnt[(g, s)]
        if log:
            K.log(f"[tpl] align round {r}: {len(ST)} groups with stacks")
    out, info = {}, []
    nt = int(round(TPL_LEN * fs))
    for g in groups:
        T = {}
        for s, st in ST.get(g, {}).items():
            env = envelope(st, fs)
            k = int(np.argmax(env))
            i0 = k - int(TPL_PRE * fs)
            if i0 < 0 or i0 + nt > nL:
                continue
            med = float(np.median(env))
            mad = float(1.4826 * np.median(np.abs(env - med)))
            h = (env[k] - med) / mad if mad > 0 else 0.0
            if h < MIN_PEAK_MAD:
                continue
            T[s] = (st[:, i0:i0 + nt].copy(), i0 / fs - LONG_PRE, float(h))
        if len(T) >= 4:
            out[g] = T
        info.append(dict(grp=g, n_sta_stack=len(ST.get(g, {})), n_sta=len(T)))
    if log:
        K.log(f"[tpl] {len(out)}/{len(groups)} groups keep >=4 stations")
    return out, pd.DataFrame(info)


def matched_filter(days, stas, TPL, dates, fs=K.FS, nmad=NMAD, dev="cuda", batch=24):
    """Network matched filter over `dates`.  Returns DataFrame(grp, t, cc, nmad)."""
    import torch
    torch.backends.cuda.matmul.allow_tf32 = True
    nt = int(round(TPL_LEN * fs))
    gs = sorted(TPL)
    out = []
    for date in dates:
        D = days.day(date)
        have = [s for s in stas if s in D]
        if len(have) < 4:
            continue
        N = min(D[s][1].shape[1] for s in have)
        NV = N - nt + 1
        acc = torch.zeros((len(gs), NV), device=dev)
        cnt = torch.zeros((len(gs), 1), device=dev)
        Xf, DEN, T0 = {}, {}, {}
        for s in have:
            t0, X = D[s]
            X = torch.from_numpy(np.ascontiguousarray(X[:, :N], dtype=np.float32)).to(dev)
            Xf[s] = torch.fft.rfft(X, n=N)
            xd = X.double()
            cs = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev,
                                                     dtype=torch.float64), xd * xd], 1), 1)
            e = (cs[:, nt:N + 1] - cs[:, :NV]).sum(0)
            # DATA QC: a window inside a gap has ~zero energy; 1/sqrt(0) turned the normalised cc
            # into 1e13 and every "detection" on such a day was a gap edge.  Gaps give cc = 0.
            emed = e.median()
            DEN[s] = torch.where(e > 1e-4 * emed, e.sqrt(),
                                 torch.full_like(e, float("inf"))).float()
            T0[s] = float(t0)
            del X, xd, cs, e
        for b0 in range(0, len(gs), batch):
            sub = gs[b0:b0 + batch]
            for s in have:
                idx = [i for i, g in enumerate(sub) if s in TPL[g]]
                if not idx:
                    continue
                M = np.zeros((len(idx), 3, N), np.float32)
                offs, nrm = [], []
                for r, i in enumerate(idx):
                    tpl, off, _ = TPL[sub[i]][s]
                    M[r, :, :nt] = tpl
                    offs.append(off)
                    nrm.append(np.sqrt((tpl ** 2).sum()))
                Tt = torch.from_numpy(M).to(dev)
                Tf = torch.fft.rfft(Tt, n=N)
                P = (Xf[s][None] * Tf.conj()).sum(1)
                del Tf, Tt
                cc = torch.fft.irfft(P, n=N)[:, :NV]
                del P
                nrmt = torch.tensor(nrm, device=dev, dtype=torch.float32)[:, None]
                cc = cc / (nrmt * DEN[s][None]).clamp_min(1e-30)
                for r, i in enumerate(idx):
                    k = int(round(offs[r] * fs))
                    a = max(-k, 0)
                    b = min(NV - k, NV)
                    if b <= a:
                        continue
                    acc[b0 + i, a:b] += cc[r, a + k:b + k]
                    cnt[b0 + i] += 1
                del cc
            torch.cuda.empty_cache()
        acc = acc / cnt.clamp_min(1.0)
        med = acc.median(dim=1, keepdim=True).values
        mad = (acc - med).abs().median(dim=1, keepdim=True).values
        thr = med + nmad * mad.clamp_min(1e-9)
        dec = int(round(DECL_S * fs))
        hit = acc >= thr
        vals = torch.where(hit, acc, torch.full((), -1e9, device=dev))
        pooled = torch.nn.functional.max_pool1d(vals[None], kernel_size=2 * dec + 1,
                                                stride=1, padding=dec)[0]
        keep = hit & (vals >= pooled - 1e-9)
        ij = torch.nonzero(keep, as_tuple=False)
        if ij.numel():
            gi = ij[:, 0].cpu().numpy()
            ti = ij[:, 1].cpu().numpy()
            v = acc[ij[:, 0], ij[:, 1]].cpu().numpy()
            nmv = ((acc[ij[:, 0], ij[:, 1]] - med[ij[:, 0], 0]) /
                   mad[ij[:, 0], 0].clamp_min(1e-9)).cpu().numpy()
            t0ref = T0[have[0]]
            out.append(pd.DataFrame(dict(grp=[gs[i] for i in gi],
                                         t=t0ref + ti / fs, cc=v, nmad=nmv, date=date)))
        del acc, hit, vals, pooled, keep, ij, Xf, DEN
        torch.cuda.empty_cache()
        K.log(f"[mf] {date} {len(have)} sta -> {len(out[-1]) if out else 0}")
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(
        columns=["grp", "t", "cc", "nmad", "date"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nodes", default=f"{K.OUT}/g20k_nodes.csv")
    ap.add_argument("--summary", default=f"{K.OUT}/g20k_summary.csv")
    ap.add_argument("--max-groups", type=int, default=0)
    ap.add_argument("--iters", type=int, default=3)
    ap.add_argument("--dates", default="")
    ap.add_argument("--stas", default="")
    ap.add_argument("--nmad-first", type=float, default=NMAD,
                    help="threshold of the FIRST pass (Bostock iterated to convergence; a looser\n                    first pass gathers stack members, later passes use 8xMAD)")
    ap.add_argument("--out", default=f"{K.OUT}/grow")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    stas = a.stas.split(",") if a.stas else K.CLEAN11
    dates = a.dates.split(",") if a.dates else K.bed_days()
    NO = pd.read_csv(a.nodes)
    G = pd.read_csv(a.summary).sort_values(["n_nodes", "nmad_max"], ascending=False)
    gsel = [int(x) for x in (G.grp.head(a.max_groups) if a.max_groups else G.grp)]
    groups = {g: NO[NO.grp == g].t.to_numpy() for g in gsel}
    K.log(f"{len(groups)} groups, sizes {min(len(v) for v in groups.values())}-"
          f"{max(len(v) for v in groups.values())}")
    dc = K.DayCache(stas)
    for it in range(a.iters + 1):
        TPL, info = build_templates(dc, stas, groups)
        if not TPL:
            K.log("no templates")
            break
        np.savez_compressed(f"{a.out}/tpl_it{it}.npz",
                            **{f"{g}_{s}": TPL[g][s][0].astype(np.float32)
                               for g in TPL for s in TPL[g]},
                            meta=np.array(json.dumps(
                                {f"{g}_{s}": [TPL[g][s][1], TPL[g][s][2]]
                                 for g in TPL for s in TPL[g]})))
        D = matched_filter(dc, stas, TPL, dates,
                           nmad=(a.nmad_first if it == 0 else NMAD))
        D.to_parquet(f"{a.out}/det_it{it}.parquet", index=False)
        n = D.groupby("grp").size()
        K.log(f"[iter {it}] {len(TPL)} templates -> {len(D)} detections; "
              f"per group median {n.median():.0f} max {n.max()}")
        if it == a.iters:
            break
        groups = {}
        for g, gg in D.groupby("grp"):
            t = gg.nlargest(MAX_STACK, "cc").t.to_numpy()
            if len(t) >= 3:
                groups[int(g)] = np.sort(t)
    K.log("done")


if __name__ == "__main__":
    main()
