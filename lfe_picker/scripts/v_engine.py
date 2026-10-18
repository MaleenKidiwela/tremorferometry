#!/usr/bin/env python
"""P4-0b -- the GPU matched-filter engine (one day at a time).  Prefix `v_`.

Implements the network statistic described in `v_common`'s docstring.  It is the frozen runner's
measurement core (`o_run_month.py`: rFFT full-day normalised correlation, FLOAT64 prefix-sum
denominator, station series = mean over 3 channels, per-(template, day) median/MAD network
normalisation) with Merlin's three P4 corrections:

  * per-station-day lowhalf MAD z BEFORE the network combination (loud-day proof);
  * per-SITE mean of the member stations, then a hard CLIP at z = 8 (anti-dominance);
  * EQUAL-WEIGHT sum over the sites observable at that sample (no SNR weights), with the
    coincidence clause >= 2 sites at per-site z >= 3 carried as its own accumulator.

Nothing in `o_*`/`s_*`/`g_*`/`p_*` is imported for its side effects and nothing is edited.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v_common as V                                                        # noqa: E402
import s_common as S                                                        # noqa: E402

FS = V.FS


def wait_gpu(need_gb=12.0, timeout_s=7200.0):
    import torch
    if not torch.cuda.is_available():
        return
    t0 = said = 0
    t0 = time.time()
    while True:
        torch.cuda.empty_cache()
        free = torch.cuda.mem_get_info()[0] / 2 ** 30
        if free >= need_gb or time.time() - t0 > timeout_s:
            return
        if not said:
            V.log(f"[gpu] only {free:.1f} GB free, waiting for {need_gb} GB ...")
            said = 1
        time.sleep(20)


class DayData:
    """One day's 20 Hz ZNE traces for a station set, pre-FFT'd, with the rolling-energy denominator."""

    def __init__(self, date, stas, n, dev="cuda", idx=None):
        import torch
        self.dev = dev
        self.n = int(n)
        self.stas = []
        self.t0 = {}
        self.DH = {}
        self.E = {}          # sqrt of the rolling 3-C-per-channel energy, (3, NV) float32
        self.valid = {}      # (NV,) bool
        self.NV = {}
        y, doy = int(date.year), int(date.dayofyear)
        for s in stas:
            v = S.load_day_zne(s, y, doy, idx)
            if v is None:
                continue
            t0, x = v
            N = x.shape[1]
            if N < 20 * self.n:
                continue
            NV = N - self.n
            NF = int(2 ** np.ceil(np.log2(N + self.n)))
            X = torch.from_numpy(np.ascontiguousarray(x)).to(dev)
            self.DH[s] = torch.fft.rfft(X, n=NF)
            xd = X.double()
            cs = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64),
                                         xd * xd], 1), 1)
            ew = (cs[:, self.n:N + 1] - cs[:, :NV + 1])[:, :NV]
            e3 = ew.sum(0)
            self.valid[s] = ((e3 > 1e-6 * e3.median()) & (ew.min(0).values > 0))
            floor = 1e-10 * ew.median()
            self.E[s] = ew.clamp_min(floor).sqrt().float()
            self.NV[s] = NV
            self.t0[s] = t0
            self.NF = NF
            self.stas.append(s)
            del X, xd, cs, ew, e3
        self.sites = sorted({V.SITE_OF[s] for s in self.stas})

    def station_z(self, s, W):
        """(B, NV) lowhalf-MAD z of the per-station 3-C mean normalised cc.  W is (B,3,n) torch."""
        import torch
        NF, NV = self.NF, self.NV[s]
        WH = torch.fft.rfft(W, n=NF)
        num = torch.fft.irfft(WH.conj() * self.DH[s][None], n=NF)[:, :, :NV]
        cc = (num / self.E[s][None]).mean(1)
        del WH, num
        m = cc.median(dim=-1, keepdim=True).values
        q = torch.kthvalue(cc, max(1, NV // 4), dim=-1, keepdim=True).values
        sc = (1.4826 * (m - q)).clamp_min(1e-9)
        return (cc - m) / sc

    def free(self):
        import torch
        self.DH.clear(); self.E.clear(); self.valid.clear()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


def run_day(date, stas, W, off, thr, dd=None, idx=None, block=48, dev="cuda",
            floor=V.T_FLOOR, max_peaks=400):
    """Grow every template of a bank over ONE day.

    W    (B, S, 3, n) float32, nan for a (template, station) that has no template
    off  (B, S) seconds: template start MINUS origin time
    thr  (B,) per-template net_mad threshold (use `floor` where unknown)

    Returns a DataFrame tpl, date, t (ORIGIN axis, epoch s), net_mad, n_site_z3, n_sites.
    Peaks are declustered at V.DECLUSTER_S by the greedy max rule."""
    import torch
    B, NS = W.shape[0], W.shape[1]
    n = W.shape[3]
    dd = DayData(date, stas, n, dev=dev, idx=idx) if dd is None else dd
    if len(dd.sites) < V.MIN_SITES:
        return pd.DataFrame(columns=["tpl", "date", "t", "net_mad", "n_site_z3", "n_sites"]), dd
    have = {s: i for i, s in enumerate(stas)}
    A = [s for s in dd.stas if s in have]
    t0max = max(dd.t0[s] for s in A)
    oo = off[:, [have[s] for s in A]]
    if not np.isfinite(oo).any():
        return pd.DataFrame(columns=["tpl", "date", "t", "net_mad", "n_site_z3", "n_sites"]), dd
    offmin = float(np.nanmin(oo))
    tg0 = t0max + max(0.0, -offmin) + 1.0
    K = {}
    NG = 10 ** 9
    for s in A:
        o = off[:, have[s]]
        k = np.rint((tg0 + np.where(np.isfinite(o), o, 0.0) - dd.t0[s]) * FS).astype(np.int64)
        k = np.where(np.isfinite(o), k, 0)
        k = np.clip(k, 0, max(dd.NV[s] - 2, 0))
        K[s] = k
        NG = min(NG, int(dd.NV[s] - k.max() - 1))
    NG = int(min(NG, 86400 * FS))
    if NG <= 1000:
        return pd.DataFrame(columns=["tpl", "date", "t", "net_mad", "n_site_z3", "n_sites"]), dd
    thr = np.asarray(thr, float)
    rows = []
    b0 = 0
    blk = int(block)
    while b0 < B:
        bs = list(range(b0, min(b0 + blk, B)))
        nb = len(bs)
        try:
            acc = torch.zeros((nb, NG), device=dev)
            nz3 = torch.zeros((nb, NG), device=dev)
            nsv = torch.zeros((nb, NG), device=dev)
            for q in dd.sites:
                mem = [s for s in V.SITE_MEMBERS[q] if s in A]
                sbuf = torch.zeros((nb, NG), device=dev)
                cnt = torch.zeros((nb, NG), device=dev)
                for s in mem:
                    si = have[s]
                    wblk = W[bs, si]
                    ok = np.isfinite(wblk).all(axis=(1, 2))
                    if not ok.any():
                        continue
                    ridx = np.nonzero(ok)[0]
                    Wt = torch.from_numpy(np.ascontiguousarray(wblk[ridx])).to(dev)
                    Z = dd.station_z(s, Wt)
                    del Wt
                    vmask = dd.valid[s].float()
                    kk = K[s][bs][ridx]
                    for kv in np.unique(kk):
                        rr = torch.from_numpy(np.nonzero(kk == kv)[0]).to(dev)
                        gr = torch.from_numpy(ridx[(kk == kv)]).to(dev)
                        vm = vmask[int(kv):int(kv) + NG]
                        sbuf.index_add_(0, gr, Z.index_select(0, rr)[:, int(kv):int(kv) + NG] * vm[None])
                        cnt.index_add_(0, gr, vm[None].expand(len(gr), NG).contiguous())
                    del Z
                okc = cnt > 0
                sz = torch.where(okc, sbuf / cnt.clamp_min(1.0), torch.zeros((), device=dev))
                sz = sz.clamp(max=V.SITE_Z_CLIP)
                acc += torch.where(okc, sz, torch.zeros((), device=dev))
                nz3 += ((sz >= V.SITE_Z_COINC) & okc).float()
                nsv += okc.float()
                del sbuf, cnt, okc, sz
        except torch.OutOfMemoryError:
            torch.cuda.empty_cache()
            if blk == 1:
                raise
            blk = max(1, blk // 2)
            V.log(f"[oom] block -> {blk}")
            continue
        m = acc.median(dim=-1, keepdim=True).values
        qq = torch.kthvalue(acc, max(1, NG // 4), dim=-1, keepdim=True).values
        sc = (1.4826 * (m - qq)).clamp_min(1e-9)
        nm = (acc - m) / sc
        tb = torch.from_numpy(np.maximum(thr[bs], floor)).to(dev).float()[:, None]
        hit = (nm >= tb) & (nz3 >= 2) & (nsv >= V.MIN_SITES)
        ij = torch.nonzero(hit, as_tuple=False)
        if len(ij):
            ii = ij[:, 0]; jj = ij[:, 1]
            vals = nm[ii, jj].cpu().numpy()
            nz = nz3[ii, jj].cpu().numpy()
            ns = nsv[ii, jj].cpu().numpy()
            ii = ii.cpu().numpy(); jj = jj.cpu().numpy()
            for r in np.unique(ii):
                sel = ii == r
                t = tg0 + jj[sel] / FS
                v = vals[sel]
                o = np.argsort(-v)[:20000]      # guard: a pathological template cannot stall
                kept_t, kept = [], []
                for p in o:
                    if all(abs(t[p] - u) >= V.DECLUSTER_S for u in kept_t):
                        kept_t.append(t[p]); kept.append(p)
                        if len(kept) >= max_peaks:
                            break
                kept = np.asarray(kept, int)
                rows.append(pd.DataFrame(dict(tpl=bs[int(r)], date=date.strftime("%Y-%m-%d"),
                                              t=t[kept], net_mad=v[kept],
                                              n_site_z3=nz[sel][kept], n_sites=ns[sel][kept])))
        del acc, nz3, nsv, nm, hit, ij
        b0 += nb
        if dev == "cuda":
            torch.cuda.empty_cache()
    out = (pd.concat(rows, ignore_index=True) if rows else
           pd.DataFrame(columns=["tpl", "date", "t", "net_mad", "n_site_z3", "n_sites"]))
    return out, dd


def grow(bank_path, days, out_path, thr=None, block=48, dev="cuda", floor=V.T_FLOOR,
         resume=True, tag=""):
    """Run a bank over a list of days (list of pandas Timestamps), writing one parquet.

    Idempotent: per-day part files under `<out_path>.parts/`; a day already written is skipped."""
    tids, stas, W, off, meta = V.load_bank(bank_path)
    B = len(tids)
    thr = np.full(B, floor) if thr is None else np.asarray(thr, float)
    parts = out_path + ".parts"
    os.makedirs(parts, exist_ok=True)
    idx = V.cache_idx()
    t00 = time.time()
    for i, d in enumerate(days):
        p = f"{parts}/{d.strftime('%Y%m%d')}.parquet"
        if resume and os.path.exists(p):
            continue
        wait_gpu()
        avail = [s for s in stas if (s, int(d.year), int(d.dayofyear)) in idx]
        if len({V.SITE_OF[s] for s in avail}) < V.MIN_SITES:
            pd.DataFrame(columns=["tpl", "date", "t", "net_mad", "n_site_z3", "n_sites"]).to_parquet(p)
            continue
        H, dd = run_day(d, stas, W, off, thr, idx=idx, block=block, dev=dev, floor=floor)
        dd.free()
        H.to_parquet(p, index=False)
        V.log(f"[grow{tag}] {d.date()} {len(avail)} sta / "
              f"{len({V.SITE_OF[s] for s in avail})} sites -> {len(H)} hits "
              f"({i+1}/{len(days)}, {time.time()-t00:.0f}s)")
    fs = sorted(__import__("glob").glob(f"{parts}/*.parquet"))
    A = [pd.read_parquet(f) for f in fs]
    A = [x for x in A if len(x)]
    R = pd.concat(A, ignore_index=True) if A else pd.DataFrame(
        columns=["tpl", "date", "t", "net_mad", "n_site_z3", "n_sites"])
    R["tid"] = [tids[int(i)] for i in R.tpl] if len(R) else []
    R.to_parquet(out_path, index=False)
    V.log(f"[grow{tag}] {len(days)} days, {len(R)} hits -> {out_path}")
    return R
