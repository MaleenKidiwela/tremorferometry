#!/usr/bin/env python
"""ARM B -- Chestler & Creager's OWN family-finding chain, run as literally as this bed allows.
Prefix `o_`.

Arm A (the rest of the queue) keeps Bostock's blind all-pairs autocorrelation and only changes
the site treatment.  Arm B instead copies the Olympic paper's chain end to end
(quoted in results/WA_2010_ARRAY_RESULTS.md sec.0):

  seeds   per ARRAY BEAM, peaks of the 5-8 Hz envelope above 18 x MAD; a +-1 s boxcar series per
          array; the boxcar series STACKED over arrays; peaks of that stack (>= MIN_ARR arrays,
          2 s decluster) are the CANDIDATE LFEs.  This is Chestler sec.2.1 eq.[1]-[2] with the
          one deviation forced by the cache (5-8 Hz instead of 5-9 Hz).
  grow    a 15-s, 2-8 Hz, 3-component raw window at each candidate time on THREE STATIONS PER
          ARRAY plus the single stations, matched-filtered through the tremor days; the network
          combination is Chestler's **2nd-root stack**, the threshold **45 x MAD** of that day's
          own function; the best 150 windows are stacked into a template; 3 iterations.
  dedupe  R_com > 0.5 (common detections / detections of the smaller family) removes the family
          with the lower median cc -- Chestler sec.2.2, which took her 57 families to 43.

Location and screening are then the SAME n_ machinery arm A uses, because this project has no
3-D Olympic model and no P picks on LFE stacks (sec.0).

  O_BED=2010 python scripts/o_chestler.py --stage seeds --out <dir>
  O_BED=2010 python scripts/o_chestler.py --stage grow  --out <dir> --max-seeds 1500
  O_BED=2010 python scripts/o_chestler.py --stage dedupe --out <dir>
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
import o_common as O                                                        # noqa: E402

PEAK_MAD = 18.0        # Chestler: peaks larger than 18 x MAD on each array beam
BOX_S = 1.0            # ... boxcar B(t - t_l) = 1 for |t - t_l| < 1 s
MIN_ARR = 3            # a candidate must be seen by at least this many arrays
SEED_DECL_S = 2.0
TPL_LEN = 15.0         # Chestler: 15-second windows
BAND = (2.0, 8.0)      # Chestler: band-passed between 2 and 8 Hz
NROOT = 2.0            # ... nth-root stack where n = 2
NMAD = 45.0            # ... exceeds 45 times the median absolute deviation
MAX_STACK = 150        # ... up to 150 of the best-correlating windows
ITERS = 3              # ... after three iterations
DECL_S = 4.0
RCOM = 0.5             # Chestler sec.2.2
STA_PER_ARRAY = 3      # "using a set of 3-stations from each array"


def hp2(X, fs=20.0):
    """The cache is 1-8 Hz; Chestler's band is 2-8 Hz -- one 2 Hz high-pass, nothing else."""
    from scipy.signal import butter, sosfiltfilt
    return sosfiltfilt(butter(4, BAND[0] / (fs / 2), btype="high", output="sos"),
                       X, axis=-1).astype(np.float32)


# ---------------------------------------------------------------- seeds -----------------------
def stage_seeds(a):
    from scipy.signal import butter, sosfiltfilt, hilbert
    K = O.activate()
    A = sorted(O.arrays())
    dates = O.tremor_dates()
    fs = K.FS
    sos = butter(4, [5.0 / (fs / 2), 7.999 / (fs / 2)], btype="band", output="sos")
    dec = int(round(SEED_DECL_S * fs))
    rows = []
    for date in dates:
        stack = None
        n_arr = 0
        for site in A:
            v = K.load_day(site, date)
            if v is None:
                continue
            t0, X = v
            z = sosfiltfilt(sos, X[0].astype(np.float64))
            e = np.abs(hilbert(z))
            k = max(int(round(0.4 * fs)), 1)
            e = np.convolve(e, np.ones(k) / k, mode="same")
            med = float(np.median(e))
            mad = float(np.median(np.abs(e - med)))
            if mad <= 0:
                continue
            thr = med + PEAK_MAD * mad
            hit = np.where(e >= thr)[0]
            if not hit.size:
                continue
            # keep the strongest sample per 2 s
            keep = []
            i = 0
            while i < hit.size:
                j = i
                while j + 1 < hit.size and hit[j + 1] - hit[i] < dec:
                    j += 1
                seg = hit[i:j + 1]
                keep.append(seg[int(np.argmax(e[seg]))])
                i = j + 1
            b = np.zeros(X.shape[1], np.float32)
            w = int(round(BOX_S * fs))
            for p in keep:
                b[max(p - w, 0):p + w] = 1.0                    # Chestler's A_l = 1 option
            stack = b if stack is None else stack + b
            n_arr += 1
        if stack is None or n_arr < MIN_ARR:
            continue
        cand = np.where(stack >= MIN_ARR)[0]
        if not cand.size:
            continue
        i = 0
        t0d = pd.Timestamp(date).timestamp()
        while i < cand.size:
            j = i
            while j + 1 < cand.size and cand[j + 1] - cand[i] < dec:
                j += 1
            seg = cand[i:j + 1]
            p = seg[int(np.argmax(stack[seg]))]
            rows.append(dict(date=date, t=t0d + p / fs, n_arrays=float(stack[p]),
                             n_arrays_live=n_arr))
            i = j + 1
        O.log(f"[seeds] {date}: {n_arr} array beams -> {sum(1 for r in rows if r['date']==date)}"
              f" candidate LFEs")
    S = pd.DataFrame(rows)
    os.makedirs(a.out, exist_ok=True)
    S.to_csv(f"{a.out}/seeds.csv", index=False)
    O.log(f"[seeds] {len(S)} candidates on {S.date.nunique() if len(S) else 0} days "
          f"-> {a.out}/seeds.csv")


def chestler_sites():
    """Chestler's correlation set: 3 stations per array (the 3 closest to the array centroid)
    plus every single station."""
    import s_common as S
    ST = O.member_table()
    out = []
    for r in O.sites().itertuples():
        mem = str(r.members).split(",")
        if r.kind == "SINGLE":
            out += mem
        else:
            d = [float(S.epi_km(r.lat, r.lon, float(ST.lat[m]), float(ST.lon[m]))) for m in mem]
            out += [m for _, m in sorted(zip(d, mem))[:STA_PER_ARRAY]]
    return sorted(out)


# ---------------------------------------------------------------- nth-root matched filter -----
def mf_nroot(days, stas, TPL, dates, fs=20.0, nmad=NMAD, dev="cuda", batch=24, tpl_len=TPL_LEN):
    """k_grow.matched_filter with Chestler's 2nd-root network stack and 45 x MAD.  Every other
    line (FFT cc, gap QC, 4-s declustering, per-day median/MAD) is k_grow's."""
    import torch
    nt = int(round(tpl_len * fs))
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
            DEN[s] = torch.where(e > 1e-4 * e.median(), e.sqrt(),
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
                    n_ = min(tpl.shape[1], nt)
                    M[r, :, :n_] = tpl[:, :n_]
                    offs.append(off)
                    nrm.append(np.sqrt((tpl[:, :nt] ** 2).sum()))
                Tt = torch.from_numpy(M).to(dev)
                Tf = torch.fft.rfft(Tt, n=N)
                P = (Xf[s][None] * Tf.conj()).sum(1)
                del Tf, Tt
                cc = torch.fft.irfft(P, n=N)[:, :NV]
                del P
                nt_ = torch.tensor(nrm, device=dev, dtype=torch.float32)[:, None]
                cc = cc / (nt_ * DEN[s][None]).clamp_min(1e-30)
                cc = torch.sign(cc) * cc.abs().pow(1.0 / NROOT)        # nth-root stack, n = 2
                for r, i in enumerate(idx):
                    k = int(round(offs[r] * fs))
                    lo, hi = max(-k, 0), min(NV - k, NV)
                    if hi <= lo:
                        continue
                    acc[b0 + i, lo:hi] += cc[r, lo + k:hi + k]
                    cnt[b0 + i] += 1
                del cc
            torch.cuda.empty_cache()
        acc = acc / cnt.clamp_min(1.0)
        acc = torch.sign(acc) * acc.abs().pow(NROOT)                   # ... raised back
        med = acc.median(dim=1, keepdim=True).values
        mad = (acc - med).abs().median(dim=1, keepdim=True).values
        thr = med + nmad * mad.clamp_min(1e-12)
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
            v = acc[ij[:, 0], ij[:, 1]]
            nv = ((v - med[ij[:, 0], 0]) / mad[ij[:, 0], 0].clamp_min(1e-12)).cpu().numpy()
            out.append(pd.DataFrame(dict(grp=[gs[i] for i in gi],
                                         t=T0[have[0]] + ti / fs,
                                         cc=v.cpu().numpy(), nmad=nv, date=date)))
        O.log(f"[mf-nroot] {date} {len(have)} sta -> {len(out[-1]) if out else 0}")
        del acc, hit, vals, pooled, keep, ij, Xf, DEN
        torch.cuda.empty_cache()
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(
        columns=["grp", "t", "cc", "nmad", "date"])


class BandCache:
    """k_common.DayCache with Chestler's 2-8 Hz band."""

    def __init__(self, stas):
        import o_beam as OB
        self.stas, self.OB = list(stas), OB
        self.date, self.D = None, {}

    def day(self, date):
        if date != self.date:
            self.D = {}
            for s in self.stas:
                v = self.OB.load_member(s, date)
                if v is not None:
                    self.D[s] = (v[0], hp2(v[1].astype(np.float64)))
            self.date = date
        return self.D

    def cut(self, sta, t_abs, len_s, date=None):
        date = (pd.Timestamp(float(t_abs), unit="s").strftime("%Y-%m-%d") if date is None
                else date)
        v = self.day(date).get(sta)
        if v is None:
            return None
        t0, X = v
        n = int(round(len_s * 20.0))
        i0 = int(round((t_abs - t0) * 20.0))
        if i0 < 0 or i0 + n > X.shape[1]:
            return None
        seg = X[:, i0:i0 + n].astype(np.float64)
        return seg if np.isfinite(seg).all() else None


def stage_grow(a):
    import k_grow as KG
    O.activate()
    S = pd.read_csv(f"{a.out}/seeds.csv")
    if a.max_seeds and len(S) > a.max_seeds:
        S = S.sort_values(["n_arrays", "t"], ascending=[False, True]).head(a.max_seeds)
    S = S.sort_values("t").reset_index(drop=True)
    stas = chestler_sites()
    O.log(f"[grow] {len(S)} seeds, {len(stas)} correlation stations "
          f"(3 per array + singles): {stas}")
    dc = BandCache(stas)
    dates = O.tremor_dates()
    # iteration 0: the RAW 15-s window at the candidate time is the template (Chestler)
    TPL = {}
    for r in S.itertuples():
        date = pd.Timestamp(float(r.t), unit="s").strftime("%Y-%m-%d")
        T = {}
        for s in stas:
            c = dc.cut(s, float(r.t) - TPL_LEN / 2, TPL_LEN, date)
            if c is None or (c ** 2).sum() <= 0:
                continue
            T[s] = (c.astype(np.float32), -TPL_LEN / 2, 99.0)
        if len(T) >= 4:
            TPL[int(r.Index)] = T
    O.log(f"[grow] iteration 0: {len(TPL)} seed templates with >= 4 stations")
    groups = None
    for it in range(ITERS + 1):
        if it > 0:
            TPL, _ = KG.build_templates(dc, stas, groups, log=True)
            if not TPL:
                break
        np.savez_compressed(f"{a.out}/tpl_it{it}.npz",
                            **{f"{g}_{s}": TPL[g][s][0].astype(np.float32)
                               for g in TPL for s in TPL[g]},
                            meta=np.array(json.dumps(
                                {f"{g}_{s}": [float(TPL[g][s][1]), float(TPL[g][s][2])]
                                 for g in TPL for s in TPL[g]})))
        D = mf_nroot(dc, stas, TPL, dates,
                     tpl_len=(TPL_LEN if it == 0 else KG.TPL_LEN))
        D.to_parquet(f"{a.out}/det_it{it}.parquet", index=False)
        n = D.groupby("grp").size() if len(D) else pd.Series(dtype=int)
        O.log(f"[grow] iter {it}: {len(TPL)} templates -> {len(D)} detections; "
              f"per group median {n.median() if len(n) else 0}")
        if it == ITERS or not len(D):
            break
        groups = {int(g): np.sort(gg.nlargest(MAX_STACK, "cc").t.to_numpy())
                  for g, gg in D.groupby("grp") if len(gg) >= 3}
        if not groups:
            break
    O.log("[grow] done")


def stage_dedupe(a):
    """Chestler sec.2.2: R_com > 0.5 -> drop the family with the lower median cc."""
    it = a.iter
    D = pd.read_parquet(f"{a.out}/det_it{it}.parquet")
    if not len(D):
        pd.DataFrame(columns=["grp"]).to_csv(f"{a.out}/keep.csv", index=False)
        return
    D["k"] = np.round(D.t / DECL_S).astype(np.int64)
    sets = {int(g): set(gg.k) for g, gg in D.groupby("grp")}
    med = D.groupby("grp").cc.median().to_dict()
    order = sorted(sets, key=lambda g: (-med[g], -len(sets[g])))
    keep = []
    for g in order:
        drop = False
        for h in keep:
            inter = len(sets[g] & sets[h])
            if inter / max(min(len(sets[g]), len(sets[h])), 1) > RCOM:
                drop = True
                break
        if not drop:
            keep.append(g)
    pd.DataFrame(dict(grp=sorted(keep))).to_csv(f"{a.out}/keep.csv", index=False)
    json.dump(dict(iteration=it, n_families_before=len(sets), n_families_after=len(keep),
                   rcom=RCOM, n_detections=int(len(D))),
              open(f"{a.out}/dedupe.json", "w"), indent=1)
    O.log(f"[dedupe] {len(sets)} -> {len(keep)} families at R_com > {RCOM}")
    # prune the detection file and the template file to the kept families
    D[D.grp.isin(set(keep))].drop(columns=["k"]).to_parquet(f"{a.out}/det_final.parquet",
                                                            index=False)
    z = np.load(f"{a.out}/tpl_it{it}.npz", allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    ks = [k for k in meta if int(k.split("_", 1)[0]) in set(keep)]
    np.savez_compressed(f"{a.out}/tpl_final.npz", **{k: z[k] for k in ks},
                        meta=np.array(json.dumps({k: meta[k] for k in ks})))


def stage_fa(a):
    """The quiet-day false-alarm run for arm B, through the SAME 2nd-root / 45 x MAD filter the
    discovery used -- otherwise clause (a) would compare two different statistics."""
    import k_grow as KG
    O.activate()
    z = np.load(f"{a.out}/tpl_final.npz", allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    TPL = {}
    for k, v in meta.items():
        g, s = k.split("_", 1)
        TPL.setdefault(int(g), {})[s] = (np.asarray(z[k], np.float32), float(v[0]), float(v[1]))
    stas = chestler_sites()
    dc = BandCache(stas)
    D = mf_nroot(dc, stas, TPL, O.quiet_days(), tpl_len=KG.TPL_LEN)
    D.to_parquet(f"{a.out}/det_fa.parquet", index=False)
    O.log(f"[fa] {len(TPL)} templates x {len(O.quiet_days())} quiet days -> {len(D)} detections "
          f"({len(D) / max(len(TPL) * len(O.quiet_days()), 1):.3f} per template-day)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["seeds", "grow", "dedupe", "fa"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-seeds", type=int, default=1500)
    ap.add_argument("--iter", type=int, default=ITERS)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    {"seeds": stage_seeds, "grow": stage_grow, "dedupe": stage_dedupe,
     "fa": stage_fa}[a.stage](a)


if __name__ == "__main__":
    main()
