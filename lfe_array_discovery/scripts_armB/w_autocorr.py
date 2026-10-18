#!/usr/bin/env python
"""Stage 1 -- NETWORK AUTOCORRELATION (Brown 2008 / Bostock 2012), GPU.  Prefix `w_`.

Merlin's statistic, verbatim (PLAN.md sec.10, review Q1/(a)):

  * frozen chain, 20 Hz, 1-8 Hz, 3-C ZNE, one hour at a time;
  * 10-s windows at 0.5-s hop, cut at the SAME ABSOLUTE TIME at every station, NO MOVEOUT MODEL
    (same-source events share their inter-station moveout, so the network correlation peaks at a
    common lag without knowing the cell);
  * per station: the normalised 3-C mean cc of window i against the WHOLE HOUR at full 1-sample
    lag resolution, converted to a lower-half-MAD z over the hour, the median and the MAD taken
    over VALID samples only (prerequisite P3);
  * per site: the mean of the member stations' z, CLIPPED at 8;
  * network: EQUAL-WEIGHT sum over the sites observable at that lag, re-normalised over lag by its
    own median / lower-half MAD (`net_mad`);
  * |lag| < 10 s excluded;
  * hit = net_mad >= T AND >= 2 sites with site-z >= 3 AND >= 4 sites observable (P4);
  * peaks declustered at 10 s.

T is NOT set here: peaks are recorded above `--t-record` and the FA-only threshold is derived from
the quiet hours by w_cluster.threshold() (prerequisite P1 -- no 8.0 floor, engineering floor 5.0).

  python scripts/w_autocorr.py --hours quiet   --out results/step1/ac_quiet
  python scripts/w_autocorr.py --hours tremor  --out results/step1/ac_tremor --t-record 6.0
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402

FS = W.FS


def _robust(x, sub, dim=-1):
    """(median, 1.4826*(median - Q25)) of x along dim, estimated on every `sub`-th column."""
    import torch
    v = x[:, ::sub] if sub > 1 else x
    m = v.median(dim=dim, keepdim=True).values
    q = torch.kthvalue(v, max(1, v.shape[dim] // 4), dim=dim, keepdim=True).values
    return m, (1.4826 * (m - q)).clamp_min(1e-9)


def autocorr_hour(loader, date, hour, mem, of, order, t_record, dev="cuda",
                  row_block=1024, max_peaks=W.AC_MAX_PEAKS_PER_WINDOW, sub=W.AC_Z_SUB,
                  info=None):
    """One station-hour block -> DataFrame(date, hour, t_i, t_j, net_mad, n_site_z3, n_sites)."""
    import torch
    n = int(round(W.AC_WIN_S * FS))
    hop = int(round(W.AC_HOP_S * FS))
    minlag = int(round(W.AC_MIN_LAG_S * FS))
    dec = int(round(W.AC_DECLUSTER_S * FS))
    det = [s for q in order for s in mem[q]]

    D = {}
    t0 = None
    veto_health, veto_flat = [], []
    for s in det:
        # F1: a station-hour the health census vetoed (flat / H-V rms ratio / 1-8 Hz tone) never
        # enters the sum.  The hour still runs on the remaining stations.
        if W.is_vetoed(s, date, hour):
            veto_health.append(s)
            continue
        v = loader.hour(s, date, hour)
        if v is None:
            continue
        ts, x = v
        if t0 is None:
            t0 = ts
        # half a sample: DayLoader.hour already cuts every station at the same ABSOLUTE time, so
        # the residual is only the sub-sample grid offset of the frozen cache's per-file t0
        # (0 / 8.4 / 15 ms on this bed).  A strict equality here is what dropped whole sites.
        if abs(ts - t0) > 0.5 / FS or x.shape[1] < int(0.5 * W.AC_BLOCK_S * FS):
            continue
        D[s] = x
    if not D:
        return None
    NH = min(x.shape[1] for x in D.values())
    NV = NH - n + 1
    NW = (NH - n) // hop + 1
    have_sites = sorted({of[s] for s in D})
    if len(have_sites) < W.MIN_SITES:
        return None

    A, HT, E, VAL, RE = {}, {}, {}, {}, {}
    for s, x in list(D.items()):
        X = torch.from_numpy(np.ascontiguousarray(x[:, :NH])).to(dev)
        xd = X.double()
        cs = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64),
                                     xd * xd], 1), 1)
        ew = (cs[:, n:NH + 1] - cs[:, :NV])[:, :NV]
        e3 = ew.sum(0)
        emed = float(e3.median())
        # F6: a station-hour more than half flat has ew.median() == 0 -> E == 0 -> cc = inf ->
        # z = NaN, and `NaN * 0` defeats the valid mask, so the WHOLE hour recorded zero pairs
        # (26/299 tremor hours in 1c).  The station is dropped from the hour instead.
        if W.DROP_ZERO_MEDIAN and not (emed > 0):
            veto_flat.append(s)
            del D[s], X, xd, cs, ew, e3
            continue
        a = X.unfold(-1, n, hop)[:, :NW].contiguous()                       # (3, NW, n)
        r = a.pow(2).sum(-1, keepdim=True).sqrt()
        A[s] = torch.where(r > 0, a / r.clamp_min(1e-30), torch.zeros_like(a))
        HT[s] = [X[c].unfold(0, n, 1)[:NV].contiguous() for c in range(3)]  # (NV, n) each
        VAL[s] = ((e3 > 1e-6 * e3.median()) & (ew.min(0).values > 0))
        E[s] = ew.clamp_min(1e-10 * ew.median()).sqrt().float()
        if W.WIN_TRANSIENT_FACTOR > 0:
            RE[s] = (e3 / max(emed, 1e-30)).float()                         # (NV,) hour-normalised
        del X, xd, cs, ew, e3, a, r
    if not D:
        return None
    NW = min(NW, min(A[s].shape[1] for s in A))
    have_sites = sorted({of[s] for s in D})
    if len(have_sites) < W.MIN_SITES:
        return None
    vidx = {s: (torch.nonzero(VAL[s], as_tuple=False).squeeze(1)
                if int(VAL[s].sum()) >= 1000 else torch.arange(NV, device=dev)) for s in D}

    # ---- WINDOW-TRANSIENT VETO (identical in quiet and tremor hours).  Each site is normalised
    # by its OWN hour median; a 10-s window whose NETWORK-MEDIAN energy exceeds the factor is
    # excluded both as row i and as lag position j.
    badcol = None
    if W.WIN_TRANSIENT_FACTOR > 0 and len(RE) and len(have_sites) >= W.MIN_SITES:
        per_site = []
        for q in have_sites:
            mm = [RE[s] for s in mem[q] if s in RE]
            if mm:
                per_site.append(torch.stack(mm, 0).median(0).values if len(mm) > 1 else mm[0])
        if len(per_site) >= W.MIN_SITES:
            netE = torch.stack(per_site, 0).median(0).values
            badcol = (netE > W.WIN_TRANSIENT_FACTOR)
        del per_site
    RE.clear()

    rows = []
    jj_all = torch.arange(NV, device=dev, dtype=torch.int32)
    for b0 in range(0, NW, row_block):
        b1 = min(b0 + row_block, NW)
        nb = b1 - b0
        acc = torch.zeros((nb, NV), device=dev)
        nz3 = torch.zeros((nb, NV), device=dev)
        nsv = torch.zeros((nb, NV), device=dev)
        for q in have_sites:
            sbuf = torch.zeros((nb, NV), device=dev)
            cnt = torch.zeros((nb, NV), device=dev)
            for s in mem[q]:
                if s not in D:
                    continue
                cc = torch.zeros((nb, NV), device=dev)
                for c in range(3):
                    cc += (A[s][c][b0:b1] @ HT[s][c].T) / E[s][c][None]
                cc /= 3.0
                m, sc = _robust(cc.index_select(1, vidx[s]), sub)
                z = (cc - m) / sc
                del cc
                # F6: a non-finite z (dead channel, zero-energy window) is NEVER summed and NEVER
                # counted -- `NaN * 0` is NaN, which propagated into net_mad and voided the hour.
                fin = torch.isfinite(z)
                z = torch.where(fin, z, torch.zeros((), device=dev))
                vm = VAL[s].float()[None] * fin.float()
                sbuf += z * vm
                cnt += vm
                del z, fin
            okc = cnt > 0
            sz = torch.where(okc, sbuf / cnt.clamp_min(1.0), torch.zeros((), device=dev))
            sz = sz.clamp(max=W.AC_SITE_Z_CLIP)
            acc += torch.where(okc, sz, torch.zeros((), device=dev))
            nz3 += ((sz >= W.AC_SITE_Z_COINC) & okc).float()
            nsv += okc.float()
            del sbuf, cnt, okc, sz
        anyv = (nsv >= W.MIN_SITES)
        vcols = torch.nonzero(anyv.any(0), as_tuple=False).squeeze(1)
        if vcols.numel() < 1000:
            del acc, nz3, nsv, anyv
            continue
        m, sc = _robust(acc.index_select(1, vcols), sub)
        nm = (acc - m) / sc
        del acc
        ti = torch.arange(b0, b1, device=dev, dtype=torch.int32)[:, None] * hop
        lagmask = (jj_all[None, :] - ti).abs() >= minlag
        if badcol is not None:
            rowbad = badcol[ti[:, 0].long().clamp(max=NV - 1)]
            lagmask = lagmask & (~badcol)[None, :] & (~rowbad)[:, None]
            del rowbad
        hit = (nm >= t_record) & (nz3 >= W.AC_MIN_SITE_Z_SITES) & anyv & lagmask
        vals = torch.where(hit, nm, torch.full((), -1e9, device=dev))
        pooled = torch.nn.functional.max_pool1d(vals[None], kernel_size=2 * dec + 1,
                                                stride=1, padding=dec)[0]
        keep = hit & (vals >= pooled - 1e-6)
        ij = torch.nonzero(keep, as_tuple=False)
        if ij.numel():
            ii, jv = ij[:, 0], ij[:, 1]
            v = nm[ii, jv]
            if ij.shape[0] > nb * max_peaks:
                o = torch.argsort(v, descending=True)[: nb * max_peaks]
                ii, jv, v = ii[o], jv[o], v[o]
            rows.append(pd.DataFrame(dict(
                i=(ii + b0).cpu().numpy().astype(np.int32),
                j=jv.cpu().numpy().astype(np.int32),
                net_mad=v.cpu().numpy().astype(np.float32),
                n_site_z3=nz3[ii, jv].cpu().numpy().astype(np.int8),
                n_sites=nsv[ii, jv].cpu().numpy().astype(np.int8))))
        del nm, nz3, nsv, hit, vals, pooled, keep, ij, lagmask, anyv
        torch.cuda.empty_cache()
    for s in list(D):
        del A[s], HT[s], E[s], VAL[s]
    torch.cuda.empty_cache()
    if info is not None:
        info.update(date=date, hour=int(hour), n_hour_sites=len(have_sites),
                    n_sta=len(D), n_veto_health=len(veto_health), n_veto_flat=len(veto_flat),
                    veto_health=",".join(veto_health), veto_flat=",".join(veto_flat),
                    n_windows_vetoed=(int(badcol.sum()) if badcol is not None else 0),
                    n_windows=int(NV))
    if not rows:
        return pd.DataFrame(columns=["date", "hour", "t_i", "t_j", "net_mad", "n_site_z3",
                                     "n_sites"])
    R = pd.concat(rows, ignore_index=True)
    R["t_i"] = t0 + R.i.to_numpy(np.float64) * hop / FS
    R["t_j"] = t0 + R.j.to_numpy(np.float64) / FS
    R["date"] = date
    R["hour"] = int(hour)
    R["n_hour_sites"] = len(have_sites)
    R["max_grid_residual_s"] = float(max((abs(v) for k, v in loader.last_residual.items()
                                          if k in D), default=0.0))
    return R[["date", "hour", "t_i", "t_j", "net_mad", "n_site_z3", "n_sites", "n_hour_sites"]]


def run(hours, out, t_record, shift=None, limit=0, dev="cuda"):
    mem, of, order, det = W.sites()
    loader = W.DayLoader(shift=shift, maxn=len(det) + 2)
    parts = out + ".parts"
    os.makedirs(parts, exist_ok=True)
    H = hours if not limit else hours.head(limit)
    t00 = time.time()
    tim = []
    for k, r in enumerate(H.itertuples()):
        p = f"{parts}/{r.date}_{int(r.hour):02d}.parquet"
        if os.path.exists(p):
            continue
        t1 = time.time()
        inf = {}
        R = autocorr_hour(loader, r.date, int(r.hour), mem, of, order, t_record, dev=dev,
                          info=inf)
        if inf:
            pd.DataFrame([inf]).to_csv(f"{parts}/{r.date}_{int(r.hour):02d}_info.csv", index=False)
        if R is None:
            R = pd.DataFrame(columns=["date", "hour", "t_i", "t_j", "net_mad", "n_site_z3",
                                      "n_sites", "n_hour_sites"])
        R.to_parquet(p, index=False)
        tim.append(time.time() - t1)
        W.log(f"[ac] {r.date} {int(r.hour):02d}h -> {len(R)} pairs "
              f"({tim[-1]:.1f}s; {k+1}/{len(H)}; {time.time()-t00:.0f}s total)")
    ifs = sorted(__import__("glob").glob(f"{parts}/*_info.csv"))
    if ifs:
        pd.concat([pd.read_csv(f) for f in ifs], ignore_index=True).to_csv(
            out + "_hourinfo.csv", index=False)
    fs = sorted(__import__("glob").glob(f"{parts}/*.parquet"))
    A = [pd.read_parquet(f) for f in fs]
    A = [x for x in A if len(x)]
    R = (pd.concat(A, ignore_index=True) if A else
         pd.DataFrame(columns=["date", "hour", "t_i", "t_j", "net_mad", "n_site_z3", "n_sites",
                               "n_hour_sites"]))
    R.to_parquet(out + ".parquet", index=False)
    json.dump(dict(n_hours=int(len(fs)), n_pairs=int(len(R)), t_record=float(t_record),
                   repair_flags=W.repair_flags(),
                   median_s_per_hour=float(np.median(tim)) if tim else None,
                   total_s=float(time.time() - t00), shift=shift or {}),
              open(out + "_meta.json", "w"), indent=1)
    W.log(f"[ac] {len(fs)} hours, {len(R)} pairs -> {out}.parquet")
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hours", default="tremor", choices=["tremor", "quiet", "file"])
    ap.add_argument("--hours-file", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--t-record", type=float, default=W.T_ENG_FLOOR)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--shift", default="")
    a = ap.parse_args()
    if a.hours == "quiet":
        H, nd = W.quiet_hours()
    elif a.hours == "tremor":
        H = W.candidate_hours()
        if W.AC_HOUR_BUDGET:
            H = H.sort_values("n_epi", ascending=False).head(W.AC_HOUR_BUDGET)
            H = H.sort_values(["date", "hour"]).reset_index(drop=True)
    else:
        H = pd.read_csv(a.hours_file)
    shift = json.loads(a.shift) if a.shift else None
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    run(H, a.out, a.t_record, shift=shift, limit=a.limit)


if __name__ == "__main__":
    main()
