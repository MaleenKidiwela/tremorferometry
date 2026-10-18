#!/usr/bin/env python
"""STAGE 4 -- Royer & Bostock (2014) sec.3.2 WAVEFORM-CORRELATION CLUSTER ANALYSIS, literally.

  "We proceeded to combine all channels corresponding to a given detection into a single 'super'
   trace and cross-correlated all such traces against one another to determine maximum correlation
   coefficients and corresponding lags.  The correlation coefficients were used to populate a
   similarity matrix employed within a hierarchical cluster analysis, allowing the detection
   waveforms to be grouped into clusters.  Waveforms for all channels (stations/components)
   available for detections within a given cluster were shifted and stacked to produce an initial
   LFE template.  This procedure resulted in 224, 70 and 54 initial templates."

This REPLACES the union-find connected components of `k_group.py:51-65,96-99`, which chains: A~B and
B~C put A and C in one family even when cc(A,C) is 0, so it produces a few giant components instead
of 224 clusters.  Hierarchical clustering with a cut does not chain.

Stage 4a  the detection CATALOGUE.  sec.3.1's culled catalogue is a list of times, not of pairs
          (overlapping detections are excluded), so the window pairs from ab_auto are expanded to
          their constituent window times and the overlap cull is applied on the time axis.
          Reported against Royer's 4915 (NW) / 3306 (SW).
Stage 4b  the super traces, the all-pairs lagged cross-correlation, the similarity matrix.
Stage 4c  hierarchical clustering.  The paper is SILENT on linkage and cut, so BOTH are SWEPT and
          the whole sweep is a table; the headline configuration is chosen by the mean silhouette
          (a data-internal criterion that never looks at 224/70) and the cut whose cluster count is
          nearest the published one is reported beside it as a separate column, never substituted.
Stage 4d  shift-and-stack each cluster into an initial template.

  python scripts/ab_cluster.py --sub NW --auto results/ab_royer/NW/auto_window \
      --out results/ab_royer/NW/clusters
-> <out>/detections.parquet   t, cc, nmad, n_ch, snr, date, hour, n_pair
-> <out>/simmat.npz           S (Nd,Nd) float32 max cc, L (Nd,Nd) int16 lag samples, nch (Nd,Nd)
-> <out>/sweep.csv            linkage x cut -> n_clusters, n_clusters_ge{2,3,4,5}, silhouette
-> <out>/clusters.csv         one row per detection: cluster id under the headline configuration
-> <out>/tpl_it0.npz          initial templates, one (3, n) array per (cluster, station)
-> <out>/summary.json         counts beside Royer's published yields
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

LAG_MAX_S = 2.0          # half-width of the super-trace lag search (SILENT in the paper)
# 2026-09-28 (user-approved Merlin step 3): --lag-max / --lag-mode.  sec.3.2 "cross-correlated all such
# traces against one another to determine maximum correlation coefficients and corresponding lags"
# implies the FULL window.  Two lag modes:
#   data    (default, the old behaviour): the member's window slides into the CONTINUOUS record around it
#           (super traces cut win + 2*lag_max), every shifted window re-normalised.  Fine for small lags,
#           WRONG for |lag| ~ window: at 15 s it compares the reference with a different, adjacent window.
#   zeropad the literal cross-correlation of the two 15 s super traces: each channel's detection window is
#           unit-normalised ONCE, zero-padded by lag_max each side (linear, not circular), and the cc at lag
#           l is the overlap sum (partial overlap at large |l| lowers cc, as in a normalised xcorr of two
#           finite traces).  Channels never slide into one another (per-channel blocks, channel-mean cc).
#           Use with --lag-max 15 (= A.WIN_S) for the full window.
LAG_MODES = ["data", "zeropad"]
# SILENT in the paper.  AUDIT FIX A10: the old sweep stopped at distance 0.75 (similarity 0.25), but
# the network-mean cc of 15 s tremor windows lives at 0.1-0.3 (smoke test: pair cc median 0.13,
# similarity p99 0.10), so every configuration was all singletons.  The sweep now spans the whole
# distance axis 1 - S in [0.10, 0.98].
CUTS = np.round(np.arange(0.10, 0.981, 0.02), 3)
LINKAGES = ["average", "complete", "single", "weighted"]
MIN_MEMBERS = 2          # SILENT in the paper; the headline min size is chosen with the cut
MIN_SIZES = [2, 3, 4, 5, 6, 8, 10, 15, 20, 30]


# ------------------------------------------------------------------ 4a detection catalogue -----
def _strict_time_cull(P, w):
    """Greedy, strongest nmad first: keep a window PAIR only if neither of its windows lies within
    one window length of a window of an already kept pair (the pair-level 2-D suppression of
    ab_auto lets many pairs share one window in strong tremor)."""
    import bisect
    T, keep = [], []
    for r in P.sort_values("nmad", ascending=False).itertuples():
        ok = True
        for t in (r.t_i, r.t_j):
            k = bisect.bisect_left(T, t - w + 1e-9)
            if k < len(T) and T[k] < t + w - 1e-9:
                ok = False
                break
        if ok:
            keep.append(r.Index)
            bisect.insort(T, r.t_i)
            bisect.insort(T, r.t_j)
    return P.loc[keep]


def detections(auto_dir, decl_s=A.WIN_S, cull="strict"):
    """sec.3.1 "This catalogue was then culled to exclude overlapping detections".

    AUDIT FIX A14 (the SW 2.4x over-count).  A detection is a window PAIR; "overlapping" is read in
    TIME: with `cull="strict"` (headline) a pair is dropped when either of its windows overlaps a
    window of a stronger kept pair.  The old reading (`cull="pairnms"`: ab_auto's 2-D suppression,
    then a per-time cull) kept every pair that re-used a window with a partner far away, which is
    what inflated SW (strong 2007 tremor under the array: many partners per window).  Both counts
    are written to summary.json.  The catalogue is the two window times of every kept pair."""
    files = sorted(glob.glob(f"{auto_dir}/*.parquet"))
    if not files:
        raise SystemExit(f"no autocorrelation output in {auto_dir}")
    P = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    K = _strict_time_cull(P, decl_s)
    counts = dict(n_pairs_pairnms=int(len(P)), n_pairs_strict=int(len(K)))
    Q = K if cull == "strict" else P
    rows = pd.concat([
        Q[["t_i", "cc", "nmad", "n_ch", "snr_i", "date", "hour"]]
        .rename(columns={"t_i": "t", "snr_i": "snr"}),
        Q[["t_j", "cc", "nmad", "n_ch", "snr_j", "date", "hour"]]
        .rename(columns={"t_j": "t", "snr_j": "snr"})], ignore_index=True)
    rows["t"] = np.round(rows.t.to_numpy(float) * A.FS) / A.FS
    g = (rows.groupby("t")
             .agg(cc=("cc", "max"), nmad=("nmad", "max"), n_ch=("n_ch", "max"),
                  snr=("snr", "median"), date=("date", "first"), hour=("hour", "first"),
                  n_pair=("cc", "size")).reset_index().sort_values("t"))
    order = np.argsort(-g.nmad.to_numpy())
    t = g.t.to_numpy()
    taken = np.zeros(len(g), bool)
    keep = []
    for i in order:
        if taken[i]:
            continue
        keep.append(i)
        taken |= np.abs(t - t[i]) < decl_s
    D = g.iloc[sorted(keep)].reset_index(drop=True)
    counts["n_times"] = int(len(D))
    return D, counts


# ------------------------------------------------------------------ 4b super traces ------------
def super_traces(dc, stas, D, win_s=A.WIN_S, lag_max_s=LAG_MAX_S, lag_mode="data"):
    """(P (Nd, C, n+2L) float32 padded windows, V (Nd, C) bool channel validity, chans).
    lag_mode 'data': padding is real record; 'zeropad': only the window is cut, padding stays zero."""
    n = int(round(win_s * A.FS))
    L = int(round(lag_max_s * A.FS))
    chans = [(s, c) for s in stas for c in range(3)]
    Nd = len(D)
    P = np.zeros((Nd, len(chans), n + 2 * L), np.float32)
    V = np.zeros((Nd, len(chans)), bool)
    rc = A.RefCache(dc, n)                   # 2026-09-25 gap guard (was `(x ** 2).sum() <= 0`)
    for k, r in enumerate(D.itertuples()):
        for s in stas:
            if lag_mode == "zeropad":
                seg = dc.cut(s, r.t, win_s, r.date)
            else:
                seg = dc.cut(s, r.t - lag_max_s, win_s + 2 * lag_max_s, r.date)
            if seg is None:
                continue
            okc = A.segment_channel_ok(seg, rc.get(s, r.date), n)
            for c in range(3):
                x = seg[c]
                if not okc[c]:
                    continue
                if lag_mode == "zeropad":
                    P[k, chans.index((s, c)), L:L + n] = x[:n].astype(np.float32)
                else:
                    P[k, chans.index((s, c))] = x.astype(np.float32)
                V[k, chans.index((s, c))] = True
        if (k + 1) % 500 == 0:
            A.log(f"  [super] {k+1}/{Nd}")
    return P, V, chans


def simmat(P, V, win_s=A.WIN_S, lag_max_s=LAG_MAX_S, dev=None, block=512):
    """Similarity matrix: max over lag of the channel-mean normalised cc between super traces.

    For each lag the per-channel blocks are unit-normalised, so one GEMM over the concatenated
    channel axis IS the sum of per-channel normalised cc; the divisor is the number of channels
    valid in BOTH detections (a second small GEMM), never a constant."""
    import torch
    dev = dev or A.torch_dev()
    Nd, C, npad = P.shape
    n = int(round(win_s * A.FS))
    L = int(round(lag_max_s * A.FS))
    Vt = torch.from_numpy(V.astype(np.float32)).to(dev)
    NCH = Vt @ Vt.T                                         # (Nd, Nd) shared-channel count

    def unit(block_arr, i0):
        """(nb, C*n) unit-normalised per channel block at sample offset i0."""
        x = torch.from_numpy(np.ascontiguousarray(block_arr[:, :, i0:i0 + n])).to(dev)
        r = x.pow(2).sum(-1, keepdim=True).sqrt()
        x = torch.where(r > 0, x / r.clamp_min(1e-30), torch.zeros_like(x))
        return x.reshape(x.shape[0], -1)

    S = torch.full((Nd, Nd), -2.0, device=dev)
    LAG = torch.zeros((Nd, Nd), dtype=torch.int16, device=dev)
    ref = torch.cat([unit(P[b0:b0 + block], L) for b0 in range(0, Nd, block)], 0)  # (Nd, C*n)
    for l in range(-L, L + 1):
        M = torch.cat([unit(P[b0:b0 + block], L + l) for b0 in range(0, Nd, block)], 0)
        G = (ref @ M.T) / NCH.clamp_min(1.0)
        upd = G > S
        S = torch.where(upd, G, S)
        LAG = torch.where(upd, torch.full_like(LAG, l), LAG)
        del M, G, upd
    S = torch.where(NCH > 0, S, torch.zeros_like(S))
    S = 0.5 * (S + S.T)                                     # symmetrise; lags stay signed
    out = (S.cpu().numpy().astype(np.float32), LAG.cpu().numpy(),
           NCH.cpu().numpy().astype(np.int16))
    if dev == "cuda":
        torch.cuda.empty_cache()
    return out


def simmat_zeropad(P, V, win_s=A.WIN_S, lag_max_s=LAG_MAX_S, dev=None):
    """lag_mode 'zeropad': linear (non-circular) normalised cross-correlation of the super traces.

    P holds each channel's window at [L, L+n) and zeros elsewhere.  Each channel window is
    unit-normalised once; at lag l the reference (window part) is dotted with the member slice
    [L+l, L+l+n) -- member sample i+l against reference sample i, zeros outside the member's window --
    so G is the sum over shared channels of per-channel normalised cc, divided by the shared-channel
    count.  Same sign convention as simmat: LAG[ref, k] = l means member k's matching window starts
    at t_k + l/fs.  Exactly symmetric (cc_ij(l) = cc_ji(-l)); all arrays live on the GPU."""
    import torch
    dev = dev or A.torch_dev()
    Nd, C, npad = P.shape
    n = int(round(win_s * A.FS))
    L = int(round(lag_max_s * A.FS))
    assert npad == n + 2 * L
    Vt = torch.from_numpy(V.astype(np.float32)).to(dev)
    NCH = Vt @ Vt.T
    X = torch.from_numpy(P).to(dev)
    r = X[:, :, L:L + n].pow(2).sum(-1, keepdim=True).sqrt()
    X = torch.where(r > 0, X / r.clamp_min(1e-30), torch.zeros_like(X))
    ref = X[:, :, L:L + n].reshape(Nd, -1).contiguous()
    inv = 1.0 / NCH.clamp_min(1.0)
    S = torch.full((Nd, Nd), -2.0, device=dev)
    LAG = torch.zeros((Nd, Nd), dtype=torch.int16, device=dev)
    for l in range(-L, L + 1):
        M = X[:, :, L + l:L + l + n].reshape(Nd, -1)
        G = (ref @ M.T) * inv
        upd = G > S
        S = torch.where(upd, G, S)
        LAG = torch.where(upd, torch.full_like(LAG, l), LAG)
        del M, G, upd
        if (l + L) % 200 == 0:
            A.log(f"  [simmat zeropad] lag {l}/{L}")
    S = torch.where(NCH > 0, S, torch.zeros_like(S))
    S = 0.5 * (S + S.T)
    out = (S.cpu().numpy().astype(np.float32), LAG.cpu().numpy(),
           NCH.cpu().numpy().astype(np.int16))
    del X, ref, S, LAG
    if dev == "cuda":
        torch.cuda.empty_cache()
    return out


# ------------------------------------------------------------------ 4c clustering --------------
def cluster_sweep(S, cuts=CUTS, linkages=LINKAGES):
    """Hierarchical cluster analysis (Rowe et al. 2002) over a sweep of linkage and cut."""
    from scipy.cluster.hierarchy import fcluster, linkage as lk
    from scipy.spatial.distance import squareform
    Dm = 1.0 - np.clip(S, -1.0, 1.0)
    np.fill_diagonal(Dm, 0.0)
    Dm = 0.5 * (Dm + Dm.T)
    cond = squareform(Dm, checks=False)
    rows, labels = [], {}
    Dm32 = Dm.astype(np.float32)                 # silhouette memory (Nd^2) at half the cost
    for meth in linkages:
        Z = lk(cond, method=meth)
        for t in cuts:
            lab = fcluster(Z, t=float(t), criterion="distance")
            sizes = pd.Series(lab).value_counts()
            sil = silhouette(Dm32, lab)
            r = dict(linkage=meth, cut=float(t), n_clusters=int(sizes.size),
                     n_singleton=int((sizes == 1).sum()), max_size=int(sizes.max()),
                     silhouette=sil)
            for m in MIN_SIZES:
                r[f"n_ge{m}"] = int((sizes >= m).sum())
            rows.append(r)
            labels[(meth, float(t))] = lab
    return pd.DataFrame(rows), labels


def silhouette(Dm, lab):
    """Mean silhouette on a precomputed distance matrix; NaN if fewer than two usable clusters.

    AUDIT FIX A3: vectorised with a sparse one-hot GEMM (the per-point Python loop was O(Nd x K)
    per configuration).  Same definition as before: points in singleton clusters are skipped."""
    from scipy import sparse
    lab = np.asarray(lab)
    u, inv = np.unique(lab, return_inverse=True)
    if u.size < 2:
        return float("nan")
    n = lab.size
    sizes = np.bincount(inv).astype(float)
    H = sparse.csr_matrix((np.ones(n), (np.arange(n), inv)), shape=(n, u.size))
    Ssum = np.asarray((H.T @ Dm.T).T)                      # (n, K) sum of distances to each cluster
    own = sizes[inv]
    ok = own >= 2
    if not ok.any():
        return float("nan")
    a = Ssum[np.arange(n), inv] / np.maximum(own - 1, 1)
    M = Ssum / sizes[None, :]
    M[np.arange(n), inv] = np.inf
    b = M.min(1)
    den = np.maximum(a, b)
    sv = np.where(den > 0, (b - a) / np.where(den > 0, den, 1.0), np.nan)
    sv = sv[ok & (den > 0)]
    return float(np.mean(sv)) if sv.size else float("nan")


# ------------------------------------------------------------------ 4d stack -------------------
def stack_clusters(dc, stas, D, lab, LAG, S, min_members=MIN_MEMBERS, win_s=A.WIN_S,
                   long_pre=12.0, long_len=44.0):
    """Shift and stack every channel available in a cluster -> {cid: {sta: (3, nL)}}.

    The stack window is longer than the detection window so the template can later be cut (the
    network moveout across the subarray is several seconds).  Loops DATE-outer so each day of the
    25-station pool is loaded once (same sums as the cluster-outer loop)."""
    nL = int(round(long_len * A.FS))
    t = D.t.to_numpy(float)
    dates = D.date.to_numpy()
    jobs, info = [], []
    for cid in np.unique(lab):
        mem = np.flatnonzero(lab == cid)
        if mem.size < min_members:
            info.append(dict(cid=int(cid), n_mem=int(mem.size), status="too small"))
            continue
        # the medoid (highest mean similarity to the rest of the cluster) is the reference; every
        # member is then shifted by its lag against that reference before stacking
        ref = int(mem[int(np.argmax(S[np.ix_(mem, mem)].mean(1)))])
        for k in mem:
            # AUDIT FIX A2: LAG[ref, k] = l means member k's window matching ref starts at
            # t_k + l/fs (simmat compares ref at offset L with member k at offset L + l)
            jobs.append((str(dates[k]), int(cid), t[k] + float(LAG[ref, k]) / A.FS))
        info.append(dict(cid=int(cid), n_mem=int(mem.size), status="ok"))
    acc, cnt = {}, {}
    rc = A.RefCache(dc, int(round(win_s * A.FS)))     # 2026-09-25 gap guard: dead channels -> 0
    for date, cid, tk in sorted(jobs):
        for s in stas:
            seg = dc.cut(s, tk - long_pre, long_len, date)
            if seg is None:
                continue
            okc = rc.get(s, date) > 0
            if not okc.any():
                continue
            if not okc.all():
                seg = np.where(okc[:, None], seg, 0.0)
            r = np.sqrt((seg ** 2).mean())
            if r <= 0:
                continue
            key = (cid, s)
            if key not in acc:
                acc[key] = np.zeros((3, nL))
                cnt[key] = 0
            acc[key] += seg[:, :nL] / r
            cnt[key] += 1
    out = {}
    for (cid, s), v in acc.items():
        if cnt[(cid, s)] >= 2:
            out.setdefault(cid, {})[s] = v / cnt[(cid, s)]
    info = pd.DataFrame(info)
    info["n_sta"] = [len(out.get(c, {})) for c in info.cid]
    return out, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", required=True, choices=["NW", "SW", "NVI", "MENDO"])
    ap.add_argument("--auto", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--lag-max", type=float, default=LAG_MAX_S,
                    help="half-width (s) of the super-trace lag search; 15 = the full window")
    ap.add_argument("--lag-mode", default="data", choices=LAG_MODES,
                    help="data = slide into the continuous record (old); zeropad = literal xcorr of the "
                         "finite windows (use with --lag-max 15)")
    ap.add_argument("--simmat-only", action="store_true",
                    help="stop after writing simmat.npz + detections.parquet")
    ap.add_argument("--min-members", type=int, default=MIN_MEMBERS)
    ap.add_argument("--max-det", type=int, default=0, help="0 = all; else the strongest N by nmad")
    ap.add_argument("--cull", default="pairnms", choices=["strict", "pairnms"])   # v3: pair-overlap cull
    ap.add_argument("--sweep-only", action="store_true")
    ap.add_argument("--null-sweep", default="", help="sweep_long.csv of the null (headline rule)")
    ap.add_argument("--headline-from", default="", help="summary.json whose headline to copy")
    a = ap.parse_args()
    auto = a.auto or f"{A.RES}/{a.sub}/auto_window"
    out = a.out or f"{A.RES}/{a.sub}/clusters"
    os.makedirs(out, exist_ok=True)
    stas = A.SUBARRAY[a.sub]["stas"]
    pub = A.SUBARRAY[a.sub]

    D, counts = detections(auto, cull=a.cull)
    n_pairs = counts["n_pairs_strict"] if a.cull == "strict" else counts["n_pairs_pairnms"]
    if a.max_det and len(D) > a.max_det:
        D = D.nlargest(a.max_det, "nmad").sort_values("t").reset_index(drop=True)
    D.to_parquet(f"{out}/detections.parquet", index=False)
    A.log(f"[4a] pairs: {counts['n_pairs_pairnms']} (2-D pair suppression) / "
          f"{counts['n_pairs_strict']} (strict time-overlap cull, cull={a.cull}) -> "
          f"{len(D)} detection times (Royer: {pub['n_det_pub']} initial detections)")

    dc = A.DayCache(stas)
    P, V, chans = super_traces(dc, stas, D, lag_max_s=a.lag_max, lag_mode=a.lag_mode)
    A.log(f"[4b] super traces {P.shape} (lag_mode {a.lag_mode}, lag_max {a.lag_max} s); "
          f"mean channels per detection {V.sum(1).mean():.1f} / {len(chans)}")
    if a.lag_mode == "zeropad":
        S, LAG, NCH = simmat_zeropad(P, V, lag_max_s=a.lag_max)
    else:
        S, LAG, NCH = simmat(P, V, lag_max_s=a.lag_max)
    del P
    np.savez_compressed(f"{out}/simmat.npz", S=S, L=LAG, nch=NCH,
                        t=D.t.to_numpy(), chans=np.array([f"{s}.{c}" for s, c in chans]))
    A.log(f"[4b] similarity matrix: median {np.median(S):.3f}, p99 {np.percentile(S,99):.3f}")
    if a.simmat_only:
        return
    off = S[~np.eye(S.shape[0], dtype=bool)]
    s_med = float(np.median(off))
    s_mad = float(np.median(np.abs(off - s_med)))
    cut8 = float(np.round(1.0 - (s_med + A.NMAD * s_mad), 4))
    SW, labels = cluster_sweep(S, cuts=np.unique(np.r_[CUTS, cut8]))
    # long table: linkage x cut x min size, with the published count and the distance to it
    L = SW.melt(id_vars=["linkage", "cut", "n_clusters", "n_singleton", "max_size", "silhouette"],
                value_vars=[f"n_ge{m}" for m in MIN_SIZES], var_name="min_size",
                value_name="n_templates")
    L["min_size"] = L.min_size.str[4:].astype(int)
    L["published"] = pub["n_tpl0_pub"]
    L["d_to_published"] = (L.n_templates - pub["n_tpl0_pub"]).abs()
    L["nearest_published"] = L.d_to_published == L.d_to_published.min()
    SW.to_csv(f"{out}/sweep.csv", index=False)
    if a.sweep_only:
        L.to_csv(f"{out}/sweep_long.csv", index=False)
        json.dump(dict(subarray=a.sub, sweep_only=True, **counts,
                       similarity_median=s_med, similarity_mad=s_mad),
                  open(f"{out}/sweep_summary.json", "w"), indent=1)
        A.log(f"[4c] sweep only -> {out}/sweep_long.csv")
        return
    # HEADLINE (AUDIT FIX A15).  Linkage, cut and minimum cluster size are SILENT (sec.3.2; Rowe et
    # al. 2002 and Bostock et al. 2012 give no numbers).  Headline, in order of precedence:
    #   --headline-from <summary.json>  (the null MUST use the real run's configuration);
    #   --null-sweep <sweep_long.csv>   the configuration with the largest REAL-minus-NULL template
    #                                   count (the site-shift null of the same stage) -- a
    #                                   data-internal choice that never looks at 224/70;
    #   otherwise average linkage at the median + 8 MAD similarity cut, min size 2 (old A11).
    # HEADLINE (AUDIT FIX A16, supersedes A15).  Rowe et al. (2002) -- the method Royer cites --
    # is a "dendrogram-based, hierarchical PAIR-GROUP algorithm" in which "CC coefficients are used
    # to fuse event pairs with high waveform similarity into clusters, and clustering is considered
    # complete when a user-defined CC coefficient cutoff threshold is reached" (as summarised in
    # USGS PP1769 ch.5, which used Rowe's code).  Pair-group with arithmetic averaging = UPGMA =
    # scipy 'average'.  Royer does not give the cutoff: we use the paper's own significance rule,
    # similarity > median + 8 MAD (sec.3.1/3.3).  Rowe's clustering implies no minimum size and
    # Royer says each cluster was stacked into a template, so min size = 2 (a single detection is
    # not a cluster).  The null-1 sweep, the real-minus-null excess and the null/real ratio are
    # REPORTED per (linkage, cut, min size) -- A15's "max excess" choice was dropped: it always
    # favours pairs and its merge compared different cut grids (a bug that assigned null = 0).
    how = "Rowe pair-group (UPGMA) at median + 8 MAD, min size 2"
    if a.headline_from:
        h = json.load(open(a.headline_from))["headline"]
        hl, hm = h["linkage"], int(h["min_size"])
        hc = cut8 if h.get("cut_rule") == "8mad" else float(h["cut"])
        how = f"copied rule from {a.headline_from}"
    else:
        hl, hc, hm = "average", cut8, 2
    if a.null_sweep and os.path.exists(a.null_sweep):
        N = pd.read_csv(a.null_sweep)[["linkage", "cut", "min_size", "n_templates"]]
        N["cut"] = N.cut.round(2)
        N = N.rename(columns={"n_templates": "n_null"}).drop_duplicates(
            ["linkage", "cut", "min_size"])
        L["cut_r"] = L.cut.round(2)
        L = L.merge(N.rename(columns={"cut": "cut_r"}), on=["linkage", "cut_r", "min_size"],
                    how="left").drop(columns="cut_r")
        L["excess_over_null"] = L.n_templates - L.n_null
        L["null_over_real"] = L.n_null / L.n_templates.clip(lower=1)
    L.to_csv(f"{out}/sweep_long.csv", index=False)
    key = min(labels, key=lambda k: (k[0] != hl, abs(k[1] - hc)))
    lab = labels[key]
    hc = key[1]
    sizes = pd.Series(lab).value_counts()
    n_head = int((sizes >= hm).sum())
    near = L.sort_values("d_to_published").iloc[0]
    A.log(f"[4c] similarity median {s_med:.4f} MAD {s_mad:.4f} (8xMAD cut {cut8})")
    A.log(f"[4c] headline ({how}): {hl} cut {hc} min size {hm} -> {n_head} clusters "
          f"(Royer: {pub['n_tpl0_pub']}); nearest-published config (reported only): "
          f"{near.linkage} cut {near.cut} min size {near.min_size} -> {near.n_templates}")
    D2 = D.copy()
    D2["cluster"] = lab
    D2["cluster_nearpub"] = labels[min(labels, key=lambda k: (k[0] != near.linkage,
                                                               abs(k[1] - near.cut)))]
    D2.to_csv(f"{out}/clusters.csv", index=False)
    a.min_members = hm
    head = dict(linkage=hl, cut=hc, min_size=hm, n_clusters=n_head, how=how, cut_rule="8mad")

    # sec.3.2: "Waveforms for all channels (stations/components) available for detections within a
    # given cluster were shifted and stacked" -> stack over the whole 25-station pool (A.POOL)
    pool = list(A.POOL[a.sub])
    TPL, info = stack_clusters(A.DayCache(pool), pool, D, lab, LAG, S, min_members=a.min_members)
    info.to_csv(f"{out}/cluster_stacks.csv", index=False)
    np.savez_compressed(f"{out}/tpl_it0.npz",
                        **{f"{g}_{s}": TPL[g][s].astype(np.float32) for g in TPL for s in TPL[g]},
                        meta=np.array(json.dumps(dict(long_pre=12.0, long_len=44.0, fs=A.FS))))
    json.dump(dict(subarray=a.sub, n_pairs=int(n_pairs), n_detections=int(len(D)),
                   n_detections_published=pub["n_det_pub"],
                   headline=head, **counts, cull=a.cull,
                   similarity_median=s_med, similarity_mad=s_mad, cut_8mad=cut8,
                   nearest_published=dict(linkage=str(near.linkage), cut=float(near.cut),
                                          min_size=int(near.min_size),
                                          n_clusters=int(near.n_templates)),
                   n_initial_templates=int(len(TPL)),
                   n_initial_templates_published=pub["n_tpl0_pub"],
                   min_members=a.min_members, lag_max_s=a.lag_max,
                   lag_mode=a.lag_mode),
              open(f"{out}/summary.json", "w"), indent=1)
    A.log(f"[4d] {len(TPL)} initial templates with stacks (Royer: {pub['n_tpl0_pub']})")


if __name__ == "__main__":
    main()
