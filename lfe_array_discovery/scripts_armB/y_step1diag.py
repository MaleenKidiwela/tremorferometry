#!/usr/bin/env python
"""STEP-1 POST-MORTEM DIAGNOSTICS D4-D7 (agent G, 2026-09-20).  Prefix `y_`.

DIAGNOSTIC ONLY.  Nothing here changes a frozen verdict, and no engine constant is retuned:
every measurement either reads the frozen run's own output or re-measures it with the frozen
constants, except where a variant is explicitly labelled an ABLATION.

  D4  Why T = 12.50 on the 2005 bed.  Per quiet hour and per tremor hour: pairs at every T on the
      0.25 grid, largest connected component, |lag| distribution, and -- on the GPU -- the
      per-SITE z of the pairs the engine actually reported, so a site-coherent-noise or
      single-station explanation can be told from a real network one.
  D5  True-pair recoverability.  For KEPT families with >= 50 Bostock events in the bed window,
      same-family same-hour event pairs, measured EXACTLY as `w_autocorr.autocorr_hour` measures a
      pair: per-station 3-C cc of the 10-s window against the whole hour, lower-half-MAD z over
      valid lags, per-site mean clipped at 8, equal-weight site sum, per-row median/lower-half-MAD
      renormalisation -> net_mad at the true lag.
  D6  Ablations (window length, per-site clip, site weighting, MAD normalisation, band/sample
      rate) on the true pairs and on 4 quiet hours.
  D7  Hour-flag coverage: the fraction of Bostock's 2005 events that fall inside the hours the
      array-energy flag selected.

The 30 hidden families are excluded everywhere (`s_common.holdouts`).
"""
from __future__ import annotations

import argparse
import collections
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import w_cluster as C                                                       # noqa: E402
import s_common as S                                                        # noqa: E402

FS = W.FS
T_GRID = np.arange(W.T_ENG_FLOOR, 16.0 + 1e-9, W.T_GRID)


# =================================================================== shared GPU machinery =======
class Hour:
    """One prepared station-hour block, built EXACTLY as `w_autocorr.autocorr_hour` builds it."""

    def __init__(self, loader, date, hour, mem, of, order, win_s=W.AC_WIN_S, dev="cuda"):
        import torch
        self.ok = False
        self.date, self.hour, self.dev = date, int(hour), dev
        self.n = int(round(win_s * FS))
        self.mem, self.of, self.order = mem, of, order
        det = [s for q in order for s in mem[q]]
        D, t0 = {}, None
        for s in det:
            v = loader.hour(s, date, hour)
            if v is None:
                continue
            ts, x = v
            if t0 is None:
                t0 = ts
            if abs(ts - t0) > 0.5 / FS or x.shape[1] < int(0.5 * W.AC_BLOCK_S * FS):
                continue
            D[s] = x
        if not D:
            return
        self.t0 = t0
        self.NH = NH = min(x.shape[1] for x in D.values())
        self.NV = NV = NH - self.n + 1
        self.have = sorted({of[s] for s in D})
        self.D = D
        self.HT, self.E, self.VAL, self.vidx = {}, {}, {}, {}
        for s, x in D.items():
            X = torch.from_numpy(np.ascontiguousarray(x[:, :NH])).to(dev)
            self.HT[s] = [X[c].unfold(0, self.n, 1)[:NV].contiguous() for c in range(3)]
            xd = X.double()
            cs = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64),
                                         xd * xd], 1), 1)
            ew = (cs[:, self.n:NH + 1] - cs[:, :NV])[:, :NV]
            e3 = ew.sum(0)
            self.VAL[s] = ((e3 > 1e-6 * e3.median()) & (ew.min(0).values > 0))
            self.E[s] = ew.clamp_min(1e-10 * ew.median()).sqrt().float()
            self.vidx[s] = (torch.nonzero(self.VAL[s], as_tuple=False).squeeze(1)
                            if int(self.VAL[s].sum()) >= 1000
                            else torch.arange(NV, device=dev))
            del X, xd, cs, ew, e3
        self.ok = len(self.have) >= W.MIN_SITES

    def close(self):
        import torch
        for s in list(getattr(self, "HT", {})):
            del self.HT[s], self.E[s], self.VAL[s], self.vidx[s]
        self.HT = self.E = self.VAL = self.vidx = self.D = None
        torch.cuda.empty_cache()

    # -------------------------------------------------------------- the statistic ---------------
    def _robust(self, x, sub=W.AC_Z_SUB, mad="lowhalf"):
        import torch
        v = x[..., ::sub] if sub > 1 else x
        m = v.median(dim=-1, keepdim=True).values
        if mad == "full":
            sc = 1.4826 * (v - m).abs().median(dim=-1, keepdim=True).values
        else:
            q = torch.kthvalue(v, max(1, v.shape[-1] // 4), dim=-1, keepdim=True).values
            sc = 1.4826 * (m - q)
        return m, sc.clamp_min(1e-9)

    def station_z(self, sta, i0, sub=W.AC_Z_SUB, mad="lowhalf"):
        """(NV,) lower-half-MAD z of one window of one station against its whole hour."""
        import torch
        a = torch.stack([self.HT[sta][c][i0] for c in range(3)], 0)[:, None, :]   # (3,1,n)
        r = a.pow(2).sum(-1, keepdim=True).sqrt()
        a = torch.where(r > 0, a / r.clamp_min(1e-30), torch.zeros_like(a))
        cc = torch.zeros((1, self.NV), device=self.dev)
        for c in range(3):
            cc += (a[c] @ self.HT[sta][c].T) / self.E[sta][c][None]
        cc /= 3.0
        m, sc = self._robust(cc.index_select(1, self.vidx[sta]), sub, mad)
        return ((cc - m) / sc)[0]

    def moveout_rows(self, i0_by_sta, clip=W.AC_SITE_Z_CLIP, sub=W.AC_Z_SUB, mad="lowhalf",
                     weights=None):
        """ORACLE MOVEOUT-AWARE pairing (DIAGNOSTIC, not the frozen statistic): every station is
        cut at its OWN predicted S instead of at one common absolute time, and the per-site z are
        summed on the common LAG axis (inter-event interval), where a same-source pair aligns at
        every station by construction.  Returns (net_mad (nd,), site_z {site:(nd,)}, nz3, nsv,
        d0) where lag index d = d0 + k samples."""
        import torch
        Zs = {}
        for s, i0 in i0_by_sta.items():
            if i0 < 0 or i0 + self.n > self.NH:
                continue
            Zs[s] = self.station_z(s, int(i0), sub, mad)
        if not Zs:
            return None
        lo = -min(int(i0_by_sta[s]) for s in Zs)
        hi = self.NV - max(int(i0_by_sta[s]) for s in Zs)
        if hi - lo < 1000:
            return None
        nd = hi - lo
        acc = torch.zeros(nd, device=self.dev)
        nz3 = torch.zeros(nd, device=self.dev)
        nsv = torch.zeros(nd, device=self.dev)
        SZ = {}
        for q in self.order:
            mm = [s for s in self.mem[q] if s in Zs]
            if not mm:
                continue
            sz = torch.stack([Zs[s][int(i0_by_sta[s]) + lo: int(i0_by_sta[s]) + hi]
                              for s in mm], 0).mean(0)
            if clip is not None:
                sz = sz.clamp(max=clip)
            wq = 1.0 if weights is None else float(weights.get(q, 1.0))
            acc += sz * wq
            nz3 += (sz >= W.AC_SITE_Z_COINC).float()
            nsv += 1.0
            SZ[q] = sz
        m, sc = self._robust(acc[None], sub, mad)
        return (acc - m[0]) / sc[0], SZ, nz3, nsv, lo

    def rows(self, i0s, clip=W.AC_SITE_Z_CLIP, sub=W.AC_Z_SUB, mad="lowhalf", weights=None,
             keep_site=True):
        """(net_mad (nb,NV), site_z {site:(nb,NV)}, nz3 (nb,NV), nsv (nb,NV)) for window rows i0s
        (sample indices of the window START).  Defaults reproduce `w_autocorr` exactly."""
        import torch
        dev, n, NV = self.dev, self.n, self.NV
        idx = torch.as_tensor(np.asarray(i0s, dtype=np.int64), device=dev)
        nb = len(idx)

        def _rob(x):
            v = x[:, ::sub] if sub > 1 else x
            m = v.median(dim=-1, keepdim=True).values
            if mad == "full":
                sc = 1.4826 * (v - m).abs().median(dim=-1, keepdim=True).values
            else:
                q = torch.kthvalue(v, max(1, v.shape[-1] // 4), dim=-1, keepdim=True).values
                sc = 1.4826 * (m - q)
            return m, sc.clamp_min(1e-9)

        acc = torch.zeros((nb, NV), device=dev)
        nz3 = torch.zeros((nb, NV), device=dev)
        nsv = torch.zeros((nb, NV), device=dev)
        SZ = {}
        for q in self.have:
            sbuf = torch.zeros((nb, NV), device=dev)
            cnt = torch.zeros((nb, NV), device=dev)
            for s in self.mem[q]:
                if s not in self.D:
                    continue
                a = torch.stack([self.HT[s][c][idx] for c in range(3)], 0)       # (3, nb, n)
                r = a.pow(2).sum(-1, keepdim=True).sqrt()
                a = torch.where(r > 0, a / r.clamp_min(1e-30), torch.zeros_like(a))
                cc = torch.zeros((nb, NV), device=dev)
                for c in range(3):
                    cc += (a[c] @ self.HT[s][c].T) / self.E[s][c][None]
                cc /= 3.0
                m, sc = _rob(cc.index_select(1, self.vidx[s]))
                z = (cc - m) / sc
                vm = self.VAL[s].float()[None]
                sbuf += z * vm
                cnt += vm.expand(nb, NV)
                del cc, z, a, r
            okc = cnt > 0
            sz = torch.where(okc, sbuf / cnt.clamp_min(1.0), torch.zeros((), device=dev))
            if clip is not None:
                sz = sz.clamp(max=clip)
            wq = 1.0 if weights is None else float(weights.get(q, 1.0))
            acc += torch.where(okc, sz * wq, torch.zeros((), device=dev))
            nz3 += ((sz >= W.AC_SITE_Z_COINC) & okc).float()
            nsv += okc.float()
            if keep_site:
                SZ[q] = sz
            else:
                del sz
            del sbuf, cnt, okc
        anyv = (nsv >= W.MIN_SITES)
        vcols = torch.nonzero(anyv.any(0), as_tuple=False).squeeze(1)
        if vcols.numel() < 1000:
            vcols = torch.arange(NV, device=dev)
        m, sc = _rob(acc.index_select(1, vcols))
        nm = (acc - m) / sc
        return nm, SZ, nz3, nsv


def bed_ctx(year=None):
    mem, of, order, det = W.sites(year)
    ST = W.station_table()
    return mem, of, order, det, ST


def kept_bostock(year, window_days=None):
    import p_common as P
    hid, _, _ = S.holdouts()
    B = P.bostock()
    B = B[(B.year == year) & (~B.bkey.isin(hid))].copy()
    B["date"] = B.time.dt.strftime("%Y-%m-%d")
    B["hh"] = B.time.dt.hour.astype(int)
    if window_days is not None:
        B = B[B.date.isin(set(window_days))]
    return B.reset_index(drop=True)


# =================================================================== D4 =========================
def comp_sizes(P, merge_s=W.NODE_MERGE_S):
    """Connected-component sizes of the pair graph, same node/edge rule as `w_cluster` (nodes =
    window times merged at `merge_s`, edges = pairs, components computed day by day) but with a
    sparse-graph solver instead of the Python union-find, so it can be called thousands of times.
    Verified against `w_cluster.components_of` in the queue's `selftest` stage."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    out = []
    if not len(P):
        return out
    for d, g in P.groupby("date"):
        ti = g.t_i.to_numpy(float)
        tj = g.t_j.to_numpy(float)
        t = np.unique(np.concatenate([ti, tj]))
        grp = np.zeros(len(t), int)
        if len(t) > 1:
            grp[1:] = np.cumsum(np.diff(t) > merge_s)
        ng = grp[-1] + 1
        ia = grp[np.searchsorted(t, ti)]
        ib = grp[np.searchsorted(t, tj)]
        A = coo_matrix((np.ones(len(ia)), (ia, ib)), shape=(ng, ng))
        _, lab = connected_components(A, directed=False)
        out += list(np.bincount(lab))
    return out


def fast_threshold(P, n_hours, grid=None):
    """`w_cluster.threshold`'s FA-only rule, evaluated with `comp_sizes`."""
    grid = np.arange(W.T_ENG_FLOOR, W.T_MAX + 1e-9, W.T_GRID) if grid is None else grid
    for T in grid:
        Q = P[P.net_mad >= T]
        fa = len(Q) / max(n_hours, 1)
        sz = comp_sizes(Q) if len(Q) else []
        nc = int(sum(1 for x in sz if x >= W.QUIET_MAX_COMPONENT))
        if fa <= W.QUIET_FA_PAIRS_PER_HOUR and nc == 0:
            return float(T), ("engineering_floor" if T <= W.T_ENG_FLOOR + 1e-9 else "fa_rule")
    return float(W.T_MAX), "saturated"


def per_hour_table(P, tag):
    """pairs / largest component / lag summary per hour, at every T on the grid."""
    rows = []
    if not len(P):
        return pd.DataFrame(columns=["tag", "date", "hour", "T", "n_pairs", "max_component",
                                     "n_comp_ge5", "med_abs_lag_s", "frac_lag_lt_20s",
                                     "n_hour_sites"])
    P = P.copy()
    P["alag"] = (P.t_j - P.t_i).abs()
    for (d, h), g in P.groupby(["date", "hour"]):
        for T in T_GRID:
            Q = g[g.net_mad >= T]
            if not len(Q):
                rows.append(dict(tag=tag, date=d, hour=int(h), T=float(T), n_pairs=0,
                                 max_component=0, n_comp_ge5=0, med_abs_lag_s=np.nan,
                                 frac_lag_lt_20s=np.nan,
                                 n_hour_sites=int(g.n_hour_sites.max())))
                continue
            sz = comp_sizes(Q)
            rows.append(dict(tag=tag, date=d, hour=int(h), T=float(T), n_pairs=int(len(Q)),
                             max_component=int(max(sz)) if sz else 0,
                             n_comp_ge5=int(sum(1 for x in sz if x >= 5)),
                             med_abs_lag_s=float(Q.alag.median()),
                             frac_lag_lt_20s=float((Q.alag < 20.0).mean()),
                             n_hour_sites=int(g.n_hour_sites.max())))
    return pd.DataFrame(rows)


def d4(out, o_run):
    mem, of, order, det, ST = bed_ctx()
    Q = pd.read_parquet(f"{o_run}/ac_quiet.parquet")
    A = pd.read_parquet(f"{o_run}/ac.parquet")
    nq = json.load(open(f"{o_run}/ac_quiet_meta.json"))["n_hours"]
    nt = json.load(open(f"{o_run}/ac_meta.json"))["n_hours"]

    TQ = per_hour_table(Q, "quiet")
    TA = per_hour_table(A, "tremor")
    pd.concat([TQ, TA], ignore_index=True).to_csv(out + "_D4_per_hour.csv", index=False)

    # pooled ratio vs T
    rows = []
    for T in T_GRID:
        q = int((Q.net_mad >= T).sum())
        a = int((A.net_mad >= T).sum())
        qh = TQ[(TQ["T"].round(4) == round(float(T), 4))]
        rows.append(dict(T=float(T), quiet_pairs=q, quiet_per_hour=q / max(nq, 1),
                         tremor_pairs=a, tremor_per_hour=a / max(nt, 1),
                         ratio=(a / max(nt, 1)) / max(q / max(nq, 1), 1e-12),
                         quiet_hours_with_pairs=int((qh.n_pairs > 0).sum()) if len(qh) else 0,
                         quiet_max_component=int(qh.max_component.max()) if len(qh) else 0))
    R = pd.DataFrame(rows)
    R.to_csv(out + "_D4_ratio_vs_T.csv", index=False)

    # leave-one-quiet-hour-out T (DIAGNOSTIC -- the frozen T is not changed)
    loo = []
    hours = sorted({(d, int(h)) for d, h in zip(Q.date, Q.hour)})
    allh = pd.read_csv(f"{W.RES}/step1c_quiet_hours.csv") if os.path.exists(
        f"{W.RES}/step1c_quiet_hours.csv") else None
    for (d, h) in hours:
        Q2 = Q[~((Q.date == d) & (Q.hour == h))]
        T2, b2 = fast_threshold(Q2, nq - 1)
        loo.append(dict(dropped_date=d, dropped_hour=h,
                        n_pairs_dropped=int(((Q.date == d) & (Q.hour == h)).sum()),
                        T_without=float(T2), binds=b2))
    L = pd.DataFrame(loo).sort_values("T_without")
    L.to_csv(out + "_D4_leave_one_hour_out_T.csv", index=False)

    # GPU site attribution on the hours that matter
    worst_q = (Q[Q.net_mad >= 8.0].groupby(["date", "hour"]).size().sort_values(ascending=False)
               .head(3).index.tolist())
    busy_q = (Q.groupby(["date", "hour"]).size().sort_values(ascending=False)
              .head(3).index.tolist())
    busy_t = (A.groupby(["date", "hour"]).size().sort_values(ascending=False)
              .head(4).index.tolist())
    sel = []
    for lab, lst in (("quiet_high_netmad", worst_q), ("quiet_busy", busy_q),
                     ("tremor_busy", busy_t)):
        for (d, h) in lst:
            if (d, int(h)) not in [(x[1], x[2]) for x in sel]:
                sel.append((lab, d, int(h)))
    loader = W.DayLoader(maxn=len(det) + 4)
    att = []
    for lab, d, h in sel:
        src = Q if lab.startswith("quiet") else A
        g = src[(src.date == d) & (src.hour == h)].sort_values("net_mad", ascending=False).head(60)
        if not len(g):
            continue
        H = Hour(loader, d, h, mem, of, order)
        if not H.ok:
            H.close()
            continue
        i0 = np.round((g.t_i.to_numpy(float) - H.t0) * FS).astype(int)
        jj = np.round((g.t_j.to_numpy(float) - H.t0) * FS).astype(int)
        keep = (i0 >= 0) & (i0 + H.n <= H.NH) & (jj >= 0) & (jj < H.NV)
        i0, jj, gk = i0[keep], jj[keep], g[keep]
        if not len(i0):
            H.close()
            continue
        nm, SZ, nz3, nsv = H.rows(i0)
        ar = np.arange(len(i0))
        rec = dict(net_mad=nm[ar, jj].detach().cpu().numpy(),
                   n_site_z3=nz3[ar, jj].detach().cpu().numpy(),
                   n_sites=nsv[ar, jj].detach().cpu().numpy())
        for q in H.have:
            rec[f"z_{q}"] = SZ[q][ar, jj].detach().cpu().numpy()
        Rt = pd.DataFrame(rec)
        Rt["label"] = lab
        Rt["date"] = d
        Rt["hour"] = h
        Rt["stored_net_mad"] = gk.net_mad.to_numpy()
        Rt["abs_lag_s"] = np.abs(gk.t_j.to_numpy() - gk.t_i.to_numpy())
        att.append(Rt)
        del nm, SZ, nz3, nsv
        H.close()
    AT = pd.concat(att, ignore_index=True) if att else pd.DataFrame()
    AT.to_csv(out + "_D4_site_attribution.csv", index=False)

    summ = dict(n_quiet_hours=int(nq), n_tremor_hours=int(nt),
                T_frozen=float(json.load(open(f"{o_run}/T.json"))["T"]),
                quiet_pairs_total=int(len(Q)), tremor_pairs_total=int(len(A)),
                ratio_crosses_one_at_T=float(R[R.ratio < 1.0]["T"].min())
                if (R.ratio < 1.0).any() else None,
                worst_quiet_hour=[f"{d} {h:02d}" for d, h in worst_q],
                T_loo_min=float(L.T_without.min()), T_loo_max=float(L.T_without.max()),
                T_without_worst_hour=float(L.iloc[0].T_without),
                worst_hour_dropped=f"{L.iloc[0].dropped_date} {int(L.iloc[0].dropped_hour):02d}")
    json.dump(summ, open(out + "_D4.json", "w"), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))


# =================================================================== D4b (NaN voiding) =========
def dnan(out, o_run):
    """D4b -- how many hours the engine VOIDED.

    `w_autocorr` normalises each station's cc by the per-window energy
    `E = ew.clamp_min(1e-10 * ew.median()).sqrt()`.  If a station-hour is more than half flat
    (a gap filled with zeros, a dead channel), `ew.median()` is 0, the clamp is a no-op, E is 0 and
    the cc is inf -> the station z is NaN -> the site mean is NaN -> the equal-weight site sum is
    NaN at EVERY lag -> `hit = (nm >= t_record) & ...` is False everywhere and the WHOLE HOUR
    silently records zero pairs.  This measures how often that happened on the bed.
    """
    mem, of, order, det, ST = bed_ctx()
    H1 = pd.read_csv(f"{W.RES}/step1c_candidate_hours.csv") if W.BED_MODE == "c" \
        else W.candidate_hours()
    H2, _ = W.quiet_hours()
    hours = ([("tremor", r.date, int(r.hour)) for r in H1.itertuples()]
             + [("quiet", r.date, int(r.hour)) for r in H2.itertuples()])
    got = {}
    for tag, sub in (("tremor", "ac"), ("quiet", "ac_quiet")):
        for p in glob.glob(f"{o_run}/{sub}.parts/*.parquet"):
            b = os.path.basename(p)[:-8]
            try:
                got[(tag, b[:10], int(b[11:]))] = int(len(pd.read_parquet(p)))
            except Exception:
                pass
    loader = W.DayLoader(maxn=len(det) + 4)
    n = int(round(W.AC_WIN_S * FS))
    rows = []
    for tag, d, h in hours:
        bad_sta, seen = [], []
        for s in det:
            v = loader.hour(s, d, h)
            if v is None:
                continue
            x = v[1].astype(np.float64)
            seen.append(s)
            NH = x.shape[1]
            if NH < n + 1:
                bad_sta.append(s)
                continue
            # exactly `w_autocorr`'s per-component window energy and its GLOBAL median, which is
            # what the clamp uses: E = ew.clamp_min(1e-10 * ew.median()).sqrt()
            cs = np.concatenate([np.zeros((3, 1)), np.cumsum(x * x, 1)], 1)
            ew = cs[:, n:NH + 1] - cs[:, :NH - n + 1]
            if np.median(ew) <= 0:
                bad_sta.append(s)
        # A single degenerate MEMBER poisons its site: the member z is NaN at every lag, and
        # `sbuf += z * vm` makes NaN * 0 = NaN, so the site mean -- and then the whole network sum
        # -- is NaN everywhere.  So a site is voided if ANY present member is degenerate.
        bad_sites = sorted({of[s] for s in bad_sta})
        rows.append(dict(tag=tag, date=d, hour=h, n_sta_seen=len(seen),
                         n_sta_degenerate=len(bad_sta), degenerate_sta=",".join(bad_sta),
                         n_sites_voided=len(bad_sites), voided_sites=",".join(bad_sites),
                         hour_voided=bool(bad_sites),
                         pairs_recorded=got.get((tag, d, h))))
    R = pd.DataFrame(rows)
    R.to_csv(out + "_D4b_voided_hours.csv", index=False)
    summ = {}
    for tag, g in R.groupby("tag"):
        z = g[g.pairs_recorded == 0] if "pairs_recorded" in g else g.iloc[:0]
        summ[tag] = dict(
            n_hours=int(len(g)),
            n_hours_voided=int(g.hour_voided.sum()),
            frac_hours_voided=float(g.hour_voided.mean()),
            n_hours_zero_pairs=int((g.pairs_recorded == 0).sum()),
            n_voided_and_zero=int((g.hour_voided & (g.pairs_recorded == 0)).sum()),
            n_zero_not_voided=int(((~g.hour_voided) & (g.pairs_recorded == 0)).sum()),
            n_voided_not_zero=int((g.hour_voided & (g.pairs_recorded > 0)).sum()),
            top_offenders=(g[g.hour_voided].voided_sites.value_counts().head(5).to_dict()
                           if g.hour_voided.any() else {}))
    json.dump(summ, open(out + "_D4b.json", "w"), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))


# =================================================================== D5 =========================
def _true_pairs(B, min_events, dt_lo=15.0, dt_hi=1200.0, max_per_family=40):
    """[(bkey, date, hour, ta, tb)] same-family pairs inside one hour, kept families only."""
    g = B.groupby("bkey").size()
    fams = sorted(int(b) for b in g[g >= min_events].index)
    out = []
    for bk in fams:
        gg = B[B.bkey == bk].sort_values("t_origin")
        t = gg.t_origin.to_numpy(float)
        d = gg.date.to_numpy()
        h = gg.hh.to_numpy()
        got = 0
        for i in range(len(t) - 1):
            for j in range(i + 1, min(i + 12, len(t))):
                if d[i] != d[j] or h[i] != h[j]:
                    continue
                dt = t[j] - t[i]
                if dt_lo <= dt <= dt_hi:
                    out.append((bk, str(d[i]), int(h[i]), float(t[i]), float(t[j])))
                    got += 1
                    break
            if got >= max_per_family:
                break
    return out, fams


def _pair_rows(H, tt, ta, tb, scan_pre=10.0, scan_post=2.0, **kw):
    """Measure one true pair the way the engine would: scan the engine's own 0.5-s window grid over
    the alignments that can contain the family's S arrivals, take the best net_mad at the true lag.
    Returns (best dict) or None."""
    tmin = ta + min(tt.values())
    lo = tmin - scan_pre
    hi = tmin + scan_post
    k0 = int(np.ceil((lo - H.t0) * FS / (W.AC_HOP_S * FS)))
    k1 = int(np.floor((hi - H.t0) * FS / (W.AC_HOP_S * FS)))
    hop = int(round(W.AC_HOP_S * FS))
    i0s = np.array([k * hop for k in range(max(k0, 0), k1 + 1)], dtype=int)
    i0s = i0s[(i0s >= 0) & (i0s + H.n <= H.NH)]
    if not len(i0s):
        return None
    dj = int(round((tb - ta) * FS))
    jj = i0s + dj
    ok = (jj >= 0) & (jj < H.NV)
    i0s, jj = i0s[ok], jj[ok]
    if not len(i0s):
        return None
    nm, SZ, nz3, nsv = H.rows(i0s, **kw)
    ar = np.arange(len(i0s))
    v = nm[ar, jj].detach().cpu().numpy()
    # A station whose hour is (nearly) all gap gives E = 0 -> inf cc -> NaN net_mad on that row.
    # `w_autocorr` drops such rows implicitly (NaN >= T is False), so the diagnostic does the same
    # rather than letting np.argmax return the first NaN.
    nan_frac = float(np.mean(~np.isfinite(v)))
    v = np.where(np.isfinite(v), v, -np.inf)
    if not np.isfinite(v).any():
        del nm, SZ, nz3, nsv
        return None
    b = int(np.argmax(v))
    nmm = nm.detach().cpu().numpy()
    wj = int(round(0.5 * FS))                       # +-0.5 s, for catalogue-time jitter
    row = nmm[b]
    jt = int(jj[b])
    _rng = np.random.default_rng(20260920 + jt)
    row = np.where(np.isfinite(row), row, -np.inf)
    nullw = [float(row[max(k - wj, 0):k + wj + 1].max())
             for k in _rng.integers(wj + 1, max(len(row) - wj - 1, wj + 2), 20)]
    rec = dict(net_mad_win=float(row[max(jt - wj, 0):jt + wj + 1].max()),
               net_mad_win_null=float(np.median(nullw)),
               net_mad_win_null_max=float(np.max(nullw)),
               net_mad=float(v[b]), n_site_z3=float(nz3[ar, jj].detach().cpu().numpy()[b]),
               n_sites=float(nsv[ar, jj].detach().cpu().numpy()[b]),
               n_rows=int(len(i0s)), align_s=float(H.t0 + i0s[b] / FS - tmin),
               nan_row_frac=nan_frac,
               net_mad_max_any_lag=float(np.nanmax(nmm[b])),
               net_mad_self_lag=float(nmm[b, i0s[b]]) if i0s[b] < H.NV else np.nan)
    for q in H.have:
        rec[f"z_{q}"] = float(SZ[q][ar, jj].detach().cpu().numpy()[b])
    del nm, SZ, nz3, nsv
    return rec


def _pair_moveout(H, tt, ta, tb, **kw):
    """Same true pair, ORACLE moveout-aware cut (diagnostic)."""
    i0 = {s: int(round((ta + tt[s] - 1.0 - H.t0) * FS)) for s in tt}
    i0 = {s: v for s, v in i0.items() if 0 <= v and v + H.n <= H.NH}
    if len({H.of[s] for s in i0}) < W.MIN_SITES:
        return None
    r = H.moveout_rows(i0, **kw)
    if r is None:
        return None
    nm, SZ, nz3, nsv, lo = r
    d = int(round((tb - ta) * FS)) - lo
    if d < 0 or d >= nm.shape[0]:
        return None
    _nm = nm.detach().cpu().numpy()
    _nm = np.where(np.isfinite(_nm), _nm, -np.inf)
    if not np.isfinite(_nm).any():
        del nm, SZ, nz3, nsv
        return None
    wj = int(round(0.5 * FS))
    _rng = np.random.default_rng(20260921 + d)
    nullw = [float(_nm[max(k - wj, 0):k + wj + 1].max())
             for k in _rng.integers(wj + 1, max(len(_nm) - wj - 1, wj + 2), 20)]
    out = dict(mo_net_mad_win=float(_nm[max(d - wj, 0):d + wj + 1].max()),
               mo_net_mad_win_null=float(np.median(nullw)),
               mo_net_mad_win_null_max=float(np.max(nullw)),
               mo_net_mad=float(_nm[d]), mo_n_site_z3=float(nz3[d].item()),
               mo_n_sites=float(nsv[d].item()),
               mo_net_mad_max_any_lag=float(_nm.max()),
               mo_rank=float((_nm >= _nm[d]).sum()))
    for q, z in SZ.items():
        out[f"moz_{q}"] = float(z[d].item())
    del nm, SZ, nz3, nsv
    return out


def d5(out, min_events=50, max_hours=40, max_pairs_per_hour=40, budget=1200):
    mem, of, order, det, ST = bed_ctx()
    import p_common as P
    FAM = P.families()
    FT = S.fam_tt(det, ST=ST)
    days = {d.strftime("%Y-%m-%d") for d in W.bed_days()}
    B = kept_bostock(W.STEP1_YEAR, days)
    pairs, fams = _true_pairs(B, min_events)
    byh = collections.defaultdict(list)
    for (bk, d, h, ta, tb) in pairs:
        byh[(d, h)].append((bk, ta, tb))
    hours = sorted(byh, key=lambda k: -len(byh[k]))[:max_hours]
    loader = W.DayLoader(maxn=len(det) + 4)
    recs, n = [], 0
    for (d, h) in hours:
        if n >= budget:
            break
        H = Hour(loader, d, h, mem, of, order)
        if not H.ok:
            H.close()
            continue
        pl = byh[(d, h)][:max_pairs_per_hour]
        for (bk, ta, tb) in pl:
            if n >= budget:
                break
            tt = {s: FT[(bk, s)] for s in det if s in H.D}
            if len(tt) < W.MIN_SITES:
                continue
            r = _pair_rows(H, tt, ta, tb)
            if r is None:
                continue
            mo = _pair_moveout(H, tt, ta, tb)
            if mo:
                r.update(mo)
            r.update(bkey=int(bk), date=d, hour=int(h), dt=float(tb - ta),
                     n_hour_sites=len(H.have))
            recs.append(r)
            n += 1
        H.close()
        W.log(f"[D5] {d} {h:02d}h -> {n} pairs measured")
    D = pd.DataFrame(recs)
    D.to_csv(out + "_D5_true_pairs.csv", index=False)

    # per-site long form vs distance
    rows = []
    for r in recs:
        for q in order:
            k = f"z_{q}"
            if k not in r or not np.isfinite(r[k]):
                continue
            mems = [s for s in mem[q] if s in ST.index]
            dd = float(np.median([S.epi_km(FAM.lat[r["bkey"]], FAM.lon[r["bkey"]],
                                           ST.lat[s], ST.lon[s]) for s in mems])) if mems else np.nan
            rows.append(dict(bkey=r["bkey"], site=q, dist_km=dd, z=float(r[k]),
                             z_moveout=float(r.get(f"moz_{q}", np.nan)),
                             net_mad=r["net_mad"], mo_net_mad=r.get("mo_net_mad", np.nan)))
    SZ = pd.DataFrame(rows)
    SZ.to_csv(out + "_D5_site_z_vs_distance.csv", index=False)
    bysite = (SZ.groupby("site").agg(dist_km=("dist_km", "median"), n=("z", "size"),
                                     z_med=("z", "median"), z_p90=("z", lambda x: x.quantile(0.9)),
                                     frac_z_ge3=("z", lambda x: float((x >= 3).mean())),
                                     zmo_med=("z_moveout", "median"),
                                     zmo_p90=("z_moveout", lambda x: x.quantile(0.9)),
                                     frac_zmo_ge3=("z_moveout", lambda x: float((x >= 3).mean())))
              .sort_values("dist_km"))
    bysite.to_csv(out + "_D5_by_site.csv")
    summ = dict(n_families_ge_min=len(fams), min_events=int(min_events),
                n_pairs_measured=int(len(D)), n_hours=len(hours),
                median_net_mad=float(D.net_mad.median()) if len(D) else None,
                p90_net_mad=float(D.net_mad.quantile(0.9)) if len(D) else None,
                max_net_mad=float(D.net_mad.max()) if len(D) else None,
                frac_ge_12p5=float((D.net_mad >= 12.5).mean()) if len(D) else None,
                frac_ge_7p75=float((D.net_mad >= 7.75).mean()) if len(D) else None,
                frac_ge_5p0=float((D.net_mad >= 5.0).mean()) if len(D) else None,
                frac_2sites_z3=float((D.n_site_z3 >= 2).mean()) if len(D) else None,
                median_net_mad_self_lag=float(D.net_mad_self_lag.median())
                if len(D) and "net_mad_self_lag" in D else None,
                moveout_median_net_mad=float(D.mo_net_mad.median())
                if len(D) and "mo_net_mad" in D else None,
                moveout_p90_net_mad=float(D.mo_net_mad.quantile(0.9))
                if len(D) and "mo_net_mad" in D else None,
                moveout_max_net_mad=float(D.mo_net_mad.max())
                if len(D) and "mo_net_mad" in D else None,
                moveout_frac_ge_12p5=float((D.mo_net_mad >= 12.5).mean())
                if len(D) and "mo_net_mad" in D else None,
                moveout_frac_ge_7p75=float((D.mo_net_mad >= 7.75).mean())
                if len(D) and "mo_net_mad" in D else None,
                moveout_frac_ge_5p0=float((D.mo_net_mad >= 5.0).mean())
                if len(D) and "mo_net_mad" in D else None,
                moveout_frac_2sites_z3=float((D.mo_n_site_z3 >= 2).mean())
                if len(D) and "mo_n_site_z3" in D else None,
                median_net_mad_win=float(D.net_mad_win.median())
                if len(D) and "net_mad_win" in D else None,
                median_net_mad_win_null=float(D.net_mad_win_null.median())
                if len(D) and "net_mad_win_null" in D else None,
                frac_win_gt_null_max=float((D.net_mad_win > D.net_mad_win_null_max).mean())
                if len(D) and "net_mad_win" in D else None,
                moveout_median_net_mad_win=float(D.mo_net_mad_win.median())
                if len(D) and "mo_net_mad_win" in D else None,
                moveout_median_net_mad_win_null=float(D.mo_net_mad_win_null.median())
                if len(D) and "mo_net_mad_win_null" in D else None,
                moveout_frac_win_gt_null_max=float((D.mo_net_mad_win
                                                    > D.mo_net_mad_win_null_max).mean())
                if len(D) and "mo_net_mad_win" in D else None,
                moveout_frac_win_ge_12p5=float((D.mo_net_mad_win >= 12.5).mean())
                if len(D) and "mo_net_mad_win" in D else None,
                moveout_frac_win_ge_7p75=float((D.mo_net_mad_win >= 7.75).mean())
                if len(D) and "mo_net_mad_win" in D else None,
                moveout_frac_win_ge_5p0=float((D.mo_net_mad_win >= 5.0).mean())
                if len(D) and "mo_net_mad_win" in D else None,
                frac_win_ge_12p5=float((D.net_mad_win >= 12.5).mean())
                if len(D) and "net_mad_win" in D else None,
                frac_win_ge_7p75=float((D.net_mad_win >= 7.75).mean())
                if len(D) and "net_mad_win" in D else None,
                frac_win_ge_5p0=float((D.net_mad_win >= 5.0).mean())
                if len(D) and "net_mad_win" in D else None,
                frac_hit_at_frozen_T=float(((D.net_mad >= 12.5) & (D.n_site_z3 >= 2)
                                            & (D.n_sites >= W.MIN_SITES)).mean()) if len(D) else None,
                by_site=json.loads(bysite.reset_index().to_json(orient="records")))
    json.dump(summ, open(out + "_D5.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in summ.items() if k != "by_site"}, indent=1, default=str))
    print(bysite.to_string())


# =================================================================== D6 =========================
VARIANTS = [
    ("frozen",        dict(win_s=10.0, kw={})),
    ("win6",          dict(win_s=6.0, kw={})),
    ("win8",          dict(win_s=8.0, kw={})),
    ("win15",         dict(win_s=15.0, kw={})),
    ("noclip",        dict(win_s=10.0, kw=dict(clip=None))),
    ("mad_full",      dict(win_s=10.0, kw=dict(mad="full"))),
    ("dist_weighted", dict(win_s=10.0, kw=dict(weights="dist"))),
]


def d6(out, n_hours=8, pairs_per_hour=15, quiet_hours=4):
    mem, of, order, det, ST = bed_ctx()
    import p_common as P
    FAM = P.families()
    FT = S.fam_tt(det, ST=ST)
    days = {d.strftime("%Y-%m-%d") for d in W.bed_days()}
    B = kept_bostock(W.STEP1_YEAR, days)
    pairs, _ = _true_pairs(B, 50)
    byh = collections.defaultdict(list)
    for (bk, d, h, ta, tb) in pairs:
        byh[(d, h)].append((bk, ta, tb))
    hrs = sorted(byh, key=lambda k: -len(byh[k]))[:n_hours]
    loader = W.DayLoader(maxn=len(det) + 4)

    def wmap(bk):
        w = {}
        for q in order:
            mm = [s for s in mem[q] if s in ST.index]
            d = np.median([S.epi_km(FAM.lat[bk], FAM.lon[bk], ST.lat[s], ST.lon[s])
                           for s in mm]) if mm else np.nan
            w[q] = 1.0 / max(float(d), 1.0)
        tot = sum(w.values())
        return {q: v * len(order) / tot for q, v in w.items()}

    recs = []
    for (d, h) in hrs:
        for name, cfg in VARIANTS:
            H = Hour(loader, d, h, mem, of, order, win_s=cfg["win_s"])
            if not H.ok:
                H.close()
                continue
            for (bk, ta, tb) in byh[(d, h)][:pairs_per_hour]:
                tt = {s: FT[(bk, s)] for s in det if s in H.D}
                if len(tt) < W.MIN_SITES:
                    continue
                kw = dict(cfg["kw"])
                if kw.get("weights") == "dist":
                    kw["weights"] = wmap(bk)
                r = _pair_rows(H, tt, ta, tb, **kw)
                if r is None:
                    continue
                recs.append(dict(variant=name, bkey=int(bk), date=d, hour=int(h),
                                 net_mad=r["net_mad"], n_site_z3=r["n_site_z3"],
                                 n_sites=r["n_sites"]))
            H.close()
        W.log(f"[D6/true] {d} {h:02d}h done ({len(recs)} rows)")
    TP = pd.DataFrame(recs)
    TP.to_csv(out + "_D6_true_pair_variants.csv", index=False)
    TS = (TP.groupby("variant").agg(n=("net_mad", "size"), med=("net_mad", "median"),
                                    p90=("net_mad", lambda x: x.quantile(0.9)),
                                    mx=("net_mad", "max"),
                                    frac_ge5=("net_mad", lambda x: float((x >= 5).mean())),
                                    frac_ge7p75=("net_mad", lambda x: float((x >= 7.75).mean())),
                                    frac_ge12p5=("net_mad", lambda x: float((x >= 12.5).mean())))
          if len(TP) else pd.DataFrame())
    TS.to_csv(out + "_D6_true_pair_summary.csv")

    # ---- quiet-hour false-alarm side of the same ablations
    QH, _ = W.quiet_hours()
    QH = QH.head(quiet_hours)
    qrec = []
    for r in QH.itertuples():
        for name, cfg in VARIANTS:
            H = Hour(loader, r.date, int(r.hour), mem, of, order, win_s=cfg["win_s"])
            if not H.ok:
                H.close()
                continue
            kw = dict(cfg["kw"])
            if kw.get("weights") == "dist":
                kw["weights"] = wmap(int(FAM.index[0]))
            hop = int(round(W.AC_HOP_S * FS))
            NW = (H.NH - H.n) // hop + 1
            minlag = int(round(W.AC_MIN_LAG_S * FS))
            best = []
            blk = 256
            import torch
            for b0 in range(0, NW, blk):
                i0s = np.arange(b0, min(b0 + blk, NW)) * hop
                i0s = i0s[i0s + H.n <= H.NH]
                if not len(i0s):
                    continue
                nm, SZ, nz3, nsv = H.rows(i0s, keep_site=False, **kw)
                jj = torch.arange(H.NV, device=H.dev)[None, :]
                ti = torch.as_tensor(i0s, device=H.dev)[:, None]
                ok = ((jj - ti).abs() >= minlag) & (nz3 >= W.AC_MIN_SITE_Z_SITES) \
                     & (nsv >= W.MIN_SITES)
                v = torch.where(ok, nm, torch.full((), -1e9, device=H.dev))
                best.append(v.max(dim=1).values.detach().cpu().numpy())
                del nm, SZ, nz3, nsv, ok, v
                torch.cuda.empty_cache()
            bb = np.concatenate(best) if best else np.array([])
            for T in (5.0, 6.0, 7.0, 7.75, 8.0, 10.0, 12.5):
                qrec.append(dict(variant=name, date=r.date, hour=int(r.hour), T=float(T),
                                 n_windows_over_T=int((bb >= T).sum()), n_windows=int(len(bb)),
                                 max_net_mad=float(bb.max()) if len(bb) else np.nan))
            H.close()
        W.log(f"[D6/quiet] {r.date} {int(r.hour):02d}h done")
    QR = pd.DataFrame(qrec)
    QR.to_csv(out + "_D6_quiet_variants.csv", index=False)
    summ = dict(true_pair_summary=json.loads(TS.reset_index().to_json(orient="records"))
                if len(TS) else [],
                quiet_summary=json.loads(QR.groupby(["variant", "T"]).n_windows_over_T.mean()
                                         .reset_index().to_json(orient="records"))
                if len(QR) else [])
    json.dump(summ, open(out + "_D6.json", "w"), indent=1, default=str)
    print(TS.to_string() if len(TS) else "no true-pair rows")


# =================================================================== D6 native rate =============
def d6_native(out, sta=None, n_pairs=120):
    """Is the frozen 20 Hz / 1-8 Hz cache itself the limitation?  Recompute ONE station-hour from
    the raw 100 Hz mseed already on disk (no download) at 2-8 Hz and at 1-8 Hz, and compare the
    per-station true-pair z with the frozen cache's."""
    import obspy
    import torch
    from scipy.signal import butter, sosfiltfilt
    mem, of, order, det, ST = bed_ctx()
    FT = S.fam_tt(det, ST=ST)
    days = {d.strftime("%Y-%m-%d") for d in W.bed_days()}
    B = kept_bostock(W.STEP1_YEAR, days)
    pairs, _ = _true_pairs(B, 50)
    byh = collections.defaultdict(list)
    for (bk, d, h, ta, tb) in pairs:
        byh[(d, h)].append((bk, ta, tb))
    RAW = f"{W.STUDY}/traces/bed2005_polaris"
    cands = [sta] if sta else ["TWKB", "GLBC", "SHVB", "TSJB", "MGCB"]
    out_rows = []
    chosen = None
    for (d, h) in sorted(byh, key=lambda k: -len(byh[k])):
        ts = pd.Timestamp(d)
        for s in cands:
            if s not in det:
                continue
            hit = glob.glob(f"{RAW}/*.{s}/{ts.year}/{ts.dayofyear:03d}.mseed")
            if hit:
                chosen = (s, d, int(h), hit[0])
                break
        if chosen:
            break
    if chosen is None:
        json.dump(dict(status="no raw station-day on disk for any true-pair hour"),
                  open(out + "_D6native.json", "w"), indent=1)
        print("no raw station-day available")
        return
    s, d, h, path = chosen
    W.log(f"[D6native] {s} {d} {h:02d}h  <- {path}")
    st = obspy.read(path)
    st.merge(method=1, fill_value=0)
    t0 = obspy.UTCDateTime(pd.Timestamp(d).timestamp() + h * 3600)
    st.trim(t0, t0 + 3600, pad=True, fill_value=0)
    chan = {tr.stats.channel[-1]: tr for tr in st}
    if not all(c in chan for c in "ZNE"):
        json.dump(dict(status=f"{s} not 3-C in the raw file",
                       channels=sorted(chan)), open(out + "_D6native.json", "w"), indent=1)
        print("not 3-C")
        return
    fs_n = float(chan["Z"].stats.sampling_rate)
    x = np.stack([chan[c].data.astype(np.float64) for c in "ZNE"], 0)
    x = x - x.mean(1, keepdims=True)

    def band(lo, hi, fsx, arr):
        sos = butter(4, [lo / (fsx / 2), hi / (fsx / 2)], btype="band", output="sos")
        return np.ascontiguousarray(sosfiltfilt(sos, arr, axis=-1).astype(np.float32))

    loader = W.DayLoader(maxn=4)
    v = loader.hour(s, d, h)
    if v is None:
        json.dump(dict(status="frozen cache hour missing"), open(out + "_D6native.json", "w"),
                  indent=1)
        return
    tcache, xc = v

    def measure(arr, fsx, tag):
        n = int(round(W.AC_WIN_S * fsx))
        NH = arr.shape[1]
        NV = NH - n + 1
        X = torch.from_numpy(np.ascontiguousarray(arr)).to("cuda")
        HT = [X[c].unfold(0, n, 1)[:NV].contiguous() for c in range(3)]
        xd = X.double()
        cs = torch.cumsum(torch.cat([torch.zeros((3, 1), device="cuda", dtype=torch.float64),
                                     xd * xd], 1), 1)
        ew = (cs[:, n:NH + 1] - cs[:, :NV])[:, :NV]
        E = ew.clamp_min(1e-10 * ew.median()).sqrt().float()
        rows = []
        rng = np.random.default_rng(20260920)
        base = float(t0.timestamp)   # F7: obspy UTCDateTime.timestamp is a PROPERTY, not a call
        for (bk, ta, tb) in byh[(d, h)][:n_pairs]:
            tt = FT[(bk, s)]
            ia = int(round((ta + tt - 1.0 - base) * fsx))
            jb = int(round((tb + tt - 1.0 - base) * fsx))
            if ia < 0 or ia + n > NH or jb < 0 or jb >= NV:
                continue
            a = X[:, ia:ia + n].float()
            r = a.pow(2).sum(-1, keepdim=True).sqrt()
            a = torch.where(r > 0, a / r.clamp_min(1e-30), torch.zeros_like(a))
            cc = torch.zeros(NV, device="cuda")
            for c in range(3):
                cc += (a[c][None] @ HT[c].T)[0] / E[c]
            cc /= 3.0
            m = cc.median()
            q = torch.kthvalue(cc, max(1, NV // 4)).values
            sc = (1.4826 * (m - q)).clamp_min(1e-9)
            z = ((cc - m) / sc).detach().cpu().numpy()
            w = int(round(1.0 * fsx))
            zt = float(z[max(jb - w, 0):jb + w + 1].max())
            zn = [float(z[jr - w:jr + w + 1].max())
                  for jr in rng.integers(w + 1, NV - w - 1, 10)]
            rows.append(dict(tag=tag, fs=fsx, bkey=int(bk), z_true=zt,
                             z_null_median=float(np.median(zn)), z_null_max=float(max(zn)),
                             z_self_ceiling=float(1.0 / sc)))
            del cc, z, a, r
        del X, HT, E, xd, cs, ew
        torch.cuda.empty_cache()
        return rows

    out_rows += measure(xc[:, :int(3600 * FS)], FS, "frozen_20Hz_1-8Hz")
    out_rows += measure(band(2.0, 8.0, fs_n, x), fs_n, f"native_{int(fs_n)}Hz_2-8Hz")
    out_rows += measure(band(1.0, 8.0, fs_n, x), fs_n, f"native_{int(fs_n)}Hz_1-8Hz")
    R = pd.DataFrame(out_rows)
    R.to_csv(out + "_D6_native_rate.csv", index=False)
    R["above_null_max"] = (R.z_true > R.z_null_max).astype(float)
    G = R.groupby("tag").agg(n=("z_true", "size"), z_true_med=("z_true", "median"),
                             z_null_med=("z_null_median", "median"),
                             frac_gt_null_max=("above_null_max", "mean"),
                             ceiling=("z_self_ceiling", "median"))
    G.to_csv(out + "_D6_native_summary.csv")
    json.dump(dict(station=s, date=d, hour=h, fs_native=fs_n,
                   summary=json.loads(G.reset_index().to_json(orient="records"))),
              open(out + "_D6native.json", "w"), indent=1, default=str)
    print(G.to_string())


# =================================================================== D7 =========================
def d7(out, o_run):
    mem, of, order, det, ST = bed_ctx()
    days = sorted({d.strftime("%Y-%m-%d") for d in W.bed_days()})
    cfg = W.bed_cfg()
    win = [d for d in days if cfg.get("day0", "") <= d <= cfg.get("day1", "9999")]
    B = kept_bostock(W.STEP1_YEAR)                                  # all of the year, kept only
    import p_common as P
    hid, _, _ = S.holdouts()
    ALL = P.bostock()
    ALL = ALL[ALL.year == W.STEP1_YEAR].copy()
    ALL["date"] = ALL.time.dt.strftime("%Y-%m-%d")
    ALL["hh"] = ALL.time.dt.hour.astype(int)

    H = pd.read_csv(f"{W.RES}/step1c_candidate_hours.csv")
    sel = {(str(d), int(h)) for d, h in zip(H.date, H.hour)}
    ran = set()
    for p in glob.glob(f"{o_run}/ac.parts/*.parquet"):
        b = os.path.basename(p)[:-8]
        ran.add((b[:10], int(b[11:])))

    def frac(df, keys):
        if not len(df):
            return 0.0, 0
        m = [(d, h) in keys for d, h in zip(df.date, df.hh)]
        return float(np.mean(m)), int(np.sum(m))

    rows = []
    for lab, df in (("all_2005", ALL), ("all_2005_kept", ALL[~ALL.bkey.isin(hid)]),
                    ("window", ALL[ALL.date.isin(win)]),
                    ("window_kept", ALL[ALL.date.isin(win) & ~ALL.bkey.isin(hid)])):
        f1, n1 = frac(df, sel)
        f2, n2 = frac(df, ran)
        rows.append(dict(subset=lab, n_events=int(len(df)),
                         frac_in_selected_hours=f1, n_in_selected_hours=n1,
                         frac_in_hours_actually_run=f2, n_in_hours_actually_run=n2))
    R = pd.DataFrame(rows)
    R.to_csv(out + "_D7_hour_coverage.csv", index=False)

    Wd = ALL[ALL.date.isin(win) & ~ALL.bkey.isin(hid)]
    per = Wd.groupby(["date", "hh"]).size().rename("n_events").reset_index()
    per["selected"] = [(d, h) in sel for d, h in zip(per.date, per.hh)]
    per["ran"] = [(d, h) in ran for d, h in zip(per.date, per.hh)]
    per.to_csv(out + "_D7_events_per_hour.csv", index=False)
    summ = dict(n_selected_hours=len(sel), n_hours_run=len(ran),
                n_window_hours_total=int(len(win) * 24),
                events_window_kept=int(len(Wd)),
                frac_window_kept_in_selected=float(R[R.subset == "window_kept"]
                                                   .frac_in_selected_hours.iloc[0]),
                frac_window_kept_in_run=float(R[R.subset == "window_kept"]
                                              .frac_in_hours_actually_run.iloc[0]),
                median_events_in_selected_hour=float(per[per.selected].n_events.median())
                if per.selected.any() else None,
                median_events_in_unselected_hour=float(per[~per.selected].n_events.median())
                if (~per.selected).any() else None,
                starved=bool(R[R.subset == "window_kept"].frac_in_selected_hours.iloc[0] < 0.60))
    json.dump(summ, open(out + "_D7.json", "w"), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))
    print(R.to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["d4", "d4b", "d5", "d6", "d6native", "d7"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--run", default="results/step1c")
    ap.add_argument("--min-events", type=int, default=50)
    ap.add_argument("--budget", type=int, default=1200)
    ap.add_argument("--sta", default="")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    if a.cmd == "d4":
        d4(a.out, a.run)
    elif a.cmd == "d4b":
        dnan(a.out, a.run)
    elif a.cmd == "d5":
        d5(a.out, min_events=a.min_events, budget=a.budget)
    elif a.cmd == "d6":
        d6(a.out)
    elif a.cmd == "d6native":
        d6_native(a.out, sta=a.sta or None)
    else:
        d7(a.out, a.run)


if __name__ == "__main__":
    main()
