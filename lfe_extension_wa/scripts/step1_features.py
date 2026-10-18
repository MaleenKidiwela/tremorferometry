#!/usr/bin/env python
"""STEP 1 (Merlin, 2026-09-30) — per-candidate FEATURES for the per-repeat probability (local FDR).

Candidates = every T1 "ALL templates" matched-filter peak (record tier, >= 8.0 MAD; nothing pre-cut) of
results/mf_mt1_all (series data, control, control2) plus two further independent site-rolled scrambled
controls from results/mf_step1c34 (same run, different seed; its control / control2 -> c3 / c4; its data
series is a duplicate of mt1_all's and is only used as a check).  Series labels: data, c1, c2, c3, c4.

Per candidate (all computable after the array era; no truth-based feature):
  net_cc, net_mad                       from the MF (net_re = recomputed network cc, a check)
  n_sta, n_sites, wsum                  family-day network (stations with a template of the day's node)
  per SITE (within-site w-weighted mean of the station cc at the candidate's aligned position; for a
  scrambled control the site is read at its rolled position, exactly as the MF built the control):
     site_cc_max, site_cc_2nd, site_z_max, site_z_2nd, site_z_min, n_z2, n_z3, n_z5, dom_top
     z = (site cc - median) / (1.4826 MAD) of THAT site's series for this family and day, from 1000
     random valid origin samples of the same day (the site's own daily noise level)
     dom_top = share of the network sum carried by the largest single-site contribution
  HF:  hfr = E(10-19 Hz)/E(2-8 Hz) over the 12-s window at the family-day's HF station (wa_hf rule, all
       templates: clear first, then w), native rate; hfr_rel = log10(hfr) - median log10 of 300 random
       windows of the same HF station-day.  Scrambled control: read at the HF station's rolled position.

  python scripts/step1_features.py --workers 12 [--days 2010-08-20,2010-07-15]
  -> results/step1_prob/feat/feat_<date>.parquet
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, zlib
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

# env overrides (STEP 2 group refit, 2026-09-30): defaults = the STEP 1 run
OD = os.environ.get("STEP1_FEATDIR", f"{W.RES}/step1_prob/feat")
LOG = f"{W.STUDY}/logs/step1.log"
MF1, MF2 = os.environ.get("STEP1_MF1", "mt1_all"), os.environ.get("STEP1_MF2", "step1c34")
TPLTAG = os.environ.get("STEP1_TPL", "build")
NRAND, NRAND_HF = 1000, 300
_G = {}


def _init():
    import wa_mf
    T, Q = wa_mf.load_tpl([TPLTAG], all_templates=True)
    Q = Q.copy()
    Q["clear"] = Q.clear.fillna(False).astype(bool)
    Q = Q.sort_values(["clear", "w"], ascending=[False, False])
    S = W.load_stations()
    _G.update(T=T, Q=Q, SITE=dict(zip(S.netsta, S.site)))


def guard_np(x, live, m, n=W.NT, floor=1e-3):
    """numpy copy of wa_mf.guard: robust energy floor + cache mask, dilated by the template length.
    Returns bad (NV,) and the window energy per channel ew (3, NV) float64."""
    NV = x.shape[1] - n + 1
    xd = x.astype(np.float64)
    EN = np.concatenate([np.zeros((3, 1)), np.cumsum(xd * xd, 1)], 1)
    ew = (EN[:, n:] - EN[:, :-n])[:, :NV]
    del EN, xd
    ewl = ew[live]
    med_ch = np.median(ewl, 1)
    e3 = ewl.sum(0); med3 = np.median(e3)
    ref = np.maximum(med_ch, 1e-2 * med3)
    bad = (ewl < floor * ref[:, None]).any(0) | (e3 < floor * med3) | (ewl.min(0) <= 0)
    if med3 <= 0:
        bad[:] = True
    cmk = np.concatenate([[0], np.cumsum(~m)])
    bad |= (cmk[n:n + NV] - cmk[:NV]) > 0
    if bad.any():
        cb = np.concatenate([[0], np.cumsum(bad)])
        lo = np.clip(np.arange(NV) - n, 0, NV); hi = np.clip(np.arange(NV) + n + 1, 0, NV)
        bad = (cb[hi] - cb[lo]) > 0
    return bad, ew


def hf_station_day(ns, date, node, pos_list):
    """E(10-19)/E(2-8) at the native rate for 12-s windows starting at pos_list (seconds after midnight,
    already including the template offset).  Returns (ratios array aligned with pos_list, ref median log10)."""
    from obspy import read, UTCDateTime
    out = np.full(len(pos_list), np.nan)
    p = W.raw_path(ns, date)
    if p is None:
        return out, np.nan
    band = node.split(".")[3]
    try:
        st = read(p).select(channel=band + "?")
    except Exception:
        return out, np.nan
    locs = sorted({tr.stats.location for tr in st})
    if not locs:
        return out, np.nan
    st = st.select(location="" if "" in locs else locs[0])
    for tr in st:
        fs = tr.stats.sampling_rate
        if abs(fs - round(fs)) / fs < 1e-4:
            tr.stats.sampling_rate = float(round(fs))
    D0 = UTCDateTime(date)
    E28 = E1019 = None; fs0 = None
    for ch in sorted({tr.stats.channel for tr in st}):
        s1 = st.select(channel=ch).copy()
        s1.merge(method=1, fill_value=0)
        tr = s1[0]
        fs0 = tr.stats.sampling_rate
        if fs0 < 40:
            return out, np.nan
        x = tr.data.astype(np.float64); x -= np.median(x)
        i0 = int(round((tr.stats.starttime - D0) * fs0))
        N = int(86400 * fs0)
        y = np.zeros(N)
        a, b = max(i0, 0), min(i0 + len(x), N)
        if b <= a:
            continue
        y[a:b] = x[a - i0:b - i0]
        e1 = W.bp(y, fs0, 2.0, 8.0) ** 2
        e2 = W.bp(y, fs0, W.HF_BAND[0], W.HF_BAND[1]) ** 2
        E28 = e1 if E28 is None else E28 + e1
        E1019 = e2 if E1019 is None else E1019 + e2
    if E28 is None:
        return out, np.nan
    c1 = np.concatenate([[0], np.cumsum(E28)]); c2 = np.concatenate([[0], np.cumsum(E1019)])
    nw = int(W.TPL_S * fs0)
    j = np.round(np.asarray(pos_list, float) * fs0).astype(np.int64)
    ok = (j >= 0) & (j + nw <= len(E28))
    s28 = np.where(ok, c1[np.clip(j + nw, 0, len(c1) - 1)] - c1[np.clip(j, 0, len(c1) - 1)], np.nan)
    s19 = np.where(ok, c2[np.clip(j + nw, 0, len(c2) - 1)] - c2[np.clip(j, 0, len(c2) - 1)], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(s28 > 0, s19 / s28, np.nan)
    # reference: random windows with non-zero energy
    rng = np.random.default_rng(int(pd.Timestamp(date).strftime("%Y%m%d")) + zlib.crc32(ns.encode()) % 1000)
    jr = rng.integers(0, len(E28) - nw, 4 * NRAND_HF)
    r28 = c1[jr + nw] - c1[jr]; r19 = c2[jr + nw] - c2[jr]
    k = (r28 > 0) & (r19 > 0)
    ref = np.nan
    if k.sum() >= 50:
        ref = float(np.median(np.log10(r19[k][:NRAND_HF] / r28[k][:NRAND_HF])))
    return out, ref


def day_job(ds):
    t0 = time.time()
    T, Q, SITE = _G["T"], _G["Q"], _G["SITE"]
    n = W.NT
    NV = W.NDAY - n + 1
    D0 = pd.Timestamp(ds).timestamp()
    J1 = json.load(open(f"{W.RES}/mf_{MF1}/run_{ds}.json"))
    p2 = f"{W.RES}/mf_{MF2}/run_{ds}.json"
    J2 = json.load(open(p2)) if os.path.exists(p2) else None
    SH = {"data": {}, "c1": J1["shift_site"], "c2": J1["shift_site2"]}
    DE = pd.read_parquet(f"{W.RES}/mf_{MF1}/det_{ds}.parquet")
    DE["series"] = DE.series.map({"data": "data", "control": "c1", "control2": "c2"})
    parts = [DE]
    chk = {}
    if J2 is not None:
        SH["c3"] = J2["shift_site"]; SH["c4"] = J2["shift_site2"]
        D2 = pd.read_parquet(f"{W.RES}/mf_{MF2}/det_{ds}.parquet")
        dd = D2[D2.series == "data"]
        a = DE[DE.series == "data"]
        chk = dict(n_data_mt1=len(a), n_data_c34run=len(dd),
                   data_identical=bool(len(a) == len(dd) and np.allclose(np.sort(a.t.values), np.sort(dd.t.values))))
        D2 = D2[D2.series != "data"].copy()
        D2["series"] = D2.series.map({"control": "c3", "control2": "c4"})
        parts.append(D2)
    C = pd.concat(parts, ignore_index=True)
    C["i"] = np.round((C.t.to_numpy() - D0) * W.FS).astype(np.int64)
    FD = pd.read_parquet(f"{W.RES}/mf_{MF1}/famday_{ds}.parquet")
    FD = FD[FD.observable]
    nodes = J1["nodes"]
    fam_sta = {int(f): s.split(",") for f, s in zip(FD.fid, FD.stations)}
    rng = np.random.default_rng(1000 + int(pd.Timestamp(ds).strftime("%Y%m%d")))
    R = rng.integers(0, NV, NRAND)                            # random ORIGIN samples (same for every family)
    Cg = {f: g for f, g in C.groupby("fid")}
    # station cc arrays: cand[f][ns] (n_cand_f,), rand[f][ns] (NRAND,)
    CC, RC, NCH = {}, {}, {}
    sta_fams = {}
    for f, ss in fam_sta.items():
        for ns in ss:
            sta_fams.setdefault(ns, []).append(f)
    for ns, fl in sta_fams.items():
        p = W.cache_path(ns, ds)
        z = np.load(p)
        node = str(z["node"])
        if node != nodes.get(ns):
            continue
        x = z["x"]; m = z["m"]; live = z["live"]
        bad, ew = guard_np(x, live, m)
        efl = 1e-10 * np.median(ew[live])
        site = SITE[ns]
        for f in fl:
            key = (f, ns, node)
            if key not in T or f not in Cg:
                continue
            tp = T[key]
            ch = tp["live"] & live
            chi = np.nonzero(ch)[0]
            if not len(chi):
                continue
            d = int(round(tp["off"] * W.FS))
            g = Cg[f]
            sh = np.array([int(SH[s].get(site, 0)) % NV if s != "data" else 0 for s in g.series], np.int64)
            j = ((g.i.to_numpy() - sh) % NV) + d
            jr = R + d
            allj = np.concatenate([j, jr])
            okj = (allj >= 0) & (allj < NV)
            jj = np.clip(allj, 0, NV - 1)
            idx = jj[:, None] + np.arange(n)[None]
            tpl = tp["tpl"][chi]                                   # (nch, n)
            num = np.einsum("cnk,ck->cn", x[chi[:, None, None], idx[None]], tpl)     # (nch, N)
            E = np.sqrt(np.maximum(ew[chi[:, None], jj[None]], efl)).astype(np.float32)
            cc = (num / E).mean(0)
            cc = np.where(okj & ~bad[jj], cc, np.nan)
            CC.setdefault(f, {})[ns] = cc[:len(j)]
            RC.setdefault(f, {})[ns] = cc[len(j):]
            NCH.setdefault(f, {})[ns] = len(chi)
        del x, ew, bad, z
    # ---------------- per candidate aggregates
    out = []
    for f, g in Cg.items():
        if f not in CC:
            continue
        sts = [ns for ns in fam_sta[f] if ns in CC[f]]
        node_of = {ns: nodes[ns] for ns in sts}
        wm = {ns: float(NCH[f][ns]) * T[(f, ns, node_of[ns])]["w"] for ns in sts}   # wa_mf: n_ch(live in both) x w
        groups = {}
        for ns in sts:
            groups.setdefault(SITE[ns], []).append(ns)
        sites = sorted(groups)
        Ws = np.array([np.sqrt(sum(wm[s] ** 2 for s in groups[q])) for q in sites])
        Wt = Ws.sum()
        ncand = len(g)
        SC = np.zeros((ncand, len(sites))); SR = np.zeros((NRAND, len(sites)))
        for k, q in enumerate(sites):
            mem = groups[q]; ws = np.array([wm[s] for s in mem])
            SC[:, k] = (np.stack([CC[f][s] for s in mem], 1) * ws).sum(1) / ws.sum()
            SR[:, k] = (np.stack([RC[f][s] for s in mem], 1) * ws).sum(1) / ws.sum()
        net_re = (SC * Ws).sum(1) / Wt
        med = np.nanmedian(SR, 0)
        mad = 1.4826 * np.nanmedian(np.abs(SR - med), 0)
        mad = np.where(mad > 0, mad, np.nan)
        Z = (SC - med) / mad
        contrib = SC * Ws / Wt
        Zs = np.sort(np.nan_to_num(Z, nan=-99), 1)[:, ::-1]
        Cs = np.sort(np.nan_to_num(SC, nan=-99), 1)[:, ::-1]
        ns_ = len(sites)
        with np.errstate(invalid="ignore", divide="ignore"):
            dom = np.nanmax(contrib, 1) / np.where(net_re > 0, net_re, np.nan)
        r = pd.DataFrame(dict(
            fid=f, t=g.t.to_numpy(), series=g.series.to_numpy(), net_cc=g.net_cc.to_numpy(), net_mad=g.net_mad.to_numpy(),
            net_re=net_re, n_sta=len(sts), n_sites=ns_, wsum=float(sum(wm.values())),
            site_cc_max=Cs[:, 0], site_cc_2nd=Cs[:, 1] if ns_ > 1 else np.nan,
            site_z_max=Zs[:, 0], site_z_2nd=Zs[:, 1] if ns_ > 1 else np.nan,
            site_z_min=np.nanmin(Z, 1),
            n_z2=(Z >= 2).sum(1), n_z3=(Z >= 3).sum(1), n_z5=(Z >= 5).sum(1), dom_top=dom,
            site_noise_mad_med=float(np.nanmedian(mad))))
        r.loc[r.site_cc_max <= -98, ["site_cc_max", "site_z_max"]] = np.nan
        r.loc[r.site_cc_2nd <= -98, "site_cc_2nd"] = np.nan
        r.loc[r.site_z_2nd <= -98, "site_z_2nd"] = np.nan
        out.append(r)
    F = pd.concat(out, ignore_index=True) if out else pd.DataFrame()
    # ---------------- HF ratio at the family-day's HF station
    F["hfr"] = np.nan; F["hfr_rel"] = np.nan; F["hf_sta"] = ""
    if len(F):
        QG = dict(tuple(Q.groupby("fid")))
        hfsta = {}
        for f, ss in fam_sta.items():
            qf = QG.get(f)
            if qf is None:
                continue
            c = qf[qf.netsta.isin(ss)]
            c = c[[nodes.get(a) == b for a, b in zip(c.netsta, c.node)]]
            if len(c):
                hfsta[f] = (c.iloc[0].netsta, c.iloc[0].node, float(c.iloc[0].offset_s))
        F["hf_sta"] = F.fid.map(lambda f: hfsta.get(f, ("",))[0])
        for ns, gi in F.groupby("hf_sta").groups.items():
            if not ns:
                continue
            gg = F.loc[gi]
            node = nodes[ns]
            site = SITE[ns]
            off = gg.fid.map(lambda f: hfsta[f][2]).to_numpy()
            i = np.round((gg.t.to_numpy() - D0) * W.FS).astype(np.int64)
            sh = np.array([int(SH[s].get(site, 0)) % NV if s != "data" else 0 for s in gg.series], np.int64)
            pos = ((i - sh) % NV) / W.FS + off
            hr, ref = hf_station_day(ns, ds, node, pos)
            F.loc[gi, "hfr"] = hr
            with np.errstate(divide="ignore", invalid="ignore"):
                F.loc[gi, "hfr_rel"] = np.log10(hr) - ref
    F["date"] = ds
    os.makedirs(OD, exist_ok=True)
    F.to_parquet(f"{OD}/feat_{ds}.parquet.tmp", index=False)
    os.replace(f"{OD}/feat_{ds}.parquet.tmp", f"{OD}/feat_{ds}.parquet")
    d = F[F.series == "data"]
    err = float(np.nanmax(np.abs(d.net_re - d.net_cc))) if len(d) else np.nan
    errc = float(np.nanmax(np.abs(F.net_re - F.net_cc))) if len(F) else np.nan
    return ds, len(F), err, errc, chk, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--days", default="")
    ap.add_argument("--windows", default="V2008,V2010,A2011")
    a = ap.parse_args()
    os.makedirs(OD, exist_ok=True)
    days = [d for w in a.windows.split(",") for d, _ in W.window_days(w)]
    if a.days:
        days = a.days.split(",")
    days = [d for d in days if os.path.exists(f"{W.RES}/mf_{MF1}/det_{d}.parquet")
            and not os.path.exists(f"{OD}/feat_{d}.parquet")]
    W.log(LOG, f"[feat] {len(days)} days, {a.workers} workers")
    with ProcessPoolExecutor(a.workers, initializer=_init) as ex:
        fut = {ex.submit(day_job, d): d for d in days}
        for fu in as_completed(fut):
            ds, nr, err, errc, chk, secs = fu.result()
            W.log(LOG, f"[feat] {ds}: {nr} candidates, max |net_re - net_cc| data {err:.2e} all {errc:.2e}, {chk}, {secs:.0f}s")


if __name__ == "__main__":
    main()
