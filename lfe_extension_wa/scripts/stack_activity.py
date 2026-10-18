#!/usr/bin/env python
"""Family ACTIVITY by stacking (engine) — lfe_extension_wa/results/stack_activity/.

Question: when is each WA family active?  Single detections on the permanent network cannot be told from
the scrambled control for most families; here we stack many low-threshold candidate times per family and
station, and ask whether the stack looks like the station's template.  Nothing is gated: every quantity is
stored and the scrambled control and random times get exactly the same treatment.

Per day, per family f (every family with templates on >= 2 sites that day):
  templates   ALL (family, station-epoch) stacks of >= 30 array-era events (wa_templates.py; clear flag kept
              as a column), first template tag listed wins, a station-day uses the template of ITS node only.
  network     frozen 1-8 Hz chain (wa_mf.py core): per-station normalised cc (mean over live channels),
              energy-floor / mask guard, template offset to the ORIGIN axis, site merge B
              (w = n_channels x SNR weight).  Deviation from wa_mf: a sample is averaged over the stations /
              sites VALID there (valid-weighted mean) instead of requiring every station valid (with ~34
              stations per family the all-valid rule would lose most of the day).
  candidates  LEAVE-ONE-SITE-OUT: for a station of site q the candidate times come from the network of the
              OTHER sites (valid where the valid other-site weight >= 50 % of their total), local maxima
              declustered +-4 s, TOP-K per day with K = 2 per scored hour (48 per full day; 23 on the
              test half of a validation day).  Fixed before any scoring.  So a station's stack is built at
              times chosen without its own data.
  stack       at each candidate time the station's 2-8 Hz waveform (Butterworth-4 zero-phase, |H|^2 applied
              in the frequency domain, designed at fs = 20 Hz = the cache rate) is cut over the template
              window (12 s from origin + template offset), each window unit-max normalised over its 3
              components (template-building rule), windows touching invalid data skipped; the day's sum S
              (3 x 240), the summed per-channel energy e and the count are stored -> any window W of days
              is the sum of its days.
  series      data; control1..3 = SITE-ROLLED scrambled control (each site's cc series rolled by its own
              shift in [300, 3000] s, independent draws per control, candidates picked on the rolled
              network, the station's data stacked at the rolled positions); random1..6 = K uniform random
              times from the same valid, scored samples of the same day (the random-times reference).
  scored      --scored test: origin AND station-data times in the TEST half (2-h blocks, 60-s guards; the
              templates were built on the build half); --scored all.

  python scripts/stack_activity.py --tag valA --windows V2008,V2010 --tpl build --scored test
  -> results/stack_activity/run_<tag>/day_<date>.npz  (rows x 10 series stacks), fam_<date>.parquet,
     cand_<date>.parquet (full-network top-K of data / control1, diagnostic), run_<date>.json
"""
from __future__ import annotations
import os, sys, json, time, argparse
from concurrent.futures import ThreadPoolExecutor
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

OUT = f"{W.RES}/stack_activity"
LOG = f"{W.STUDY}/logs/stack_activity.log"
SERIES = ["data", "control1", "control2", "control3"] + [f"random{i}" for i in range(1, 7)]
N_CTL, N_RAND = 3, 6
DECL = int(4.0 * W.FS)          # +-4 s declustering
NL = int(W.LONG_LEN * W.FS)
IPRE = int(W.LONG_PRE * W.FS)


def load_templates(tags):
    """all >= 30-event templates: 1-8 Hz (MF) + 2-8 Hz (stack correlation) from the long stack."""
    rows, T18, T28 = [], {}, {}
    for tag in tags:
        q = pd.read_csv(f"{W.RES}/templates/tpl_qc_{tag}.csv", dtype={"live": str})
        z = np.load(f"{W.RES}/templates/tpl_{tag}.npz")
        zs = np.load(f"{W.RES}/templates/stacks_{tag}.npz")
        kt = {k: i for i, k in enumerate(z["keys"])}  # noqa
        ks = {k: i for i, k in enumerate(zs["keys"])}
        stk = zs["stacks"]
        tpl_arr = z["tpl"]                      # NpzFile re-reads on every [] access: load once
        q = q[q.n_events >= W.MIN_EVENTS].copy()
        q["live"] = q.live.astype(str).str.zfill(3)
        for r in q.itertuples():
            key = (int(r.fid), r.netsta, r.node)
            k = f"{r.fid}|{r.netsta}|{r.node}"
            if key in T18 or k not in kt or k not in ks or not np.isfinite(r.offset_s):
                continue
            live = np.array([c == "1" for c in r.live])
            j0 = IPRE + int(round(r.offset_s * W.FS))
            st = stk[ks[k]].astype(np.float64)
            t28 = W.bp(st, W.FS, 2.0, 8.0)[:, j0:j0 + W.NT].copy()
            t28[~live] = 0
            nr = np.sqrt((t28 ** 2).sum(1))
            ok = live & (nr > 0)
            t28[ok] /= nr[ok, None]
            t28[~ok] = 0
            T18[key] = tpl_arr[kt[k]].astype(np.float32)
            T28[key] = t28.astype(np.float32)
            w = float(r.w) if np.isfinite(r.w) else 1.0
            rows.append(dict(fid=int(r.fid), netsta=r.netsta, node=r.node, tag=tag, clear=bool(r.clear),
                             n_events=int(r.n_events), imp=r.imp, dur=r.dur, snr_db=r.snr_db, w=w,
                             offset_s=float(r.offset_s), live=r.live))
    return pd.DataFrame(rows), T18, T28


def load_day(stations, ds):
    D = {}
    for ns in stations:
        p = W.cache_path(ns, ds)
        if os.path.exists(p):
            try:
                z = np.load(p)
                D[ns] = dict(x=z["x"], m=z["m"], live=z["live"], node=str(z["node"]))
            except Exception:
                pass
    return D


def main():
    import torch
    import torch.nn.functional as Fnn
    from scipy.signal import butter, sosfreqz
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--windows", default="")
    ap.add_argument("--date-range", default="", help="YYYY-MM-DD:YYYY-MM-DD (continuous series)")
    ap.add_argument("--days", default="", help="comma list (smoke)")
    ap.add_argument("--tpl", required=True)
    ap.add_argument("--stations", default=f"{W.RES}/setup/stations.csv")
    ap.add_argument("--families", default=f"{W.RES}/setup/families.csv")
    ap.add_argument("--fids", default="", help="comma list or csv with fid")
    ap.add_argument("--only-stations", default="", help="comma list or csv with netsta")
    ap.add_argument("--scored", choices=["test", "all"], required=True)
    ap.add_argument("--kph", type=float, default=2.0, help="candidates per scored hour (top-K per day)")
    ap.add_argument("--energy-floor", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=20260930)
    ap.add_argument("--sta-chunk", type=int, default=12)
    a = ap.parse_args()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    od = f"{OUT}/run_{a.tag}"
    os.makedirs(od, exist_ok=True)
    S = pd.read_csv(a.stations)
    if a.only_stations:
        keep = (set(pd.read_csv(a.only_stations).netsta) if a.only_stations.endswith(".csv")
                else set(a.only_stations.split(",")))
        S = S[S.netsta.isin(keep)]
    SITE = dict(zip(S.netsta, S.site))
    F = pd.read_csv(a.families)
    if a.fids:
        fk = set(pd.read_csv(a.fids).fid) if a.fids.endswith(".csv") else {int(x) for x in a.fids.split(",")}
        F = F[F.fid.isin(fk)]
    Q, T18, T28 = load_templates(a.tpl.split(","))
    Q = Q[Q.netsta.isin(SITE) & Q.fid.isin(set(F.fid))].reset_index(drop=True)
    Q["site"] = Q.netsta.map(SITE)
    Q.to_csv(f"{od}/templates_used.csv", index=False)
    TQ = {(r.fid, r.netsta, r.node): r for r in Q.itertuples()}
    stations = sorted(Q.netsta.unique())
    if a.date_range:
        d0, d1 = a.date_range.split(":")
        days = [(d.strftime("%Y-%m-%d"), "cont") for d in pd.date_range(d0, d1, freq="D")]
    else:
        days = [(d, k) for w in a.windows.split(",") for d, k in W.window_days(w)]
    if a.days:
        ks = set(a.days.split(","))
        days = [x for x in days if x[0] in ks]
    n = W.NT
    NV = W.NDAY - n + 1
    NF = int(2 ** np.ceil(np.log2(W.NDAY + n)))
    sos = butter(4, [2.0, 8.0], btype="band", fs=W.FS, output="sos")
    _, h = sosfreqz(sos, worN=np.fft.rfftfreq(NF, 1.0 / W.FS), fs=W.FS)
    H2 = torch.from_numpy((np.abs(h) ** 2).astype(np.complex64)).to(dev)
    ar = torch.arange(n, device=dev)
    NB = (NV + DECL - 1) // DECL
    json.dump(dict(vars(a), n_templates=len(Q), n_clear=int(Q.clear.sum()), n_fam=int(Q.fid.nunique()),
                   stations=stations, series=SERIES, decluster_s=4.0, NF=NF),
              open(f"{od}/config.json", "w"), indent=1)
    W.log(LOG, f"[sa:{a.tag}] {len(days)} days, {len(Q)} templates ({int(Q.clear.sum())} clear) on "
          f"{len(stations)} stations, {Q.fid.nunique()} families, scored={a.scored}, kph={a.kph}, dev {dev}")
    todo = [(ds, k) for ds, k in days if not os.path.exists(f"{od}/day_{ds}.npz")]
    ex = ThreadPoolExecutor(1)
    fut = ex.submit(load_day, stations, todo[0][0]) if todo else None
    for it, (ds, kind) in enumerate(todo):
        t_day = time.time()
        Dn = fut.result()
        if it + 1 < len(todo):
            fut = ex.submit(load_day, stations, todo[it + 1][0])
        D0 = pd.Timestamp(ds).timestamp()
        tt = D0 + np.arange(NV) / W.FS
        if a.scored == "test":
            scored_np = W.half_of(tt) == 1
        else:
            scored_np = np.ones(NV, bool)
        scored = torch.from_numpy(scored_np).to(dev)
        K = int(round(a.kph * scored_np.sum() / (3600 * W.FS)))
        # ---------- stations of the day (GPU)
        P, gi = {}, {}
        X28_all = torch.empty((max(len(Dn), 1), 3, W.NDAY), device=dev)
        BADU = torch.empty((max(len(Dn), 1), NV), dtype=torch.bool, device=dev)
        for ns, d in Dn.items():
            X = torch.from_numpy(d["x"]).to(dev)
            live = d["live"]
            DH = torch.fft.rfft(X, n=NF)
            gi[ns] = len(gi)
            X28_all[gi[ns]] = torch.fft.irfft(DH * H2[None], n=NF)[:, :W.NDAY]
            xd = X.double()
            EN = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64), xd * xd], 1), 1)
            ew = (EN[:, n:] - EN[:, :-n])[:, :NV]
            del EN, xd
            lv = torch.from_numpy(live).to(dev)
            ewl = ew[lv]
            med_ch = ewl.median(dim=1).values
            e3 = ewl.sum(0)
            med3 = e3.median()
            ref_ch = torch.maximum(med_ch, 1e-2 * med3)
            bad = ((ewl < a.energy_floor * ref_ch[:, None]).any(0) | (e3 < a.energy_floor * med3)
                   | (ewl.min(0).values <= 0))
            if float(med3) <= 0:
                bad = torch.ones_like(bad)
            mk = torch.from_numpy(~d["m"]).to(dev)
            cmk = torch.cat([torch.zeros(1, device=dev), torch.cumsum(mk.float(), 0)])
            bad |= (cmk[n:n + NV] - cmk[:NV]) > 0
            BADU[gi[ns]] = bad
            if bool(bad.any()):
                bad = Fnn.max_pool1d(bad.float()[None, None], kernel_size=2 * n + 1, stride=1, padding=n)[0, 0] > 0
            efl = 1e-10 * ewl.median()
            E = ew.clamp_min(efl).sqrt().float()
            del ew, ewl, e3, X
            P[ns] = dict(DH=DH, E=E, bad_dil=bad, live=live, node=d["node"])
        del Dn
        rng = np.random.default_rng(a.seed + int(pd.Timestamp(ds).strftime("%Y%m%d")))
        sites_today = sorted({SITE[ns] for ns in P})
        shifts = [{q: int(rng.uniform(*W.ROLL_S) * W.FS) for q in sites_today} for _ in range(N_CTL)]
        g = torch.Generator(device=dev)
        g.manual_seed(a.seed + int(pd.Timestamp(ds).strftime("%Y%m%d")))
        out_rows, out_S, out_e, out_c, fam_rows, cand_rows = [], [], [], [], [], []
        for f in sorted(Q.fid.unique()):
            ss = []
            for ns in P:
                key = (f, ns, P[ns]["node"])
                if key in TQ:
                    tl = np.array([c == "1" for c in TQ[key].live])
                    ch = tl & P[ns]["live"]
                    if ch.any():
                        ss.append((ns, key, ch))
            sites = sorted({SITE[ns] for ns, _, _ in ss})
            fr = dict(fid=f, date=ds, kind=kind, n_sta=len(ss), n_sites=len(sites), K=K,
                      n_sta_clear=sum(bool(TQ[k].clear) for _, k, _ in ss),
                      n_sites_clear=len({SITE[ns] for ns, k, _ in ss if TQ[k].clear}),
                      stations=",".join(ns for ns, _, _ in ss))
            if len(sites) < 2:
                fam_rows.append(fr)
                continue
            si = {q: i for i, q in enumerate(sites)}
            nq = len(sites)
            wm = {ns: float(ch.sum()) * float(TQ[k].w) for ns, k, ch in ss}
            Wq = torch.zeros(nq, device=dev, dtype=torch.float32)
            for ns, _, _ in ss:
                Wq[si[SITE[ns]]] += wm[ns] ** 2
            Wq = Wq.sqrt()
            snum = torch.zeros((nq, NV), device=dev)
            sden = torch.zeros((nq, NV), device=dev)
            for i0 in range(0, len(ss), a.sta_chunk):
                cs = ss[i0:i0 + a.sta_chunk]
                Tt = torch.from_numpy(np.stack([T18[k] for _, k, _ in cs])).to(dev)
                WH = torch.fft.rfft(Tt, n=NF)
                DHb = torch.stack([P[ns]["DH"] for ns, _, _ in cs])
                num = torch.fft.irfft(WH.conj() * DHb, n=NF)[:, :, :NV]
                del WH, DHb, Tt
                chm = torch.from_numpy(np.stack([ch for _, _, ch in cs]).astype(np.float32)).to(dev)
                for j, (ns, k, ch) in enumerate(cs):
                    cc = (num[j] / P[ns]["E"] * chm[j][:, None]).sum(0) / chm[j].sum()
                    cc[P[ns]["bad_dil"]] = 0.0
                    dd = int(round(TQ[k].offset_s * W.FS))
                    a0, b0 = max(0, -dd), min(NV, NV - dd)
                    q = si[SITE[ns]]
                    tv = torch.zeros(NV, device=dev)
                    tv[a0:b0] = (~P[ns]["bad_dil"][a0 + dd:b0 + dd]).float()
                    tmp = torch.zeros(NV, device=dev)
                    tmp[a0:b0] = cc[a0 + dd:b0 + dd]
                    snum[q] += wm[ns] * tv * tmp
                    sden[q] += wm[ns] * tv
                    del cc, tv, tmp
                del num
            vq = sden > 0
            mq = torch.where(vq, snum / sden.clamp_min(1e-20), torch.zeros_like(snum))
            del snum, sden
            Wtot = Wq.sum()
            need = 0.5 * (Wtot - Wq)[:, None]

            def lo_net(m, v):
                Wv = Wq[:, None] * v.float()
                wm_ = Wv * m
                num_all = wm_.sum(0); den_all = Wv.sum(0)
                den = den_all[None] - Wv
                net = (num_all[None] - wm_) / den.clamp_min(1e-20)
                ok = (den >= need) & scored[None]
                full = num_all / den_all.clamp_min(1e-20)
                okf = (den_all >= 0.5 * Wtot) & scored
                return net, ok, full, okf

            def pick(x, ok, k):
                # declustering +-4 s: maxima of 4-s blocks that are >= the left and > the right neighbouring
                # block maxima (so every pick is the maximum within >= +-4 s and picks are >= 4 s apart)
                nr = x.shape[0]
                xp = torch.full((nr, NB * DECL), -float("inf"), device=dev)
                xp[:, :NV] = torch.where(ok, x, torch.full_like(x, -float("inf")))
                bm, ba = xp.view(nr, NB, DECL).max(2)
                del xp
                ninf = torch.full((nr, 1), -float("inf"), device=dev)
                left = torch.cat([ninf, bm[:, :-1]], 1); right = torch.cat([bm[:, 1:], ninf], 1)
                bm2 = torch.where((bm >= left) & (bm > right), bm, torch.full_like(bm, -float("inf")))
                v, ib = torch.topk(bm2, min(k, NB), dim=1)
                return ib * DECL + ba.gather(1, ib), torch.isfinite(v), v

            idx_s, good_s = [], []                 # per series: (nq, K) data-origin index, valid
            net, ok, full, okf = lo_net(mq, vq)
            ix, gd, vv = pick(net, ok, K)
            idx_s.append(ix); good_s.append(gd)
            fx, fg, fv = pick(full[None], okf[None], K)
            tsel = fx[0][fg[0]].cpu().numpy()
            for t_, v_ in zip(tsel, fv[0][fg[0]].cpu().numpy()):
                cand_rows.append((f, "data", D0 + t_ / W.FS, float(v_)))
            fr["cand_mean_data"] = float(vv[gd].mean()) if bool(gd.any()) else np.nan
            ok_data = ok
            del net, full, okf
            for r in range(N_CTL):
                shl = [shifts[r][q] for q in sites]
                shq = torch.tensor(shl, device=dev)
                mr = torch.stack([torch.roll(mq[i], shl[i]) for i in range(nq)])
                vr = torch.stack([torch.roll(vq[i], shl[i]) for i in range(nq)])
                net, okr, full, okf = lo_net(mr, vr)
                ix, gd, vv = pick(net, okr, K)
                idx_s.append(torch.remainder(ix - shq[:, None], NV)); good_s.append(gd)
                if r == 0:
                    fx, fg, fv = pick(full[None], okf[None], K)
                    tsel = fx[0][fg[0]].cpu().numpy()
                    for t_, v_ in zip(tsel, fv[0][fg[0]].cpu().numpy()):
                        cand_rows.append((f, "control1", D0 + t_ / W.FS, float(v_)))
                fr[f"cand_mean_control{r+1}"] = float(vv[gd].mean()) if bool(gd.any()) else np.nan
                del mr, vr, net, okr, full, okf
            # random times: uniform draws, the first K that fall on valid scored samples of the leave-out network
            M = 8 * K
            cr = torch.randint(0, NV, (nq, N_RAND * M), device=dev, generator=g)
            okc = ok_data.gather(1, cr).view(nq, N_RAND, M)
            key = torch.where(okc, (M - torch.arange(M, device=dev)).expand(nq, N_RAND, M), torch.full_like(cr.view(nq, N_RAND, M), -1))
            kv, pos = key.topk(K, dim=2)
            rix = cr.view(nq, N_RAND, M).gather(2, pos)
            for r in range(N_RAND):
                idx_s.append(rix[:, r]); good_s.append(kv[:, r] > 0)
            del cr, okc, key, ok_data, mq, vq
            IX = torch.stack(idx_s, 1)             # (nq, 10, K)
            GD = torch.stack(good_s, 1)
            # ---------- stacks per station
            qs = torch.tensor([si[SITE[ns]] for ns, _, _ in ss], device=dev)
            gs = torch.tensor([gi[ns] for ns, _, _ in ss], device=dev)
            dds = torch.tensor([int(round(TQ[k].offset_s * W.FS)) for _, k, _ in ss], device=dev)
            Iv = IX[qs]                                                # (nsta, 10, K)
            J = Iv + dds[:, None, None]
            okw = GD[qs] & (J >= 0) & (J < NV)
            Jc = J.clamp(0, NV - 1)
            okw &= (~BADU[gs[:, None, None], Jc]) & scored[Iv]
            base = (gs[:, None] * 3 + torch.arange(3, device=dev)[None]) * W.NDAY      # (nsta, 3)
            flat = base[:, None, None, :, None] + Jc[..., None, None] + ar          # (nsta, 10, K, 3, n)
            Xw = X28_all.view(-1)[flat]
            del flat
            mx = Xw.abs().amax(dim=(3, 4))
            okw &= mx > 0
            Xw = Xw / torch.where(okw, mx, torch.ones_like(mx))[..., None, None]
            Xw = Xw * okw[..., None, None].float()
            out_S.append(Xw.sum(2).half().cpu().numpy())               # (nsta, 10, 3, n)
            out_e.append((Xw ** 2).sum((2, 4)).cpu().numpy().astype(np.float32))
            out_c.append(okw.sum(2).cpu().numpy().astype(np.int16))
            out_rows += [(f, ns, k[2], bool(TQ[k].clear)) for ns, k, _ in ss]
            del Xw
            fam_rows.append(fr)
            del IX, GD
        del P, X28_all, BADU
        R = pd.DataFrame(out_rows, columns=["fid", "netsta", "node", "clear"])
        np.savez(f"{od}/day_{ds}.tmp.npz", fid=R.fid.to_numpy(), netsta=R.netsta.to_numpy().astype(str),
                 node=R.node.to_numpy().astype(str), clear=R.clear.to_numpy(),
                 S=np.concatenate(out_S) if out_S else np.zeros((0, len(SERIES), 3, n), np.float16),
                 e=np.concatenate(out_e) if out_e else np.zeros((0, len(SERIES), 3), np.float32),
                 cnt=np.concatenate(out_c) if out_c else np.zeros((0, len(SERIES)), np.int16),
                 series=np.array(SERIES), K=K, date=ds, kind=kind)
        pd.DataFrame(fam_rows).to_parquet(f"{od}/fam_{ds}.parquet", index=False)
        pd.DataFrame(cand_rows, columns=["fid", "series", "t", "val"]).to_parquet(f"{od}/cand_{ds}.parquet", index=False)
        os.replace(f"{od}/day_{ds}.tmp.npz", f"{od}/day_{ds}.npz")
        nf = sum(1 for x in fam_rows if x["n_sites"] >= 2)
        json.dump(dict(date=ds, kind=kind, n_sta=len(sites_today), n_fam=nf, n_rows=len(R), K=K,
                       secs=time.time() - t_day, shifts=[{k_: v for k_, v in s.items()} for s in shifts]),
                  open(f"{od}/run_{ds}.json", "w"), indent=1)
        if dev == "cuda":
            torch.cuda.empty_cache()
        W.log(LOG, f"[sa:{a.tag}] {ds} {kind}: {len(sites_today)} sites, {nf} families, {len(R)} rows, K {K}, "
              f"{time.time()-t_day:.0f}s")
    ex.shutdown()


if __name__ == "__main__":
    main()
