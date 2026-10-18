#!/usr/bin/env python
"""SH_CONF stage 3 (GPU): PNET S probability at every candidate's expected S time, + a Mendocino fine-tune.

PNET = lfe_picker P1b PhaseNet-style U-Net trained on Bostock SVI LFEs (pnet_w05.pt; input contract
lfe_picker/scripts/p_scan.cf_day: 60-s windows, 20 Hz, 1-8 Hz Butterworth-4 zero-phase chain, Z/N/E, per-window 3-C
peak normalisation).  The lfe_south cache20hz IS that chain (t3_cache.py), recorded order (alphabetical channels):
  * horizontals rotated to N/E with the FDSN azimuths of the day's epoch (results/sh_ext2026/station_epochs_raw.csv;
    PB boreholes EH1/EH2, BK.BRIC 7/97 deg before 2024-05-21, the rest ~0/90);
  * vertical-only stations (NC.KPR, NC.KBU EHZ) enter as (Z, 0, 0) -- the SVI model never saw this input.
Expected S at a station = candidate origin + tt_s + best_tau of that family's template at that station (the MF's
pick); a scrambled-control candidate is read at its site's ROLLED position (sh_conf_common.s_index), exactly as the
control MF built it.  Station value = max CF within +-0.5 s; site value = max over the site's stations.
Features (suffix _w = off-the-shelf, _f = Mendocino fine-tune): pnet_max, pnet_2nd, n_sites_pnet03, n_sites_pnet05,
pnet_mean, pnet_rel_max (site value - that site's median at NRAND random origins of the same family-day), per-station
values pn<sfx>_<NET.STA>.

Fine-tune (the WA recipe, lfe_extension_wa/scripts/pnet_labels.py + pnet_finetune.py): start = SVI weights; windows
ONLY on v3b template-BUILD days (never a fit / evaluation day: every build day of any family is excluded from the
confidence fit and evaluation) and ONLY at the 14 training stations (PNET_HO = BK.JEWT, NC.KBU, PB.B933 never used).
Positives = every Shelly (2025) event of that day at stations where that family's template S pick has pick_q >= 3
(weight 1 if pick_q >= 10, else 0.5), Gaussian sigma 0.5 s at every labelled S in the window, loss weight only within
+-3 s of labelled S; negatives = random windows of the same station-days, weight 1 except +-3 s around ANY catalogued
event's S (any family).  Masked BCE, 1/2 positive + 1/2 negative batches, horizontals zeroed with p = 0.3 on 3-C
windows, polarity flip; internal validation = every 10th build date; best internal AUC checkpoint kept (step 0 =
the SVI weights is a candidate).

  python scripts/sh_conf_pnet.py labels --workers 8
  python scripts/sh_conf_pnet.py finetune --steps 6000
  python scripts/sh_conf_pnet.py feat --shard 0 --nshard 3
"""
import os, sys, json, time, argparse, warnings, traceback
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sh_conf_common as C
PICKER = "/home/jovyan/tremorferometry/lfe_picker/scripts"
PNET_W05 = "/home/jovyan/tremorferometry/lfe_picker/results/p1/models/pnet_w05.pt"
PD = f"{C.OUT}/pnet"
NW = 1200; SIG = 10.0; MASK = 60
PN = ["pnet_max", "pnet_2nd", "n_sites_pnet03", "n_sites_pnet05", "pnet_mean", "pnet_rel_max"]
_AZ = {}


def az_table():
    if not _AZ:
        E = pd.read_csv(f"{C.EXT}/station_epochs_raw.csv")
        E["t0"] = pd.to_datetime(E.start)
        en = E.end.fillna("").astype(str)
        E["t1"] = pd.to_datetime(en.where(~(en.eq("") | en.str.startswith("3000") | en.eq("nan")), "2100-01-01T00:00:00"))
        _AZ["E"] = E
    return _AZ["E"]


def azimuth(ns, cha, date):
    E = az_table(); d = pd.Timestamp(date) + pd.Timedelta(hours=12)
    g = E[(E.netsta == ns) & (E.cha == cha) & (E.t0 <= d) & (E.t1 > d)]
    return float(g.az.iloc[0]) if len(g) else np.nan


def load_zne(ns, date):
    """cache (C, N) -> (3, N) float32 Z, N, E; status.  (Z, 0, 0) for vertical-only stations."""
    net, sta = ns.split(".")
    p = C.npz_path(net, sta, date)
    if not os.path.exists(p):
        return None, "nofile"
    try:
        z = np.load(p); x = z["x"].astype(np.float32); t0 = float(z["t0"]); ch = [str(c) for c in z["chans"]]
    except Exception as e:
        return None, "error"
    if x.ndim != 2 or x.shape[1] < 0.5 * 86400 * C.FS or not np.any(x):
        return None, "short"
    zi = [k for k, c in enumerate(ch) if c.endswith("Z")]
    hi = [k for k, c in enumerate(ch) if not c.endswith("Z")]
    y = np.zeros((3, x.shape[1]), np.float32)
    if zi: y[0] = x[zi[0]]
    nh = 0
    if len(hi) == 2 and np.any(x[hi[0]]) and np.any(x[hi[1]]):
        a1, a2 = azimuth(ns, ch[hi[0]], date), azimuth(ns, ch[hi[1]], date)
        if not (np.isfinite(a1) and np.isfinite(a2)):
            names = (ch[hi[0]][-1], ch[hi[1]][-1])
            a1, a2 = {"N": 0.0, "E": 90.0}.get(names[0], np.nan), {"N": 0.0, "E": 90.0}.get(names[1], np.nan)
        if np.isfinite(a1) and np.isfinite(a2):
            r1, r2 = np.radians(a1), np.radians(a2)
            A = np.array([[np.cos(r1), np.sin(r1)], [np.cos(r2), np.sin(r2)]])
            Ri = np.linalg.inv(A).astype(np.float32)
            y[1:3] = Ri @ x[hi]
            nh = 2
    return dict(x=y, t0=t0, n_horiz=nh, chans=ch), "ok"


def maxfilt(v, half=10):
    from scipy.ndimage import maximum_filter1d
    return maximum_filter1d(v, size=2 * half + 1, mode="constant", cval=0.0)


# ------------------------------------------------------------------------------------------------ labels
_G = {}
def _init_labels():
    M, _ = C.load_tpl(C.TPL["v3b"])
    S = C.his_catalogue(); S["date"] = S.t.dt.strftime("%Y-%m-%d")
    _G.update(M=M, TR={d: g for d, g in S.groupby("date")})


def label_job(args):
    ns, ds, cap_pos, n_neg, seed = args
    z, st = load_zne(ns, ds)
    if z is None:
        return ns, ds, st, None
    M = _G["M"]; x = z["x"]; t0 = z["t0"]
    rng = np.random.default_rng(seed)
    tg = _G["TR"].get(ds)
    lab_t, lab_w, msk_t = [], [], []
    m = M[M.netsta == ns].set_index("fam")
    if tg is not None:
        for f, g in tg.groupby("fam"):
            if f not in m.index: continue
            r = m.loc[f]
            s = (g.t_epoch.to_numpy() + float(r.tt_s) + float(r.best_tau) - t0) * C.FS
            msk_t.append(s)
            q = float(r.pick_q) if np.isfinite(r.pick_q) else 0.0
            if q >= 3:
                lab_t.append(s); lab_w.append(np.full(len(s), 1.0 if q >= 10 else 0.5))
    lab_t = np.concatenate(lab_t) if lab_t else np.zeros(0); lab_w = np.concatenate(lab_w) if lab_w else np.zeros(0)
    msk_t = np.sort(np.concatenate(msk_t)) if msk_t else np.zeros(0)
    o = np.argsort(lab_t); lab_t, lab_w = lab_t[o], lab_w[o]
    N = x.shape[1]
    valid = np.abs(x).sum(0) > 0
    cv = np.concatenate([[0], np.cumsum(~valid)])
    X, Y, W, K = [], [], [], []
    tt = np.arange(NW)

    def add(s0, kind):
        s0 = int(s0)
        if s0 < 0 or s0 + NW > N or cv[s0 + NW] - cv[s0] > NW // 10:
            return False
        seg = x[:, s0:s0 + NW]; mx = np.abs(seg).max()
        if not mx > 0: return False
        y = np.zeros(NW, np.float32); w = np.zeros(NW, np.float32)
        if kind == 0:
            lo, hi = np.searchsorted(lab_t, s0 - MASK), np.searchsorted(lab_t, s0 + NW + MASK)
            for tl, wl in zip(lab_t[lo:hi], lab_w[lo:hi]):
                rr = tl - s0
                y = np.maximum(y, np.exp(-0.5 * ((tt - rr) / SIG) ** 2)); mm = np.abs(tt - rr) <= MASK
                w[mm] = np.maximum(w[mm], wl)
        else:
            w[:] = 1.0
            ml, mh = np.searchsorted(msk_t, s0 - MASK), np.searchsorted(msk_t, s0 + NW + MASK)
            for tm in msk_t[ml:mh]:
                w[np.abs(tt - (tm - s0)) <= MASK] = 0.0
        w[~valid[s0:s0 + NW]] = 0.0
        X.append((seg / mx).astype(np.float16)); Y.append((y * 255).round().astype(np.uint8))
        W.append((w * 255).round().astype(np.uint8)); K.append(kind)
        return True
    if len(lab_t):
        npos = 0
        for i in rng.permutation(len(lab_t))[:4 * cap_pos]:
            if npos >= cap_pos: break
            if add(lab_t[i] - rng.uniform(5, 55) * C.FS, 0): npos += 1
    k = 0
    for s0 in rng.integers(0, N - NW, 6 * n_neg):
        if k >= n_neg: break
        if add(s0, 1): k += 1
    if not X:
        return ns, ds, "ok", None
    return ns, ds, "ok", (np.stack(X), np.stack(Y), np.stack(W), np.array(K, np.int8), z["n_horiz"])


def labels(a):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    BF = C.build_fam_days(); bdays = sorted(set().union(*BF.values()))
    if a.days: bdays = [d for d in bdays if d in set(a.days.split(","))]
    stas = [s for s in C.STAS if s not in C.PNET_HO]
    jobs = [(ns, d, a.cap_pos, a.n_neg, 7 + i * 31 + j) for i, d in enumerate(bdays) for j, ns in enumerate(stas)]
    C.log(f"[labels] {len(jobs)} station-days on {len(bdays)} build days, {len(stas)} training stations (held out {C.PNET_HO})")
    X, Y, Wt, meta, status = [], [], [], [], []
    t0 = time.time()
    with ProcessPoolExecutor(a.workers, initializer=_init_labels) as ex:
        fut = [ex.submit(label_job, j) for j in jobs]
        for i, fu in enumerate(as_completed(fut)):
            ns, ds, st, r = fu.result()
            status.append(("build", ds, ns, st, ""))
            if r is None: continue
            x, y, w, k, nh = r
            X.append(x); Y.append(y); Wt.append(w); meta.append(pd.DataFrame(dict(netsta=ns, date=ds, kind=k, n_horiz=nh)))
            if (i + 1) % 1000 == 0: C.log(f"[labels] {i + 1}/{len(jobs)} {time.time() - t0:.0f}s")
    C.write_status("pnetlabels", "0", status)
    od = f"{PD}/train"; os.makedirs(od, exist_ok=True)
    Mt = pd.concat(meta, ignore_index=True)
    np.save(f"{od}/X.npy", np.concatenate(X)); np.save(f"{od}/Y.npy", np.concatenate(Y)); np.save(f"{od}/Wt.npy", np.concatenate(Wt))
    Mt.to_parquet(f"{od}/meta.parquet", index=False)
    C.log(f"[labels] {len(Mt)} windows {Mt.groupby(['kind', 'n_horiz']).size().to_dict()}; stations {sorted(Mt.netsta.unique())}; "
          f"dates {Mt.date.nunique()} ({Mt.date.min()}..{Mt.date.max()}); {time.time() - t0:.0f}s")
    C.guard("pnetlabels")


# ------------------------------------------------------------------------------------------------ fine-tune
def masked_bce(p, y, w):
    import torch
    p = p.clamp(1e-6, 1 - 1e-6)
    b = -(y * torch.log(p) + (1 - y) * torch.log(1 - p)) * w
    return (b.sum(-1) / w.sum(-1).clamp_min(1e-6)).mean()


def auc(pos, neg):
    from scipy.stats import rankdata
    pos = np.asarray(pos, float); neg = np.asarray(neg, float)
    pos = pos[np.isfinite(pos)]; neg = neg[np.isfinite(neg)]
    if len(pos) < 5 or len(neg) < 5: return np.nan
    r = rankdata(np.concatenate([pos, neg]))
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def finetune(a):
    import torch
    sys.path.insert(0, PICKER); import p_scan
    torch.manual_seed(0); rng = np.random.default_rng(0); dev = "cuda"
    od = f"{PD}/train"
    X = np.load(f"{od}/X.npy", mmap_mode="r"); Y = np.load(f"{od}/Y.npy"); Wt = np.load(f"{od}/Wt.npy")
    Mt = pd.read_parquet(f"{od}/meta.parquet")
    assert not Mt.netsta.isin(C.PNET_HO).any(), "held-out station in the training windows"
    dates = sorted(Mt.date.unique()); vd = set(dates[5::10])
    val = Mt.date.isin(vd).to_numpy() if len(dates) >= 10 else (rng.random(len(Mt)) < 0.1)   # tiny smoke sets: random 10 %
    kind = Mt.kind.to_numpy(); nh = Mt.n_horiz.to_numpy()
    idx = {k: np.nonzero((kind == k) & ~val)[0] for k in (0, 1)}
    vidx = np.nonzero(val)[0]
    if len(vidx) > 30000: vidx = np.sort(rng.choice(vidx, 30000, replace=False))
    C.log(f"[ft] windows {len(Mt)}; train pos/neg {[len(idx[k]) for k in (0, 1)]}; val {len(vidx)} ({len(vd)} dates)")
    Xg = torch.from_numpy(np.ascontiguousarray(X)).to(dev); Yg = torch.from_numpy(Y).to(dev); Wg = torch.from_numpy(Wt).to(dev)
    NHg = torch.from_numpy(nh).to(dev)
    model, meta = p_scan.load_model(PNET_W05); model.train()
    opt = torch.optim.Adam(model.parameters(), lr=a.lr)
    nb = [a.batch // 2, a.batch - a.batch // 2]

    def batch_of(ii, aug):
        x = Xg[ii].float()
        if aug:
            zo = (torch.rand(len(ii), device=dev) < 0.3) & (NHg[ii] > 0); x[zo, 1:] = 0
            sg = torch.where(torch.rand(len(ii), device=dev) < 0.5, -1.0, 1.0); x = x * sg[:, None, None]
        x = x / x.abs().amax(dim=(1, 2), keepdim=True).clamp_min(1e-12)
        return x, Yg[ii].float() / 255.0, Wg[ii].float() / 255.0

    @torch.no_grad()
    def validate():
        model.eval(); pos, neg, vl = [], [], 0.0
        for i0 in range(0, len(vidx), 1024):
            ii = torch.from_numpy(vidx[i0:i0 + 1024]).to(dev); x, y, w = batch_of(ii, False)
            with torch.amp.autocast("cuda", dtype=torch.bfloat16):
                p = model(x).float()
            vl += float(masked_bce(p, y, w)) * len(ii)
            mp = torch.nn.functional.max_pool1d(p[:, None], 21, 1, 10)[:, 0]
            k = torch.from_numpy(kind[vidx[i0:i0 + 1024]]).to(dev); ip = k == 0
            if ip.any(): pos.append(mp[ip].gather(1, y[ip].argmax(1)[:, None])[:, 0].cpu().numpy())
            if (~ip).any(): neg.append(mp[~ip][:, 600].cpu().numpy())
        model.train()
        pos = np.concatenate(pos) if pos else np.zeros(0); neg = np.concatenate(neg) if neg else np.zeros(0)
        return vl / len(vidx), auc(pos, neg), float(np.median(pos)) if len(pos) else np.nan, float(np.quantile(neg, 0.99)) if len(neg) else np.nan

    os.makedirs(f"{PD}/models", exist_ok=True); out = f"{PD}/models/ft1.pt"
    vl, va, pm, n99 = validate(); best = va; hist = [dict(step=0, val_loss=vl, val_auc=va, pos_med=pm, neg_p99=n99)]
    C.log(f"[ft] step 0 (SVI weights) val loss {vl:.4f} auc {va:.4f} pos_med {pm:.3f} neg_p99 {n99:.3f}")
    torch.save(dict(arch="pnet", fac=2.0, win_s=60.0, state=model.state_dict(), step=0, val_auc=va, held_out=C.PNET_HO,
                    train_dates=[str(d) for d in dates]), out)
    t0 = time.time(); run = 0.0
    for step in range(1, a.steps + 1):
        ii = np.concatenate([rng.choice(idx[k], nb[j]) for j, k in enumerate((0, 1))])
        x, y, w = batch_of(torch.from_numpy(ii).to(dev), True)
        with torch.amp.autocast("cuda", dtype=torch.bfloat16):
            p = model(x)
        loss = masked_bce(p.float(), y, w); opt.zero_grad(set_to_none=True); loss.backward(); opt.step()
        run += float(loss.detach())
        if step % a.val_every == 0 or step == a.steps:
            vl, va, pm, n99 = validate()
            hist.append(dict(step=step, val_loss=vl, val_auc=va, pos_med=pm, neg_p99=n99, train_loss=run / a.val_every)); run = 0.0
            C.log(f"[ft] step {step} val loss {vl:.4f} auc {va:.4f} pos_med {pm:.3f} neg_p99 {n99:.3f} {time.time() - t0:.0f}s")
            if va > best:
                best = va
                torch.save(dict(arch="pnet", fac=2.0, win_s=60.0, state=model.state_dict(), step=step, val_auc=va, held_out=C.PNET_HO,
                                train_dates=[str(d) for d in dates]), out)
            pd.DataFrame(hist).to_csv(f"{PD}/models/train_ft1.csv", index=False)
    z = torch.load(out, weights_only=False)
    C.log(f"[ft] done {time.time() - t0:.0f}s; kept step {z['step']} (internal val AUC {z['val_auc']:.4f})")


# ------------------------------------------------------------------------------------------------ features
def feat(a):
    import torch
    from concurrent.futures import ThreadPoolExecutor
    sys.path.insert(0, PICKER); import p_scan
    models = [("w", PNET_W05)]
    if os.path.exists(f"{PD}/models/ft1.pt") and not a.no_ft:
        models.append(("f", f"{PD}/models/ft1.pt"))
    MOD = {s: p_scan.load_model(p) for s, p in models}
    C.log(f"[pnet-feat {a.shard}/{a.nshard}] models {[(s, p) for s, p in models]}")
    D = pd.read_csv(f"{C.OUT}/day_table.csv", dtype={"date": str})
    D = D[D.run.isin(a.runs.split(","))]
    if a.days: D = D[D.date.isin(a.days.split(","))]
    D = D.reset_index(drop=True); D = D[D.index % a.nshard == a.shard]
    tag = "".join(s for s, _ in models)
    todo = [r for r in D.itertuples() if os.path.exists(f"{C.OUT}/feat/{r.run}/{r.date}.json")
            and not os.path.exists(f"{PD}/feat_{tag}/{r.run}/{r.date}.parquet")]
    SITE = C.sites(); TPLS = {}
    status = []; t00 = time.time()
    ex = ThreadPoolExecutor(4)
    for n, r in enumerate(todo):
        try:
            meta = json.load(open(f"{C.OUT}/feat/{r.run}/{r.date}.json"))
            F = pd.read_parquet(f"{C.OUT}/feat/{r.run}/{r.date}.parquet", columns=["series", "k", "fam", "i"]) \
                if os.path.getsize(f"{C.OUT}/feat/{r.run}/{r.date}.parquet") > 0 else pd.DataFrame()
            if r.tpl not in TPLS: TPLS[r.tpl] = C.load_tpl(r.tpl)[0]
            M = TPLS[r.tpl]
            stns = list(meta["stations"])
            futs = {ns: ex.submit(load_zne, ns, r.date) for ns in stns}
            tref = pd.Timestamp(r.date).timestamp(); L = meta["L"]; pre_n = meta["pre_n"]; PRE = pre_n / C.FS
            CF = {}
            for ns in stns:
                z, st = futs[ns].result()
                status.append((r.run, r.date, ns, st, ""))
                if z is None: continue
                for s, (m, mt) in MOD.items():
                    CF[(s, ns)] = maxfilt(p_scan.cf_day(m, mt, z["x"]), 10)
                del z
            R = C.rand_origins(r.date, meta["NS"])
            out = F[["series", "k", "fam"]].copy() if len(F) else pd.DataFrame(columns=["series", "k", "fam"])
            cols = {}
            for f, gi in (F.groupby("fam").groups.items() if len(F) else []):
                g = F.loc[gi]
                sub = M[(M.fam == f) & M.netsta.isin(stns)]
                sts = [ns for ns in sub.netsta if (models[0][0], ns) in CF]
                if not sts: continue
                sites = sorted({SITE[s] for s in sts}); qi = {q: k for k, q in enumerate(sites)}
                for s, _ in models:
                    V = np.full((len(g), len(sites)), -1.0); VR = np.full((len(R), len(sites)), -1.0)
                    for rr in sub.itertuples():
                        ns = rr.netsta
                        if (s, ns) not in CF: continue
                        info = meta["stations"][ns]; nw = info["N"] - L + 1
                        sh = int(round((float(rr.tt_s) + float(rr.best_tau) - PRE - (info["t0"] - tref)) * C.FS))
                        shift = int(meta["shifts"][str(SITE[ns])])
                        cf = CF[(s, ns)]
                        jd = np.full(len(g), -1, np.int64)
                        for ser, sft in (("data", 0), ("c1", shift)):
                            mm = (g.series == ser).to_numpy()
                            if mm.any(): jd[mm] = C.s_index(g.i.to_numpy()[mm], sh, sft, nw, pre_n)
                        v = np.where((jd >= 0) & (jd < len(cf)), cf[np.clip(jd, 0, len(cf) - 1)], np.nan)
                        jr = C.s_index(R, sh, 0, nw, pre_n)
                        vr = np.where((jr >= 0) & (jr < len(cf)), cf[np.clip(jr, 0, len(cf) - 1)], np.nan)
                        cols.setdefault(f"pn{s}_{ns}", pd.Series(np.nan, index=F.index, dtype=np.float32)).loc[gi] = v.astype(np.float32)
                        k = qi[SITE[ns]]
                        V[:, k] = np.fmax(V[:, k], np.nan_to_num(v, nan=-1.0)); VR[:, k] = np.fmax(VR[:, k], np.nan_to_num(vr, nan=-1.0))
                    V[V < 0] = np.nan; VR[VR < 0] = np.nan
                    Vs = np.sort(np.nan_to_num(V, nan=-1.0), 1)[:, ::-1]
                    with np.errstate(invalid="ignore", all="ignore"):
                        med = np.nanmedian(VR, 0) if np.isfinite(VR).any() else np.full(len(sites), np.nan)
                        rel = V - med[None]
                        vals = dict(pnet_max=np.where(Vs[:, 0] >= 0, Vs[:, 0], np.nan),
                                    pnet_2nd=np.where(Vs[:, 1] >= 0, Vs[:, 1], np.nan) if len(sites) > 1 else np.full(len(g), np.nan),
                                    n_sites_pnet03=(np.nan_to_num(V, nan=0) > 0.3).sum(1).astype(float),
                                    n_sites_pnet05=(np.nan_to_num(V, nan=0) > 0.5).sum(1).astype(float),
                                    pnet_mean=np.nanmean(V, 1), pnet_rel_max=np.nanmax(np.where(np.isfinite(rel), rel, -9), 1))
                    vals["pnet_rel_max"] = np.where(vals["pnet_rel_max"] <= -8, np.nan, vals["pnet_rel_max"])
                    for c, v in vals.items():
                        cols.setdefault(f"{c}_{s}", pd.Series(np.nan, index=F.index, dtype=float)).loc[gi] = v
            if cols:
                out = pd.concat([out, pd.DataFrame(cols)], axis=1)
            od = f"{PD}/feat_{tag}/{r.run}"; os.makedirs(od, exist_ok=True)
            out.to_parquet(f"{od}/{r.date}.parquet.tmp", index=False); os.replace(f"{od}/{r.date}.parquet.tmp", f"{od}/{r.date}.parquet")
            del CF
            torch.cuda.empty_cache()
            if n % 10 == 0 or n < 3:
                C.log(f"  [{a.shard}] {r.run} {r.date}: {len(out)} rows, pnet_max_w defined "
                      f"{out['pnet_max_w'].notna().mean() if 'pnet_max_w' in out else 0:.2f}; {n + 1}/{len(todo)} {time.time() - t00:.0f}s")
        except Exception as e:
            C.log(f"  ERROR {r.run} {r.date}: {e}\n{traceback.format_exc()[-800:]}")
            status += [(r.run, r.date, ns, "error", f"day job: {str(e)[:100]}") for ns in C.STAS]
        if n % 20 == 0: C.write_status("pnet", str(a.shard), status)
    C.write_status("pnet", str(a.shard), status)
    C.log(f"[pnet-feat {a.shard}] done {len(todo)} run-days {time.time() - t00:.0f}s")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what")
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--cap-pos", type=int, default=60)
    ap.add_argument("--n-neg", type=int, default=40); ap.add_argument("--days", default="")
    ap.add_argument("--steps", type=int, default=6000); ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--lr", type=float, default=5e-5); ap.add_argument("--val-every", type=int, default=500)
    ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--runs", default="v3b,fwd"); ap.add_argument("--no-ft", action="store_true")
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    dict(labels=labels, finetune=finetune, feat=feat)[a.what](a)
