#!/usr/bin/env python
"""DUC_FIN PNET: S probability at every candidate's expected S time (off-the-shelf SVI model + a fine-tune on HIS events).

PNET = lfe_picker P1b PhaseNet-style U-Net trained on Bostock SVI LFEs (pnet_w05.pt; 60-s windows, 20 Hz, 1-8 Hz
Butterworth-4 zero-phase, Z/N/E, per-window 3-C peak normalisation; p_scan.cf_day).  Input here = cache_duc (his chain:
response removed, 1.5-9 Hz, 20 Hz) + a 1-8 Hz Butterworth-4 zero-phase pass (so 1.5-8 Hz), horizontals rotated to N/E
with the epoch azimuths of results/duc_fix/station_epochs.csv, vertical-only stations (his SHZ / EHZ, 1-component)
entered as (Z, 0, 0).  Expected S = detection time + duc_fin_common.arrival_offsets offset of that family-station
(template peak; CVM1.7 moveout for a station without his template); scrambled-control rows at the site's rolled
(original-time) position.  Station value = max CF within +-0.5 s; site value = max over the site's stations.
Features (suffix _w off-the-shelf, _f fine-tuned): pnet_max, pnet_2nd, n_sites_pnet03, n_sites_pnet05, pnet_mean,
pnet_rel_max, per-station pn<sfx>_<NET.STA>.

Fine-tune (the WA / SH_CONF recipe): start = SVI weights; windows ONLY on the 'pnet_train' days (duc_fin_common.day_table:
every 10th 2008-2010 day outside the held-out months; never fitted or evaluated by the confidence model) and ONLY at the
stations not in PNET_HO (PB.B935, NC.KHMB, NC.LBP = held out).  Positives = his catalogue events that pass his rule (cc x
nch >= threshold) at the family's template stations, Gaussian sigma 0.5 s at every labelled S, loss weight within +-3 s;
negatives = random windows of the same station-days, weight 1 except +-3 s around the S of ANY of his detections (all 66
families, raw).  Masked BCE, 1/2 positive + 1/2 negative batches, horizontals zeroed with p = 0.3, polarity flip; internal
validation = every 10th training date; best internal-AUC checkpoint kept (step 0 = SVI weights is a candidate).

  python scripts/duc_fin_pnet.py labels --workers 8
  python scripts/duc_fin_pnet.py finetune --steps 6000
  python scripts/duc_fin_pnet.py feat --shard 0 --nshard 3
"""
import os, sys, json, time, argparse, warnings, traceback
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duc_fin_common as F
import duc_fix_common as DC
PICKER = "/home/jovyan/tremorferometry/lfe_picker/scripts"
PNET_W05 = "/home/jovyan/tremorferometry/lfe_picker/results/p1/models/pnet_w05.pt"
PD = f"{F.OUT}/pnet"
FS = DC.FS; NDAY = int(round(86400 * FS))
NW = 1200; SIG = 10.0; MASK = 60
_G = {}
_SOS = None


def sos18():
    global _SOS
    if _SOS is None:
        from scipy.signal import butter
        _SOS = butter(4, [1.0 / (FS / 2), 8.0 / (FS / 2)], btype="band", output="sos")
    return _SOS


def load_zne(ns, date, look):
    """cache_duc (C, N) -> (3, N) float32 Z, N, E (1-8 Hz pass); (Z, 0, 0) for 1-component stations; status."""
    from scipy.signal import sosfiltfilt
    net, sta = ns.split(".")
    ep = look(sta, pd.Timestamp(date))
    if ep is None: return None, "nodata"
    p = DC.cache_path(net, sta, date)
    if not os.path.exists(p):
        if "plan" not in _G: _G["plan"] = F.plan_status()
        return None, ("missing" if _G["plan"].get((ns, str(pd.Timestamp(date).date()))) in ("built", "cached") else "nodata")
    v = DC.load_cache(net, sta, date)
    if v is None: return None, "error"
    t0, x, val, ch = v
    if x.ndim != 2 or x.shape[1] < 0.5 * 86400 * FS or not np.any(x): return None, "short"
    x = np.where(val, x, 0.0)
    try:
        x = sosfiltfilt(sos18(), x, axis=-1)
    except Exception:
        return None, "error"
    x = np.where(val, x, 0.0).astype(np.float32)
    zi = [k for k, c in enumerate(ch) if c.endswith("Z")]; hi = [k for k, c in enumerate(ch) if not c.endswith("Z")]
    y = np.zeros((3, x.shape[1]), np.float32)
    if zi: y[0] = x[zi[0]]
    nh = 0
    if len(hi) == 2 and isinstance(ep.az, str) and "/" in ep.az and np.any(x[hi[0]]) and np.any(x[hi[1]]):
        a1, a2 = [float(q) for q in ep.az.split("/")]
        r1, r2 = np.radians(a1), np.radians(a2)
        A = np.array([[np.cos(r1), np.sin(r1)], [np.cos(r2), np.sin(r2)]])
        y[1:3] = (np.linalg.inv(A) @ x[hi]).astype(np.float32); nh = 2
    return dict(x=y, t0=t0, n_horiz=nh, chans=ch), "ok"


def maxfilt(v, half=10):
    from scipy.ndimage import maximum_filter1d
    return maximum_filter1d(v, size=2 * half + 1, mode="constant", cval=0.0)


def offsets():
    O = pd.read_csv(f"{F.OUT}/arrival_offsets.csv", dtype={"fam": str})
    return O


# ------------------------------------------------------------------------------------------------ labels
def _init_labels():
    O = offsets(); O = O[O.offset_source == "template"]
    H = F.his_catalogue(); H["date"] = H.t.dt.strftime("%Y-%m-%d")
    _G.update(off={(r.fam, r.netsta): r.offset_s for r in O.itertuples()}, TR={d: g for d, g in H.groupby("date")},
              look=DC.epoch_lookup())


def label_job(args):
    ns, ds, cap_pos, n_neg, seed = args
    z, st = load_zne(ns, ds, _G["look"])
    if z is None: return ns, ds, st, None
    x = z["x"]; t0 = z["t0"]; rng = np.random.default_rng(seed)
    tg = _G["TR"].get(ds)
    lab_t, msk_t = [], []
    if tg is not None:
        for f, g in tg.groupby("fam"):
            o = _G["off"].get((f, ns))
            if o is None: continue
            s = (g.t_epoch.to_numpy() + o - t0) * FS
            msk_t.append(s)
            lab_t.append(s[g.hp.to_numpy()])
    lab_t = np.sort(np.concatenate(lab_t)) if lab_t else np.zeros(0)
    msk_t = np.sort(np.concatenate(msk_t)) if msk_t else np.zeros(0)
    N = x.shape[1]; valid = np.abs(x).sum(0) > 0; cv = np.concatenate([[0], np.cumsum(~valid)])
    X, Y, W, K = [], [], [], []; tt = np.arange(NW)
    def add(s0, kind):
        s0 = int(s0)
        if s0 < 0 or s0 + NW > N or cv[s0 + NW] - cv[s0] > NW // 10: return False
        seg = x[:, s0:s0 + NW]; mx = np.abs(seg).max()
        if not mx > 0: return False
        y = np.zeros(NW, np.float32); w = np.zeros(NW, np.float32)
        if kind == 0:
            lo, hi = np.searchsorted(lab_t, s0 - MASK), np.searchsorted(lab_t, s0 + NW + MASK)
            for tl in lab_t[lo:hi]:
                rr = tl - s0
                y = np.maximum(y, np.exp(-0.5 * ((tt - rr) / SIG) ** 2)); w[np.abs(tt - rr) <= MASK] = 1.0
        else:
            w[:] = 1.0
            ml, mh = np.searchsorted(msk_t, s0 - MASK), np.searchsorted(msk_t, s0 + NW + MASK)
            for tm in msk_t[ml:mh]: w[np.abs(tt - (tm - s0)) <= MASK] = 0.0
        w[~valid[s0:s0 + NW]] = 0.0
        X.append((seg / mx).astype(np.float16)); Y.append((y * 255).round().astype(np.uint8))
        W.append((w * 255).round().astype(np.uint8)); K.append(kind)
        return True
    if len(lab_t):
        npos = 0
        for i in rng.permutation(len(lab_t))[:4 * cap_pos]:
            if npos >= cap_pos: break
            if add(lab_t[i] - rng.uniform(5, 55) * FS, 0): npos += 1
    k = 0
    for s0 in rng.integers(0, N - NW, 6 * n_neg):
        if k >= n_neg: break
        if add(s0, 1): k += 1
    if not X: return ns, ds, "ok", None
    return ns, ds, "ok", (np.stack(X), np.stack(Y), np.stack(W), np.array(K, np.int8), z["n_horiz"])


def labels(a):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    DT = F.day_table(); bdays = list(DT[DT.role == "pnet_train"].date)
    if a.days: bdays = [d for d in bdays if d in set(a.days.split(","))] or a.days.split(",")
    stas = [s for s in F.stations() if s not in F.PNET_HO]
    jobs = [(ns, d, a.cap_pos, a.n_neg, 7 + i * 31 + j) for i, d in enumerate(bdays) for j, ns in enumerate(stas)]
    F.log(f"[labels] {len(jobs)} station-days on {len(bdays)} pnet_train days, {len(stas)} stations (held out {F.PNET_HO})")
    X, Y, Wt, meta, status = [], [], [], [], []; t0 = time.time()
    with ProcessPoolExecutor(a.workers, initializer=_init_labels) as ex:
        fut = [ex.submit(label_job, j) for j in jobs]
        for i, fu in enumerate(as_completed(fut)):
            ns, ds, st, r = fu.result(); status.append((ds, ns, st, ""))
            if r is None: continue
            x, y, w, k, nh = r
            X.append(x); Y.append(y); Wt.append(w); meta.append(pd.DataFrame(dict(netsta=ns, date=ds, kind=k, n_horiz=nh)))
            if (i + 1) % 1000 == 0: F.log(f"[labels] {i + 1}/{len(jobs)} {time.time() - t0:.0f}s")
    F.write_status("pnetlabels", "0", status)
    od = f"{PD}/train"; os.makedirs(od, exist_ok=True)
    Mt = pd.concat(meta, ignore_index=True)
    np.save(f"{od}/X.npy", np.concatenate(X)); np.save(f"{od}/Y.npy", np.concatenate(Y)); np.save(f"{od}/Wt.npy", np.concatenate(Wt))
    Mt.to_parquet(f"{od}/meta.parquet", index=False)
    F.log(f"[labels] {len(Mt)} windows {Mt.groupby(['kind', 'n_horiz']).size().to_dict()}; stations {Mt.netsta.nunique()}; "
          f"dates {Mt.date.nunique()}; {time.time() - t0:.0f}s")
    F.guard("pnetlabels")


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
    assert not Mt.netsta.isin(F.PNET_HO).any(), "held-out station in the training windows"
    dates = sorted(Mt.date.unique()); vd = set(dates[5::10])
    val = Mt.date.isin(vd).to_numpy() if len(dates) >= 10 else (rng.random(len(Mt)) < 0.1)
    kind = Mt.kind.to_numpy(); nh = Mt.n_horiz.to_numpy()
    idx = {k: np.nonzero((kind == k) & ~val)[0] for k in (0, 1)}
    vidx = np.nonzero(val)[0]
    if len(vidx) > 30000: vidx = np.sort(rng.choice(vidx, 30000, replace=False))
    F.log(f"[ft] windows {len(Mt)}; train pos/neg {[len(idx[k]) for k in (0, 1)]}; val {len(vidx)} ({len(vd)} dates)")
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
        return vl / max(len(vidx), 1), auc(pos, neg), float(np.median(pos)) if len(pos) else np.nan, float(np.quantile(neg, 0.99)) if len(neg) else np.nan

    os.makedirs(f"{PD}/models", exist_ok=True); out = f"{PD}/models/ft1.pt"
    vl, va, pm, n99 = validate(); best = va
    hist = [dict(step=0, val_loss=vl, val_auc=va, pos_med=pm, neg_p99=n99, train_loss=np.nan)]
    F.log(f"[ft] step 0 (SVI weights) val loss {vl:.4f} auc {va:.4f} pos_med {pm:.3f} neg_p99 {n99:.3f}")
    torch.save(dict(arch="pnet", fac=2.0, win_s=60.0, state=model.state_dict(), step=0, val_auc=va, held_out=F.PNET_HO,
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
            nstep = a.val_every if step % a.val_every == 0 else step % a.val_every
            hist.append(dict(step=step, val_loss=vl, val_auc=va, pos_med=pm, neg_p99=n99, train_loss=run / nstep)); run = 0.0
            F.log(f"[ft] step {step} train loss {hist[-1]['train_loss']:.4f} val loss {vl:.4f} auc {va:.4f} pos_med {pm:.3f} neg_p99 {n99:.3f} {time.time() - t0:.0f}s")
            if va > best:
                best = va
                torch.save(dict(arch="pnet", fac=2.0, win_s=60.0, state=model.state_dict(), step=step, val_auc=va, held_out=F.PNET_HO,
                                train_dates=[str(d) for d in dates]), out)
            pd.DataFrame(hist).to_csv(f"{PD}/models/train_ft1.csv", index=False)
    pd.DataFrame(hist).to_csv(f"{PD}/models/train_ft1.csv", index=False)
    z = torch.load(out, weights_only=False)
    F.log(f"[ft] done {time.time() - t0:.0f}s; kept step {z['step']} (internal val AUC {z['val_auc']:.4f})")
    learning_curve()


def learning_curve():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    p = f"{PD}/models/train_ft1.csv"
    if not os.path.exists(p): return
    Hc = pd.read_csv(p); z = None
    try:
        import torch; z = torch.load(f"{PD}/models/ft1.pt", weights_only=False, map_location="cpu")
    except Exception: pass
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot(Hc.step, Hc.train_loss, "o-", label="train loss (mean over the last 500 steps)")
    ax[0].plot(Hc.step, Hc.val_loss, "s-", label="validation loss")
    ax[0].set_xlabel("step"); ax[0].set_ylabel("masked BCE"); ax[0].legend(fontsize=8)
    ax[1].plot(Hc.step, Hc.val_auc, "o-", color="C2")
    if z is not None: ax[1].axvline(z["step"], color="k", ls=":", label=f"kept step {z['step']}"); ax[1].legend(fontsize=8)
    ax[1].set_xlabel("step"); ax[1].set_ylabel("validation AUC (pos at S vs random)")
    fig.suptitle("PNET fine-tune on Ducellier's catalogue events (pnet_train days, 3 stations held out)", fontsize=10)
    fig.tight_layout(); os.makedirs(f"{F.OUT}/figures", exist_ok=True)
    fig.savefig(f"{F.OUT}/figures/pnet_learning_curve.png", dpi=130); plt.close(fig)


# ------------------------------------------------------------------------------------------------ features
def feat(a):
    import torch
    from concurrent.futures import ThreadPoolExecutor
    sys.path.insert(0, PICKER); import p_scan
    models = [("w", PNET_W05)]
    if os.path.exists(f"{PD}/models/ft1.pt") and not a.no_ft: models.append(("f", f"{PD}/models/ft1.pt"))
    MOD = {s: p_scan.load_model(p) for s, p in models}
    tag = "".join(s for s, _ in models)
    F.log(f"[pnet-feat {a.shard}/{a.nshard}] models {models}")
    days = a.days.split(",") if a.days else list(F.day_table().date)
    days = [d for j, d in enumerate(days) if j % a.nshard == a.shard]
    todo = [d for d in days if os.path.exists(f"{F.OUT}/feat/{d}.json") and not os.path.exists(f"{PD}/feat_{tag}/{d}.parquet")]
    O = offsets(); OFF = {(r.fam, r.netsta): r.offset_s for r in O.itertuples()}
    SITE = F.sites(); look = DC.epoch_lookup(); famst = {r.family: r.stations.split(",") for r in DC.families().itertuples()}
    no = F.netof()
    status = []; t00 = time.time(); ex = ThreadPoolExecutor(4)
    for n, ds in enumerate(todo):
        try:
            meta = json.load(open(f"{F.OUT}/feat/{ds}.json"))
            pq = f"{F.OUT}/feat/{ds}.parquet"
            Fr = pd.read_parquet(pq, columns=["series", "k", "fam", "g", "t_epoch"]) if os.path.getsize(pq) > 0 else pd.DataFrame()
            stns = list(meta["stations"])
            futs = {ns: ex.submit(load_zne, ns, ds, look) for ns in stns}
            tref = pd.Timestamp(ds).timestamp(); CF = {}
            for ns in stns:
                z, st = futs[ns].result(); status.append((ds, ns, st, ""))
                if z is None: continue
                if abs(z["t0"] - tref) > 1e-3: status[-1] = (ds, ns, "error", "cache t0 != midnight"); continue
                for s, (m, mt) in MOD.items(): CF[(s, ns)] = maxfilt(p_scan.cf_day(m, mt, z["x"]), 10)
                del z
            R = F.rand_origins(ds, NDAY)
            out = Fr[["series", "k"]].copy() if len(Fr) else pd.DataFrame(columns=["series", "k"])
            cols = {}
            for f, gi in (Fr.groupby("fam").groups.items() if len(Fr) else []):
                g = Fr.loc[gi]; f = str(f)
                sts = [f"{no[s]}.{s}" for s in famst.get(f, []) if (models[0][0], f"{no[s]}.{s}") in CF]
                if not sts: continue
                sites = sorted({SITE[s] for s in sts}); qi = {q: k for k, q in enumerate(sites)}
                mc = (g.series == "c1").to_numpy()
                for s, _ in models:
                    V = np.full((len(g), len(sites)), -1.0); VR = np.full((len(R), len(sites)), -1.0)
                    for ns in sts:
                        cf = CF[(s, ns)]; off = OFF.get((f, ns))
                        if off is None: continue
                        tdet = g.t_epoch.to_numpy().copy()
                        if mc.any():
                            sh = int(meta["shifts"][str(SITE[ns])])
                            tdet[mc] = tref + ((g.g.to_numpy()[mc] % NDAY + sh) % NDAY) / FS
                        j = np.round((tdet + off - tref) * FS).astype(np.int64)
                        v = np.where((j >= 0) & (j < len(cf)), cf[np.clip(j, 0, len(cf) - 1)], np.nan)
                        jr = np.round(R + off * FS).astype(np.int64)
                        vr = np.where((jr >= 0) & (jr < len(cf)), cf[np.clip(jr, 0, len(cf) - 1)], np.nan)
                        cols.setdefault(f"pn{s}_{ns}", pd.Series(np.nan, index=Fr.index, dtype=np.float32)).loc[gi] = v.astype(np.float32)
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
                        cols.setdefault(f"{c}_{s}", pd.Series(np.nan, index=Fr.index, dtype=float)).loc[gi] = v
            if cols: out = pd.concat([out, pd.DataFrame(cols)], axis=1)
            od = f"{PD}/feat_{tag}"; os.makedirs(od, exist_ok=True)
            out.to_parquet(f"{od}/{ds}.parquet.tmp", index=False); os.replace(f"{od}/{ds}.parquet.tmp", f"{od}/{ds}.parquet")
            del CF; torch.cuda.empty_cache()
            if n % 20 == 0 or n < 3:
                F.log(f"  [{a.shard}] {ds}: {len(out)} rows, pnet_max_w defined {out['pnet_max_w'].notna().mean() if 'pnet_max_w' in out else 0:.2f}; "
                      f"{n + 1}/{len(todo)} {time.time() - t00:.0f}s; gpu {torch.cuda.max_memory_allocated() / 1e9:.1f} GB")
        except Exception as e:
            F.log(f"  ERROR {ds}: {e}\n{traceback.format_exc()[-800:]}")
            status += [(ds, ns, "error", f"day job: {str(e)[:100]}") for ns in F.stations()]
        if n % 20 == 0: F.write_status("pnet", str(a.shard), status)
    F.write_status("pnet", str(a.shard), status)
    F.log(f"[pnet-feat {a.shard}] done {len(todo)} days {time.time() - t00:.0f}s")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what")
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--cap-pos", type=int, default=60)
    ap.add_argument("--n-neg", type=int, default=40); ap.add_argument("--days", default="")
    ap.add_argument("--steps", type=int, default=6000); ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--lr", type=float, default=5e-5); ap.add_argument("--val-every", type=int, default=500)
    ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--no-ft", action="store_true")
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    dict(labels=labels, finetune=finetune, feat=feat, curve=lambda a: learning_curve())[a.what](a)
