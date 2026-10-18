#!/usr/bin/env python
"""Stage B: fine-tune PNET (lfe_picker P1b pnet_w05, trained on Bostock SVI LFEs) on the WA permanent stations.

Start = the SVI weights.  Data = results/pnet/train (pnet_labels.py; BUILD half of V2008 + V2010, training stations
only).  Loss = lfe_picker's per-window-normalised masked BCE (p1b_train.masked_bce).  Batches: 1/2 positive windows,
1/4 quiet-day windows, 1/4 scrambled-control windows.  Augmentation: horizontals zeroed with p = 0.3 on 3-C windows
(the WA network has many 1-C EHZ stations; the SVI model returns ~0 on (Z, 0, 0)), then re-peak-normalised; random
polarity flip.  Internal validation = every 10th date (still BUILD half; the test half / A2011 / held-out stations are
never used here): window AUC of the CF max within +-0.5 s of the first labelled S (positive windows) vs at the centre
of negative windows; the checkpoint with the best internal AUC is kept (step 0 = the SVI weights is a candidate).

  python scripts/pnet_finetune.py --steps 12000 --tag ft1  -> results/pnet/models/<tag>.pt (+ train_<tag>.csv)
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import pnet_common as PC
sys.path.insert(0, PC.PICKER)


def masked_bce(p, y, w):
    import torch
    p = p.clamp(1e-6, 1 - 1e-6)
    bce = -(y * torch.log(p) + (1 - y) * torch.log(1 - p)) * w
    den = w.sum(-1).clamp_min(1e-6)
    return (bce.sum(-1) / den).mean()


def main():
    import torch
    import p_scan
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=f"{PC.OUT}/train")
    ap.add_argument("--tag", default="ft1")
    ap.add_argument("--init", default=PC.PNET_W05)
    ap.add_argument("--steps", type=int, default=12000)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--lr", type=float, default=5e-5)
    ap.add_argument("--p-zonly", type=float, default=0.3)
    ap.add_argument("--val-every", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    torch.manual_seed(a.seed)
    rng = np.random.default_rng(a.seed)
    dev = "cuda"
    X = np.load(f"{a.data}/X.npy", mmap_mode="r")
    Y = np.load(f"{a.data}/Y.npy"); Wt = np.load(f"{a.data}/Wt.npy")
    M = pd.read_parquet(f"{a.data}/meta.parquet")
    dates = sorted(M.date.unique())
    vdates = set(dates[5::10])
    M["val"] = M.date.isin(vdates) if len(dates) >= 10 else rng.random(len(M)) < 0.1   # tiny smoke sets: random 10 %
    kind = M.kind.to_numpy(); val = M.val.to_numpy(); nh = M.n_horiz.to_numpy()
    idx = {k: np.nonzero((kind == k) & ~val)[0] for k in (0, 1, 2)}
    vidx = np.nonzero(val)[0]
    if len(vidx) > 30000:
        vidx = np.sort(rng.choice(vidx, 30000, replace=False))
    W.log(PC.LOG, f"[ft:{a.tag}] windows {len(M)}; train pos/quiet/ctl {[len(idx[k]) for k in (0, 1, 2)]}; val {len(vidx)} "
          f"({len(vdates)} dates)")
    Xg = torch.from_numpy(np.ascontiguousarray(X)).to(dev)            # float16 on the GPU (few GB)
    Yg = torch.from_numpy(Y).to(dev); Wg = torch.from_numpy(Wt).to(dev)
    NHg = torch.from_numpy(nh).to(dev)
    model, meta = p_scan.load_model(a.init)
    model.train()
    opt = torch.optim.Adam(model.parameters(), lr=a.lr)
    nb = [a.batch // 2, a.batch // 4, a.batch - a.batch // 2 - a.batch // 4]

    def batch_of(ii, aug):
        x = Xg[ii].float()
        if aug:
            zo = (torch.rand(len(ii), device=dev) < a.p_zonly) & (NHg[ii] > 0)
            x[zo, 1:] = 0
            sg = torch.where(torch.rand(len(ii), device=dev) < 0.5, -1.0, 1.0)
            x = x * sg[:, None, None]
        x = x / x.abs().amax(dim=(1, 2), keepdim=True).clamp_min(1e-12)
        return x, Yg[ii].float() / 255.0, Wg[ii].float() / 255.0

    @torch.no_grad()
    def validate():
        model.eval()
        pos, neg, vl = [], [], 0.0
        for i0 in range(0, len(vidx), 1024):
            ii = torch.from_numpy(vidx[i0:i0 + 1024]).to(dev)
            x, y, w = batch_of(ii, False)
            with torch.amp.autocast("cuda", dtype=torch.bfloat16):
                p = model(x).float()
            vl += float(masked_bce(p, y, w)) * len(ii)
            mp = torch.nn.functional.max_pool1d(p[:, None], 21, 1, 10)[:, 0]
            k = torch.from_numpy(kind[vidx[i0:i0 + 1024]]).to(dev)
            ispos = k == 0
            if ispos.any():
                ya = y[ispos].argmax(1)                       # a labelled S (the largest target)
                pos.append(mp[ispos].gather(1, ya[:, None])[:, 0].cpu().numpy())
            if (~ispos).any():
                neg.append(mp[~ispos][:, 600].cpu().numpy())
        model.train()
        pos = np.concatenate(pos) if pos else np.zeros(0); neg = np.concatenate(neg) if neg else np.zeros(0)
        return vl / len(vidx), PC.auc(pos, neg), float(np.median(pos)) if len(pos) else np.nan, \
            float(np.quantile(neg, 0.99)) if len(neg) else np.nan

    os.makedirs(f"{PC.OUT}/models", exist_ok=True)
    out = f"{PC.OUT}/models/{a.tag}.pt"
    hist = []
    vl, vauc, pm, n99 = validate()
    best = vauc
    hist.append(dict(step=0, val_loss=vl, val_auc=vauc, pos_med=pm, neg_p99=n99, train_loss=np.nan))
    W.log(PC.LOG, f"[ft:{a.tag}] step 0 (SVI weights) val loss {vl:.4f} auc {vauc:.4f} pos_med {pm:.3f} neg_p99 {n99:.3f}")
    torch.save(dict(arch="pnet", fac=2.0, win_s=60.0, state=model.state_dict(), cfg=vars(a), step=0, val_auc=vauc), out)
    t0 = time.time(); run = 0.0; tl = np.nan
    for step in range(1, a.steps + 1):
        ii = np.concatenate([rng.choice(idx[k], nb[j]) for j, k in enumerate((0, 1, 2))])
        x, y, w = batch_of(torch.from_numpy(ii).to(dev), True)
        with torch.amp.autocast("cuda", dtype=torch.bfloat16):
            p = model(x)
        loss = masked_bce(p.float(), y, w)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        run += float(loss.detach())
        if step % 200 == 0:
            W.log(PC.LOG, f"[ft:{a.tag}] step {step} loss {run/200:.4f} {step*a.batch/(time.time()-t0):.0f} win/s")
            tl = run / 200; run = 0.0
        if step % a.val_every == 0 or step == a.steps:
            vl, vauc, pm, n99 = validate()
            hist.append(dict(step=step, val_loss=vl, val_auc=vauc, pos_med=pm, neg_p99=n99, train_loss=tl))
            W.log(PC.LOG, f"[ft:{a.tag}] [val] step {step} loss {vl:.4f} auc {vauc:.4f} pos_med {pm:.3f} neg_p99 {n99:.3f}")
            if vauc > best:
                best = vauc
                torch.save(dict(arch="pnet", fac=2.0, win_s=60.0, state=model.state_dict(), cfg=vars(a), step=step,
                                val_auc=vauc), out)
            pd.DataFrame(hist).to_csv(f"{PC.OUT}/models/train_{a.tag}.csv", index=False)
    z = torch.load(out, weights_only=False)
    W.log(PC.LOG, f"[ft:{a.tag}] done {time.time()-t0:.0f}s; kept step {z['step']} (internal val AUC {z['val_auc']:.4f})")


if __name__ == "__main__":
    main()
