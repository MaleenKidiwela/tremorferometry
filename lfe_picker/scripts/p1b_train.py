#!/usr/bin/env python
"""P1b — train one of the two baselines on the P0 labels.

  python scripts/p1b_train.py --arch pnet --steps 60000 --batch 256 --tag pnet_a
  python scripts/p1b_train.py --arch lin  --steps 60000 --batch 256 --tag lin_a

Training scope (frozen by the referee protocol): years 2008/2009/2010/2012, stations S15 minus
B926/NLLB/LRIV, families minus the 30 hidden.  2011 is never touched here; 2013 is forbidden.
Loss = per-window-normalised weighted BCE, so a positive window (supervised only within +-3 s of a
labelled S) counts as much as a fully supervised negative window.
"""
from __future__ import annotations
import os, sys, json, time, argparse
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P
import p_dataset as D
import p_models as M

OUT = f"{P.RES}/p1"


def masked_bce(p, y, w):
    p = p.clamp(1e-6, 1 - 1e-6)
    bce = -(y * torch.log(p) + (1 - y) * torch.log(1 - p)) * w
    den = w.sum(-1).clamp_min(1e-6)
    return (bce.sum(-1) / den).mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", choices=["lin", "pnet"], required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--steps", type=int, default=60000)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--w-pred", type=float, default=0.5, help="weight of `predicted` S labels")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--fac", type=float, default=2.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--val-every", type=int, default=2000)
    a = ap.parse_args()
    os.makedirs(f"{OUT}/models", exist_ok=True)
    torch.manual_seed(a.seed)

    win_s = 15.0 if a.arch == "lin" else 60.0
    stas = [s for s in P.S15 if s not in P.HELD_OUT_STATIONS]
    dev = "cuda"
    model = (M.LinUNet(a.fac) if a.arch == "lin" else M.PhaseNetS()).to(dev)
    nparam = sum(q.numel() for q in model.parameters())
    opt = torch.optim.Adam(model.parameters(), lr=a.lr)
    scaler = torch.amp.GradScaler("cuda")
    print(f"[cfg] arch={a.arch} params={nparam/1e6:.2f}M win={win_s}s batch={a.batch} "
          f"steps={a.steps} w_pred={a.w_pred}", flush=True)

    tw = D.TorchWindows(P.TRAIN_YEARS, stas, a.batch, win_s=win_s, seed=a.seed + 1,
                        w_pred=a.w_pred, workers=a.workers)
    it = iter(tw)
    # a small fixed validation batch stream (same distribution, different seed/draws)
    vw = D.TorchWindows(P.DEV_YEARS, stas, 512, win_s=win_s, seed=a.seed + 777,
                        w_pred=a.w_pred, workers=2)
    vit = iter(vw)

    hist = []
    t0 = time.time()
    run = 0.0
    for step in range(1, a.steps + 1):
        X, Y, W, _U = next(it)
        X, Y, W = X.to(dev, non_blocking=True), Y.to(dev, non_blocking=True), W.to(dev, non_blocking=True)
        with torch.amp.autocast("cuda", dtype=torch.bfloat16):
            if a.arch == "lin":
                Xu = M.upsample2(X)
                Yu = torch.nn.functional.interpolate(Y[:, None], scale_factor=2, mode="linear",
                                                     align_corners=False)[:, 0]
                Wu = torch.nn.functional.interpolate(W[:, None], scale_factor=2,
                                                     mode="nearest")[:, 0]
                p = model(M.zen2inp(Xu))
                loss = masked_bce(p.float(), Yu, Wu)
            else:
                p = model(X)
                loss = masked_bce(p.float(), Y, W)
        opt.zero_grad(set_to_none=True)
        scaler.scale(loss).backward()
        scaler.step(opt)
        scaler.update()
        run += float(loss.detach())
        if step % 200 == 0:
            print(f"  step {step:6d}  loss {run/200:.4f}  {step*a.batch/(time.time()-t0):.0f} win/s",
                  flush=True)
            run = 0.0
        if step % a.val_every == 0 or step == a.steps:
            model.eval()
            with torch.no_grad(), torch.amp.autocast("cuda", dtype=torch.bfloat16):
                vl, vp, vn = 0.0, [], []
                for _ in range(4):
                    X, Y, W, _U = next(vit)
                    X, Y, W = X.to(dev), Y.to(dev), W.to(dev)
                    if a.arch == "lin":
                        Xu = M.upsample2(X)
                        Yv = torch.nn.functional.interpolate(Y[:, None], scale_factor=2,
                                                             mode="linear", align_corners=False)[:, 0]
                        Wv = torch.nn.functional.interpolate(W[:, None], scale_factor=2,
                                                             mode="nearest")[:, 0]
                        q = model(M.zen2inp(Xu)).float()
                    else:
                        Yv, Wv, q = Y, W, model(X).float()
                    vl += float(masked_bce(q, Yv, Wv))
                    ispos = (Yv.max(-1).values > 0.5)
                    vp.append(q[ispos].max(-1).values.cpu().numpy())
                    vn.append(q[~ispos].max(-1).values.cpu().numpy())
            vp, vn = np.concatenate(vp), np.concatenate(vn)
            # window-level AUC of max-CF, positives vs certain/hard negatives (dev years)
            allv = np.concatenate([vp, vn]); lab = np.r_[np.ones(len(vp)), np.zeros(len(vn))]
            o = np.argsort(allv); r = np.empty(len(allv)); r[o] = np.arange(len(allv))
            auc = float((r[lab == 1].sum() - len(vp) * (len(vp) - 1) / 2) / (len(vp) * len(vn)))
            hist.append(dict(step=step, val_loss=vl / 4, win_auc=auc,
                             pos_med=float(np.median(vp)), neg_p99=float(np.quantile(vn, 0.99))))
            print(f"  [val] step {step} loss {vl/4:.4f} win_auc {auc:.4f} "
                  f"pos_med {np.median(vp):.3f} neg_p99 {np.quantile(vn,0.99):.3f}", flush=True)
            model.train()
            torch.save(dict(arch=a.arch, fac=a.fac, win_s=win_s, state=model.state_dict(),
                            cfg=vars(a), nparam=nparam),
                       f"{OUT}/models/{a.tag}.pt")
            import pandas as pd
            pd.DataFrame(hist).to_csv(f"{OUT}/train_{a.tag}.csv", index=False)
    print(f"[done] {time.time()-t0:.0f}s -> {OUT}/models/{a.tag}.pt")


if __name__ == "__main__":
    main()
