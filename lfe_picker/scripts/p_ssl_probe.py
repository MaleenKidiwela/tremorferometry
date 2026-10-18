#!/usr/bin/env python
"""The GATE on SSL pretraining (Merlin P2): cross-station prediction, then a linear probe.

Pretext: mask one target station's 60-s 3-C window and reconstruct its ENVELOPE from the other
stations' windows of the same absolute time.  Constraints Merlin set:
  * training data 2014-2016 ONLY (never 2010-2013, so no referee leakage);
  * the target station must be >= 15 km from EVERY input station (so intra-site coherence -- the
    PGC and LRIV clusters, his F4 -- cannot be what is learned);
  * catalogued-earthquake windows excluded.
Probe: freeze the encoder, fit a LOGISTIC REGRESSION on its pooled features to separate 2012
Bostock event windows from in-tremor non-event windows of the same days.
  PASS if AUC >= 0.85 -> pretraining is allowed in P2.  Otherwise pretraining is DROPPED and said so.

  python scripts/p_ssl_probe.py --build      # the 2014-2016 day store (subsampled)
  python scripts/p_ssl_probe.py --train --steps 20000
  python scripts/p_ssl_probe.py --probe
"""
from __future__ import annotations
import os, sys, json, time, argparse
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P
import p_models as M

OUT = f"{P.RES}/p2"
TAG = "_ssl"
FS = P.FS
WIN = 1200
MIN_KM = 15.0


def pairs_ok():
    ST = P.station_table()
    d = {}
    for s in P.S15:
        d[s] = [q for q in P.S15 if q != s and
                P.epi_km(ST.lat[s], ST.lon[s], ST.lat[q], ST.lon[q]) >= MIN_KM]
    return d


class Enc(nn.Module):
    """encoder shared by the pretext and the probe: (B,3,1200) -> (B,C,75) features."""

    def __init__(self, ch=(24, 48, 96, 192)):
        super().__init__()
        self.inp = M.Block(3, ch[0])
        self.down = nn.ModuleList([M.Block(ch[i], ch[i + 1]) for i in range(len(ch) - 1)])
        self.out_ch = ch[-1]

    def forward(self, x):
        h = self.inp(x)
        for d in self.down:
            h = d(F.max_pool1d(h, 2))
        return h                                  # (B, 192, 150)


class Pretext(nn.Module):
    def __init__(self, nin=6):
        super().__init__()
        self.enc = Enc()
        self.mix = nn.Conv1d(self.enc.out_ch, self.enc.out_ch, 1)
        self.head = nn.Sequential(nn.Conv1d(self.enc.out_ch, 96, 5, padding=2), nn.ReLU(),
                                  nn.Conv1d(96, 1, 1))

    def forward(self, xs):
        """xs (B, N, 3, 1200) inputs from N stations -> predicted target ENVELOPE (B, 150)."""
        B, N = xs.shape[:2]
        h = self.enc(xs.reshape(B * N, 3, WIN)).reshape(B, N, self.enc.out_ch, -1)
        h = self.mix(h.mean(1))
        return self.head(h).squeeze(1)


def envelope(x, n=8):
    e = (x ** 2).sum(1, keepdim=True).sqrt()
    return F.avg_pool1d(e, n, n).squeeze(1)


def build(a):
    import subprocess
    cmd = [sys.executable, f"{P.HERE}/p_build_mm.py", "--years", "2014,2015,2016",
           "--tag", TAG, "--stride", str(a.stride), "--workers", str(a.workers)]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def sample_batch(store, days_by_doy, ok, rng, B=64, nin=4):
    XS, YT = [], []
    keys = list(days_by_doy)
    while len(XS) < B:
        k = keys[rng.integers(len(keys))]
        avail = days_by_doy[k]
        tgt = avail[rng.integers(len(avail))]
        inp = [s for s in avail if s in ok[tgt]]
        if len(inp) < nin:
            continue
        inp = list(rng.permutation(inp)[:nin])
        i_t = store.idx(tgt, *k)
        n = int(store.nvalid[i_t])
        j = int(rng.integers(0, max(1, n - WIN)))
        t0 = store.t0[i_t] + j / FS
        xs, bad = [], False
        for s in inp:
            i = store.idx(s, *k)
            js = int(round((t0 - store.t0[i]) * FS))
            if js < 0 or js + WIN > int(store.nvalid[i]):
                bad = True
                break
            w = np.asarray(store.get_i(i)[:, js:js + WIN], np.float32)
            m = np.abs(w).max()
            if m <= 0:
                bad = True
                break
            xs.append(w / m)
        if bad:
            continue
        y = np.asarray(store.get_i(i_t)[:, j:j + WIN], np.float32)
        m = np.abs(y).max()
        if m <= 0:
            continue
        XS.append(np.stack(xs))
        YT.append(y / m)
    return np.stack(XS), np.stack(YT)


class RamStore:
    """The SSL memmap lives on NFS and random 7 kB reads from it run at ~1 % CPU.  Pull a fixed
    subset of whole days into RAM once (sequential) and sample from there."""

    def __init__(self, store, keys):
        self.t0, self.nvalid, self.key = {}, {}, {}
        self.X = np.empty((len(keys), 3, P.NDAY), np.float16)
        for n, k in enumerate(keys):
            i = store.idx(*k)
            self.X[n] = store.get_i(i)
            self.t0[k] = store.t0[i]
            self.nvalid[k] = store.nvalid[i]
            self.key[k] = n

    def idx(self, sta, year, doy):
        return (sta, year, doy)

    def get_i(self, k):
        return self.X[self.key[k]]


def train(a):
    store = P.DayStore(tag=TAG)
    ok = pairs_ok()
    by = {}
    for r in store.di.itertuples():
        by.setdefault((r.year, r.doy), []).append(r.sta)
    by = {k: v for k, v in by.items() if len(v) >= 6}
    days = sorted(by)[::max(1, len(by) // a.ndays)][:a.ndays]
    by = {k: by[k] for k in days}
    keys = [(s, y, d) for (y, d) in days for s in by[(y, d)]]
    t_ram = time.time()
    store = RamStore(store, keys)
    print(f"[ssl] {len(keys)} station-days over {len(by)} days pulled into RAM "
          f"({len(keys)*3*P.NDAY*2/1e9:.1f} GB, {time.time()-t_ram:.0f}s)", flush=True)
    dev = "cuda"
    m = Pretext().to(dev)
    opt = torch.optim.Adam(m.parameters(), 3e-4)
    rng = np.random.default_rng(0)
    t0 = time.time()
    run = 0.0
    for step in range(1, a.steps + 1):
        X, Y = sample_batch(store, by, ok, rng, B=a.batch, nin=4)
        X = torch.from_numpy(X).to(dev)
        Y = torch.from_numpy(Y).to(dev)
        with torch.amp.autocast("cuda", dtype=torch.bfloat16):
            p = m(X).float()
            tgt = envelope(Y, WIN // p.shape[-1])
            tgt = tgt / tgt.mean(-1, keepdim=True).clamp_min(1e-9)
            loss = F.smooth_l1_loss(p, torch.log1p(tgt))
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        run += float(loss.detach())
        if step % 200 == 0:
            print(f"  step {step} loss {run/200:.4f} {step*a.batch/(time.time()-t0):.0f}/s",
                  flush=True)
            run = 0.0
    torch.save(dict(state=m.state_dict(), cfg=vars(a)), f"{OUT}/models/ssl_pretext.pt")
    print(f"[done] -> {OUT}/models/ssl_pretext.pt")


@torch.no_grad()
def probe(a):
    """2012 Bostock event windows vs IN-TREMOR non-event windows of the same station-days."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold
    dev = "cuda"
    m = Pretext().to(dev)
    m.load_state_dict(torch.load(f"{OUT}/models/ssl_pretext.pt", map_location=dev,
                                 weights_only=False)["state"])
    m.eval()
    store = P.DayStore()
    L = pd.read_parquet(f"{P.RES}/labels_v0.parquet")
    ps = L[(L.kind == "pos_S") & (L.year == 2012) & (~L.hidden_family)]
    ets = ps.groupby("doy").size()
    ets = set(ets[ets >= 2000].index)
    ps = ps[ps.doy.isin(ets)]
    rng = np.random.default_rng(5)
    Xs, ys, gs = [], [], []
    for (sta, doy), g in ps.groupby(["sta", "doy"]):
        if not store.has(sta, 2012, doy):
            continue
        i = store.idx(sta, 2012, doy)
        n = int(store.nvalid[i])
        tev = np.sort(g.t.values)
        pick = rng.choice(len(tev), min(60, len(tev)), replace=False)
        for k in pick:
            j = int(round((tev[k] - 30.0 - store.t0[i]) * FS))
            if j < 0 or j + WIN > n:
                continue
            Xs.append(np.asarray(store.get_i(i)[:, j:j + WIN], np.float32))
            ys.append(1); gs.append(f"{sta}_{doy}")
        # in-tremor non-event windows of the SAME station-day, >= 60 s from any event
        got = 0
        for _ in range(400):
            if got >= 60:
                break
            j = int(rng.integers(0, max(1, n - WIN)))
            t0 = store.t0[i] + (j + WIN / 2) / FS
            if np.min(np.abs(tev - t0)) < 60.0:
                continue
            Xs.append(np.asarray(store.get_i(i)[:, j:j + WIN], np.float32))
            ys.append(0); gs.append(f"{sta}_{doy}")
            got += 1
    X = np.stack(Xs)
    X = X / np.abs(X).max((1, 2), keepdims=True).clip(1e-9)
    y = np.array(ys)
    g = np.array(gs)
    feats = []
    for i0 in range(0, len(X), 512):
        b = torch.from_numpy(X[i0:i0 + 512]).to(dev)
        with torch.amp.autocast("cuda", dtype=torch.bfloat16):
            h = m.enc(b).float()
        feats.append(torch.cat([h.mean(-1), h.amax(-1)], 1).cpu().numpy())
    Fm = np.concatenate(feats)
    aucs = []
    for tr, te in GroupKFold(n_splits=4).split(Fm, y, g):
        lr = LogisticRegression(max_iter=2000, C=1.0)
        mu, sd = Fm[tr].mean(0), Fm[tr].std(0) + 1e-9
        lr.fit((Fm[tr] - mu) / sd, y[tr])
        aucs.append(roc_auc_score(y[te], lr.predict_proba((Fm[te] - mu) / sd)[:, 1]))
    auc = float(np.mean(aucs))
    res = dict(n=int(len(y)), n_pos=int(y.sum()), auc_mean=auc, auc_folds=[float(x) for x in aucs],
               bar=0.85, pass_=bool(auc >= 0.85),
               decision=("pretraining ALLOWED" if auc >= 0.85 else
                         "pretraining DROPPED (Merlin's gate not met)"))
    json.dump(res, open(f"{OUT}/ssl_probe.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--train", action="store_true")
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--stride", type=int, default=4)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--batch", type=int, default=48)
    ap.add_argument("--ndays", type=int, default=45)
    a = ap.parse_args()
    os.makedirs(f"{OUT}/models", exist_ok=True)
    if a.build:
        build(a)
    if a.train:
        train(a)
    if a.probe:
        probe(a)
