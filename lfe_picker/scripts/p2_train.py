#!/usr/bin/env python
"""P2 — the discovery CF: the better P1b architecture, retrained with a station embedding under
dropout and with REAL-WAVEFORM TRANSPLANTS at random Slab2 cells.

Transplants (Merlin P2, "training extras"): take a LOUD Bostock event's actual 3-C windows on its own
day across the training stations, and add them into a quiet day at the moveout of a RANDOM Slab2
grid cell, amplitude 0.2-1 of the donor (annealed loud -> quiet), at most --transplant-frac of the
samples.  The donor waveform is real; only its position on the interface (hence its inter-station
timing) is new, which is exactly what exercises the CF and the back-projection away from the 130
known patches without teaching a position prior.  The injected fraction is logged.

Station identity enters only as an EMBEDDING that is dropped (replaced by an "unknown station"
token) with probability --sta-drop, so the weights stay station-agnostic; whole-site drops are an
evaluation-time ablation of the back-projection, not of the CF.

  python scripts/p2_train.py --tag p2core --steps 150000 --init results/p1/models/pnet_w05.pt
"""
from __future__ import annotations
import os, sys, json, time, argparse
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P
import p_dataset as D
import p_models as M

OUT = f"{P.RES}/p2"
FS = P.FS


class Transplanter:
    """Plants real loud-event 3-C snippets into quiet-day windows at a random cell's moveout."""

    def __init__(self, store, bank, stas, snip=(-4.0, 8.0), mag_q=0.6, seed=0,
                 tt_npz=f"{OUT}/grid_tt.npz"):
        self.st, self.b, self.stas = store, bank, stas
        self.rng = np.random.default_rng(seed)
        self.pre, self.post = snip
        z = np.load(tt_npz)
        gs = [str(s) for s in z["stas"]]
        self.cell_tt = {s: z["tt"][:, gs.index(s)] for s in stas}
        self.ncell = z["tt"].shape[0]
        self.lat, self.lon = z["lat"], z["lon"]
        B = P.bostock()
        B = B[B.year.isin(bank.years) & (~B.hidden)]
        thr = B.mag.quantile(mag_q)
        self.don = B[B.mag >= thr].reset_index(drop=True)
        self.tt = P.fam_tt_labels(stas)
        self.negdays = [k for k in bank.negdays if store.has(*k)]
        n = int(round((self.post - self.pre) * FS))
        self.nsnip = n
        w = np.ones(n, np.float32)
        e = int(0.1 * n)
        w[:e] = 0.5 * (1 - np.cos(np.pi * np.arange(e) / e))
        w[-e:] = w[:e][::-1]
        self.taper = w

    def _snip(self, sta, year, doy, t_s):
        if not self.st.has(sta, year, doy):
            return None
        i = self.st.idx(sta, year, doy)
        j = self.st.sample_of(i, t_s + self.pre)
        if j < 0 or j + self.nsnip > int(self.st.nvalid[i]):
            return None
        x = np.asarray(self.st.get_i(i)[:, j:j + self.nsnip], np.float32)
        if np.abs(x).max() <= 0:
            return None
        return x * self.taper

    def draw(self, amp, win_s, nmax=6):
        """-> list of (x[3,nw], y[nw], w[nw], sta) for one planted event, or []."""
        nw = int(round(win_s * FS))
        for _ in range(12):
            e = self.don.iloc[self.rng.integers(len(self.don))]
            cell = int(self.rng.integers(self.ncell))
            k = self.negdays[self.rng.integers(len(self.negdays))]
            _, qy, qd = k
            i_q = self.st.idx(*k)
            nq = int(self.st.nvalid[i_q])
            t_q = self.st.t0[i_q] + self.rng.uniform(600, max(601, nq / FS - 600))
            g = self.b.guard.get((qy, qd), np.zeros(0))
            if g.size and ((g >= t_q - 200) & (g <= t_q + 200)).any():
                continue
            order = self.rng.permutation(self.stas)[:nmax]
            out = []
            for sta in order:
                if not self.st.has(sta, qy, qd):
                    continue
                tp = self.tt.get((int(e.bkey), sta))
                if tp is None:
                    continue
                sn = self._snip(sta, int(e.year), int(e.doy), float(e.t_origin) + tp[0])
                if sn is None:
                    continue
                t_arr = t_q + float(self.cell_tt[sta][cell])
                i_r = self.st.idx(sta, qy, qd)
                off = self.rng.uniform(4.0, win_s - 4.0)
                i0 = int(round((t_arr - off - self.st.t0[i_r]) * FS))
                if i0 < 0 or i0 + nw > int(self.st.nvalid[i_r]):
                    continue
                x = np.asarray(self.st.get_i(i_r)[:, i0:i0 + nw], np.float32)
                js = int(round((t_arr + self.pre - self.st.t0[i_r]) * FS)) - i0
                if js < 0 or js + self.nsnip > nw:
                    continue
                x = x.copy()
                x[:, js:js + self.nsnip] += amp * sn
                m = np.abs(x).max()
                if m <= 0:
                    continue
                tt_ = np.arange(nw) / FS
                y = np.exp(-0.5 * ((tt_ - off) / D.SIG) ** 2).astype(np.float32)
                w = (np.abs(tt_ - off) <= D.SUPP).astype(np.float32)
                out.append((x / m, y, w, sta))
            if len(out) >= 2:
                return out, cell, int(e.bkey), amp
        return [], -1, -1, amp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="p2core")
    ap.add_argument("--arch", default="pnet", choices=["pnet", "lin"])
    ap.add_argument("--init", default="", help="warm start from a P1b checkpoint")
    ap.add_argument("--steps", type=int, default=150000)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--w-pred", type=float, default=0.5)
    ap.add_argument("--sta-drop", type=float, default=0.3)
    ap.add_argument("--transplant-frac", type=float, default=0.25)
    ap.add_argument("--nnpu-pi", type=float, default=0.0,
                    help="class prior of the UNLABELLED ETS-day tremor; 0 = drop it (Merlin F6 "
                         "allows either). > 0 adds the non-negative PU risk "
                         "max(0, R_u^- - pi * R_p^-).")
    ap.add_argument("--nnpu-lam", type=float, default=1.0)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--val-every", type=int, default=2500)
    a = ap.parse_args()
    assert a.transplant_frac <= 0.30, "Merlin caps transplants at 30 % of samples"
    os.makedirs(f"{OUT}/models", exist_ok=True)
    torch.manual_seed(a.seed)
    win_s = 15.0 if a.arch == "lin" else 60.0
    stas = [s for s in P.S15 if s not in P.HELD_OUT_STATIONS]
    sid_of = {s: i for i, s in enumerate(P.S15)}      # held-out stations keep an id they never see
    dev = "cuda"
    model = (M.LinUNet() if a.arch == "lin" else M.PhaseNetS(nsta=len(P.S15))).to(dev)
    if a.init:
        z = torch.load(a.init, map_location="cpu", weights_only=False)
        miss = model.load_state_dict(z["state"], strict=False)
        print(f"[init] warm start from {a.init}: missing={list(miss.missing_keys)}")
    opt = torch.optim.Adam(model.parameters(), lr=a.lr)
    scaler = torch.amp.GradScaler("cuda")

    # ---- batches (base windows + transplants) are built in worker processes ----
    from torch.utils.data import DataLoader, IterableDataset

    class Mixed(IterableDataset):
        def __iter__(self):
            info = torch.utils.data.get_worker_info()
            wid = 0 if info is None else info.id
            rng = np.random.default_rng(a.seed * 977 + wid)
            store = P.DayStore()
            bank = D.LabelBank(P.TRAIN_YEARS, stas, w_pred=a.w_pred)
            smp = D.WindowSampler(bank, store, win_s=win_s, rng=rng)
            tp = Transplanter(store, bank, stas, seed=a.seed * 31 + wid)
            ntp = int(round(a.transplant_frac * a.batch))
            step = 0
            while True:
                step += 1
                X, Y, W, U, S = smp.batch(a.batch)
                amp = float(np.clip((1.0 - 0.8 * min(1.0, step / (a.steps / max(1, a.workers))))
                                    * float(np.exp(rng.uniform(-0.4, 0.2))), 0.2, 1.0))
                got, k = [], 0
                while len(got) < ntp and k < 10:
                    o, _c, _b, _a = tp.draw(amp, win_s, nmax=12)
                    got += o
                    k += 1
                for j, (x, y, w, sta) in enumerate(got[:ntp]):
                    X[j], Y[j], W[j], S[j] = x, y, w, sta
                    U[j] = 0.0                    # a transplant window is fully accounted for
                sid = np.array([sid_of[s] for s in S], np.int64)
                sid[rng.random(len(sid)) < a.sta_drop] = len(P.S15)
                yield (torch.from_numpy(X), torch.from_numpy(Y), torch.from_numpy(W),
                       torch.from_numpy(U), torch.from_numpy(sid),
                       torch.tensor(min(ntp, len(got))), torch.tensor(amp))

    dl = DataLoader(Mixed(), batch_size=None, num_workers=a.workers,
                    persistent_workers=True, prefetch_factor=4)
    it = iter(dl)
    print(f"[tp] transplants in {a.workers} workers, cap {a.transplant_frac:.0%} of samples",
          flush=True)

    ntp_win, hist, t0, run, amp = 0, [], time.time(), 0.0, 1.0
    for step in range(1, a.steps + 1):
        X, Y, W, U, sid, ntp_i, amp_t = next(it)
        ntp_win += int(ntp_i)
        amp = float(amp_t)
        X, Y, W, U, sid = X.to(dev), Y.to(dev), W.to(dev), U.to(dev), sid.to(dev)
        with torch.amp.autocast("cuda", dtype=torch.bfloat16):
            p = model(X, sid) if a.arch == "pnet" else model(M.zen2inp(M.upsample2(X)))
            q = p.float().clamp(1e-6, 1 - 1e-6)
            bce = -(Y * torch.log(q) + (1 - Y) * torch.log(1 - q)) * W
            loss = (bce.sum(-1) / W.sum(-1).clamp_min(1e-6)).mean()
            if a.nnpu_pi > 0:
                # non-negative PU risk on the ETS-day tremor that F6 forbids calling negative:
                #   max(0, E_u[-log(1-p)] - pi * E_p[-log(1-p)])
                nl = -torch.log(1 - q)
                pw = W * (Y > 0.5)                       # the labelled-positive support
                Ru = (nl * U).sum() / U.sum().clamp_min(1.0)
                Rp = (nl * pw).sum() / pw.sum().clamp_min(1.0)
                loss = loss + a.nnpu_lam * torch.clamp(Ru - a.nnpu_pi * Rp, min=0.0)
        opt.zero_grad(set_to_none=True)
        scaler.scale(loss).backward()
        scaler.step(opt)
        scaler.update()
        run += float(loss.detach())
        if step % 200 == 0:
            print(f"  step {step:6d} loss {run/200:.4f} amp {amp:.2f} "
                  f"tp_win {ntp_win} ({ntp_win/(step*a.batch):.1%})  "
                  f"{step*a.batch/(time.time()-t0):.0f} win/s", flush=True)
            run = 0.0
        if step % a.val_every == 0 or step == a.steps:
            torch.save(dict(arch=a.arch, fac=2.0, win_s=win_s, nsta=len(P.S15),
                            state=model.state_dict(), cfg=vars(a),
                            transplant_windows=ntp_win, sid_of=sid_of),
                       f"{OUT}/models/{a.tag}.pt")
            hist.append(dict(step=step, loss=float(loss.detach()), amp=amp,
                             tp_frac=ntp_win / (step * a.batch)))
            pd.DataFrame(hist).to_csv(f"{OUT}/train_{a.tag}.csv", index=False)
    print(f"[done] {time.time()-t0:.0f}s; transplanted {ntp_win} windows "
          f"({ntp_win/(a.steps*a.batch):.1%} of samples)")


if __name__ == "__main__":
    main()
