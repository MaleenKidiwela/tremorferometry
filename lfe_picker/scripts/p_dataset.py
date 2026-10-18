#!/usr/bin/env python
"""Window sampler for P1b / P2 — the labelling policy in one place.

POSITIVES     every pos_S row of labels_v0 in the training years, at the training stations, for the
              KEPT families.  Target = Gaussian, sigma 0.5 s, at every labelled S inside the window
              (ETS bursts put several in one 60-s window; all of them are labelled).
              Sample weight 1.0 for `picked` provenance, W_PRED (default 0.5) for `predicted`.
NEGATIVES     (a) certain: windows on PNSN-quiet station-days, with +-120 s around EVERY Bostock
              event of ANY family (hidden included) excluded, target 0 at weight 1 everywhere;
              (b) hard: catalogued EQ/blast windows, HF-veto quiet-day transients, T_ext quiet-day
              detections -- target 0 at weight 1 everywhere.
UNLABELLED    everything else.  Concretely (Merlin F6): inside a POSITIVE window the loss is applied
              ONLY within +-3 s of a labelled S; the rest of the window is tremor background of
              unknown content and carries weight 0.  Windows within +-3 s of a HIDDEN family's event
              carry weight 0 there too, so the holdout leaks neither as a positive nor as a negative.

Everything reads the ZNE memmap day store (p_build_mm.py); windows are peak-normalised over the
three components at sampling time.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

FS = P.FS
SIG = 0.5                 # target Gaussian sigma, s (labels_v0)
SUPP = 3.0                # supervised half-width around a labelled S, s
GUARD = 120.0             # exclusion around any Bostock event on a certain-negative day, s


class LabelBank:
    """Per station-day event lists + the negative day inventory, filtered to a training scope."""

    def __init__(self, years, stas, drop_hidden=True, w_pred=0.5, labels=None):
        self.years, self.stas = set(years), set(stas)
        self.w_pred = w_pred
        L = labels if labels is not None else pd.read_parquet(f"{P.RES}/labels_v0.parquet")
        self.raw = L
        ps = L[(L.kind == "pos_S") & L.year.isin(self.years) & L.sta.isin(self.stas)].copy()
        ps["w"] = np.where(ps.provenance.values == "picked", 1.0, w_pred)
        self.kept = ps[~ps.hidden_family] if drop_hidden else ps
        self.hid = ps[ps.hidden_family]
        # event times per station-day
        self.ev = {}
        for (s, y, d), g in self.kept.groupby(["sta", "year", "doy"], sort=False):
            o = np.argsort(g.t.values)
            self.ev[(s, y, d)] = (g.t.values[o], g.w.values[o], g.fam.values[o])
        self.evh = {k: np.sort(g.t.values)
                    for k, g in self.hid.groupby(["sta", "year", "doy"], sort=False)}
        # certain-negative station-days
        nd = L[(L.kind == "neg_certain_day") & L.year.isin(self.years) & L.sta.isin(self.stas)]
        self.negdays = [(r.sta, r.year, r.doy) for r in nd.itertuples()]
        # all Bostock origin-derived S times per (year, doy) for the guard (any family, any station)
        allpos = L[(L.kind == "pos_S") & L.year.isin(self.years)]
        self.guard = {k: np.sort(g.t.values) for k, g in allpos.groupby(["year", "doy"], sort=False)}
        # hard negatives: station-specific (eq) and network-level (hfveto / text -> all stations)
        hn = L[L.kind.isin(["hard_neg_eq", "hard_neg_hfveto", "hard_neg_text"])
               & L.year.isin(self.years)].copy()
        hn = hn[np.isfinite(hn.t.values)]
        self.hard = hn[["kind", "sta", "year", "doy", "t"]].reset_index(drop=True)
        self.posdays = sorted(self.ev.keys())

    def summary(self):
        return dict(pos_rows=int(len(self.kept)), hidden_rows=int(len(self.hid)),
                    pos_station_days=len(self.ev), neg_station_days=len(self.negdays),
                    hard_rows=int(len(self.hard)),
                    picked=int((self.kept.provenance == "picked").sum()),
                    predicted=int((self.kept.provenance == "predicted").sum()))


class WindowSampler:
    """Draws (x[3,nw], y[nw], w[nw]) windows; x is peak-normalised, y the S target, w the mask."""

    def __init__(self, bank: LabelBank, store: P.DayStore, win_s=60.0, rng=None,
                 pos_frac=0.5, hard_frac=0.15, edge_s=4.0):
        self.b, self.st = bank, store
        self.nw = int(round(win_s * FS))
        self.win_s = win_s
        self.rng = rng if rng is not None else np.random.default_rng(0)
        self.pos_frac, self.hard_frac, self.edge_s = pos_frac, hard_frac, edge_s
        self.posdays = [k for k in bank.posdays if store.has(*k)]
        self.negdays = [k for k in bank.negdays if store.has(*k)]
        h = bank.hard
        self.hard_net = h[h.kind != "hard_neg_eq"][["year", "doy", "t"]].to_numpy()
        he = h[h.kind == "hard_neg_eq"]
        self.hard_sta = he[["sta", "year", "doy", "t"]].to_numpy() if len(he) else np.zeros((0, 4))
        self.allstas = sorted({k[0] for k in self.posdays} | {k[0] for k in self.negdays})

    # -------------------------------------------------------------- helpers
    def _cut(self, sta, year, doy, i0):
        i = self.st.idx(sta, year, doy)
        n = int(self.st.nvalid[i])
        if i0 < 0 or i0 + self.nw > n:
            return None, None
        x = np.asarray(self.st.get_i(i)[:, i0:i0 + self.nw], dtype=np.float32)
        if not np.isfinite(x).all():
            return None, None
        return x, self.st.t0[i] + i0 / FS

    @staticmethod
    def _norm(x):
        m = np.abs(x).max()
        return x / m if m > 0 else x

    def _target(self, t_win, ev_t, ev_w, hid_t):
        """-> y (target), w (labelled weight), u (UNLABELLED mask: tremor of unknown content).

        u = 1 exactly where the loss is not applied AND the sample is not part of a hidden family's
        exclusion zone.  Those are the samples nnPU may use as unlabelled (Merlin F6 allows
        "excluded OR nnPU with a class prior"); with --nnpu-pi 0 they are simply dropped."""
        nw = self.nw
        tt = t_win + np.arange(nw) / FS
        y = np.zeros(nw, np.float32)
        w = np.zeros(nw, np.float32)
        u = np.ones(nw, np.float32)
        for te, we in zip(ev_t, ev_w):
            d = tt - te
            m = np.abs(d) <= SUPP
            if not m.any():
                continue
            y[m] = np.maximum(y[m], np.exp(-0.5 * (d[m] / SIG) ** 2))
            w[m] = np.maximum(w[m], we)
            u[m] = 0.0
        for te in hid_t:                        # hidden families: never a label, never a negative
            m = np.abs(tt - te) <= SUPP
            w[m] = 0.0
            u[m] = 0.0
        return y, w, u

    # -------------------------------------------------------------- draws
    def draw_pos(self):
        for _ in range(20):
            k = self.posdays[self.rng.integers(len(self.posdays))]
            t, wv, _f = self.b.ev[k]
            j = self.rng.integers(len(t))
            i = self.st.idx(*k)
            # place that S uniformly inside the window, edge_s from either end
            off = self.rng.uniform(self.edge_s, self.win_s - self.edge_s)
            i0 = int(round((t[j] - off - self.st.t0[i]) * FS))
            x, t_win = self._cut(*k, i0)
            if x is None:
                continue
            lo, hi = t_win - SUPP, t_win + self.win_s + SUPP
            sel = (t >= lo) & (t <= hi)
            hid = self.b.evh.get(k, np.zeros(0))
            hid = hid[(hid >= lo) & (hid <= hi)]
            y, w, u = self._target(t_win, t[sel], wv[sel], hid)
            return self._norm(x), y, w, u, k[0], t_win
        return None

    def draw_neg(self):
        for _ in range(20):
            if len(self.negdays) == 0:
                break
            k = self.negdays[self.rng.integers(len(self.negdays))]
            i = self.st.idx(*k)
            n = int(self.st.nvalid[i])
            i0 = int(self.rng.integers(0, max(1, n - self.nw)))
            x, t_win = self._cut(*k, i0)
            if x is None:
                continue
            g = self.b.guard.get((k[1], k[2]), np.zeros(0))
            if g.size:
                a, b = t_win - GUARD, t_win + self.win_s + GUARD
                if ((g >= a) & (g <= b)).any():
                    continue
            return (self._norm(x), np.zeros(self.nw, np.float32),
                    np.ones(self.nw, np.float32), np.zeros(self.nw, np.float32), k[0], t_win)
        return None

    def draw_hard(self):
        for _ in range(20):
            if self.rng.random() < 0.5 and len(self.hard_sta):
                r = self.hard_sta[self.rng.integers(len(self.hard_sta))]
                sta, year, doy, t = str(r[0]), int(r[1]), int(r[2]), float(r[3])
            elif len(self.hard_net):
                r = self.hard_net[self.rng.integers(len(self.hard_net))]
                year, doy, t = int(r[0]), int(r[1]), float(r[2])
                sta = self.allstas[self.rng.integers(len(self.allstas))]
            else:
                break
            if sta not in self.b.stas or not self.st.has(sta, year, doy):
                continue
            i = self.st.idx(sta, year, doy)
            off = self.rng.uniform(self.edge_s, self.win_s - self.edge_s)
            i0 = int(round((t - off - self.st.t0[i]) * FS))
            x, t_win = self._cut(sta, year, doy, i0)
            if x is None:
                continue
            g = self.b.guard.get((year, doy), np.zeros(0))
            if g.size:
                a, b = t_win - GUARD, t_win + self.win_s + GUARD
                if ((g >= a) & (g <= b)).any():
                    continue                     # never teach "0" next to a real event
            return (self._norm(x), np.zeros(self.nw, np.float32),
                    np.ones(self.nw, np.float32), np.zeros(self.nw, np.float32), sta, t_win)
        return None

    def batch(self, n):
        X = np.zeros((n, 3, self.nw), np.float32)
        Y = np.zeros((n, self.nw), np.float32)
        W = np.zeros((n, self.nw), np.float32)
        U = np.zeros((n, self.nw), np.float32)
        S = []
        k = 0
        while k < n:
            q = self.rng.random()
            r = self.draw_pos() if q < self.pos_frac else (
                self.draw_hard() if q < self.pos_frac + self.hard_frac else self.draw_neg())
            if r is None:
                r = self.draw_pos()
            if r is None:
                continue
            X[k], Y[k], W[k], U[k] = r[0], r[1], r[2], r[3]
            S.append(r[4])
            k += 1
        return X, Y, W, U, S


# ---------------------------------------------------------------- torch plumbing ---------------
class TorchWindows:
    """Infinite iterable of batches, built in worker processes so the GPU never waits."""

    def __init__(self, years, stas, batch, win_s=60.0, seed=0, w_pred=0.5,
                 pos_frac=0.5, hard_frac=0.15, workers=8):
        self.args = dict(years=years, stas=stas, batch=batch, win_s=win_s, seed=seed,
                         w_pred=w_pred, pos_frac=pos_frac, hard_frac=hard_frac)
        self.workers = workers

    def __iter__(self):
        import torch
        from torch.utils.data import DataLoader, IterableDataset
        A = self.args

        class DS(IterableDataset):
            def __iter__(self):
                info = torch.utils.data.get_worker_info()
                wid = 0 if info is None else info.id
                bank = LabelBank(A["years"], A["stas"], w_pred=A["w_pred"])
                store = P.DayStore()
                smp = WindowSampler(bank, store, win_s=A["win_s"],
                                    rng=np.random.default_rng(A["seed"] * 1000 + wid),
                                    pos_frac=A["pos_frac"], hard_frac=A["hard_frac"])
                while True:
                    X, Y, W, U, _ = smp.batch(A["batch"])
                    yield (torch.from_numpy(X), torch.from_numpy(Y), torch.from_numpy(W),
                           torch.from_numpy(U))

        dl = DataLoader(DS(), batch_size=None, num_workers=self.workers,
                        persistent_workers=self.workers > 0, prefetch_factor=4 if self.workers else None)
        return iter(dl)
