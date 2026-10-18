#!/usr/bin/env python
"""The two P1b baselines, in PyTorch.

LIN  — Lin et al. (2024) `unet_tools.make_large_unet(fac, sr=40, ncomps=3)` re-implemented layer for
       layer: 15-s windows at HIS rate, 6 input channels (log|amp| and sign per component, in his
       Z,E,N order), the Conv1D 21/15/11 encoder, the Flatten->Dense(75)->Reshape bottleneck (which
       is what makes his model non-convolutional and fixes its input length), skip concatenations,
       a per-sample sigmoid, binary cross-entropy, Adam 1e-4.  Our data is the frozen 1-8 Hz chain
       at 20 Hz, so windows are polyphase-upsampled 20 -> 40 Hz to enter his architecture; the band
       is unchanged and 20 Hz already oversamples it, so nothing is invented.

PNET — a PhaseNet-style fully convolutional U-Net native to our chain: 60-s windows at 20 Hz, 3 raw
       peak-normalised components in, 4 stride-2 levels, GroupNorm + residual blocks, two sigmoid
       output channels (S, and P as an auxiliary head).  Fully convolutional, so it runs over a
       whole station-day in one pass.
"""
from __future__ import annotations
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

EPS = 1e-6


# ------------------------------------------------------------------ Lin -----------------------
def zen2inp(x):
    """(B,3,n) peak-normalised Z,N,E  ->  (B,6,n) Lin features in his Z,E,N channel order."""
    z, n_, e = x[:, 0], x[:, 1], x[:, 2]
    out = []
    for c in (z, e, n_):
        out += [torch.log(c.abs() + EPS), torch.sign(c)]
    return torch.stack(out, 1)


class LinUNet(nn.Module):
    """make_large_unet(fac, sr=40, ncomps=3): input (B, 6, 600), output (B, 600) in [0,1]."""

    def __init__(self, fac=2.0, n=600):
        super().__init__()
        self.n = n
        f1, f2, f3 = int(fac * 32), int(fac * 64), int(fac * 128)
        self.c1 = nn.Conv1d(6, f1, 21, padding=10)
        self.c2 = nn.Conv1d(f1, f2, 15, padding=7)
        self.c3 = nn.Conv1d(f2, f3, 11, padding=5)
        nb = n // 8
        self.dense = nn.Linear(f3 * nb, nb)
        self.u3 = nn.Conv1d(1, f3, 11, padding=5)
        self.u2 = nn.Conv1d(f3 + f3, f2, 15, padding=7)
        self.u1 = nn.Conv1d(f2 + f2, f1, 21, padding=10)
        self.out = nn.Conv1d(f1 + f1, 1, 21, padding=10)

    def forward(self, x):
        l1 = F.relu(self.c1(x))
        l2 = F.relu(self.c2(F.max_pool1d(l1, 2)))
        l3 = F.relu(self.c3(F.max_pool1d(l2, 2)))
        b = F.relu(self.dense(F.max_pool1d(l3, 2).flatten(1))).unsqueeze(1)
        u = F.interpolate(F.relu(self.u3(b)), scale_factor=2, mode="nearest")
        u = F.interpolate(F.relu(self.u2(torch.cat([u, l3], 1))), scale_factor=2, mode="nearest")
        u = F.interpolate(F.relu(self.u1(torch.cat([u, l2], 1))), scale_factor=2, mode="nearest")
        return torch.sigmoid(self.out(torch.cat([u, l1], 1))).squeeze(1)


# ------------------------------------------------------------------ PhaseNet-style ------------
class Block(nn.Module):
    def __init__(self, cin, cout, k=9):
        super().__init__()
        self.a = nn.Conv1d(cin, cout, k, padding=k // 2)
        self.na = nn.GroupNorm(8, cout)
        self.b = nn.Conv1d(cout, cout, k, padding=k // 2)
        self.nb = nn.GroupNorm(8, cout)
        self.sk = nn.Conv1d(cin, cout, 1) if cin != cout else nn.Identity()

    def forward(self, x):
        h = F.relu(self.na(self.a(x)))
        h = self.nb(self.b(h))
        return F.relu(h + self.sk(x))


class PhaseNetS(nn.Module):
    """(B, 3, n) -> (B, n) sigmoid S characteristic function.

    P is NOT predicted: Merlin F5(c) rules our P labels unusable (residuals +0.5..+2.5 s beyond
    60 km, Lin's own sigma 4.2 s), so an auxiliary P head would be supervised by noise.  Optional
    station embedding is added to the encoder trunk (used under station dropout in P2)."""

    def __init__(self, ch=(24, 48, 96, 192), nsta=0, emb=16):
        super().__init__()
        c0 = ch[0]
        self.inp = Block(3, c0)
        self.down = nn.ModuleList([Block(ch[i], ch[i + 1]) for i in range(len(ch) - 1)])
        self.mid = Block(ch[-1], ch[-1])
        self.up = nn.ModuleList([Block(ch[i + 1] + ch[i], ch[i]) for i in reversed(range(len(ch) - 1))])
        self.head = nn.Conv1d(c0, 1, 1)
        self.nsta = nsta
        if nsta:
            self.emb = nn.Embedding(nsta + 1, emb)     # index nsta == "unknown station" (dropout)
            self.proj = nn.Conv1d(emb, c0, 1)

    def forward(self, x, sid=None):
        h = self.inp(x)
        if self.nsta and sid is not None:
            h = h + self.proj(self.emb(sid).unsqueeze(-1))
        skips = []
        for d in self.down:
            skips.append(h)
            h = d(F.max_pool1d(h, 2))
        h = self.mid(h)
        for u, s in zip(self.up, reversed(skips)):
            h = F.interpolate(h, size=s.shape[-1], mode="nearest")
            h = u(torch.cat([h, s], 1))
        return torch.sigmoid(self.head(h)).squeeze(1)


# ------------------------------------------------------------------ helpers -------------------
def upsample2(x):
    """20 -> 40 Hz polyphase upsampling on the GPU (windowed-sinc, up=2)."""
    if not hasattr(upsample2, "k"):
        from scipy.signal import firwin
        k = firwin(31, 0.5).astype(np.float32) * 2.0
        upsample2.k = torch.from_numpy(k)
    k = upsample2.k.to(x.device, x.dtype)
    B, C, n = x.shape
    z = torch.zeros(B, C, 2 * n, device=x.device, dtype=x.dtype)
    z[:, :, ::2] = x
    return F.conv1d(z.reshape(B * C, 1, 2 * n), k.view(1, 1, -1), padding=len(k) // 2
                    ).reshape(B, C, 2 * n)[:, :, :2 * n]


def make(name, **kw):
    if name == "lin":
        return LinUNet(**kw)
    if name == "pnet":
        return PhaseNetS(**kw)
    raise ValueError(name)
