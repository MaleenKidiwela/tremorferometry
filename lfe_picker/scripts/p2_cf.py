#!/usr/bin/env python
"""Write the per-station-day characteristic function of a P2 model to disk (float16 .npy per day).

  python scripts/p2_cf.py --model results/p2/models/p2core.pt --year 2012 --out results/p2/cf_p2core_2012

--sid known   use each station's own embedding (training stations only)
--sid unknown use the "unknown station" token everywhere (the default: keeps the CF station-agnostic
              and is the only legal choice at B926/NLLB/LRIV, whose embeddings were never trained)
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P
import p_models as M
import p_scan as S

FS = P.FS


@torch.no_grad()
def cf_day_sid(model, meta, x, sid, dev="cuda", batch=1024):
    nw = int(round(meta["win_s"] * FS))
    step = nw // 3
    N = x.shape[1]
    starts = np.arange(0, max(1, N - nw + 1), step, dtype=np.int64)
    if starts[-1] != N - nw:
        starts = np.append(starts, N - nw)
    out = np.zeros(N, np.float32)
    cov = np.zeros(N, np.bool_)
    X = torch.from_numpy(np.ascontiguousarray(x)).to(dev)
    for i0 in range(0, len(starts), batch):
        sl = starts[i0:i0 + batch]
        idx = torch.from_numpy(sl).to(dev)[:, None] + torch.arange(nw, device=dev)[None, :]
        w = X[:, idx.reshape(-1)].reshape(3, len(sl), nw).permute(1, 0, 2).float()
        w = w / w.abs().amax(dim=(1, 2), keepdim=True).clamp_min(1e-12)
        sids = torch.full((len(sl),), sid, dtype=torch.long, device=dev)
        with torch.amp.autocast("cuda", dtype=torch.bfloat16):
            p = model(w, sids).float().cpu().numpy()
        a, b = nw // 3, 2 * (nw // 3)
        for k, s in enumerate(sl):
            lo, hi = (0 if s == 0 else a), (nw if s == starts[-1] else b)
            out[s + lo:s + hi] = np.maximum(out[s + lo:s + hi], p[k, lo:hi])
            cov[s + lo:s + hi] = True
    out[~cov] = 0.0
    out[(np.asarray(x, np.float32) == 0).all(0)] = 0.0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sid", default="unknown", choices=["unknown", "known"])
    ap.add_argument("--doys", default="")
    a = ap.parse_args()
    assert a.year != 2013, "2013 is sealed"
    os.makedirs(a.out, exist_ok=True)
    z = torch.load(a.model, map_location="cuda", weights_only=False)
    model = M.PhaseNetS(nsta=z.get("nsta", 0))
    model.load_state_dict(z["state"])
    model.cuda().eval()
    sid_of = z.get("sid_of", {s: i for i, s in enumerate(P.S15)})
    unk = z.get("nsta", len(P.S15))
    st = P.DayStore()
    di = st.di[st.di.year == a.year]
    if a.doys:
        di = di[di.doy.isin([int(x) for x in a.doys.split(",")])]
    print(f"[cf] {len(di)} station-days, sid={a.sid}", flush=True)
    t0 = time.time()
    for k, r in enumerate(di.itertuples()):
        f = f"{a.out}/{r.sta}_{r.year}_{r.doy}.npy"
        if os.path.exists(f):
            continue
        i = st.idx(r.sta, r.year, r.doy)
        n = int(st.nvalid[i])
        x = np.asarray(st.get_i(i)[:, :n], np.float32)
        sid = unk if (a.sid == "unknown" or r.sta in P.HELD_OUT_STATIONS) else sid_of[r.sta]
        cf = cf_day_sid(model, z, x, sid)
        np.save(f, cf.astype(np.float16))
        if (k + 1) % 50 == 0:
            print(f"  {k+1}/{len(di)} {time.time()-t0:.0f}s", flush=True)
    print(f"[done] {time.time()-t0:.0f}s -> {a.out}")


if __name__ == "__main__":
    main()
