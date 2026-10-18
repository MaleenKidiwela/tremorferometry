#!/usr/bin/env python
"""P3D STEP 4a -- the CF timing KILL TEST on the dense network (GPU).

At every hp event of the TEST episodes, read the frozen P2 characteristic function
(`results/p2/models/p2nnpu.pt`) at each station and record the CF peak time inside +-5 s of the
Crosson-predicted S.  The arm closes if the median |residual| exceeds 1.5 s at more than half of
the sites whose median hypocentral distance to the family field is <= 60 km.

Station identity: every station is run through the UNKNOWN-station embedding (index `nsta`), which
is what `p2_cf.py --sid unknown` does.  The code path is `p2_cf.cf_day_sid(model, meta, x, sid)`
with `sid` passed explicitly -- it NEVER indexes by S15 position, so a new station is legal here
(verified by reading p2_cf.main: `sid_of` is consulted only under `--sid known`).

Days are read from the 20 Hz caches through `s_common` (ZNE-rotated with the StationXML azimuths),
not through `p_common.DayStore`, whose index is filtered to S15.

  python scripts/s_step4_cf.py --stage cf     # writes results/p3d/cf/<sta>_<year>_<doy>.npy
  python scripts/s_step4_cf.py --stage 4a     # writes results/p3d/step4a_residuals.csv
"""
from __future__ import annotations
import os, sys, time, json, argparse, warnings
import numpy as np
import pandas as pd
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                        # noqa: E402

CFDIR = f"{S.P3D}/cf"
MODEL = f"{S.PICK}/results/p2/models/p2nnpu.pt"
YEARS = list(range(2019, 2027))


def target_days(episodes, quiet_months):
    W = S.sse_windows()
    days = set()
    for r in W[W.id.isin(episodes)].itertuples():
        for t in pd.date_range(r.start, r.end, freq="D"):
            days.add((int(t.year), int(t.dayofyear)))
    for m in quiet_months:
        for t in pd.date_range(f"{m}-01", pd.Timestamp(f"{m}-01") + pd.offsets.MonthEnd(0)):
            days.add((int(t.year), int(t.dayofyear)))
    return sorted(days)


@torch.no_grad()
def cf_at_events(model, x, idx_ev, sid, nw, dev="cuda", half=100, batch=1024):
    """The frozen CF (p2_cf.cf_day_sid) evaluated ONLY on the window grid that covers +-half samples
    around each event, by the SAME rule.  cf_day_sid tiles the day with windows every nw//3 samples
    and an interior window starting at s writes only samples [s+nw//3, s+2nw//3); consecutive
    windows tile exactly, so the covering window of sample j is s = step*(floor(j/step) - 1).
    Evaluating only those windows reproduces the full-day CF sample-for-sample in the interior."""
    N = x.shape[1]
    step = nw // 3
    last = N - nw
    need = set()
    for i in idx_ev:
        for j in range(int(i) - half, int(i) + half + 1, max(1, step // 4)):
            need.add(int(np.clip(step * (j // step - 1), 0, last)))
        need.add(int(np.clip(step * ((int(i) + half) // step - 1), 0, last)))
    need = np.array(sorted(need), dtype=np.int64)
    if not len(need):
        return None
    X = torch.from_numpy(np.ascontiguousarray(x)).to(dev)
    P = np.zeros((len(need), nw), np.float32)
    for k0 in range(0, len(need), batch):
        sl = need[k0:k0 + batch]
        ix = torch.from_numpy(sl).to(dev)[:, None] + torch.arange(nw, device=dev)[None, :]
        w = X[:, ix.reshape(-1)].reshape(3, len(sl), nw).permute(1, 0, 2).float()
        w = w / w.abs().amax(dim=(1, 2), keepdim=True).clamp_min(1e-12)
        sids = torch.full((len(sl),), sid, dtype=torch.long, device=dev)
        with torch.amp.autocast("cuda", dtype=torch.bfloat16):
            P[k0:k0 + len(sl)] = model(w, sids).float().cpu().numpy()
    del X
    pos = {int(v): k for k, v in enumerate(need)}
    zero = (np.asarray(x, np.float32) == 0).all(0)
    out = np.zeros((len(idx_ev), 2 * half + 1), np.float32)
    for e, i in enumerate(idx_ev):
        for m, j in enumerate(range(int(i) - half, int(i) + half + 1)):
            c = int(np.clip(step * (j // step - 1), 0, last))
            v = 0.0
            if c in pos and 0 <= j - c < nw:
                lo = 0 if c == 0 else step
                hi = nw if c == last else 2 * step
                if lo <= j - c < hi:
                    v = float(P[pos[c], j - c])
            out[e, m] = 0.0 if zero[j] else v
    return out


def stage_4a(a):
    import p_models as M
    ST = (S.station_table() if a.variant == "v1"
          else pd.read_csv(f"{S.P3D}/stations_p3d_v2.csv").set_index("sta"))
    SI = S.site_table(ST)
    stas = [s for s in ST.index if ST["mode"][s] == "3C"]
    if a.stations:
        stas = [s for s in a.stations.split(",") if s in stas]
    z = torch.load(MODEL, map_location="cuda", weights_only=False)
    model = M.PhaseNetS(nsta=z.get("nsta", 0))
    model.load_state_dict(z["state"])
    model.cuda().eval()
    unk = int(z.get("nsta", 15))
    nw = int(round(float(z["win_s"]) * S.FS))
    print(f"[4a] model win_s={z['win_s']} nw={nw}, unknown-station sid={unk} "
          f"(sid passed explicitly; never indexed by S15 position)", flush=True)
    tts = S.fam_tt(stas, ST)
    H = S.hp_events(0.8, episodes=a.episodes.split(","))
    yr, doy = S.yd_of(H.t.values)
    H = H.assign(year=yr, doy=doy)
    idx = S.cache_index(years=YEARS, stas=stas)
    half = int(round(5.0 * S.FS))
    rows, t0 = [], time.time()
    for s in stas:
        ns = 0
        for (y, d), g in H.groupby(["year", "doy"]):
            v = S.load_day_zne(s, int(y), int(d), idx)
            if v is None:
                continue
            tt0, x = v
            x = np.asarray(x, np.float32)
            ts = g.t.to_numpy() + np.array([tts[(int(b), s)] for b in g.bkey])
            iev = np.rint((ts - tt0) * S.FS).astype(np.int64)
            ok = (iev - half - nw > 0) & (iev + half + nw < x.shape[1])
            if not ok.any():
                continue
            cf = cf_at_events(model, x, iev[ok], unk, nw, half=half)
            if cf is None:
                continue
            for e, r in enumerate(g[ok].itertuples()):
                seg = cf[e]
                if seg.max() <= 0:
                    continue
                j = int(np.argmax(seg))
                rows.append(dict(sta=s, site=ST.site[s], bkey=int(r.bkey), t=r.t,
                                 episode=r.episode, resid_s=(j - half) / S.FS,
                                 cf_peak=float(seg[j]), cf_edge=float(seg[[0, -1]].mean())))
                ns += 1
        print(f"  [4a] {s}: {ns} events  ({time.time()-t0:.0f}s)", flush=True)
    D = pd.DataFrame(rows)
    if not len(D):
        print("[4a] no events -- nothing written (check the roster/variant and the cache)")
        return
    D.to_csv(f"{S.P3D}/step4a_residuals{a.out_tag}.csv", index=False)
    le60 = set(SI[SI.geom_le60].site)
    su = []
    for site, g in D.groupby("site"):
        su.append(dict(site=site, n=len(g), le60=site in le60,
                       med_abs_resid_s=float(np.median(np.abs(g.resid_s))),
                       med_resid_s=float(np.median(g.resid_s)),
                       frac_within_1s=float((np.abs(g.resid_s) <= 1.0).mean()),
                       frac_within_0p75s=float((np.abs(g.resid_s) <= 0.75).mean()),
                       med_cf_peak=float(np.median(g.cf_peak)),
                       med_cf_edge=float(np.median(g.cf_edge))))
    SU = pd.DataFrame(su).sort_values("med_abs_resid_s")
    SU.to_csv(f"{S.P3D}/step4a_summary{a.out_tag}.csv", index=False)
    print(SU.to_string(index=False))
    k = SU[SU.le60]
    bad = int((k.med_abs_resid_s > 1.5).sum())
    verdict = "CLOSE" if len(k) and bad > len(k) / 2 else "SURVIVES"
    print(f"[4a] {bad}/{len(k)} of the <=60 km sites have median |residual| > 1.5 s -> {verdict}; "
          f"a CF with NO timing information peaks uniformly over +-5 s -> chance median = 2.50 s")
    json.dump(dict(verdict=verdict, n_le60=int(len(k)), n_bad=bad,
                   chance_median_abs_resid_s=2.5,
                   sites={r.site: r.med_abs_resid_s for r in SU.itertuples()}),
              open(f"{S.P3D}/step4a_verdict{a.out_tag}.json", "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["4a"], default="4a")
    ap.add_argument("--episodes", default=",".join(S.TEST_EPISODES))
    ap.add_argument("--stations", default="")
    ap.add_argument("--variant", default="v1", choices=["v1", "v2"])
    ap.add_argument("--out-tag", default="")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    stage_4a(a)


if __name__ == "__main__":
    main()
