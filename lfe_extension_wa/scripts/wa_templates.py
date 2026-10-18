#!/usr/bin/env python
"""Stage 3 (and trial re-templating): PERMANENT-station templates from each family's OWN detections.

For every family and every permanent station-epoch (bridge node) that has data at >= MIN_EVENTS (30,
e5r_templates.py) of the family's detections, stack the long window [t-20, t+32] s of the 20 Hz
1-8 Hz cache, each event normalised to unit maximum across its 3 components (the Royer-chain restack
rule, ab_grow.restack/unit_max3), all of the event window inside valid data.

  clear arrival  the catalogue's rule on the stack, 2-8 Hz Butterworth-4 designed at fs = 20 Hz
                 (NOT the 40-sps design that produced the 1-4 Hz slip), over the catalogue's window
                 [t-12, t+32]: imp = peak/median of the component-summed Hilbert envelope > 5 and
                 time above half max < 3 s.  Only clear templates are used by the matched filter;
                 every (family, node) keeps its row in the QC table.
  S pick         argmax of the 0.5-s-smoothed 3-C envelope of the 1-8 Hz stack in [t-5, t+20] s
  template       12 s from S - 0.5 s, each live channel L2-normalised (frozen chain)
  weight         S-SNR in dB (4 s from S-0.5 vs the 4 s 15 s earlier, e5r_templates.s_snr_db),
                 w = 10^(snr/20) capped at 3 x median, / median (e5r_templates.py)

Default (validation) events: the BUILD half (2-h blocks, e4 --half) of the families' array-era
detections inside the downloaded windows.  --events <parquet> (fid, t) re-templates from other events
(the trial's new station epochs, from pass-1 permanent detections).

  python scripts/wa_templates.py --tag build
  python scripts/wa_templates.py --tag t2012new --events results/trial/pass1_events.parquet \
         --only-new-nodes results/trial/new_nodes.csv
"""
from __future__ import annotations
import os, sys, argparse, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

LOG = f"{W.STUDY}/logs/templates.log"
NL = int(W.LONG_LEN * W.FS)          # 1040
IPRE = int(W.LONG_PRE * W.FS)        # 400 = index of t in the long window


def stack_day(args):
    date, ev, stations, only_nodes = args
    D0 = pd.Timestamp(date).timestamp()
    out = {}
    for ns in stations:
        p = W.cache_path(ns, date)
        if not os.path.exists(p):
            continue
        z = np.load(p)
        node = str(z["node"])
        if only_nodes is not None and node not in only_nodes:
            continue
        X, m, live = z["x"], z["m"], z["live"]
        i0 = np.round((ev.t.to_numpy() - D0) * W.FS).astype(int) - IPRE
        ok = (i0 >= 0) & (i0 + NL <= W.NDAY)
        if not ok.any():
            continue
        cm = np.concatenate([[0], np.cumsum(~m)])
        ii = i0[ok]
        ok2 = (cm[ii + NL] - cm[ii]) == 0
        fids = ev.fid.to_numpy()[ok][ok2]
        ii = ii[ok2]
        if len(ii) == 0:
            continue
        idx = ii[:, None] + np.arange(NL)[None]
        seg = X[:, idx].transpose(1, 0, 2).astype(np.float64)      # (nev, 3, NL)
        mx = np.abs(seg).max(axis=(1, 2))
        good = mx > 0
        seg = seg[good] / mx[good, None, None]
        fids = fids[good]
        for f in np.unique(fids):
            s = seg[fids == f]
            out[(int(f), ns, node)] = (s.sum(0).astype(np.float32), len(s), live.copy())
    return date, out


def env_smooth(st, smooth_s=0.5):
    from scipy.signal import hilbert
    env = np.sqrt((np.abs(hilbert(st, axis=-1)) ** 2).sum(0))
    k = int(round(smooth_s * W.FS))
    return np.convolve(env, np.ones(k) / k, mode="same")


def s_snr_db(st, s_idx, win_s=4.0, noise_off_s=15.0):
    n = int(round(win_s * W.FS)); a = max(s_idx - int(0.5 * W.FS), 0); sig = st[:, a:a + n]
    b = max(a - int(noise_off_s * W.FS), 0); noi = st[:, b:b + n]
    if sig.size == 0 or noi.size == 0 or (noi ** 2).mean() <= 0:
        return np.nan
    return float(10 * np.log10((sig ** 2).mean() / (noi ** 2).mean()))


def make_templates(acc, stations):
    site = dict(zip(stations.netsta, stations.site))
    rows, TPL = [], {}
    for (f, ns, node), (s, n, live) in acc.items():
        st = s / n
        chl = np.array([bool(live[c]) and float((st[c] ** 2).sum()) > 0 for c in range(3)])
        row = dict(fid=f, netsta=ns, node=node, site=site.get(ns, ns), n_events=n, nch=int(chl.sum()),
                   live="".join("1" if c else "0" for c in chl))
        if n < W.MIN_EVENTS or not chl.any():
            row.update(clear=False, imp=np.nan, dur=np.nan, why="n_events<30" if n < W.MIN_EVENTS else "dead")
            rows.append(row); continue
        st28 = W.bp(st.astype(np.float64), W.FS, 2.0, 8.0)
        i_a = int((W.LONG_PRE - W.CLEAR_PRE) * W.FS)
        imp, dur = W.env_rule(st28[chl][:, i_a:i_a + int(W.CLEAR_LEN * W.FS)], W.FS)
        clear = bool(np.isfinite(imp) and imp > 5 and dur < 3)
        env = env_smooth(st[chl].astype(np.float64))
        a_, b_ = IPRE + int(-5 * W.FS), IPRE + int(20 * W.FS)
        k = a_ + int(np.argmax(env[a_:b_]))
        j0 = k - int(W.PRE * W.FS)
        tpl = st[:, j0:j0 + W.NT].astype(np.float64).copy()
        tpl[~chl] = 0.0
        nrm = np.sqrt((tpl ** 2).sum(1))
        chl2 = chl & (nrm > 0)
        tpl[chl2] /= nrm[chl2, None]
        snr = s_snr_db(st[chl].astype(np.float64), k)
        row.update(clear=clear, imp=imp, dur=dur, why="" if clear else "not clear",
                   s_off=(k - IPRE) / W.FS, offset_s=(j0 - IPRE) / W.FS, snr_db=snr)
        rows.append(row)
        TPL[(f, ns, node)] = tpl.astype(np.float32)
    Q = pd.DataFrame(rows)
    return Q, TPL


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="build")
    ap.add_argument("--events", default="", help="parquet with fid,t (default: truth build half)")
    ap.add_argument("--only-new-nodes", default="", help="csv with a `node` column: stack only these nodes")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--weights-ref", default="", help="QC csv whose median raw SNR weight normalises w")
    a = ap.parse_args()
    S = W.load_stations()
    if a.events:
        EV = pd.read_parquet(a.events)[["fid", "t"]]
    else:
        T = pd.read_parquet(f"{W.SETUP}/truth.parquet")
        EV = T[T.half == 0][["fid", "t"]]
    EV = EV.copy()
    EV["date"] = pd.to_datetime(EV.t, unit="s").dt.strftime("%Y-%m-%d")
    only = set(pd.read_csv(a.only_new_nodes).node) if a.only_new_nodes else None
    W.log(LOG, f"[tpl:{a.tag}] {len(EV)} events on {EV.date.nunique()} days, {EV.fid.nunique()} families"
          + (f", nodes {sorted(only)}" if only else ""))
    acc = {}
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        fut = [ex.submit(stack_day, (d, g[["fid", "t"]], list(S.netsta), only)) for d, g in EV.groupby("date")]
        for i, f in enumerate(as_completed(fut)):
            d, out = f.result()
            for k, (s, n, live) in out.items():
                if k in acc:
                    s0, n0, l0 = acc[k]; acc[k] = (s0 + s, n0 + n, l0 | live)
                else:
                    acc[k] = (s.astype(np.float64), n, live)
            if (i + 1) % 10 == 0:
                W.log(LOG, f"[tpl:{a.tag}] {i+1}/{len(fut)} days {time.time()-t0:.0f}s")
    Q, TPL = make_templates(acc, S)
    for c_ in ["fid", "netsta", "node", "site", "n_events", "nch", "live", "clear", "imp", "dur", "why", "s_off",
               "offset_s", "snr_db"]:
        if c_ not in Q:                  # e.g. every stack < 30 events: no S pick / SNR columns
            Q[c_] = pd.Series(dtype=float if c_ not in ("netsta", "node", "site", "live", "why") else object)
    Q["clear"] = Q.clear.fillna(False).astype(bool)
    # SNR weights (e5r rule) over the usable (clear) templates
    use = Q.clear.fillna(False).astype(bool)
    wraw = 10 ** (Q.snr_db / 20.0)
    ref = float(np.median(wraw[use])) if use.any() else 1.0
    if a.weights_ref:
        R = pd.read_csv(a.weights_ref)
        R = R[R.clear.fillna(False).astype(bool)]
        ref = float(np.median(10 ** (R.snr_db / 20.0)))
    Q["w"] = np.minimum(wraw, 3 * ref) / ref
    od = f"{W.RES}/templates"
    os.makedirs(od, exist_ok=True)
    Q.to_csv(f"{od}/tpl_qc_{a.tag}.csv", index=False)
    keys = [k for k in TPL]
    np.savez_compressed(f"{od}/stacks_{a.tag}.npz",
                        keys=np.array([f"{f}|{ns}|{node}" for f, ns, node in acc]),
                        stacks=np.stack([(v[0] / v[1]).astype(np.float32) for v in acc.values()]) if acc else np.zeros((0, 3, NL)))
    np.savez_compressed(f"{od}/tpl_{a.tag}.npz", keys=np.array([f"{f}|{ns}|{node}" for f, ns, node in keys]),
                        tpl=np.stack([TPL[k] for k in keys]) if keys else np.zeros((0, 3, W.NT), np.float32))
    c = Q[use]
    fam_sites = c.groupby("fid").site.nunique()
    W.log(LOG, f"[tpl:{a.tag}] (fam,node) stacks {len(Q)}; clear {int(use.sum())}; families with >=1 clear "
          f"{c.fid.nunique()}, >=2 clear sites {(fam_sites >= 2).sum()}; {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
