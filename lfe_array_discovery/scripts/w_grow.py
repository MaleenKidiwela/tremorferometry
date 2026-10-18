#!/usr/bin/env python
"""Stage 3 -- HANDOFF (align_stack -> pass-0 templates), GROWTH, SURVIVAL, SCRAMBLE NULL.

Thin wrapper over the P4-0b growth engine (`v_engine`), with Merlin's four prerequisites installed
in process by `w_common.patch_v` (NO 8.0 floor, data-driven sites, valid-only station_z, MIN_SITES
4, per-day availability never intersected).  Nothing in lfe_picker is edited.

  handoff   per seed cluster: cut HANDOFF_CUT_LEN_S windows at every station at the cluster member
            times (the SAME absolute time at every station -- the moveout is embedded in the stack,
            which is all growth needs), align_stack (2 iterations, +- 1.5 s), find the stack's
            ENERGY ONSET, cut a 20-s pass-0 template from onset - 2 s.  `off` = template start
            minus the cluster time axis.
  scramble  the TEMPLATE-SCRAMBLE null: the same off, but each station's template is cut from a
            DIFFERENT random tremor window.  Grown identically; real >= 3 x scramble.
  quiet     FA-only T per template on the quiet days (<= 0.1 false hits per template-day, grid
            0.25, engineering floor 5.0); reports where the floor binds.
  grow      growth over every bed day with per-day availability.
  survive   >= 20 independent hits (60 s), >= 3 distinct days, <= 20 % within +- 10 min of the seed
            hour, and scramble ratio >= 3.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402
import s_common as S                                                        # noqa: E402


def _cap(days):
    """W_MAX_DAYS=<n> truncates every day list -- SMOKE TESTS ONLY; never set in the queue."""
    k = int(os.environ.get("W_MAX_DAYS", "0") or 0)
    return days[:k] if k else days


def _cut(loader, sta, t_abs, n):
    # F1: a window inside a vetoed station-hour (flat / H-V ratio / 1-8 Hz tone) is never stacked
    # into a template -- otherwise SHVB's 4.76 Hz tone would be built into every 1c' template.
    if W.vetoed_at(sta, t_abs):
        return None
    ts = pd.Timestamp(float(t_abs), unit="s")
    v = loader.get(sta, int(ts.year), int(ts.dayofyear))
    if v is None:
        return None
    t0, X = v
    i0 = int(round((t_abs - t0) * W.FS))
    if i0 < 0 or i0 + n > X.shape[1]:
        return None
    seg = X[:, i0:i0 + n].astype(np.float64)
    if not np.isfinite(seg).all() or (seg ** 2).sum(1).min() <= 0:
        return None
    r = np.sqrt((seg ** 2).sum(-1, keepdims=True))
    return (seg / np.maximum(r, 1e-30)).astype(np.float32)


def onset_index(st, nom):
    """Energy onset of a 3-C stack, searched in HANDOFF_ONSET_SEARCH around the cluster time."""
    env = S.envelope(st)
    a = max(nom + int(round(W.HANDOFF_ONSET_SEARCH[0] * W.FS)), 0)
    b = min(nom + int(round(W.HANDOFF_ONSET_SEARCH[1] * W.FS)), len(env))
    return int(a + np.argmax(env[a:b]))


def handoff(members, seeds, out_bank, shift=None, scramble=False, seed_rng=20260920,
            tremor_times=None):
    """Seed clusters -> pass-0 template bank."""
    V, E = W.patch_v(shift=shift)
    loader = E._W_LOADER
    mem, of, order, det = W.sites()
    n_cut = int(round(W.HANDOFF_CUT_LEN_S * W.FS))
    nom = int(round(W.HANDOFF_CUT_PRE_S * W.FS))
    n = int(round(W.TPL0_LEN_S * W.FS))
    rng = np.random.default_rng(seed_rng)
    cids = list(seeds.cid)
    Wb = np.full((len(cids), len(det), 3, n), np.nan, np.float32)
    off = np.full((len(cids), len(det)), np.nan)
    diag, tids = [], []
    for bi, cid in enumerate(cids):
        t = np.sort(members[members.cid == cid].t.to_numpy(float))
        if len(t) > W.MAX_STACK:      # bounded stacking cost; evenly spaced over the cluster
            t = t[np.linspace(0, len(t) - 1, W.MAX_STACK).round().astype(int)]
        got = 0
        for si, s in enumerate(det):
            segs = [w for w in (_cut(loader, s, tt - W.HANDOFF_CUT_PRE_S, n_cut) for tt in t)
                    if w is not None]
            if len(segs) < W.HANDOFF_MIN_SEGMENTS:
                continue
            st, _, _ = S.align_stack(np.stack(segs), iters=W.ALIGN_ITERS, lag=W.ALIGN_LAG)
            i_on = onset_index(st, nom)
            a = i_on - int(round(W.TPL0_PRE_S * W.FS))
            if a < 0 or a + n > st.shape[1]:
                continue
            o = -W.HANDOFF_CUT_PRE_S + a / W.FS
            if scramble:
                # TEMPLATE-SCRAMBLE NULL: same off, a DIFFERENT random tremor window per station
                for _ in range(20):
                    tr = float(rng.choice(tremor_times))
                    w = _cut(loader, s, tr, n)
                    if w is not None:
                        Wb[bi, si] = w
                        break
                else:
                    continue
            else:
                Wb[bi, si] = V.norm_channels(st[:, a:a + n])
            off[bi, si] = o
            got += 1
            diag.append(dict(cid=cid, sta=s, n_stacked=len(segs), onset_s=float((i_on - nom) / W.FS),
                             s_snr_db=float(S.s_snr_db(st, i_on))))
        tids.append(cid)
    ok = np.isfinite(Wb).all(axis=(2, 3))
    nsite = np.array([len({of[det[i]] for i in np.nonzero(ok[b])[0]}) for b in range(len(cids))])
    sel = np.nonzero(nsite >= W.MIN_SITES)[0]
    pd.DataFrame(diag).to_csv(out_bank.replace(".npz", "_diag.csv"), index=False)
    V.save_bank(out_bank, [tids[i] for i in sel], det, Wb[sel], off[sel],
                dict(kind="scramble0" if scramble else "pass0", len_s=W.TPL0_LEN_S,
                     pre_s=W.TPL0_PRE_S, n_in=len(cids),
                     n_site_median=float(np.median(nsite)) if len(nsite) else 0.0))
    W.log(f"[handoff{'/scramble' if scramble else ''}] {len(sel)}/{len(cids)} clusters with "
          f">= {W.MIN_SITES} sites -> {out_bank}")
    return out_bank


def tremor_window_times(n=4000, seed=20260920):
    """Random window start times on the CANDIDATE (tremor) hours, for the scramble null."""
    H = W.candidate_hours()
    rng = np.random.default_rng(seed)
    if not len(H):
        return np.zeros(0)
    t0 = np.array([pd.Timestamp(r.date).timestamp() + r.hour * 3600 for r in H.itertuples()])
    k = rng.integers(0, len(t0), n)
    return t0[k] + rng.uniform(60.0, 3400.0, n)


def quiet_T(bank, out, shift=None, block=32):
    V, E = W.patch_v(shift=shift)
    tids, stas, Wt, off, meta = V.load_bank(bank)
    if len(tids) == 0:
        pd.DataFrame(columns=W.THR_COLS).to_csv(out + "_thresholds.csv", index=False)
        pd.DataFrame(columns=W.HIT_COLS + ["tid"]).to_parquet(out + "_quiethits.parquet",
                                                              index=False)
        W.log("[quiet] empty bank -- 0 templates, empty thresholds written")
        return pd.DataFrame(columns=W.THR_COLS)
    Q = W.quiet_days()
    days = _cap([pd.Timestamp(d) for d in Q[Q.quiet].date])
    H = E.grow(bank, days, out + "_quiethits.parquet", thr=None, floor=W.T_ENG_FLOOR,
               block=block, tag="/quiet")
    nday = len(days)
    grid = np.arange(W.T_ENG_FLOOR, W.T_MAX + 1e-9, W.T_GRID)
    rows = []
    for bi, tid in enumerate(tids):
        v = H[H.tpl == bi].net_mad.to_numpy() if len(H) else np.zeros(0)
        fa = np.array([(v >= t).sum() / max(nday, 1) for t in grid])
        ok = np.nonzero(fa <= W.FA_PER_TPL_DAY)[0]
        T = float(grid[ok[0]]) if len(ok) else float(W.T_MAX)
        rows.append(dict(tpl=bi, tid=tid, T_seed=T, n_quiet_hits=int(len(v)), quiet_days=nday,
                         fa_at_T=float(fa[ok[0]]) if len(ok) else float(fa[-1]),
                         floor_binds=bool(len(ok) and T <= W.T_ENG_FLOOR + 1e-9),
                         saturated=bool(not len(ok))))
    Tt = pd.DataFrame(rows, columns=W.THR_COLS) if rows else pd.DataFrame(columns=W.THR_COLS)
    Tt.to_csv(out + "_thresholds.csv", index=False)
    if len(Tt):
        W.log(f"[quiet] {len(Tt)} templates / {nday} quiet days; T median {Tt.T_seed.median():.2f} "
              f"[{Tt.T_seed.min():.2f}-{Tt.T_seed.max():.2f}]; engineering floor binds for "
              f"{int(Tt.floor_binds.sum())}/{len(Tt)}; saturated {int(Tt.saturated.sum())}")
    return Tt


def grow(bank, out, thr_csv="", shift=None, block=32):
    V, E = W.patch_v(shift=shift)
    tids, stas, Wt, off, meta = V.load_bank(bank)
    if len(tids) == 0:
        pd.DataFrame(columns=W.HIT_COLS + ["tid"]).to_parquet(out + "_hits.parquet", index=False)
        W.log("[grow] empty bank -- 0 templates, empty hits written")
        return
    Tt = W.read_csv_safe(thr_csv or (out + "_thresholds.csv"), W.THR_COLS)
    T = (Tt.set_index("tpl").T_seed.reindex(range(len(tids))).fillna(W.T_ENG_FLOOR).to_numpy()
         if len(Tt) else np.full(len(tids), W.T_ENG_FLOOR))
    days = _cap(W.bed_days())
    E.grow(bank, days, out + "_hits.parquet", thr=T, floor=W.T_ENG_FLOOR, block=block, tag="/grow")


def survive(bank, out, seeds_csv, scramble_hits="", scramble_bank=""):
    V, E = W.patch_v()
    tids, stas, Wt, off, meta = V.load_bank(bank)
    if len(tids) == 0:
        pd.DataFrame(columns=W.SURV_COLS).to_csv(out + "_survival.csv", index=False)
        pd.DataFrame(columns=W.HIT_COLS + ["tid"]).to_parquet(out + "_indep.parquet", index=False)
        W.log("[survive] empty bank -- 0 templates, 0 survivors")
        return pd.DataFrame(columns=W.SURV_COLS)
    H = (pd.read_parquet(out + "_hits.parquet") if os.path.exists(out + "_hits.parquet")
         else pd.DataFrame(columns=W.HIT_COLS + ["tid"]))
    _sd = W.read_csv_safe(seeds_csv)
    Sd = _sd.set_index("cid") if len(_sd) else None
    # F5 (Merlin's FINAL ruling): survival counts ETS-WINDOW days ONLY.  The FA-calibration quiet
    # days are still GROWN over -- their independent hits are the empirical false-alarm count and
    # are reported in n_indep_quiet -- but they can no longer supply a "distinct day" or an
    # independent hit to the gate.  Flag off = the frozen 1c' behaviour, byte for byte.
    ETS = W.ets_date_set() if W.SURV_ETS_ONLY else None
    NM = None
    if scramble_hits and os.path.exists(scramble_hits) and scramble_bank:
        NH = pd.read_parquet(scramble_hits)
        if ETS is not None and len(NH):
            NH = NH[NH.date.astype(str).isin(ETS)]
        ntids, *_ = V.load_bank(scramble_bank)
        # a template whose twin was dropped upstream has n_scramble = 0, not NaN, when the flag is
        # on (the twin produced nothing); frozen default leaves it NaN, which vetoes the template.
        NM = {t: 0 for t in tids} if W.SCRAMBLE_MISSING_ZERO else {}
        for bi, g in NH.groupby("tpl"):
            NM[ntids[int(bi)]] = int(W.independent(g.sort_values("t").t.to_numpy()).sum())
    rows, indep = [], []
    for bi, tid in enumerate(tids):
        g = H[H.tpl == bi].sort_values("t") if len(H) else H
        if not len(g):
            rows.append(dict(tpl=bi, tid=tid, n_hits=0, n_indep=0, n_days=0, frac_near_seed=np.nan,
                             n_scramble=NM.get(tid, np.nan) if NM else np.nan,
                             scramble_ratio=np.nan, survives=False, n_indep_quiet=0))
            continue
        m = W.independent(g.t.to_numpy())
        gi = g[m].copy()
        gi["tid"] = tid
        indep.append(gi)                     # UNCHANGED: every independent hit, all bed days
        gc = gi[gi.date.astype(str).isin(ETS)] if ETS is not None else gi
        n_quiet = int(len(gi) - len(gc))
        ts = float(Sd.t_seed[tid]) if (Sd is not None and tid in Sd.index) else np.nan
        if W.NEAR_SEED_MODE == "either" and Sd is not None and tid in Sd.index:
            # a 2-node pair seed has two seed windows; Merlin's clause is "<= 20 % within
            # +- 10 min of EITHER seed window".
            t0s = float(Sd.t0[tid]) if "t0" in Sd.columns else ts
            t1s = float(Sd.t1[tid]) if "t1" in Sd.columns else ts
            tv = gc.t.to_numpy()
            near = (float(np.mean((np.abs(tv - t0s) <= W.NEAR_SEED_S)
                                  | (np.abs(tv - t1s) <= W.NEAR_SEED_S))) if len(tv) else np.nan)
        else:
            near = (float(np.mean(np.abs(gc.t.to_numpy() - ts) <= W.NEAR_SEED_S))
                    if (np.isfinite(ts) and len(gc)) else np.nan)
        nm = NM.get(tid, np.nan) if NM else np.nan
        ratio = (len(gc) / max(nm, 1)) if (NM and np.isfinite(nm)) else np.nan
        surv = bool(len(gc) >= W.SURV_MIN_INDEP and gc.date.nunique() >= W.SURV_MIN_DAYS
                    and (not np.isfinite(near) or near <= W.SURV_MAX_NEAR_SEED)
                    and (NM is None or (np.isfinite(ratio) and ratio >= W.SCRAMBLE_RATIO)))
        rows.append(dict(tpl=bi, tid=tid, n_hits=int(len(g)), n_indep=int(len(gc)),
                         n_days=int(gc.date.nunique()), frac_near_seed=near, n_scramble=nm,
                         scramble_ratio=ratio, survives=surv, n_indep_quiet=n_quiet))
    R = (pd.DataFrame(rows, columns=W.SURV_COLS + ["n_indep_quiet"]) if rows
         else pd.DataFrame(columns=W.SURV_COLS + ["n_indep_quiet"]))
    R.to_csv(out + "_survival.csv", index=False)
    if indep:
        pd.concat(indep, ignore_index=True).to_parquet(out + "_indep.parquet", index=False)
    else:
        pd.DataFrame(columns=["tpl", "date", "t", "net_mad", "n_site_z3", "n_sites",
                              "tid"]).to_parquet(out + "_indep.parquet", index=False)
    W.log(f"[survive] {len(R)} templates -> {int(R.survives.sum()) if len(R) else 0} survivors")
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["handoff", "scramble", "quiet", "grow", "survive"])
    ap.add_argument("--clusters", default="")
    ap.add_argument("--seeds", default="")
    ap.add_argument("--bank", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--thr", default="")
    ap.add_argument("--scramble-hits", default="")
    ap.add_argument("--scramble-bank", default="")
    ap.add_argument("--shift", default="")
    ap.add_argument("--block", type=int, default=32)
    a = ap.parse_args()
    shift = json.loads(a.shift) if a.shift else None
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    if a.cmd in ("handoff", "scramble"):
        if os.path.exists(a.out):
            W.log(f"[handoff] {a.out} exists -- skip")
            return
        M = pd.read_parquet(a.clusters)
        Sd = W.read_csv_safe(a.seeds)
        if not len(Sd):
            import v_common as V
            W.patch_v(shift=shift)
            V.save_bank(a.out, [], W.sites()[3], np.zeros((0, len(W.sites()[3]), 3, 1), np.float32),
                        np.zeros((0, len(W.sites()[3]))), dict(kind="empty"))
            W.log("[handoff] no seed clusters -- empty bank")
            return
        handoff(M, Sd, a.out, shift=shift, scramble=(a.cmd == "scramble"),
                tremor_times=tremor_window_times() if a.cmd == "scramble" else None)
    elif a.cmd == "quiet":
        quiet_T(a.bank, a.out, shift=shift, block=a.block)
    elif a.cmd == "grow":
        grow(a.bank, a.out, a.thr, shift=shift, block=a.block)
    else:
        survive(a.bank, a.out, a.seeds, a.scramble_hits, a.scramble_bank)


if __name__ == "__main__":
    main()
