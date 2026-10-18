#!/usr/bin/env python
"""P3D transfer-stack engine (Merlin Q1).  Used by Step 1 (co-located machinery controls) and
Step 2 (transfer stacks for every new 3-C station + LRIV-HH).

One (family, station) leg =
  events   hp (p_real >= 0.8) of that family on the SSE-window days of the given episodes;
           N >= 30, else the labelled fallback (cat_core & p_real >= 0.5, N >= 60), else insufficient
  place    Crosson 1-D S at the family hypocentre, a5 correction (0 for a new station), +-20 s window
  stack    e5r_templates.cut_windows (rms-normalised) + align_stack (2 iterations, +-1.5 s)
  pick     S on the 0.5-s-smoothed 3-C envelope of the ALIGNED stack, +-4 s of nominal
  template 12 s from S-2 s, each component L2-normalised, frozen npz format

ACCEPTANCE (existence, not identity; alignment-free; every freedom charged to the null)
  (a) S-window excess 10log10 E(S-1..S+5)/E(S-12..S-6) of the UNALIGNED stack, maximised over a
      +-1.5 s constant offset, vs 200 matched RANDOM-TIME stacks put through the identical offset
      search: true > p99 AND true >= +3 dB
  (b) 100 WRONG-CELL moveout nulls (predicted S differing by >= 5 s, drawn from the P2 Slab2 grid):
      true > p95 of the shuffled
  (c) matched-N SAME-EVENTS control at the nearest existing template station (PGC/VGZ/B926): if the
      control itself fails (a), the leg is UNTESTABLE, not FAIL
  (d) self-consistency: an alignment-free re-stack at the template's own placement correlates
      >= 0.7 with the template

DIAGNOSTICS: single-event leave-one-out registration rho (median LOO cc of each event window against
its own stack over S-1..S+7), the PGC-dominance half-split control, and the N-matched S-SNR that
feeds the weight rule w_new = w_ref * 10^((SNR_new - SNR_ref)/20).
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                        # noqa: E402

HALF_W = 20.0                     # +-20 s cut window (frozen)
WIDE = 50.0                       # one wide cut per event; every offset is then a numpy slice
NW = int(round(2 * HALF_W * S.FS))
NOM = int(round(HALF_W * S.FS))   # nominal S sample inside the window
TPL_PRE = 2.0                     # 12-s template starts at S - 2 s (Merlin Q1)
TPL_LEN = 12.0
NTPL = int(round(TPL_LEN * S.FS))
OFF_MAX = int(round(1.5 * S.FS))  # +-1.5 s constant-offset search (charged to every null)
N_RAND = 200
N_CELL = 100
RNG_SEED = 20260919


NWIDE = int(round(2 * WIDE * S.FS))
O0 = int(round((WIDE - HALF_W) * S.FS))       # offset of the nominal +-20 s window inside the wide cut


def win_at(Wd, ds_samples, norm=True):
    """The +-20 s window at a time offset of `ds_samples` inside the wide cut, rms-renormalised
    per event exactly as the frozen cut_windows does."""
    a = O0 + int(ds_samples)
    if a < 0 or a + NW > Wd.shape[2]:
        return None
    seg = Wd[:, :, a:a + NW].astype(np.float64)
    if norm:
        r = np.sqrt((seg ** 2).mean(axis=(1, 2)))
        ok = np.isfinite(r) & (r > 0)
        if ok.sum() < 5:
            return None
        seg = seg[ok] / r[ok][:, None, None]
    return seg


def _stack(dc, sta, starts):
    W, keep = S.cut_windows_abs(dc, sta, starts, NW)
    if W is None or len(W) < 5:
        return None
    return W, keep


def excess_offsearch(st, nom=NOM, off=OFF_MAX):
    """max over a +-off-sample constant offset of the alignment-free S-window excess (dB)."""
    v = [S.s_excess_db(st, nom + o) for o in range(-off, off + 1, 2)]
    v = [x for x in v if np.isfinite(x)]
    return (max(v) if v else np.nan)


def excess_at_offsets(Wd, offsets_s):
    """Alignment-free S excess (with the +-1.5 s search) of the UNALIGNED stack at each time offset,
    all from one wide cut.  `offsets_s` in seconds."""
    out = []
    for d in np.atleast_1d(offsets_s):
        seg = win_at(Wd, int(round(float(d) * S.FS)))
        out.append(np.nan if seg is None else excess_offsearch(seg.mean(0)))
    return np.asarray(out, float)


def loo_registration(W, st, s_idx, pre=1.0, post=7.0):
    """median leave-one-out per-event cc of the individual windows against their own stack."""
    a = max(s_idx - int(round(pre * S.FS)), 0)
    b = min(s_idx + int(round(post * S.FS)), W.shape[2])
    x = W[:, :, a:b].reshape(len(W), -1)
    x = x - x.mean(1, keepdims=True)
    N = len(x)
    if N < 3:
        return np.nan
    tot = x.sum(0)
    ccs = []
    for i in range(N):
        r = (tot - x[i]) / (N - 1)
        r = r - r.mean()
        d = np.linalg.norm(x[i]) * np.linalg.norm(r)
        ccs.append(float(x[i] @ r / d) if d > 0 else np.nan)
    return float(np.nanmedian(ccs))


def cc3(a, b, lag=0):
    """3-C correlation coefficient of two (3, n) arrays, max over +-lag samples."""
    n = min(a.shape[1], b.shape[1])
    a = a[:, :n] - a[:, :n].mean(1, keepdims=True)
    best = -2.0
    for k in range(-lag, lag + 1):
        bb = np.roll(b[:, :n], k, axis=1)
        bb = bb - bb.mean(1, keepdims=True)
        d = np.linalg.norm(a) * np.linalg.norm(bb)
        if d > 0:
            best = max(best, float((a * bb).sum() / d))
    return best


def wrong_cell_offsets(bkey, sta, ST, rng, n=N_CELL, min_ds=5.0, max_ds=25.0):
    """Offsets (s) of the predicted S at `sta` for cells of the P2 Slab2 grid whose predicted S
    differs from the family's own by >= min_ds; the moveout null of Merlin (b)."""
    import p_vel1d as V
    F = S.families()
    key = (bkey, sta)
    if key not in wrong_cell_offsets._cache_tt:
        wrong_cell_offsets._cache_tt[key] = float(
            V.tt(float(S.epi_km(F.lat[bkey], F.lon[bkey], ST.lat[sta], ST.lon[sta])),
                 float(F.dep[bkey]), "S"))
    t_true = wrong_cell_offsets._cache_tt[key]
    if sta not in wrong_cell_offsets._cache_grid:
        import p_common as P
        g, _, _, _ = P.build_grid()
        d = S.epi_km(g.lat.values, g.lon.values, ST.lat[sta], ST.lon[sta])
        tt = np.array([V.tt(float(dd), float(zz), "S") for dd, zz in zip(d, g.dep.values)])
        wrong_cell_offsets._cache_grid[sta] = tt
    ds = wrong_cell_offsets._cache_grid[sta] - t_true
    ok = ds[(np.abs(ds) >= min_ds) & (np.abs(ds) <= max_ds)]
    if len(ok) == 0:
        ok = np.concatenate([rng.uniform(min_ds, max_ds, 200), -rng.uniform(min_ds, max_ds, 200)])
    return rng.choice(ok, size=min(n, len(ok)), replace=len(ok) < n)


wrong_cell_offsets._cache_tt = {}
wrong_cell_offsets._cache_grid = {}


def random_time_starts(dc, sta, days, n, rng):
    """n random window starts on the SAME station-days the real events use (matched-time null)."""
    out = []
    for _ in range(n):
        y, doy = days[rng.integers(len(days))]
        v = dc.get(sta, y, doy)
        if v is None:
            continue
        t0, X = v
        lim = X.shape[1] - NW - 1
        if lim <= 0:
            continue
        out.append(t0 + rng.integers(0, lim) / S.FS)
    return np.asarray(out, float)


def sp_seconds(bkey, sta, ST):
    """Crosson S-P at this (family, station) -- the P arrival is what contaminates the moveout null."""
    import p_vel1d as V
    F = S.families()
    d = float(S.epi_km(F.lat[bkey], F.lon[bkey], ST.lat[sta], ST.lon[sta]))
    z = float(F.dep[bkey])
    return float(V.tt(d, z, "S") - V.tt(d, z, "P"))


def leg(dc, sta, bkey, t_events, ttS, ST, rng, want_nulls=True, n_rand=N_RAND,
        n_cell=N_CELL, ref_template=None, keep_stack=True):
    """Build one (family, station) transfer leg and score it.  Returns a dict (never raises)."""
    out = dict(sta=sta, bkey=int(bkey), n_events=int(len(t_events)))
    t_events = np.asarray(t_events, float)
    Wd, keep = S.cut_fast(dc, sta, t_events + ttS - WIDE, NWIDE, norm=False)
    if Wd is None or len(Wd) < 5:
        out["status"] = "no_data"
        return out
    W = win_at(Wd, 0)
    if W is None:
        out["status"] = "no_data"
        return out
    out["n_stacked"] = int(len(W))
    st_un = W.mean(0)
    st_al, W_al, shifts = S.align_stack(W, iters=2, lag=OFF_MAX)
    s_idx = S.pick_s(st_al, NOM, search_s=4.0)
    out["s_pick_s"] = (s_idx - NOM) / S.FS
    out["snr_db"] = S.s_snr_db(st_al, s_idx)
    out["excess_db"] = excess_offsearch(st_un)
    out["rho_loo"] = loo_registration(W_al, st_al, s_idx)
    out["align_mad_s"] = float(1.4826 * np.median(np.abs(shifts - np.median(shifts))) / S.FS)

    # 12-s template from S - 2 s, each component L2-normalised (frozen format)
    a = s_idx - int(round(TPL_PRE * S.FS))
    tpl = np.zeros((3, NTPL))
    if 0 <= a and a + NTPL <= st_al.shape[1]:
        t = st_al[:, a:a + NTPL].copy()
        nrm = np.linalg.norm(t, axis=1, keepdims=True)
        tpl = t / np.where(nrm > 0, nrm, 1.0)
    out["tpl"] = tpl
    if keep_stack:
        out["stack_pm20"] = st_al.astype(np.float32)
        out["stack_pm20_unaligned"] = st_un.astype(np.float32)

    # (d) self-consistency: alignment-free re-stack at the template's OWN placement
    a2 = O0 + s_idx - int(round(TPL_PRE * S.FS))
    if 0 <= a2 and a2 + NTPL <= Wd.shape[2]:
        seg = Wd[:, :, a2:a2 + NTPL].astype(np.float64)
        r = np.sqrt((seg ** 2).mean(axis=(1, 2)))
        ok = np.isfinite(r) & (r > 0)
        out["cc_selfconsist"] = cc3(tpl, (seg[ok] / r[ok][:, None, None]).mean(0), lag=0)
    else:
        out["cc_selfconsist"] = np.nan
    out["pass_d"] = bool(np.isfinite(out["cc_selfconsist"]) and out["cc_selfconsist"] >= 0.7)

    if not want_nulls:
        out["status"] = "stacked"
        return out

    # (a) matched RANDOM-TIME stacks on the same station-days, identical +-1.5 s offset search
    days = sorted(set(zip(*S.yd_of(t_events[keep]))))
    N = len(W)
    # all n_rand x N random windows are cut in ONE day-sorted pass, then split into n_rand stacks
    rs = random_time_starts(dc, sta, days, n_rand * N, rng)
    nulls = []
    if len(rs):
        Wn, kn = S.cut_fast(dc, sta, rs, NW)
        if Wn is not None:
            grp = (kn // N)
            for g in range(n_rand):
                m = grp == g
                if m.sum() < max(5, N // 2):
                    continue
                nulls.append(excess_offsearch(Wn[m].astype(np.float64).mean(0)))
    nulls = np.asarray([x for x in nulls if np.isfinite(x)], float)
    out["n_null_rand"] = int(len(nulls))
    out["null_rand_p99"] = float(np.percentile(nulls, 99)) if len(nulls) >= 20 else np.nan
    out["null_rand_med"] = float(np.median(nulls)) if len(nulls) else np.nan
    out["pass_a"] = bool(np.isfinite(out["excess_db"]) and np.isfinite(out["null_rand_p99"])
                         and out["excess_db"] > out["null_rand_p99"]
                         and out["excess_db"] >= out["null_rand_med"] + 3.0)

    # (b) wrong-cell moveout nulls, all sliced out of the SAME wide cut
    offs = wrong_cell_offsets(int(bkey), sta, ST, rng, n=n_cell, max_ds=WIDE - HALF_W - 2.0)
    cn = excess_at_offsets(Wd, offs)
    fin = np.isfinite(cn)
    out["n_null_cell"] = int(fin.sum())
    out["null_cell_p95"] = float(np.percentile(cn[fin], 95)) if fin.sum() >= 20 else np.nan
    out["pass_b"] = bool(np.isfinite(out["excess_db"]) and np.isfinite(out["null_cell_p95"])
                         and out["excess_db"] > out["null_cell_p95"])
    # (b') the SAME null with the P-contaminated offsets removed.  A wrong cell whose predicted S
    # sits (S-P) s EARLY puts the family's own P arrival inside the S window, so the "wrong" cell
    # still scores a large excess: at a SINGLE station the moveout null degenerates.  Reported as a
    # DIAGNOSTIC alongside the literal (b), never instead of it.
    sp = sp_seconds(int(bkey), sta, ST)
    okm = fin & (np.abs(np.asarray(offs, float) + sp) > 4.0)
    out["sp_s"] = sp
    out["n_null_cell_noP"] = int(okm.sum())
    out["null_cell_noP_p95"] = float(np.percentile(cn[okm], 95)) if okm.sum() >= 20 else np.nan
    out["pass_b_noP"] = bool(np.isfinite(out["excess_db"]) and np.isfinite(out["null_cell_noP_p95"])
                             and out["excess_db"] > out["null_cell_noP_p95"])
    out["status"] = "scored"
    if ref_template is not None:
        out["cc_vs_ref_template"] = cc3(tpl, ref_template, lag=int(round(1.0 * S.FS)))
    return out
