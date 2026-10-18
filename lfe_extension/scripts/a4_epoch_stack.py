#!/usr/bin/env python
"""Track A-2 STEP 4 -- epoch transfer and IDENTITY test (Merlin plan-review step 4).

For one (station, epoch boundary):
  * reference events = the detections of `a4_run_selection.sh` for that boundary, i.e. selected by the
    network with THAT STATION EXCLUDED from the sum, at the frozen T_core (8.75 MAD), with the frozen
    two-site coincidence gate, on the tremor days of the SSE window;
  * the station's 3-C 20 Hz / 1-8 Hz data is cut at the PREDICTED S arrival (t_detection + ttS from the
    Crosson 1-D model to the family hypocentre), rms-normalised, cc-aligned (ONE shift per station for
    all three components) and linearly stacked -- e5r_templates.cut_windows / align_stack, unchanged;
  * the new stack is compared with the PREVIOUS-EPOCH (2012, from-scratch) stack of the same
    (family, station) over S-1..S+7 s: signed 3-C cc maximised over +-2 s, the lag at the extremum
    (= the S-pick shift dS), and the per-event registration (median LEAVE-ONE-OUT cc of the individual
    event windows against the stack, over the same 8 s).

The same script is run on the *_old pass (2012 tremor days, same withheld station) -- the CONTROL that
says what this procedure costs by itself, before any instrument change.  The verdict uses the
pre-registered absolute bars (|cc| >= 0.6, |dS| <= 0.5 s, registration >= 0.8x the previous epoch) with
the previous-epoch value taken from that control, so the two sides are measured identically.

Outputs (per station): appends to results/a4_epoch_identity.csv and writes
results/a4_templates_<STA>_<epoch>.npz (accepted families only, --extra-tpl format for Track C).
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e5r_templates as T                                      # noqa: E402
import d3_vel1d as V                                           # noqa: E402

REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
EXT = os.path.dirname(HERE)
FS, HALF_W, PRE_TPL, PRE_TT = 20.0, 20.0, 1.0, 0.5
IDX_NOM = int(round(HALF_W * FS))
NCMP = 160                          # S-1 .. S+7 s at 20 Hz = the comparison window
LAG_MAX = int(round(2.0 * FS))      # +-2 s search for the identity cc


def cc3(A, B):
    """3-C zero-lag cc of two (3, n) blocks, components summed (one common normalisation)."""
    A = A - A.mean(1, keepdims=True); B = B - B.mean(1, keepdims=True)
    d = np.sqrt((A ** 2).sum() * (B ** 2).sum())
    return float((A * B).sum() / d) if d > 0 else np.nan


def best_cc(ref, new, i_ref, lag_max=LAG_MAX, n=NCMP):
    """Signed 3-C cc of `new` against `ref[:, i_ref:i_ref+n]`, maximised in |cc| over +-lag_max
    samples; returns (cc, lag_seconds)."""
    A = ref[:, i_ref:i_ref + n]
    if A.shape[1] < n:
        return np.nan, np.nan
    out = []
    for k in range(-lag_max, lag_max + 1):
        j = i_ref + k
        if j < 0 or j + n > new.shape[1]:
            continue
        out.append((cc3(A, new[:, j:j + n]), k / FS))
    if not out:
        return np.nan, np.nan
    v = max(out, key=lambda z: abs(z[0]) if np.isfinite(z[0]) else -1)
    return v[0], v[1]


def registration(Sal, st, i0, n=NCMP):
    """Median LEAVE-ONE-OUT per-event cc of the aligned windows against the stack over [i0, i0+n).
    LOO removes the trivial N-dependence of correlating an event with a stack containing it."""
    N = len(Sal)
    if N < 5 or i0 < 0 or i0 + n > st.shape[1]:
        return np.nan
    tot = st * N
    ccs = []
    for i in range(N):
        loo = (tot[:, i0:i0 + n] - Sal[i][:, i0:i0 + n]) / (N - 1)
        ccs.append(cc3(loo, Sal[i][:, i0:i0 + n]))
    return float(np.nanmedian(ccs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sta", required=True)
    ap.add_argument("--netsta", required=True)
    ap.add_argument("--det", required=True, help="e5r_detections_<tag>.parquet of the withheld-station pass")
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--epoch", required=True, help="label, e.g. 2018_HH or 2012_BH")
    ap.add_argument("--boundary", required=True)
    ap.add_argument("--leg", default="new", choices=["new", "old"])
    ap.add_argument("--tmad", type=float, default=8.75)
    ap.add_argument("--min-events", type=int, default=20)
    ap.add_argument("--max-events", type=int, default=400)
    ap.add_argument("--out", default="")
    ap.add_argument("--match-n", default="", help="a4_epoch_identity_<key>_new.csv: cap each family's "
                    "event count at the number the NEW leg used, so cc is compared at matched N "
                    "(cc between two stacks rises with N on either side)")
    a = ap.parse_args()

    T.CACHES = [f"{EXT}/cache20hz", f"{REP}/cache20hz"]
    T.YEAR = a.year
    stab = pd.read_csv(f"{REP}/results/e5r_stations.csv").set_index("sta")
    lat_s, lon_s = float(stab.loc[a.sta, "lat"]), float(stab.loc[a.sta, "lon"])
    F = pd.read_csv(f"{REP}/results/e5r_families_2012.csv").set_index("bkey")
    Z = np.load(f"{REP}/results/e5r_templates_2012.npz", allow_pickle=True)
    QC = pd.read_csv(f"{REP}/results/e5r_template_qc_2012.csv")
    QC = QC[(QC.sta == a.sta) & QC.usable & QC.template].set_index("bkey")

    D = pd.read_parquet(a.det)
    D = D[(D.net_mad >= a.tmad) & (D.L == 12.0)]
    NMATCH = {}
    if a.match_n:
        m = pd.read_csv(a.match_n)
        NMATCH = {int(r.family): int(r.n_events) for r in m[m.status == "ok"].itertuples()}
        print(f"[{a.sta}/{a.leg}] matched-N control: caps from {os.path.basename(a.match_n)} for "
              f"{len(NMATCH)} families (median {int(np.median(list(NMATCH.values()))) if NMATCH else 0})",
              flush=True)
    print(f"[{a.sta}/{a.leg}] {len(D)} detections >= {a.tmad} MAD over {D.fam.nunique()} families, "
          f"{D.doy.nunique()} days, year {a.year}", flush=True)

    days = T.load_days(a.netsta, a.year)
    if not days:
        raise SystemExit(f"no 20 Hz cache for {a.netsta} {a.year}")
    print(f"[{a.sta}/{a.leg}] {len(days)} cached days", flush=True)
    rng = np.random.default_rng(7)
    nwin = int(round(2 * HALF_W * FS))
    rows, TPL, SP, STK = [], {}, {}, {}
    for bkey in sorted(QC.index):
        if bkey not in F.index:
            continue
        g = D[D.fam == float(bkey)]
        if len(g) < a.min_events:
            rows.append(dict(sta=a.sta, leg=a.leg, epoch=a.epoch, boundary=a.boundary, family=int(bkey),
                             n_events=len(g), status="too_few_detections"))
            continue
        cap = min(a.max_events, NMATCH.get(int(bkey), a.max_events)) if NMATCH else a.max_events
        if NMATCH and int(bkey) not in NMATCH:
            continue                      # the new leg never stacked this family: nothing to match
        g = g.sort_values("net_mad", ascending=False).head(cap)
        fr = F.loc[bkey]
        dx = T.epi_km(float(fr.lat), float(fr.lon), lat_s, lon_s)
        ttP, ttS = float(V.tt(dx, float(fr.dep), "P")), float(V.tt(dx, float(fr.dep), "S"))
        S, keep = T.cut_windows(days, g.t.to_numpy() + ttS - HALF_W, nwin)
        if S is None or len(S) < a.min_events:
            rows.append(dict(sta=a.sta, leg=a.leg, epoch=a.epoch, boundary=a.boundary, family=int(bkey),
                             n_events=0 if S is None else len(S), status="no_data"))
            continue
        st, Sal, sh = T.align_stack(S, iters=2)
        sp_idx = T.pick_s(st, IDX_NOM)
        snr = T.s_snr_db(st, sp_idx)
        tmin, tmax = float(g.t.min()), float(g.t.max())
        Sn, _ = T.cut_windows(days, rng.uniform(tmin, tmax, len(S)) + ttS - HALF_W, nwin)
        snr_n = np.nan
        if Sn is not None and len(Sn) >= 10:
            stn, _, _ = T.align_stack(Sn, iters=2)
            snr_n = T.s_snr_db(stn, T.pick_s(stn, IDX_NOM))

        ref = np.asarray(Z[f"{a.sta}|{int(bkey)}"], float)      # the 2012 (previous-epoch) stack
        i_ref = IDX_NOM + int(round(float(QC.loc[bkey, "s_pick_s"]) * FS)) - int(round(PRE_TPL * FS))
        cc, dS = best_cc(ref, st, i_ref)
        reg = registration(Sal, st, i_ref + int(round(dS * FS)))
        # the 2012 stack's own registration cannot be recomputed here (its event windows were not
        # kept), so the ratio is formed against the *_old control leg in a4_epoch_report.py.
        d = {}
        i0 = sp_idx - int(round(PRE_TPL * FS))
        ok = i0 >= 0 and i0 + int(12 * FS) <= st.shape[1]
        if ok:
            seg = st[:, i0:i0 + int(12 * FS)].astype(np.float32)
            nn = np.linalg.norm(seg, axis=1)
            ok = bool(np.all(nn > 0))
            if ok:
                d[12.0] = (seg / nn[:, None]).astype(np.float32)
        if ok:
            key = (float(bkey), a.sta)
            TPL[key] = d
            SP[key] = (ttS - ttP) + (sp_idx - IDX_NOM) / FS - (PRE_TPL - PRE_TT)
            STK[f"{a.sta}|{int(bkey)}"] = st.astype(np.float32)
        rows.append(dict(sta=a.sta, leg=a.leg, epoch=a.epoch, boundary=a.boundary, family=int(bkey),
                         n_events=int(len(S)), n_det=int(len(g)), epi_km=dx,
                         cc=cc, dS=dS, registration=reg, snr_db=snr, snr_null_db=snr_n,
                         s_pick_s=(sp_idx - IDX_NOM) / FS, s_pick_2012=float(QC.loc[bkey, "s_pick_s"]),
                         align_shift_s=float(np.median(sh)) / FS, template=int(ok), status="ok"))
    R = pd.DataFrame(rows)
    out = a.out or f"{EXT}/results/a4_epoch_identity_{a.sta}_{a.leg}.csv"
    R.to_csv(out, index=False)
    np.savez_compressed(f"{EXT}/results/a4_stacks_{a.sta}_{a.epoch}_{a.leg}.npz",
                        tpl=np.array(TPL, dtype=object), sp=np.array(SP, dtype=object), **STK)
    o = R[R.status == "ok"]
    print(f"[{a.sta}/{a.leg}] {len(o)} families stacked; |cc| >= 0.6: {int((o.cc.abs() >= 0.6).sum())}; "
          f"cc median {o.cc.median():.3f}; |dS| <= 0.5 s: {int((o.dS.abs() <= 0.5).sum())}; "
          f"dS median {o.dS.median():+.2f} s; registration median {o.registration.median():.3f}", flush=True)
    print(f"[write] {out}", flush=True)


if __name__ == "__main__":
    main()
