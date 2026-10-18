#!/usr/bin/env python
"""Stage 6 -- LOCATION (g_locate + four controls) and GATES 0-2 (g_gates.evaluate).

Location (Merlin Q4, reduced to what this 6-site bed can support):
  per site, the independent hits are stacked (align_stack) at the template's own per-station
  offset; S is picked on the 3-C envelope within +- 3 s of the template onset and accepted at
  peak/noise-rms >= 3; the site S times go into `g_locate.locate_interface` on the 2-D Slab2 grid
  (centre 48.50 / -123.75, half-width 62 km), with the free-depth 3-D run as the sanity check.
  Controls: 200 shuffled-pick nulls (rms ratio < 0.5), 50-realisation synthetic floor (<= 3 km),
  leave-one-site-out (<= 5 km), free depth within +- 10 km of Slab2 and not at a bound.

The located cell converts the cluster time axis to the ORIGIN axis
  t_origin = t_cluster + median_sites(tS_site - ttS(cell -> site) - corr_site)
so the proposal handed to `g_gates` is on the axis the gates expect; `offset_search_s = 30` is
still passed (the tested value) and `offset_fit_s` is read back.

Gate 3 is NOT run: it is UNTESTABLE inside one ETS window (Merlin F6).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402
import s_common as S                                                        # noqa: E402


def _pick(st, nom, search_s=W.LOC_PICK_SEARCH_S):
    """(index, peak/noise-rms ratio) of the S onset on the 3-C envelope."""
    env = S.envelope(st)
    w = int(round(search_s * W.FS))
    a, b = max(nom - w, 0), min(nom + w, len(env))
    i = int(a + np.argmax(env[a:b]))
    n0 = max(i - int(round(14 * W.FS)), 0)
    n1 = max(i - int(round(6 * W.FS)), 1)
    noi = env[n0:n1]
    r = float(env[i] / np.sqrt((noi ** 2).mean())) if len(noi) and (noi ** 2).mean() > 0 else np.nan
    return i, r


def site_picks(loader, tids, stas, Wt, off, meta, H, bi, tid, mem, of, order):
    """{site: (tS relative to the cluster axis, snr, n_stacked)} for one template."""
    import w_grow as GR
    pre = float(meta.get("pre_s", W.TPL1_PRE_S))
    g = H[H.tid == tid].sort_values("net_mad", ascending=False).head(W.MAX_STACK)
    t = g.t.to_numpy(float)
    n_cut = int(round(W.HANDOFF_CUT_LEN_S * W.FS))
    nom = int(round(W.HANDOFF_CUT_PRE_S * W.FS))
    per = {}
    for si, s in enumerate(stas):
        o = off[bi, si]
        if not np.isfinite(o):
            continue
        onset = o + pre                       # onset time relative to the cluster axis
        segs = [w for w in (GR._cut(loader, s, tt + onset - W.HANDOFF_CUT_PRE_S, n_cut)
                            for tt in t) if w is not None]
        if len(segs) < 10:
            continue
        st, _, _ = S.align_stack(np.stack(segs), iters=W.ALIGN_ITERS, lag=W.ALIGN_LAG)
        i_s, snr = _pick(st, nom)
        per[s] = (onset + (i_s - nom) / W.FS, snr, len(segs))
    out = {}
    for q in order:
        v = [per[s] for s in mem[q] if s in per and np.isfinite(per[s][1])]
        if not v:
            continue
        out[q] = (float(np.median([x[0] for x in v])), float(np.median([x[1] for x in v])),
                  int(np.median([x[2] for x in v])))
    return out


def locate_one(grid, site_ll, picks, order, rng):
    """locate_interface + the four controls.  Returns a dict of statistics."""
    import g_locate as GL
    sel_sites = [q for q in order if q in picks and picks[q][1] >= W.LOC_MIN_PICK_SNR]
    if len(sel_sites) < 3:
        return dict(n_picked=len(sel_sites), located=False, reason="picks<3")
    tS = np.full(len(grid.sites), np.nan)
    tP = np.full(len(grid.sites), np.nan)
    for q in sel_sites:
        tS[grid.sites.index(q)] = picks[q][0]
    sel = np.array([grid.sites.index(q) for q in sel_sites], int)
    r = GL.locate_interface(grid, tS, tP, sel)
    if r is None:
        return dict(n_picked=len(sel_sites), located=False, reason="no_solution")
    rf = GL.locate_free(grid, tS, tP, sel)
    # shuffled-pick null
    null = []
    for _ in range(W.LOC_N_PICK_NULL):
        tt = tS.copy()
        v = tt[sel]
        tt[sel] = rng.permutation(v)
        q = GL.locate_interface(grid, tt, tP, sel)
        if q is not None:
            null.append(q["rms_s"])
    ratio = float(r["rms_s"] / np.median(null)) if null else np.nan
    # synthetic floor
    med, p90, _ = GL.synthetic_floor(grid, site_ll, sel_sites, r["lat"], r["lon"], r["z"],
                                     [False] * len(sel_sites), n_real=W.LOC_N_SYNTH,
                                     seed=int(rng.integers(1 << 30)))
    # leave-one-site-out
    loso = []
    for k in range(len(sel)):
        s2 = np.delete(sel, k)
        if len(s2) < 3:
            continue
        q = GL.locate_interface(grid, tS, tP, s2)
        if q is not None:
            loso.append(GL.G.epi_km(q["lat"], q["lon"], r["lat"], r["lon"]))
    gap = GL.azimuth_gap(r["lat"], r["lon"], sel_sites, site_ll)
    d60 = sum(1 for q in sel_sites if GL.G.epi_km(r["lat"], r["lon"], *site_ll[q]) <= 60.0)
    dz = abs(rf["z"] - r["z"]) if rf else np.nan
    located = bool(d60 >= W.LOC_MIN_SITES and gap <= W.LOC_MAX_GAP_DEG
                   and np.isfinite(ratio) and ratio < W.LOC_NULL_RATIO
                   and np.isfinite(med) and med <= W.LOC_SYNTH_BAR_KM
                   and (not loso or max(loso) <= W.LOC_LOSO_BAR_KM)
                   and np.isfinite(dz) and dz <= W.LOC_DEPTH_BAR_KM
                   and not r["at_bound"])
    return dict(lat=r["lat"], lon=r["lon"], dep=r["z"], rms_s=r["rms_s"], at_bound=r["at_bound"],
                n_picked=len(sel_sites), n_sites_60km=d60, gap_deg=gap,
                shuffled_ratio=ratio, synth_med_km=med, synth_p90_km=p90,
                loso_max_km=float(max(loso)) if loso else np.nan,
                free_lat=rf["lat"] if rf else np.nan, free_lon=rf["lon"] if rf else np.nan,
                free_dep=rf["z"] if rf else np.nan, depth_diff_km=dz,
                sites=",".join(sel_sites), located=located)


def run(bank, indep, dedupe_csv, survival_csv, out, gates=True, shift=None):
    V, E = W.patch_v(shift=shift)
    import g_common as G
    import g_locate as GL
    mem, of, order, det = W.sites()
    tids, stas, Wt, off, meta = V.load_bank(bank)
    H = (pd.read_parquet(indep) if os.path.exists(indep)
         else pd.DataFrame(columns=W.HIT_COLS + ["tid"]))
    keep = list(tids)
    D = W.read_csv_safe(dedupe_csv, W.DEDUPE_COLS)
    if len(D):
        keep = [t for t in tids if t in set(D[D.kept.astype(bool)].tid)]
    SV = W.read_csv_safe(survival_csv, W.SURV_COLS)
    if len(SV):
        keep = [t for t in keep if t in set(SV[SV.survives.astype(bool)].tid)]
    W.log(f"[locate] {len(keep)}/{len(tids)} templates enter location")
    if not keep:
        pd.DataFrame(columns=["tid", "n_picked", "located", "reason"]).to_csv(
            out + "_locations.csv", index=False)
        pd.DataFrame(columns=["pid", "verdict"]).to_csv(out + "_gates.csv", index=False)
        pd.DataFrame(columns=["cid", "t", "lat", "lon"]).to_parquet(out + "_clusters.parquet",
                                                                    index=False)
        json.dump(dict(n_templates=0, n_located=0, n_gate_pass=0), open(out + ".json", "w"),
                  indent=1)
        return
    ST = W.station_table()
    site_ll = {q: (float(ST.lat[mem[q][0]]), float(ST.lon[mem[q][0]])) for q in order}
    grid = GL.Grid(order, site_ll, W.GRID_LAT0, W.GRID_LON0, halfwidth=W.GRID_HALFWIDTH)
    rng = np.random.default_rng(W.SITESHIFT_SEED)
    loader = E._W_LOADER
    rows, clus = [], []
    for tid in keep:
        bi = tids.index(tid)
        try:
            pk = site_picks(loader, tids, stas, Wt, off, meta, H, bi, tid, mem, of, order)
            r = locate_one(grid, site_ll, pk, order, rng)
        except Exception as e:                                              # noqa: BLE001
            W.log(f"[locate] {tid} ERROR {e}")
            pk, r = {}, dict(n_picked=0, located=False, reason=f"error:{str(e)[:120]}")
        r["tid"] = tid
        g = H[H.tid == tid]
        r["n_members"] = int(len(g))
        for q in order:
            r[f"tS_{q}"] = pk[q][0] if q in pk else np.nan
            r[f"snr_{q}"] = pk[q][1] if q in pk else np.nan
        if r.get("located"):
            tt = V.cell_tt(r["lat"], r["lon"], r["dep"], stas=det, ST=ST)
            dd = []
            for q in order:
                if q in pk and pk[q][1] >= W.LOC_MIN_PICK_SNR:
                    dd.append(pk[q][0] - float(np.median([tt[s] for s in mem[q]])))
            r["axis_offset_s"] = float(np.median(dd)) if dd else np.nan
            torg = g.t.to_numpy(float) + r["axis_offset_s"]
            clus.append(pd.DataFrame(dict(cid=tid, t=torg, lat=r["lat"], lon=r["lon"])))
        rows.append(r)
    L = pd.DataFrame(rows)
    L.to_csv(out + "_locations.csv", index=False)
    C = (pd.concat(clus, ignore_index=True) if clus
         else pd.DataFrame(columns=["cid", "t", "lat", "lon"]))
    C.to_parquet(out + "_clusters.parquet", index=False)
    W.log(f"[locate] {int(L.located.sum()) if 'located' in L else 0} located")
    ver = []
    if gates and len(C):
        import g_gates as GG
        ctx = GG.GateContext(GG.GateConfig(withheld=tuple(W.WITHHELD)), verbose=False)
        days = np.array([np.datetime64(d, "D") for d in W.availability().date])
        for cid, g in C.groupby("cid"):
            p = GG.Proposal(pid=str(cid), times=g.t.to_numpy(float), lat=float(g.lat.iloc[0]),
                            lon=float(g.lon.iloc[0]), dep=np.nan,
                            offset_search_s=W.GATES_OFFSET_SEARCH_S, meta=dict(drop_bkey=-1))
            try:
                res = GG.evaluate(ctx, p, gates=W.GATES, evaluated_days=days)
                ver.append(GG.flatten(res))
            except Exception as e:                                          # noqa: BLE001
                ver.append(dict(pid=str(cid), verdict="ERROR", error=str(e)[:200]))
            W.log(f"[gates] {cid} -> {ver[-1].get('verdict')}")
    Gv = pd.DataFrame(ver) if ver else pd.DataFrame(columns=["pid", "verdict"])
    Gv.to_csv(out + "_gates.csv", index=False)
    npass = int((Gv.verdict == "PASS").sum()) if len(Gv) else 0
    json.dump(dict(n_templates=len(keep), n_located=int(L.located.sum()) if "located" in L else 0,
                   n_gate_pass=npass), open(out + ".json", "w"), indent=1)
    W.log(f"[gates] {npass} PASS of {len(Gv)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bank", required=True)
    ap.add_argument("--indep", required=True)
    ap.add_argument("--dedupe", default="")
    ap.add_argument("--survival", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-gates", action="store_true")
    ap.add_argument("--shift", default="")
    a = ap.parse_args()
    shift = json.loads(a.shift) if a.shift else None
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    run(a.bank, a.indep, a.dedupe, a.survival, a.out, gates=not a.no_gates, shift=shift)


if __name__ == "__main__":
    main()
