#!/usr/bin/env python
"""P3D STEP 2 -- transfer stacks on the CALIBRATION episodes for every new 3-C station + LRIV-HH.

Protocol: Merlin Q1 (see s_transfer.py docstring).  Per (family, station):
  events   hp (p_real >= 0.8) on the SSE-window days of the station's own calibration episodes
           (S01-S13; BUCKS/TURTL/ODUC build on S18), N >= 30; labelled fallback (cat_core &
           p_real >= 0.5) to N >= 60; else `insufficient`
  score    (a) random-time existence null, (b) wrong-cell moveout null (+ the P-guarded diagnostic
           (b')), (c) matched-N SAME-EVENTS control at the nearest existing template station ->
           FAIL becomes UNTESTABLE, (d) self-consistency
  weights  w_new = w_ref * 10^((S-SNR_new - S-SNR_ref)/20) with BOTH S-SNRs measured on the SAME
           events at the SAME N (the reference leg is built on exactly the same event list), then
           the frozen 3x-median clip and median normalisation
  extras   single-event leave-one-out registration rho at the new site AND at PGC/VGZ on the same
           events (THE number the discovery arm hinges on); PGC-dominance half-split control

Outputs results/p3d/step2_legs.csv, step2_summary.csv, step2_rho.csv,
        bank_transfer_v1/{bank_full.npz, bank_production.npz, bank_*.csv}
"""
from __future__ import annotations
import os, sys, json, time, argparse, warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                        # noqa: E402
import s_transfer as T                                                      # noqa: E402

YEARS = list(range(2019, 2027))
REF_CHOICES = ["PGC", "VGZ", "B926"]
MAX_EV = 300


def episodes_for(sta):
    """Per-station CALIBRATION episodes (pre-registered split, incl. the PQ amendment):
    pre-2023 stations S01-S13; JRDRB/ALBHB S08-S13; OHOH S12-S13;
    BUCKS/TURTL/ODUC/LHLYB/DAOB/ODEER build on S18 alone."""
    return S.cal_episodes_for(sta)


def _one_station(args):
    sta, fams, max_ev, do_nulls, variant = args
    warnings.filterwarnings("ignore")
    ST = (S.station_table() if variant == "v1"
          else pd.read_csv(f"{S.P3D}/stations_p3d_v2.csv").set_index("sta"))
    F = S.families()
    eps = episodes_for(sta)
    H = S.hp_events(0.8, episodes=eps)
    HF = S.hp_events(0.5, episodes=eps, column="cat_core")
    stas = [sta] + REF_CHOICES
    tts = S.fam_tt(stas, ST)
    idx = S.cache_index(years=YEARS, stas=stas)
    have = {s: {(y, d) for (ss, y, d) in idx if ss == s} for s in stas}
    for DF in (H, HF):
        yy, dd = S.yd_of(DF.t.values)
        DF["ok_day"] = np.array([(y, d) in have[sta] for y, d in zip(yy, dd)])
    dc = S.DayCache(maxn=2200, years=YEARS, stas=stas)
    rng = np.random.default_rng(T.RNG_SEED + int(__import__("hashlib").sha256(sta.encode()).hexdigest()[:6], 16))
    # PGC-site z per detection for the dominance half-split control
    PS = persite_z(eps, "PGC")
    rows, tpls = [], {}
    for b in fams:
        # reference station = nearest of PGC/VGZ/B926 with a 2012 template for this family
        rec = dict(sta=sta, bkey=int(b), site=ST.site[sta])
        e = H[(H.bkey == b) & H.ok_day]
        kind = "hp"
        if len(e) < 30:
            e2 = HF[(HF.bkey == b) & HF.ok_day]
            if len(e2) >= 60:
                e, kind = e2, "fallback"
            else:
                rec.update(status="insufficient", n_events=int(max(len(e), len(e2))), kind="none")
                rows.append(rec)
                continue
        if len(e) > max_ev:
            e = e.sort_values("p_real", ascending=False).head(max_ev)
        e = e.sort_values("t")
        t_ev = e.t.to_numpy()
        rec.update(kind=kind, n_events=len(t_ev))
        L = T.leg(dc, sta, b, t_ev, tts[(b, sta)], ST, rng, want_nulls=do_nulls)
        for k, v in L.items():
            if not isinstance(v, np.ndarray):
                rec[k] = v
        if L.get("status") == "no_data":
            rows.append(rec)
            continue
        tpls[(int(b), sta)] = dict(tpl=L["tpl"], stack=L["stack_pm20"], s_pick_s=L["s_pick_s"],
                                   n=L.get("n_stacked", 0), kind=kind)
        # ---- reference legs on the SAME events (control (c), matched-N weights, rho comparison)
        refs = sorted(REF_CHOICES, key=lambda s: float(S.hyp_km(F.lat[b], F.lon[b], F.dep[b],
                                                                ST.lat[s], ST.lon[s])))
        for rs in refs:
            if not have.get(rs):
                continue
            RL = T.leg(dc, rs, b, t_ev, tts[(b, rs)], ST, rng,
                       want_nulls=(do_nulls and rs == refs[0]), keep_stack=False)
            if RL.get("status") == "no_data":
                continue
            rec[f"ref_{rs}_snr_db"] = RL.get("snr_db", np.nan)
            rec[f"ref_{rs}_rho"] = RL.get("rho_loo", np.nan)
            rec[f"ref_{rs}_excess_db"] = RL.get("excess_db", np.nan)
            rec[f"ref_{rs}_n"] = RL.get("n_stacked", 0)
            if rs == refs[0]:
                rec["ref_sta"] = rs
                rec["ref_pass_a"] = RL.get("pass_a", False)
                rec["ref_pass_b"] = RL.get("pass_b", False)
                rec["ref_pass_b_noP"] = RL.get("pass_b_noP", False)
                rec["ref_pass_d"] = RL.get("pass_d", False)
        # ---- PGC-dominance half-split control
        rec.update(half_split(dc, sta, b, e, PS, tts[(b, sta)], ST, rng))
        rows.append(rec)
    return rows, tpls


def persite_z(episodes, site):
    """{rounded detection time: site z} for the PGC site over the episodes' months."""
    W = S.sse_windows()
    W = W[W.id.isin(episodes)]
    months = sorted({m.strftime("%Y-%m") for r in W.itertuples()
                     for m in pd.date_range(r.start, r.end, freq="D")})
    out = {}
    for m in months:
        p = f"{S.EXT}/results/trackE/persite_{m}.parquet"
        if not os.path.exists(p):
            continue
        d = pd.read_parquet(p, columns=["fam", "t", "site", "z_site", "ctrl"])
        d = d[(d.site == site) & (d.ctrl == 0)]
        for fam, t, z in zip(d.fam.values, d.t.values, d.z_site.values):
            out[(int(fam), round(float(t), 2))] = float(z)
    return out


def half_split(dc, sta, b, e, PS, ttS, ST, rng):
    """Split the event set at the median PGC-site z and require the two new-station stacks to agree
    above the station's own cross-family chance level (Merlin Q1 PGC-dominance control)."""
    z = np.array([PS.get((int(b), round(float(t), 2)), np.nan) for t in e.t.values])
    ok = np.isfinite(z)
    if ok.sum() < 30:
        return dict(hs_status="no_z", hs_n=int(ok.sum()))
    tt = e.t.to_numpy()[ok]
    zz = z[ok]
    med = np.median(zz)
    lo, hi = tt[zz <= med], tt[zz > med]
    if min(len(lo), len(hi)) < 12:
        return dict(hs_status="too_few", hs_n=int(ok.sum()))
    Ls = [T.leg(dc, sta, b, x, ttS, ST, rng, want_nulls=False, keep_stack=True) for x in (lo, hi)]
    if any(L.get("status") == "no_data" for L in Ls):
        return dict(hs_status="no_data")
    cc = T.cc3(Ls[0]["tpl"], Ls[1]["tpl"], lag=int(round(1.0 * S.FS)))
    return dict(hs_status="ok", hs_n=int(ok.sum()), hs_n_lo=len(lo), hs_n_hi=len(hi),
                hs_cc=cc, hs_excess_lo=Ls[0]["excess_db"], hs_excess_hi=Ls[1]["excess_db"],
                hs_rho_lo=Ls[0]["rho_loo"], hs_rho_hi=Ls[1]["rho_loo"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", default="")
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--max-events", type=int, default=MAX_EV)
    ap.add_argument("--no-nulls", action="store_true")
    ap.add_argument("--tag", default="")
    ap.add_argument("--variant", default="v1", choices=["v1", "v2"])
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    ST = (S.station_table() if a.variant == "v1"
          else pd.read_csv(f"{S.P3D}/stations_p3d_v2.csv").set_index("sta"))
    targets = [s for s in ST.index if ST.transfer_target[s] and ST["mode"][s] == "3C"]
    if a.stations:
        targets = [s for s in a.stations.split(",") if s in targets]
    fams = sorted(int(b) for b in S.families().index)
    print(f"[step2] {len(targets)} stations x {len(fams)} families, workers={a.workers}", flush=True)
    jobs = [(s, fams, a.max_events, not a.no_nulls, a.variant) for s in targets]
    t0 = time.time()
    allrows, alltpl = [], {}
    import multiprocessing as mp
    with mp.Pool(min(a.workers, len(jobs))) as pool:
        for rows, tpls in pool.imap_unordered(_one_station, jobs):
            allrows += rows
            alltpl.update(tpls)
            done = {r["sta"] for r in rows}
            print(f"  done {done} ({len(rows)} legs)  {time.time()-t0:.0f}s", flush=True)
    D = pd.DataFrame(allrows)
    # ---- verdicts
    D["accept"] = (D.get("pass_a", False).fillna(False) & D.get("pass_b", False).fillna(False)
                   & D.get("pass_d", False).fillna(False))
    D["accept_noP"] = (D.get("pass_a", False).fillna(False)
                       & D.get("pass_b_noP", False).fillna(False)
                       & D.get("pass_d", False).fillna(False))
    ctl_fail = ~(D.get("ref_pass_a", False).fillna(False) & D.get("ref_pass_d", False).fillna(False))
    D["verdict"] = np.where(D.status.eq("insufficient"), "insufficient",
                   np.where(D.status.eq("no_data"), "no_data",
                   np.where(ctl_fail, "untestable",
                   np.where(D.accept_noP, "accept", "reject"))))
    D.to_csv(f"{S.P3D}/step2_legs{a.tag}.csv", index=False)

    # ---- N-matched weights
    W2012 = pd.read_csv(f"{S.IN}/e5r_weights_snr_2012.csv")
    wref = {(int(r.fam), r.sta): float(r.w) for r in W2012.itertuples()}
    w = []
    for r in D.itertuples():
        rs = getattr(r, "ref_sta", None)
        s_new = getattr(r, "snr_db", np.nan)
        s_ref = getattr(r, f"ref_{rs}_snr_db", np.nan) if isinstance(rs, str) else np.nan
        w0 = wref.get((int(r.bkey), rs), np.nan) if isinstance(rs, str) else np.nan
        w.append(w0 * 10 ** ((s_new - s_ref) / 20.0)
                 if np.isfinite(s_new) and np.isfinite(s_ref) and np.isfinite(w0) else np.nan)
    D["w_raw"] = w
    med = float(np.nanmedian(D.w_raw[D.verdict == "accept"])) if (D.verdict == "accept").any() else np.nan
    D["w_new"] = np.minimum(D.w_raw, 3.0 * med) / med
    D.to_csv(f"{S.P3D}/step2_legs{a.tag}.csv", index=False)

    # ---- the bank, in the frozen npz format; hidden families INCLUDED / EXCLUDED
    hid, iso, _ = S.holdouts()
    bd = f"{S.P3D}/bank_transfer_v1{a.tag}"   # --tag keeps a partial run from clobbering the main bank
    os.makedirs(bd, exist_ok=True)
    acc = D[D.verdict == "accept"]
    for tag, drop in (("full", set()), ("production", hid)):
        tpl, sp, extra, keys = {}, {}, {}, []
        for r in acc.itertuples():
            k = (int(r.bkey), r.sta)
            if k[0] in drop or k not in alltpl:
                continue
            v = alltpl[k]
            tpl[(float(k[0]), k[1])] = {12.0: np.asarray(v["tpl"], np.float32)}
            sp[(float(k[0]), k[1])] = float(v["s_pick_s"])
            extra[f"{k[1]}|{k[0]}"] = np.asarray(v["stack"], np.float32)
            keys.append(dict(bkey=k[0], sta=k[1], site=r.site, n=v["n"], kind=v["kind"],
                             w_new=r.w_new, snr_db=r.snr_db, rho_loo=r.rho_loo,
                             excess_db=r.excess_db, hidden=k[0] in hid))
        np.savez_compressed(f"{bd}/bank_{tag}.npz",
                            tpl=np.array(tpl, dtype=object), sp=np.array(sp, dtype=object), **extra)
        pd.DataFrame(keys).to_csv(f"{bd}/bank_{tag}_keys.csv", index=False)
        print(f"[step2] bank_{tag}: {len(keys)} (family, station) templates, "
              f"{sum(k['hidden'] for k in keys)} hidden", flush=True)

    # ---- per-station summary + the rho table
    su, rh = [], []
    for sta, g in D.groupby("sta"):
        gg = g[g.status.isin(["scored", "stacked"])]
        su.append(dict(sta=sta, site=ST.site[sta], attempted=len(g),
                       insufficient=int((g.verdict == "insufficient").sum()),
                       no_data=int((g.verdict == "no_data").sum()),
                       untestable=int((g.verdict == "untestable").sum()),
                       accepted=int((g.verdict == "accept").sum()),
                       rejected=int((g.verdict == "reject").sum()),
                       accept_literal_b=int(g.accept.fillna(False).sum()),
                       med_snr_db=float(gg.snr_db.median()) if len(gg) else np.nan,
                       med_excess_db=float(gg.excess_db.median()) if len(gg) else np.nan,
                       med_rho=float(gg.rho_loo.median()) if len(gg) else np.nan,
                       med_hs_cc=float(gg.hs_cc.median()) if "hs_cc" in gg else np.nan))
        for ref in REF_CHOICES:
            c = f"ref_{ref}_rho"
            if c in gg and gg[c].notna().sum() >= 5:
                m = gg[gg[c].notna() & gg.rho_loo.notna()]
                rh.append(dict(sta=sta, ref=ref, n=len(m),
                               rho_new=float(m.rho_loo.median()),
                               rho_ref=float(m[c].median()),
                               ratio=float((m.rho_loo / m[c]).median()),
                               snr_new=float(m.snr_db.median()),
                               snr_ref=float(m[f"ref_{ref}_snr_db"].median())))
    SU = pd.DataFrame(su).sort_values("accepted", ascending=False)
    RH = pd.DataFrame(rh)
    SU.to_csv(f"{S.P3D}/step2_summary{a.tag}.csv", index=False)
    RH.to_csv(f"{S.P3D}/step2_rho{a.tag}.csv", index=False)
    print(SU.to_string(index=False))
    print(RH.to_string(index=False))
    n_symb = int((D[(D.sta == "SYMB") & (D.verdict == "accept")]).shape[0])
    print(f"[step2] ACCEPT gate: (family, SYMB) accepted = {n_symb} "
          f"(proceed >= 30; FAIL if < 10 while the reference legs pass)")


if __name__ == "__main__":
    main()
