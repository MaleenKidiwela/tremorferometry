#!/usr/bin/env python
"""Step 5b -- INJECTION CONTROL: plant northern high-purity LFEs into Washington station-days.

Merlin F4/Q5(b).  Gate 1 has no Washington control table, so a non-PASS is UNTESTABLE rather than
FAIL.  The injection supplies one -- and, with it, the location-error distribution of the whole
chain at a KNOWN cell.  It is EXPLICITLY NOT a sensitivity ceiling: `notes/2026-09-11_Notes.md`
E3 step 3 records that injected stacks recover 4-11x off real at matched s/n, so "the injected
family was recovered at amplitude A" says nothing about what amplitude a real family needs.

DESIGN (frozen in results/freeze_wa_v1.json before anything is planted)
  * 4 synthetic families at 4 RANDOM active cells of the sub-region (seeded, not chosen by
    coverage, not chosen near the Chestler coordinates);
  * each cell planted at 3 AMPLITUDE LEVELS -- 0.5x, 1.0x and 2.0x the donor event's OWN
    signal-to-noise at its own northern station -- so the 12 planted repeaters span the north's
    observed registration;
  * ~40 DISTINCT northern hp events per synthetic family (`lfe_extension/results/trackE/
    realness_v1_1.parquet`, hp == True), drawn from distinct northern families, with the 30 hidden
    families of `lfe_picker/results/holdouts.json` EXCLUDED as donors;
  * DONOR CLASS MATCHED TO RECIPIENT CLASS: a Washington borehole receives a northern borehole
    snippet, a Washington surface station receives a northern surface snippet, and within the class
    the donor station is the one whose own epicentral distance is closest to the recipient's
    distance from the planted cell;
  * inserted into ADMITTED **and** WITHHELD station-days at the planted cell's moveout
    (Crosson 1-D + the frozen station correction), additively, on the 20 Hz 1-8 Hz cache;
  * the injected copy is a SEPARATE cache root (`results/step5/inject/cache20hz`): every bed day is
    a SYMLINK to the real cache and only the injected station-days are real files, so the real
    cache is never modified and the same chain can be pointed at either root with `X_CACHE`.

    python scripts/x_inject.py build   --out results/step5/inject
    python scripts/x_inject.py score   --root results/step5/inject --lg results/step5/inject/lg
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
import x_common as X                                                        # noqa: E402

EXT = "/home/jovyan/tremorferometry/lfe_extension"
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
DONOR_CACHES = [f"{EXT}/cache20hz", f"{REP}/cache20hz", f"{X.PICKER}/cache20hz"]
N_FAMILIES = 4
N_EVENTS = 40
AMPS = (0.5, 1.0, 2.0)
N_DAYS = 6
SNIP_PRE_S, SNIP_LEN_S = 2.0, 12.0          # 12 s from the donor's S - 2 s (the template geometry)
SEARCH_S = 3.0                              # envelope-peak search about the predicted donor S
MIN_SEP_S = 120.0                           # planted events are >= 2 min apart across ALL families
SEED = 20260920
BOREHOLE_NETS = ("PB",)


def donor_class(sta, net):
    return "borehole" if str(net) in BOREHOLE_NETS and str(sta).startswith("B") else "surface"


def north_stations():
    """Donor stations with 20 Hz day files on the northern caches, with lat/lon and class."""
    import glob
    import g_common as G
    S = G.stations()
    have = {}
    for c in DONOR_CACHES:
        for p in glob.glob(f"{c}/*/*/*.npz"):
            parts = p.split(os.sep)
            sta = parts[-3].split(".")[1]
            have.setdefault(sta, {})[(int(parts[-2]), int(parts[-1][:3]))] = p
    rows = []
    for s in sorted(have):
        if s not in S.index:
            continue
        rows.append(dict(sta=s, netsta=str(S.loc[s, "netsta"]), lat=float(S.loc[s, "lat"]),
                         lon=float(S.loc[s, "lon"]),
                         klass=donor_class(s, str(S.loc[s, "netsta"]).split(".")[0]),
                         n_days=len(have[s])))
    return pd.DataFrame(rows), have


def load_north_day(path):
    z = np.load(path)
    x = z["x"]
    if x.shape[0] != 3:
        return None
    return float(z["t0"]), x


def snippet(path, t_pred, fs=20.0):
    """(3, n) donor snippet, 12 s from the envelope peak - 2 s within +-3 s of the predicted S,
    plus its own signal-to-noise (rms in the snippet / rms of the 12 s before it)."""
    import s_common as S
    v = load_north_day(path)
    if v is None:
        return None, np.nan
    t0, Xd = v
    n = int(round(SNIP_LEN_S * fs))
    w = int(round(SEARCH_S * fs))
    i_pred = int(round((t_pred - t0) * fs))
    a = i_pred - w - int(round(SNIP_PRE_S * fs))
    b = i_pred + w + n
    if a - n < 0 or b + n > Xd.shape[1]:
        return None, np.nan
    env = S.envelope(np.asarray(Xd[:, i_pred - w:i_pred + w], float))
    i_on = i_pred - w + int(np.argmax(env))
    i0 = i_on - int(round(SNIP_PRE_S * fs))
    if i0 - n < 0 or i0 + n > Xd.shape[1]:
        return None, np.nan
    seg = np.asarray(Xd[:, i0:i0 + n], float)
    noi = np.asarray(Xd[:, i0 - n:i0], float)
    r_s = float(np.sqrt((seg ** 2).mean()))
    r_n = float(np.sqrt((noi ** 2).mean()))
    if not np.isfinite(r_s) or r_s <= 0 or not np.isfinite(r_n) or r_n <= 0:
        return None, np.nan
    return seg / r_s, r_s / r_n


def build(out_root, n_families=N_FAMILIES, n_events=N_EVENTS,
          amps=AMPS, n_days=N_DAYS, only_stations=None):
    import glob
    import g_common as G
    import s_common as S
    W = X.activate()
    C = X.cfg()
    rng = np.random.default_rng(SEED)
    os.makedirs(f"{out_root}/cache20hz", exist_ok=True)
    cells = pd.read_csv(f"{X.OUT}/active_cells.csv")
    pick = rng.choice(len(cells), size=n_families, replace=False)
    fam_cells = cells.iloc[np.sort(pick)].reset_index(drop=True)
    X.log(f"[inject] planted cells:\n{fam_cells.to_string(index=False)}")

    # ---- injection days: the N_DAYS window days with the most candidate hours
    H = W.candidate_hours()
    cnt = H.groupby("date").size().sort_values(ascending=False)
    days = sorted(cnt.head(n_days).index)
    hours = H[H.date.isin(days)].reset_index(drop=True)
    X.log(f"[inject] {len(days)} injection days {days}, {len(hours)} candidate hours")

    # ---- donors: northern hp events, hidden families excluded
    hold = json.load(open(f"{X.PICKER}/results/holdouts.json"))
    hidden = {int(h["bkey"]) for h in hold["hidden_families"]}
    R = pd.read_parquet(f"{EXT}/results/trackE/realness_v1_1.parquet",
                        columns=["family_id", "t", "date", "hp", "p_real", "net_mad"])
    R = R[R.hp & ~R.family_id.isin(hidden)].reset_index(drop=True)
    F = pd.read_csv(f"{EXT}/results/trackD/family_table_v1_1.csv")
    F = F.set_index("family_id")[["lat_bostock", "lon_bostock", "depth_bostock"]]
    NS, have = north_stations()
    X.log(f"[inject] {len(R)} northern hp donor events over {R.family_id.nunique()} families "
          f"({len(hidden)} hidden excluded); {len(NS)} donor stations "
          f"({NS.klass.value_counts().to_dict()})")

    # ---- planted event times: >= MIN_SEP_S apart across every family
    n_plant = n_families * len(amps)
    times = []
    while len(times) < n_plant * n_events:
        r = hours.iloc[int(rng.integers(0, len(hours)))]
        t = pd.Timestamp(r.date).timestamp() + int(r.hour) * 3600 + float(rng.uniform(60, 3400))
        if all(abs(t - u) >= MIN_SEP_S for u in times[-400:]):
            times.append(t)
    times = np.sort(np.array(times))
    rng.shuffle(times)

    ST = X.wa_station_table()
    corr = X.station_corrections()
    mem = C["site_members"]
    recips = sorted(set([s for q in C["site_order"] for s in mem[q]]
                        + list(C["withheld_stations"])))
    recips = [s for s in recips if s in ST.index]
    if only_stations:
        recips = [s for s in recips if s in set(only_stations)]

    # ---- plan
    plan, k = [], 0
    for fi in range(n_families):
        cell = fam_cells.iloc[fi]
        don = (R.sample(n=min(n_events * len(amps) * 4, len(R)), random_state=SEED + fi)
                 .drop_duplicates("family_id", keep="first"))
        if len(don) < n_events:
            don = R.sample(n=n_events * 3, random_state=SEED + 100 + fi).head(n_events)
        don = don.head(n_events * len(amps)).reset_index(drop=True)
        for ai, amp in enumerate(amps):
            sid = f"SYN{fi}A{ai}"
            dd = don.iloc[ai * n_events:(ai + 1) * n_events]
            if len(dd) < n_events:
                dd = don.sample(n=n_events, replace=len(don) < n_events,
                                random_state=SEED + 7 * fi + ai)
            for j, d in enumerate(dd.itertuples()):
                plan.append(dict(sid=sid, fam=fi, amp=amp, lat=float(cell.lat),
                                 lon=float(cell.lon), dep=float(cell.dep),
                                 t_plant=float(times[k]), donor_family=int(d.family_id),
                                 donor_t=float(d.t), donor_date=str(d.date)))
                k += 1
    P = pd.DataFrame(plan)
    X.log(f"[inject] {len(P)} planted events = {n_plant} synthetic families x {n_events}")

    # ---- plant, station-day by station-day
    idx = W.cache_idx()
    os.makedirs(f"{out_root}/cache20hz", exist_ok=True)
    src_days = {}
    for (sta, yr, doy), p in idx.items():
        src_days.setdefault((sta, yr, doy), p)
    inj_doys = {int(pd.Timestamp(d).dayofyear) for d in days}
    wrote, rows = 0, []
    for sta in recips:
        ns = f"{ST.net[sta]}.{sta}"
        os.makedirs(f"{out_root}/cache20hz/{ns}/{X.YEAR}", exist_ok=True)
        kl = donor_class(sta, str(ST.net[sta]))
        pool = NS[NS.klass == kl]
        if not len(pool):
            pool = NS
        for (s2, yr, doy), p in sorted(src_days.items()):
            if s2 != sta:
                continue
            dst = f"{out_root}/cache20hz/{ns}/{yr}/{doy:03d}.npz"
            if doy not in inj_doys:
                if not os.path.lexists(dst):
                    os.symlink(os.path.abspath(p), dst)
                continue
            if os.path.exists(dst) and not os.path.islink(dst):
                continue
            z = np.load(p)
            t0, Xw, order, fs = float(z["t0"]), np.array(z["x"], dtype=np.float32), z["order"], z["fs"]
            day_lo = t0
            day_hi = t0 + Xw.shape[1] / 20.0
            sub = P[(P.t_plant >= day_lo) & (P.t_plant < day_hi - 60)]
            dsta_lat, dsta_lon = float(ST.lat[sta]), float(ST.lon[sta])
            nplant = 0
            for r in sub.itertuples():
                dist = float(S.epi_km(r.lat, r.lon, dsta_lat, dsta_lon))
                ttS, _ = G.ttSP(dist, r.dep)
                ttS += corr.get(sta, (0.0, 0.0))[0]
                if r.donor_family not in F.index:
                    continue
                flat, flon, fdep = (float(F.lat_bostock[r.donor_family]),
                                    float(F.lon_bostock[r.donor_family]),
                                    float(F.depth_bostock[r.donor_family]))
                dd = np.array([float(S.epi_km(flat, flon, a, b))
                               for a, b in zip(pool.lat, pool.lon)])
                cand = pool.iloc[np.argsort(np.abs(dd - dist))]
                dtime = pd.Timestamp(r.donor_date)
                key = (int(dtime.year), int(dtime.dayofyear))
                seg, snr = None, np.nan
                for c2, dist_d in zip(cand.itertuples(), np.sort(np.abs(dd - dist))):
                    pth = have.get(c2.sta, {}).get(key)
                    if pth is None:
                        continue
                    ttd, _ = G.ttSP(float(S.epi_km(flat, flon, c2.lat, c2.lon)), fdep)
                    seg, snr = snippet(pth, r.donor_t + ttd)
                    if seg is not None:
                        dsta = c2.sta
                        break
                if seg is None:
                    continue
                i0 = int(round((r.t_plant + ttS - SNIP_PRE_S - t0) * 20.0))
                n = seg.shape[1]
                if i0 < n or i0 + 2 * n > Xw.shape[1]:
                    continue
                noi = Xw[:, i0 - n:i0].astype(np.float64)
                rn = float(np.sqrt((noi ** 2).mean()))
                if not np.isfinite(rn) or rn <= 0 or not np.isfinite(snr):
                    continue
                scale = r.amp * snr * rn                # target rms = amp x donor SNR x local noise
                Xw[:, i0:i0 + n] += (scale * seg).astype(np.float32)
                nplant += 1
                rows.append(dict(sid=r.sid, sta=sta, donor_sta=dsta, donor_family=r.donor_family,
                                 t_plant=r.t_plant, t_s=r.t_plant + ttS, dist_km=dist,
                                 amp=r.amp, donor_snr=float(snr), local_noise_rms=rn,
                                 planted_rms=float(scale)))
            if os.path.islink(dst):
                os.unlink(dst)
            np.savez(dst, x=Xw, t0=t0, order=order, fs=fs)
            wrote += 1
            if wrote % 50 == 0:
                X.log(f"[inject] {wrote} injected station-days, {len(rows)} planted snippets")
    D = pd.DataFrame(rows)
    os.makedirs(out_root, exist_ok=True)
    P.to_csv(f"{out_root}/plan.csv", index=False)
    D.to_csv(f"{out_root}/planted.csv", index=False)
    hours.to_csv(f"{out_root}/hours.csv", index=False)
    truth = (P.groupby("sid").agg(fam=("fam", "first"), amp=("amp", "first"),
                                  lat=("lat", "first"), lon=("lon", "first"),
                                  dep=("dep", "first"), n_planned=("t_plant", "size"))
             .reset_index())
    if len(D):
        got = D.groupby("sid").agg(n_stations=("sta", "nunique"),
                                   n_events=("t_plant", "nunique")).reset_index()
        truth = truth.merge(got, on="sid", how="left")
    truth.to_csv(f"{out_root}/truth.csv", index=False)
    json.dump(dict(n_synthetic_families=int(len(truth)), n_events_each=n_events,
                   amplitudes=list(amps), n_days=len(days), days=days,
                   n_injected_station_days=int(wrote), n_planted_snippets=int(len(D)),
                   seed=SEED, cells=fam_cells.to_dict("records"),
                   note=("mechanics + location-error distribution ONLY; never a sensitivity "
                         "ceiling (E3: injected stacks recover 4-11x off real)")),
              open(f"{out_root}/meta.json", "w"), indent=1)
    X.log(f"[inject] wrote {wrote} injected station-days, {len(D)} planted snippets, "
          f"{len(truth)} synthetic families -> {out_root}")


def score(root, lg):
    """Recovery of the planted repeaters: match located/gated proposals to the planted cells."""
    import s_common as S
    T = pd.read_csv(f"{root}/truth.csv")
    D = pd.read_csv(f"{root}/planted.csv") if os.path.exists(f"{root}/planted.csv") \
        else pd.DataFrame()
    P = pd.read_csv(f"{root}/plan.csv")
    L = pd.read_csv(f"{lg}_locations.csv") if os.path.exists(f"{lg}_locations.csv") \
        else pd.DataFrame()
    Gv = pd.read_csv(f"{lg}_gates.csv") if os.path.exists(f"{lg}_gates.csv") else pd.DataFrame()
    C = (pd.read_parquet(f"{lg}_clusters.parquet") if os.path.exists(f"{lg}_clusters.parquet")
         else pd.DataFrame(columns=["cid", "t", "lat", "lon"]))
    rows = []
    for t in T.itertuples():
        pl = P[P.sid == t.sid].t_plant.to_numpy(float)
        best = dict(sid=t.sid, amp=t.amp, lat=t.lat, lon=t.lon, n_planned=int(t.n_planned),
                    recovered=False, cid=None, n_matched=0, frac_matched=0.0,
                    loc_err_km=np.nan, gate_verdict=None, g1_sites_pass=np.nan)
        for cid, g in (C.groupby("cid") if len(C) else []):
            tt = g.t.to_numpy(float)
            m = int(sum(1 for u in pl if np.min(np.abs(tt - u)) <= 2.0)) if len(tt) else 0
            if m > best["n_matched"]:
                gv = Gv[Gv.pid.astype(str) == str(cid)] if len(Gv) else pd.DataFrame()
                best.update(cid=str(cid), n_matched=m, frac_matched=m / max(len(pl), 1),
                            loc_err_km=float(S.epi_km(float(g.lat.iloc[0]), float(g.lon.iloc[0]),
                                                      t.lat, t.lon)),
                            gate_verdict=(str(gv.verdict.iloc[0]) if len(gv) else None),
                            g1_sites_pass=(float(gv.g1_sites_pass.iloc[0])
                                           if len(gv) and "g1_sites_pass" in gv else np.nan),
                            recovered=bool(m >= 0.3 * len(pl)))
        rows.append(best)
    Rr = pd.DataFrame(rows)
    Rr.to_csv(f"{root}/recovery.csv", index=False)
    # Gate-1 control table for the real run: one row per withheld station per recovered family
    ctl = []
    for r in Rr[Rr.recovered].itertuples():
        full = f"{lg}_gate_full_{r.cid}.json"
        if not os.path.exists(full):
            continue
        j = json.load(open(full))
        for p in j.get("gate1", {}).get("per_station", []):
            ctl.append(dict(sta=p.get("sta"), N=int(p.get("n_used", 0)),
                            dist_km=float(p.get("dist_km", np.nan)),
                            excess_db=float(p.get("excess_db", np.nan)) if p.get("excess_db")
                            is not None else np.nan,
                            verdict=p.get("verdict"), sid=r.sid, amp=r.amp))
    pd.DataFrame(ctl).to_csv(f"{root}/gate1_controls_wa.csv", index=False)
    X.log(f"[inject/score] {int(Rr.recovered.sum())}/{len(Rr)} planted repeaters recovered; "
          f"median location error "
          f"{float(Rr[Rr.recovered].loc_err_km.median()) if Rr.recovered.any() else float('nan'):.2f} km; "
          f"{len(ctl)} Gate-1 control rows")
    print(Rr.to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "score"])
    ap.add_argument("--out", default=f"{X.OUT}/inject")
    ap.add_argument("--root", default=f"{X.OUT}/inject")
    ap.add_argument("--lg", default=f"{X.OUT}/inject/lg")
    ap.add_argument("--n-families", type=int, default=N_FAMILIES)
    ap.add_argument("--n-events", type=int, default=N_EVENTS)
    ap.add_argument("--n-days", type=int, default=N_DAYS)
    ap.add_argument("--amps", default=",".join(str(x) for x in AMPS))
    ap.add_argument("--only-stations", default="")
    a = ap.parse_args()
    if a.cmd == "build":
        build(a.out, n_families=a.n_families, n_events=a.n_events,
              amps=tuple(float(x) for x in a.amps.split(",")), n_days=a.n_days,
              only_stations=[s for s in a.only_stations.split(",") if s])
    else:
        score(a.root, a.lg)


if __name__ == "__main__":
    main()
