#!/usr/bin/env python
"""Track D, D1(d) — IS THE WINTER QUIET-DAY EXCESS REAL LFE ACTIVITY?  (coordinator request 2026-09-17)

Track A step 3b asked, on the 2012 referee, whether quiet-day detections are false alarms or real
inter-ETS repeats, by stacking at a station EXCLUDED from the network sum.  This applies the same idea,
the same configuration and the same null to production months, split by season, so that the winter /
summer quiet-day ratio (~3x, section 1c) can be attributed:

  * real inter-ETS activity  -> the winter quiet-day detections register at the withheld stations more
    strongly than the summer ones;
  * noise                    -> both seasons sit at their nulls and the winter excess is false alarms.

TWO statistics, because the production quiet-day rate is ~4x lower than 2012's and no single family has
enough quiet-day detections for a coherent per-family stack:

  (1) POOLED INCOHERENT (primary).  For every quiet-day detection, cut the withheld station's 3-C trace
      around that family's own predicted S (Crosson tt + the frozen template's S-pick refinement) and
      form r = E(S-2s .. S+2s) / E(pre-arrival 4 s of the same 40-s window).  The mean of log10 r over
      thousands of detections is an INCOHERENT stack: if a fraction p of the detections carry a real S
      arrival at the predicted time, the mean rises above the matched random-time null by ~p x (a/sigma)^2
      in linear units.  It pools across families (each window is normalised by its own pre-arrival
      noise), which is what makes it usable at 2-7 detections per family.
  (2) COHERENT PER FAMILY (secondary, reported where N allows).  Track A step 3b's own statistic: the
      linear stack aligned on the predicted S, its S-window SNR against a matched random-time stack, and
      its cc against the family's 2012 template.

Withheld network: CN.VGZ + CN.SNB + PB.B926 excluded from the network sum (`d1d_run_withheld.sh`), the
frozen 12-s 2012 templates, the frozen HF veto, quiet family-days taken from Track C's own `quiet_flag`,
and every null drawn from the SAME quiet days so that data coverage is matched.

THREE arms, because a raised S-window energy at a withheld station is not by itself proof of an LFE:
  real      the detection times with that family's OWN predicted S offset;
  time-null matched RANDOM times on the same quiet days with the same offset -- tests "is there anything
            there at all";
  mo-null   the SAME detection times with ANOTHER family's offset at that station (offsets permuted
            across families) -- tests whether the energy sits at the FAMILY-SPECIFIC predicted arrival.
            A broadband regional transient that made the network fire would raise the real arm AND the
            moveout-shuffled arm; a real LFE raises only the real arm."""
from __future__ import annotations
import argparse, glob, json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e5r_templates as T  # noqa: E402
import d3_vel1d as v1d  # noqa: E402

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
OUT = f"{EXT}/results"
TC = f"{EXT}/results/trackC"
TD = f"{EXT}/results/trackD"
CACHES = [f"{EXT}/cache20hz", "/home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz"]
FS = T.FS
WITHHELD = ["VGZ", "SNB", "B926"]
HF_CUT = 3.6384
SIG_HALF = 2.0          # +-2 s about the predicted S
NOI_S = 4.0             # the first 4 s of the 40-s window = 20..16 s BEFORE the predicted S


def energies(W, idx):
    """per window: (S-window energy, pre-arrival energy), summed over the three components."""
    ns = int(round(SIG_HALF * FS)); nn = int(round(NOI_S * FS))
    es = (W[:, :, idx - ns:idx + ns] ** 2).sum(axis=(1, 2)) / (2 * ns)
    en = (W[:, :, 0:nn] ** 2).sum(axis=(1, 2)) / nn
    return es, en


def snr_db(st, idx):
    ns = int(round(SIG_HALF * FS)); nn = int(round(NOI_S * FS))
    sig = st[:, idx - ns:idx + ns]; noi = st[:, 0:nn]
    return float(10 * np.log10(max((sig ** 2).mean(), 1e-30) / max((noi ** 2).mean(), 1e-30)))


def cc_tpl(st, idx, tpl, lag_s=1.0):
    nL = tpl.shape[1]; i0 = idx - int(round(1.0 * FS)); lag = int(round(lag_s * FS))
    best = -9.0
    for k in range(-lag, lag + 1):
        a = i0 + k
        if a < 0 or a + nL > st.shape[1]:
            continue
        seg = st[:, a:a + nL]
        nn = np.linalg.norm(seg, axis=1)
        if np.any(nn <= 0):
            continue
        best = max(best, float((seg / nn[:, None] * tpl).sum() / tpl.shape[0]))
    return best


def load_month_days(netsta, ym):
    p = pd.Period(ym, freq="M"); y = p.year
    D = {}
    for d in range(p.start_time.dayofyear, p.end_time.dayofyear + 1):
        f = next((q for q in (f"{c}/{netsta}/{y}/{d:03d}.npz" for c in CACHES) if os.path.exists(q)), None)
        if f is None:
            continue
        z = np.load(f); D[d] = (float(z["t0"]), z["x"])
    return D


def season_months(a):
    return {"DJF": [m for m in a.djf.split(",") if m], "JJA": [m for m in a.jja.split(",") if m]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--djf", default="2015-02,2014-02,2014-12")
    ap.add_argument("--jja", default="2014-08,2015-08,2014-07")
    ap.add_argument("--tiers", default="T_core:8.75,T_ext:8.25")
    ap.add_argument("--draws", type=int, default=40)
    ap.add_argument("--mo-draws", type=int, default=20)
    ap.add_argument("--mo-min-sep-s", type=float, default=5.0,
                    help="the shuffled family's S offset must differ by at least this, so a real arrival "
                         "cannot leak into the +-2 s signal window of the moveout-shuffled arm")
    ap.add_argument("--min-n-fam", type=int, default=15)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--no-hf", action="store_true")
    ap.add_argument("--out", default=f"{TD}/season_stack_test.csv")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    S = pd.read_csv(f"{IN}/e5r_stations.csv").set_index("sta")
    F = pd.read_csv(f"{IN}/e5r_families_2012.csv").set_index("bkey")
    qc = pd.read_csv(f"{IN}/e5r_template_qc_2012.csv"); qc = qc[qc.usable & qc.template]
    z = np.load(f"{IN}/e5r_templates_2012.npz", allow_pickle=True); TPL = z["tpl"].item()
    pick_sta = {}
    for f, g in qc[qc.sta.isin(WITHHELD)].groupby("bkey"):
        r = g.sort_values("snr_db", ascending=False).iloc[0]
        pick_sta[int(f)] = (r.sta, float(r.s_pick_s))
    n = int(round(2 * T.HALF_W * FS)); idx = int(round(T.HALF_W * FS))
    TIERS = [(k, float(v)) for k, v in (x.split(":") for x in a.tiers.split(","))]
    SEA = season_months(a)

    pooled, fam_rows, STK = [], [], {}
    for sea, months in SEA.items():
        for tier, Tt in TIERS:
            # ---- gather, per (station, month), the detection times and their family S offsets
            agg = {}          # sta -> dict(real=[log10 r], null=[[log10 r] per draw])
            nq_tot = ndet_tot = 0
            fam_pack = {}     # (sta, fam) -> list of window arrays
            for ym in months:
                fp = f"{TC}/family_days_{ym}.parquet"
                dp = f"{OUT}/e5r_detections_d1d_wh_{ym}.parquet"
                if not (os.path.exists(fp) and os.path.exists(dp)):
                    continue
                fd = pd.read_parquet(fp, columns=["fam", "date", "observable", "quiet_flag"])
                fd = fd[fd.observable & fd.quiet_flag][["fam", "date"]]
                d = pd.read_parquet(dp); d = d[d.L == 12.0].copy()
                hp = f"{OUT}/a1_hf_ratio_d1d_wh_{ym}.parquet"
                if not a.no_hf and os.path.exists(hp):
                    H = pd.read_parquet(hp); H = H[H.L == 12.0][["doy", "fam", "t", "hfr_PGC"]]
                    d = d.merge(H, on=["doy", "fam", "t"], how="left")
                    v = d.hfr_PGC.to_numpy(); d = d[~(np.isfinite(v) & (v > HF_CUT))]
                d["date"] = [str(pd.Timestamp(t_, unit="s").date()) for t_ in d.t]
                d = d.merge(fd, on=["fam", "date"], how="inner")
                d = d[d.net_mad >= Tt]
                nq_tot += len(fd); ndet_tot += len(d)
                for sta in WITHHELD:
                    fams = [f_ for f_, (s_, _) in pick_sta.items() if s_ == sta]
                    dd = d[d.fam.astype(int).isin(fams)]
                    if not len(dd):
                        continue
                    D = load_month_days(S.loc[sta, "netsta"], ym)
                    if not D:
                        continue
                    qd = fd[fd.fam.astype(int).isin(fams)]
                    doys = sorted({pd.Timestamp(x).dayofyear for x in qd.date} & set(D))
                    A_ = agg.setdefault(sta, dict(real=[], null=[[] for _ in range(a.draws)],
                                                  mo=[[] for _ in range(a.mo_draws)]))
                    offs_all = {}
                    for f2 in fams:
                        fr2 = F.loc[f2]
                        dx2 = T.epi_km(fr2.lat, fr2.lon, float(S.loc[sta, "lat"]), float(S.loc[sta, "lon"]))
                        offs_all[f2] = float(v1d.tt(dx2, fr2.dep, "S")) + pick_sta[f2][1] - T.HALF_W
                    for f_, g in dd.groupby(dd.fam.astype(int)):
                        fr = F.loc[f_]; pick = pick_sta[f_][1]
                        dx = T.epi_km(fr.lat, fr.lon, float(S.loc[sta, "lat"]), float(S.loc[sta, "lon"]))
                        off = float(v1d.tt(dx, fr.dep, "S")) + pick - T.HALF_W
                        tt = np.sort(g.t.to_numpy())
                        W, _ = T.cut_windows(D, tt + off, n)
                        if W is not None and len(W):
                            es, en = energies(W, idx)
                            k = np.isfinite(es) & np.isfinite(en) & (en > 0) & (es > 0)
                            A_["real"].append(np.log10(es[k] / en[k]))
                            fam_pack.setdefault((sta, f_), []).append(W)
                        # moveout-shuffled arm: same times, another family's offset at this station
                        cand = [o for f2, o in offs_all.items()
                                if f2 != f_ and abs(o - off) > a.mo_min_sep_s]
                        for kmo in range(a.mo_draws):
                            if not cand:
                                break
                            Wm, _ = T.cut_windows(D, tt + cand[rng.integers(len(cand))], n)
                            if Wm is not None and len(Wm):
                                es, en = energies(Wm, idx)
                                k = np.isfinite(es) & np.isfinite(en) & (en > 0) & (es > 0)
                                A_["mo"][kmo].append(np.log10(es[k] / en[k]))
                        if not doys:
                            continue
                        for kdr in range(a.draws):
                            dsel = rng.choice(doys, len(tt))
                            t0s = np.array([D[int(x)][0] for x in dsel]) + rng.uniform(0, 86000, len(tt))
                            Wn, _ = T.cut_windows(D, t0s + off, n)
                            if Wn is not None and len(Wn):
                                es, en = energies(Wn, idx)
                                k = np.isfinite(es) & np.isfinite(en) & (en > 0) & (es > 0)
                                A_["null"][kdr].append(np.log10(es[k] / en[k]))
                    del D
            # ---- (1) pooled incoherent, per station and over all three
            allr, alln, allmo = [], [[] for _ in range(a.draws)], [[] for _ in range(a.mo_draws)]
            for sta, A_ in agg.items():
                r = np.concatenate(A_["real"]) if A_["real"] else np.array([])
                nl = [np.concatenate(x) if x else np.array([]) for x in A_["null"]]
                mo = [np.concatenate(x) if x else np.array([]) for x in A_["mo"]]
                if len(r) < 20:
                    continue
                nm = np.array([v.mean() for v in nl if len(v) >= 20])
                mm = np.array([v.mean() for v in mo if len(v) >= 20])
                if not len(nm):
                    continue
                exc = 10 * (r.mean() - nm.mean())         # dB, log10 -> 10*
                se = 10 * r.std() / np.sqrt(len(r))
                pooled.append(dict(season=sea, tier=tier, T=Tt, sta=sta, level="station",
                                   n_quiet_family_days=nq_tot, n_det=int(len(r)),
                                   mean_log10r=float(r.mean()), null_mean_log10r=float(nm.mean()),
                                   null_sd_log10r=float(nm.std()), excess_db=float(exc),
                                   excess_db_se=float(se), z=float(exc / max(10 * nm.std(), 1e-9)),
                                   p_null=float((nm >= r.mean()).sum() + 1) / (len(nm) + 1),
                                   mo_mean_log10r=float(mm.mean()) if len(mm) else np.nan,
                                   excess_vs_mo_db=float(10 * (r.mean() - mm.mean())) if len(mm) else np.nan,
                                   p_mo=(float((mm >= r.mean()).sum() + 1) / (len(mm) + 1)) if len(mm) else np.nan))
                allr.append(r)
                for kdr, v in enumerate(nl):
                    if len(v):
                        alln[kdr].append(v)
                for kmo, v in enumerate(mo):
                    if len(v):
                        allmo[kmo].append(v)
            if allr:
                r = np.concatenate(allr)
                nm = np.array([np.concatenate(x).mean() for x in alln if x])
                mm = np.array([np.concatenate(x).mean() for x in allmo if x])
                exc = 10 * (r.mean() - nm.mean())
                pooled.append(dict(season=sea, tier=tier, T=Tt, sta="ALL", level="pooled",
                                   n_quiet_family_days=nq_tot, n_det=int(len(r)),
                                   mean_log10r=float(r.mean()), null_mean_log10r=float(nm.mean()),
                                   null_sd_log10r=float(nm.std()), excess_db=float(exc),
                                   excess_db_se=float(10 * r.std() / np.sqrt(len(r))),
                                   z=float(exc / max(10 * nm.std(), 1e-9)),
                                   p_null=float((nm >= r.mean()).sum() + 1) / (len(nm) + 1),
                                   mo_mean_log10r=float(mm.mean()) if len(mm) else np.nan,
                                   excess_vs_mo_db=float(10 * (r.mean() - mm.mean())) if len(mm) else np.nan,
                                   p_mo=(float((mm >= r.mean()).sum() + 1) / (len(mm) + 1)) if len(mm) else np.nan))
                print(f"[{sea} {tier}] pooled incoherent: {len(r)} detections at the withheld stations, "
                      f"excess {exc:+.3f} dB (null sd {10*nm.std():.3f})", flush=True)
            # ---- (2) coherent per family, where N allows
            for (sta, f_), parts in fam_pack.items():
                Wc = np.concatenate(parts, 0)
                if len(Wc) < a.min_n_fam:
                    continue
                st = Wc.mean(0)
                tp = TPL.get((float(f_), sta), {}).get(8.0)
                fam_rows.append(dict(season=sea, tier=tier, T=Tt, sta=sta, fam=f_, n_stacked=len(Wc),
                                     snr_db=snr_db(st, idx),
                                     cc_template=cc_tpl(st, idx, tp) if tp is not None else np.nan))
                STK[f"{sea}|{tier}|{sta}|{f_}"] = st.astype(np.float32)
    P = pd.DataFrame(pooled)
    if not len(P):
        raise SystemExit("[d1d] nothing to test — run scripts/d1d_run_withheld.sh first")
    P.to_csv(a.out, index=False)
    if fam_rows:
        pd.DataFrame(fam_rows).to_csv(a.out.replace(".csv", "_byfamily.csv"), index=False)
        np.savez_compressed(a.out.replace(".csv", "_stacks.npz"), **STK)
    SU = P[P.level == "pooled"].copy()
    # a real-repeat fraction, scaled off Track A's 2012 ETS positive control on the same stations
    ctl = pd.read_csv(f"{OUT}/a3b_quietday_stack_test.csv")
    ctl = ctl[(ctl.set == "ets_control") & ctl.excess_db.notna()]
    SU["ets_control_excess_db_2012"] = float(ctl.excess_db.median()) if len(ctl) else np.nan
    SU.to_csv(a.out.replace(".csv", "_summary.csv"), index=False)
    pd.set_option("display.width", 250)
    print("\n[pooled incoherent withheld-station test, by season]")
    print(SU[["season", "tier", "n_quiet_family_days", "n_det", "mean_log10r", "null_mean_log10r",
              "excess_db", "excess_db_se", "z", "p_null", "mo_mean_log10r", "excess_vs_mo_db",
              "p_mo"]].round(4).to_string(index=False))
    print("\n[per station]")
    print(P[P.level == "station"][["season", "tier", "sta", "n_det", "excess_db", "excess_db_se", "z",
                                   "p_null", "excess_vs_mo_db", "p_mo"]].round(4).to_string(index=False))
    json.dump(dict(seasons=SEA, tiers=TIERS, draws=a.draws, sig_half_s=SIG_HALF, noise_s=NOI_S,
                   statistic="mean log10 of (S-window energy / pre-arrival energy) at a station withheld "
                             "from the network sum; null = matched random times on the same quiet days"),
              open(a.out.replace(".csv", "_meta.json"), "w"), indent=1)
    print(f"[write] {a.out} (+ _summary.csv, _byfamily.csv, _stacks.npz, _meta.json)")


if __name__ == "__main__":
    main()
