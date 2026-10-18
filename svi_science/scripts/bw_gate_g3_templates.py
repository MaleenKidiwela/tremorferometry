#!/usr/bin/env python
"""bw_gate_g3_templates -- G3: era templates for the POLARIS + LZB stations from the 2004-07 ETS (BUILD blocks only).

Build events are re-timed onto the frozen axis first (bw_gate_retime.py: the 2004 catalogue seconds are truncated).
The forward template recipe is called unchanged: lfe_extension/scripts/e5r_templates.do_station (predicted S, +-20 s,
rms-normalised, cc-aligned, S pick, S-SNR vs a same-N random-time stack, 8/12-s templates from S-1 s, 3-C L2-normalised,
placement sp, amplitude reference on the quiet days).  Usable = S-SNR >= null + 3 dB, >= 30 build events.  Weights
w = min(10^(SNR/20), 3 M) / M with M = the median of the 2012 usable pairs (so era weights sit on the frozen scale).

Time axis: Bostock's 2004 event times are ~0.3-0.6 s off the frozen axis (bw_gate_g2_lagdiag.py).  Per family the
build-block event times are shifted by delta_f = median per-event lag (|lag| <= 3 s, cc >= 0.3, >= 10 events) of the
NEVER-changed PGC frozen 12-s template at PGC; families without enough take the median delta over families.

Out: OUT/tables/g3_era_templates.npz (tpl/sp, --extra-tpl format), g3_era_weights.csv, g3_era_ampref.csv,
     g3_template_qc.csv, g3_delta_pgc.csv, g3_carriers.csv (per family: frozen CN4 carriers + era carriers)
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_gate_g2_epoch as E
G = E.G
TAB = E.TAB
QUIET = range(156, 170)        # 2004-06-04 .. 06-17


def frozen_offsets():
    R, C = G.patch_forward()
    fz = C.frozen(); L = float(fz["template_length_s"]); tpath, _ = C.resolve_templates(fz, "")
    tpl, sp, F, wgt, mv, ampref = C.load_frozen_tables(fz, tpath, L)
    return tpl, C.placement(tpl, sp, F, mv, C.station_table(), L)


def delta_pgc(ev, year=2004, caches=None):
    """per-family lag of the frozen PGC template at PGC around t_origin + frozen placement (detector axis)."""
    S = E.setup(); T = S["T"]
    T.CACHES = caches or [G.CACHE_BW]; T.YEAR = year
    D = T.load_days("CN.PGC", year)
    tp, off = frozen_offsets()
    QC = S["QC"]; QC = QC[(QC.sta == "PGC") & QC.usable & QC.template].set_index("bkey")
    rows = []
    for f in sorted(ev.fam.unique()):
        if (float(f), "PGC") not in tp:
            continue
        w = np.asarray(tp[(float(f), "PGC")][12.0], float); n = w.shape[1]
        if (float(f), "PGC") not in off:
            continue
        g = ev[ev.fam == f].sort_values("mag", ascending=False).head(400)
        Sw, _ = T.cut_windows(D, g.t_origin.to_numpy() + off[(float(f), "PGC")] - 3.0, n + 120)   # frozen placement, +-3 s scan
        if Sw is None:
            continue
        from numpy.lib.stride_tricks import sliding_window_view as swv
        lags = []
        for s in Sw:
            V3 = swv(s, n, axis=1)                         # (3, 121, n)
            num = np.einsum("ckn,cn->k", V3, w)
            den = np.sqrt(np.einsum("ckn,ckn->k", V3, V3)) * np.sqrt((w * w).sum())
            c = num / np.maximum(den, 1e-20); k = int(np.argmax(c))
            if c[k] >= 0.3:
                lags.append((k - 60) / 20.0)
        rows.append(dict(fam=int(f), n_ev=len(Sw), n_good=len(lags), delta=float(np.median(lags)) if len(lags) >= 10 else np.nan,
                         delta_mad=float(np.median(np.abs(np.array(lags) - np.median(lags)))) if len(lags) >= 10 else np.nan))
    Dl = pd.DataFrame(rows)
    med = float(Dl.delta.median())
    Dl["delta_used"] = Dl.delta.fillna(med); Dl["fallback"] = Dl.delta.isna()
    return Dl, med


def main():
    S = E.setup(); T = S["T"]
    ev = E.events("2004-07", build_blocks=True)
    Dl, med = delta_pgc(ev)                                   # diagnostic: per-family PGC lag at the catalogue times
    Dl.to_csv(f"{TAB}/g3_delta_pgc.csv", index=False)
    G.log(f"[g3 tpl] PGC delta at catalogue times: median {med:+.3f} s over {Dl.delta.notna().sum()} families")
    # Deviation (recorded in GATES.md): 2004 catalogue seconds are TRUNCATED to integers, so instead of one shift per
    # family every build event is re-timed onto the frozen axis by the CN4 frozen templates (bw_gate_retime.py).
    import bw_gate_retime as RT
    ev = RT.retime(ev, 2004, [G.CACHE_BW], list(G.CN4))
    ev[["fam", "t_cat", "t_origin", "lag", "net_cc"]].to_parquet(f"{TAB}/g3_build_events_retimed.parquet", index=False)
    G.log(f"[g3 tpl] re-timed {len(ev)} build events: lag median {ev.lag.median():+.2f} s, IQR {ev.lag.quantile(.25):+.2f}..{ev.lag.quantile(.75):+.2f}")
    ev = ev.assign(bkey=ev.fam.astype(int))
    F = pd.read_csv(f"{G.REP}/results/e5r_families_2012.csv")
    cnt = ev.groupby("bkey").size(); F = F.assign(n_build=F.bkey.map(cnt).fillna(0).astype(int))
    fam_tbl = F[F.n_build >= 30][["bkey", "sb_fam_id", "lat", "lon", "dep", "n_build"]].to_dict("records")
    EP = G.epochs_table().drop_duplicates("sta").set_index("sta")
    T.YEAR = 2004; T.CACHES = [G.CACHE_BW]; T.QUIET_DOYS = QUIET; T.LENS = (8.0, 12.0)
    jobs = [(s, G.NETSTA[s], float(EP.loc[s, "lat"]), float(EP.loc[s, "lon"]), "era", fam_tbl, ev[["bkey", "t_origin"]], {}, 500 + i)
            for i, s in enumerate(G.ERA)]
    from multiprocessing import Pool
    with Pool(min(12, len(jobs))) as pool:
        out = pool.map(T.do_station, jobs)
    rows, TPL, SP, REF = [], {}, {}, []
    for sta, rw, Tt, P, K, ref in out:
        rows += rw; TPL.update(Tt); SP.update(P)
        for L, v in ref.items():
            REF.append(dict(sta=sta, L=L, e3_ref=v))
    Q = pd.DataFrame(rows); Q.to_csv(f"{TAB}/g3_template_qc.csv", index=False)
    pd.DataFrame(REF).to_csv(f"{TAB}/g3_era_ampref.csv", index=False)
    use = Q[Q.usable & Q.template]
    TPLu = {k: v for k, v in TPL.items() if ((use.bkey == int(k[0])) & (use.sta == k[1])).any()}
    np.savez_compressed(f"{TAB}/g3_era_templates.npz", tpl=np.array(TPLu, dtype=object), sp=np.array({k: SP[k] for k in TPLu}, dtype=object))
    W12 = pd.read_csv(f"{G.REP}/results/e5r_weights_snr_2012.csv")
    M = float(np.median(10 ** (W12.snr_db / 20.0)))
    w = np.minimum(10 ** (use.snr_db / 20.0), 3 * M) / M
    pd.DataFrame(dict(fam=use.bkey.astype(float), sta=use.sta, w=w.to_numpy(), snr_db=use.snr_db, n_events=use.n_events)) \
        .to_csv(f"{TAB}/g3_era_weights.csv", index=False)
    # carriers per family
    QC12 = S["QC"]; QC12 = QC12[QC12.usable & QC12.template & QC12.sta.isin(G.CN4)]
    fz = QC12.groupby("bkey").sta.apply(lambda s: ",".join(sorted(s)))
    fe = use.groupby("bkey").sta.apply(lambda s: ",".join(sorted(s)))
    C = pd.DataFrame(dict(fam=range(130)))
    C["frozen_cn4"] = C.fam.map(fz).fillna(""); C["era"] = C.fam.map(fe).fillna("")
    C["n_frozen_cn4"] = C.frozen_cn4.str.count(",") + (C.frozen_cn4 != ""); C["n_era"] = C.era.str.count(",") + (C.era != "")
    C["n_total"] = C.n_frozen_cn4 + C.n_era
    C["n_build_events"] = C.fam.map(cnt).fillna(0).astype(int)
    C.to_csv(f"{TAB}/g3_carriers.csv", index=False)
    G.log(f"[g3 tpl] usable era templates {len(TPLu)}; per station {use.groupby('sta').size().to_dict()}; "
          f"families with >= 1 era template {int((C.n_era > 0).sum())}; median carriers frozen CN4 {C.n_frozen_cn4.median()}, "
          f"CN4 + era {C.n_total.median()}; families with >= 2 sites CN4 {int((C.n_frozen_cn4 >= 2).sum())}, CN4+era {int((C.n_total >= 2).sum())}")
    print(Q.groupby("sta").agg(n=("bkey", "size"), usable=("usable", "sum"), snr=("snr_db", "median"), snr_null=("snr_null_db", "median"),
                               spick=("s_pick_s", "median"), epi=("epi_km", "median")).round(2).to_string())


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "check":      # sanity: the same delta on 2012 test blocks / 2008 should be ~0
        for w, yr in (("2012-09", 2012), ("2008-05", 2008), ("2004-07", 2004)):
            Dl, med = delta_pgc(E.events(w, test_blocks=(w == "2012-09")), yr, E.CACHES[yr])
            Dl.to_csv(f"{TAB}/g3_delta_pgc_check_{w}.csv", index=False)
            G.log(f"[g3 check] {w}: median family delta {med:+.3f} s, n {int(Dl.delta.notna().sum())}, "
                  f"family-delta MAD {(Dl.delta - med).abs().median():.3f}, median within-family MAD {Dl.delta_mad.median():.3f}")
    else:
        main()
