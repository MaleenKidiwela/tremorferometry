#!/usr/bin/env python
"""bw_gate_g0_score -- G0 pilot verdict: download completeness per station (CN from NRCan, PB from EarthScope), request
statistics (the NRCan error mode), sizes, cache status, and the end-to-end detector pass on the pilot days
(guard + candidates per day, data vs scrambled control, ETS vs quiet).  The 2008 detector recall is NOT scored
(2003-2009 truth stays unspent for the one look per network model, SCOPE R7).

Out: OUT/tables/g0_download_by_station.csv, g0_requests.csv, g0_detector_days.csv, g0_decision.json, figures/g0_pilot.png
"""
from __future__ import annotations
import os, sys, glob, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_gate_common as G

TAB = f"{G.OUT}/tables"
PILOT = [str(x.date()) for x in pd.date_range("2008-05-01", "2008-06-04")] + [str(x.date()) for x in pd.date_range("2008-08-01", "2008-08-31")]
ETS = ("2008-05-10", "2008-06-02")


def epoch_days(sta, days):
    out = []
    for d in days:
        for net, s, band, fs, a, b in G.BW_EPOCHS:
            if s == sta and a <= d < b:
                out.append(d); break
    return out


def main():
    L = pd.read_csv(G.LEDGER, dtype=str).fillna("")
    L = L[~L.tag.isin(["test", "cachecheck"])]
    first = L[(L.chan_log != "resumed") & ~L.tag.str.startswith("repair")].drop_duplicates(["net", "sta", "date"], keep="first")
    rep = L[L.chan_log.str.startswith("repair")]
    L = L.drop_duplicates(["net", "sta", "date"], keep="last")          # final state on disk (resumed / repair rows verify it)
    win = lambda d: ("pilot_ets" if "2008-05-01" <= d <= "2008-06-04" else "pilot_quiet" if d.startswith("2008-08")
                     else "g2_2004" if d.startswith("2004") else "g2_2005")
    L["tag"] = L.date.map(win)
    first = first.assign(requests=pd.to_numeric(first.requests, errors="coerce"))
    L["covf"] = pd.to_numeric(L.coverage_frac, errors="coerce").fillna(0)
    L["nch"] = pd.to_numeric(L.n_channels, errors="coerce").fillna(0)
    L["ok3"] = (L.status == "data") & (L.nch >= 3) & (L.covf >= 0.5)
    L["requests"] = pd.to_numeric(L.requests, errors="coerce")
    L["bytes"] = pd.to_numeric(L.bytes, errors="coerce").fillna(0)
    rows = []
    for tag, g0 in L.groupby("tag"):
        for (net, sta), g in g0.groupby(["net", "sta"]):
            days = sorted(g.date)
            ep = epoch_days(sta, days)
            rows.append(dict(tag=tag, net=net, sta=sta, epoch_days=len(ep), rows=len(g), ok3=int(g.ok3.sum()),
                             frac_ok3=g.ok3.sum() / max(len(ep), 1), partial=int((g.status == "partial").sum()),
                             nodata=int((g.status == "nodata").sum()), error=int((g.status == "error").sum()),
                             GB=g.bytes.sum() / 1e9,
                             req_per_day=first[(first.sta == sta) & first.date.isin(days)].requests.mean(),
                             secs_per_day=pd.to_numeric(first[(first.sta == sta) & first.date.isin(days)].secs, errors="coerce").mean()))
    S = pd.DataFrame(rows); S.to_csv(f"{TAB}/g0_download_by_station.csv", index=False)
    # NRCan per-request success rate (from the chan_log of the first attempts)
    rq = []
    for r in first[first.net == "CN"].itertuples():
        for tok in str(r.chan_log).split():
            if ":ok@" in tok:
                k = int(tok.split("@")[1]); rq.append(dict(sta=r.sta, date=r.date, tries=k, ok=True))
            elif ":fail@" in tok:
                rq.append(dict(sta=r.sta, date=r.date, tries=int(tok.split("@")[1].split("(")[0]), ok=False))
    Q = pd.DataFrame(rq); Q.to_csv(f"{TAB}/g0_requests.csv", index=False)
    p_req = float(Q.ok.sum() / Q.tries.sum()) if len(Q) else np.nan
    rep_ok = int(rep.chan_log.str.contains(":ok@").sum()); rep_n = len(rep)
    pil = S[S.tag.isin(["pilot_ets", "pilot_quiet"])]
    cn = pil[pil.net == "CN"]
    cn_frac = float(cn.ok3.sum() / max(cn.epoch_days.sum(), 1))
    # detector
    fd = sorted(glob.glob(f"{G.OUT}/mf/feat/g0_pilot_B3/*.json"))
    drows = []
    for f in fd:
        m = json.load(open(f))
        d = m["date"]
        p = f"{G.OUT}/mf/feat/g0_pilot_B3/{d}.parquet"
        n8 = dict(data=0, c1=0)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            try:
                A = pd.read_parquet(p, columns=["series", "net_mad"])
                for s in ("data", "c1"):
                    n8[s] = int(((A.series == s) & (A.net_mad >= 8.75)).sum())
            except Exception:
                pass
        drows.append(dict(date=d, n_sta=len(m.get("stations", [])), stations=",".join(m.get("stations", [])),
                          n_data_ge8=m.get("n_data", 0), n_c1_ge8=m.get("n_c1", 0), n_data_ge875=n8["data"], n_c1_ge875=n8["c1"],
                          ets=ETS[0] <= d <= ETS[1], secs=m.get("secs")))
    D = pd.DataFrame(drows); D.to_csv(f"{TAB}/g0_detector_days.csv", index=False)
    st = pd.concat([pd.read_csv(f, dtype=str) for f in glob.glob(f"{G.OUT}/mf/status/mf_g0_pilot_B3_s*.csv")], ignore_index=True) \
        if glob.glob(f"{G.OUT}/mf/status/mf_g0_pilot_B3_s*.csv") else pd.DataFrame(columns=["status"])
    st = st.drop_duplicates(["date", "netsta"], keep="last") if len(st) else st
    nbad = int(st.status.isin(["error", "missing"]).sum()) if len(st) else -1
    nok = int((st.status == "ok").sum()) if len(st) else 0
    err_frac = nbad / max(nok + nbad, 1) if nbad >= 0 else np.nan
    dec = dict(cn_pilot_ok3_station_days=int(cn.ok3.sum()), cn_pilot_epoch_days=int(cn.epoch_days.sum()), cn_frac_ok3=cn_frac,
               cn_by_station=cn[["tag", "sta", "ok3", "epoch_days", "frac_ok3"]].to_dict("records"),
               nrcan_per_request_success=p_req, nrcan_requests_per_station_day=float(first[first.net == "CN"].requests.mean()),
               nrcan_channel_fails_after_20=int((~Q.ok).sum()) if len(Q) else 0, repair_days=rep_n,
               repair_days_recovered=int((rep.n_channels.astype(str) == "3").sum()), repair_log=rep[["sta", "date", "chan_log"]].values.tolist(),
               cn_secs_per_station_day=float(pd.to_numeric(first[first.net == "CN"].secs, errors="coerce").mean()),
               pb_pilot_frac_ok3=float(pil[pil.net == "PB"].ok3.sum() / max(pil[pil.net == "PB"].epoch_days.sum(), 1)),
               GB_all=float(L.bytes.sum() / 1e9), GB_by_tag=L.groupby("tag").bytes.sum().div(1e9).round(2).to_dict(),
               detector_days=len(D), detector_days_expected=len(PILOT), detector_station_day_err_frac=err_frac,
               detector_status=st.status.value_counts().to_dict() if len(st) else {},
               ets_data_per_day_ge875=float(D[D.ets].n_data_ge875.mean()) if len(D) else np.nan,
               ets_c1_per_day_ge875=float(D[D.ets].n_c1_ge875.mean()) if len(D) else np.nan,
               quiet_data_per_day_ge875=float(D[D.date.str.startswith("2008-08")].n_data_ge875.mean()) if len(D) else np.nan,
               quiet_c1_per_day_ge875=float(D[D.date.str.startswith("2008-08")].n_c1_ge875.mean()) if len(D) else np.nan)
    dec["pass_download"] = bool(cn_frac >= 0.95)
    dec["pass_end_to_end"] = bool(len(D) == len(PILOT) and np.isfinite(err_frac) and err_frac <= 0.03)
    dec["pass"] = dec["pass_download"] and dec["pass_end_to_end"]
    json.dump(dec, open(f"{TAB}/g0_decision.json", "w"), indent=1, default=float)
    G.log(json.dumps({k: v for k, v in dec.items() if k != "cn_by_station"}, default=float))
    if len(D):
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=False)
        x = pd.to_datetime(D.date)
        for k, (lo, hi) in enumerate((("2008-05-01", "2008-06-04"), ("2008-08-01", "2008-08-31"))):
            m = (D.date >= lo) & (D.date <= hi)
            ax[k].bar(x[m] - pd.Timedelta(hours=4), D.n_data_ge875[m], width=0.4, color="#2471a3", label="data, >= 8.75 MAD")
            ax[k].bar(x[m] + pd.Timedelta(hours=4), D.n_c1_ge875[m], width=0.4, color="0.6", label="scrambled control, >= 8.75 MAD")
            ax[k].set_ylabel("candidates / day"); ax[k].legend(fontsize=8)
            if k == 0:
                ax[k].axvspan(pd.Timestamp(ETS[0]), pd.Timestamp(ETS[1]) + pd.Timedelta(days=1), color="#f5cba7", alpha=.4, zorder=0)
                ax[k].set_title("G0 pilot, B3 network (= forward S14): 2008 ETS (shaded = Bostock episode) and quiet 2008-08")
        plt.tight_layout(); plt.savefig(f"{G.OUT}/figures/g0_pilot.png", dpi=130); plt.close()


if __name__ == "__main__":
    main()
