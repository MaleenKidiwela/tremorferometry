#!/usr/bin/env python
"""bw_scope_volume -- SVI backward scope, step 5: download / cache / GPU / storage estimates, 2003-01-01..2013-12-31.

Station-days per station-year = days inside the station's metadata epochs (inventory_channels.csv, velocity bands)
x the delivered fraction measured by the 10-min samples (availability_station_year.csv) = lower estimate;
the *_upper columns assume every metadata-epoch day exists (NRCan misses ~2/3 of single requests, so the sample
fraction under-states CN/PO availability).
MB per raw station-day = median of our own past downloads of that station (lfe_extension / lfe_picker / 2005 bed
manifests); fallback 10 MB at 40 Hz, 20 MB at 100 Hz.  Cache = 20.74 MB per station-day (measured, 3-C float32 20 Hz npz).
Days already on disk (ondisk_station_days_pre2014.csv) are subtracted from the download.
GPU = forward production matched filter (data + scrambled control + features): 4641 days in 19,115 s wall on 3 shards
(lfe_extension/results/refine/logs/prodmf_*.log) = 4.12 s per network-day.
Two station plans:
  A  frozen detector: the 15 template stations only (+ none)                      -> what the forward run would see
  B  A + every non-template station within 60 km that delivers (POLARIS/C8/LZB/PFB/YOUB/...) -> needed only if
     era templates are built
Out: volume_station_year.csv, volume_summary.csv
"""
import numpy as np, pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/backward_scope"
TPL = ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B001", "B003", "B004", "B005", "B006", "B007", "B009", "B010", "B011", "B926"]
CACHE_MB = 20.74
GPU_S_PER_DAY = 19115.0 / 4641


def mb_per_day():
    out = {}
    for p in [f"{ROOT}/lfe_extension/results/download_manifest.csv", f"{ROOT}/lfe_picker/results/download_manifest.csv",
              f"{ROOT}/lfe_array_discovery/results/step1c_download_manifest.csv"]:
        m = pd.read_csv(p, low_memory=False)
        m = m[m.status == "data"].copy(); m["bytes"] = pd.to_numeric(m.bytes, errors="coerce")
        for (s, fs), v in m.groupby(["sta", "fs"]).bytes.median().items():
            out.setdefault((s, float(fs)), v / 1e6)
    return out


def main():
    A = pd.read_csv(f"{OUT}/availability_station_year.csv")
    C = pd.read_csv(f"{OUT}/inventory_channels.csv")
    C = C[C.cha.str[:2].isin(["HH", "BH", "EH", "SH"]) & C.cha.str.endswith(("Z", "3"))]
    OD = pd.read_csv(f"{OUT}/ondisk_station_days_pre2014.csv", index_col=0)
    MB = mb_per_day()
    rows = []
    for r in A.itertuples():
        y0, y1 = pd.Timestamp(f"{r.year}-01-01"), pd.Timestamp(f"{r.year + 1}-01-01")
        c = C[(C.net == r.net) & (C.sta == r.sta)]
        if c.empty:
            continue
        # best band per day = highest fs epoch; count days in union of epochs
        days = set()
        fs_main = 0.0
        for e in c.itertuples():
            a = max(pd.Timestamp(e.start[:10]), y0); b = min(pd.Timestamp(e.end[:10]), y1)
            if b > a:
                days |= set(pd.date_range(a, b - pd.Timedelta(days=1)).date)
                fs_main = max(fs_main, float(e.fs))
        nd = len(days)
        frac = r.n_delivered / max(r.n_in_epoch, 1)
        exp_days = nd * frac
        mb = MB.get((r.sta, fs_main), 10.0 if fs_main <= 50 else 20.0)
        ondisk = float(OD.loc[f"{r.net}.{r.sta}", str(r.year)]) if (f"{r.net}.{r.sta}" in OD.index and str(r.year) in OD.columns) else 0.0
        rows.append(dict(net=r.net, sta=r.sta, year=r.year, template=r.sta in TPL, dist_med_km=r.dist_med_km,
                         epoch_days=nd, delivered_frac=round(frac, 2), exp_station_days=round(exp_days),
                         fs=fs_main, mb_per_day=round(mb, 1), ondisk_days=ondisk,
                         raw_GB=round(exp_days * mb / 1e3, 1), raw_GB_upper=round(nd * mb / 1e3, 1),
                         dl_GB_upper=round(max(nd - ondisk, 0) * mb / 1e3, 1), cache_GB_upper=round(nd * CACHE_MB / 1e3, 1), dl_GB=round(max(exp_days - ondisk, 0) * mb / 1e3, 1),
                         cache_GB=round(exp_days * CACHE_MB / 1e3, 1)))
    V = pd.DataFrame(rows)
    V.to_csv(f"{OUT}/volume_station_year.csv", index=False)
    near = V.dist_med_km <= 80          # non-template stations worth an era template (median family distance)
    summ = []
    for nm, sel in [("A_template_stations", V.template), ("B_plus_nontemplate_le80km_median", V.template | near)]:
        W = V[sel]
        summ.append(dict(plan=nm, n_stations=W.sta.nunique(), station_days=int(W.exp_station_days.sum()),
                         raw_TB=round(W.raw_GB.sum() / 1e3, 2), download_TB=round(W.dl_GB.sum() / 1e3, 2),
                         cache_TB=round(W.cache_GB.sum() / 1e3, 2), station_days_upper=int(W.epoch_days.sum()),
                         raw_TB_upper=round(W.raw_GB_upper.sum() / 1e3, 2), download_TB_upper=round(W.dl_GB_upper.sum() / 1e3, 2),
                         cache_TB_upper=round(W.cache_GB_upper.sum() / 1e3, 2)))
    S = pd.DataFrame(summ)
    days = (pd.Timestamp("2014-01-01") - pd.Timestamp("2003-01-01")).days
    S["network_days"] = days
    S["gpu_h_mf_data_plus_control"] = round(days * GPU_S_PER_DAY / 3600, 1)
    S.to_csv(f"{OUT}/volume_summary.csv", index=False)
    print(V.groupby(["template", "year"])[["exp_station_days", "raw_GB", "dl_GB", "cache_GB"]].sum().to_string())
    print(S.to_string())


if __name__ == "__main__":
    main()
