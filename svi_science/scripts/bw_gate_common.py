#!/usr/bin/env python
"""bw_gate_common -- shared definitions for the SVI BACKWARD-extension gates (G0 / G1 / G2).

Nothing of the forward pipeline is edited: the forward modules (lfe_extension/scripts/svi_ref_common.py,
c_common.py, svi_ref_mf.py) are imported and PATCHED IN MEMORY from here:
  * extra station sets (CN4 / B1 / B2 / B3) are added to svi_ref_common.SSETS;
  * when BW_EPOCHS=1, c_common.epochs() also returns the pre-2010-08 station epochs (forward table starts
    2010-08-01); c_common.manifest() then also returns the backward download ledger; c_common.trace_path()
    also looks in the backward trace root (for the PGC HF ratio).
"Scrambled control" = the same detector on the same day with each station site's data shifted in time.
"""
from __future__ import annotations
import os, sys, glob
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
LFE = f"{ROOT}/lfe_extension"
LFE_SCRIPTS = f"{LFE}/scripts"
REP = f"{ROOT}/lfe_catalog_repair"
SVI = f"{ROOT}/svi_science"
OUT = f"{SVI}/results/backward_gates"
SCOPE = f"{SVI}/results/backward_scope"
REFINE = f"{LFE}/results/refine"
TRACES_BW = f"{LFE}/traces_backward"            # new raw root for 2003-2013 (NET.STA/YYYY/DOY.mseed)
CACHE_BW = f"{LFE}/cache20hz_backward"          # 20 Hz / 1-8 Hz caches of the new raw (forward chain)
LEDGER = os.environ.get("BW_LEDGER", f"{OUT}/tables/download_ledger.csv")
# forward-chain 20 Hz caches already on disk (same e4_cache20 chain; read-only), searched after CACHE_BW
FWD_CACHES = [f"{LFE}/cache20hz_r1", f"{LFE}/cache20hz", f"{REP}/cache20hz", f"{ROOT}/lfe_picker/cache20hz"]
BOST = f"{REP}/results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv"

CN4 = ["PGC", "VGZ", "NLLB", "SNB"]
OLY = ["B001", "B003", "B004", "B005", "B006", "B007"]
PGCB = ["B009", "B010", "B011"]
B1 = CN4 + OLY
B2 = CN4 + OLY + PGCB
B3 = B2 + ["B926"]                                  # = forward S14
NETSTA = {"PGC": "CN.PGC", "VGZ": "CN.VGZ", "NLLB": "CN.NLLB", "SNB": "CN.SNB", "LRIV": "UW.LRIV",
          **{s: f"PB.{s}" for s in OLY + PGCB + ["B926"]}}

# pre-2010-08-01 band epochs of the 15 template stations (FDSN metadata, backward_scope/inventory_channels.csv);
# end EXCLUSIVE and capped at 2010-08-01, where results/station_epochs.csv (forward) takes over.
# VGZ before 2003-02-27 (EH S-13 short period) is left out on purpose (SCOPE section 3: not worth a template).
BW_EPOCHS = [
    ("CN", "PGC", "BH", 40, "2000-06-28", "2010-08-01"),
    ("CN", "NLLB", "BH", 40, "2003-01-01", "2010-08-01"),
    ("CN", "SNB", "BH", 40, "2000-05-25", "2010-08-01"),
    ("CN", "VGZ", "BH", 40, "2003-02-27", "2009-07-31"),
    ("CN", "VGZ", "HH", 100, "2009-07-31", "2010-08-01"),
    ("UW", "LRIV", "BH", 40, "2009-08-27", "2010-08-01"),
    ("PB", "B001", "EH", 100, "2005-06-30", "2010-08-01"),
    ("PB", "B003", "EH", 100, "2005-06-26", "2010-08-01"),
    ("PB", "B004", "EH", 100, "2005-06-15", "2010-08-01"),
    ("PB", "B005", "EH", 100, "2005-07-21", "2010-08-01"),
    ("PB", "B006", "EH", 100, "2005-07-28", "2010-08-01"),
    ("PB", "B007", "EH", 100, "2005-07-22", "2010-08-01"),
    ("PB", "B009", "EH", 100, "2005-11-02", "2010-08-01"),
    ("PB", "B010", "EH", 100, "2005-11-02", "2010-08-01"),
    ("PB", "B011", "EH", 100, "2005-11-02", "2010-08-01"),
    ("PB", "B926", "EH", 100, "2007-11-13", "2010-08-01"),
]
# G3: POLARIS + LZB era stations (FDSN metadata epochs, end exclusive; KELB -> KLNB rename 2004-06-03: KLNB only)
ERA = ["SHVB", "MGCB", "TWBB", "SHDB", "KHVB", "TWKB", "TSJB", "KLNB", "SILB", "SSIB", "TWGB", "LZB", "KELB"]   # KELB added 2026-10-09 (user): = KLNB site before 2004-06-03
BW_EPOCHS += [
    ("C8", "SHVB", "BH", 40, "2003-02-27", "2004-09-30"), ("C8", "MGCB", "HH", 100, "2003-08-12", "2006-05-13"),
    ("C8", "TWBB", "HH", 100, "2003-07-31", "2006-06-20"), ("C8", "SHDB", "BH", 40, "2003-02-28", "2005-01-11"),
    ("C8", "KHVB", "BH", 40, "2003-02-28", "2004-08-05"), ("PO", "TWKB", "HH", 100, "2003-08-01", "2006-05-14"),
    ("PO", "TSJB", "HH", 100, "2003-08-01", "2006-06-19"), ("PO", "KLNB", "HH", 100, "2004-06-04", "2006-06-22"),
    ("PO", "SILB", "HH", 100, "2003-08-03", "2006-07-04"), ("PO", "SSIB", "HH", 100, "2003-07-25", "2006-06-16"),
    ("PO", "TWGB", "HH", 100, "2003-08-02", "2005-05-10"), ("CN", "LZB", "BH", 40, "2001-01-10", "2006-07-01"),   # LZB used only as an era station (end capped 2026-10-09)
    ("PO", "KELB", "HH", 100, "2002-09-12", "2004-06-04")]   # same site as KLNB (48.661 -123.570, CMG3ESP/Trident; gain 6.17e7 to 2003-07-31, 3.07e8 after; KLNB 3.97e8)
NETSTA.update({s_: f"{n_}.{s_}" for n_, s_, *_ in BW_EPOCHS})
# instrument epochs inside a band (sensor/digitiser), for G2 and the tpl_epoch column
INSTR = {
    "VGZ": [("2003-02-27", "2007-02-28", "BH CMG3ESP/GD1"), ("2007-02-28", "2009-07-31", "BH CMG3T/GD2_v21"),
            ("2009-07-31", "2017-03-03", "HH CMG3ESP/GD2_v21 (2012 build epoch)")],
    "SNB": [("2000-05-25", "2006-04-20", "BH CMG40T/GD1"), ("2006-04-20", "2008-12-11", "BH CMG3ESP/GD2_v12"),
            ("2008-12-11", "2011-03-17", "BH CMG3ESP/GD1"), ("2011-03-17", "2017-06-08", "HH CMG3ESP/GD2_v21 (2012 build epoch)")],
}


def log(m):
    import time
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def epochs_table():
    rows = []
    for net, sta, band, fs, a, b in BW_EPOCHS:
        lat = lon = np.nan
        rows.append(dict(net=net, sta=sta, band=band, fs=fs, loc="", start=a, end=b, lat=lat, lon=lon, n_comp=3))
    B = pd.DataFrame(rows)
    F = pd.read_csv(f"{LFE}/results/station_epochs.csv")
    ll = F.groupby("sta")[["lat", "lon"]].first()
    I = pd.read_csv(f"{SCOPE}/inventory_stations.csv").groupby("sta")[["lat", "lon"]].first()
    B["lat"] = B.sta.map(ll.lat).fillna(B.sta.map(I.lat)); B["lon"] = B.sta.map(ll.lon).fillna(B.sta.map(I.lon))
    return pd.concat([B, F], ignore_index=True)


def patch_forward(backward=False):
    """import the forward modules, add the station sets; with backward=True also extend epochs / manifest / traces."""
    if LFE_SCRIPTS not in sys.path:
        sys.path.insert(0, LFE_SCRIPTS)
    import svi_ref_common as R
    import c_common as C
    R.SSETS.update(CN4=list(CN4), B1=list(B1), B2=list(B2), B3=list(B3))
    for base, mem in (("CN4", CN4), ("B1", B1), ("B3", B3)):          # withheld-station sets for G2 (A-2 protocol)
        for x in CN4:
            R.SSETS[f"{base}no{x}"] = [s_ for s_ in mem if s_ != x]
    R.SSETS["CN4ERA"] = list(CN4) + list(ERA)
    R.SSETS["S15ERA"] = list(R.S15) + list(ERA)            # production era-augmented detector (frozen network + era)
    if os.environ.get("BW_SITELIST_EXT"):                  # G3: era stations appended to the scrambled-control site list
        R.SITE_LIST = list(R.SITE_LIST) + [s_ for s_ in sorted(ERA) if s_ not in R.SITE_LIST]
    era = os.environ.get("BW_ERA_TPL", "")
    if era:                                                 # G3: merge the era templates / weights / amplitude refs
        _lft = C.load_frozen_tables

        def load_frozen_tables(fz, tpl_path, L):
            tpl, sp, F, wgt, mv, ampref = _lft(fz, tpl_path, L)
            z = np.load(era, allow_pickle=True)
            et, es = z["tpl"].item(), z["sp"].item()
            for k, v in et.items():
                if L in v:
                    tpl[k] = v; sp[k] = es[k]
            ov = os.environ.get("BW_TPL_OVERRIDE", "")          # P3: e.g. SNB era templates replacing the 2012 ones
            if ov:
                zo = np.load(ov, allow_pickle=True); ot, osp = zo["tpl"].item(), zo["sp"].item()
                for k, v in ot.items():
                    if L in v:
                        tpl[k] = v; sp[k] = osp[k]
                Wo = ov.replace(".npz", "_weights.csv")
                if os.path.exists(Wo):
                    wgt.update({(float(r.fam), str(r.sta)): float(r.w) for r in pd.read_csv(Wo).itertuples()})
            W = pd.read_csv(era.replace("g3_era_templates.npz", "g3_era_weights.csv"))
            wgt.update({(float(r.fam), str(r.sta)): float(r.w) for r in W.itertuples()})
            A = pd.read_csv(era.replace("g3_era_templates.npz", "g3_era_ampref.csv"))
            ampref.update({(str(r.sta), float(r.L)): float(r.e3_ref) for r in A.itertuples()})
            return tpl, sp, F, wgt, mv, ampref
        C.load_frozen_tables = load_frozen_tables
        _st = C.station_table

        def station_table():
            s0 = _st()
            E = epochs_table().drop_duplicates("sta").set_index("sta")
            add = pd.DataFrame([dict(sta=s_, net=E.loc[s_, "net"], netsta=NETSTA[s_], lat=E.loc[s_, "lat"], lon=E.loc[s_, "lon"])
                                for s_ in ERA if s_ not in s0.index]).set_index("sta")
            return pd.concat([s0, add])
        C.station_table = station_table
    if backward:
        def epochs():
            e = epochs_table()
            e["start"] = pd.to_datetime(e.start); e["end"] = pd.to_datetime(e.end)
            e["netsta"] = e.net + "." + e.sta
            e["epoch_id"] = e.netsta + ":" + e.band + ":" + e.start.dt.strftime("%Y-%m-%d")
            return e
        C.epochs = epochs
        _man = C.manifest

        def manifest():
            m = _man()
            lfs = [x for x in os.environ.get("BW_LEDGERS", LEDGER).split(",") if x and os.path.exists(x)]
            if lfs:
                L = pd.concat([pd.read_csv(x, dtype=str) for x in lfs], ignore_index=True).fillna("")
                L = L[L.status.isin(["data", "partial", "nodata", "error"])]
                L = L.drop_duplicates(["net", "sta", "date"], keep="last")
                add = pd.DataFrame(dict(net=L.net, sta=L.sta, band=L.band, date=L.date,
                                        doy=pd.to_datetime(L.date).dt.dayofyear.astype("Int64"),
                                        status=np.where(L.status == "data", "data", "nodata"),
                                        coverage_frac=pd.to_numeric(L.coverage_frac, errors="coerce"), note=""))
                add["netsta"] = add.net + "." + add.sta
                m = pd.concat([m, add], ignore_index=True)
            return m
        C.manifest = manifest
        _tp = C.trace_path

        def trace_path(netsta, year, doy, note=""):
            p = f"{TRACES_BW}/{netsta}/{year}/{doy:03d}.mseed"
            if os.path.exists(p) and os.path.getsize(p) > 0:
                return p
            return _tp(netsta, year, doy, note)
        C.trace_path = trace_path
    return R, C


def bostock_all():
    """Bostock (2015) events of ALL years mapped to the 2012 family keys, on the detector's origin axis (frozen 2012
    per-family offsets) -- svi_ref_common.bostock() logic without its sealed-year filter (2013 is never read here)."""
    R, C = patch_forward()
    FAM = pd.read_csv(R.FAMS12); OFF = pd.read_csv(R.OFFS12)
    omap = {int(r.bkey): float(r.offset_s) for r in OFF.itertuples() if bool(r.accept)}
    B = pd.read_csv(BOST)
    B = B[B.year != 2013]
    B["time"] = pd.to_datetime(dict(year=B.year, month=B.month, day=B.day, hour=B.hour, minute=B["min"],
                                    second=B.sec.astype(int)), errors="coerce") + pd.to_timedelta(B.sec % 1.0, unit="s")
    k = {(round(r.lat, 4), round(r.lon, 4), round(r.dep, 4)): int(r.bkey) for r in FAM.itertuples()}
    B["fam"] = [k.get((round(a, 4), round(b, 4), round(c, 4)), -1) for a, b, c in zip(B.lat, B.lon, B.dep)]
    B = B[B.fam >= 0].copy()
    B["t_origin"] = R.to_epoch(B.time) + B.fam.map(omap)
    B = B.dropna(subset=["t_origin"]).sort_values("t_origin").reset_index(drop=True)
    B["date"] = pd.to_datetime(B.t_origin, unit="s").dt.strftime("%Y-%m-%d")
    B["year"] = B.date.str[:4].astype(int)
    return B[["t_origin", "fam", "date", "year", "mag", "lat", "lon", "dep"]]
