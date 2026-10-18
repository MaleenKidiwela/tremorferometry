#!/usr/bin/env python
"""P3D wrapper around the FROZEN matched filter (`o_run_month.py`), adding the new stations.

`o_run_month.py` and `o_common.py` are NOT edited: every change below is a monkey-patch applied in
this process only and recorded in the run json by the tag.  `o_run_month` derives its station list
from the TEMPLATE BANK (`STAS = C.order_stations({k[1] for k in tpl})`), so a bank that contains
transferred (family, new-station) templates brings the new stations into the network automatically
once the supporting tables know about them:

  station_table   + the ten new stations and the Z-only helpers (lat/lon for the placement)
  epochs          + results/p3d/station_epochs.csv (netsta, band, alive window)
  SITES           + the CO-LOCATED merges B926<-CLRS and B001<-OSQM (Merlin flaw 3): without this
                  the two co-located pairs would count as four independent sites in the site sum
  weights         + the N-matched transfer weights from the bank's key table
  amp_ref         + a P3D quiet-month amplitude reference for the new stations (rel_amp only)
  epoch policy    optional --undrop for a station the frozen policy excludes (LRIV after 2019-05-30)

  python scripts/s_run_month.py 2024-08 --bank results/p3d/bank_transfer_v1/bank_production.npz \
      --drop-sites B926q,GOBB,LOPEZ,TURTL --tag dense_nowithheld --out-dir results/p3d/trackS
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402
import o_common as C                                                       # noqa: E402
import o_run_month as R                                                    # noqa: E402

_DROP = set()
_BANK = ""
_UNDROP = set()
_orig_load = C.load_frozen_tables
_orig_ep = C.epochs
_orig_tbl = C.station_table
_orig_pol = C.epoch_policy
_orig_share = C.share_ok_matrix


def _station_table():
    t = _orig_tbl()
    n = pd.read_csv(f"{S.P3D}/stations_p3d.csv").set_index("sta")
    n = n[~n.index.isin(t.index)]
    add = pd.DataFrame(dict(net=n.net, netsta=n.netsta, lat=n.lat, lon=n.lon, elev_m=n.elev_m,
                            chan=n.band.astype(str) + "?", fs_native=n.fs_native, klass=n.klass),
                       index=n.index)
    return pd.concat([t, add])


def _epochs():
    e = _orig_ep()
    p = pd.read_csv(f"{S.P3D}/station_epochs.csv")
    p = p[["net", "sta", "band", "fs", "loc", "start", "end", "lat", "lon", "n_comp"]]
    p["start"] = pd.to_datetime(p.start)
    p["end"] = pd.to_datetime(p.end)
    p["netsta"] = p.net + "." + p.sta
    p["epoch_id"] = p.netsta + ":" + p.band + ":" + p.start.dt.strftime("%Y-%m-%d")
    p = p[~p.sta.isin(set(e.sta))]
    return pd.concat([e, p], ignore_index=True)


def _policy():
    P = dict(_orig_pol())
    for k in list(P):
        if k[0] in _UNDROP:
            P[k] = dict(P[k], drop=False, label=f"p3d_transfer_v1", npz="")
    return P


def _ampref_extra(L):
    p = f"{S.P3D}/amp_ref_p3d.csv"
    if not os.path.exists(p):
        return {}
    A = pd.read_csv(p)
    return {(str(r.sta), float(r.L)): float(r.e3_ref) for r in A.itertuples()}


def _load(fz, tpl_path, L):
    tpl, sp, F, wgt, mv, ampref = _orig_load(fz, tpl_path, L)
    ampref.update(_ampref_extra(L))
    if _BANK:
        z = np.load(_BANK, allow_pickle=True)
        bt = z["tpl"].item()
        bs = z["sp"].item()
        n_add = 0
        for k, v in bt.items():
            if L in v:
                tpl[k] = v
                sp[k] = bs[k]
                n_add += 1
        kf = _BANK.replace(".npz", "_keys.csv")
        if os.path.exists(kf):
            K = pd.read_csv(kf)
            for r in K.itertuples():
                wgt[(float(r.bkey), str(r.sta))] = float(r.w_new)
        print(f"[s_run_month] bank {os.path.basename(_BANK)}: +{n_add} (family, station) templates",
              flush=True)
    if _DROP:
        drop_sta = set()
        ST = pd.read_csv(f"{S.P3D}/stations_p3d.csv").set_index("sta")
        for s in ST.index:
            if ST.site[s] in _DROP or s in _DROP:
                drop_sta.add(s)
        tpl = {k: v for k, v in tpl.items() if k[1] not in drop_sta}
        sp = {k: v for k, v in sp.items() if k[1] not in drop_sta}
        wgt = {k: v for k, v in wgt.items() if k[1] not in drop_sta}
        mv = {k: v for k, v in mv.items() if k[1] not in drop_sta}
        print(f"[s_run_month] withheld {sorted(_DROP)} -> stations {sorted(drop_sta)}; "
              f"{len(tpl)} (family, station) templates remain", flush=True)
    return tpl, sp, F, wgt, mv, ampref


def _share(stas=None):
    return _orig_share(stas)


def main():
    global _BANK
    argv = sys.argv[1:]
    out = []
    i = 0
    while i < len(argv):
        if argv[i] == "--drop-sites":
            _DROP.update(x for x in argv[i + 1].split(",") if x)
            i += 2
        elif argv[i] == "--bank":
            _BANK = argv[i + 1]
            i += 2
        elif argv[i] == "--undrop":
            _UNDROP.update(x for x in argv[i + 1].split(",") if x)
            i += 2
        else:
            out.append(argv[i])
            i += 1
    # the co-located merge, applied to the frozen site map
    C.SITES = dict(C.SITES)
    C.SITES["B926"] = ["B926", "CLRS"]
    C.SITES["B001"] = ["B001", "OSQM"]
    C.STAS_ORDER = S.STAS_FROZEN + S.NEW_3C
    C.load_frozen_tables = _load
    C.epochs = _epochs
    C.station_table = _station_table
    C.epoch_policy = _policy
    C.share_ok_matrix = _share
    for k, v in (("load_frozen_tables", _load), ("epochs", _epochs), ("station_table", _station_table),
                 ("epoch_policy", _policy), ("share_ok_matrix", _share)):
        setattr(R.C, k, v)
    R.C.SITES = C.SITES
    R.C.STAS_ORDER = C.STAS_ORDER
    sys.argv = [sys.argv[0]] + out
    R.main()


if __name__ == "__main__":
    main()
