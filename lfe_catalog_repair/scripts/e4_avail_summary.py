#!/usr/bin/env python
"""E4 step 6 summary — who was recording 3-C at >= 40 Hz, and what would each download cost?

Built from results/e4_station_epochs_full.csv (an unrestricted FDSN channel-metadata query over
1990-2026 for every station of interest).  Station CODES collide across networks -- TWGB and TWKB
also exist in the Taiwan network TW, LCBC in Colombia's CM, MGCB in ZR (Chile), OCP in CI and UU,
B011 in SN and ZF -- so only the CASCADIA network codes (PO, C8, CN, PB, UW) and a coordinate box
around the margin are accepted. Without that filter the query reports POLARIS stations as "open".
"""
from __future__ import annotations
import numpy as np, pandas as pd

RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
TEMPLATE21 = ["GLBC", "GOWB", "JRBC", "KELB", "KHVB", "LCBC", "LZB", "MGCB", "PFB", "PGC",
              "SHVB", "SILB", "SNB", "SOKB", "SSIB", "TSJB", "TWBB", "TWGB", "TWKB", "VGZ", "YOUB"]
CURRENT = ["PGC", "LZB", "VGZ", "SNB", "YOUB", "PFB", "SOKB", "NLLB", "B011", "B926", "LRIV",
           "OCP", "STW", "CLRS", "SYMB", "PTRF"]
NETS_OK = {"PO", "C8", "CN", "PB", "UW", "UO"}
WIN = [("2010ETS", "2010-08-01", "2010-09-30"), ("2012ETS", "2012-08-15", "2012-10-15"),
       ("2019", "2019-01-01", "2019-12-31"), ("2021", "2021-01-01", "2021-12-31"),
       ("2024", "2024-08-01", "2024-09-30"), ("2026", "2026-01-01", "2026-09-01")]


def main():
    D = pd.read_csv(f"{RES}/e4_station_epochs_full.csv", keep_default_na=False)
    D["fs"] = pd.to_numeric(D.fs, errors="coerce")
    D = D[D.net.isin(NETS_OK) & (D.lat > 39) & (D.lat < 51) & (D.lon > -128) & (D.lon < -119)]
    D = D[(D.fs >= 40) & D.chan.str[-1].isin(list("ZNE12"))].copy()
    D["e"] = D.end.replace("open", "2030-01-01")
    D["band"] = D.chan.str[:2]
    rows = []
    for sta, h in D.groupby("sta"):
        r = dict(sta=sta, nets=";".join(sorted(h.net.unique())), first=h.start.min(),
                 last=h.e.max(), lat=float(h.lat.iloc[0]), lon=float(h.lon.iloc[0]))
        for lab, t0, t1 in WIN:
            z = h[(h.start <= t1) & (h.e >= t0)]
            n3 = 0
            if len(z):
                cnt = z.groupby(["net", "loc", "band"]).chan.nunique()
                n3 = int((cnt >= 3).any())
            r[lab] = n3
            r[lab + "_1c"] = int(len(z) > 0)
        rows.append(r)
    G = pd.DataFrame(rows).sort_values("sta")
    G.to_csv(f"{RES}/e4_station_availability.csv", index=False)
    pd.set_option("display.width", 220)
    print(G[["sta", "nets", "first", "last"] + [w[0] for w in WIN]].to_string(index=False))
    for roster, lst in (("21 SeisBench template stations", TEMPLATE21),
                        ("current ring + boreholes", CURRENT)):
        print(f"\n### {roster}")
        for lab, _, _ in WIN:
            h = G[G.sta.isin(lst)]
            on = sorted(h[h[lab] == 1].sta)
            print(f"  {lab}: {len(on)}/{len(lst)} with 3-C >= 40 Hz -> {', '.join(on)}")

    def gb(nsta, days, fs, nc=3):
        return nsta * days * nc * (0.42 * 24 * fs / 100.0) / 1024.0
    print("\n=== candidate downloads (NOT executed; volumes at the native rate and at 40 Hz) ===")
    n10 = int(G[G.sta.isin(TEMPLATE21)]["2010ETS"].sum())
    n12 = int(G[G.sta.isin(TEMPLATE21)]["2012ETS"].sum())
    print(f"(a) the 21 template stations, 2010-08-01..09-30 (61 d): only {n10} are 3-C >= 40 Hz "
          f"-> {gb(n10, 61, 100):.1f} GB native / {gb(n10, 61, 40):.1f} GB at 40 Hz, and they are "
          f"EXACTLY the outer ring we already have on disk (+ SOKB).")
    print(f"(b) southern VI, 2012-08-15..10-15 (62 d): {n12} template stations "
          f"({gb(n12, 62, 100):.1f} / {gb(n12, 62, 40):.1f} GB) plus PB B011/B926, CN NLLB and "
          f"UW LRIV = 11 stations -> {gb(11, 62, 100):.1f} GB native, {gb(11, 62, 40):.1f} GB at "
          f"40 Hz.")
    print(f"(c) POLARIS Aug-Sep 2010: IMPOSSIBLE. Every inner POLARIS/C8 station "
          f"(GLBC GOWB JRBC KELB KHVB LCBC MGCB SHVB SILB SSIB TSJB TWBB TWGB TWKB) closed between "
          f"2004-06 and 2007-02; none was recording in 2010 or 2012. Their data exist only for "
          f"~2002-2007, which covers Muenchmeyer's 2003-2005 windows, not ours.")
    print(f"[write] {RES}/e4_station_availability.csv")


if __name__ == "__main__":
    main()
